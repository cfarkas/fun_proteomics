#!/usr/bin/env python3
"""Rebuild DIA delivery from verified completed searches, with all-pairs models."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import sqlite3
import subprocess
import sys

import pandas as pd
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'pipeline_lib'), str(ROOT/'modules/diapasef')]
from common import filter_report, load_report, read_method
from matrix_catalog import export_matrix_formats
from quantification import file_provenance, load_matrix
from quantification_config import differential_options, validate_design_rows
from quantification_analysis import write_dashboard


def stamp():
    return datetime.now(timezone.utc).isoformat()


def same_file_content(left, right):
    a, b = file_provenance(left), file_provenance(right)
    if a['sha256'] != b['sha256'] or a['bytes'] != b['bytes']:
        raise ValueError(f'Integrity mismatch: {left} differs from {right}')
    return a


def cached_report(report):
    manifest = report.with_name('report.manifest.txt')
    if not manifest.is_file():
        raise ValueError('Missing engine manifest; cannot verify completed search: ' + str(report))
    candidates = []
    for line in manifest.read_text().splitlines():
        for label, value in json.loads(line).items():
            if 'Main' in label and '.parquet' in label:
                candidates.append(Path(value))
    if len(candidates) != 1 or not candidates[0].is_file():
        raise ValueError('Original cached report is missing or ambiguous: ' + str(manifest))
    original = candidates[0].resolve()
    task = original.parent.parent
    if not (task/'.exitcode').is_file() or (task/'.exitcode').read_text().strip() != '0':
        raise ValueError('Cached search did not complete successfully: ' + str(task))
    log = report.with_name('engine.log')
    if not log.is_file() or 'Finished' not in log.read_text():
        raise ValueError('Engine completion evidence missing: ' + str(report))
    same_file_content(report, original)
    return task, original


def audit_report(report, samples, method, fasta):
    task, original = cached_report(report)
    method_record = same_file_content(method, task/'method.toml')
    fasta_record = same_file_content(fasta, task/'search.fasta')
    cached_samples = pd.read_csv(task/'samples.tsv', sep='\t', dtype=str, keep_default_na=False)
    keys = ['sample_id', 'sample_type', 'run']
    if set(map(tuple, samples[keys].values)) != set(map(tuple, cached_samples[keys].values)):
        raise ValueError('Requested samples differ from completed search: ' + str(report))
    # Decode every Parquet column/row group, not just its footer.
    parquet = pq.ParquetFile(report)
    decoded = sum(batch.num_rows for batch in parquet.iter_batches(batch_size=32768))
    if not decoded or decoded != parquet.metadata.num_rows:
        raise ValueError('Incomplete or empty Parquet report: ' + str(report))
    frame = load_report(report)
    if set(frame.Run.dropna()) != set(samples.run):
        raise ValueError('Report runs do not match requested samples: ' + str(report))
    filtered = filter_report(frame, read_method(method))
    if filtered.empty:
        raise ValueError('No positive protein evidence passes the configured confidence filters')
    return dict(report=file_provenance(report), cached_report=file_provenance(original),
                task=str(task), method=method_record, fasta=fasta_record,
                decoded_rows=decoded, passing_precursor_rows=len(filtered), samples=len(samples),
                checks=['successful task exit', 'engine completion', 'published/cache SHA256 equality',
                        'method and FASTA SHA256 equality', 'sample/run identity', 'complete Parquet decoding',
                        'required confidence fields and positive quantities'])


def audit_raw(samples, hash_raw=False):
    records = []
    for row in samples.itertuples(index=False):
        raw = Path(row.raw_path).resolve()
        if raw.suffix.lower() != '.d':
            raise ValueError('Expected a Bruker acquisition: ' + str(raw))
        files = []
        for name in ('analysis.tdf', 'analysis.tdf_bin'):
            path = raw/name
            if not path.is_file() or not path.stat().st_size:
                raise ValueError('Missing/empty acquisition component: ' + str(path))
            before = path.stat()
            record = file_provenance(path) if hash_raw else dict(path=str(path), bytes=before.st_size)
            after = path.stat()
            if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
                raise ValueError('Acquisition changed during integrity check: ' + str(path))
            record['mtime_ns'] = after.st_mtime_ns
            files.append(record)
        with sqlite3.connect((raw/'analysis.tdf').as_uri()+'?mode=ro', uri=True) as db:
            if db.execute('PRAGMA quick_check').fetchall() != [('ok',)]:
                raise ValueError('Bruker metadata SQLite integrity failed: ' + str(raw))
            types = {value[0] for value in db.execute('SELECT DISTINCT MsMsType FROM Frames')}
            if types != {0, 9}:
                raise ValueError('Not a DIA-PASEF acquisition: ' + str(raw))
        records.append(dict(sample_id=row.sample_id, files=files, sqlite_quick_check='ok', frame_types=sorted(types)))
        print(f'RAW VERIFIED {len(records)}/{len(samples)}: {row.sample_id}', flush=True)
    return records


def write_output_manifest(output):
    rows = [file_provenance(path) for path in sorted(output.rglob('*')) if path.is_file() and path.name != 'OUTPUT_SHA256SUMS']
    (output/'OUTPUT_SHA256SUMS').write_text(''.join(row['sha256']+'  '+Path(row['path']).relative_to(output).as_posix()+'\n' for row in rows))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-run', required=True, help='Completed run containing prepared/ and search/ (or analysis_details/search/)')
    parser.add_argument('--outdir', required=True, help='New empty result directory')
    parser.add_argument('--samples', help='Prepared TSV with sample IDs, conditions, raw paths and run names')
    parser.add_argument('--method', help='Defaults to this repository method.toml; must match the completed task')
    parser.add_argument('--service-matrix', action='append', default=[], metavar='COHORT=PATH')
    parser.add_argument('--independent-samples', action='store_true', help='Explicitly confirm one independent biological sample per row; use sample IDs as replicate IDs')
    parser.add_argument('--hash-raw', action='store_true', help='SHA256 every raw component in addition to SQLite/size checks')
    parser.add_argument('--fdr', type=float, default=.05)
    parser.add_argument('--min-replicates', type=int, default=2)
    parser.add_argument('--rscript', default=os.environ.get('PROTEOMICS_RSCRIPT', 'Rscript'))
    args = parser.parse_args()
    source, output = Path(args.source_run).resolve(), Path(args.outdir).resolve()
    if output == source or (output.exists() and any(output.iterdir())):
        parser.error('Use a new or empty output directory; completed results are not overwritten')
    sheet = Path(args.samples).resolve() if args.samples else source/'prepared/samples.tsv'
    samples = pd.read_csv(sheet, sep='\t', dtype=str, keep_default_na=False)
    required = {'sample_id','sample_type','condition','replicate','raw_path','run','library_reference'}
    if not required.issubset(samples):
        raise ValueError('Sample sheet needs: ' + ', '.join(sorted(required)))
    if samples.duplicated(['sample_type','sample_id']).any() or samples.run.duplicated().any() or samples.raw_path.duplicated().any():
        raise ValueError('Duplicate sample IDs, acquisitions or runs')
    if any(not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*', value) for value in samples.sample_type):
        raise ValueError('Invalid cohort label')
    for key in ('donor', 'biological_sample'):
        if key in samples and (samples[key].eq('').any() or samples.duplicated(['sample_type',key]).any()):
            raise ValueError('Missing/repeated biological identity requires a tailored design: ' + key)
    if args.independent_samples:
        if 'source_replicate' not in samples:
            samples['source_replicate'] = samples.replicate
        samples['replicate'] = samples.sample_id
    options = differential_options(dict(enabled=True, comparisons='all_pairwise', fdr=args.fdr,
                                        min_replicates=args.min_replicates, rscript=args.rscript, gene_level=True))
    services = {}
    for value in args.service_matrix:
        group, path = value.split('=', 1)
        if group in services or group not in set(samples.sample_type):
            raise ValueError('Duplicate or unknown service cohort: ' + group)
        services[group] = str(Path(path).resolve())
    method = Path(args.method).resolve() if args.method else ROOT/'modules/diapasef/method.toml'
    fasta = source/'prepared/search.fasta'
    reports, audits = {}, {}
    for group, selected in samples.groupby('sample_type', sort=True):
        validate_design_rows(selected.to_dict('records'), options)
        if group in services:
            load_matrix(services[group], selected)
        candidates = [source/base/group/'cohort/report.parquet' for base in ('analysis_details/search', 'search')]
        candidates = [p for p in candidates if p.is_file()]
        if len(candidates) != 1:
            raise ValueError('Expected exactly one completed report for ' + group)
        reports[group] = candidates[0]
        audits[group] = audit_report(candidates[0], selected, method, fasta)
        print('SEARCH VERIFIED: '+group+'; reuse '+str(candidates[0]), flush=True)
    output.mkdir(parents=True, exist_ok=True)
    audit = output/'audit'; audit.mkdir()
    logs = output/'logs'; logs.mkdir()
    inputs = output/'inputs'; inputs.mkdir()
    state = dict(status='checking_raw_integrity', started=stamp(), source_run=str(source), repository=str(ROOT),
                 command=[sys.executable, *sys.argv], independent_samples_confirmed=args.independent_samples,
                 differential=options, searches=audits,
                 raw_hash_note='New integrity baseline; no historical raw checksums were supplied. Search report copies are compared against completed task outputs.')
    state_path = audit/'run.json'
    state_path.write_text(json.dumps(state, indent=2)+'\n')
    try:
        raw = audit_raw(samples, args.hash_raw)
        (audit/'raw_integrity.json').write_text(json.dumps(raw, indent=2)+'\n')
        (audit/'search_reuse.json').write_text(json.dumps(audits, indent=2)+'\n')
        samples.to_csv(inputs/'samples.tsv', sep='\t', index=False)
        export_matrix_formats(inputs/'samples.tsv')
        shutil.copy2(method, inputs/'method.toml')
        shutil.copy2(fasta, inputs/'search.fasta')
        state['input_provenance'] = [file_provenance(sheet), file_provenance(method), file_provenance(fasta)] + [file_provenance(p) for p in services.values()]
        state['code_provenance'] = [file_provenance(p) for folder in ('pipeline_lib', 'scripts', 'modules/diapasef')
                                    for p in sorted((ROOT/folder).glob('*')) if p.suffix in ('.py','.R','.nf','.toml')]
        state.update(status='running_reporting', raw_files_checked=sum(len(r['files']) for r in raw))
        state_path.write_text(json.dumps(state, indent=2)+'\n')
        env = os.environ.copy()
        env['R_LIBS'] = os.pathsep.join(filter(None, [str(ROOT/'.runtime/R-library'), env.get('R_LIBS')]))
        env['MPLCONFIGDIR'] = str(output/'analysis_details/.matplotlib')
        groups = []
        for group, selected in samples.groupby('sample_type', sort=True):
            group_sheet = inputs/(group+'.tsv'); selected.to_csv(group_sheet, sep='\t', index=False)
            destination = output/'analysis_details/groups'/group
            command = [sys.executable, str(ROOT/'modules/diapasef/05_count_proteins.py'), '--samples', str(group_sheet),
                       '--report', str(reports[group]), '--method', str(inputs/'method.toml'), '--fasta', str(inputs/'search.fasta'),
                       '--outdir', str(destination), '--differential', json.dumps(options)]
            if group in services:
                command += ['--service-matrix', services[group]]
            with (logs/(group+'_reporting.log')).open('w') as handle:
                subprocess.run(command, cwd=ROOT, env=env, stdout=handle, stderr=subprocess.STDOUT, check=True)
            groups.append(str(destination))
            print('REPORTING COMPLETE: '+group, flush=True)
        with (logs/'collect_results.log').open('w') as handle:
            subprocess.run([sys.executable, str(ROOT/'modules/diapasef/06_collect_results.py'), '--samples', str(inputs/'samples.tsv'),
                            '--groups', *groups], cwd=output, env=env, stdout=handle, stderr=subprocess.STDOUT, check=True)
        final = output/'final_results'
        (final/'FINAL_RESULTS_COMPLETE').unlink(missing_ok=True)
        # Preserve completed optional evidence with provenance; do not imply a new phosphosite model.
        historical = source/'phosphorylation'
        optional_records = []
        if (historical/'final/REPORT.md').is_file():
            target = final/'optional_analyses/phosphorylation'
            for name in ('final', 'groups'):
                if (historical/name).is_dir():
                    shutil.copytree(historical/name, target/name)
                    for path in sorted((historical/name).rglob('*')):
                        if path.is_file():
                            optional_records.append(same_file_content(path, target/name/path.relative_to(historical/name)))
            for path in sorted(target.rglob('*.tsv')):
                export_matrix_formats(path)
            (target/'README.md').write_text('Preserved phosphorylation evidence from the completed source run. Checksums are recorded in audit/optional_evidence.json at the run root. No phosphosite differential model was fitted; protein/gene all-pairs results are in differential_expression/. See final/REPORT.md for original thresholds and limitations.\n')
            (audit/'optional_evidence.json').write_text(json.dumps(optional_records, indent=2)+'\n')
        qc = final/'quality_control/reuse'; qc.mkdir(parents=True)
        shutil.copy2(audit/'search_reuse.json', qc/'search_reuse.json')
        (final/'README.md').write_text((final/'README.md').read_text()+
            '\nThis run reused completed DIA-NN searches after SHA256 comparison with the original task outputs, full Parquet decoding and sample/method/FASTA checks. Raw integrity and source/code checksums are in ../audit/. Rows are analyzed as independent biological samples; when explicitly confirmed with --independent-samples, the original replicate label is retained as source_replicate. Protein and gene limma-trend analyses include every unordered pair of conditions. Optional historical phosphorylation evidence is preserved separately and has no newly fitted site-level comparisons.\n')
        write_dashboard(final, 'DIA-PASEF — all pairwise comparisons')
        for group, report in reports.items():
            if file_provenance(report)['sha256'] != audits[group]['report']['sha256']:
                raise ValueError('Source report changed during analysis: ' + str(report))
        (final/'FINAL_RESULTS_COMPLETE').write_text('Verified DIA reuse, quantification and all-pairs protein/gene models completed.\n')
        write_output_manifest(final)
        state.update(status='completed', ended=stamp(), comparisons=len(pd.read_csv(final/'results_index.csv')),
                     final_results=str(final), raw_searches_reused=len(reports), raw_searches_rerun=0)
        state_path.write_text(json.dumps(state, indent=2)+'\n')
        (output/'START_HERE.html').write_text('<!doctype html><meta charset="utf-8"><meta http-equiv="refresh" content="0;url=final_results/START_HERE.html"><a href="final_results/START_HERE.html">Open organized results</a>')
        print('COMPLETED: '+str(final/'START_HERE.html'), flush=True)
    except Exception as error:
        state.update(status='failed', ended=stamp(), error=str(error))
        state_path.write_text(json.dumps(state, indent=2)+'\n')
        raise


if __name__ == '__main__':
    main()
