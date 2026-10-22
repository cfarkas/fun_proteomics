#!/usr/bin/env python3
"""Portable DDA-PASEF input validation, tool preparation and FragPipe execution."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import sqlite3
import stat
import subprocess
import sys
import urllib.request
import zipfile


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def save(path, data):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(data, indent=2, ensure_ascii=False) + '\n')


def read(path):
    return json.loads(Path(path).read_text())


def resolve(value, base):
    path = Path(value).expanduser()
    result = str((path if path.is_absolute() else Path(base) / path).resolve())
    if any(c in result for c in '\n\r\t'):
        raise ValueError('Paths cannot contain tabs or line breaks')
    return result


def validate(args):
    import jsonschema
    cfg = read(args.config)
    from optional_analyses import migrate_options
    migrate_options(cfg)
    jsonschema.Draft202012Validator(read(args.schema)).validate(cfg)
    from optional_analyses import validate_options
    validate_options(cfg)
    if 'proteomics' in cfg:
        choice=cfg['proteomics']
        if choice.get('confirmed') is not True or choice.get('acquisition')!='DDA-PASEF' or choice.get('analysis') not in ('global-label-free','phosphoproteomics','immunopeptidomics'):
            raise ValueError('A confirmed, supported DDA-PASEF workflow is required')
    base = args.base
    cfg.setdefault('resources', {'cpus': 32, 'memory_gib': 128})
    cfg.setdefault('organism', None)
    cfg['outdir'] = resolve(cfg['outdir'], base)
    ids, raw_paths, audit = set(), set(), []
    for sample in cfg['samples']:
        sid = sample['sample_id']
        if sid in ids:
            raise ValueError(f'Duplicate sample_id: {sid}')
        ids.add(sid)
        sample['raw_path'] = resolve(sample['raw_path'], base)
        raw = Path(sample['raw_path'])
        if raw in raw_paths:
            raise ValueError(f'Raw directory used twice: {raw}')
        raw_paths.add(raw)
        if not raw.is_dir() or raw.suffix.lower() != '.d':
            raise ValueError(f'Expected native Bruker .d directory: {raw}')
        for name in ['analysis.tdf', 'analysis.tdf_bin']:
            if not (raw / name).is_file():
                raise ValueError(f'Missing {raw / name}')
        # URI encoding also handles spaces, # and ? in original acquisition names.
        with sqlite3.connect((raw / 'analysis.tdf').as_uri() + '?mode=ro', uri=True) as db:
            frame_types = dict(db.execute('SELECT MsMsType, COUNT(*) FROM Frames GROUP BY MsMsType'))
            if set(frame_types) != {0, 8}:
                raise ValueError(f'{sid}: expected DDA-PASEF frame types 0 and 8; found {frame_types}')
            precursors = db.execute('SELECT COUNT(*) FROM Precursors').fetchone()[0]
        audit.append({'sample_id': sid, 'frame_counts': frame_types, 'precursors': precursors,
                      'analysis_tdf_sha256': digest(raw / 'analysis.tdf'),
                      'raw_files': {str(p.relative_to(raw)): {'size': p.stat().st_size,
                          'mtime_ns': p.stat().st_mtime_ns} for p in sorted(raw.rglob('*')) if p.is_file()}})
    for key, value in cfg['tools'].items():
        cfg['tools'][key] = resolve(value, base)
        if not Path(cfg['tools'][key]).is_file():
            raise ValueError(f'Missing tool archive: {cfg["tools"][key]}')
    database = cfg['database']
    for key in ['fasta', 'contaminants_zip']:
        if key in database:
            database[key] = resolve(database[key], base)
            if not Path(database[key]).is_file():
                raise ValueError(f'Missing database input: {database[key]}')
    if digest(database['fasta']) != database['sha256']:
        raise ValueError('FASTA SHA256 does not match the configuration')
    if not database['contains_decoys'] and 'contaminants_zip' not in database:
        raise ValueError('A target-only FASTA requires contaminants_zip and contaminants_sha256')
    if 'contaminants_zip' in database and digest(database['contaminants_zip']) != database.get('contaminants_sha256'):
        raise ValueError('Contaminant archive SHA256 does not match the configuration')
    for cohort, spec in cfg.get('service_quantification',{}).items():
        if cohort not in {s['cohort'] for s in cfg['samples']}: raise ValueError('Unknown benchmark cohort: '+cohort)
        spec['matrix']=resolve(spec['matrix'],base)
        from quantification import load_matrix
        import pandas as pd
        samples=spec.get('samples') or selected_samples(cfg,cohort)
        load_matrix(spec['matrix'],pd.DataFrame(samples),sheet=spec.get('sheet'))
        spec['sha256']=digest(spec['matrix'])
    cfg['raw_audit'] = audit
    save('validated.json', cfg)
    save('input_audit.json', {'samples': audit, 'config_sha256': digest(args.config)})
    print(f'Validated {len(ids)} native DDA-PASEF samples in {len({s["cohort"] for s in cfg["samples"]})} cohorts')


def safe_extract(archive, target):
    target = Path(target).resolve()
    target.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive) as z:
        for info in z.infolist():
            dest = (target / info.filename).resolve()
            if not dest.is_relative_to(target) or stat.S_ISLNK(info.external_attr >> 16):
                raise ValueError(f'Unsafe archive member: {info.filename}')
        z.extractall(target)
        for info in z.infolist():
            mode = (info.external_attr >> 16) & 0o777
            if mode:
                (target / info.filename).chmod(mode)


def prepare_tools(args):
    cfg, versions = read(args.plan), read(args.versions)
    root = Path('toolset').resolve()
    root.mkdir(exist_ok=True)
    provenance = {}
    for name, key in [('FragPipe', 'fragpipe_zip'), ('MSFragger', 'msfragger_zip'),
                      ('IonQuant', 'ionquant_zip'), ('diaTracer', 'diatracer_zip')]:
        spec = versions[name]
        archive = cfg['tools'].get(key)
        if archive is None:
            if name != 'FragPipe':
                raise ValueError(f'Register and provide {key}; see README mandatory third parties')
            archive = Path('FragPipe-24.0-linux.zip')
            urllib.request.urlretrieve(spec['url'], archive)
        observed = digest(archive)
        if observed != spec['archive_sha256']:
            raise ValueError(f'{name} archive differs from pinned {spec["version"]}; expected {spec["archive_sha256"]}, got {observed}')
        safe_extract(archive, root)
        provenance[name] = {'version': spec['version'], 'archive_sha256': observed}
    required = ['fragpipe-24.0/bin/fragpipe', 'MSFragger-4.4.1/MSFragger-4.4.1.jar',
                'IonQuant-1.11.20/IonQuant-1.11.20.jar', 'diatracer-2.2.1/diaTracer-2.2.1.jar',
                'fragpipe-24.0/tools/diann/1.8.2_beta_8/linux/diann-1.8.1.8']
    for relative in required:
        path = root / relative
        if not path.is_file():
            raise ValueError(f'Unexpected archive layout: missing {relative}')
        if relative.endswith(('bin/fragpipe', 'diann-1.8.1.8')):
            path.chmod(path.stat().st_mode | 0o111)
    executable = root / required[-1]
    link = subprocess.run(['ldd', str(executable)], text=True, capture_output=True)
    if 'not found' in link.stdout:
        raise RuntimeError('DIA-NN missing system libraries:\n' + link.stdout)
    provenance['bundled'] = {'MSBooster': '1.4.14', 'Percolator': '3.7.1',
                              'Philosopher': '5.1.3-RC9', 'DIA-NN': '1.8.2_beta_8'}
    provenance['runtime'] = {'python': sys.version, 'platform': platform.platform(),
        'java': subprocess.run(['java', '-version'], capture_output=True, text=True, check=True).stderr}
    save(root / 'software_versions.json', provenance)


def fasta_records(text):
    header, seq = None, []
    for line in text.splitlines():
        if line.startswith('>'):
            if header is not None:
                yield header, ''.join(seq)
            header, seq = line[1:], []
        elif line.strip():
            if header is None:
                raise ValueError('Sequence before first FASTA header')
            seq.append(line.strip())
    if header is not None:
        yield header, ''.join(seq)


def accession(header):
    fields = header.split()[0].split('|')
    return fields[1] if len(fields) > 1 else fields[0]


def build_database(args):
    db = read(args.plan)['database']
    proteins = list(fasta_records(Path(db['fasta']).read_text()))
    if not proteins or any(not s or not re.fullmatch('[A-Z]+', s) for _, s in proteins):
        raise ValueError('FASTA must contain nonempty uppercase amino-acid sequences')
    names = [h.split()[0] for h, _ in proteins]
    if len(names) != len(set(names)):
        raise ValueError('Duplicate FASTA identifiers')
    if db['contains_decoys']:
        targets = {h.split()[0] for h, _ in proteins if not h.startswith('rev_')}
        decoys = {h.split()[0][4:] for h, _ in proteins if h.startswith('rev_')}
        if not targets or targets != decoys:
            raise ValueError('Prebuilt FASTA must have exactly one rev_ decoy per target identifier')
        if not any(h.startswith('contam_') for h in targets):
            raise ValueError('Prebuilt FASTA must include contam_ entries')
        shutil.copyfile(db['fasta'], 'search.fasta')
        method = 'Prebuilt target+contaminant+decoy FASTA; copied without modification'
    else:
        if any(h.startswith(('rev_', 'contam_')) for h, _ in proteins):
            raise ValueError('Target-only FASTA already contains rev_ or contam_ entries')
        accessions = {accession(h) for h, _ in proteins}
        with zipfile.ZipFile(db['contaminants_zip']) as z:
            contaminants = list(fasta_records(z.read('crap-gpmdb.fas').decode()))
        proteins += [('contam_' + h, s) for h, s in contaminants if accession(h) not in accessions]
        decoys = []
        for header, sequence in proteins:
            segments = re.findall(r'.*?[KR]|.+$', sequence)
            sequence = ''.join(s[0] + s[1:-1][::-1] + s[-1] if len(s) > 1 else s for s in segments)
            decoys.append(('rev_' + header, sequence))
        proteins += decoys
        Path('search.fasta').write_text(''.join(f'>{h}\n{s}\n' for h, s in proteins))
        method = 'KR segment internal reversal; preserve first and last residue; prefix rev_'
    save('database_provenance.json', {**db, 'search_fasta_sha256': digest('search.fasta'),
         'entries': len(proteins), 'decoy_method': method})


def selected_samples(plan, cohort):
    samples = [s for s in plan['samples'] if s['cohort'] == cohort]
    if not samples:
        raise ValueError(f'Unknown or empty cohort: {cohort}')
    return samples


def search_properties(args, cfg):
    choice = cfg.get('proteomics', {})
    if choice.get('analysis') == 'phosphoproteomics' and digest(args.template) != choice.get('module_sha256'):
        raise ValueError('Phosphoproteomics module changed; confirm the installed module again.')
    hla = choice.get('analysis') == 'immunopeptidomics'
    if hla and (args.mode != 'DDA' or cfg.get('modes') != ['DDA']):
        raise ValueError('The immunopeptidomics preset supports DDA searches only.')
    properties = dict(line.split('=', 1) for line in Path(args.template).read_text().splitlines()
                      if line and not line.startswith('#') and '=' in line)
    if hla:
        from immunopeptidomics import validate_workflow
        validate_workflow(properties)
    properties.update({'database.db-path': str(Path(args.database).resolve()), 'database.decoy-tag': 'rev_',
        'workflow.description': f'{args.cohort} {args.mode}; ' + ('HLA-enriched peptides, MBR off, 1% PSM FDR' if hla else 'cohort-specific inference and MBR'),
        'workflow.saved-with-ver': '24.0', 'workflow.threads': str(args.cpus), 'workflow.ram': str(args.memory),
        'msfragger.output_report_topN_dda_plus': '1' if hla else '5', 'percolator.keep-tsv-files': 'true',
        'phi-report.filter': '--psm 0.01 --pep 1 --ion 1 --prot=1' if hla else '--sequential --prot 0.01 --picked',
        'tab-run.delete_temp_files': 'false'})
    return properties


def search(args):
    cfg = read(args.plan)
    properties = search_properties(args, cfg)
    samples = selected_samples(cfg, args.cohort)
    root = Path('run').resolve()
    root.mkdir(exist_ok=True)
    (root/'SEARCH_COMPLETE').unlink(missing_ok=True)
    raw_root = Path('raw').resolve()
    raw_root.mkdir(exist_ok=True)
    manifest = []
    for sample in samples:
        source = Path(sample['raw_path'])
        target = raw_root / (sample['sample_id'] + source.suffix)
        experiment = sample.get('biological_sample') or sample['sample_id']
        if args.dry_run:
            manifest.append('\t'.join([str(source), experiment, '', args.mode]))
            continue
        marker = raw_root / (sample['sample_id'] + '.copy_complete.json')
        if not marker.exists():
            if target.exists():
                raise ValueError(f'Incomplete raw copy requires inspection before retry: {target}')
            subprocess.run(['cp', '--reflink=auto', '-a', '-L', str(source), str(target)], check=True)
            original = {str(p.relative_to(source)): p.stat().st_size for p in source.rglob('*') if p.is_file()} if source.is_dir() else {'sha256': digest(source)}
            copied = {str(p.relative_to(target)): p.stat().st_size for p in target.rglob('*') if p.is_file()} if target.is_dir() else {'sha256': digest(target)}
            if original != copied:
                raise ValueError(f'Raw staging verification failed for {sample["sample_id"]}')
            save(marker, {'source': str(source), 'files': original, 'copy_complete': True})
        # Vendor caches in a previously completed staged copy may change on retry.
        # They remain confined to this task's independent acquisition copy.
        # A unique experiment name and blank replicate retain the exact sample ID.
        manifest.append('\t'.join([str(target), experiment, '', args.mode]))
    (root / 'fragpipe.manifest').write_text('\n'.join(manifest) + '\n')
    (root / 'fragpipe.workflow').write_text('\n'.join(f'{k}={v}' for k, v in properties.items()) + '\n')
    tools = Path(args.tools).resolve()  # FragPipe discovery does not descend symlink directories.
    fragpipe = tools / 'fragpipe-24.0/bin/fragpipe'
    command = [str(fragpipe), '--headless', '--workflow', str(root / 'fragpipe.workflow'),
               '--manifest', str(root / 'fragpipe.manifest'), '--workdir', str(root),
               '--threads', str(args.cpus), '--ram', str(args.memory),
               '--config-tools-folder', str(tools), '--config-diann',
               str(tools / 'fragpipe-24.0/tools/diann/1.8.2_beta_8/linux/diann-1.8.1.8'),
               '--config-python', sys.executable]
    if args.dry_run:
        command.append('--dry-run')
    if hasattr(os, 'sched_getaffinity'):
        allowed = sorted(os.sched_getaffinity(0))
        if len(allowed) < args.cpus:
            raise ValueError(f'Only {len(allowed)} CPUs available; requested {args.cpus}')
        os.sched_setaffinity(0, allowed[:args.cpus])
    env = dict(os.environ, OMP_NUM_THREADS=str(args.cpus), OPENBLAS_NUM_THREADS=str(args.cpus), MKL_NUM_THREADS=str(args.cpus))
    save(root / 'run_provenance.json', {'cohort': args.cohort, 'mode': args.mode, 'samples': samples,
        'database_sha256': digest(args.database), 'template_sha256': digest(args.template),
        'command': command, 'cpus': args.cpus, 'memory_gib': args.memory,
        'dry_run': args.dry_run, 'stub': False, 'analysis': cfg.get('proteomics',{}).get('analysis','global-label-free'),
        'immunopeptidomics': cfg.get('immunopeptidomics')})
    shutil.copyfile(tools / 'software_versions.json', root / 'software_versions.json')
    with (root / 'pipeline.log').open('w') as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, env=env, check=True)
    if not args.dry_run and not (root / 'combined_protein.tsv').is_file():
        raise RuntimeError('FragPipe exited without combined_protein.tsv; inspect pipeline.log')
    (root / ('DRY_RUN_COMPLETE' if args.dry_run else 'SEARCH_COMPLETE')).write_text('complete\n')
    if not args.dry_run:
        from search_checkpoint import write_checkpoint
        write_checkpoint(root, cfg, args.cohort, args.mode, properties, digest(args.database),
                         read(tools/'software_versions.json'))


def stub_search(args):
    samples = selected_samples(read(args.plan), args.cohort)
    root = Path('run')
    root.mkdir(exist_ok=True)
    rows = []
    for i in range(3):
        row = {'Protein': f'sp|TEST{i}|SYNTHETIC', 'Protein ID': f'TEST{i}', 'Gene': f'TEST{i}',
               'Description': 'SYNTHETIC STUB - NOT SCIENTIFIC RESULTS', 'Organism': 'Homo sapiens',
               'Indistinguishable Proteins': ''}
        for j, sample in enumerate(samples):
            sid = sample['sample_id']
            row[sid + ' Spectral Count'] = 2 if i == 0 or (args.mode == 'DDA+' and i == 1) else 0
            row[sid + ' Intensity'] = 100 if i < 2 else 0
            if len(samples) > 1:
                row[sid + ' MaxLFQ Intensity'] = 100 if i < 2 else 0
        rows.append(row)
    with (root / 'combined_protein.tsv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    save(root / 'run_provenance.json', {'cohort': args.cohort, 'mode': args.mode, 'samples': samples, 'stub': True})
    (root / 'STUB_NOT_SCIENTIFIC_RESULTS').write_text('Synthetic workflow test only.\n')
    (root / 'SEARCH_COMPLETE').write_text('stub\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('validate')
    for name in ['config', 'schema', 'base']:
        p.add_argument('--' + name, required=True)
    p.set_defaults(func=validate)
    p = sub.add_parser('tools')
    p.add_argument('--plan', required=True)
    p.add_argument('--versions', required=True)
    p.set_defaults(func=prepare_tools)
    p = sub.add_parser('database')
    p.add_argument('--plan', required=True)
    p.set_defaults(func=build_database)
    for name, func in [('search', search), ('stub-search', stub_search)]:
        p = sub.add_parser(name)
        p.add_argument('--plan', required=True)
        p.add_argument('--cohort', required=True)
        p.add_argument('--mode', choices=['DDA', 'DDA+'], required=True)
        p.add_argument('--tools')
        p.add_argument('--template')
        p.add_argument('--database')
        p.add_argument('--cpus', type=int, default=32)
        p.add_argument('--memory', type=int, default=128)
        p.add_argument('--dry-run', action='store_true')
        p.set_defaults(func=func)
    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
