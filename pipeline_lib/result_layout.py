"""Deliver browsable scientific results while keeping model records in analysis_details."""
import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from matrix_catalog import build_matrix_catalog

FEATURES = {'protein': 'Proteins', 'gene': 'Genes (proteomics)', 'peptide': 'Peptides',
            'modification_gene': 'Modified peptides per gene', 'modification_protein': 'Modified peptides per protein',
            'modification': 'Supporting modified-peptide PSMs', 'isoform': 'Isoform-discriminating evidence'}
INDEX_FIELDS = ['section', 'analysis', 'cohort', 'mode', 'feature_type', 'method', 'test_type', 'comparison',
                'reported_features', 'tested_features', 'significant_features', 'status', 'fdr',
                'full_results', 'significant_results']
GENE_FEATURES = {'gene', 'modification_gene'}


def ordered_results(rows):
    return sorted(rows, key=lambda row: (row['section'] != 'Main results',
        row['feature_type'] != 'protein', row['analysis'], row['cohort'], row['mode'],
        row['comparison'], row['feature_type'], row['method']))


def read_rows(path):
    if not Path(path).is_file():
        return []
    with Path(path).open(newline='') as handle:
        return list(csv.DictReader(handle))


def write_rows(path, fields, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024*1024), b''):
            value.update(chunk)
    return value.hexdigest()


def segment(value):
    value = re.sub(r'[^A-Za-z0-9_.-]', '_', str(value)).strip('._')
    if not value:
        raise ValueError('Empty result path component')
    return value


def safe_source(base, relative):
    path = (Path(base)/relative).resolve()
    if not path.is_relative_to(Path(base).resolve()):
        raise ValueError('Result index points outside its results directory: '+str(relative))
    return path


def publish_results(out, preserve_original=False, search_summary=None):
    """Rearrange existing outputs only: never search spectra, refit or modify result tables."""
    out = Path(out).resolve()
    if not (out/'FINAL_RESULTS_COMPLETE').is_file():
        raise ValueError('Result reorganization requires a completed final_results directory.')
    if (out/'main_results').is_dir() and (out/'results_index.csv').is_file():
        result = upgrade_gene_layout(out, preserve_original)
        build_matrix_catalog(out)
        return result
    job = json.loads((out/'analysis_settings.json').read_text())
    audit_target = out.parent/'analysis_details/reporting/result_audit'
    work = Path(tempfile.mkdtemp(prefix='.result_layout_', dir=out.parent))
    final, audit = work/'final_results', work/'audit'
    final.mkdir(); audit.mkdir()
    index, files, manifest = [], [], []
    copied, destinations = set(), {}

    def copy(source, relative, technical=False, title='', category=''):
        source = Path(source).resolve()
        if not source.is_file():
            raise ValueError('Missing result file: '+str(source))
        destination = (audit if technical else final)/relative
        key = str(destination)
        if key in destinations:
            if destinations[key] != str(source):
                raise ValueError('Ambiguous result destination: '+str(relative))
            return str(relative)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
        checksum = digest(source)
        if digest(destination) != checksum:
            raise ValueError('Result copy verification failed: '+str(relative))
        destinations[key] = str(source); copied.add(source)
        manifest.append({'area': 'audit' if technical else 'results', 'file': str(relative),
                         'source': os.path.relpath(source, out), 'sha256': checksum, 'bytes': source.stat().st_size})
        if not technical and title:
            files.append({'category': category, 'title': title, 'path': str(relative)})
        return str(relative)

    def copy_tree(source, relative, technical=False, title_prefix='', category=''):
        if not Path(source).is_dir():
            return
        for path in sorted(Path(source).rglob('*')):
            if path.is_file() and not any(part.startswith('mass_shift') or part=='open_search' for part in path.relative_to(source).parts):
                copy(path, Path(relative)/path.relative_to(source), technical,
                     title=(title_prefix+' · '+path.name if title_prefix else ''), category=category)

    def status(tested, significant, reported):
        if not reported:
            return 'No eligible evidence'
        if not tested:
            return 'Insufficient evidence'
        return 'Significant results' if significant else 'No significant results'

    def differential(source, section, analysis, relative):
        for row in read_rows(source/'differential/results_index.csv'):
            if row['feature_type'] not in FEATURES:
                continue
            test_type = row['test_type']
            comparison = row['comparison']
            if comparison == Path(row['all_results']).name:
                comparison = Path(comparison).stem.removeprefix(row['method']+'_')
            contrast = 'overall_condition_test' if test_type=='omnibus' else segment(comparison)
            cutoff = job.get('differential',{}).get('fdr',.05) if section=='Main results' else job.get('optional_analyses',{}).get('statistics',{}).get('fdr',.05)
            gene = section=='Main results' and row['feature_type'] in GENE_FEATURES
            if gene and job.get('optional_analyses',{}).get('gene_counts') is False:
                continue
            folder = Path('optional_analyses/gene_counts/comparisons' if gene else relative)/segment(row['cohort'])/('ddaplus' if row['mode']=='DDA+' else 'dda')/contrast/segment(row['feature_type'])
            full = copy(safe_source(source/'differential', row['all_results']), folder/(segment(row['method'])+'_all_results.csv'))
            sig = copy(safe_source(source/'differential', row['significant_results']), folder/(segment(row['method'])+'_significant_results.csv'))
            tested, reported, significant = (int(row[key]) for key in ('tested_features','reported_features','significant_features'))
            index.append(dict(section='Optional analyses' if gene else section,analysis='Gene counts and comparisons' if gene else analysis,cohort=row['cohort'],mode=row['mode'],
                feature_type=row['feature_type'],method=row['method'],test_type=test_type,
                comparison='Any condition difference' if test_type=='omnibus' else comparison,
                reported_features=reported,tested_features=tested,significant_features=significant,
                status=status(tested,significant,reported),fdr=cutoff,full_results=full,significant_results=sig))

    try:
        differential(out, 'Main results', 'Differential expression', 'main_results/differential_expression')
        optional = out/'optional_analyses'
        statistics = optional/'statistics'
        differential(statistics, 'Optional analyses', 'Peptide and isoform comparisons', 'optional_analyses/statistics/comparisons')
        # Counts are grouped by cohort and search, avoiding blank columns from other cohorts.
        for group in job.get('groups', []):
            cohort, mode = group['cohort'], group['mode']
            slug = 'ddaplus' if mode=='DDA+' else 'dda'
            raw = (out/group['directory']).resolve()
            prefix = Path('main_results/counts')/segment(cohort)/slug
            for kind, label in FEATURES.items():
                gene = kind in GENE_FEATURES
                if gene and job.get('optional_analyses',{}).get('gene_counts') is False:
                    continue
                path = raw/(kind+'_counts.csv')
                if path.is_file():
                    target = Path('optional_analyses/gene_counts/counts')/segment(cohort)/slug if gene else prefix
                    copy(path, target/'raw'/path.name, title=cohort+' · '+mode+' · '+label+' · raw counts', category='Gene count matrices' if gene else 'Main count matrices')
            for filename in ('protein_intensity.tsv','protein_maxlfq_intensity.tsv'):
                path=raw/filename
                if path.is_file():
                    copy(path,Path('matrices')/segment(cohort)/slug/filename,
                         title=cohort+' · '+mode+' · '+filename.removesuffix('.tsv').replace('_',' '),
                         category='Protein abundance matrices')
            for name, source in [('normalized',out/'matrices'/cohort/slug), ('zscore',out.parent/'zscore'/cohort/slug)]:
                for path in sorted(source.glob('*.csv')):
                    gene = path.name.startswith(('gene_', 'modification_gene_'))
                    if gene and job.get('optional_analyses',{}).get('gene_counts') is False:
                        continue
                    target = Path('optional_analyses/gene_counts/counts')/segment(cohort)/slug if gene else prefix
                    copy(path,target/name/path.name,title=cohort+' · '+mode+' · '+name+' · '+path.name,
                         category='Gene count matrices' if gene else 'Main count matrices')
        for name in ['sample_annotations.csv','analysis_summary.csv','comparisons.csv','method_notes.csv']:
            if (out/name).is_file():
                copy(out/name,Path('main_results/quality_control')/name,title=name.replace('_',' ').removesuffix('.csv'),category='Study and coverage')
        for name in ['sample_detection_counts.tsv','cohort_detection_counts.tsv','dda_vs_ddaplus.tsv']:
            path = (Path(search_summary) if search_summary else out.parent/'analysis_details/final')/name
            if path.is_file():
                copy(path,Path('main_results/quality_control')/name,title=name.replace('_',' ').removesuffix('.tsv'),category='Study and coverage')
        condition_path = out/'differential/condition_columns.csv'
        if condition_path.is_file():
            copy(condition_path,'main_results/quality_control/condition_columns.csv',title='Condition names and table column labels',category='Study and coverage')
        for branch, label in [('peptide_isoforms','Peptide and isoform evidence')]:
            base = optional/branch
            for name in ['peptide_evidence.csv','peptide_condition_summary.csv','isoform_evidence.csv']:
                for path in sorted(base.glob('*/*/'+name)):
                    cohort, slug = path.relative_to(base).parts[:2]
                    copy(path,Path('optional_analyses')/branch/path.relative_to(base),
                         title=cohort+' · '+('DDA+' if slug=='ddaplus' else 'DDA')+' · '+name.replace('_',' ').removesuffix('.csv'),category=label)
        for path in sorted((out/'immunopeptidomics').glob('*/*/*')):
            if path.is_file() and path.name in ('peptides.csv','sample_summary.csv','length_distribution.csv','binding_input.txt'):
                copy(path,Path('optional_analyses/immunopeptidomics')/path.relative_to(out/'immunopeptidomics'),
                     title=str(path.relative_to(out/'immunopeptidomics')),category='HLA peptide evidence')
        # Optional count models use their own normalization and evidence definitions.
        if (statistics/'analysis_settings.json').is_file():
            stats_job = json.loads((statistics/'analysis_settings.json').read_text())
            for group in stats_job.get('groups', []):
                cohort, mode = group['cohort'], group['mode'];slug='ddaplus' if mode=='DDA+' else 'dda'
                raw = (statistics/group['directory']).resolve()
                prefix = Path('optional_analyses/statistics/counts')/segment(cohort)/slug
                for kind in ('peptide','isoform'):
                    if (raw/(kind+'_counts.csv')).is_file():
                        copy(raw/(kind+'_counts.csv'),prefix/'raw'/(kind+'_counts.csv'),
                             title=cohort+' · '+mode+' · '+FEATURES[kind]+' · optional counts',category='Optional count matrices')
                for name in ('matrices','zscore'):
                    copy_tree(statistics/name/cohort/slug,prefix/name,title_prefix=cohort+' · '+mode+' · optional '+name,category='Optional count matrices')
        drugs = optional/'drug_interaction'
        drug_settings = json.loads((drugs/'settings.json').read_text()) if (drugs/'settings.json').is_file() else {}
        for row in read_rows(drugs/'analysis_summary.csv'):
            kind, method = row['feature_type'], row.get('method') or 'edgeR'
            slug = 'ddaplus' if row['mode']=='DDA+' else 'dda'
            default = Path(row['cohort'])/slug/kind/'interaction.csv'
            source = safe_source(drugs, row.get('results_file') or default)
            folder = Path('optional_analyses/drug_interaction')/segment(row['cohort'])/slug/segment(kind)
            full = copy(source,folder/(segment(method)+'_all_results.csv'))
            sig = copy(source.with_name('significant.csv'),folder/(segment(method)+'_significant_results.csv'))
            tested, reported, significant = (int(row[key]) for key in ('tested_features','input_features','significant_features'))
            index.append(dict(section='Optional analyses',analysis='Drug interaction',cohort=row['cohort'],mode=row['mode'],
                feature_type=kind,method=method,test_type='interaction',comparison='Combination − drug A − drug B + control',
                reported_features=reported,tested_features=tested,significant_features=significant,
                status=status(tested,significant,reported),fdr=drug_settings.get('fdr',.05),full_results=full,significant_results=sig))
        copy_tree(out/'benchmark',Path('benchmark'),title_prefix='Service benchmark',category='Service benchmarks')
        # The selected final files are copied unchanged. Remaining supported outputs are model records.
        for path in sorted(out.rglob('*')):
            if not path.is_file() or path.resolve() in copied:
                continue
            relative = path.relative_to(out)
            if 'open_search' in relative.parts or any(part.startswith('mass_shift') for part in relative.parts):
                continue
            if path.name in ('START_HERE.html','FINAL_RESULTS_COMPLETE'):
                continue
            if relative.parts[:2]==('optional_analyses','statistics'):
                target = Path('optional_statistics')/Path(*relative.parts[2:])
            elif relative.parts[:2]==('optional_analyses','drug_interaction'):
                target = Path('drug_interaction')/Path(*relative.parts[2:])
            elif relative.parts[:2]==('optional_analyses','peptide_isoforms'):
                target = Path('peptide_isoforms')/Path(*relative.parts[2:])
            else:
                target = Path('main_statistics')/relative
            copy(path,target,technical=True)
        index = ordered_results(index)
        write_rows(final/'results_index.csv',INDEX_FIELDS,index)
        write_rows(audit/'file_manifest.csv',['area','file','source','sha256','bytes'],manifest)
        # Human-readable summaries expose no pooled significance across tests.
        overview = {'created': datetime.now(timezone.utc).isoformat(), 'analyses': job.get('analysis','proteomics'),
            'samples': read_rows(out/'sample_annotations.csv'), 'groups': [{k:g[k] for k in ('cohort','mode')} for g in job.get('groups',[])],
            'coverage': read_rows(out/'analysis_summary.csv'), 'results': index, 'files': files, 'feature_names': FEATURES,
            'drug_settings': drug_settings, 'main_settings': job.get('differential',{}),
            'optional_settings': {key:([kind for kind in value if kind in ('peptide','isoform')] if key=='features' else value)
                                 for key,value in job.get('optional_analyses',{}).get('statistics',{}).items()},
            'audit_path':'../analysis_details/reporting/result_audit/README.md'}
        template = Path(__file__).with_name('results_dashboard.html').read_text()
        payload = json.dumps(overview,ensure_ascii=False).replace('<','\\u003c').replace('&','\\u0026')
        (final/'START_HERE.html').write_text(template.replace('__RESULT_DATA__',payload))
        (final/'README.md').write_text('''# Proteomics results

Open **START_HERE.html** for the study overview, searchable main and optional test tables, and count matrices. The report works offline. **results_index.csv** provides the same test catalogue with file paths.

- **main_results/differential_expression/**: primary protein comparisons first, followed by supporting peptide and selected protein-modification tests. Each method has an all-results table and a significant-results table.
- **main_results/counts/**: per-cohort/search raw counts, normalized CPM, log2 CPM and visualization z-scores.
- **main_results/quality_control/**: sample assignments, coverage and statistical testing summaries.
- **matrices/**: direct links to CSV/TSV counts and available protein intensity/MaxLFQ matrices.
- **benchmark/**: matched service detection counts, gains/losses, complete-group and accession audits, and sample mappings.
- **optional_analyses/gene_counts/**: selected gene counts, normalized matrices, z-scores and comparisons, including selected gene-modification tests. Gene comparisons use DEPs settings and separate models/FDR; they do not change protein counts or protein statistics.
- **optional_analyses/statistics/**: separate peptide and isoform count comparisons and their matrices.
- **optional_analyses/drug_interaction/**: edgeR, limma-voom and/or Monte Carlo interaction results as selected.
- **optional_analyses/peptide_isoforms/**: sequence mapping, condition summaries and isoform evidence.
- **optional_analyses/immunopeptidomics/**: HLA peptide evidence when selected.

Only generated analyses appear. **No significant results** is a valid outcome. Full tables retain nonsignificant results; rows without estimable tests have blank p-values/FDR and an explicit status. **No eligible evidence** and **Insufficient evidence** are distinguished from a completed test with no discoveries. Empty significant tables are retained.

FDR is Benjamini–Hochberg adjusted separately for each method, cohort, search mode, feature type and comparison. Do not pool discoveries or select the smallest p-value across methods. DDA and DDA+ are separate searches of the samples, never biological replicates. Positive pairwise logFC means comparison minus reference. Overall condition tests have no single fold-change direction. Interaction log2FC is combination minus A minus B plus control; direction describes proteomic nonadditivity. Matching doses and times are recorded in interaction tables. Monte Carlo uncertainty concerns the simulated tail probability; its fitted nuisance parameters remain fixed.

Protein and peptide counts are spectral counts. Gene counts use unambiguous peptide assignments. Modification-by-gene/protein counts measure distinct modified peptide sequences per modification type, not occupancy. Isoform tests require discriminating peptides and alternative isoforms in the searched FASTA. HLA evidence does not measure source-protein abundance.

Model inputs, filtering masks, design matrices, normalization factors, software versions, settings and the verified file manifest remain under **../analysis_details/reporting/**. The layout step copies result tables without changing their statistics. Spectrum-search files and logs remain in the run's analysis_details and logs directories.
''')
        (audit/'README.md').write_text('''# Analysis records

The user-facing results are in ../../../final_results/START_HERE.html. This directory contains settings and supporting records moved out of the final result folders. Original table contents are preserved.

- main_statistics/: run settings, R session, detailed result indices and supporting records.
- optional_statistics/: optional count-model inputs, exclusion/filter reports, designs, normalization factors and settings.
- drug_interaction/: four-condition designs, fitted null models, voom weights, simulation parameters and feature status.
- peptide_isoforms/: FASTA mapping settings and provenance.
- file_manifest.csv: original file locations, delivered locations and SHA-256 checksums verified during copying.

Main model inputs remain in ../<cohort>/<mode>/. Historical settings describe the run that produced the results. Files with only headers and analyses with no significant features are valid outputs. Archived layouts, when requested, are in analysis_details/previous_result_layout_*.
''')
        (audit/'layout.json').write_text(json.dumps({'version':2,'tests':len(index),'created':overview['created'],'refitted':False},indent=2)+'\n')
        (final/'FINAL_RESULTS_COMPLETE').write_text('Final result tables and supporting files organized and verified.\n')
        build_matrix_catalog(final)
        # Both destination sets are ready before the existing results are replaced.
        stamp = datetime.now().strftime('%Y%m%d_%H%M%S_%f')
        old = out.parent/'analysis_details'/('previous_result_layout_'+stamp) if preserve_original else work/'original'
        old.parent.mkdir(parents=True,exist_ok=True)
        old_audit = old.with_name(old.name+'_audit') if preserve_original else work/'old_audit'
        audit_target.parent.mkdir(parents=True,exist_ok=True)
        out.rename(old)
        try:
            if audit_target.exists():audit_target.rename(old_audit)
            audit.rename(audit_target)
            final.rename(out)
        except Exception:
            if audit_target.exists() and not audit.exists():audit_target.rename(audit)
            if old_audit.exists():old_audit.rename(audit_target)
            old.rename(out)
            raise
        if preserve_original and (out.parent/'zscore').is_dir():
            (out.parent/'zscore').rename(old.with_name(old.name+'_zscore'))
        return {'tests':len(index),'result_files':sum(row['area']=='results' for row in manifest),
                'archive':str(old) if preserve_original else None,'audit':str(audit_target)}
    finally:
        shutil.rmtree(work,ignore_errors=True)


def upgrade_gene_layout(out, preserve_original=False):
    """Upgrade an already delivered report without changing any statistical CSV."""
    out = Path(out).resolve()
    index = read_rows(out/'results_index.csv')
    if not any(r['section']=='Main results' and r['feature_type'] in GENE_FEATURES for r in index) and not any(p.name.startswith(('gene_','modification_gene_')) for p in (out/'main_results/counts').rglob('*.csv')):
        return {'already_organized': True, 'tests':len(index)}
    audit = out.parent/'analysis_details/reporting/result_audit'
    work = Path(tempfile.mkdtemp(prefix='.result_layout_',dir=out.parent))
    stage = work/'final_results'
    try:
        shutil.copytree(out,stage)
        html = (out/'START_HERE.html').read_text()
        payload = json.loads(re.search(r'<script id="result-data" type="application/json">(.*?)</script>',html,re.S)[1])
        moved = {}
        def move(old, new):
            if old in moved:return
            source, target = safe_source(stage,old), safe_source(stage,new)
            if target.exists():raise ValueError('Result path already exists: '+new)
            target.parent.mkdir(parents=True,exist_ok=True)
            if digest(source)!=digest(safe_source(out,old)):raise ValueError('Result changed during layout upgrade')
            source.rename(target);moved[old]=new
        for row in index:
            if row['section']=='Main results' and row['feature_type'] in GENE_FEATURES:
                row.update(section='Optional analyses',analysis='Gene counts and comparisons')
                for key in ('full_results','significant_results'):
                    new=row[key].replace('main_results/differential_expression/','optional_analyses/gene_counts/comparisons/',1)
                    move(row[key],new);row[key]=new
        for path in (stage/'main_results/counts').rglob('*.csv'):
            if path.name.startswith(('gene_','modification_gene_')):
                old=str(path.relative_to(stage));move(old,old.replace('main_results/counts/','optional_analyses/gene_counts/counts/',1))
        for entry in payload['files']:
            if entry['path'] in moved:
                entry.update(path=moved[entry['path']],category='Gene count matrices')
        index=ordered_results(index);payload['results']=index
        payload['created']=datetime.now(timezone.utc).isoformat()
        for path in sorted(stage.rglob('*'),key=lambda p:len(p.parts),reverse=True):
            if path.is_dir() and not any(path.iterdir()):path.rmdir()
        write_rows(stage/'results_index.csv',INDEX_FIELDS,index)
        encoded=json.dumps(payload,ensure_ascii=False).replace('<','\\u003c').replace('&','\\u0026')
        (stage/'START_HERE.html').write_text(Path(__file__).with_name('results_dashboard.html').read_text().replace('__RESULT_DATA__',encoded))
        readme=(stage/'README.md').read_text()
        readme=readme.replace('- **main_results/differential_expression/**: named condition comparisons and overall condition tests',
                              '- **main_results/differential_expression/**: primary protein results and supporting peptide/protein-modification tests')
        (stage/'README.md').write_text(readme+'\nProtein analyses are primary. Gene counts, matrices and comparisons (including gene modifications) are under **optional_analyses/gene_counts/**. Only their locations changed; statistics were not recomputed.\n')
        manifest=read_rows(audit/'file_manifest.csv')
        for row in manifest:
            if row['area']=='results' and row['file'] in moved:row['file']=moved[row['file']]
            target=(audit if row['area']=='audit' else stage)/row['file']
            if digest(target)!=row['sha256']:raise ValueError('File integrity check failed: '+str(target))
        stamp=datetime.now().strftime('%Y%m%d_%H%M%S_%f')
        archive=out.parent/'analysis_details'/('previous_result_layout_'+stamp) if preserve_original else work/'original'
        archive.parent.mkdir(parents=True,exist_ok=True)
        out.rename(archive)
        try:stage.rename(out)
        except Exception:
            archive.rename(out);raise
        write_rows(audit/'file_manifest.csv',['area','file','source','sha256','bytes'],manifest)
        (audit/'layout.json').write_text(json.dumps({'version':2,'tests':len(index),'created':payload['created'],'refitted':False},indent=2)+'\n')
        return {'tests':len(index),'moved_files':len(moved),'archive':str(archive) if preserve_original else None}
    finally:
        shutil.rmtree(work,ignore_errors=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('out',help='Completed final_results directory')
    parser.add_argument('--preserve-original',action='store_true')
    args=parser.parse_args()
    print(json.dumps(publish_results(args.out,args.preserve_original),indent=2))
