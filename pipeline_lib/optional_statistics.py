"""Prepare optional count evidence for the existing differential model and reports."""
import collections
import csv
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile


CATEGORIES = {'peptide': 'peptides', 'isoform': 'isoforms'}
COUNT_UNITS = {'peptide': 'peptide spectral counts',
               'isoform': 'spectral counts from isoform-discriminating peptides'}


def validate_statistics(plan):
    options = plan.get('optional_analyses', {})
    settings = options.get('statistics', {})
    if not isinstance(settings, dict) or not isinstance(settings.get('enabled', False), bool):
        raise ValueError('Optional statistics enabled must be true or false.')
    if not settings.get('enabled'):
        return
    features, methods = settings.get('features', []), settings.get('methods', [])
    if (not isinstance(features, list) or not features or any(f not in CATEGORIES for f in features)
            or len(set(features)) != len(features)):
        raise ValueError('Select valid evidence types for optional statistical comparisons.')
    if (not isinstance(methods, list) or not methods or any(m not in ('edgeR', 'limma') for m in methods)
            or len(set(methods)) != len(methods)):
        raise ValueError('Optional statistical methods must include edgeR and/or limma.')
    hla = plan.get('proteomics', {}).get('analysis') == 'immunopeptidomics'
    allowed = {'peptide': options.get('peptide_isoforms') or hla,
               'isoform': options.get('peptide_isoforms')}
    if any(not allowed[f] for f in features):
        raise ValueError('Enable the parent optional analysis for each selected statistical evidence type.')
    fdr = settings.get('fdr', 0.05)
    if isinstance(fdr, bool) or not isinstance(fdr, (int, float)) or not math.isfinite(fdr) or not 0 < fdr <= 1:
        raise ValueError('Optional-analysis FDR must be greater than 0 and at most 1.')
    if not isinstance(settings.get('filter_low_counts', False), bool):
        raise ValueError('Optional-analysis count filtering must be true or false.')
    if settings.get('filter_low_counts'):
        for key, fallback in [('min_count', 1), ('min_total_count', 3)]:
            value = settings.get(key, fallback)
            if (isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value)
                    or value < 1 or (key == 'min_total_count' and value != int(value))):
                raise ValueError('Optional-analysis minimum counts must be positive; total counts must be integers.')
    comparison, control = settings.get('comparisons', 'all_pairwise'), settings.get('control', '')
    if comparison not in ('control', 'all_pairwise') or (comparison == 'control' and not control):
        raise ValueError('Choose a reference condition or all pairs for optional statistics.')
    cohorts = collections.defaultdict(lambda: collections.defaultdict(list))
    for sample in plan.get('samples', []):
        cohorts[sample['cohort']][sample['condition']].append(sample['replicate'])
    if not cohorts:
        raise ValueError('Optional statistics require samples.')
    for cohort, groups in cohorts.items():
        if len(groups) < 2 or any(len(reps) < 2 or len(set(reps)) != len(reps) for reps in groups.values()):
            raise ValueError('Optional statistics require two conditions and two independent biological replicates per condition: '+cohort)
        if comparison == 'control' and control not in groups:
            raise ValueError('Optional-statistics reference condition is missing from cohort '+cohort+': '+control)


def read_table(path, delimiter=','):
    with Path(path).open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        return reader.fieldnames or [], list(reader)


def write_table(path, fields, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def count(value):
    try:
        parsed = int(value)
    except (ValueError, TypeError):
        raise ValueError('Optional statistics require nonnegative integer counts: '+str(value)) from None
    if parsed < 0 or str(parsed) != str(value).strip():
        raise ValueError('Optional statistics require nonnegative integer counts: '+str(value))
    return parsed


def isoform_counts(folder, ids, exclusions):
    source = Path(folder)/'isoform_evidence.csv'
    fields, rows = read_table(source)
    if not set(['Accession', 'Family_entries_in_FASTA', 'Isoform_discriminating_peptides']+ids) <= set(fields):
        raise ValueError('Isoform evidence lacks sample counts or discrimination annotations.')
    accepted, rejected = [], []
    for row in rows:
        if count(row['Family_entries_in_FASTA']) < 2:
            reason = 'No alternative isoform sequence in the searched FASTA'
        elif count(row['Isoform_discriminating_peptides']) == 0:
            reason = 'No observed peptide uniquely distinguishes this isoform'
        else:
            reason = ''
        if reason:
            rejected.append({**row, 'Exclusion_reason': reason})
        else:
            accepted.append({'FeatureID': row['Accession'], **row, **{sid: count(row[sid]) for sid in ids}})
    write_table(exclusions, fields+['Exclusion_reason'], rejected)
    return ['FeatureID']+fields, accepted, source


def build_statistics(plan, final_out, out, rscript, env):
    from pipeline import digest
    from result_reports import organize
    settings = dict(plan['optional_analyses']['statistics'])
    selected = settings['features']
    categories = {kind: CATEGORIES[kind] for kind in selected}
    manifest, provenance = [], []
    combined = {kind: [] for kind in selected}
    annotations = {}
    all_ids = [sample['sample_id'] for sample in plan['samples']]
    for cohort in sorted({sample['cohort'] for sample in plan['samples']}):
        samples = [sample for sample in plan['samples'] if sample['cohort'] == cohort]
        ids = [sample['sample_id'] for sample in samples]
        for mode in plan['modes']:
            slug = 'ddaplus' if mode == 'DDA+' else 'dda'
            raw = out/'inputs'/cohort/slug
            for kind in selected:
                if kind == 'isoform':
                    fields, rows, source = isoform_counts(final_out/'optional_analyses/peptide_isoforms'/cohort/slug, ids, raw/'excluded_isoforms.csv')
                else:
                    source = final_out.parent/'analysis_details/reporting'/cohort/slug/'peptide_counts.csv'
                    fields, rows = read_table(source)
                if not set(['FeatureID']+ids) <= set(fields):
                    raise ValueError('Missing sample count columns in optional evidence: '+str(source))
                if len({row['FeatureID'] for row in rows}) != len(rows):
                    raise ValueError('Duplicate optional feature IDs in '+str(source))
                for row in rows:
                    for sid in ids:
                        row[sid] = count(row[sid])
                write_table(raw/(kind+'_counts.csv'), fields, rows)
                ann = [field for field in fields if field not in ids and field != 'FeatureID']
                annotations[kind] = list(dict.fromkeys(annotations.get(kind, [])+ann))
                combined[kind].extend({'Cohort': cohort, 'Mode': mode, **row, **{sid: row.get(sid, '') for sid in all_ids}} for row in rows)
                provenance.append({'cohort': cohort, 'mode': mode, 'feature_type': kind,
                                   'source': os.path.relpath(source, out), 'sha256': digest(source), 'features': len(rows)})
            write_table(raw/'samples.csv', ['sample_id', 'condition', 'replicate', 'cohort'],
                        [{key: sample[key] for key in ('sample_id', 'condition', 'replicate', 'cohort')} for sample in samples])
            manifest.append({'cohort': cohort, 'mode': mode, 'directory': str(raw.resolve())})
    for kind in selected:
        write_table(out/(kind+'_counts.csv'), ['Cohort', 'Mode', 'FeatureID']+annotations[kind]+all_ids, combined[kind])
    # All optional evidence uses sparse-count modeling and preserves full-library offsets.
    options = {**settings, 'modification_min_count': settings.get('min_count', 1),
               'modification_min_total_count': settings.get('min_total_count', 3)}
    options.pop('features')
    job = {'groups': manifest, 'outdir': str(out.resolve()), 'analysis': 'optional-count-statistics',
           'parent_analysis': plan.get('proteomics', {}).get('analysis', 'global-label-free'),
           'differential': options, 'feature_categories': categories,
           'sparse_kinds': selected, 'count_units': {kind: COUNT_UNITS[kind] for kind in selected},
           'retain_untested_rows': True, 'zscore_directory': str((out/'zscore').resolve())}
    jobfile = out/'analysis_settings.json'
    jobfile.write_text(json.dumps(job, indent=2)+'\n')
    (out/'input_sources.json').write_text(json.dumps(provenance, indent=2)+'\n')
    subprocess.run([rscript, str(Path(__file__).with_name('differential.R')), str(jobfile)], check=True, env=env)
    organize(out)
    job['outdir'] = '.'
    job['zscore_directory'] = 'zscore'
    for group in job['groups']:
        group['directory'] = os.path.relpath(group['directory'], out)
    jobfile.write_text(json.dumps(job, indent=2)+'\n')
    (out/'README.md').write_text('''# Statistical comparisons for optional analyses

Open `START_HERE.html` for the test index. `analysis_settings.json` records the selected evidence, methods, reference/all-pairs comparison, FDR and filtering thresholds. `input_sources.json` records the exact source table hashes. These settings run independently of the main DEPs switch. Drug-combination interaction has its own four-condition model and FDR in `../drug_interaction/` when selected.

The engine fits a joint condition model separately for each feature type, cohort and search mode. edgeR uses its small-count quasi-likelihood GLM with robust moderation; optional limma-voom provides a comparison method. TMMwsp normalization and library sizes calculated before count filtering are retained. The edgeR overall condition test is included; selected methods also produce named pairwise comparisons. Samples are independent biological replicates. There is no inferred paired-subject or batch adjustment. DDA and DDA+ are separate analyses, never biological replicates.

`differential/` contains full and significant result tables, a test index and condition-name mapping. Pairwise logFC is comparison minus reference. `PValue` is unadjusted; `FDR` is Benjamini–Hochberg adjusted within each feature type, method, contrast, cohort and search mode. The FDR cutoff is recorded in settings and the HTML guide. This tests condition differences; it is separate from peptide, protein, isoform or modification identification confidence. Inspect effect sizes, condition means and detection counts alongside p-values.

Peptide features use primary-search spectral counts. Isoform features sum only peptides unique across the searched target FASTA to an accession whose family contains alternative isoforms, with I/L equivalence. Accessions lacking alternatives or discriminating peptides are listed in `inputs/<cohort>/<mode>/excluded_isoforms.csv`; none receive inferred isoform-specific counts or p-values. Canonical-only FASTA files can therefore produce empty isoform tests. These are comparisons of isoform-associated evidence, not tests of within-gene isoform usage or switching. HLA counts describe peptide evidence rather than source-protein abundance.

Low coverage is not a run failure. Full tables retain nonsignificant results. All-zero, filtered or unestimable rows have explicit test_status and blank p-values/FDR; they are not evidence of no effect. No pseudocounts or imputation enter the count models. Empty significant subsets and header-only tests with no eligible evidence are valid outcomes. Filtering reports, design/contrast matrices and normalization factors are under `inputs/`; normalized matrices and visualization-only z-scores have their own folders.

Methods: https://bioconductor.org/packages/release/bioc/html/edgeR.html
''')


def run_statistics(plan, final_out, rscript, env):
    validate_statistics(plan)
    if not plan.get('optional_analyses', {}).get('statistics', {}).get('enabled'):
        return
    final_out = Path(final_out).resolve()
    parent = final_out/'optional_analyses'
    parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix='.statistics_', dir=parent))
    build_statistics(plan, final_out, stage, rscript, env)
    target = parent/'statistics'
    backup = stage.with_name(stage.name+'_previous')
    if target.exists():
        target.rename(backup)
    try:
        stage.rename(target)
    except Exception:
        if backup.exists():
            backup.rename(target)
        raise
    if backup.exists():
        shutil.rmtree(backup)
