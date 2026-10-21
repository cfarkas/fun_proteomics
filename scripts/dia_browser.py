"""Translate browser studies into independent native DIA-NN input files.

The browser study remains the editable record. DDA settings are retained there
but never translated into DIA scientific settings.
"""
import copy
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

from diapasef_config import normalize_config, path_value
from fasta_download import read_records, accession

DIA_FIELDS = {'confirmed', 'canonical_fasta', 'contaminant_fasta', 'install_cache',
              'python', 'raw_root', 'library_reference_ids', 'differential', 'service_matrices'}
NOTICE = ('DIA-NN uses the bundled DIA method and produces identification and quantification '
          'protein and gene matrices, optional limma intensity comparisons and service benchmarks. '
          'DDA search and count-model settings remain saved separately.')


def is_dia_study(cfg):
    return cfg.get('proteomics', {}).get('acquisition') == 'DIA-PASEF'


def python312():
    candidates = [shutil.which('python3.12'), '/usr/bin/python3.12',
                  '/usr/bin/python3', shutil.which('python3')]
    for candidate in dict.fromkeys(p for p in candidates if p):
        try:
            version = subprocess.check_output([candidate, '-c',
                'import sys; print("%d.%d" % sys.version_info[:2])'], text=True,
                stderr=subprocess.DEVNULL, timeout=10).strip()
            if version == '3.12':
                return candidate
        except (OSError, subprocess.SubprocessError):
            pass
    raise ValueError('DIA-NN requires Python 3.12. Set the DIA Python executable path.')


def _targets(study, options, root):
    database = study.get('database', {})
    custom = bool(options.get('canonical_fasta'))
    fasta = path_value(options.get('canonical_fasta') or database.get('fasta'), root, 'file')
    digest = hashlib.sha256(Path(fasta).read_bytes()).hexdigest()
    if not custom and (not database.get('sha256') or digest != database['sha256']):
        raise ValueError('The saved FASTA does not match its SHA256 checksum.')
    proteins = read_records(fasta)
    decoys = [(h, s) for h, s in proteins if h.startswith('rev_')]
    targets = [(h, s) for h, s in proteins if not h.startswith('rev_')]
    if custom and decoys:
        raise ValueError('The DIA target-only FASTA must not contain decoys.')
    if not custom and database.get('contains_decoys') is not bool(decoys):
        raise ValueError('FASTA decoy metadata disagrees with its contents; only the pipeline rev_ prefix is supported.')
    if decoys and {h.split()[0] for h, _ in targets} != {h.split()[0][4:] for h, _ in decoys}:
        raise ValueError('The saved FASTA requires exactly one rev_ decoy per target identifier.')
    if not targets:
        raise ValueError('FASTA contains no target proteins.')
    for header, _ in targets:
        if any(re.match(r'(?i)^(?:rev(?:erse)?|decoy|xxx)[_|-]', part)
               for part in header.split()[0].split('|')):
            raise ValueError('Unrecognized decoy prefix in target FASTA; use a verified target-only FASTA.')
    # Preserve existing contaminant targets; add an explicitly configured local
    # contaminant archive using the same parsing rules as the DDA workflow.
    if not custom and database.get('contaminants_zip'):
        contamination = path_value(database['contaminants_zip'], root, 'file')
        observed = hashlib.sha256(Path(contamination).read_bytes()).hexdigest()
        if database.get('contaminants_sha256') and observed != database['contaminants_sha256']:
            raise ValueError('Contaminant archive does not match its saved SHA256 checksum.')
        seen = {accession(h) for h, _ in targets}
        sequences = {sequence for _, sequence in targets}
        for header, sequence in read_records(contamination, zipped=True):
            if header.startswith('rev_'):
                raise ValueError('Contaminant input must not contain decoys.')
            header = header if header.startswith('contam_') else 'contam_' + header
            if accession(header) not in seen and sequence not in sequences:
                targets.append((header, sequence))
                seen.add(accession(header))
                sequences.add(sequence)
    return targets, {'source_fasta': fasta, 'source_sha256': digest,
                     'removed_rev_decoys': len(decoys), 'target_entries': len(targets)}


def browser_inputs(study, root):
    """Validate browser-only choices and construct inputs without modifying study."""
    root = Path(root).resolve()
    choice = study.get('proteomics', {})
    if not is_dia_study(study) or choice.get('confirmed') is not True:
        raise ValueError('Confirm DIA-PASEF in the Proteomics step before preparing DIA inputs.')
    if choice.get('analysis') != 'global-label-free':
        raise ValueError('The DIA browser adapter supports global label-free proteomics only.')
    options = study.get('dia', {})
    if not isinstance(options, dict) or set(options) - DIA_FIELDS:
        raise ValueError('Unknown DIA browser settings: ' + ', '.join(sorted(set(options) - DIA_FIELDS)) if isinstance(options, dict) else 'DIA settings must be an object')
    if options.get('confirmed') is not True:
        raise ValueError('Confirm the DIA method and supported outputs in the DIA settings. ' + NOTICE)
    rows = []
    if not study.get('samples'):
        raise ValueError('Select samples before preparing DIA inputs.')
    raw_base = Path(path_value(options.get('raw_root', str(root)), root, 'directory'))
    for sample in study['samples']:
        row = {key: sample[key] for key in ('sample_id', 'condition', 'replicate')}
        if any(not str(row[key]).strip() or any(c in str(row[key]) for c in '\t\r\n') for key in row):
            raise ValueError('Sample names and annotations must be nonempty and cannot contain tabs or line breaks.')
        row['sample_type'] = sample['cohort']
        if sample.get('service_sample'): row['service_sample'] = sample['service_sample']
        if row['sample_type'] == 'unassigned':
            raise ValueError('Assign a cohort for ' + str(row['sample_id']) + ' before running.')
        row['raw_path'] = path_value(sample['raw_path'], raw_base, 'directory')
        if 'library_reference' in sample and not isinstance(sample['library_reference'], bool):
            raise ValueError('library_reference must be true or false.')
        rows.append(row)
    raw_root = Path(os.path.commonpath([str(Path(row['raw_path']).parent) for row in rows]))
    if options.get('raw_root'):
        raw_root = raw_base
        if any(not Path(row['raw_path']).is_relative_to(raw_root) for row in rows):
            raise ValueError('Every raw acquisition must be inside the selected DIA raw root.')
    if raw_root == Path(raw_root.anchor):
        raise ValueError('Select DIA samples under one shared raw-data directory.')
    explicit = options.get('library_reference_ids')
    per_sample = None
    if explicit is None and any('library_reference' in sample for sample in study['samples']):
        per_sample = {(sample['cohort'], sample['sample_id']) for sample in study['samples'] if sample.get('library_reference') is True}
    if explicit is not None and (not isinstance(explicit, list) or any(not isinstance(v, str) for v in explicit)):
        raise ValueError('DIA library_reference_ids must be a list of sample names.')
    if explicit is not None and set(explicit) - {r['sample_id'] for r in rows}:
        raise ValueError('DIA library reference names do not match selected samples.')
    selected = set(explicit or [])
    if explicit is None and per_sample is None:
        for group in sorted({r['sample_type'] for r in rows}):
            candidates = sorted((r for r in rows if r['sample_type'] == group), key=lambda r: (r['raw_path'], r['sample_id']))
            selected.update((group, row['sample_id']) for row in candidates[:20])
    for row in rows:
        flag = (row['sample_type'], row['sample_id']) in per_sample if per_sample is not None else row['sample_id'] in selected if explicit is not None else (row['sample_type'], row['sample_id']) in selected
        row['library_reference'] = str(flag).lower()
    targets, provenance = _targets(study, options, root)
    resources = study.get('resources', {})
    cfg = {'acquisition': 'DIA-PASEF', 'raw_root': str(raw_root),
           'install_cache': path_value(options.get('install_cache', str(root / '.runtime/dia-cache')), root),
           'outdir': path_value(study.get('outdir'), root), 'cpus': resources.get('cpus', 32),
           'memory_gib': resources.get('memory_gib', 96), 'expected_samples': len(rows),
           'python': options.get('python') or python312(),
           'method_config': str(root / 'modules/diapasef/method.toml'),
           'requirements_lock': str(root / 'modules/diapasef/requirements-lock.txt')}
    cfg['differential'] = copy.deepcopy(options.get('differential', {'enabled': False}))
    cfg['service_matrices'] = {group: path_value(path, root, 'file') for group, path in options.get('service_matrices', {}).items()}
    if options.get('contaminant_fasta'):
        cfg['contaminant_fasta'] = path_value(options['contaminant_fasta'], root, 'file')
    return cfg, rows, targets, provenance


def _write_inputs(cfg, rows, path, targets, provenance, logical_path=None):
    logical_path = logical_path or path
    fasta = path.with_suffix('.targets.fasta')
    with fasta.open('x') as handle:
        for header, sequence in targets:
            handle.write(f'>{header}\n{sequence}\n')
    with path.with_suffix('.fasta-provenance.json').open('x') as handle:
        json.dump(dict(provenance, target_sha256=hashlib.sha256(fasta.read_bytes()).hexdigest()), handle, indent=2)
        handle.write('\n')
    sheet = path.with_suffix('.samples.tsv')
    with sheet.open('x', newline='') as handle:
        writer = csv.DictWriter(handle, delimiter='\t', fieldnames=['sample_id', 'sample_type', 'condition', 'replicate', 'raw_path', 'library_reference', 'service_sample'])
        writer.writeheader()
        writer.writerows(rows)
    native = dict(cfg, samplesheet=str(logical_path.with_suffix('.samples.tsv')), canonical_fasta=str(logical_path.with_suffix('.targets.fasta')))
    with path.open('x') as handle:
        json.dump(native, handle, indent=2)
        handle.write('\n')
    return native


def validate_browser_study(study, root):
    cfg, rows, targets, provenance = browser_inputs(study, root)
    # Validation writes temporary manifests only; it never changes the study,
    # starts Nextflow, downloads tools, or creates the eventual output folder.
    with tempfile.TemporaryDirectory(prefix='pasef-dia-check-') as directory:
        folder = Path(directory)
        validation = copy.deepcopy(cfg)
        if not Path(validation['install_cache']).exists():
            validation['install_cache'] = str(folder)
        path = folder / 'dia.json'
        _write_inputs(validation, rows, path, targets, provenance)
        try:
            normalized = normalize_config(path, Path(root) / 'modules/diapasef')
        except subprocess.SubprocessError as error:
            raise ValueError('The selected DIA Python executable could not run its version check. Select Python 3.12.') from error
    return {'samples': len(rows), 'groups': normalized['groups'], 'notice': NOTICE,
            'library_references': {group: [r['sample_id'] for r in rows if r['sample_type'] == group and r['library_reference'] == 'true'] for group in normalized['groups']},
            'fasta': provenance}


def prepare_browser_study(study_path, root):
    """Write content-addressed native inputs, preserving browser JSON verbatim."""
    study_path = Path(study_path).resolve()
    study = json.loads(study_path.read_text())
    cfg, rows, targets, provenance = browser_inputs(study, root)
    validate_browser_study(study, root)
    Path(cfg['install_cache']).mkdir(parents=True, exist_ok=True)
    fingerprint = hashlib.sha256(json.dumps([cfg, rows, targets, provenance], sort_keys=True).encode()).hexdigest()
    parent = Path(root).resolve() / 'ui_configs/dia'
    parent.mkdir(parents=True, exist_ok=True)
    destination = parent / fingerprint
    path = destination / 'dia.json'
    with tempfile.TemporaryDirectory(prefix='.prepare-', dir=parent) as directory:
        bundle = Path(directory) / 'bundle'
        bundle.mkdir()
        _write_inputs(cfg, rows, bundle / 'dia.json', targets, provenance, logical_path=path)
        if destination.exists():
            for expected in bundle.iterdir():
                actual = destination / expected.name
                if not actual.is_file() or actual.read_bytes() != expected.read_bytes():
                    raise ValueError('Saved DIA inputs changed; restore the original generated input bundle before resuming.')
        else:
            bundle.rename(destination)
    return path
