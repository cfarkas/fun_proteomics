"""Validate biological grouping before raw BAF search."""
from pathlib import Path
import re
import shutil
import subprocess
import sys


def biological_samples(frame, subject=None):
    rows = frame.copy()
    if subject and subject not in rows: raise ValueError('Enter the donor/subject ID for each BAF acquisition')
    if 'biological_sample' not in rows: rows['biological_sample'] = rows.sample_id
    if rows.biological_sample.eq('').any(): raise ValueError('Assign a biological sample to every technical injection')
    keep = ['condition', *(['cohort'] if 'cohort' in rows else []), *([subject] if subject else []),
            *(['service_sample'] if 'service_sample' in rows else [])]
    for name, group in rows.groupby('biological_sample', sort=False):
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', name): raise ValueError('Biological sample IDs need letters, numbers, dot, dash or underscore')
        if any(group[key].nunique(dropna=False) != 1 for key in keep):
            raise ValueError('Technical injections disagree on condition/donor/service mapping for '+name)
    result = rows.drop_duplicates('biological_sample').copy()
    result['sample_id'] = result.biological_sample
    result['replicate'] = result[subject] if subject else [str(i+1) for i in range(len(result))]
    return result.drop(columns=['raw_path','biological_sample'], errors='ignore')


def converter_python(cfg, workflow):
    path = cfg.get('python') or str(workflow/'.runtime/baf-python/bin/python')
    if not Path(path).is_file():
        raise ValueError('BAF reader missing. Run fun_proteomics install --baf, or set python to a runtime with pyBaf2Sql and psims.')
    probe=subprocess.run([path,'-c','import psims; from pyBaf2Sql.init_baf2sql import init_baf2sql_api; init_baf2sql_api()'],capture_output=True,text=True,timeout=30)
    if probe.returncode: raise ValueError('BAF runtime cannot load psims/Baf2Sql: '+probe.stderr[-1000:])
    return path


def check(cfg, frame, workflow):
    from .study import imports
    imports(workflow)
    from baf import source_files
    from quantification_config import differential_options, validate_design_rows
    from quantification import load_matrix
    from pipeline import digest, fasta_records
    if frame.cohort.nunique()!=1: raise ValueError('Use one BAF project per cohort')
    for row in frame.to_dict('records'):
        if Path(cfg['project']).resolve().is_relative_to(Path(row['raw_path']).resolve()):
            raise ValueError('The output project must be outside every source acquisition directory')
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', row['sample_id']): raise ValueError('Acquisition IDs need letters, numbers, dot, dash or underscore')
        for p in source_files(row['raw_path']):
            if not p.stat().st_size: raise ValueError('Empty BAF input: '+str(p))
    bio = biological_samples(frame, cfg.get('subject_column'))
    options = differential_options(dict(enabled=cfg['differential'], comparisons='all_pairwise',
        fdr=cfg['fdr'], min_replicates=cfg['min_replicates'], normalization=cfg['normalization'],
        **({'subject_column':cfg['subject_column']} if cfg.get('subject_column') else {})))
    validate_design_rows(bio.to_dict('records'), options)
    if not cfg.get('fasta') or not Path(cfg['fasta']).is_file(): raise ValueError('Select a target/decoy search FASTA')
    records = list(fasta_records(Path(cfg['fasta']).read_text()))
    targets = {h.split()[0] for h,_ in records if not h.startswith('rev_')}
    decoys = {h.split()[0][4:] for h,_ in records if h.startswith('rev_')}
    if not targets or targets != decoys or not any(h.startswith('contam_') for h in targets):
        raise ValueError('BAF requires one matching rev_ decoy per target and contam_ entries')
    if len({h.split()[0] for h,_ in records}) != len(records): raise ValueError('Duplicate FASTA identifiers')
    from .continuation import read_json
    versions = read_json(workflow/'assets/tool_versions.json')
    for tool,key in [('MSFragger','msfragger_zip'),('IonQuant','ionquant_zip'),('diaTracer','diatracer_zip')]:
        if not cfg.get(key) or not Path(cfg[key]).is_file(): raise ValueError('Select registered '+key)
        if digest(cfg[key]) != versions[tool]['archive_sha256']: raise ValueError('Unexpected '+tool+' archive version')
    converter_python(cfg, workflow)
    if not (workflow/'.runtime/env/bin/java').is_file() and not shutil.which('java'):
        raise ValueError('Java is missing; run fun_proteomics install --conda before starting conversion')
    if not 1 <= int(cfg['cpus']) <= 512 or not 4 <= int(cfg['memory_gib']) <= 4096: raise ValueError('Invalid CPU/RAM allocation')
    if cfg.get('service_matrix'):
        from .assignments import service_samples_for
        load_matrix(cfg['service_matrix'], service_samples_for(cfg,bio,workflow), sheet=cfg.get('service_sheet'))
    conditions = bio.condition.nunique(); pairs = conditions*(conditions-1)//2 if cfg['differential'] else 0
    return dict(mode='baf', samples=len(bio), acquisitions=len(frame), comparisons=pairs,
        output_folder=cfg['project'], results_folder=str(Path(cfg['project'])/'results/final_results'), warnings=[],
        summary=f'{len(frame)} BAF acquisitions → {len(bio)} biological samples; {pairs} condition pairs. Technical injections share one FragPipe experiment; donor pairing is preserved when specified.')
