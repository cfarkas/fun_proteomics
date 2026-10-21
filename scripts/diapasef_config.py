#!/usr/bin/env python3
"""Validate the DIA-NN branch before scheduling spectrum searches."""
import argparse
import csv
import json
import os
from pathlib import Path
import re
import shutil
import sqlite3
import subprocess
import tomllib
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pipeline_lib"))
from quantification_config import audit_service, differential_options, validate_design_rows

FIELDS = {
    'acquisition', 'samplesheet', 'raw_root', 'canonical_fasta',
    'contaminant_fasta', 'method_config', 'requirements_lock', 'install_cache',
    'differential', 'service_matrices', 'expected_samples', 'cpus', 'memory_gib', 'python', 'outdir',
}


def integer(value, name, minimum=1):
    if isinstance(value, bool) or not re.fullmatch(r'[0-9]+', str(value)) or int(value) < minimum:
        raise ValueError(f'{name} must be an integer >= {minimum}')
    return int(value)


def path_value(value, base, kind=None):
    if not isinstance(value, str) or not value.strip() or any(c in value for c in '\t\n\r'):
        raise ValueError('Paths must be nonempty strings without tabs or line breaks')
    path = (base / Path(value).expanduser()).resolve()
    if kind == 'file' and not path.is_file():
        raise ValueError(f'Missing file: {path}')
    if kind == 'directory' and not path.is_dir():
        raise ValueError(f'Missing directory: {path}')
    return str(path)


def audit_samples(sheet, raw_root, expected, services):
    with Path(sheet).open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        required = {'sample_id', 'sample_type', 'raw_path', 'library_reference'}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f'Sample table requires {sorted(required)}')
        rows = list(reader)
    if not rows or (expected and len(rows) != expected):
        raise ValueError(f'Expected {expected or "at least one"} samples; found {len(rows)}')
    seen_ids, seen_paths, seen_runs, groups = set(), set(), set(), {}
    for row in rows:
        if any(not row.get(key) for key in required):
            raise ValueError('Missing required sample value')
        group, sid = row['sample_type'], row['sample_id']
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*', group):
            raise ValueError(f'Invalid sample type: {group}')
        if (group, sid) in seen_ids:
            raise ValueError(f'Duplicate sample ID within type: {group}/{sid}')
        seen_ids.add((group, sid))
        flag = row['library_reference'].lower()
        if flag not in ('true', 'false'):
            raise ValueError('library_reference must be true or false')
        raw = Path(path_value(row['raw_path'], Path(raw_root), 'directory'))
        if raw.suffix.lower() != '.d' or 'BLANK' in raw.name.upper():
            raise ValueError(f'Expected a biological Bruker .d acquisition: {raw}')
        if raw in seen_paths or (group, raw.name) in seen_runs:
            raise ValueError(f'Duplicate raw acquisition: {raw}')
        seen_paths.add(raw)
        seen_runs.add((group, raw.name))
        for filename in ('analysis.tdf', 'analysis.tdf_bin'):
            component = raw / filename
            if not component.is_file() or not component.stat().st_size:
                raise ValueError(f'Incomplete raw acquisition: {component}')
        with sqlite3.connect((raw / 'analysis.tdf').as_uri() + '?mode=ro', uri=True) as connection:
            types = {r[0] for r in connection.execute('SELECT DISTINCT MsMsType FROM Frames')}
        if types != {0, 9}:
            raise ValueError(f'{group}/{sid}: expected DIA-PASEF frame types 0 and 9; found {sorted(types)}')
        entry = groups.setdefault(group, {'samples': 0, 'library_references': 0, 'runs': []})
        entry['samples'] += 1
        entry['library_references'] += flag == 'true'
        entry['runs'].append(raw.name)
    unknown = set(services) - set(groups)
    if unknown:
        raise ValueError(f'Service matrices configured for absent types: {sorted(unknown)}')
    for group, info in groups.items():
        if not info['library_references']:
            raise ValueError(f'At least one library reference is required for {group}')
        if group in services:
            audit_service(services[group], [row for row in rows if row['sample_type'] == group])
        del info['runs']
    return groups


def normalize_config(config, module_dir, threads=None, memory=None, outdir=None):
    config = Path(config).resolve()
    cfg = json.loads(config.read_text())
    if not isinstance(cfg, dict):
        raise ValueError('DIA configuration must be a JSON object')
    if set(cfg) - FIELDS:
        raise ValueError(f'Unknown DIA configuration fields: {sorted(set(cfg) - FIELDS)}')
    if cfg.get('acquisition', 'DIA-PASEF') != 'DIA-PASEF':
        raise ValueError('The DIA-NN branch requires acquisition DIA-PASEF')
    base, module_dir = config.parent, Path(module_dir).resolve()
    for name, kind in [('samplesheet', 'file'), ('raw_root', 'directory'),
                       ('canonical_fasta', 'file'), ('install_cache', 'directory')]:
        cfg[name] = path_value(cfg.get(name), base, kind)
    for name, default in [('method_config', 'method.toml'), ('requirements_lock', 'requirements-lock.txt')]:
        cfg[name] = path_value(cfg.get(name, str(module_dir / default)), base, 'file')
    cfg['contaminant_fasta'] = path_value(cfg['contaminant_fasta'], base, 'file') if cfg.get('contaminant_fasta') else None
    services = cfg.get('service_matrices', {})
    if not isinstance(services, dict):
        raise ValueError('service_matrices must map sample types to matrix paths')
    cfg['service_matrices'] = {key: path_value(value, base, 'file') for key, value in services.items()}
    cfg['outdir'] = path_value(outdir, Path.cwd()) if outdir else path_value(cfg.get('outdir'), base)
    output, raw_root = Path(cfg['outdir']), Path(cfg['raw_root'])
    if output == raw_root or output.is_relative_to(raw_root) or raw_root.is_relative_to(output):
        raise ValueError('The output directory and raw root must be separate directory trees')
    cfg['cpus'] = integer(threads if threads is not None else cfg.get('cpus', 32), 'cpus')
    available = len(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else os.cpu_count()
    if cfg['cpus'] > available:
        raise ValueError(f'Requested {cfg["cpus"]} CPUs; only {available} available')
    cfg['memory_gib'] = integer(memory if memory is not None else cfg.get('memory_gib', 96), 'memory_gib', 4)
    cfg['expected_samples'] = integer(cfg.get('expected_samples', 0), 'expected_samples', 0)
    python = cfg.get('python') or shutil.which('python3.12') or shutil.which('python3')
    cfg['python'] = shutil.which(python) if python else None
    if not cfg['python']:
        raise ValueError('Python 3.12 is required; configure python with its executable path')
    version = subprocess.check_output([cfg['python'], '-c', 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")'], text=True).strip()
    if version != '3.12':
        raise ValueError(f'DIA-NN environment requires Python 3.12; selected interpreter is {version}')
    with Path(cfg['method_config']).open('rb') as handle:
        method = tomllib.load(handle)
    for name in ('precursor_q', 'global_protein_q', 'run_protein_q'):
        if not 0 < method['search'][name] <= 0.01:
            raise ValueError(f'{name} must be > 0 and <= 0.01')
    cfg['differential'] = differential_options(cfg.get('differential'))
    cfg['groups'] = audit_samples(cfg['samplesheet'], cfg['raw_root'], cfg['expected_samples'], cfg['service_matrices'])
    with Path(cfg['samplesheet']).open(newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    for group in cfg['groups']:
        validate_design_rows([row for row in rows if row['sample_type'] == group], cfg['differential'])
    if cfg['differential']['enabled']:
        selected_r = shutil.which(cfg['differential']['rscript'])
        if not selected_r:
            raise ValueError('DIA differential analysis requires Rscript with limma and jsonlite')
        subprocess.run([selected_r, '-e', 'stopifnot(requireNamespace("limma",quietly=TRUE),requireNamespace("jsonlite",quietly=TRUE))'], check=True, capture_output=True)
        cfg['differential']['rscript'] = selected_r
    cfg['expected_samples'] = sum(value['samples'] for value in cfg['groups'].values())
    return cfg


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True)
    parser.add_argument('--module-dir', required=True)
    parser.add_argument('--threads')
    parser.add_argument('--memory')
    parser.add_argument('--outdir')
    args = parser.parse_args()
    print(json.dumps(normalize_config(args.config, args.module_dir, args.threads, args.memory, args.outdir)))


if __name__ == '__main__':
    main()
