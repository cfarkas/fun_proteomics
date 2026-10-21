"""Shared scientific operations for the numbered Nextflow tasks."""

import os
from pathlib import Path
import shutil
import subprocess
import tomllib

import numpy as np
import pandas as pd
import pyarrow.parquet as pq


def read_method(path):
    with Path(path).open('rb') as handle:
        method = tomllib.load(handle)
    for key in ['precursor_q', 'run_protein_q', 'global_protein_q']:
        if not 0 < method['search'][key] <= .01:
            raise ValueError(f'{key} must be in (0, 0.01]')
    return method


def read_samples(path):
    table = pd.read_csv(path, sep='\t', dtype={'sample_id': str, 'sample_type': str}, keep_default_na=False)
    table['library_reference'] = table.library_reference.astype(str).str.lower().eq('true')
    return table


def modifications(method):
    settings = method['method']
    args = ['--var-mods', str(settings['max_variable_modifications'])]
    if method['search']['fixed_carbamidomethyl_c']:
        args += ['--fixed-mod', 'UniMod:4,57.021464,C']
    if settings['methionine_oxidation']:
        args += ['--var-mod', 'UniMod:35,15.994915,M']
    if settings['protein_nterm_acetylation']:
        args += ['--var-mod', 'UniMod:1,42.010565,*n']
    return args


def run_engine(environment, arguments, output, threads, method, required):
    """Run one foreground engine; Nextflow owns task completion and resume."""
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    if shutil.disk_usage(output).free / 1e9 < method['compute']['minimum_free_gb']:
        raise RuntimeError('Insufficient free space for the configured minimum_free_gb')
    binaries = list((Path(environment).resolve() / 'software/diann-2.7.0').rglob('diann-linux'))
    if len(binaries) != 1:
        raise RuntimeError('The environment must contain exactly one DIA-NN executable')
    binary = binaries[0]
    command = [str(binary), '--parent', str(binary.parent), '--threads', str(threads), *map(str, arguments)]
    env = {**os.environ, 'OMP_NUM_THREADS': str(threads), 'OPENBLAS_NUM_THREADS': '1',
           'LD_LIBRARY_PATH': str(binary.parent) + ':' + os.environ.get('LD_LIBRARY_PATH', '')}
    with (output / 'engine.log').open('a') as log:
        subprocess.run(command, cwd=output, env=env, stdin=subprocess.DEVNULL,
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    for name in required:
        path = output / name
        if not path.is_file() or not path.stat().st_size:
            raise RuntimeError(f'Missing engine output: {path}')


def search_arguments(method, fasta, library, samples, output, build_library=False):
    output = Path(output).resolve()
    (output / 'cache').mkdir(parents=True, exist_ok=True)
    search = method['search']
    args = ['--lib', Path(library).resolve(), '--fasta', Path(fasta).resolve()]
    for raw in samples.raw_path:
        args += ['--f', raw]
    args += ['--out', output / 'report.parquet', '--matrices',
             '--qvalue', search['precursor_q'], '--matrix-qvalue', search['global_protein_q'],
             '--matrix-spec-q', search['run_protein_q'],
             '--mass-acc', search['ms2_ppm'], '--mass-acc-ms1', search['ms1_ppm'],
             '--window', method['method']['scan_window'],
             '--read-threads', method['compute']['read_threads'],
             '--temp', output / 'cache', '--use-quant', '--export-quant', *modifications(method)]
    if build_library or method['method']['force_ms2']:
        args += ['--force-ms2']
    if build_library:
        args += ['--out-lib', output / 'empirical.parquet', '--gen-spec-lib',
                 '--out-lib-qvalue', search['precursor_q']]
    else:
        args += ['--no-prot-inf']
    return args


def validate_report(path, samples):
    observed = set(pd.read_parquet(path, columns=['Run']).Run)
    if observed != set(samples.run):
        raise ValueError('Engine report does not contain exactly the requested runs')


def canonical_group(value):
    return ';'.join(sorted(set(str(value).split(';'))))


def filter_report(report, method):
    """Apply run, global and empirical-library confidence without relaxing cutoffs."""
    search = method['search']
    limits = {'Q.Value': search['precursor_q'], 'Global.Q.Value': search['precursor_q'],
              'Lib.Q.Value': search['precursor_q'], 'PG.Q.Value': search['run_protein_q'],
              'Global.PG.Q.Value': search['global_protein_q'], 'Lib.PG.Q.Value': search['global_protein_q']}
    if method['method']['methionine_oxidation']:
        limits['Lib.Peptidoform.Q.Value'] = search['precursor_q']
    missing = set(limits) - set(report.columns)
    if missing:
        raise ValueError(f'Missing required confidence fields: {sorted(missing)}')
    keep = (np.isfinite(report['PG.MaxLFQ']) & report['PG.MaxLFQ'].gt(0)
            & report['Protein.Group'].notna()
            & ~report['Protein.Group'].str.contains(r'(?:^|[;|])(?:contam|rev|reverse|decoy|xxx)[_|-]', case=False, na=True))
    for column, cutoff in limits.items():
        keep &= np.isfinite(report[column]) & report[column].between(0, cutoff)
    if 'Decoy' in report:
        keep &= report.Decoy.eq(0)
    return report.loc[keep].copy()


def load_report(path):
    schema = pq.ParquetFile(path).schema_arrow.names
    columns = ['Run', 'Protein.Group', 'Protein.Names', 'Genes', 'Precursor.Id',
               'Stripped.Sequence', 'Proteotypic', 'Q.Value', 'PG.Q.Value', 'Global.Q.Value',
               'Global.PG.Q.Value', 'Lib.Q.Value', 'Lib.PG.Q.Value', 'Lib.Peptidoform.Q.Value',
               'PG.MaxLFQ', 'Precursor.Normalised', 'Precursor.Quantity', 'Decoy', 'RT', 'IM']
    columns += [name for name in schema if name.startswith('Fr.') and name.endswith('.Quantity')]
    return pd.read_parquet(path, columns=[name for name in columns if name in schema])
