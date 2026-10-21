#!/usr/bin/env python3
"""Build an empirical library only from this sample type's declared references."""

import argparse
from pathlib import Path
import pandas as pd
from common import read_method, read_samples, run_engine, search_arguments, validate_report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['environment', 'method', 'fasta', 'samples', 'prediction']:
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--threads', type=int, required=True)
    args = parser.parse_args()
    method = read_method(args.method)
    samples = read_samples(args.samples)
    selected = samples.loc[samples.library_reference].copy()
    if selected.empty or samples.sample_type.nunique() != 1:
        raise ValueError('A library requires references from exactly one sample type')
    output = Path('library').resolve()
    predicted = Path(args.prediction) / 'database.predicted.speclib'
    command = search_arguments(method, args.fasta, predicted, selected, output, build_library=True)
    run_engine(args.environment, command, output, args.threads, method,
               ['report.parquet', 'report.pg_matrix.tsv', 'empirical.parquet'])
    validate_report(output / 'report.parquet', selected)
    entries = pd.read_parquet(output / 'empirical.parquet', columns=['Precursor.Id', 'Decoy'])
    if not entries.Decoy.eq(0).any() or not entries.Decoy.eq(1).any():
        raise ValueError('Empirical library must preserve targets and decoys')
    print(f'Library complete for {samples.sample_type.iloc[0]}: {len(selected)} reference runs')


if __name__ == '__main__':
    main()
