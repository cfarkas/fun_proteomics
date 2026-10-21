#!/usr/bin/env python3
"""Quantify all requested runs of one sample type against its empirical library."""

import argparse
from pathlib import Path
from common import read_method, read_samples, run_engine, search_arguments, validate_report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['environment', 'method', 'fasta', 'samples', 'library']:
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--threads', type=int, required=True)
    args = parser.parse_args()
    method, samples = read_method(args.method), read_samples(args.samples)
    if samples.sample_type.nunique() != 1:
        raise ValueError('Each quantification task processes one sample type')
    output = Path('cohort').resolve()
    command = search_arguments(method, args.fasta, Path(args.library) / 'empirical.parquet', samples, output)
    run_engine(args.environment, command, output, args.threads, method, ['report.parquet', 'report.pg_matrix.tsv'])
    validate_report(output / 'report.parquet', samples)


if __name__ == '__main__':
    main()
