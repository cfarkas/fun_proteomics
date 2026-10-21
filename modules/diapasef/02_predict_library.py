#!/usr/bin/env python3
"""Predict a spectral library from the supplied target/contaminant FASTA."""

import argparse
from pathlib import Path
from common import modifications, read_method, run_engine


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--environment', required=True)
    parser.add_argument('--method', required=True)
    parser.add_argument('--fasta', required=True)
    parser.add_argument('--threads', type=int, required=True)
    args = parser.parse_args()
    method = read_method(args.method)
    search = method['search']
    output = Path('predicted').resolve()
    command = ['--fasta', Path(args.fasta).resolve(), '--fasta-search', '--predictor',
               '--out-lib', output / 'database.parquet', '--cut', 'K*,R*,!*P',
               '--missed-cleavages', method['method']['missed_cleavages'],
               '--min-pep-len', search['peptide_length_min'], '--max-pep-len', search['peptide_length_max'],
               '--min-pr-charge', search['charge_min'], '--max-pr-charge', search['charge_max'],
               '--min-pr-mz', search['precursor_mz_min'], '--max-pr-mz', search['precursor_mz_max'],
               *modifications(method)]
    run_engine(args.environment, command, output, args.threads, method, ['database.predicted.speclib'])


if __name__ == '__main__':
    main()
