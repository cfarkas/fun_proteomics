#!/usr/bin/env python3
"""Validate arbitrary sample IDs/types and prepare independent group manifests."""

import argparse
import json
from pathlib import Path
import re
import sqlite3

import pandas as pd


def read_fasta(path):
    header, sequence = None, []
    with Path(path).open() as handle:
        for line in handle:
            line = line.strip()
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(sequence)
                header, sequence = line[1:], []
            elif line:
                sequence.append(line)
    if header is not None:
        yield header, ''.join(sequence)


def accession(header):
    token = header.split()[0]
    return token.split('|')[1] if token.count('|') >= 2 else token


def build_fasta(targets, contaminants, destination):
    records = list(read_fasta(targets))
    if not records:
        raise ValueError('Target FASTA is empty')
    ids, sequences = {accession(h) for h, _ in records}, {s for _, s in records}
    if contaminants:
        for index, (header, sequence) in enumerate(read_fasta(contaminants)):
            if accession(header) not in ids and sequence not in sequences:
                identifier = f'CONTAM_{index:05d}'
                records.append((f'sp|{identifier}|{identifier} {header}', sequence))
                sequences.add(sequence)
    with Path(destination).open('w') as handle:
        for header, sequence in records:
            if not sequence or not re.fullmatch('[A-Z]+', sequence):
                raise ValueError(f'Invalid FASTA sequence: {header}')
            handle.write(f'>{header}\n')
            for start in range(0, len(sequence), 70):
                handle.write(sequence[start:start + 70] + '\n')


def read_samplesheet(path, raw_root, expected_samples=0, check_metadata=True):
    table = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
    required = ['sample_id', 'sample_type', 'raw_path', 'library_reference']
    if table.empty or any(name not in table for name in required):
        raise ValueError(f'A nonempty sample table requires columns: {required}')
    if table[required].eq('').any().any():
        raise ValueError('Missing required sample-table value')
    if table[['sample_type', 'sample_id']].duplicated().any():
        raise ValueError('Sample IDs must be unique within each sample_type')
    if not table.sample_type.str.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*').all():
        raise ValueError('sample_type must use letters, digits, underscores or hyphens')
    if table.sample_id.str.contains(r'[\t\r\n]').any():
        raise ValueError('Sample IDs cannot contain tabs or line breaks')
    flags = table.library_reference.str.lower()
    if not flags.isin(['true', 'false']).all():
        raise ValueError('library_reference must be true or false')
    table['library_reference'] = flags.eq('true')
    paths, dates = [], []
    for value in table.raw_path:
        raw = (Path(raw_root) / value).resolve()
        if not raw.is_dir() or raw.suffix != '.d' or 'BLANK' in raw.name.upper():
            raise ValueError(f'Missing or invalid biological Bruker raw: {raw}')
        for name in ['analysis.tdf', 'analysis.tdf_bin']:
            if not (raw / name).is_file() or not (raw / name).stat().st_size:
                raise ValueError(f'Incomplete raw directory: {raw}')
        date = ''
        if check_metadata:
            with sqlite3.connect((raw / 'analysis.tdf').as_uri() + '?mode=ro', uri=True) as connection:
                meta = dict(connection.execute('SELECT Key, Value FROM GlobalMetadata'))
                frame_types = {row[0] for row in connection.execute('SELECT DISTINCT MsMsType FROM Frames')}
            if frame_types != {0, 9}:
                raise ValueError(f'Expected only MS1 and diaPASEF frames in {raw}; found {sorted(frame_types)}')
            date = meta.get('AcquisitionDateTime', '')
        paths.append(str(raw))
        dates.append(date)
    table['raw_path'] = paths
    table['run'] = [Path(path).stem for path in paths]
    table['acquisition_datetime'] = dates
    if not table.raw_path.is_unique or table[['sample_type', 'run']].duplicated().any():
        raise ValueError('Duplicate raw path or duplicate run basename within a sample type')
    if expected_samples and len(table) != expected_samples:
        raise ValueError(f'Expected {expected_samples} samples, received {len(table)}')
    for name, group in table.groupby('sample_type'):
        if not group.library_reference.any():
            raise ValueError(f'At least one library_reference=true is required for {name}')
    return table.sort_values(['sample_type', 'acquisition_datetime', 'run'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--samplesheet', required=True)
    parser.add_argument('--raw-root', required=True)
    parser.add_argument('--canonical-fasta', required=True)
    parser.add_argument('--contaminant-fasta')
    parser.add_argument('--expected-samples', type=int, default=0)
    parser.add_argument('--service-types', default='[]')
    args = parser.parse_args()
    samples = read_samplesheet(args.samplesheet, args.raw_root, args.expected_samples)
    unknown = set(json.loads(args.service_types)) - set(samples.sample_type)
    if unknown:
        raise ValueError(f'Service matrices configured for absent sample types: {sorted(unknown)}')
    output = Path('prepared')
    (output / 'groups').mkdir(parents=True, exist_ok=True)
    samples.to_csv(output / 'samples.tsv', sep='\t', index=False)
    for sample_type, group in samples.groupby('sample_type', sort=True):
        group.to_csv(output / 'groups' / f'{sample_type}.tsv', sep='\t', index=False)
        print(f'{sample_type}: {len(group)} samples; {int(group.library_reference.sum())} library references')
    build_fasta(args.canonical_fasta, args.contaminant_fasta, output / 'search.fasta')


if __name__ == '__main__':
    main()
