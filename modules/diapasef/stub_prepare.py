"""Write routing fixtures only; deliberately produce no scientific measurements."""
import csv
from pathlib import Path
import re

output = Path('prepared')
(output / 'groups').mkdir(parents=True, exist_ok=True)
with Path('samplesheet.tsv').open() as handle:
    reader = csv.DictReader(handle, delimiter='\t')
    columns, rows = reader.fieldnames, list(reader)
for group in sorted({row['sample_type'] for row in rows}):
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*', group):
        raise ValueError('Invalid stub sample type')
    with (output / 'groups' / f'{group}.tsv').open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t')
        writer.writeheader()
        writer.writerows(row for row in rows if row['sample_type'] == group)
(output / 'samples.tsv').write_text(Path('samplesheet.tsv').read_text())
(output / 'search.fasta').write_text('STUB ONLY - not a real FASTA\n')
