#!/usr/bin/env python3
"""Check native main.nf DIA routing and resume using synthetic acquisitions."""
import csv
import os
from pathlib import Path
import subprocess
from test_diapasef_config import make_fixture

ROOT = Path(__file__).resolve().parents[1]
fixture = ROOT / '.test-runtime/dia'
config = make_fixture(fixture)
nextflow = os.environ.get('NEXTFLOW_BINARY', 'nextflow')
command = [nextflow, '-log', str(fixture / 'nextflow.log'), 'run', str(ROOT / 'main.nf'),
           '--dia_config', str(config), '-profile', 'test', '-stub-run', '-ansi-log', 'false',
           '-work-dir', str(fixture / 'work')]
subprocess.run([*command, '-with-trace', 'trace.tsv'], cwd=fixture, check=True)
with (fixture / 'trace.tsv').open() as handle:
    first = list(csv.DictReader(handle, delimiter='\t'))
assert len(first) == 10, first
assert all(row['status'] == 'COMPLETED' for row in first), first
assert all('DIA_PASEF:DIA_' in row['name'] for row in first), first
for group in ('serum', 'cell_pellet'):
    for stage in ('libraries', 'search', 'groups'):
        assert (fixture / 'results' / 'analysis_details' / stage / group).is_dir()
assert (fixture / 'results/final_results/STUB.txt').is_file()
assert not (fixture / 'results/final_results/sample_counts.tsv').exists()
subprocess.run([*command, '-resume', '-with-trace', 'resume_trace.tsv'], cwd=fixture, check=True)
with (fixture / 'resume_trace.tsv').open() as handle:
    resumed = list(csv.DictReader(handle, delimiter='\t'))
assert len(resumed) == 10 and all(row['status'] == 'CACHED' for row in resumed), resumed
print('PASS: main.nf DIA branch routes two types through ten native tasks; all ten resume from cache.')
