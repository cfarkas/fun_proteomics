#!/usr/bin/env python3
"""Check native main.nf DIA routing and resume using synthetic acquisitions."""
import csv
import json
import os
from pathlib import Path
import subprocess
from test_diapasef_config import make_fixture

ROOT = Path(__file__).resolve().parents[1]
fixture = ROOT / '.test-runtime/dia'
config = make_fixture(fixture)
# Keep a private lock copy so the final check changes a real scientific input.
lock = fixture / 'requirements-lock.txt'
lock.write_bytes((ROOT / 'modules/diapasef/requirements-lock.txt').read_bytes())
settings = json.loads(config.read_text())
settings['requirements_lock'] = str(lock)
config.write_text(json.dumps(settings))
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
# A shared download cache may gain unrelated environments between invocations.
unrelated = fixture / 'cache/unrelated R library'
unrelated.mkdir(exist_ok=True)
(unrelated / 'installed.txt').write_text('Not an analysis input\n')
subprocess.run([*command, '-resume', '-with-trace', 'cache_growth_trace.tsv'], cwd=fixture, check=True)
with (fixture / 'cache_growth_trace.tsv').open() as handle:
    cache_growth = list(csv.DictReader(handle, delimiter='\t'))
assert len(cache_growth) == 10 and all(row['status'] == 'CACHED' for row in cache_growth), cache_growth
# Actual dependency-lock changes must still invalidate the environment task.
lock.write_text(lock.read_text() + '\n# Changed lock input for cache invalidation test\n')
subprocess.run([*command, '-resume', '-with-trace', 'changed_lock_trace.tsv'], cwd=fixture, check=True)
with (fixture / 'changed_lock_trace.tsv').open() as handle:
    changed = list(csv.DictReader(handle, delimiter='\t'))
assert any(row['name'].endswith(':DIA_ENVIRONMENT') and row['status'] == 'COMPLETED' for row in changed), changed
print('PASS: ten native DIA tasks resume after unrelated installation-cache growth; dependency-lock changes invalidate the environment.')
