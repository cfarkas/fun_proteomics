#!/usr/bin/env bash
set -euo pipefail
repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
fixture="$repo/.test-runtime"
python "$repo/tests/make_fixture.py" "$fixture"
cd "$fixture"
nextflow run "$repo/main.nf" --config "$fixture/test.local.json" -stub-run -profile test -ansi-log false -with-trace trace.tsv
python - "$fixture/results/analysis_details/final" <<'PY'
import csv,json,sys
from pathlib import Path
out=Path(sys.argv[1])
with (out/'sample_detection_counts.tsv').open() as h:
    rows=list(csv.DictReader(h,delimiter='\t'))
assert len(rows)==6, rows
assert all(r['maxlfq_protein_groups']=='' for r in rows if r['cohort']=='serum')
with (out/'dda_vs_ddaplus.tsv').open() as h:
    comparisons=list(csv.DictReader(h,delimiter='\t'))
assert len(comparisons)==3
assert all(r['direct_protein_groups_delta']=='1' for r in comparisons)
assert not (out/'provenance.json').exists()
for cohort, mode in [('cell', 'dda'), ('cell', 'ddaplus'), ('serum', 'dda'), ('serum', 'ddaplus')]:
    folder = out.parent / cohort / mode
    assert (folder / 'run/run_provenance.json').is_file()
    assert not (folder / f'tables_{cohort}_{mode}/report_provenance.json').exists()
assert (out/'STUB_COMPLETE_NOT_SCIENTIFIC_RESULTS').is_file()
assert (out/'updated_count_tables.xlsx').stat().st_size>1000
print('PASS: two modes, separate cell/serum cohorts, six sample reports, paired deltas and Excel export')
PY
nextflow run "$repo/main.nf" --config_file "$fixture/test.local.json" -stub-run -profile test -ansi-log false -resume -with-trace resume_trace.tsv
python - "$fixture/resume_trace.tsv" <<'PY'
import csv,sys
with open(sys.argv[1]) as h:
    rows=list(csv.DictReader(h,delimiter='\t'))
assert len(rows)==16, rows
assert all(r['status']=='CACHED' for r in rows if r['name']!='VALIDATE_INPUTS' and not r['name'].startswith('CHECK_SEARCH')), rows
assert next(r for r in rows if r['name']=='VALIDATE_INPUTS')['status']=='COMPLETED'
print('PASS: resume reruns validation/search verification and reuses all 11 eligible tasks')
PY
# Supplying both aliases must still fail before any analysis starts.
if nextflow run "$repo/main.nf" --config "$fixture/test.local.json" --config_file "$fixture/test.local.json" -ansi-log false > "$fixture/alias_conflict.log" 2>&1; then
    echo 'FAIL: conflicting configuration aliases were accepted' >&2
    exit 1
fi
python - "$fixture/alias_conflict.log" <<'CHECK_ALIAS'
from pathlib import Path
import sys
assert 'Use --config OR --config_file, not both' in Path(sys.argv[1]).read_text()
CHECK_ALIAS
echo 'PASS: both config aliases work independently and conflict only when both are supplied'
