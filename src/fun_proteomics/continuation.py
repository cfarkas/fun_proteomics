"""Discover run folders and verify completed artifacts before choosing reuse."""
from datetime import datetime, timezone
import csv
import json
import re
from pathlib import Path
import shlex
from .project import digest


def now():
    return datetime.now(timezone.utc).isoformat()


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name+'.tmp')
    temporary.write_text(json.dumps(value, indent=2)+'\n')
    temporary.replace(path)


def run_root(folder):
    root = Path(folder).expanduser().resolve()
    if not root.is_dir():
        raise ValueError('Run folder does not exist: '+str(root))
    return root.parent if root.name == 'final_results' else root


def verify_delivery(root):
    """Return invalid/missing relative files; require a complete hash inventory."""
    final = Path(root)/'final_results'
    manifest = final/'OUTPUT_SHA256SUMS'
    if not manifest.is_file() or not (final/'FINAL_RESULTS_COMPLETE').is_file():
        return ['completion marker or output checksum manifest missing']
    failures = []
    entries = manifest.read_text().splitlines()
    if not entries:
        return ['output checksum manifest is empty']
    seen = set()
    for line in entries:
        expected, relative = line.split('  ', 1)
        path = final/relative
        if not re.fullmatch(r'[0-9a-f]{64}', expected) or Path(relative).is_absolute() or '..' in Path(relative).parts or not path.resolve().is_relative_to(final.resolve()):
            raise ValueError('Invalid delivery checksum entry: '+relative)
        if relative in seen:
            raise ValueError('Duplicate delivery checksum entry: '+relative)
        seen.add(relative)
        if not path.is_file() or digest(path) != expected:
            failures.append(relative)
    actual = {p.relative_to(final).as_posix() for p in final.rglob('*') if p.is_file() and p.name != 'OUTPUT_SHA256SUMS'}
    failures.extend(sorted(actual-seen))
    return failures


def verify_recorded_inputs(root, hash_raw=False):
    """A changed acquisition invalidates reuse of the old search report."""
    root = Path(root)
    record_path = root/'audit/run.json'
    count_records = [root/'final_results/input_provenance.json', root/'analysis_details/reporting/result_audit/main_statistics/input_provenance.json']
    for count_record in count_records:
        for item in (read_json(count_record) if count_record.is_file() else []):
            path = Path(item['path'])
            if not path.is_file() or digest(path) != item['sha256']:
                raise ValueError('Recorded count analysis input changed or disappeared: '+str(path))
    if not record_path.is_file():
        return
    record = read_json(record_path)
    for item in record.get('input_provenance', []):
        path = Path(item['path'])
        if not path.is_file() or digest(path) != item['sha256']:
            raise ValueError('Recorded analysis input changed or disappeared: '+str(path))
    for entry in record.get('searches', {}).values():
        for key in ('report', 'cached_report'):
            item = entry[key]
            if not Path(item['path']).is_file() or digest(item['path']) != item['sha256']:
                raise ValueError('Cached search report integrity failed: '+item['path'])
    raw_path = root/'audit/raw_integrity.json'
    if raw_path.is_file():
        for entry in read_json(raw_path):
            for item in entry['files']:
                path = Path(item['path'])
                if not path.is_file():
                    raise ValueError('Raw acquisition component missing: '+str(path))
                stat = path.stat()
                changed = stat.st_size != item['bytes'] or stat.st_mtime_ns != item['mtime_ns']
                if changed or hash_raw:
                    if not item.get('sha256'):
                        raise ValueError('No recorded raw SHA256 baseline; cannot verify acquisition identity: '+str(path))
                    if digest(path) != item['sha256']:
                        raise ValueError('Raw acquisition changed; regenerate its search before continuing: '+str(path))


def option_values(argv, flag):
    values = []
    for i, value in enumerate(argv):
        if value == flag and i+1 < len(argv):
            values.append(argv[i+1])
        elif value.startswith(flag+'='):
            values.append(value.split('=',1)[1])
    return values


def dia_settings(root):
    root = Path(root)
    audit = root/'audit/run.json'
    if audit.is_file():
        record = read_json(audit)
        if record.get('source_run') and record.get('searches'):
            source = Path(record['source_run']).resolve()
            if not source.is_dir():
                raise ValueError('Original reusable search folder is unavailable: '+str(source))
            argv = record.get('command', [])
            samples = root/'inputs/samples.tsv'
            raw_audit = root/'audit/raw_integrity.json'
            hashed_raw = raw_audit.is_file() and any(item.get('sha256') for entry in read_json(raw_audit) for item in entry.get('files', []))
            return dict(source_run=str(source), samples=str(samples) if samples.is_file() else None,
                        independent=record.get('independent_samples_confirmed', False),
                        services=option_values(argv, '--service-matrix'),
                        method=str(root/'inputs/method.toml'),
                        r_library=str(Path(record.get('repository', ''))/'.runtime/R-library'),
                        differential=record.get('differential', {}), hash_raw=hashed_raw)
    reports = list((root/'analysis_details/search').glob('*/cohort/report.parquet')) + list((root/'search').glob('*/cohort/report.parquet'))
    if reports and (root/'prepared/samples.tsv').is_file() and (root/'prepared/search.fasta').is_file():
        with (root/'prepared/samples.tsv').open(newline='') as handle:
            groups = {row['sample_type'] for row in csv.DictReader(handle, delimiter='\t')}
        if len(reports) != len(groups) or {p.parents[1].name for p in reports} != groups:
            return None
        return dict(source_run=str(root), samples=None, independent=False, services=[], method=None,
                    r_library=None, differential={})
    return None


def _study_identity(argv, launch):
    """Match study paths, including quoted paths, independently of executable names."""
    launch = Path(launch).resolve()
    result = {}
    for flag in ('--config', '--config_file', '--dia_config', '--outdir'):
        values = option_values(argv, flag)
        if values:
            path = Path(values[-1]).expanduser()
            result[flag] = str((path if path.is_absolute() else launch/path).resolve())
    if 'run' in argv and argv.index('run')+1 < len(argv):
        path = Path(argv[argv.index('run')+1])
        result['workflow'] = str((path if path.is_absolute() else launch/path).resolve())
    return result


def nextflow_session(launch, command=None):
    history = Path(launch)/'.nextflow/history'
    if not history.is_file():
        return None
    expected = _study_identity(command, launch) if command else None
    for line in reversed(history.read_text().splitlines()):
        fields = line.split('\t')
        if len(fields) >= 7 and (Path(launch)/'.nextflow/cache'/fields[5]).is_dir():
            argv = shlex.split(fields[6])
            if any(flag in argv for flag in ('--help','--interactively','--validate_only','-stub-run','-preview')):
                continue
            if expected is not None and _study_identity(argv, launch) != expected:
                continue
            return fields[5]
    return None


def latest_legacy_nextflow(root):
    records = sorted((Path(root)/'logs').glob('run_*.json'), key=lambda p:p.stat().st_mtime, reverse=True)
    for path in records:
        value = read_json(path)
        if value.get('kind') == 'analysis' and value.get('command') and value.get('root'):
            return value
    return None


def fresh_output(root):
    return Path(root)/'continuations'/datetime.now(timezone.utc).strftime('%Y-%m-%d_%H-%M-%S_%f')
