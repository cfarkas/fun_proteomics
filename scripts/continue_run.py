"""Select the latest real Nextflow session for this study, including completed runs."""
import argparse
import json
from pathlib import Path
import shlex
import uuid


def option(argv, names, default=None):
    for i, arg in enumerate(argv):
        if arg in names and i+1 < len(argv):
            return argv[i+1]
        if '=' in arg and arg.split('=', 1)[0] in names:
            return arg.split('=', 1)[1]
    return default


def absolute(value, base):
    path = Path(value).expanduser()
    return (path if path.is_absolute() else Path(base)/path).resolve()


def select_session(project, launch, config, outdir=None):
    project, launch = Path(project).resolve(), Path(launch).resolve()
    config = absolute(config, launch)
    empty = {'session': None, 'workdir': None}
    try:
        study = json.loads(config.read_text())
        output = absolute(outdir, launch) if outdir else absolute(study['outdir'], config.parent)
        lines = (launch/'.nextflow/history').read_text().splitlines()
    except (OSError, ValueError, KeyError):
        return empty
    excluded = {'--help', '--interactively', '--validate_only', '--dry_run', '--dia_config', '-preview', '-stub', '-stub-run'}
    for line in reversed(lines):
        fields = line.split('\t')
        if len(fields) < 7:
            continue
        try:
            argv = shlex.split(fields[6])
            if any(arg.split('=')[0] in excluded for arg in argv):
                continue
            run = argv.index('run')
            if absolute(argv[run+1], launch) != project/'main.nf':
                continue
            prior = absolute(option(argv, ('--config', '--config_file'), project/'config.local.json'), launch)
            previous = json.loads(prior.read_text())
            override = option(argv, ('--outdir',))
            previous_output = absolute(override, launch) if override else absolute(previous['outdir'], prior.parent)
            if previous_output != output:
                continue
            session = str(uuid.UUID(fields[5]))
            if not (launch/'.nextflow/cache'/session).is_dir():
                continue
            work = option(argv, ('-w', '-work-dir'))
            return {'session': session, 'workdir': str(absolute(work, launch)) if work else None,
                    'name': fields[2]}
        except (OSError, ValueError, KeyError, IndexError):
            continue
    return empty


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project', required=True)
    parser.add_argument('--launch', required=True)
    parser.add_argument('--config', required=True)
    parser.add_argument('--outdir')
    args = parser.parse_args()
    print(json.dumps(select_session(args.project, args.launch, args.config, args.outdir)))
