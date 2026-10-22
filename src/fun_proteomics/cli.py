"""A public Python entry point; raw workflows remain native Nextflow runs."""
import argparse
import fcntl
import importlib.util
import json
import os
from pathlib import Path
import shlex
import shutil
import signal
import subprocess
import sys
from . import __version__
from .project import digest, resolve_project
from .continuation import (dia_settings, fresh_output, latest_legacy_nextflow, nextflow_session,
                           now, option_values, read_json, run_root, verify_delivery, verify_recorded_inputs, write_json)


def parser():
    p = argparse.ArgumentParser(prog='fun_proteomics', description='Set up a study, check it, run it, and open clearly organized results.', epilog='Start: fun_proteomics setup --project /data/my-study. SSH: add --terminal. Resume: --continue /data/my-study')
    p.add_argument('--interactive', action='store_true', help='Shortcut for setup; follow with --project STUDY')
    p.add_argument('--version', action='version', version='fun_proteomics '+__version__)
    p.add_argument('--workflow-dir', '--project', dest='project', help='Reuse a workflow checkout and its runtime (--project is the legacy alias)')
    p.add_argument('--cache-dir', help='Writable location for the packaged workflow')
    p.add_argument('--rscript', help='Rscript executable for native statistics')
    p.add_argument('--r-library', help='Existing R package library')
    p.add_argument('--nextflow', help='Nextflow executable')
    p.add_argument('--print-command', action='store_true', help='Inspect the launch plan without starting work or changing results')
    p.add_argument('--continue', dest='continue_folder', metavar='FOLDER', help='Detect and continue an existing run folder')
    p.add_argument('--outdir', help='Explicit new output folder for native continuation')
    p.add_argument('--service-matrix', action='append', default=[], metavar='COHORT=PATH')
    p.add_argument('--independent-samples', action='store_true', help='Confirm one independent biological sample per row')
    p.add_argument('--hash-raw', action='store_true', help='Hash all raw components during verification')
    sub = p.add_subparsers(dest='command')
    from .study import add_commands
    add_commands(sub)
    run = sub.add_parser('run', help='Run or continue a saved study')
    group = run.add_mutually_exclusive_group(required=True)
    group.add_argument('--config', help='Saved PROJECT/config/run.json (also accepts legacy DDA JSON)')
    group.add_argument('--dia-config', help='DIA study JSON')
    run.add_argument('--outdir', dest='run_outdir', help='Required only for legacy raw-workflow JSON')
    run.add_argument('--validate-only', action='store_true')
    ui = sub.add_parser('ui', help='Open guided setup; --advanced opens the original raw editor')
    ui.add_argument('--config')
    ui.add_argument('--port', type=int, default=0)
    ui.add_argument('--advanced', action='store_true')
    ui.add_argument('--no-browser', action='store_true')
    doctor = sub.add_parser('doctor', help='Check tools and explain what is missing')
    doctor.add_argument('--json', action='store_true')
    doctor.add_argument('--mode', choices=['matrix', 'dia', 'dda', 'baf'], default='matrix')
    install = sub.add_parser('install', help='Install analysis tools in the writable workflow cache')
    install.add_argument('--baf', action='store_true', help='Install the separate native BAF reader and mzML writer')
    install.add_argument('--conda', action='store_true', help='Install pinned Python/Java/Nextflow and R tools')
    install.add_argument('--r-only', action='store_true', help='Install analysis packages using an existing Rscript')
    sub.add_parser('workflow-path', help='Print the workflow directory for direct nextflow run')
    for command, help_text in [('quantify','Import a quantification matrix using native Python'),
                               ('counts','Reanalyze delivered DDA counts using native Python'),
                               ('matrices','Export CSV/TSV matrix companions and index'),
                               ('raw-baf','Search native BAF acquisitions and quantify MS1 LFQ'),
                               ('benchmark','Benchmark any matched pipeline and service matrices')]:
        child = sub.add_parser(command, help=help_text, add_help=False)
        child.add_argument('arguments', nargs=argparse.REMAINDER)
    return p


def environment(project, args, inherited_r_library=None, inherited_rscript=None):
    env = os.environ.copy()
    rscript = args.rscript or inherited_rscript or env.get('PROTEOMICS_RSCRIPT')
    if not rscript:
        settings=project/'.runtime/runtime-settings.json'
        installed=read_json(settings).get('rscript') if settings.is_file() else None
        local=project/'.runtime/r-env/bin/Rscript'
        rscript=installed if installed and Path(installed).is_file() else str(local) if local.is_file() else None
    if not rscript:
        rscript = '/usr/bin/Rscript' if Path('/usr/bin/Rscript').is_file() else shutil.which('Rscript')
    if rscript:
        env['PROTEOMICS_RSCRIPT'] = rscript
    library = args.r_library or inherited_r_library
    if not library and (project/'.runtime/R-library').is_dir():
        library = str(project/'.runtime/R-library')
    if library:
        if not Path(library).is_dir():
            if args.r_library:
                raise ValueError('R library does not exist: '+library)
        else:
            env['R_LIBS'] = os.pathsep.join(filter(None,[str(Path(library).resolve()),env.get('R_LIBS')]))
    java=project/'.runtime/env'
    if (java/'bin/java').is_file():
        env['JAVA_HOME']=str(java)
        env['PATH']=str(java/'bin')+os.pathsep+env.get('PATH','')
    env.setdefault('OPENBLAS_NUM_THREADS','1')
    env.setdefault('OMP_NUM_THREADS','1')
    return env


def nextflow_binary(project, args, env):
    path = args.nextflow or str(project/'.runtime/env/bin/nextflow')
    binary = shutil.which(path) or (shutil.which('nextflow') if not args.nextflow else None)
    if not binary:
        raise ValueError('Nextflow is not installed. Install Nextflow/Java or use --nextflow /path/to/nextflow. See doctor.')
    java = project/'.runtime/env'
    if (java/'bin/java').is_file():
        env['JAVA_HOME'] = str(java)
        env['PATH'] = str(java/'bin')+os.pathsep+env.get('PATH','')
    env.setdefault('NXF_VER','26.04.2')
    env.setdefault('NXF_DISABLE_CHECK_LATEST','true')
    return str(Path(binary).absolute())


def process_identity(pid):
    try:
        stat = Path('/proc')/str(pid)/'stat'
        return Path('/proc/sys/kernel/random/boot_id').read_text().strip()+':'+stat.read_text().split(') ',1)[1].split()[19]
    except (OSError, IndexError):
        return None


def absolute_options(arguments, flags):
    arguments = list(arguments)
    for flag in flags:
        for index, token in enumerate(arguments):
            if token == flag and index+1 < len(arguments):
                arguments[index+1] = str(Path(arguments[index+1]).expanduser().resolve())
            elif token.startswith(flag+'='):
                arguments[index] = flag+'='+str(Path(token.split('=',1)[1]).expanduser().resolve())
    return arguments


def replace_option(arguments, flag, value):
    arguments = list(arguments)
    for index, token in enumerate(arguments):
        if token == flag:
            arguments[index+1] = str(value)
            return arguments
        if token.startswith(flag+'='):
            arguments[index] = flag+'='+str(value)
            return arguments
    raise ValueError('Required argument missing: '+flag)


def recorded_files(paths):
    files = set()
    for value in paths:
        path = Path(value)
        if path.is_dir():
            files.update(p.resolve() for p in path.rglob('*') if p.is_file())
        elif path.is_file():
            files.add(path.resolve())
        else:
            raise ValueError('Input does not exist: '+str(path))
    return [dict(path=str(path), sha256=digest(path)) for path in sorted(files)]


def execute(command, cwd, env, state_path, record, inspect=False):
    print(shlex.join(command), flush=True)
    if inspect:
        return 0
    state_path = Path(state_path)
    state_path.parent.mkdir(parents=True, exist_ok=True)
    Path(cwd).mkdir(parents=True, exist_ok=True)
    with (state_path.parent/'run.lock').open('a') as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise ValueError('Another launcher is already using this run folder') from error
        record = dict(record, command=command, launch_dir=str(Path(cwd).resolve()), status='running', started=now(), package_version=__version__)
        write_json(state_path, record)
        process = None
        try:
            process = subprocess.Popen(command, cwd=cwd, env=env, start_new_session=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            record.update(pid=process.pid, process_identity=process_identity(process.pid))
            write_json(state_path, record)
            with (state_path.parent/'analysis.log').open('a') as log:
                for line in process.stdout:
                    print(line,end='',flush=True);log.write(line);log.flush()
            code = process.wait()
        except KeyboardInterrupt:
            if process is not None and process.poll() is None:
                os.killpg(process.pid, signal.SIGINT)
                try:
                    process.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGTERM)
                    try:
                        process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                        process.wait()
            code = 130
        except Exception as error:
            record.update(status='failed', ended=now(), error=str(error))
            write_json(state_path, record)
            raise
        finally:
            if process is not None and process.stdout is not None: process.stdout.close()
        record.update(status='completed' if code == 0 else 'failed', exit_code=code, ended=now())
        if record['engine'] == 'nextflow':
            record['session'] = nextflow_session(cwd, command)
        final = Path(record['outdir'])/'final_results'
        if code == 0 and record['engine'] in ('native-quantification','native-benchmark') and (final/'START_HERE.html').is_file():
            (final/'FINAL_RESULTS_COMPLETE').write_text('Native quantification and reporting completed.\n')
        if code == 0 and (final/'FINAL_RESULTS_COMPLETE').is_file():
            final = Path(record['outdir'])/'final_results'
            paths = sorted(p for p in final.rglob('*') if p.is_file() and p.name != 'OUTPUT_SHA256SUMS')
            (final/'OUTPUT_SHA256SUMS').write_text(''.join(digest(p)+'  '+p.relative_to(final).as_posix()+'\n' for p in paths))
        write_json(state_path, record)
        return code


def continue_run(args, project):
    folder = run_root(args.continue_folder)
    state_path = folder/'.fun_proteomics/run.json'
    saved = read_json(state_path) if state_path.is_file() else None
    if saved and saved.get('status') == 'running' and saved.get('process_identity') and process_identity(saved.get('pid')) == saved['process_identity']:
        raise ValueError('This run is still active (PID '+str(saved['pid'])+'); no second run was started')
    active_output = Path(saved['outdir']) if saved else folder
    for item in (saved or {}).get('input_provenance', []):
        if not Path(item['path']).is_file() or digest(item['path']) != item['sha256']:
            raise ValueError('Analysis input changed or disappeared: '+item['path']+'; start a new run with the updated inputs')
    if saved and saved.get('engine')=='native-baf':
        if args.outdir or args.service_matrix or args.independent_samples:
            raise ValueError('BAF continuation uses its saved settings; create a new project to change them')
        command=[sys.executable,str(project/'scripts/run_baf.py'),*saved['command'][2:]]
        return execute(command,project,environment(project,args,saved.get('r_library'),saved.get('rscript')),state_path,saved,args.print_command)
    verify_recorded_inputs(active_output, args.hash_raw)
    invalid = verify_delivery(active_output)
    if not args.outdir and not args.service_matrix and not invalid and (not saved or saved.get('engine') != 'nextflow'):
        if saved and saved.get('config') and digest(saved['config']) != saved.get('config_sha256'):
            invalid = ['study configuration changed']
        else:
            print('Completed outputs verified. No analysis needs repeating.')
            print('Results: '+str(active_output/'final_results/START_HERE.html'))
            return 0
    if saved and saved.get('engine') in ('native-quantification', 'native-counts', 'native-benchmark') and saved.get('command'):
        if args.service_matrix or args.independent_samples:
            raise ValueError('This continuation reuses its recorded arguments; use quantify or counts to change the analysis')
        output = Path(args.outdir).expanduser().resolve() if args.outdir else fresh_output(folder)
        if output == folder or (output.exists() and any(output.iterdir())):
            raise ValueError('Use a new output folder; previous results are preserved')
        script = {'native-quantification':'analyze_quantification.py','native-counts':'reanalyze_counts.py','native-benchmark':'benchmark_service.py'}[saved['engine']]
        command = [sys.executable, str(project/'scripts'/script), *saved['command'][2:]]
        command = replace_option(command, '--outdir', output/'final_results' if saved['engine'] in ('native-quantification','native-benchmark') else output)
        env = environment(project,args,saved.get('r_library'),saved.get('rscript'))
        record = dict(saved, outdir=str(output), continued_from=str(folder), project=str(project),
                      r_library=env.get('R_LIBS','').split(os.pathsep)[0] or None, rscript=env.get('PROTEOMICS_RSCRIPT'))
        return execute(command, project, env, state_path, record, args.print_command)
    settings = None if saved and saved.get('engine') == 'nextflow' else (dia_settings(active_output) or dia_settings(folder))
    if settings:
        output = Path(args.outdir).expanduser().resolve() if args.outdir else fresh_output(folder)
        if output == folder or (output.exists() and any(output.iterdir())):
            raise ValueError('Use a new output folder; previous results are preserved')
        env = environment(project, args, (saved or {}).get('r_library') or settings.get('r_library'), (saved or {}).get('rscript'))
        command = [sys.executable, str(project/'scripts/reanalyze_dia.py'), '--source-run', settings['source_run'], '--outdir', str(output)]
        for flag, value in [('--samples',settings.get('samples')), ('--method',settings.get('method')),
                            ('--rscript',env.get('PROTEOMICS_RSCRIPT'))]:
            if value:
                command += [flag,value]
        if args.independent_samples or settings['independent']:
            command.append('--independent-samples')
        if args.hash_raw or settings.get('hash_raw'):
            command.append('--hash-raw')
        for service in args.service_matrix or settings['services']:
            name, path = service.split('=',1)
            command += ['--service-matrix',name+'='+str(Path(path).expanduser().resolve())]
        options = settings.get('differential', {})
        for name in ('fdr', 'min_replicates'):
            if name in options:
                command += ['--'+name.replace('_','-'),str(options[name])]
        record = dict(engine='native-dia', source_run=settings['source_run'], outdir=str(output), continued_from=str(folder), project=str(project),
                      r_library=env.get('R_LIBS','').split(os.pathsep)[0] or None, rscript=env.get('PROTEOMICS_RSCRIPT'))
        return execute(command, project, env, state_path, record, args.print_command)
    if saved and saved['engine'] == 'nextflow':
        env = environment(project,args)
        command = list(saved['command'])
        command[0] = nextflow_binary(Path(saved.get('project', project)),args,env)
        if '-resume' in command:
            position = command.index('-resume'); del command[position:position+2]
        session = saved.get('session') or nextflow_session(saved['launch_dir'], saved['command'])
        if session and not (Path(saved['launch_dir'])/'.nextflow/cache'/session).is_dir():
            session = None
        if not session:
            raise ValueError('Nextflow cache/session missing; start a fresh run with run --config or --dia-config')
        command += ['-resume', session]
        if args.outdir or args.service_matrix or args.independent_samples:
            raise ValueError('Nextflow continuation uses its recorded configuration and output folder; edit its study JSON to change analysis options')
        if not args.print_command:
            final = active_output/'final_results'
            # Remove no data: quarantine damaged published files before cached tasks republish them.
            for relative in invalid:
                path = final/relative
                if path.is_file() and path.resolve().is_relative_to(final.resolve()):
                    target = state_path.parent/'quarantine'/now().replace(':','-')/relative
                    target.parent.mkdir(parents=True,exist_ok=True);path.rename(target)
        if saved.get('config'):
            saved['config_sha256'] = digest(saved['config'])
        return execute(command, saved['launch_dir'],env,state_path,saved,args.print_command)
    legacy = latest_legacy_nextflow(folder)
    if legacy:
        if args.outdir or args.service_matrix or args.independent_samples:
            raise ValueError('Nextflow continuation uses its recorded configuration and output folder; edit its study JSON to change analysis options')
        launch = Path(legacy['root'])
        session = nextflow_session(launch, legacy['command'])
        if not session:
            raise ValueError('Original Nextflow session/work cache is unavailable: '+str(launch))
        env = environment(launch,args)
        command = list(legacy['command']);command[0]=nextflow_binary(launch,args,env)
        if '-resume' in command:
            position=command.index('-resume'); del command[position:position+2]
        command += ['-resume',session]
        record=dict(engine='nextflow',outdir=str(folder),project=str(launch),session=session)
        return execute(command,launch,env,state_path,record,args.print_command)
    final = folder/'final_results'
    if (final/'main_results/counts').is_dir():
        if args.service_matrix or args.independent_samples:
            raise ValueError('Count continuation uses the delivered sample metadata; use quantify for service quantification')
        output = Path(args.outdir).expanduser().resolve() if args.outdir else fresh_output(folder)
        command=[sys.executable,str(project/'scripts/reanalyze_counts.py'),'--results',str(final),'--outdir',str(output)]
        return execute(command,project,environment(project,args),state_path,dict(engine='native-counts',outdir=str(output),source_run=str(folder)),args.print_command)
    raise ValueError('No recognized run state or completed search found in '+str(folder)+'. Supply a run root or its final_results folder.')


def main(argv=None):
    p = parser()
    tokens = list(sys.argv[1:] if argv is None else argv)
    if '--interactive' in tokens or '--interactively' in tokens:
        index=next(i for i,v in enumerate(tokens) if v in ('--interactive','--interactively'))
        tokens=tokens[:index]+['setup']+tokens[index+1:]
    native = {'quantify', 'counts', 'matrices', 'raw-baf', 'benchmark'}
    value_flags = {'--workflow-dir','--project','--cache-dir','--rscript','--r-library','--nextflow','--continue','--outdir','--service-matrix'}
    boundary = None
    skip = False
    for index, token in enumerate(tokens):
        if skip:
            skip = False
            continue
        if token in value_flags:
            skip = True
        elif token in {'setup','check','system','status','demo','run','ui','doctor','install','workflow-path',*native}:
            if token in native: boundary = index
            break
    args = p.parse_args(tokens if boundary is None else tokens[:boundary+1])
    if boundary is not None:
        args.arguments = tokens[boundary+1:]
        if args.arguments[:1] == ['--']:
            args.arguments = args.arguments[1:]
    if args.continue_folder and args.command:
        p.error('Use --continue FOLDER or a subcommand, not both')
    if not args.continue_folder and not args.command:
        print('fun_proteomics | Proteomics analysis\n\nStart:     fun_proteomics setup --project /data/my-study\nTerminal:  add --terminal\nContinue:  fun_proteomics --continue /data/my-study\nCheck tools: fun_proteomics doctor\nGuides and videos: https://cfarkas.github.io/fun_proteomics/\nAll commands: fun_proteomics --help');return 0
    try:
        project = resolve_project(args.project, args.cache_dir)
        if args.command == 'workflow-path':
            print(project);return 0
        from . import study
        if args.command in ('setup', 'check', 'system', 'status', 'demo'):
            return study.dispatch(args, project)
        if args.command == 'doctor':
            from .requirements import doctor
            return doctor(project,args)
        if args.command == 'install':
            from .requirements import install
            return install(project,args)
        if args.continue_folder:
            config = Path(args.continue_folder).expanduser()/'config/run.json'
            if config.is_file() and study.is_study(config):
                if args.outdir or args.service_matrix or args.independent_samples:
                    raise ValueError('A guided study uses its saved settings; create a new project to change inputs or design')
                return study.run_study(config, args, project)
            return continue_run(args,project)
        env=environment(project,args)
        if args.command in ('quantify','counts','matrices','raw-baf','benchmark'):
            script={'quantify':'scripts/analyze_quantification.py','counts':'scripts/reanalyze_counts.py','matrices':'pipeline_lib/matrix_catalog.py','raw-baf':'scripts/run_baf.py','benchmark':'scripts/benchmark_service.py'}[args.command]
            command=[sys.executable,str(project/script),*args.arguments]
            if args.command == 'matrices' or any(value in args.arguments for value in ('-h', '--help')):
                if args.print_command:
                    print(shlex.join(command));return 0
                return subprocess.run(command,env=env).returncode
            paths = ('--matrix','--samples','--fasta','--service-matrix','--service-samples','--results','--config','--outdir')
            command = absolute_options(command, paths)
            destinations = option_values(command, '--outdir')
            if len(destinations) != 1:
                raise ValueError('Supply exactly one --outdir for this analysis')
            output = Path(destinations[0])
            if args.command!='raw-baf' and output.exists() and any(output.iterdir()):
                raise ValueError('Output folder is not empty; use --continue '+str(output))
            inputs = [value for flag in paths if flag != '--outdir' for value in option_values(command,flag)]
            record = dict(engine={'quantify':'native-quantification','raw-baf':'native-baf','benchmark':'native-benchmark','counts':'native-counts'}[args.command],
                          outdir=str(output), project=str(project), input_provenance=recorded_files(inputs),
                          r_library=env.get('R_LIBS','').split(os.pathsep)[0] or None, rscript=env.get('PROTEOMICS_RSCRIPT'))
            if args.command in ('quantify','benchmark'):
                command = replace_option(command, '--outdir', output/'final_results')
            return execute(command,project,env,output/'.fun_proteomics/run.json',record,args.print_command)
        if args.command == 'ui':
            if not args.advanced and (not args.config or study.is_study(Path(args.config).expanduser())):
                from .setup_web import serve
                return serve(args,project)
            args.port = args.port or 8765
            command=[nextflow_binary(project,args,env),'run',str(project/'main.nf'),'--interactively','--ui_port',str(args.port)]
            if args.config:command+=['--config',str(Path(args.config).expanduser().resolve())]
            if args.print_command:print(shlex.join(command));return 0
            return subprocess.run(command,env=env).returncode
        if args.command == 'run':
            if args.config and study.is_study(Path(args.config).expanduser()):
                if args.run_outdir:
                    raise ValueError('The saved project defines its results folder; omit --outdir')
                return study.run_study(Path(args.config).expanduser(),args,project,check_only=args.validate_only)
            if not args.run_outdir:
                raise ValueError('Legacy raw configurations require --outdir; guided studies use setup --project FOLDER')
            output=Path(args.run_outdir).expanduser().resolve()
            state=output/'.fun_proteomics/run.json'
            if output.exists() and any(output.iterdir()):
                raise ValueError('Output folder is not empty; use --continue '+str(output))
            config=Path(args.config or args.dia_config).expanduser().resolve()
            if not config.is_file():raise ValueError('Configuration does not exist: '+str(config))
            controller=output/'.fun_proteomics'
            command=[nextflow_binary(project,args,env),'-log',str(controller/'nextflow.log'),'run',str(project/'main.nf'),
                     '--config' if args.config else '--dia_config',str(config),'--outdir',str(output),'-work-dir',str(controller/'work'),'-ansi-log','false']
            if args.validate_only:command.append('--validate_only')
            return execute(command,controller/'nextflow',env,state,dict(engine='nextflow',project=str(project),outdir=str(output),
                           config=str(config),config_sha256=digest(config)),args.print_command)
    except (ValueError, OSError, KeyError, json.JSONDecodeError) as error:
        print('fun_proteomics: '+str(error),file=sys.stderr);return 2
    except (KeyboardInterrupt, EOFError):
        print('\nSetup interrupted. Saved inputs and results are preserved.',file=sys.stderr);return 130
    return 0
