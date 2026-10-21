"""Control one durable Nextflow run and expose its on-disk log to the browser."""
import collections
import fcntl
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import threading
import time
from ui_run_state import ACTIVE, process_matches, read_state, write_json
from dia_browser import is_dia_study, prepare_browser_study

class Runner:
    def __init__(self, root, command_prefix=None):
        self.root=Path(root).resolve(); self.runtime=Path(os.environ.get('PROTEOMICS_RUNTIME_DIR',str(self.root/'.runtime'))); self.command_prefix=command_prefix or [str(self.runtime/'env/bin/nextflow')]
        self.lock=threading.RLock(); self.proc=None; self.lines=collections.deque(maxlen=2000)
        self.started=None;self.ended=None;self.state='idle';self.exit_code=None;self.command=[];self.log_path=None;self.serial=0
        self.index=self.root/'.runtime/ui_run.json';self.state_path=None;self.record={};self.offset=0;self.pending=b''

    def _refresh(self):
        if not self.index.is_file():return
        index=json.loads(self.index.read_text())
        self.state_path=Path(index['state_path'])
        self.record=read_state(self.state_path)
        record=self.record
        path=Path(record['log_path'])
        if path!=self.log_path:
            self.log_path=path;self.offset=0;self.pending=b'';self.lines.clear();self.serial=0
        self.state=record['state'];self.started=record['started'];self.ended=record.get('ended')
        self.exit_code=record.get('exit_code');self.command=record['command']
        if path.is_file():
            if path.stat().st_size<self.offset:self.offset=0;self.pending=b'';self.lines.clear()
            with path.open('rb') as handle:
                handle.seek(self.offset)
                chunk=handle.read()
                self.offset=handle.tell()
            lines=(self.pending+chunk).split(b'\n');self.pending=lines.pop()
            if self.state not in ACTIVE and self.pending:lines.append(self.pending);self.pending=b''
            for line in lines:self.lines.append(line.decode('utf-8',errors='replace'));self.serial+=1
    def external_run(self):
        # Reject another controller launched manually in this project, including
        # absolute-main.nf launches from another directory.
        for entry in Path('/proc').iterdir():
            if not entry.name.isdigit():continue
            try:
                # /proc terminates argv with NUL; it is not an empty argument
                # following a bare --interactively flag.
                argv=(entry/'cmdline').read_bytes().decode(errors='replace').rstrip('\x00').split('\x00')
                args=' '.join(argv)
                # The interactive Nextflow launcher owns the web server, not a search.
                if '--interactively' in argv:
                    i=argv.index('--interactively')
                    if i+1==len(argv) or argv[i+1].startswith('-') or argv[i+1].lower()=='true':continue
                if '--interactively=true' in argv:continue
                if 'nextflow' not in args or 'run' not in argv or not any(a=='main.nf' or a.endswith('/main.nf') for a in argv):continue
                command=argv.index('run')
                if command+1>=len(argv):continue
                workflow=Path(argv[command+1])
                if not workflow.is_absolute():workflow=(entry/'cwd').resolve()/workflow
                if workflow.resolve()==self.root/'main.nf':return int(entry.name)
            except (OSError,ValueError):pass
        return None
    def resume_candidate(self):
        history=self.root/'.nextflow/history'
        empty={'available':False,'message':'No incomplete analysis detected.'}
        if not history.is_file():return empty
        managed={}
        if self.index.is_file():
            try:managed=read_state(Path(json.loads(self.index.read_text())['state_path']))
            except (KeyError,OSError,ValueError):pass
        dia_managed=managed.get('workflow')=='dia'
        if dia_managed and (managed.get('state')=='completed' or managed.get('kind')!='analysis'):return empty
        for line in reversed(history.read_text().splitlines()):
            fields=line.split('\t')
            if len(fields)<7:continue
            try:argv=shlex.split(fields[6])
            except ValueError:continue
            if dia_managed:
                native=None
                for i,arg in enumerate(argv):
                    if arg=='--dia_config' and i+1<len(argv):native=Path(argv[i+1])
                    if arg.startswith('--dia_config='):native=Path(arg.split('=',1)[1])
                if native is None:continue
                if not native.is_absolute():native=self.root/native
                if native.resolve()!=Path(managed['analysis_config']).resolve():continue
                if fields[3]=='OK':return empty
                config=Path(managed['config']);cache=self.root/'.nextflow/cache'/fields[5]
                available=config.is_file() and native.is_file() and cache.is_dir()
                return {'available':available,'name':fields[2],'session':fields[5],'config':str(config),
                        'dia_config':str(native),'message':('Incomplete DIA analysis: '+fields[2]+'. Resume reuses matching completed tasks.' if available else 'The DIA study or Nextflow cache is missing.')}
            excluded={'--interactively','--help','--validate_only','--dry_run','--report_only','--dia_config','-stub-run','-stub'}
            if any(a.split('=')[0] in excluded for a in argv):continue
            if not any(a=='main.nf' or a.endswith('/main.nf') for a in argv):continue
            if fields[3]=='OK':return empty
            config=self.root/'config.local.json'
            workdir=None
            for i,arg in enumerate(argv):
                if arg in ('--config','--config_file') and i+1<len(argv):config=Path(argv[i+1])
                if arg.startswith(('--config=','--config_file=')):config=Path(arg.split('=',1)[1])
                if arg in ('-w','-work-dir') and i+1<len(argv):workdir=argv[i+1]
            if not config.is_absolute():config=self.root/config
            cache=self.root/'.nextflow/cache'/fields[5]
            available=config.is_file() and cache.is_dir()
            return {'available':available,'name':fields[2],'session':fields[5],'config':str(config),
                    'workdir':workdir,'message':('Incomplete analysis: '+fields[2]+'. Resume reuses matching completed tasks; unfinished tasks restart.' if available else 'Previous analysis found, but its configuration or Nextflow cache is missing.')}
        return empty
    def resume_session(self):
        candidate=self.resume_candidate()
        if not candidate['available']:raise ValueError(candidate['message'])
        return candidate['session']
    def resume_previous(self):
        with self.lock:
            candidate=self.resume_candidate()
            if not candidate['available']:raise ValueError(candidate['message'])
            return self.start(candidate['config'],resume=True,dia_config=candidate.get('dia_config'),resume_session=candidate['session'])
    def start(self, config, resume=False, validate_only=False, report_only=False, dia_config=None, resume_session=None):
        self.index.parent.mkdir(parents=True,exist_ok=True)
        # Two browser servers on different ports must not launch duplicate runs.
        with self.lock, self.index.with_suffix('.lock').open('a') as handle:
            fcntl.flock(handle,fcntl.LOCK_EX)
            return self._start(config,resume,validate_only,report_only,dia_config,resume_session)

    def _start(self, config, resume=False, validate_only=False, report_only=False, dia_config=None, resume_session=None):
        with self.lock:
            self._refresh()
            if self.state in ('running','stopping'):raise ValueError('A run is already active. Stop it before starting another.')
            external=self.external_run()
            if external:raise ValueError(f'Nextflow is already running for this project (PID {external}).')
            config=Path(config)
            if not config.is_absolute():config=self.root/config
            config=config.resolve()
            study=json.loads(config.read_text())
            dia=is_dia_study(study)
            if validate_only and report_only:raise ValueError('Choose input checking or report-only analysis')
            if dia and report_only:raise ValueError('Report-only DDA statistics are unavailable for DIA-NN.')
            if dia_config is not None and not dia:raise ValueError('A DIA configuration requires a confirmed DIA browser study.')
            analysis_config=Path(dia_config).resolve() if dia_config else prepare_browser_study(config,self.root) if dia else config
            if dia and not analysis_config.is_file():raise ValueError('The saved DIA input configuration is missing.')
            output=Path(study['outdir']).expanduser()
            if not output.is_absolute():output=config.parent/output
            logs=output.resolve()/'logs';logs.mkdir(parents=True,exist_ok=True)
            run_id=str(time.time_ns())
            log_path=logs/f'run_{run_id}.log'
            self.command=self.command_prefix+['-log',str(logs/f'nextflow_{run_id}.log'),'run',str(self.root/'main.nf'),
                '--dia_config' if dia else '--config',str(analysis_config)]
            if not dia:self.command+=['--verbosity','4']
            self.command+=['-ansi-log','false','-with-trace',str(logs/f'nextflow_{run_id}_trace.tsv'),
                '-with-report',str(logs/f'nextflow_{run_id}_report.html'),'-with-timeline',str(logs/f'nextflow_{run_id}_timeline.html')]
            if resume:
                self.command+=['-resume']+([resume_session] if resume_session else []) if dia else ['--continue']
            if validate_only:self.command+=['--validate_only']
            if report_only:self.command+=['--report_only']
            env=os.environ.copy();env['PATH']=str(self.runtime/'env/bin')+':'+env['PATH']
            state_path=logs/f'run_{run_id}.json'
            record={'state':'running','root':str(self.root),'config':str(config),'analysis_config':str(analysis_config),'workflow':'dia' if dia else 'dda','command':self.command,
                    'log_path':str(log_path),'state_path':str(state_path),'stop_path':str(state_path.with_suffix('.stop')),
                    'started':time.time(),'ended':None,'exit_code':None,'pid':None,'supervisor_pid':None,
                    'kind':'validation' if validate_only else 'differential' if report_only else 'analysis'}
            write_json(state_path,record)
            write_json(self.index,{'state_path':str(state_path)})
            with log_path.open('a') as log:
                try:
                    self.proc=subprocess.Popen([sys.executable,str(Path(__file__).with_name('ui_run_state.py')),str(state_path)],
                        cwd=self.root,env=env,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
                except OSError as error:
                    record.update(state='failed',ended=time.time(),message=str(error))
                    write_json(state_path,record)
                    raise
            threading.Thread(target=self.proc.wait,daemon=True).start()
            return self.status()

    def stop(self):
        with self.lock:
            self._refresh()
            if self.state not in ('running','stopping'):
                state=self.status()
                state['message']='Analysis already completed. There is no active run to stop.' if self.state=='completed' else 'There is no active run to stop.'
                return state
            if self.state=='stopping':return self.status()
            if self.record.get('supervisor_pid') and not process_matches(self.record['supervisor_pid'],self.record.get('supervisor_identity')):
                raise ValueError('Run supervisor is unavailable. Inspect the run log before stopping the analysis process manually.')
            Path(self.record['stop_path']).touch()
            return self.status()

    def status(self):
        with self.lock:
            self._refresh()
            return {'state':self.state,'pid':self.record.get('pid'),'exit_code':self.exit_code,
                    'elapsed':round((self.ended or time.time())-self.started) if self.started else 0,
                    'resume':self.resume_candidate(),'kind':'validation' if '--validate_only' in self.command else 'differential' if '--report_only' in self.command else 'analysis','command':self.command,'log_path':str(self.log_path) if self.log_path else '',
                    'state_path':str(self.state_path) if self.state_path else '', 'config':self.record.get('config'),
                    'message':self.record.get('message',''),
                    'lines':list(self.lines),'serial':self.serial}
