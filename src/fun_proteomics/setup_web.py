"""Local guided setup with a session-bound API and resumable native launch."""
import copy
import html
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
import mimetypes
import os
from pathlib import Path
import secrets
import signal
import subprocess
import sys
import threading
from urllib.parse import unquote, urlparse
import webbrowser

from . import study
from .continuation import read_json, write_json
from .setup_terminal import from_args


class Application:
    def __init__(self,args,workflow):
        self.args=args;self.workflow=workflow;self.token=secrets.token_urlsafe(32)
        self.initial=from_args(args);self.project=self.initial.get('project')
        if getattr(args,'config',None): self.project=str(Path(args.config).expanduser().resolve().parent.parent)
        if self.project:
            self.project=str(Path(self.project).expanduser().resolve())
            self.initial['project']=self.project
            config=Path(self.project)/'config/run.json'
            if config.is_file(): self.initial=study.load(config)
        self.process=None;self.guard=threading.Lock()

    def load(self,project):
        cfg=study.load(Path(project).expanduser()/'config/run.json')
        self.project=cfg['project'];self.initial=cfg
        return {'config':cfg,'status':self.status()}

    def status(self):
        if not self.project or not (Path(self.project)/'config/run.json').is_file(): return {'state':'new','log_tail':''}
        info=study.status(self.project)
        log=Path(info['log'])
        if log.is_file():
            with log.open('rb') as handle:
                handle.seek(max(0,log.stat().st_size-12000));info['log_tail']=handle.read().decode(errors='replace')
        else: info['log_tail']=''
        record=Path(self.project)/'.fun_proteomics/browser-launch.json'
        if record.is_file():
            from .cli import process_identity
            entry=read_json(record)
            if entry.get('identity') and process_identity(entry['pid'])==entry['identity']:
                info['state']='running'
            elif entry.get('exit_code'):
                info['state']='failed';info['exit_code']=entry['exit_code']
        return info

    def start(self):
        from .cli import environment, process_identity
        with self.guard:
            if self.status()['state']=='running': raise ValueError('This project is already running')
            path=Path(self.project)/'config/run.json'
            study.check(study.load(path),self.workflow)
            folder=Path(self.project)/'logs';folder.mkdir(parents=True,exist_ok=True)
            command=[sys.executable,'-m','fun_proteomics',*study.runtime_flags(self.args),'run','--config',str(path)]
            env=environment(self.workflow,self.args)
            env['PYTHONPATH']=str(Path(__file__).resolve().parents[1])+os.pathsep+env.get('PYTHONPATH','')
            with (folder/'analysis.log').open('a') as log:
                self.process=subprocess.Popen(command,cwd=self.project,env=env,
                             stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            record=Path(self.project)/'.fun_proteomics/browser-launch.json'
            entry=dict(pid=self.process.pid,identity=process_identity(self.process.pid))
            write_json(record,entry)
            process=self.process
            def finish():
                code=process.wait()
                with self.guard:
                    if read_json(record).get('pid')==process.pid: write_json(record,dict(entry,exit_code=code))
            threading.Thread(target=finish,daemon=True).start()
            return self.status()

    def stop(self):
        from .cli import process_identity
        with self.guard:
            path=Path(self.project)/'.fun_proteomics/browser-launch.json'
            if not path.is_file(): raise ValueError('This browser has no recorded launcher to stop')
            entry=read_json(path)
            if not entry.get('identity') or process_identity(entry['pid'])!=entry['identity']:
                raise ValueError('The recorded launcher is no longer running')
            os.killpg(entry['pid'],signal.SIGINT)
            return {'message':'Stop requested. Keep this project and use Run / continue to recover.'}


def make_server(app,port=0):
    class Handler(BaseHTTPRequestHandler):
        def log_message(self,*_): pass
        def reply(self,value,status=200,mime='application/json',cookie=False):
            payload=value if isinstance(value,bytes) else value.encode() if isinstance(value,str) else json.dumps(value,allow_nan=False).encode()
            self.send_response(status);self.send_header('Content-Type',mime)
            self.send_header('Content-Length',str(len(payload)));self.send_header('Cache-Control','no-store')
            self.send_header('X-Content-Type-Options','nosniff');self.send_header('X-Frame-Options','DENY')
            self.send_header('Referrer-Policy','no-referrer')
            if cookie:self.send_header('Set-Cookie','proteomics_session='+app.token+'; HttpOnly; SameSite=Strict; Path=/')
            self.end_headers();self.wfile.write(payload)
        def trusted(self,token=True,cookie=False):
            allowed={f'127.0.0.1:{self.server.server_port}',f'localhost:{self.server.server_port}'}
            if self.headers.get('Host') not in allowed:return False
            if self.headers.get('Origin') and self.headers['Origin'] not in {'http://'+v for v in allowed}:return False
            if self.headers.get('Sec-Fetch-Site')=='cross-site':return False
            supplied=self.headers.get('X-Session-Token','')
            if cookie:
                from http.cookies import SimpleCookie
                cookies=SimpleCookie(self.headers.get('Cookie',''))
                supplied=cookies['proteomics_session'].value if 'proteomics_session' in cookies else ''
            return not token or secrets.compare_digest(supplied,app.token)
        def do_GET(self):
            route=urlparse(self.path).path
            if route=='/':
                if not self.trusted(False):return self.reply({'error':'Invalid host or origin'},403)
                return self.reply((app.workflow/'assets/setup.html').read_bytes(),mime='text/html; charset=utf-8')
            if route.startswith('/results/'):
                if not self.trusted(cookie=True):return self.reply({'error':'Open results from this setup session'},403)
                try:
                    cfg=study.load(Path(app.project)/'config/run.json');root,_=study.result_location(cfg);final=(root/'final_results').resolve()
                    path=(final/unquote(route.removeprefix('/results/'))).resolve()
                    if not path.is_relative_to(final) or not path.is_file(): raise ValueError('Result not found')
                    return self.reply(path.read_bytes(),mime=mimetypes.guess_type(path.name)[0] or 'application/octet-stream')
                except (ValueError,OSError):return self.reply({'error':'Result not found'},404)
            if not self.trusted():return self.reply({'error':'Invalid session; reopen the terminal URL'},403)
            try:
                if route=='/api/initial':return self.reply({'config':app.initial,'cwd':str(Path.cwd()),'status':app.status()},cookie=True)
                if route=='/api/status':return self.reply(app.status())
                return self.reply({'error':'Not found'},404)
            except (ValueError,OSError) as error:return self.reply({'error':str(error)},400)
        def do_POST(self):
            if not self.trusted():return self.reply({'error':'Invalid session or origin'},403)
            try:
                size=int(self.headers.get('Content-Length','0'))
                if not 0<size<=2_000_000:raise ValueError('Invalid request size')
                data=json.loads(self.rfile.read(size))
                if not isinstance(data,dict):raise ValueError('Expected an object')
                if self.path=='/api/browse':
                    path=Path(data.get('path') or Path.cwd()).expanduser().resolve()
                    if path.is_file():path=path.parent
                    if not path.is_dir():raise ValueError('Directory does not exist')
                    entries=sorted((p for p in path.iterdir() if not p.name.startswith('.')),key=lambda p:(not p.is_dir(),p.name.lower()))
                    return self.reply({'path':str(path),'parent':str(path.parent),'entries':[{'name':p.name,'path':str(p),'directory':p.is_dir()} for p in entries[:500]],'truncated':len(entries)>500})
                if self.path=='/api/scan':return self.reply({'samples':study.scan(data['path'])})
                if self.path=='/api/load':return self.reply(app.load(data['project']))
                if self.path=='/api/check':return self.reply(study.check(data['config'],app.workflow))
                if self.path=='/api/save':
                    with app.guard:
                        if app.status()['state']=='running':raise ValueError('An analysis is running; settings are locked')
                        path,report=study.save(data['config'],app.workflow)
                        app.project=str(path.parent.parent);app.initial=study.load(path)
                    return self.reply({'config':app.initial,'report':report,'path':str(path)})
                if self.path=='/api/run':return self.reply(app.start())
                if self.path=='/api/stop':return self.reply(app.stop())
                return self.reply({'error':'Not found'},404)
            except (ValueError,OSError,KeyError,TypeError) as error:return self.reply({'error':str(error)},400)
    return ThreadingHTTPServer(('127.0.0.1',port),Handler)


def serve(args,workflow):
    app=Application(args,workflow);server=make_server(app,getattr(args,'port',0))
    url=f'http://127.0.0.1:{server.server_port}/#'+app.token
    print('fun_proteomics | Guided setup\n'+url,flush=True)
    print('Keep this terminal open for setup. Launched analyses keep running if it closes.',flush=True)
    if not getattr(args,'no_browser',False):webbrowser.open(url)
    try:server.serve_forever(poll_interval=.25)
    except KeyboardInterrupt:print('\nSetup closed. Reopen the same project to view progress or continue.')
    finally:server.server_close()
    return 0
