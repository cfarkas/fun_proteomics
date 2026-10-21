#!/usr/bin/env python3
"""Local study editor with controlled Nextflow launch, stop and live logs."""
import argparse
import copy
import json
import os
import signal
import subprocess
import sys
import threading
import urllib.request
import webbrowser
from pathlib import Path
import secrets
import shlex
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlparse
import jsonschema
from analysis_modules import status as module_status, install as install_module
from ui_log_follow import follow, mirror_runner
from ui_server_lifecycle import start_server
from differential_tools import capabilities as differential_capabilities
from ui_runner import Runner
from dia_browser import is_dia_study, validate_browser_study, prepare_browser_study
from native_folders import choose_folder
from tool_manager import discover as discover_tools, install as install_tool
from sample_detection import detect, ms_format
from system_capabilities import capabilities
from proteomics_detection import detect_proteomics, confirm_config

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'pipeline_lib'))
from optional_analyses import validate_options, migrate_options
TOKEN=secrets.token_urlsafe(32)

def absolute_config(path):
    cfg=json.loads(path.read_text())
    if path.name.startswith('draft_'):return cfg
    if not is_dia_study(cfg):migrate_options(cfg)
    for s in cfg.get('samples', []):
        if isinstance(s.get('raw_path'),str) and s['raw_path']:s['raw_path']=str((path.parent/Path(s['raw_path']).expanduser()).resolve())
    for group,keys in [('database',['fasta','contaminants_zip']),('tools',list(cfg.get('tools', {}))),('dia',['canonical_fasta','contaminant_fasta','install_cache','raw_root'])]:
        for k in keys:
            if cfg.get(group, {}).get(k): cfg[group][k]=str((path.parent/Path(cfg[group][k]).expanduser()).resolve())
    if cfg.get('dia', {}).get('service_matrices'):
        cfg['dia']['service_matrices'] = {name: str((path.parent / Path(value).expanduser()).resolve()) for name, value in cfg['dia']['service_matrices'].items()}
    if cfg.get('outdir'):cfg['outdir']=str((path.parent/cfg['outdir']).resolve())
    return cfg

def check(cfg):
    if is_dia_study(cfg):
        schema=json.loads((ROOT/'assets/config.schema.json').read_text())
        schema['required']=['outdir','samples','proteomics','dia']
        schema.pop('allOf',None)
        # Retained DDA choices are inactive, so their validation must not block DIA.
        for name in ('database','tools','modes','differential','optional_analyses','immunopeptidomics'):
            schema['properties'][name]={}
        jsonschema.Draft202012Validator(schema).validate(cfg)
        raw_base=(ROOT/Path(cfg['dia'].get('raw_root',str(ROOT))).expanduser()).resolve()
        for sample in cfg['samples']:sample['raw_path']=str((raw_base/Path(sample['raw_path']).expanduser()).resolve())
        for group,keys in [('database',['fasta','contaminants_zip']),('dia',['canonical_fasta','contaminant_fasta','install_cache','raw_root'])]:
            for key in keys:
                if cfg.get(group,{}).get(key):cfg[group][key]=str((ROOT/Path(cfg[group][key]).expanduser()).resolve())
        cfg['outdir']=str((ROOT/Path(cfg['outdir']).expanduser()).resolve())
        confirm_config(cfg)
        limits=capabilities();resources=cfg.get('resources',{'cpus':32,'memory_gib':96})
        if resources['cpus']>limits['available_cpus']:raise ValueError('Threads exceed available CPUs: '+str(limits['available_cpus']))
        if resources['memory_gib']>limits['max_memory_gib']:raise ValueError('Memory exceeds currently available RAM. Return to Study settings and detect resources again.')
        validate_browser_study(cfg,ROOT)
        return cfg
    migrate_options(cfg)
    jsonschema.Draft202012Validator(json.loads((ROOT/'assets/config.schema.json').read_text())).validate(cfg)
    validate_options(cfg)
    confirm_config(cfg)
    de=cfg.get('differential',{})
    available=differential_capabilities()
    if not available['available']:raise ValueError('Final matrices require R, edgeR, limma, jsonlite and statmod. Run bash scripts/install_differential_tools.sh first.')
    if de.get('enabled'):
        missing=[m for m in de.get('methods',['edgeR','limma']) if not available['methods'].get(m)]
        if missing:raise ValueError('Install missing differential dependencies: '+', '.join(missing))
        for cohort in {s['cohort'] for s in cfg['samples']}:
            counts={}
            for sample in cfg['samples']:
                if sample['cohort']==cohort:counts[sample['condition']]=counts.get(sample['condition'],0)+1
            if len(counts)<2 or min(counts.values())<2:raise ValueError('Differential proteins require at least two conditions and two independent replicates per condition in '+cohort)
    limits=capabilities()
    resources=cfg.get('resources',{'cpus':32,'memory_gib':128})
    if resources['cpus']>limits['available_cpus']: raise ValueError('Threads exceed available CPUs: '+str(limits['available_cpus']))
    if resources['memory_gib']>limits['max_memory_gib']: raise ValueError('Memory exceeds currently available RAM. Return to Study settings and detect resources again.')
    ids=set(); paths=set()
    for s in cfg['samples']:
        if s['cohort']=='unassigned': raise ValueError('Assign a cohort for '+s['sample_id']+' before running; detection could not infer it reliably.')
        if s['sample_id'] in ids: raise ValueError('Duplicate sample ID: '+s['sample_id'])
        ids.add(s['sample_id']); raw=Path(s['raw_path']).expanduser().resolve()
        if raw in paths: raise ValueError('Same .d directory used twice: '+str(raw))
        paths.add(raw)
        if raw.suffix.lower()!='.d' or not raw.is_dir(): raise ValueError('Not a server .d folder: '+str(raw))
        for name in ['analysis.tdf','analysis.tdf_bin']:
            if not (raw/name).is_file(): raise ValueError('Missing '+str(raw/name))
        s['raw_path']=str(raw)
    for group,keys in [('database',['fasta','contaminants_zip']),('tools',list(cfg['tools']))]:
        for k in keys:
            if k in cfg[group]:
                path=Path(cfg[group][k]).expanduser()
                if not path.is_absolute() or not path.is_file(): raise ValueError(f'{k}: supply an existing absolute server path')
    if not Path(cfg['outdir']).is_absolute(): cfg['outdir']=str((ROOT/cfg['outdir']).resolve())
    return cfg

class Handler(BaseHTTPRequestHandler):
    def log_message(self,fmt,*args): pass
    def reply(self,data,status=200,mime='application/json'):
        payload=data.encode() if isinstance(data,str) else json.dumps(data).encode()
        self.send_response(status); self.send_header('Content-Type',mime); self.send_header('Content-Length',str(len(payload)))
        self.send_header('Cache-Control','no-store'); self.send_header('X-Content-Type-Options','nosniff'); self.end_headers(); self.wfile.write(payload)
    def trusted(self):
        return self.headers.get('Host') in (f'127.0.0.1:{self.server.server_port}',f'localhost:{self.server.server_port}')
    def do_GET(self):
        if not self.trusted(): return self.reply({'error':'Invalid host'},403)
        if self.path=='/': return self.reply((ROOT/'assets/config_ui.html').read_text().replace('__TOKEN__',TOKEN),mime='text/html; charset=utf-8')
        if self.path=='/api/session':
            origin=self.headers.get('Origin')
            if origin and origin not in (f'http://127.0.0.1:{self.server.server_port}',f'http://localhost:{self.server.server_port}'):
                return self.reply({'error':'Invalid origin'},403)
            if self.headers.get('Sec-Fetch-Site')=='cross-site': return self.reply({'error':'Invalid origin'},403)
            return self.reply({'token':TOKEN})
        if self.path=='/api/differential_tools': return self.reply(differential_capabilities())
        if self.path=='/api/modules': return self.reply({'phosphoproteomics':module_status()})
        if self.path=='/api/capabilities': return self.reply(capabilities())
        if self.path=='/api/status': return self.reply(self.server.runner.status())
        if self.path=='/api/config':
            state_path=self.server.config.with_suffix('.ui.json')
            ui_state=json.loads(state_path.read_text()) if state_path.is_file() else {}
            return self.reply({'config':absolute_config(self.server.config),'root':str(ROOT),'config_path':str(self.server.config),'resource_overrides':getattr(self.server,'resource_overrides',{}),'ui_state':ui_state})
        return self.reply({'error':'Not found'},404)
    def do_POST(self):
        if not self.trusted(): return self.reply({'error':'Invalid host'},403)
        if self.headers.get('X-Config-Token')!=TOKEN: return self.reply({'error':'Session changed','code':'TOKEN_STALE'},403)
        origin=self.headers.get('Origin')
        if origin and origin not in (f'http://127.0.0.1:{self.server.server_port}',f'http://localhost:{self.server.server_port}'):
            return self.reply({'error':'Invalid origin'},403)
        try:
            size=int(self.headers.get('Content-Length','0'))
            if not 0<size<=2_000_000: raise ValueError('Invalid request size')
            data=json.loads(self.rfile.read(size))
            if self.path=='/api/detect_tools':
                return self.reply(discover_tools(data.get('tools',{}),data.get('search_root','')))
            if self.path=='/api/install_tool':
                if self.server.runner.status()['state'] in ('running','stopping'): raise ValueError('Stop the analysis before installing tools.')
                return self.reply(install_tool(data.get('name'),data.get('tools',{}),data.get('force') is True))
            if self.path=='/api/install_module':
                if data.get('module')!='phosphoproteomics': raise ValueError('No supported installer for this module yet.')
                if self.server.runner.status()['state'] in ('running','stopping'): raise ValueError('Stop the analysis before installing a module.')
                return self.reply(install_module(data.get('tools',{}),force=data.get('force') is True))
            if self.path=='/api/download_fasta':
                taxon=str(data.get('taxon','')).strip()
                if not taxon.isascii() or not taxon.isdigit() or not 0<int(taxon)<=999999999:
                    raise ValueError('Enter a valid NCBI taxonomy ID')
                reviewed=data.get('reviewed',True)
                if not isinstance(reviewed,bool): raise ValueError('Reviewed must be true or false')
                folder=ROOT/'.database_cache'/'prepared'/secrets.token_hex(12)
                folder.mkdir(parents=True)
                command=[sys.executable,str(ROOT/'scripts/fasta_download.py'),'--taxon',taxon,'--cache',str(ROOT/'.database_cache'/'downloads')]
                if not reviewed: command.append('--unreviewed')
                try:
                    result=subprocess.run(command,cwd=folder,capture_output=True,text=True,timeout=900)
                except subprocess.TimeoutExpired:
                    raise ValueError('Download timed out. Check your connection and try again.')
                if result.returncode:
                    detail=(result.stderr or result.stdout).strip().splitlines()
                    raise ValueError('FASTA download failed: '+(detail[-1] if detail else 'unknown error'))
                info=json.loads((folder/'database_provenance.json').read_text())
                return self.reply({'database':{'fasta':str(folder/'search.fasta'),'sha256':info['search_fasta_sha256'],'contains_decoys':True},'entries':info['entries']})
            if self.path=='/api/detect_proteomics': return self.reply(detect_proteomics(data['samples']))
            if self.path=='/api/detect_samples':
                return self.reply(detect(data['parent'], data.get('existing', []), data.get('recursive', False) is True))
            if self.path=='/api/choose_folder':
                return self.reply(choose_folder(data.get('initial',''),data.get('target','input')))
            if self.path=='/api/focus_browser':
                try:
                    result=subprocess.run([sys.executable,str(ROOT/'scripts/desktop_focus.py')],capture_output=True,text=True,timeout=8)
                    return self.reply(json.loads(result.stdout) if result.returncode==0 else {'focused':False})
                except (subprocess.TimeoutExpired,ValueError): return self.reply({'focused':False})
            if self.path=='/api/resume': return self.reply(self.server.runner.resume_previous())
            if self.path=='/api/stop': return self.reply(self.server.runner.stop())
            if self.path=='/api/directories':
                parent=Path(data['parent']).expanduser().resolve()
                if not parent.is_dir(): raise ValueError('Server directory does not exist')
                return self.reply({'parent':str(parent),'up':str(parent.parent),'directories':[str(p) for p in sorted(parent.iterdir()) if p.is_dir() and not p.name.startswith('.')], 'ms_files':[{'path':str(p),'format':ms_format(p)} for p in sorted(parent.iterdir()) if p.is_file() and ms_format(p)]})
            if self.path=='/api/ensure_output_folder':
                path=Path(data['path']).expanduser().resolve()
                path.mkdir(parents=True,exist_ok=True)
                return self.reply({'path':str(path)})
            if self.path=='/api/mkdir':
                parent=Path(data['parent']).expanduser().resolve()
                name=data['name'].strip()
                if not name or name in ('.','..') or '/' in name or '\\' in name or any(ord(c)<32 for c in name): raise ValueError('Enter one folder name, without slashes')
                target=parent/name
                target.mkdir()  # Do not overwrite or silently reuse existing folders.
                return self.reply({'path':str(target)})
            if self.path=='/api/scan':
                parent=Path(data['parent']).expanduser().resolve()
                if not parent.is_dir(): raise ValueError('Server directory does not exist')
                dirs=sorted(p for p in parent.iterdir() if p.is_dir() and p.suffix.lower()=='.d')
                if len(dirs)>2000: raise ValueError('Too many directories; select a narrower parent')
                return self.reply({'folders':[{'name':p.name,'path':str(p),'valid':all((p/n).is_file() for n in ['analysis.tdf','analysis.tdf_bin'])} for p in dirs]})
            if self.path=='/api/draft':
                cfg=data.get('config');ui_state=data.get('ui_state',{})
                if not isinstance(cfg,dict) or not isinstance(ui_state,dict):raise ValueError('Draft configuration and UI state must be objects.')
                samples=cfg.get('samples',[])
                if not isinstance(samples,list) or len(samples)>2000 or any(not isinstance(sample,dict) for sample in samples):raise ValueError('Draft samples must be a list of at most 2000 sample objects.')
                target=ROOT/'ui_configs';target.mkdir(exist_ok=True)
                path=target/('draft_'+secrets.token_hex(6)+'.json')
                with path.open('x') as handle:json.dump(cfg,handle,indent=2);handle.write('\n')
                with path.with_suffix('.ui.json').open('x') as handle:json.dump(ui_state,handle,indent=2);handle.write('\n')
                self.server.config=path
                return self.reply({'path':str(path),'message':'Draft saved. No analysis was started.'})
            if self.path in ('/api/check','/api/save','/api/run'):
                ui_state=data.get('ui_state',{})
                if self.path!='/api/check' and not isinstance(ui_state,dict):raise ValueError('UI state must be an object.')
                cfg=check(copy.deepcopy(data['config']))
                dia=is_dia_study(cfg)
                if dia and data.get('report_only') is True:raise ValueError('Report-only DDA statistics are unavailable for DIA-NN. Choose DIA input checking or analysis.')
                preview=validate_browser_study(cfg,ROOT) if dia else None
                if self.path=='/api/check': return self.reply({'message':f"Basic checks passed for {len(cfg['samples'])} samples. Run --validate_only for acquisition, checksum and tool validation.",**({'dia':preview} if dia else {})})
                target=ROOT/'ui_configs'; target.mkdir(exist_ok=True)
                path=target/('study_'+secrets.token_hex(6)+'.json')
                with path.open('x') as f: json.dump(cfg,f,indent=2); f.write('\n')
                with path.with_suffix('.ui.json').open('x') as handle:json.dump(ui_state,handle,indent=2);handle.write('\n')
                native=prepare_browser_study(path,ROOT) if dia else None
                if self.path=='/api/run':
                    options={'dia_config':native} if dia else {}
                    state=self.server.runner.start(path,resume=data.get('resume') is True,validate_only=data.get('validate_only') is True,report_only=data.get('report_only') is True,**options)
                    self.server.config=path
                    return self.reply({'path':str(path),'run':state,**({'dia_config':str(native),'dia':preview} if dia else {})})
                command=f'nextflow run main.nf --dia_config {shlex.quote(str(native))} -ansi-log false' if dia else f'nextflow run main.nf --config {shlex.quote(str(path))} --verbosity 4 -ansi-log false'
                self.server.config=path
                resume_flag=' -resume' if dia else ' --continue'
                return self.reply({'path':str(path),'command':command+(resume_flag if data.get('resume',True) is True else ''),'validate':command+' --validate_only',**({'dia_config':str(native),'dia':preview} if dia else {})})
            return self.reply({'error':'Not found'},404)
        except jsonschema.ValidationError as e: return self.reply({'error':'.'.join(map(str,e.absolute_path))+': '+e.message},400)
        except (ValueError,KeyError,OSError) as e: return self.reply({'error':str(e)},400)

def open_browser(url):
    try:
        # Desktop opener returns promptly and uses the user's chosen browser.
        subprocess.Popen(['xdg-open',url],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    except OSError:
        if not webbrowser.open(url): print('Open this URL in your browser: '+url,flush=True)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--port',type=int,default=8765)
    p.add_argument('--config',type=Path,default=ROOT/'config.local.json')
    p.add_argument('--open-browser',action='store_true')
    p.add_argument('--detach',action='store_true',help='Keep the setup server available after the launcher terminal closes')
    p.add_argument('--no-log-mirror',action='store_true',help=argparse.SUPPRESS)
    p.add_argument('--overrides',default='{}',help='Initial resource/output settings supplied by the Nextflow launcher')
    a=p.parse_args();url=f'http://127.0.0.1:{a.port}'
    if a.detach:
        try:
            url,log_path,process=start_server(ROOT,a.config,a.port,json.loads(a.overrides))
        except (OSError,ValueError,RuntimeError) as e:
            raise SystemExit(str(e))
        print(('Configuration editor ready: ' if process else 'Interface already running: ')+url,flush=True)
        print('Closing this terminal leaves the setup page and any analysis available. Server log: '+str(log_path),flush=True)
        if a.open_browser:open_browser(url)
        return follow(url)
    try:
        server=ThreadingHTTPServer(('127.0.0.1',a.port),Handler)
    except OSError:
        try:
            with urllib.request.urlopen(url+'/api/config',timeout=3) as r: existing=json.load(r)
            if not isinstance(existing.get('root'),str) or not existing['root'] or Path(existing['root']).resolve()!=ROOT: raise ValueError('Port belongs to another application')
        except Exception as e: raise SystemExit(f'Cannot use {url}: {e}. Choose another --ui_port.')
        print('Interface already running. Reopening '+url+' with its current settings.',flush=True)
        if a.open_browser:open_browser(url)
        return follow(url)
    server.config=a.config.resolve();server.runner=Runner(ROOT)
    restored=server.runner.status()
    if restored.get('config') and Path(restored['config']).is_file() and (restored['state'] in ('running','stopping') or server.config==ROOT/'config.local.json'):
        server.config=Path(restored['config'])
    overrides=json.loads(a.overrides);server.resource_overrides=overrides
    if overrides:
        cfg=absolute_config(server.config)
        if 'outdir' in overrides:cfg['outdir']=overrides['outdir']
        cfg.setdefault('resources',{'cpus':32,'memory_gib':128})
        for k in ['cpus','memory_gib']:
            if k in overrides:cfg['resources'][k]=overrides[k]
        folder=ROOT/'ui_configs';folder.mkdir(exist_ok=True)
        path=folder/('interactive_'+secrets.token_hex(6)+'.json')
        path.write_text(json.dumps(cfg,indent=2)+'\n');server.config=path
    print('Configuration editor: '+url,flush=True)
    print('No new analysis starts until you click Run analysis. Runs and their logs survive interface closure.',flush=True)
    print('Use Stop run in the browser to stop analysis. Reopen the interface to reconnect.',flush=True)
    stopped=threading.Event()
    mirror=None
    if not a.no_log_mirror:
        mirror=threading.Thread(target=mirror_runner,args=(server.runner,stopped),daemon=True);mirror.start()
    if a.open_browser:threading.Timer(.4,open_browser,args=(url,)).start()
    def interrupted(signum,frame):raise KeyboardInterrupt
    signal.signal(signal.SIGTERM,interrupted)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print('Interface closed. Any active analysis continues; reopen the interface to follow or stop it.',flush=True)
    finally:
        stopped.set()
        if mirror:mirror.join(timeout=2)
        server.server_close()
    return 0
if __name__=='__main__': raise SystemExit(main())
