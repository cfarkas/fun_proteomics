#!/usr/bin/env python3
"""Record real synthetic setup/analysis sessions; requires Playwright and ffmpeg.

Browser clips are continuous screen recordings. The terminal clip replays an
actual PTY transcript at a readable pace. Raw DIA demonstrates setup only.
No real acquisitions, credentials, patient paths or external data are used.
"""
import argparse
import errno
import hashlib
import json
import os
from pathlib import Path
import pty
import re
import select
import shutil
import signal
import sqlite3
import subprocess
import time

ROOT=Path(__file__).resolve().parents[1]
CAPTION_JS="""([title,text])=>{let e=document.getElementById('tutorial-caption');if(!e){e=document.createElement('div');e.id='tutorial-caption';e.style.cssText='position:fixed;z-index:100000;bottom:0;left:0;right:0;background:#123c45;color:white;padding:16px 30px;min-height:106px;font:23px/1.45 system-ui;box-shadow:0 -3px 12px #123c4522';document.body.style.paddingBottom='125px';document.body.append(e);}e.replaceChildren();const t=document.createElement('div');t.textContent=title;t.style.cssText='font-size:13px;color:#afe2cd;margin-bottom:6px;letter-spacing:.4px;text-transform:uppercase';const c=document.createElement('div');c.textContent=text;e.append(t,c);}"""
TERMINAL_HTML='''<!doctype html><meta charset="utf-8"><title>fun_proteomics terminal walkthrough</title><style>*{box-sizing:border-box}body{margin:0;background:#eef4f2;color:#173c45;font:18px system-ui}header{padding:24px 30px;background:#123c45;color:white;font-size:25px}header small{float:right;font-size:14px;color:#c8e7db;margin-top:8px}main{margin:24px;border-radius:12px;background:#102831;color:#e1eee9;padding:22px;height:770px;overflow:hidden}pre{margin:0;white-space:pre-wrap;overflow-wrap:anywhere;font:18px/26px 'DejaVu Sans Mono',monospace}</style><header>fun_proteomics · Terminal setup<small>Actual synthetic session · replayed at a readable pace</small></header><main id="screen"><pre id="out"></pre></main>'''


def timestamp(seconds):
    ms=round(seconds*1000);return f'{ms//3600000:02}:{ms//60000%60:02}:{ms//1000%60:02}.{ms%1000:03}'


class Clip:
    def __init__(self,page,name,title,output):
        self.page=page;self.name=name;self.title=title;self.output=output
        self.start=time.monotonic();self.captions=[]
    def say(self,text,seconds=5):
        elapsed=time.monotonic()-self.start
        if self.captions:self.captions[-1]['end']=elapsed
        self.captions.append(dict(start=elapsed,text=text))
        self.page.evaluate(CAPTION_JS,[self.title,text]);self.page.wait_for_timeout(seconds*1000)
    def finish(self,context):
        self.captions[-1]['end']=time.monotonic()-self.start
        self.page.screenshot(path=str(self.output/(self.name+'.png')))
        video=self.page.video;context.close();source=Path(video.path())
        target=self.output/(self.name+'.mp4')
        subprocess.run(['ffmpeg','-y','-loglevel','error','-i',str(source),'-an','-c:v','libx264','-preset','fast','-crf','23','-pix_fmt','yuv420p','-movflags','+faststart',str(target)],check=True)
        subprocess.run(['ffmpeg','-y','-loglevel','error','-ss','4','-i',str(target),'-frames:v','1',str(self.output/(self.name+'.jpg'))],check=True)
        vtt='WEBVTT\n\n'+''.join(f'{timestamp(c["start"])} --> {timestamp(c["end"])}\n{c["text"]}\n\n' for c in self.captions)
        (self.output/(self.name+'.vtt')).write_text(vtt)
        return {'file':target.name,'sha256':hashlib.sha256(target.read_bytes()).hexdigest(),
                'captions':self.captions,'bytes':target.stat().st_size}


def fixture(folder):
    inputs=folder/'inputs';inputs.mkdir(parents=True)
    for src,dst in [('quantification_matrix.csv','proteins.csv'),('quantification_samples.csv','samples.csv'),('service_quantification.csv','service.csv')]:
        shutil.copy2(ROOT/'examples'/src,inputs/dst)
    return inputs


def launch(args,folder,env):
    process=subprocess.Popen([args.python,'-m','fun_proteomics','setup','--project',str(folder/'study'),'--no-browser'],
                              cwd=folder,env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,start_new_session=True)
    url=None;lines=[];deadline=time.monotonic()+30
    while time.monotonic()<deadline:
        ready,_,_=select.select([process.stdout],[],[],1)
        if ready:
            data=os.read(process.stdout.fileno(),65536).decode(errors='replace');lines.append(data)
            match=re.search(r'http://127\.0\.0\.1:\d+/#[^\s]+',''.join(lines))
            if match:url=match.group(0);break
            if not data and process.poll() is not None:break
    if not url:
        process.terminate();raise RuntimeError('Setup did not open: '+''.join(lines))
    return process,url


def browser_clip(browser,args,env,raw=False):
    name='raw-dia-setup' if raw else 'browser-analysis';folder=args.work_dir/name;folder.mkdir()
    inputs=fixture(folder);process,url=launch(args,folder,env)
    context=None
    try:
        context=browser.new_context(viewport={'width':1400,'height':1040},record_video_dir=str(args.output/'capture'),record_video_size={'width':1400,'height':1040})
        page=context.new_page();errors=[];page.on('pageerror',lambda e:errors.append(str(e)))
        page.goto(url);page.wait_for_function('document.getElementById("project").value.length>0')
        clip=Clip(page,name,'DIA raw setup · synthetic metadata · no spectrum search' if raw else 'Browser analysis · synthetic inputs · real analysis',args.output)
        clip.say('Launch with: fun_proteomics setup --project /path/to/my-study. One folder holds settings and results.',6)
        if raw:
            acquisitions=folder/'raw';acquisitions.mkdir()
            for i in range(4):
                path=acquisitions/f'Sample{i+1}.d';path.mkdir();(path/'analysis.tdf_bin').write_bytes(b'SYNTHETIC METADATA DEMO - NO SPECTRA')
                with sqlite3.connect(path/'analysis.tdf') as db:
                    db.execute('CREATE TABLE Frames (MsMsType INTEGER)');db.executemany('INSERT INTO Frames VALUES (?)',[(0,),(9,)])
            fasta=inputs/'targets.fasta';fasta.write_text('>sp|P00001|DEMO GN=GENE_A\nMPEPTIDEK\n')
            page.locator('input[name=mode][value=dia]').check();page.locator('#raw_root').fill(str(acquisitions))
            page.locator('#scan').click();page.wait_for_function('document.querySelectorAll("#raw-table tbody tr").length===4')
            clip.say('Choose DIA-PASEF and scan your raw-data folder. This demo contains synthetic acquisition metadata only.',6)
            for i,row in enumerate(page.locator('#raw-table tbody tr').all()):row.locator('input').nth(1).fill('Control' if i<2 else 'Treatment')
            clip.say('Assign a condition to every sample. Conditions are never guessed from filenames.',6)
            page.locator('#next').click();page.locator('#independent_samples').check();page.locator('#fasta').fill(str(fasta))
            clip.say('Confirm independent biological replicates. Choose the matching target-only FASTA with gene annotations.',7)
            page.locator('#check').click();page.wait_for_function('document.getElementById("preview").textContent.includes("1 condition pairs")')
            page.locator('#preview').scroll_into_view_if_needed()
            clip.say('Check inputs before running. Two conditions give one comparison at protein and available gene levels.',7)
            page.locator('#save').click();page.wait_for_function('!document.getElementById("run").disabled')
            clip.say('The configuration is saved. With real spectra and installed tools, Run / continue starts the Nextflow workflow.',7)
            clip.say('This recording stops after setup. No spectrum search or biological result is claimed for synthetic acquisitions.',6)
            assert not (folder/'study/results').exists()
        else:
            page.locator('[data-browse=matrix]').click();page.locator('#browse-path').fill(str(inputs));page.locator('#browse-go').click()
            page.locator('#entries').get_by_role('button',name='proteins.csv',exact=True).wait_for()
            clip.say('Browse files on the computer running the analysis. Select your protein quantification table.',5)
            page.locator('#entries').get_by_role('button',name='proteins.csv',exact=True).click();page.locator('#samples').fill(str(inputs/'samples.csv'))
            clip.say('Select the sample table. Sample IDs match quantity columns; condition names define the comparison groups.',6)
            page.locator('#next').click();page.locator('#independent_samples').check();page.locator('#service_matrix').fill(str(inputs/'service.csv'))
            clip.say('Keep all-versus-all comparisons and confirm independent samples. An optional service matrix adds detection gains and losses.',8)
            page.locator('#check').click();page.wait_for_function('document.getElementById("preview").textContent.includes("3 condition pairs")')
            page.locator('#preview').scroll_into_view_if_needed()
            clip.say('The check finds six samples, three condition pairs, ten protein groups and ten unambiguous genes.',6)
            page.locator('#save').click();page.wait_for_function('!document.getElementById("run").disabled')
            clip.say('Save configuration and check writes config/run.json. Now start the saved analysis.',5)
            page.locator('#run').click();page.wait_for_function('document.getElementById("state").textContent==="completed"',timeout=120000)
            clip.say('The real synthetic analysis is complete. Open results for comparisons, QC and the service benchmark.',6)
            root=url.split('/#')[0];page.goto(root+'/results/MATRICES.html');page.wait_for_load_state()
            clip.say('The matrix catalogue links directly to CSV and TSV files for proteins, genes and transformed quantities.',7)
            assert 'CSV' in page.locator('body').inner_text() and 'TSV' in page.locator('body').inner_text()
            page.goto(url);page.wait_for_function('!document.getElementById("run").disabled');page.locator('#run').click()
            page.wait_for_function('document.getElementById("log").textContent.includes("No analysis needs repeating")',timeout=60000)
            clip.say('Run / continue verifies completed outputs. Unchanged results are reused without repeating the analysis.',8)
            assert len(list((folder/'study/results/final_results/differential_expression').glob('*/*/all_results.csv')))==6
        assert not errors,errors
        return clip.finish(context)
    finally:
        if context:context.close()
        process.send_signal(signal.SIGINT)
        try:process.wait(timeout=10)
        except subprocess.TimeoutExpired:process.terminate();process.wait(timeout=5)
        process.stdout.close()


def terminal_capture(args,folder,env):
    inputs=fixture(folder);events=[];origin=time.monotonic()
    def event(text):events.append([round(time.monotonic()-origin,3),'o',text])
    prompts=[('Analysis [1]:','1'),('Protein matrix (CSV, TSV or XLSX):','inputs/proteins.csv'),
             ('Sample table (sample_id, condition; raw_path for acquisitions):','inputs/samples.csv'),
             ('Matching FASTA for missing gene symbols (optional):',''),
             ('Service matrix for a detection benchmark (optional):','inputs/service.csv'),
             ('Compare every condition against every other condition? (yes/no) [yes]:','yes'),
             ('Does each sample come from a different patient or biological replicate? (yes/no) [no]:','yes')]
    for command,shown in [(['setup','--terminal','--project','./study','--run'],'setup --terminal --project ./study --run'),(['--continue','./study'],'--continue ./study')]:
        event('$ python -m fun_proteomics '+shown+'\r\n')
        master,slave=pty.openpty();process=subprocess.Popen([args.python,'-m','fun_proteomics',*command],cwd=folder,env=env,stdin=slave,stdout=slave,stderr=slave,start_new_session=True)
        os.close(slave);buffer='';index=0;expected=prompts if command[0]=='setup' else [];deadline=time.monotonic()+120
        try:
            while time.monotonic()<deadline:
                ready,_,_=select.select([master],[],[],.2)
                if ready:
                    try:data=os.read(master,65536)
                    except OSError as error:
                        if error.errno==errno.EIO:break
                        raise
                    if not data:break
                    text=data.decode(errors='replace');event(text);buffer+=text
                    if index<len(expected) and expected[index][0] in buffer:
                        time.sleep(.5);os.write(master,(expected[index][1]+'\n').encode());buffer='';index+=1
                elif process.poll() is not None:break
            else:raise RuntimeError('Terminal recording timed out at '+buffer[-600:])
            if process.wait(timeout=10)!=0 or index!=len(expected):raise RuntimeError('Terminal setup failed: '+buffer[-1500:])
        finally:
            if process.poll() is None:process.terminate();process.wait(timeout=5)
            os.close(master)
    assert 'No analysis needs repeating' in ''.join(e[2] for e in events)
    header={'version':2,'width':118,'height':30,'title':'fun_proteomics synthetic analysis and continuation','env':{'TERM':'dumb'}}
    (args.output/'terminal-analysis.cast').write_text(json.dumps(header)+'\n'+''.join(json.dumps(e)+'\n' for e in events))
    transcript=''.join(e[2] for e in events).replace('\r\n','\n')
    (args.output/'terminal-analysis.txt').write_text(transcript)
    return events


def terminal_clip(browser,args,env):
    folder=args.work_dir/'terminal-analysis';folder.mkdir();events=terminal_capture(args,folder,env)
    context=browser.new_context(viewport={'width':1400,'height':1040},record_video_dir=str(args.output/'capture'),record_video_size={'width':1400,'height':1040})
    page=context.new_page();page.set_content(TERMINAL_HTML)
    clip=Clip(page,'terminal-analysis','Terminal analysis · actual session replay · synthetic inputs',args.output)
    clip.say('For SSH or terminal users: add --terminal and --run to the same setup command.',6)
    stage=''
    for _,_,text in events:
        page.evaluate("text=>{document.getElementById('out').appendChild(document.createTextNode(text));const s=document.getElementById('screen');s.scrollTop=s.scrollHeight;}",text.replace('\r\n','\n'))
        if 'Analysis [1]' in text:description='Choose a quantification matrix, raw DIA, raw DDA or an existing run. Enter accepts the displayed default.'
        elif 'Sample table' in text and 'samples.csv' not in text:description='Provide the protein matrix and sample table. The same inputs work in browser setup.'
        elif 'Matching FASTA' in text:description='A matching FASTA can supply missing gene symbols. This synthetic table already includes gene names.'
        elif 'Service matrix' in text:description='Optionally add service quantification to benchmark shared proteins and detection gains or losses.'
        elif 'Does each sample' in text:description='Confirm independence only when samples represent different patients or biological replicates.'
        elif 'Configuration:' in text or 'Saved:' in text:description='Setup saves the configuration, reports the planned comparisons, and starts the real analysis.'
        elif '$ python' in text and '--continue' in text:description='Later, from any working directory, continue by supplying the same project folder.'
        elif 'No analysis needs repeating' in text:description='Checksums and recorded inputs match. Completed results are reused without recomputation.'
        else:description=''
        if description and description!=stage:stage=description;clip.say(description,5)
        else:page.wait_for_timeout(750)
    clip.say('Open START_HERE.html for results and MATRICES.html for clear CSV and TSV links.',7)
    return clip.finish(context)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--python',required=True,help='Python with the built package installed')
    parser.add_argument('--work-dir',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--browser-path',help='Optional existing Chromium executable')
    parser.add_argument('--only',choices=['browser','terminal','raw'])
    args=parser.parse_args();args.work_dir=args.work_dir.resolve();args.output=args.output.resolve()
    args.work_dir.mkdir(parents=True,exist_ok=True);args.output.mkdir(parents=True,exist_ok=True)
    env=os.environ.copy();env.pop('PYTHONPATH',None)
    env.update(FUN_PROTEOMICS_CACHE=str(args.work_dir/'workflow-cache'),TERM='dumb',NO_COLOR='1',PYTHONUNBUFFERED='1')
    version=subprocess.check_output([args.python,'-m','fun_proteomics','--version'],env=env,text=True).strip()
    from playwright.sync_api import sync_playwright
    output=[]
    with sync_playwright() as p:
        browser=p.chromium.launch(headless=True,**({'executable_path':args.browser_path} if args.browser_path else {}))
        try:
            if args.only in (None,'browser'):output.append(browser_clip(browser,args,env))
            if args.only in (None,'terminal'):output.append(terminal_clip(browser,args,env))
            if args.only in (None,'raw'):output.append(browser_clip(browser,args,env,raw=True))
        finally:browser.close()
    metadata={'version':version,'input_kind':'synthetic only','recorder_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'clips':output}
    (args.output/('recording-'+(args.only or 'all')+'.json')).write_text(json.dumps(metadata,indent=2)+'\n')
    for clip in output:print('Recorded '+clip['file']+' ('+str(clip['bytes'])+' bytes)')


if __name__=='__main__':main()
