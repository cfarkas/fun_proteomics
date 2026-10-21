"""Reconnect without discarding browser edits or replaying an uncertain action."""
import json
from pathlib import Path
import shutil
import subprocess
import unittest

ROOT = Path(__file__).resolve().parents[1]


@unittest.skipUnless(shutil.which('node'), 'Node unavailable for browser logic tests')
class BrowserConnectionTests(unittest.TestCase):
    def test_disconnect_reconnect_tokens_samples_and_initial_load(self):
        html = (ROOT / 'assets/config_ui.html').read_text()
        network = html[html.index('let connectionLost='):html.index('const countThresholds=')]
        poll = html[html.index('async function poll()'):html.index('async function launch()')]
        detect = html[html.index('async function detectSamples()'):html.index("action('detect_samples'")]
        resources = html[html.index('async function detectResources('):html.index('async function validateStep(')]
        initialize = html[html.index('async function initialize()'):html.index('\ninitialize();')]
        script = r'''
const assert=require('node:assert/strict'),vm=require('node:vm');
const fields=new Map(),messages=[];
function $(id){if(!fields.has(id))fields.set(id,{value:'',checked:false,disabled:false,hidden:false,dataset:{},style:{}});return fields.get(id);}
const sample={sample_id:'custom_sample',condition:'edited_condition',raw_path:'/data/sample.d'};
const original={samples:[sample]},calls=[];
let renders=0;
const context={$,AbortSignal,structuredClone,token:'old',window:{location:{port:'8766'}},
 setTimeout:()=>1,clearTimeout:()=>{},message:(text)=>messages.push(text),
 cfg:original,runBusy:false,launchPending:false,stopPending:false,maxStep:7,RUN_STEP:7,
 lastSerial:-1,lastState:'idle',restoredRun:true,pollPending:false,
 showStep:()=>{},rememberAssignments:()=>{},renderSamples:()=>{renders++;},
 inputRevision:0,sampleAssignments:new Map([[sample.raw_path,sample]]),
 updateDEPOptions:()=>{},refreshConditionOptions:()=>{},updateImmunoOptions:()=>{},loadOptionalOptions:()=>{},
 dbFields:()=>{},invalidateStep:()=>{},phosphoModule:null,
 readLocalDraft:()=>null,restoreDraftFields:()=>{},updateDiaControls:()=>{}};
vm.createContext(context);vm.runInContext(SOURCE,context);
const response=(x,status=200)=>({ok:status===200,status,json:async()=>x});
const status={state:'running',kind:'analysis',elapsed:40,pid:23,log_path:'/run.log',exit_code:null,lines:[],serial:1};
(async()=>{
 vm.runInContext('setupInitialized=true',context);
 $('cpus').value='77';$('outdir').value='/edited/output';
 context.fetch=async(path,options)=>{calls.push([path,options]);throw TypeError('NetworkError when attempting to fetch resource.');};
 await assert.rejects(context.api('run',{config:{}}),error=>error.code==='CONNECTION_LOST'&&!error.message.includes('NetworkError'));
 assert.equal(calls.length,1,'An uncertain Run request must never be replayed');
 assert.equal($('connection-status').hidden,false);
 assert.match($('connection-status').textContent,/--ui_port 8766/);
 await context.poll();assert.equal($('run').disabled,true);assert.equal($('stop').disabled,true);
 context.fetch=async(path,options)=>{calls.push([path,options]);return response(status);};
 await context.poll();
 assert.match($('connection-status').textContent,/Connection restored/);
 assert.equal($('cpus').value,'77');assert.equal($('outdir').value,'/edited/output');
 assert.equal(context.cfg,original);assert.equal(context.cfg.samples[0],sample);
 assert.equal($('stop').disabled,false);assert.equal($('run').disabled,true);
 assert.equal(calls.filter(x=>x[0]==='/api/run').length,1);

 calls.length=0;
 context.fetch=async(path,options)=>{
  calls.push([path,options]);
  if(path==='/api/session')return response({token:'new'});
  if(options.headers['X-Config-Token']==='old')return response({code:'TOKEN_STALE'},403);
  return response({state:'running'});
 };
 await context.api('run',{config:{}});
 assert.deepEqual(calls.map(x=>x[0]),['/api/run','/api/session','/api/run']);
 assert.equal(calls[2][1].headers['X-Config-Token'],'new');

 calls.length=0;
 context.fetch=async(path)=>{calls.push(path);return {ok:true,status:200,json:async()=>{throw TypeError('body disconnected');}};};
 await assert.rejects(context.api('stop',{}),error=>error.code==='CONNECTION_LOST');
 assert.equal(calls.length,1,'A lost response body must not replay Stop');
 await assert.rejects(context.detectSamples(),error=>error.code==='CONNECTION_LOST');
 assert.equal(context.cfg.samples[0],sample,'Failed detection must preserve edited sample rows');
 assert.equal(renders,0);

 vm.runInContext('setupInitialized=false',context);context.cfg=undefined;
 const config={samples:[sample],database:{},tools:{},outdir:'/initial',modes:['DDA'],resources:{cpus:8,memory_gib:32}};
 const routes={config:{config,root:'/project',config_path:'/project/study.json',resource_overrides:{cpus:8,memory_gib:32}},
  modules:{phosphoproteomics:{}},capabilities:{available_cpus:16,max_memory_gib:64,memory_available_gib:64,recommended_threads:16,recommended_memory_gib:56,supported:true}};
 context.fetch=async(path)=>{if(path==='/api/modules')throw TypeError('offline');return response(routes[path.split('/').pop()]);};
 await context.initialize();assert.equal(context.cfg,undefined,'Partial initialization must not overwrite fields');
 context.fetch=async(path)=>response(routes[path.split('/').pop()]);
 await context.initialize();assert.equal($('cpus').value,8);assert.equal($('wizard-next').disabled,false);
 $('cpus').value='11';await context.initialize();assert.equal($('cpus').value,'11','Reconnect must not reload initialized fields');
})().catch(error=>{console.error(error);process.exitCode=1;});
'''.replace('SOURCE', json.dumps(network + '\n' + poll + '\n' + detect + '\n' + resources + '\n' + initialize))
        result = subprocess.run(['node', '-'], input=script, text=True, capture_output=True, timeout=15)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
