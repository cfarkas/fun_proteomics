"""Inactive thresholds must not block saved configurations or browser launches."""
import json
from pathlib import Path
import shutil
import subprocess
import unittest

import jsonschema

ROOT = Path(__file__).resolve().parents[1]
THRESHOLDS = ['min_count', 'min_total_count', 'modification_min_count', 'modification_min_total_count']


class CountFilterSchemaTests(unittest.TestCase):
    def setUp(self):
        schema = json.loads((ROOT / 'assets/config.schema.json').read_text())
        jsonschema.Draft202012Validator.check_schema(schema)
        self.validator = jsonschema.Draft202012Validator(schema['properties']['differential'])

    def test_off_or_omitted_filter_ignores_zero_thresholds(self):
        for setting in [{}, {'filter_low_counts': False}, {'enabled': False, 'filter_low_counts': True}]:
            self.validator.validate({'enabled': True, 'modifications': True,
                                     **dict.fromkeys(THRESHOLDS, 0), **setting})

    def test_enabled_filter_rejects_zero_only_for_active_thresholds(self):
        base = {'enabled': True, 'filter_low_counts': True, 'modifications': True,
                **dict.fromkeys(THRESHOLDS, 1)}
        self.validator.validate(base)
        for key in THRESHOLDS:
            with self.subTest(key=key), self.assertRaises(jsonschema.ValidationError):
                self.validator.validate({**base, key: 0})
        self.validator.validate({**base, 'modifications': False,
                                 'modification_min_count': 0, 'modification_min_total_count': 0})


@unittest.skipUnless(shutil.which('node'), 'Node unavailable for browser logic tests')
class BrowserControlTests(unittest.TestCase):
    def test_count_inputs_receive_desktop_focus_without_overwriting_edits(self):
        html = (ROOT / 'assets/config_ui.html').read_text()
        focus = html[html.index('// Request native desktop focus'):html.index('function message(')]
        harness = r'''
const assert=require('node:assert/strict'),vm=require('node:vm');
let pointerdown,requests=0,release;
const document={activeElement:null,addEventListener:(name,handler)=>{assert.equal(name,'pointerdown');pointerdown=handler;}};
const context={document,window:{focus:()=>{}},api:async name=>{assert.equal(name,'focus_browser');requests++;}};
vm.createContext(context);vm.runInContext(SOURCE,context);
const field=(value,step='1')=>({value,step,disabled:false,readOnly:false,focusCount:0,
 setAttribute(name,value){this[name]=value;},focus(){this.focusCount++;document.activeElement=this;}});
const click=input=>pointerdown({target:{closest:selector=>{assert.equal(selector,'input[type="number"]');return input;}}});
(async()=>{
 for(const id of ['de-min-count','de-min-total','de-mod-min-count','de-mod-min-total']){
  const input=field('7');await click(input);
  assert.equal(document.activeElement,input,id);assert.equal(input.value,'7',id);
  assert.equal(input.focusCount,2,id);assert.equal(input.inputmode,'numeric');
 }
 const fraction=field('0.05','0.01');await click(fraction);assert.equal(fraction.inputmode,'decimal');
 const previousRequests=requests,disabled=field('15');disabled.disabled=true;await click(disabled);await click(null);
 assert.equal(requests,previousRequests);assert.equal(disabled.focusCount,0);
 context.api=()=>new Promise(resolve=>{release=resolve;});
 const total=field('11'),pending=click(total),other=field('3');other.focus();release();await pending;
 assert.equal(document.activeElement,other,'A completed focus request must not steal focus back');
 assert.equal(total.value,'11');assert.equal(total.focusCount,1);
 context.api=async()=>{throw Error('No desktop available');};
 await click(total);assert.equal(document.activeElement,total);assert.equal(total.value,'11');
})().catch(error=>{console.error(error);process.exitCode=1;});
'''.replace('SOURCE', json.dumps(focus))
        result = subprocess.run(['node', '-'], input=harness, text=True, capture_output=True, timeout=15)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_uncheck_resets_invalid_fields_omits_thresholds_and_stop_tracks_state(self):
        html = (ROOT / 'assets/config_ui.html').read_text()
        collect = html[html.index('const countThresholds='):html.index('function renderSamples()')]
        update = html[html.index('function updateDEPOptions()'):html.index("\n$('de-enabled').onchange")]
        poll = html[html.index('async function poll()'):html.index('async function launch()')]
        stop = html[html.index('let stopPending=false;'):html.index(';setInterval(poll,1500)')]
        harness = r'''
const assert=require('node:assert/strict'), vm=require('node:vm');
const elements=new Map(), messages=[];
function $(id){if(!elements.has(id))elements.set(id,{value:'',checked:false,disabled:false,hidden:false,dataset:{},style:{}});return elements.get(id);}
const context={$,Number,structuredClone,cfg:{database:{},samples:[],tools:{}},conditionLabels:{},phosphoModule:null,
 effectiveOutputFolder:()=>'/output',fastaDownloadReady:false,runBusy:false,lastSerial:-1,lastState:'idle',
 launchPending:false,maxStep:7,RUN_STEP:7,restoredRun:true,pollPending:false,setupInitialized:true,showStep:()=>{},message:(text)=>messages.push(text)};
vm.createContext(context);
vm.runInContext(SOURCE,context);
$('cpus').value='8';$('memory').value='32';$('modes').value='both';$('fasta-source').value='local';
$('de-enabled').checked=true;$('de-edgeR').checked=true;$('de-control').value='CTRL';$('de-modifications').checked=true;
$('de-fdr').value='0.05';$('de-upper').value='100000';
const ids=['de-min-count','de-min-total','de-mod-min-count','de-mod-min-total'];
const keys=['min_count','min_total_count','modification_min_count','modification_min_total_count'];
for(const id of ids)$(id).value='0';
$('de-filter-counts').checked=false;context.updateDEPOptions();
assert.ok(ids.every(id=>$(id).disabled&&Number($(id).value)>=1));
assert.equal($('de-filter-options').hidden,true);
// Even stale values introduced without a change event must not be submitted.
for(const id of ids)$(id).value='0';
const off=context.collect().differential;
assert.equal(off.filter_low_counts,false);assert.ok(keys.every(key=>!(key in off)));
$('de-filter-counts').checked=true;context.updateDEPOptions();
assert.ok(ids.every(id=>!$(id).disabled));assert.throws(()=>context.collect(),/minimum counts must be at least 1/);
$('de-filter-counts').checked=false;context.updateDEPOptions();
$('de-filter-counts').checked=true;context.updateDEPOptions();
const on=context.collect().differential;assert.equal(on.filter_low_counts,true);assert.equal(on.min_count,10);
$('de-modifications').checked=false;context.updateDEPOptions();
assert.equal($('de-mod-min-count').disabled,true);assert.equal('modification_min_count' in context.collect().differential,false);
(async()=>{
 for(const state of ['idle','running','stopping','completed','failed','stopped']){
  const status={state,kind:'analysis',elapsed:0,pid:123,log_path:'/run.log',exit_code:state==='completed'?0:null,lines:[],serial:1};
  context.getApi=async()=>status;await context.poll();
  assert.equal($('stop').disabled,state!=='running',state);
 }
 context.api=async()=>({state:'completed',message:'Analysis already completed. There is no active run to stop.'});
 await $('stop').onclick();assert.match(messages.at(-1),/already completed/);
 assert.ok(!messages.at(-1).includes('Stop requested'));
})().catch(error=>{console.error(error);process.exitCode=1;});
'''.replace('SOURCE', json.dumps(collect + '\n' + update + '\n' + poll + '\n' + stop))
        result = subprocess.run(['node', '-'], input=harness, text=True, capture_output=True, timeout=15)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
