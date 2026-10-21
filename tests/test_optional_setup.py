"""Saved drug designs survive sample refreshes, and setup errors identify their controls."""
import json
from pathlib import Path
import shutil
import subprocess
import unittest

ROOT = Path(__file__).resolve().parents[1]


@unittest.skipUnless(shutil.which('node'), 'Node unavailable for browser logic tests')
class OptionalSetupTests(unittest.TestCase):
    def test_saved_design_refresh_validation_navigation_and_continue_request(self):
        html = (ROOT/'assets/config_ui.html').read_text()
        source = '\n'.join([
            html[html.index('function showConfigurationError('):html.index('let connectionLost=')],
            html[html.index('const countThresholds='):html.index('function renderSamples()')],
            html[html.index('function refreshDrugConditions()'):html.index('let previousSearchModes=')],
            html[html.index('async function launch()'):html.index("action('run',launch)")],
        ])
        harness = r'''
const assert=require('node:assert/strict'),vm=require('node:vm');
const fields=new Map(),selectIds=new Set(['interaction-cohort','interaction-control','interaction-drug-a','interaction-drug-b','interaction-combination','optional-stats-control']);
function $(id){
 if(!fields.has(id)){
  const field={value:'',checked:false,disabled:false,hidden:false,options:[],style:{},attributes:{},
   replaceChildren(...options){this.options=options;this.value=options[0]?.value||'';},add(option){this.options.push(option);},
   setAttribute(key,value){this.attributes[key]=value;},removeAttribute(key){delete this.attributes[key];},
   classList:{toggle(){}},focus(){context.focused=id;},scrollIntoView(){context.scrolled=id;}};
  if(selectIds.has(id))Object.defineProperty(field,'value',{get(){return this._value||'';},set(value){this._value=this.options.some(o=>o.value===value)?value:'';}});
  fields.set(id,field);
 }
 return fields.get(id);
}
const samples=['Control','A','B','AB'].flatMap(condition=>[1,2].map(replicate=>({sample_id:condition+replicate,cohort:'cells',condition,replicate})));
const saved={cohort:'cells',control:'Control',drug_a:'A',drug_b:'B',combination:'AB',enabled:true,matched_design:false,methods:['edgeR'],fdr:.05};
const context={draftSnapshot:()=>({ui_state:{}}),rememberDraft:()=>{},$,Number,structuredClone,cfg:{samples,optional_analyses:{peptide_isoforms:true,
 drugs:{applied:true,study_type:'synergy',interaction:saved},statistics:{enabled:true,features:['peptide'],methods:['edgeR'],comparisons:'control',control:'Control',fdr:.05}}},
 Option:function(text,value){this.text=text;this.value=value;},invalidateStep(){},
 showStep(step){context.step=step;},message(text){context.messageText=text;},
 document:{querySelectorAll:()=>[...fields.values()].filter(f=>f.attributes['aria-invalid'])},
 launchPending:false,runBusy:false,maxStep:7,RUN_STEP:7,lastSerial:0,poll:async()=>{}};
vm.createContext(context);vm.runInContext(SOURCE,context);
context.loadOptionalOptions();
function assignments(){const x=context.collectOptionalOptions().drugs.interaction;return ['cohort','control','drug_a','drug_b','combination'].map(key=>x[key]);}
const expected=['cells','Control','A','B','AB'];
assert.deepEqual(assignments(),expected);
// During interactive rediscovery there may be no selected samples yet.
context.cfg.samples=[];context.refreshDrugConditions();
assert.equal($('interaction-control').value,'');
context.cfg.samples=samples;context.refreshDrugConditions();
assert.deepEqual(assignments(),expected);
assert.equal($('optional-stats-control').value,'Control');
assert.match($('statistics-drug-conditions').textContent,/Drug A alone: A/);
// A deliberate cleared selection stays cleared and cannot be replaced by the statistics reference.
$('interaction-drug-a').value='';$('interaction-drug-a').onchange();context.refreshDrugConditions();
assert.equal($('interaction-drug-a').value,'');
let missing;try{context.collectOptionalOptions();}catch(e){missing=e;}
assert.equal(missing.field,'interaction-drug-a');assert.equal(missing.optionalTab,'drugs');
assert.match(missing.message,/Drug A alone/);assert.match(missing.message,/reference in Statistics/);
context.showConfigurationError(missing);
assert.equal(context.step,5);assert.equal($('optional-panel-drugs').hidden,false);
assert.equal(context.focused,'interaction-drug-a');assert.equal($('interaction-validation').hidden,false);
$('interaction-drug-a').value='Control';$('interaction-drug-a').onchange();
assert.throws(()=>context.collectOptionalOptions(),e=>e.field==='interaction-drug-a'&&/both use/.test(e.message));
$('interaction-drug-a').value='A';$('interaction-drug-a').onchange();assert.deepEqual(assignments(),expected);
// A different cohort does not inherit the first cohort's drug assignments.
context.cfg.samples=samples.concat(samples.map(s=>({...s,cohort:'other'})));context.refreshDrugConditions();
$('interaction-cohort').value='other';$('interaction-cohort').onchange();
assert.equal($('interaction-control').value,'');
$('interaction-cohort').value='cells';$('interaction-cohort').onchange();assert.deepEqual(assignments(),expected);
// Removed conditions remain invalid; saved settings never invent evidence.
context.cfg.samples=samples.filter(s=>s.condition!=='AB');context.refreshDrugConditions();
assert.throws(()=>context.collectOptionalOptions(),e=>e.field==='interaction-combination');
context.cfg.samples=samples;context.refreshDrugConditions();assert.deepEqual(assignments(),expected);
$('drug-study-type').value='single';context.updateDrugOptions();
assert.equal('interaction' in context.collectOptionalOptions().drugs,false);
$('drugs-applied').value='no';assert.equal(context.collectOptionalOptions().drugs.applied,false);
$('drugs-applied').value='yes';$('drug-study-type').value='synergy';context.updateDrugOptions();
// Run uses exactly the continuation choice shown in the form.
const sent=[];context.api=async(path,data)=>sent.push({path,data});context.collect=()=>({optional_analyses:context.collectOptionalOptions()});
(async()=>{
 for(const resume of [false,true]){$('continue-run').checked=resume;await context.launch();assert.equal(sent.at(-1).data.resume,resume);}
 assert.equal(sent.length,2);
 $('interaction-control').value='';$('interaction-control').onchange();
 await assert.rejects(context.launch(),e=>e.field==='interaction-control');
 assert.equal(sent.length,2,'Invalid designs must not start a run');assert.equal(context.launchPending,false);
})().catch(error=>{console.error(error);process.exitCode=1;});
'''.replace('SOURCE', json.dumps(source))
        result = subprocess.run(['node', '-'], input=harness, text=True, capture_output=True, timeout=15)
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
