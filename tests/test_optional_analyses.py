"""Optional analysis routing, sequence evidence and quantitative interaction regressions."""
import copy
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import jsonschema

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'pipeline_lib'))
import optional_analyses as optional
import pipeline


def design():
    return {'applied': True, 'study_type': 'synergy', 'interaction': {
        'enabled': True, 'cohort': 'cells', 'control': 'Control', 'drug_a': 'A', 'drug_b': 'B',
        'combination': 'AB', 'matched_design': True, 'fdr': 0.05}}


def samples():
    return [{'sample_id': f'S{g}{i}', 'cohort': 'cells', 'condition': g,
             'replicate': i, 'raw_path': '/example/'+f'S{g}{i}.d'} for g in ('Control', 'A', 'B', 'AB') for i in range(1, 5)]


class OptionalTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory(); self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.plan = {'samples': samples(), 'optional_analyses': {'drugs': design()}}

    def test_design_gate_and_independent_replicates(self):
        optional.validate_options(self.plan)
        for edit in [lambda d: d.update(applied=False), lambda d: d.update(study_type='single'),
                     lambda d: d['interaction'].update(control='A'),
                     lambda d: d['interaction'].update(cohort='missing'),
                     lambda d: d['interaction'].update(matched_design='yes')]:
            plan = copy.deepcopy(self.plan); edit(plan['optional_analyses']['drugs'])
            with self.assertRaises(ValueError): optional.validate_options(plan)
        self.plan['samples'][1]['replicate'] = 1
        with self.assertRaisesRegex(ValueError, 'replicates'): optional.validate_options(self.plan)

    def test_unconfirmed_or_omitted_dose_matching_is_valid(self):
        schema = json.loads((ROOT/'assets/config.schema.json').read_text())['properties']['optional_analyses']
        validator = jsonschema.Draft202012Validator(schema)
        for matched in (False, None):
            with self.subTest(matched=matched):
                plan = copy.deepcopy(self.plan)
                interaction = plan['optional_analyses']['drugs']['interaction']
                if matched is None:
                    del interaction['matched_design']
                else:
                    interaction['matched_design'] = matched
                validator.validate(plan['optional_analyses'])
                optional.validate_options(plan)

    def test_schema_rejects_hidden_stale_drug_settings(self):
        schema = json.loads((ROOT/'assets/config.schema.json').read_text())['properties']['optional_analyses']
        validator = jsonschema.Draft202012Validator(schema)
        validator.validate({'drugs': design()})
        validator.validate({'drugs': {'applied': False}})
        validator.validate({'drugs': {'applied': True, 'study_type': 'single'}})
        for drugs in [{'applied': False, 'study_type': 'synergy'}, {'applied': True}, {**design(), 'study_type': 'single'}]:
            self.assertFalse(validator.is_valid({'drugs': drugs}))
        with self.assertRaisesRegex(ValueError, 'HLA-enriched'):
            optional.validate_options({'optional_analyses': {'hla_enriched': True}})

    def test_retired_discovery_settings_migrate_without_running_a_search(self):
        plan = copy.deepcopy(self.plan)
        plan['optional_analyses'].update(open_search=True, peptide_isoforms=True,
            statistics={'enabled': True, 'features': ['mass_shift', 'peptide'],
                        'methods': ['edgeR'], 'comparisons': 'control', 'control': 'Control'})
        optional.validate_options(plan)
        self.assertNotIn('open_search', plan['optional_analyses'])
        self.assertEqual(plan['optional_analyses']['statistics']['features'], ['peptide'])
        plan['optional_analyses']['statistics']['features'] = ['mass_shift']
        optional.validate_options(plan)
        self.assertFalse(plan['optional_analyses']['statistics']['enabled'])
        self.assertFalse(hasattr(pipeline, 'run_open_search'))
        self.assertFalse(hasattr(optional, 'export_open_search'))
        schema = json.loads((ROOT/'assets/config.schema.json').read_text())['properties']['optional_analyses']
        jsonschema.Draft202012Validator(schema).validate(plan['optional_analyses'])

    def export_sequences(self, fasta, peptides):
        database = self.root/'database.fasta';database.write_text(fasta)
        self.plan['database'] = {'fasta': str(database)}
        ids = [s['sample_id'] for s in self.plan['samples']]
        optional.write_csv(self.root/'peptide_spectral_counts.csv',['Peptide']+ids,
                           [{'Peptide': p, **dict.fromkeys(ids, 2)} for p in peptides])
        out = self.root/'isoforms'; optional.export_isoforms(self.plan,self.root,'cells',out)
        return out

    def test_isoform_specificity_shared_sequences_il_and_decoys(self):
        fasta = '>sp|P1|ONE GN=G1\nMPEPTIDEAAAKSHAREDK\n>sp|P1-2|TWO GN=G1\nMOTHERPEPKSHAREDK\n>sp|P2|THREE GN=G2\nMILKPEPTIDE\n>sp|P3|FOUR GN=G3\nMLLKPEPTIDE\n>rev_sp|P1|ONE\nMPEPTIDEAAAK\n'
        out = self.export_sequences(fasta,['PEPTIDEAAA','SHAREDK','MILKPEPTIDE','NOTFOUND'])
        rows = {r['Peptide']: r for r in optional.read_csv(out/'peptide_evidence.csv')}
        self.assertEqual(rows['PEPTIDEAAA']['Isoform_discriminating'],'True')
        self.assertEqual(rows['SHAREDK']['Matching_proteins'],'2')
        self.assertEqual(rows['MILKPEPTIDE']['Matching_proteins'],'2')
        self.assertEqual(rows['NOTFOUND']['Mapping_status'],'unmapped')
        counts = {r['Accession']: r for r in optional.read_csv(out/'isoform_evidence.csv')}
        self.assertEqual(counts['P1']['SControl1'],'2')
        self.assertEqual(counts['P1-2']['SControl1'],'0')
        self.assertEqual(len(optional.read_csv(out/'peptide_condition_summary.csv')),16)

    def test_canonical_only_and_empty_coverage_are_valid_reports(self):
        out = self.export_sequences('>P1\nPEPTIDE\n',['PEPTIDE'])
        self.assertEqual(optional.read_csv(out/'peptide_evidence.csv')[0]['Isoform_discriminating'],'False')
        out = self.export_sequences('>P1\nPEPTIDE\n',[])
        self.assertEqual(optional.read_csv(out/'isoform_evidence.csv'),[])
        self.assertEqual(json.loads((out/'settings.json').read_text())['peptides'],0)

    def test_zero_spectral_counts_do_not_claim_isoform_detection(self):
        fasta='>P1\nPEPTIDE\n>P1-2\nANOTHER\n'
        self.export_sequences(fasta,['PEPTIDE'])
        ids=[s['sample_id'] for s in self.plan['samples']]
        optional.write_csv(self.root/'peptide_spectral_counts.csv',['Peptide']+ids,
                           [{'Peptide':'PEPTIDE',**dict.fromkeys(ids,0)}])
        out=self.root/'zero';optional.export_isoforms(self.plan,self.root,'cells',out)
        self.assertEqual(optional.read_csv(out/'isoform_evidence.csv'),[])
        self.assertEqual(optional.read_csv(out/'peptide_evidence.csv')[0]['Isoform_discriminating'],'False')


class InteractionModelTests(unittest.TestCase):
    def run_model(self, hla=False, zero=False, matched=True, methods=None, monte_carlo=None,
                  feature_count=350, extra_zero=False, empty_library=False, gene_counts=None):
        tmp = tempfile.TemporaryDirectory();self.addCleanup(tmp.cleanup);root=Path(tmp.name)
        group = root/'inputs';group.mkdir()
        sample_rows = samples();ids = [s['sample_id'] for s in sample_rows]
        optional.write_csv(group/'samples.csv',['sample_id','cohort','condition','replicate'],
                           [{k:s[k] for k in ('sample_id','cohort','condition','replicate')} for s in sample_rows])
        import random
        rng = random.Random(32)
        rows = []
        for j in range(feature_count):
            values = {}
            for sample in sample_rows:
                factor = {'Control':1,'A':2,'B':3,'AB':6}[sample['condition']] if j<3 else 1
                if j==1 and sample['condition']=='AB':factor*=4
                if j==2 and sample['condition']=='AB':factor/=4
                values[sample['sample_id']] = 0 if zero or (empty_library and sample['sample_id']==ids[0]) else round(60*factor*rng.uniform(.9,1.1))
            rows.append({'FeatureID':f'P{j}','Gene':f'G{j}',**values})
        if extra_zero:
            rows.append({'FeatureID':'AllZero','Gene':'AllZero',**dict.fromkeys(ids,0)})
        for kind in ('protein','gene','peptide'):
            optional.write_csv(group/(kind+'_counts.csv'),['FeatureID','Gene']+ids,rows)
        job={'groups':[{'cohort':'cells','mode':'DDA','directory':str(group)}], 'outdir':str(root/'results'),
             'resources':{'cpus':2},
             'differential':{'enabled':False},'optional_analyses':{'drugs':design()},
             'analysis':'immunopeptidomics' if hla else 'global-label-free'}
        if gene_counts is not None: job['optional_analyses']['gene_counts'] = gene_counts
        interaction = job['optional_analyses']['drugs']['interaction']
        if methods is not None: interaction['methods'] = methods
        if monte_carlo is not None: interaction['monte_carlo'] = monte_carlo
        if matched is None:
            del interaction['matched_design']
        else:
            interaction['matched_design'] = matched
        path=root/'job.json';path.write_text(json.dumps(job))
        env=os.environ.copy();env['R_LIBS']=str(ROOT/'.runtime/R-library')
        result=subprocess.run([os.environ.get('PROTEOMICS_RSCRIPT','Rscript'),str(ROOT/'pipeline_lib/drug_interaction.R'),str(path)],env=env,text=True,capture_output=True,timeout=120)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        return root/'results/optional_analyses/drug_interaction'

    def test_optional_genes_do_not_change_protein_interactions(self):
        outputs = [self.run_model(methods=['edgeR','limma'],gene_counts=enabled) for enabled in (True,False)]
        for method in ('', 'limma'):
            self.assertEqual((outputs[0]/'cells/dda/protein'/method/'interaction.csv').read_bytes(),
                             (outputs[1]/'cells/dda/protein'/method/'interaction.csv').read_bytes())
        self.assertTrue((outputs[0]/'cells/dda/gene').is_dir())
        self.assertFalse((outputs[1]/'cells/dda/gene').exists())

    def test_log_additive_null_positive_negative_and_hla_scope(self):
        out=self.run_model(hla=True)
        rows={r['FeatureID']:r for r in optional.read_csv(out/'cells/dda/peptide/interaction.csv')}
        self.assertLess(abs(float(rows['P0']['interaction_log2FC'])),.25)
        self.assertGreater(float(rows['P1']['interaction_log2FC']),1.7)
        self.assertLess(float(rows['P2']['interaction_log2FC']),-1.7)
        self.assertLess(float(rows['P1']['FDR']),.05)
        self.assertEqual(rows['P1']['dose_time_matching_confirmed'],'TRUE')
        self.assertEqual(rows['P1']['analysis_interpretation'],'proteomic interaction')
        self.assertFalse((out/'cells/dda/protein').exists())

    def test_unconfirmed_matching_runs_and_labels_results_without_changing_statistics(self):
        reference = optional.read_csv(self.run_model(hla=True)/'cells/dda/peptide/interaction.csv')
        for matched in (False, None):
            with self.subTest(matched=matched):
                out = self.run_model(hla=True,matched=matched)
                rows = optional.read_csv(out/'cells/dda/peptide/interaction.csv')
                for expected, row in zip(reference, rows):
                    for field in ('FeatureID','interaction_log2FC','PValue','FDR','significant','test_status'):
                        self.assertEqual(row[field],expected[field])
                self.assertEqual(len(rows),len(reference))
                for table in ('interaction.csv','significant.csv'):
                    reported = optional.read_csv(out/'cells/dda/peptide'/table)
                    self.assertTrue(reported)
                    self.assertTrue(all(r['dose_time_matching_confirmed']=='FALSE' and
                                        r['analysis_interpretation']=='exploratory combination contrast' for r in reported))
                self.assertIn('were not confirmed',(out/'README.md').read_text())
                settings = json.loads((out/'settings.json').read_text())
                self.assertIs(settings['dose_time_matching_confirmed'],False)
                self.assertEqual(settings['analysis_interpretation'],'exploratory combination contrast')
                self.assertEqual(optional.read_csv(out/'analysis_summary.csv')[0]['dose_time_matching_confirmed'],'FALSE')

    def test_zero_coverage_returns_labeled_full_tables_without_failure(self):
        out=self.run_model(zero=True,matched=False)
        for kind in ('protein','gene','peptide'):
            rows=optional.read_csv(out/f'cells/dda/{kind}/interaction.csv')
            self.assertEqual(len(rows),350)
            self.assertTrue(all(r['PValue']=='' and r['FDR']=='' and r['significant']=='FALSE' for r in rows))
            self.assertTrue(all('insufficient evidence' in r['test_status'] for r in rows))
            self.assertTrue(all(r['analysis_interpretation']=='exploratory combination contrast' for r in rows))
            self.assertEqual(optional.read_csv(out/f'cells/dda/{kind}/significant.csv'),[])


@unittest.skipUnless(shutil.which('node'),'Node unavailable')
class OptionalBrowserTests(unittest.TestCase):
    def test_drug_menu_gating_and_config_roundtrip(self):
        html=(ROOT/'assets/config_ui.html').read_text()
        collect=html[html.index('function collectOptionalStatistics'):html.index('function collect(')]
        update=html[html.index('function updateInteractionSummary'):html.index("for(const id of ['drugs-applied'")]
        harness=r'''
const assert=require('node:assert/strict'),vm=require('node:vm');
const elements=new Map();function $(id){if(!elements.has(id))elements.set(id,{value:'',checked:false,disabled:false,hidden:false});return elements.get(id);}
const context={$,cfg:{samples:SAMPLES}};vm.createContext(context);vm.runInContext(SOURCE,context);
$('interaction-enabled').checked=true;$('interaction-edgeR').checked=true;$('drug-study-type').value='synergy';$('drugs-applied').value='no';
context.updateDrugOptions();assert.equal($('drug-options').hidden,true);assert.equal($('interaction-enabled').disabled,true);
assert.equal('interaction' in context.collectOptionalOptions().drugs,false);
$('drugs-applied').value='yes';$('drug-study-type').value='';context.updateDrugOptions();assert.equal($('drug-options').hidden,false);assert.throws(()=>context.collectOptionalOptions(),/synergy or single/);
$('drug-study-type').value='single';context.updateDrugOptions();assert.equal($('synergy-options').hidden,true);assert.equal('interaction' in context.collectOptionalOptions().drugs,false);
$('drug-study-type').value='synergy';context.updateDrugOptions();assert.equal($('interaction-options').hidden,false);
for(const [id,value] of Object.entries({'cohort':'cells','control':'Control','drug-a':'A','drug-b':'B','combination':'AB','fdr':'.05'}))$('interaction-'+id).value=value;
$('interaction-matched').checked=false;context.updateDrugOptions();assert.equal($('interaction-design-note').hidden,false);
assert.equal(context.collectOptionalOptions().drugs.interaction.matched_design,false);
$('interaction-matched').checked=true;context.updateDrugOptions();assert.equal($('interaction-design-note').hidden,true);
$('immuno-enabled').value='yes';
let saved=context.collectOptionalOptions();assert.equal(saved.hla_enriched,true);assert.equal('open_search' in saved,false);assert.equal(saved.drugs.interaction.combination,'AB');
$('interaction-monte_carlo').checked=true;$('interaction-limma').checked=true;context.updateDrugOptions();
assert.equal($('interaction-mc-options').hidden,false);assert.equal($('interaction-mc-simulations').disabled,false);
assert.throws(()=>context.collectOptionalOptions(),/Monte Carlo simulations/);
for(const [key,value] of Object.entries({simulations:'999',seed:'42',workers:'2'}))$('interaction-mc-'+key).value=value;
saved=context.collectOptionalOptions().drugs.interaction;
assert.equal(JSON.stringify(saved.methods),'["edgeR","limma","monte_carlo"]');assert.equal(saved.monte_carlo.seed,42);
$('interaction-mc-simulations').value='0';assert.throws(()=>context.collectOptionalOptions(),/Monte Carlo simulations/);
$('interaction-monte_carlo').checked=false;context.updateDrugOptions();
assert.equal($('interaction-mc-seed').disabled,true);assert.equal('monte_carlo' in context.collectOptionalOptions().drugs.interaction,false);
$('interaction-edgeR').checked=false;$('interaction-limma').checked=false;
assert.throws(()=>context.collectOptionalOptions(),/at least one/);$('interaction-limma').checked=true;
$('interaction-drug-b').value='A';assert.throws(()=>context.collectOptionalOptions(),error=>error.optionalTab==='drugs'&&error.field==='interaction-drug-b'&&/different condition/.test(error.message));
$('drug-study-type').value='single';context.updateDrugOptions();saved=context.collectOptionalOptions();assert.equal('interaction' in saved.drugs,false);
$('drugs-applied').value='no';context.updateDrugOptions();assert.equal(JSON.stringify(context.collectOptionalOptions().drugs),'{{"applied":false}}'.replace('{{','{').replace('}}','}'));
'''.replace('SOURCE',json.dumps(collect+'\n'+update)).replace('SAMPLES',json.dumps(samples()))
        result=subprocess.run(['node','-'],input=harness,text=True,capture_output=True,timeout=15)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
