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

import jsonschema

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'pipeline_lib'))
sys.path.insert(0, str(ROOT/'tests'))
import immunopeptidomics as hla
import pipeline
import reporting
from final_results import assemble
from make_fixture import make_fixture


class HlaTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.plan = json.loads(make_fixture(self.root).read_text())
        self.plan.update(modes=['DDA'], proteomics={'confirmed': True, 'acquisition': 'DDA-PASEF',
                         'analysis': 'immunopeptidomics'}, immunopeptidomics={'hla_class': 'I'})
        self.plan['samples'] = [self.plan['samples'][0]]
        self.sid = self.plan['samples'][0]['sample_id']
        self.run = self.root/'run'
        (self.run/self.sid).mkdir(parents=True)
        (self.run/'SEARCH_COMPLETE').touch()
        shutil.copy2(ROOT/'assets/fragpipe_hla.workflow', self.run/'fragpipe.workflow')
        self.out = self.root/'tables_cell_dda'/'immunopeptidomics'

    def row(self, spectrum='scan1', peptide='PEPTIDEAA', **values):
        return {'Spectrum': spectrum, 'Peptide': peptide, 'Protein': 'sp|TEST1|TEST_HUMAN',
                'Gene': 'G1', 'Mapped Genes': '', 'Mapped Proteins': '', 'Assigned Modifications': '',
                'Is Decoy': 'false', 'Is Contaminant': 'false', 'Charge': '2', 'Probability': '0.998', **values}

    def psms(self, rows):
        with (self.run/self.sid/'psm.tsv').open('w', newline='') as f:
            writer = csv.DictWriter(f, list(rows[0]), delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)

    def rows(self, path):
        with Path(path).open() as handle:
            return list(csv.DictReader(handle))

    def test_search_preserves_hla_fdr_and_does_not_change_global_search(self):
        args = SimpleNamespace(template=ROOT/'assets/fragpipe_hla.workflow', database=self.plan['database']['fasta'],
                               cohort='cell', mode='DDA', cpus=1, memory=4)
        props = pipeline.search_properties(args, self.plan)
        hla.validate_workflow(props)
        self.assertEqual(props['phi-report.filter'], '--psm 0.01 --pep 1 --ion 1 --prot=1')
        self.assertEqual(props['msfragger.fragment_mass_units'], '0')
        self.assertEqual(props['msfragger.fragment_mass_tolerance'], '0.03')
        self.assertEqual(props['msfragger.precursor_mass_lower'], '-15')
        self.assertEqual(props['msfragger.precursor_mass_upper'], '15')
        self.assertNotIn('STY,true', props['msfragger.table.var-mods'])
        self.assertIn('119.0041,C,true', props['msfragger.table.var-mods'])
        self.assertEqual(props['ionquant.mbr'], '0')
        args.mode = 'DDA+'
        with self.assertRaisesRegex(ValueError, 'DDA searches only'):
            pipeline.search_properties(args, self.plan)
        args.template = ROOT/'assets/fragpipe_lfq_mbr.workflow'
        global_plan = copy.deepcopy(self.plan)
        global_plan['proteomics']['analysis'] = 'global-label-free'
        props = pipeline.search_properties(args, global_plan)
        self.assertEqual(props['ionquant.mbr'], '1')
        self.assertEqual(props['phi-report.filter'], '--sequential --prot 0.01 --picked')

    def test_schema_requires_explicit_class_and_compatible_search(self):
        validator = jsonschema.Draft202012Validator(json.loads((ROOT/'assets/config.schema.json').read_text()))
        validator.validate(self.plan)
        invalid = []
        plan = copy.deepcopy(self.plan);plan.pop('immunopeptidomics');invalid.append(plan)
        plan = copy.deepcopy(self.plan);plan['modes'] = ['DDA+'];invalid.append(plan)
        plan = copy.deepcopy(self.plan);plan['immunopeptidomics']['hla_class'] = 'unknown';invalid.append(plan)
        plan = copy.deepcopy(self.plan);plan['differential'] = {'enabled': True, 'modifications': True};invalid.append(plan)
        plan = copy.deepcopy(self.plan);plan['proteomics']['analysis'] = 'phosphoproteomics';invalid.append(plan)
        for plan in invalid:
            self.assertFalse(validator.is_valid(plan), plan)

    def test_counts_independent_of_protein_groups_and_duplicate_mappings(self):
        self.psms([self.row(), self.row(**{'Mapped Proteins': 'P2', 'Mapped Genes': 'G2'}),
                   self.row('scan2'), self.row('scan3', 'ANOTHERPEP'),
                   self.row('decoy', **{'Is Decoy': 'true'}),
                   self.row('contaminant', Protein='contam_P1')])
        # No combined_protein.tsv exists: the peptide exporter must not require it.
        hla.export(self.plan, self.run, 'cell', self.out)
        rows = {r['Peptide']: r for r in self.rows(self.out/'peptides.csv')}
        self.assertEqual(rows['PEPTIDEAA'][self.sid], '2')
        self.assertEqual(rows['PEPTIDEAA']['Genes'], 'G1; G2')
        self.assertIn('P2', rows['PEPTIDEAA']['Proteins'])
        summary = self.rows(self.out/'sample_summary.csv')[0]
        self.assertEqual(summary['retained_psms'], '3')
        self.assertEqual(summary['duplicate_mapping_rows'], '1')
        self.assertEqual(summary['decoy_or_contaminant_rows'], '2')
        self.assertEqual(len(self.rows(self.out/'length_distribution.csv')), 19)
        self.assertEqual(set((self.out/'binding_input.txt').read_text().splitlines()), set(rows))
        matrix = self.rows(self.out.parent/'peptide_spectral_counts.csv')
        self.assertEqual(next(r for r in matrix if r['Peptide']=='PEPTIDEAA')['Mapped Genes'], 'G2')

    def test_species_filter_uses_fasta_not_inferred_proteins(self):
        self.plan['organism'] = 'Homo sapiens'
        self.psms([self.row(), self.row('mouse', Protein='sp|MOUSE|OTHER')])
        hla.export(self.plan, self.run, 'cell', self.out)
        self.assertEqual(self.rows(self.out/'peptides.csv')[0][self.sid], '1')
        self.assertEqual(self.rows(self.out/'sample_summary.csv')[0]['other_organism_rows'], '1')

    def test_class_two_export_and_empty_sample(self):
        self.plan['immunopeptidomics']['hla_class'] = 'II'
        self.psms([self.row(peptide='PEPTIDE'), self.row('long', 'PEPTIDEPEPTIDEPEP')])
        hla.export(self.plan, self.run, 'cell', self.out)
        self.assertEqual((self.out/'binding_input.txt').read_text(), 'PEPTIDEPEPTIDEPEP\n')
        self.psms([self.row(**{'Is Decoy': 'true'})])
        hla.export(self.plan, self.run, 'cell', self.out)
        self.assertEqual(self.rows(self.out/'peptides.csv'), [])
        self.assertEqual(self.rows(self.out/'sample_summary.csv')[0]['retained_psms'], '0')

    def test_partial_wrong_search_and_conflicting_psms_are_rejected(self):
        self.psms([self.row(), self.row(peptide='DIFFERENT')])
        with self.assertRaisesRegex(ValueError, 'multiple assignments'):
            hla.export(self.plan, self.run, 'cell', self.out)
        shutil.copy2(ROOT/'assets/fragpipe_lfq_mbr.workflow', self.run/'fragpipe.workflow')
        with self.assertRaisesRegex(ValueError, 'workflow setting'):
            hla.export(self.plan, self.run, 'cell', self.out)
        (self.run/'SEARCH_COMPLETE').unlink()
        with self.assertRaisesRegex(ValueError, 'completed scientific search'):
            hla.export(self.plan, self.run, 'cell', self.out)

    def test_reporting_and_final_results_include_hla_evidence(self):
        self.psms([self.row(), self.row('scan2', 'ANOTHERPEP')])
        self.plan['optional_analyses'] = {'hla_enriched': True, 'peptide_isoforms': True, 'drugs': {'applied': False}}
        protein = {'Protein': 'P_NOT_IN_PSMS', 'Organism': 'Homo sapiens', 'Gene': 'G_OTHER',
                   'Indistinguishable Proteins': '', self.sid+' Spectral Count': '1', self.sid+' Intensity': '100'}
        with (self.run/'combined_protein.tsv').open('w') as f:
            writer = csv.DictWriter(f, list(protein), delimiter='\t');writer.writeheader();writer.writerow(protein)
        (self.run/'run_provenance.json').write_text(json.dumps({'stub': False, 'mode': 'DDA'}))
        self.plan['differential'] = {'enabled': False}
        planfile = self.root/'plan.json';planfile.write_text(json.dumps(self.plan))
        reporting.report(SimpleNamespace(plan=planfile, run=self.run, cohort='cell', mode='DDA', out=self.out.parent))
        out = self.root/'final_results'
        assemble(self.plan, [self.out.parent], out, os.environ.get('PROTEOMICS_RSCRIPT','Rscript'), ROOT/'.runtime/R-library')
        self.assertEqual(len(self.rows(out/'main_results/counts/cell/dda/raw/peptide_counts.csv')), 2)
        self.assertTrue((out/'optional_analyses/immunopeptidomics/cell/dda/peptides.csv').is_file())
        self.assertIn('immunopeptidomics/cell/dda/peptides.csv', (out/'START_HERE.html').read_text())
        for suffix in ('peptide_isoforms/cell/dda/isoform_evidence.csv',):
            self.assertTrue((out/'optional_analyses'/suffix).is_file())
            self.assertIn('optional_analyses/'+suffix,(out/'START_HERE.html').read_text())
        self.assertEqual(json.loads((out.parent/'analysis_details/reporting/result_audit/main_statistics/analysis_settings.json').read_text())['analysis'], 'immunopeptidomics')
        self.assertIn('1% PSM FDR', json.loads((self.out.parent/'definitions.json').read_text())['fdr'])
        self.assertTrue((out/'FINAL_RESULTS_COMPLETE').is_file())


@unittest.skipUnless(shutil.which('node'), 'Node unavailable for browser logic tests')
class HlaBrowserTests(unittest.TestCase):
    def test_tab_order_toggle_roundtrip_and_mode_restoration(self):
        html = (ROOT/'assets/config_ui.html').read_text()
        collect = html[html.index('const countThresholds='):html.index('function renderSamples()')]
        update = html[html.index('function updateDEPOptions()'):html.index("\n$('de-enabled').onchange")]
        immune = html[html.index('let previousSearchModes='):html.index('setupWizard();', html.index('let previousSearchModes='))]
        statistics = html[html.index('function updateOptionalStatistics()'):html.index("for(const id of ['optional-stats-enabled'")]
        selected = html[html.index('function selectedAnalysis()'):html.index('\n', html.index('function selectedAnalysis()'))]
        steps = html[html.index('const stepIds='):html.index('let currentStep=')]
        harness = r'''
const assert=require('node:assert/strict'),vm=require('node:vm');
const elements=new Map();function $(id){if(!elements.has(id))elements.set(id,{value:'',checked:false,disabled:false,hidden:false});return elements.get(id);}
const context={$,Number,structuredClone,cfg:{proteomics:{confirmed:true,acquisition:'DDA-PASEF',analysis:'phosphoproteomics'},database:{}},
 conditionLabels:{},phosphoModule:{installed:true,sha256:'hash'},fastaDownloadReady:false,effectiveOutputFolder:()=>'/output',invalidateStep:()=>{}};
vm.createContext(context);vm.runInContext(SOURCE,context);
assert.equal(vm.runInContext("stepIds.indexOf('optional-analyses-section')",context),5);
assert.equal(vm.runInContext('RUN_STEP',context),7);
$('modes').value='both';$('fasta-source').value='local';$('include-phospho').checked=true;
$('de-enabled').checked=true;$('de-edgeR').checked=true;$('de-control').value='Control';$('de-modifications').checked=true;
$('immuno-enabled').value='yes';$('immuno-class').value='II';context.updateImmunoOptions();
const enabled=context.collect();assert.equal(enabled.proteomics.analysis,'immunopeptidomics');assert.equal(enabled.immunopeptidomics.hla_class,'II');
assert.equal(JSON.stringify(enabled.modes),'["DDA"]');assert.equal(enabled.differential.modifications,false);
assert.equal('module_sha256' in enabled.proteomics,false);assert.equal($('modes').disabled,true);assert.equal($('de-modifications').disabled,true);
assert.match($('search-mode-summary').textContent,/HLA-enriched preset/);
assert.match($('review-search-modes').textContent,/no DDA\+ results/);
$('immuno-enabled').value='no';context.updateImmunoOptions();const off=context.collect();
assert.equal(off.proteomics.analysis,'phosphoproteomics');assert.equal(off.proteomics.module_sha256,'hash');
assert.equal('immunopeptidomics' in off,false);assert.equal($('modes').value,'both');assert.equal($('modes').disabled,false);
assert.equal($('de-modifications').disabled,false);assert.equal(off.differential.modifications,true);
assert.match($('search-mode-summary').textContent,/Two separate searches/);
// A reopened DDA-only study must explicitly say that DDA+ will not be produced.
$('modes').value='DDA';$('modes').onchange();
assert.match($('search-mode-summary').textContent,/no DDA\+ results/);
$('immuno-enabled').value='yes';context.updateImmunoOptions();
$('immuno-enabled').value='no';context.updateImmunoOptions();
assert.equal($('modes').value,'DDA');assert.match($('review-search-modes').textContent,/DDA only/);
$('modes').value='DDA+';$('modes').onchange();assert.match($('review-search-modes').textContent,/no DDA results/);
'''.replace('SOURCE', json.dumps('\n'.join([collect, update, selected, steps, statistics, immune])))
        result = subprocess.run(['node','-'], input=harness, text=True, capture_output=True, timeout=15)
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
