"""Validate optional statistical inputs, browser settings and real model outputs."""
import copy
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import jsonschema

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'pipeline_lib'))
import optional_statistics as stats


def plan():
    return {'modes': ['DDA'], 'differential': {'enabled': False},
            'samples': [{'sample_id': f'S{i}', 'condition': 'Control' if i < 4 else 'Treatment',
                         'cohort': 'cells', 'replicate': i % 4+1} for i in range(8)],
            'optional_analyses': {'peptide_isoforms': True,
                'statistics': {'enabled': True, 'features': ['peptide', 'isoform'],
                               'methods': ['edgeR', 'limma'], 'comparisons': 'control',
                               'control': 'Control', 'fdr': 0.1, 'filter_low_counts': False}}}


class OptionalStatisticsTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory(); self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.plan = plan()
        self.ids = [s['sample_id'] for s in self.plan['samples']]

    def test_validation_rejects_unsupported_evidence_and_designs(self):
        schema = json.loads((ROOT/'assets/config.schema.json').read_text())['properties']['optional_analyses']
        validator = jsonschema.Draft202012Validator(schema)
        validator.validate(self.plan['optional_analyses'])
        stats.validate_statistics(self.plan)
        for edit in [lambda p: p['optional_analyses'].update(peptide_isoforms=False),
                     lambda p: p['optional_analyses']['statistics'].update(control='Missing'),
                     lambda p: p['optional_analyses']['statistics'].update(methods=[]),
                     lambda p: p['optional_analyses']['statistics'].update(fdr=0),
                     lambda p: p['optional_analyses']['statistics'].update(filter_low_counts=True, min_total_count=0),
                     lambda p: p['samples'][1].update(replicate=1)]:
            bad = copy.deepcopy(self.plan); edit(bad)
            with self.assertRaises(ValueError): stats.validate_statistics(bad)
        self.plan['optional_analyses']['statistics'].update(filter_low_counts=False, min_count=0, min_total_count=0)
        validator.validate(self.plan['optional_analyses']);stats.validate_statistics(self.plan)
        self.plan['proteomics'] = {'analysis': 'immunopeptidomics'}
        self.plan['optional_analyses'].update(peptide_isoforms=False)
        self.plan['optional_analyses']['statistics']['features'] = ['peptide']
        stats.validate_statistics(self.plan)
        self.plan['optional_analyses']['statistics']['features'] = ['protein']
        self.assertFalse(validator.is_valid(self.plan['optional_analyses']))

    def test_isoform_comparisons_exclude_canonical_only_and_shared_evidence(self):
        fields = ['Accession', 'Family_entries_in_FASTA', 'Isoform_discriminating_peptides', 'S0', 'S1']
        rows = [dict(zip(fields, values)) for values in [('P1', '1', '0', '10', '40'),
                ('P2', '2', '0', '0', '0'), ('P2-2', '2', '2', '3', '9')]]
        stats.write_table(self.root/'isoform_evidence.csv', fields, rows)
        _, accepted, _ = stats.isoform_counts(self.root, ['S0', 'S1'], self.root/'excluded.csv')
        self.assertEqual([r['FeatureID'] for r in accepted], ['P2-2'])
        self.assertEqual(accepted[0]['S1'], 9)
        excluded = stats.read_table(self.root/'excluded.csv')[1]
        self.assertEqual(len(excluded), 2)
        self.assertIn('alternative', excluded[0]['Exclusion_reason'])

    def fixture(self, empty_sample=False, no_isoforms=False):
        out = self.root/'final_results'
        raw = self.root/'analysis_details/reporting/cells/dda'
        raw.mkdir(parents=True)
        peptides, isoforms = [], []
        for j in range(160):
            values = [1000+i*3 if j < 8 and i >= 4 else 40+(i*13+j*7) % 25 for i in range(8)]
            if j == 158: values = [0]*4+[1]*4
            if j == 159: values = [0]*8
            if empty_sample: values[0] = 0
            counts = dict(zip(self.ids, values))
            peptides.append({'FeatureID': f'pep{j}', 'Peptide': f'pep{j}', 'Gene': f'G{j}', 'Mapped Genes': '', **counts})
            isoforms.append({'Accession': f'P{j}-2', 'Family_entries_in_FASTA': 1 if no_isoforms else 2,
                            'Isoform_discriminating_peptides': 0 if no_isoforms else 1, **counts})
        stats.write_table(raw/'peptide_counts.csv', list(peptides[0]), peptides)
        stats.write_table(out/'optional_analyses/peptide_isoforms/cells/dda/isoform_evidence.csv', list(isoforms[0]), isoforms)
        return out, raw

    def run_stats(self, **kwargs):
        out, raw = self.fixture(**kwargs)
        env = {**os.environ, 'R_LIBS': str(ROOT/'.runtime/R-library')}
        stats.run_statistics(self.plan, out, os.environ.get('PROTEOMICS_RSCRIPT', 'Rscript'), env)
        return out/'optional_analyses/statistics', raw

    def test_models_effect_direction_fdr_offsets_and_retained_filtered_rows(self):
        self.plan['optional_analyses']['statistics'].update(filter_low_counts=True, min_count=5, min_total_count=20)
        out, _ = self.run_stats()
        index = stats.read_table(out/'differential/results_index.csv')[1]
        self.assertEqual({row['feature_type'] for row in index}, {'peptide', 'isoform'})
        for kind, category in stats.CATEGORIES.items():
            for method in ['edgeR', 'limma']:
                path = out/f'differential/{method}/pairwise/cells/dda/{category}/{method}_Control_vs_Treatment.csv'
                rows = stats.read_table(path)[1]
                self.assertEqual(len(rows), 160)
                up = next(row for row in rows if row['FeatureID'] in ('pep0', 'P0-2'))
                self.assertGreater(float(up['logFC']), 2)
                self.assertEqual(up['significant'].lower(), 'true')
                self.assertEqual(up['count_unit'], stats.COUNT_UNITS[kind])
                self.assertEqual(up['normalization'], 'TMMwsp')
                for row in rows:
                    if row['FeatureID'] in ('pep158', 'P158-2', 'pep159', 'P159-2'):
                        self.assertEqual(row['PValue'], '');self.assertEqual(row['FDR'], '')
                        self.assertEqual(row['significant'].lower(), 'false')
                        self.assertIn('not tested:', row['test_status'])
                tested = sorted((float(row['PValue']), float(row['FDR'])) for row in rows if row['PValue'])
                previous = 1.0
                for rank in range(len(tested), 0, -1):
                    p, adjusted = tested[rank-1];previous = min(previous, p*len(tested)/rank)
                    self.assertAlmostEqual(adjusted, previous, places=9)
            factors = stats.read_table(out/f'inputs/cells/dda/models/{category}/counts_model_normalization_factors.csv')[1]
            counts = stats.read_table(out/f'inputs/cells/dda/{kind}_counts.csv')[1]
            self.assertEqual([int(row['library_size']) for row in factors], [sum(int(row[sid]) for row in counts) for sid in self.ids])
        job = json.loads((out/'analysis_settings.json').read_text())
        self.assertTrue((out/job['groups'][0]['directory']/'samples.csv').is_file())
        self.assertTrue((out/'zscore/cells/dda/isoform_zscore.csv').is_file())
        self.assertNotIn('protein_counts.csv', (out/'START_HERE.html').read_text())

    def test_empty_library_and_no_isoform_support_finish_with_explicit_status(self):
        out, _ = self.run_stats(empty_sample=True, no_isoforms=True)
        for category in ('peptides', 'isoforms'):
            rows = stats.read_table(out/f'differential/edgeR/pairwise/cells/dda/{category}/edgeR_Control_vs_Treatment.csv')[1]
            self.assertEqual(len(rows), 0 if category == 'isoforms' else 160)
            self.assertTrue(all(row['PValue'] == '' and row['FDR'] == '' and 'insufficient evidence' in row['test_status'] for row in rows))
        self.assertEqual(len(stats.read_table(out/'inputs/cells/dda/excluded_isoforms.csv')[1]), 160)
        summary = stats.read_table(out/'analysis_summary.csv')[1]
        self.assertTrue(all(row['tested_features'] == '0' for row in summary))

    def test_final_assembly_runs_optional_statistics_with_deps_disabled(self):
        from final_results import assemble
        out, raw = self.fixture()
        self.plan['optional_analyses']['statistics']['features'] = ['peptide']
        report = self.root/'tables_cells_dda';report.mkdir()
        with (report/'protein_spectral_counts.tsv').open('w') as handle:
            handle.write('\t'.join(['Protein']+self.ids)+'\n'+'\t'.join(['P1']+['50']*8)+'\n')
        fields, rows = stats.read_table(raw/'peptide_counts.csv')
        stats.write_table(report/'peptide_spectral_counts.csv', [f for f in fields if f != 'FeatureID'],
                          [{k: v for k, v in row.items() if k != 'FeatureID'} for row in rows])
        for filename, fields in [('modified_peptide_psm_counts.csv', ['Peptide', 'Assigned Modifications']),
                                 ('gene_modification_peptide_counts.csv', ['Gene', 'Modification']),
                                 ('protein_modification_peptide_counts.csv', ['Protein', 'Modification'])]:
            stats.write_table(report/'ptm_counts'/filename, fields+self.ids, [])
        from pipeline import digest
        fasta=self.root/'targets.fasta';fasta.write_text('>P1\nPEPTIDE\n')
        self.plan['database']={'fasta':str(fasta)}
        evidence=report/'optional_analyses/peptide_isoforms'
        shutil.copytree(out/'optional_analyses/peptide_isoforms/cells/dda',evidence)
        (evidence/'settings.json').write_text(json.dumps({'database_sha256':digest(fasta),
            'peptide_counts_sha256':digest(report/'peptide_spectral_counts.csv')}))
        # Benchmark staged LFQ reports before Nextflow publishes any raw search folder.
        stats.write_table(report/'protein_maxlfq_intensity.tsv', ['Protein',*self.ids],
                          [dict(Protein='P1',**{key:100 for key in self.ids})])
        text=(report/'protein_maxlfq_intensity.tsv').read_text()
        (report/'protein_maxlfq_intensity.tsv').write_text(text.replace(',', '\t'))
        service=self.root/'service.csv'
        aliases=['provider_'+key for key in self.ids]
        stats.write_table(service,['Protein.Group',*aliases],
                          [dict({'Protein.Group':group},**{key:100 for key in aliases}) for group in ['P1','P2']])
        self.plan['service_quantification']={'cells':{'matrix':str(service),'sha256':digest(service),
            'samples':[{'sample_id':key,'service_sample':alias} for key,alias in zip(self.ids,aliases)]}}
        assemble(self.plan, [report], out, os.environ.get('PROTEOMICS_RSCRIPT', 'Rscript'), ROOT/'.runtime/R-library')
        benchmark=out/'benchmark/cells/dda'
        summary=stats.read_table(benchmark/'detection_summary.csv')[1][-1]
        self.assertEqual(summary['shared'],'1')
        self.assertEqual(summary['service_only'],'1')
        self.assertEqual(summary['net_difference'],'-1')
        self.assertTrue((benchmark/'reported_detection_summary.tsv').is_file())

        self.assertTrue((out/'FINAL_RESULTS_COMPLETE').is_file())
        self.assertFalse((out/'differential').exists())
        index=stats.read_table(out/'results_index.csv')[1]
        self.assertTrue(index)
        self.assertTrue(all(row['section']=='Optional analyses' for row in index))
        self.assertTrue(all((out/row['full_results']).is_file() for row in index))
        self.assertTrue((out/'main_results/counts/cells/dda/raw/peptide_counts.csv').is_file())

    def test_single_peptide_reports_status_without_failing(self):
        out, _ = self.fixture()
        self.plan['optional_analyses']['statistics']['features'] = ['peptide']
        path = self.root/'analysis_details/reporting/cells/dda/peptide_counts.csv'
        path.write_text('\n'.join(path.read_text().splitlines()[:2])+'\n')
        stats.run_statistics(self.plan, out, os.environ.get('PROTEOMICS_RSCRIPT', 'Rscript'),
                             {**os.environ, 'R_LIBS': str(ROOT/'.runtime/R-library')})
        for method in ['edgeR', 'limma']:
            rows = stats.read_table(out/f'optional_analyses/statistics/differential/{method}/pairwise/cells/dda/peptides/{method}_Control_vs_Treatment.csv')[1]
            self.assertEqual(len(rows), 1)
            self.assertIn('test_status', rows[0])


@unittest.skipUnless(shutil.which('node'), 'Node unavailable')
class OptionalStatisticsBrowserTests(unittest.TestCase):
    def test_saved_statistics_reload_and_missing_reference_remains_explicit(self):
        html = (ROOT/'assets/config_ui.html').read_text()
        collect = html[html.index('const countThresholds='):html.index('function renderSamples()')]
        controls = html[html.index('function refreshDrugConditions()'):html.index('let previousSearchModes=')]
        saved = plan()
        saved['optional_analyses']['statistics'].update(filter_low_counts=True, min_count=2, min_total_count=7)
        harness = r'''
const assert=require('node:assert/strict'),vm=require('node:vm');
const elements=new Map();
function $(id){if(!elements.has(id))elements.set(id,{value:'',checked:false,disabled:false,hidden:false,
 replaceChildren(){},add(){},setAttribute(){},classList:{toggle(){}},focus(){}});return elements.get(id);}
const context={$,Number,cfg:CONFIG,Option:function(text,value){this.text=text;this.value=value;},invalidateStep(){}};
vm.createContext(context);vm.runInContext(SOURCE,context);context.loadOptionalOptions();
assert.equal(JSON.stringify(context.collectOptionalOptions().statistics),JSON.stringify(context.cfg.optional_analyses.statistics));
assert.equal($('optional-stats-min-total').value,7);assert.equal($('optional-stats-min-total').disabled,false);
assert.equal($('optional-stats-isoforms').disabled,false);
context.showOptionalTab('statistics');assert.equal($('optional-panel-statistics').hidden,false);
context.cfg.samples=context.cfg.samples.map(sample=>({...sample,condition:sample.condition==='Control'?'New reference':sample.condition}));
context.refreshOptionalStatisticsConditions();assert.equal($('optional-stats-control').value,'');
assert.throws(()=>context.collectOptionalOptions(),/reference condition/);
'''.replace('SOURCE', json.dumps(collect+'\n'+controls)).replace('CONFIG', json.dumps(saved))
        result = subprocess.run(['node', '-'], input=harness, text=True, capture_output=True, timeout=15)
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)

    def test_statistical_settings_gating_and_config_roundtrip(self):
        html = (ROOT/'assets/config_ui.html').read_text()
        collect = html[html.index('const countThresholds='):html.index('function renderSamples()')]
        update = html[html.index('function updateOptionalStatistics()'):html.index("for(const id of ['optional-stats-enabled'")]
        harness = r'''
const assert=require('node:assert/strict'),vm=require('node:vm');
const elements=new Map();function $(id){if(!elements.has(id))elements.set(id,{value:'',checked:false,disabled:false,hidden:false});return elements.get(id);}
const context={$,Number};vm.createContext(context);vm.runInContext(SOURCE,context);
assert.equal(context.collectOptionalOptions().statistics.enabled,false);
$('optional-stats-enabled').checked=true;$('optional-stats-edgeR').checked=true;
$('optional-stats-control').value='__all_pairs__';$('optional-stats-fdr').value='0.01';
for(const field of ['peptides','isoforms'])$('optional-stats-'+field).checked=true;
context.updateOptionalStatistics();assert.equal($('optional-stats-peptides').disabled,true);
assert.throws(()=>context.collectOptionalOptions(),/parent analysis/);
$('immuno-enabled').value='yes';context.updateOptionalStatistics();
let saved=context.collectOptionalOptions().statistics;
assert.equal(JSON.stringify(saved.features),'["peptide"]');assert.equal(saved.fdr,0.01);assert.equal(saved.comparisons,'all_pairwise');
$('optional-stats-min-total').value='0';assert.equal('min_total_count' in context.collectOptionalOptions().statistics,false);
$('optional-stats-filter').checked=true;context.updateOptionalStatistics();
assert.throws(()=>context.collectOptionalOptions(),/minimum counts/);
$('optional-stats-min-count').value='2';$('optional-stats-min-total').value='7';
$('optional-stats-control').value='Control';$('peptide-isoforms-enabled').checked=true;
context.updateOptionalStatistics();saved=context.collectOptionalOptions().statistics;
assert.equal(saved.min_total_count,7);assert.equal(saved.control,'Control');assert.equal(saved.features.length,2);
$('peptide-isoforms-enabled').checked=false;$('immuno-enabled').value='yes';
context.updateOptionalStatistics();saved=context.collectOptionalOptions().statistics;
assert.equal(JSON.stringify(saved.features),'["peptide"]');assert.equal($('optional-stats-isoforms').disabled,true);
$('optional-stats-enabled').checked=false;context.updateOptionalStatistics();
assert.equal(JSON.stringify(context.collectOptionalOptions().statistics),'{"enabled":false}');
'''.replace('SOURCE', json.dumps(collect+'\n'+update))
        result = subprocess.run(['node', '-'], input=harness, text=True, capture_output=True, timeout=15)
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
