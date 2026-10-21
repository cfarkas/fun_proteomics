"""Regression checks for moderated and simulated four-condition interactions."""
import copy
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

import jsonschema
import test_optional_analyses as fixtures

ROOT = Path(__file__).resolve().parents[1]
METHODS = ['edgeR', 'limma', 'monte_carlo']
read_csv = fixtures.optional.read_csv


class InteractionSettingsTests(unittest.TestCase):
    def test_method_and_simulation_validation(self):
        plan = {'samples': fixtures.samples(), 'optional_analyses': {'drugs': fixtures.design()}}
        settings = plan['optional_analyses']['drugs']['interaction']
        settings.update(methods=METHODS, monte_carlo={'simulations': 199, 'seed': 0, 'workers': 2})
        schema = json.loads((ROOT/'assets/config.schema.json').read_text())['properties']['optional_analyses']
        validator = jsonschema.Draft202012Validator(schema)
        validator.validate(plan['optional_analyses'])
        fixtures.optional.validate_options(plan)
        invalid = [{'methods': []}, {'methods': ['limma', 'limma']}, {'methods': ['bayesian']},
                   {'methods': 'limma'}, {'fdr': 0}, {'fdr': True}]
        for key, values in [('simulations', [98, 1000001, 99.5, True]),
                            ('seed', [-1, 2147483647, False]), ('workers', [0, 9, 1.5])]:
            invalid.extend({'monte_carlo': {key: value}} for value in values)
        for bad in invalid:
            with self.subTest(bad=bad):
                cfg = copy.deepcopy(plan)
                cfg['optional_analyses']['drugs']['interaction'].update(bad)
                self.assertFalse(validator.is_valid(cfg['optional_analyses']))
                with self.assertRaises(ValueError): fixtures.optional.validate_options(cfg)
        # Old configurations retain edgeR without requiring new keys.
        settings.pop('methods'); settings.pop('monte_carlo')
        fixtures.optional.validate_options(plan)


class InteractionMethodsTests(unittest.TestCase):
    run_model = fixtures.InteractionModelTests.run_model

    def test_effects_null_simulation_pvalues_and_separate_fdr(self):
        out = self.run_model(hla=True, methods=METHODS, extra_zero=True,
                             monte_carlo={'simulations': 199, 'seed': 42, 'workers': 2})
        common = out/'cells/dda/peptide'
        for method in METHODS:
            directory = common if method == 'edgeR' else common/method
            rows = {r['FeatureID']: r for r in read_csv(directory/'interaction.csv')}
            self.assertEqual(len(rows), 351)
            self.assertLess(abs(float(rows['P0']['interaction_log2FC'])), .25)
            self.assertGreater(float(rows['P0']['PValue']), .05)
            self.assertGreater(float(rows['P1']['interaction_log2FC']), 1.7)
            self.assertLess(float(rows['P2']['interaction_log2FC']), -1.7)
            for feature in ('P1', 'P2'):
                self.assertLess(float(rows[feature]['PValue']), .01)
                self.assertEqual(rows[feature]['method'], method)
            self.assertEqual(rows['AllZero']['PValue'], '')
            self.assertEqual(rows['AllZero']['FDR'], '')
            self.assertEqual(rows['AllZero']['significant'], 'FALSE')
            self.assertIn('all_zero', rows['AllZero']['test_status'])
            tested = sorted((r for r in rows.values() if r['PValue']), key=lambda r: float(r['PValue']))
            expected_fdr = 1.0
            for rank in range(len(tested), 0, -1):
                row = tested[rank-1]
                expected_fdr = min(expected_fdr, float(row['PValue'])*len(tested)/rank)
                self.assertAlmostEqual(float(row['FDR']), expected_fdr, places=12)
            if method == 'limma':
                self.assertGreater(float(rows['P1']['interaction_CI95_lower']), 0)
                self.assertLess(float(rows['P2']['interaction_CI95_upper']), 0)
                self.assertEqual(len(read_csv(directory/'voom_weights.csv')), 350)
            if method == 'monte_carlo':
                for row in tested:
                    self.assertEqual(int(row['MC_valid_simulations']), 199)
                    self.assertEqual(float(row['PValue']), (int(row['MC_exceedances'])+1)/200)
                    self.assertEqual(float(row['MC_minimum_pvalue']), .005)
                    self.assertGreaterEqual(float(row['MC_tail_probability_CI95_lower']), 0)
                    self.assertLessEqual(float(row['MC_tail_probability_CI95_upper']), 1)
                    self.assertTrue(math.isfinite(float(row['MC_standard_error'])))
                # Deliberately low simulation counts must not bypass BH resolution limits.
                self.assertEqual(rows['P1']['significant'], 'FALSE')
        normalization = {r['sample_id']: float(r['effective_library_size']) for r in read_csv(common/'normalization_factors.csv')}
        null = {r['FeatureID']: r for r in read_csv(common/'monte_carlo/null_fitted_counts.csv')}
        for feature in ('P0', 'P1', 'P2'):
            means = {condition: sum(float(null[feature][s['sample_id']])/normalization[s['sample_id']]
                                   for s in fixtures.samples() if s['condition'] == condition)/4
                     for condition in ('Control', 'A', 'B', 'AB')}
            self.assertAlmostEqual(means['Control']*means['AB']/(means['A']*means['B']), 1, places=10)
        baseline = {condition: float(null['P0'][f'S{condition}1'])/normalization[f'S{condition}1']
                    for condition in ('Control', 'A', 'B')}
        self.assertGreater(baseline['A']/baseline['Control'], 1.5)
        self.assertGreater(baseline['B']/baseline['Control'], 2.5)
        summary = read_csv(out/'analysis_summary.csv')
        self.assertEqual({r['method'] for r in summary}, set(METHODS))
        self.assertTrue(all((out/r['results_file']).is_file() for r in summary))
        saved = json.loads((out/'settings.json').read_text())
        self.assertEqual(saved['methods'], METHODS)
        self.assertEqual(saved['monte_carlo']['effective_workers'], 2)
        self.assertFalse((out/'cells/dda/protein').exists())

    def test_monte_carlo_reproducible_across_workers(self):
        outputs = [self.run_model(hla=True, methods=['monte_carlo'],
                    monte_carlo={'simulations': 99, 'seed': 6, 'workers': workers}) for workers in (1, 2)]
        rows = [read_csv(out/'cells/dda/peptide/monte_carlo/interaction.csv') for out in outputs]
        self.assertEqual(rows[0], rows[1])
        self.assertFalse((outputs[0]/'cells/dda/peptide/interaction.csv').exists())

    def test_simulated_null_does_not_inflate_false_positives(self):
        # A regression check under one known NB null, not a general calibration claim.
        script = r'''
suppressPackageStartupMessages(library(edgeR))
args <- commandArgs(TRUE)
source(args[1])
set.seed(91826)
n <- 1200L
band <- rep(c("low","moderate","high"),each=400)
means <- c(runif(400,.5,3),runif(400,10,30),runif(400,60,100))
phi <- .15+.2/sqrt(means)
group <- factor(rep(c("control","drug_a","drug_b","combination"),each=4),levels=c("control","drug_a","drug_b","combination"))
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
null_design <- cbind(Intercept=1,drug_a=as.integer(group %in% c("drug_a","combination")),drug_b=as.integer(group %in% c("drug_b","combination")))
a <- rnorm(n,0,.8)*(runif(n)<.2)
b <- rnorm(n,0,.8)*(runif(n)<.2)
mu <- outer(means,exp(seq(-.35,.35,length.out=16)))*2^(outer(a,null_design[,"drug_a"])+outer(b,null_design[,"drug_b"]))
x <- matrix(rnbinom(n*16,mu=as.vector(mu),size=rep(1/phi,16)),nrow=n,dimnames=list(paste0("F",seq_len(n)),paste0("S",seq_len(16))))
d <- calcNormFactors(DGEList(x,group=group),method="TMM")
rows <- list()
for(method in c("limma","monte_carlo")) {
 tab <- if(method=="limma")interaction_limma(d,design,c(1,-1,-1,1),args[2]) else interaction_monte_carlo(d,design,null_design,c(1,-1,-1,1),args[2],199L,44L,2L)
 rows[[method]] <- data.frame(method=method,band=band,PValue=tab$PValue)
}
write.csv(do.call(rbind,rows),file.path(args[2],"null_results.csv"),row.names=FALSE)
'''
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root/'null.R'; path.write_text(script)
            result = subprocess.run([os.environ.get('PROTEOMICS_RSCRIPT', 'Rscript'), str(path),
                str(ROOT/'pipeline_lib/drug_interaction_models.R'), str(root)],
                env={**os.environ, 'R_LIBS': str(ROOT/'.runtime/R-library')}, text=True, capture_output=True, timeout=120)
            self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
            rows = read_csv(root/'null_results.csv')
        for method in ('limma', 'monte_carlo'):
            for band in ('low', 'moderate', 'high'):
                pvalues = [float(r['PValue']) for r in rows if r['method']==method and r['band']==band]
                self.assertEqual(len(pvalues), 400)
                self.assertTrue(all(math.isfinite(p) and 0 < p <= 1 for p in pvalues))
                self.assertLess(sum(p < .05 for p in pvalues)/len(pvalues), .10)

    def test_low_coverage_and_empty_inputs_report_all_methods(self):
        for opts in ({'zero': True}, {'empty_library': True}, {'feature_count': 0}, {'feature_count': 1}):
            with self.subTest(opts=opts):
                out = self.run_model(hla=True, methods=METHODS, monte_carlo={'simulations': 99}, **opts)
                summary = read_csv(out/'analysis_summary.csv')
                self.assertEqual(len(summary), 3)
                for reported in summary:
                    rows = read_csv(out/reported['results_file'])
                    self.assertEqual(len(rows), opts.get('feature_count', 350))
                    for row in rows:
                        if not row['PValue']:
                            self.assertEqual(row['FDR'], '')
                            self.assertEqual(row['significant'], 'FALSE')
                            self.assertIn('insufficient evidence', row['test_status'])
                if opts.get('feature_count') == 1:
                    self.assertIn('insufficient evidence', next(r['status'] for r in summary if r['method']=='limma'))

    def test_regeneration_removes_deselected_method_outputs(self):
        out = self.run_model(hla=True, methods=['edgeR', 'limma'])
        job_path = out.parents[2]/'job.json'
        job = json.loads(job_path.read_text())
        job['optional_analyses']['drugs']['interaction']['methods'] = ['limma']
        job_path.write_text(json.dumps(job))
        env = os.environ.copy(); env['R_LIBS'] = str(ROOT/'.runtime/R-library')
        result = subprocess.run([os.environ.get('PROTEOMICS_RSCRIPT', 'Rscript'), str(ROOT/'pipeline_lib/drug_interaction.R'),
                                 str(job_path)], env=env, text=True, capture_output=True, timeout=60)
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        self.assertFalse((out/'cells/dda/peptide/interaction.csv').exists())
        self.assertEqual([r['method'] for r in read_csv(out/'analysis_summary.csv')], ['limma'])
        self.assertFalse(list(out.parent.glob('.drug_interaction_*')))


@unittest.skipUnless(shutil.which('node'), 'Node unavailable')
class InteractionBrowserTests(unittest.TestCase):
    def test_method_settings_load_and_roundtrip(self):
        html = (ROOT/'assets/config_ui.html').read_text()
        collect = html[html.index('const countThresholds='):html.index('function renderSamples()')]
        controls = html[html.index('function refreshDrugConditions'):html.index('let previousSearchModes=null;')]
        config = {'samples': fixtures.samples(), 'optional_analyses': {'drugs': fixtures.design()}}
        config['optional_analyses']['drugs']['interaction'].update(methods=METHODS,
            monte_carlo={'simulations': 299, 'seed': 18, 'workers': 2})
        harness = r'''
const assert=require('node:assert/strict'),vm=require('node:vm');
const elements=new Map();function $(id){if(!elements.has(id))elements.set(id,{value:'',checked:false,disabled:false,hidden:false,
 replaceChildren(){},add(){},setAttribute(){},classList:{toggle(){}},focus(){}});return elements.get(id);}
const context={$,Number,cfg:CONFIG,Option:function(text,value){this.text=text;this.value=value;},invalidateStep(){}};
vm.createContext(context);vm.runInContext(SOURCE,context);context.loadOptionalOptions();
let selected=context.collectOptionalOptions().drugs.interaction;
assert.equal(JSON.stringify(selected.methods),JSON.stringify(context.cfg.optional_analyses.drugs.interaction.methods));
assert.equal(JSON.stringify(selected.monte_carlo),JSON.stringify(context.cfg.optional_analyses.drugs.interaction.monte_carlo));
assert.equal($('interaction-mc-options').hidden,false);assert.equal($('interaction-mc-seed').disabled,false);
context.cfg.optional_analyses.drugs.interaction={...context.cfg.optional_analyses.drugs.interaction,methods:undefined,monte_carlo:undefined};
context.loadOptionalOptions();selected=context.collectOptionalOptions().drugs.interaction;
assert.equal(JSON.stringify(selected.methods),'["edgeR"]');assert.equal('monte_carlo' in selected,false);
assert.equal($('interaction-mc-options').hidden,true);assert.equal($('interaction-mc-seed').disabled,true);
'''.replace('SOURCE', json.dumps(collect+'\n'+controls)).replace('CONFIG', json.dumps(config))
        result = subprocess.run(['node', '-'], input=harness, text=True, capture_output=True, timeout=15)
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
