import csv,json,os,shutil,subprocess,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
class DifferentialTests(unittest.TestCase):
 def run_model(self,groups,control='CTRL',methods=None,enabled=True,modifications=False,labels=None,sparse=False,filter_counts=None,low_counts=False,thresholds=None,expected_error=None,analysis=None,zero_library=False,feature_count=120,gene_counts=None):
  tmp=tempfile.TemporaryDirectory();self.addCleanup(tmp.cleanup);root=Path(tmp.name);folder=root/'raw';folder.mkdir()
  samples=[(f'S{i}',g) for i,g in enumerate(groups)]
  with (folder/'samples.csv').open('w') as h:
   w=csv.writer(h);w.writerow(['sample_id','condition','replicate','cohort']);w.writerows((s,g,i%3+1,'test') for i,(s,g) in enumerate(samples))
  for kind in ['protein','peptide','modification','gene','modification_gene','modification_protein']:
   with (folder/(kind+'_counts.csv')).open('w') as h:
    w=csv.writer(h);w.writerow(['FeatureID','Protein','Gene','Modification']+[s for s,g in samples])
    for j in range(feature_count):
     values=[(1000+i*3 if j<10 and g=='B' else 30+(i*13+j*7)%30) for i,(s,g) in enumerate(samples)]
     if sparse and kind.startswith('modification_'):
      values=[(j*7+i*11)%5 if (j+i)%3==0 else 0 for i in range(len(samples))]
      if j<10:values=[0 if g!='B' else 2+i%2 for i,(s,g) in enumerate(samples)]
     if j==119:values=[0]*len(samples)
     if zero_library:values[0]=0
     if low_counts and j==118:values=[0 if g=='CTRL' else 1 for i,(s,g) in enumerate(samples)]
     if low_counts and j==119:values=[0]*len(samples)
     w.writerow([f'P{j}',f'P{j}',f'G{j}','Oxidation']+values)
  job={'groups':[{'cohort':'test','mode':'DDA','directory':str(folder)}],'outdir':str(root/'final_results'),'differential':{'enabled':enabled,'methods':methods or ['edgeR','limma'],'control':control,'edgeRun_upper':100000,'modifications':modifications,'condition_labels':labels or {}}}
  if gene_counts is not None:job['optional_analyses']={'gene_counts':gene_counts}
  if analysis:job['analysis']=analysis
  if filter_counts is not None:job['differential']['filter_low_counts']=filter_counts
  if thresholds:job['differential'].update(thresholds)
  path=root/'job.json';path.write_text(json.dumps(job));env=os.environ.copy();env['R_LIBS']=os.pathsep.join(filter(None,[str(ROOT/'.runtime/R-library'),env.get('R_LIBS')]))
  script=root/'differential.R';shutil.copy2(ROOT/'pipeline_lib/differential.R',script)
  result=subprocess.run([os.environ.get('PROTEOMICS_RSCRIPT','Rscript'),str(script),str(path)],env=env,capture_output=True,text=True)
  if expected_error:
   self.assertNotEqual(result.returncode,0)
   self.assertIn(expected_error,result.stdout+result.stderr)
   return root
  self.assertEqual(result.returncode,0,result.stdout+result.stderr)
  return root
 def rows(self,path):return list(csv.DictReader(path.read_text().splitlines()))
 def test_empty_library_keeps_unestimable_tables_and_other_outputs(self):
  root=self.run_model(['CTRL']*3+['B']*3,methods=['edgeR','limma'],modifications=True,zero_library=True)
  for category in ['proteins','genes','peptides','modifications_by_gene','modifications_by_protein']:
   for method in ['edgeR','limma']:
    rows=self.rows(root/f'final_results/differential/{method}/test/dda/{method}_CTRL_vs_B_{category}.csv')
    self.assertEqual(len(rows),120)
    self.assertTrue(all(r['PValue']=='' and r['adjusted_p_value']=='' and r['significant']=='FALSE' for r in rows))
    self.assertTrue(all('insufficient evidence' in r['test_status'] for r in rows))
  self.assertTrue((root/'final_results/matrices/test/dda/peptide_normalized_cpm.csv').is_file())
 def test_one_feature_study_finishes_with_explicit_evidence_status(self):
  root=self.run_model(['CTRL']*3+['B']*3,feature_count=1,methods=['edgeR','limma'])
  for method in ['edgeR','limma']:
   rows=self.rows(root/f'final_results/differential/{method}/test/dda/{method}_CTRL_vs_B_peptides.csv')
   self.assertEqual(len(rows),1)
   self.assertIn('test_status',rows[0])
 def test_disabled_filter_ignores_zero_thresholds_for_every_method(self):
  thresholds=dict.fromkeys(['min_count','min_total_count','modification_min_count','modification_min_total_count'],0)
  root=self.run_model(['CTRL']*4+['B']*4,methods=['edgeR','edgeRun','limma'],modifications=True,low_counts=True,filter_counts=False,thresholds=thresholds)
  for category in ['proteins','genes','peptides','modifications_by_gene','modifications_by_protein']:
   for method in ['edgeR','edgeRun','limma']:
    result=self.rows(root/f'final_results/differential/{method}/test/dda/{method}_CTRL_vs_B_{category}.csv')
    ids={row['FeatureID'] for row in result}
    self.assertIn('P118',ids)
    self.assertNotIn('P119',ids)
  self.assertTrue(all(row['filter_low_counts']=='FALSE' for row in self.rows(root/'final_results/analysis_summary.csv')))
 def test_hla_differential_tests_only_peptides(self):
  root=self.run_model(['CTRL']*3+['B']*3,methods=['edgeR','edgeRun','limma'],analysis='immunopeptidomics')
  self.assertEqual({r['FeatureType'] for r in self.rows(root/'final_results/analysis_summary.csv')},{'peptide'})
  for method in ['edgeR','edgeRun','limma']:
   self.assertTrue((root/f'final_results/differential/{method}/test/dda/{method}_CTRL_vs_B_peptides.csv').is_file())
   self.assertFalse((root/f'final_results/differential/{method}/test/dda/{method}_CTRL_vs_B_proteins.csv').exists())
 def test_active_filter_still_rejects_zero_thresholds(self):
  self.run_model(['CTRL']*3+['B']*3,filter_counts=True,thresholds={'min_count':0},expected_error='when low-count filtering is enabled')
  self.run_model(['CTRL']*3+['B']*3,filter_counts=True,modifications=True,thresholds={'modification_min_count':0},expected_error='when modification filtering is enabled')
 def test_two_group_glm_and_direction_all_feature_types(self):
  root=self.run_model(['CTRL']*3+['B']*3,methods=['edgeR','edgeRun','limma'],modifications=True)
  for category in ['proteins','genes','peptides','modifications_by_gene','modifications_by_protein']:
   for method in ['edgeR','edgeRun','limma']:
    rows=self.rows(root/f'final_results/differential/{method}/test/dda/{method}_CTRL_vs_B_{category}.csv')
    self.assertGreater(float(next(r for r in rows if r['FeatureID']=='P0')['logFC']),0)
    self.assertTrue(all(0<=float(r['adjusted_p_value'])<=1 for r in rows))
  self.assertTrue((root/'final_results/differential/edgeR/test/dda/edgeR_GLM_omnibus_genes.csv').is_file())
 def test_no_control_uses_glm_and_all_pairs(self):
  root=self.run_model(['A']*3+['B']*3+['C']*3)
  self.assertEqual(len(self.rows(root/'final_results/comparisons.csv')),3)
  self.assertTrue((root/'final_results/differential/limma/test/dda/limma_A_vs_C_peptides.csv').is_file())
  rows=self.rows(root/'zscore/test/dda/protein_zscore.csv')
  self.assertAlmostEqual(sum(float(v) for k,v in rows[0].items() if k!='FeatureID'),0,places=8)
 def test_disabled_skips_testing_but_keeps_normalization(self):
  root=self.run_model(['OnlyGroup']*3,enabled=False)
  self.assertFalse((root/'final_results/differential').exists())
  for kind in ['protein','peptide','modification','gene','modification_gene','modification_protein']:
   self.assertTrue((root/f'final_results/matrices/test/dda/{kind}_normalized_cpm.csv').is_file())
   self.assertTrue((root/f'zscore/test/dda/{kind}_zscore.csv').is_file())
 def test_any_reference_named_files_and_modification_tests(self):
  root=self.run_model(['Treatment']*3+['Normal tissue']*3,control='Treatment',modifications=True,labels={'Treatment':'TENI','Normal tissue':'NILO'})
  for category in ['proteins','genes','peptides','modifications_by_gene','modifications_by_protein']:
   for method in ['edgeR','limma']:
    rows=self.rows(root/f'final_results/differential/{method}/test/dda/{method}_TENI_vs_NILO_{category}.csv')
    self.assertTrue(all(r['baseline']=='Treatment' and r['comparison']=='Normal tissue' and r['direction']=='Normal tissue minus Treatment' for r in rows))
  self.assertFalse(list(root.rglob('*contrast_0*.csv')))
 def test_sparse_modifications_use_small_count_model(self):
  root=self.run_model(['CTRL']*4+['B']*4,modifications=True,sparse=True,filter_counts=True)
  rows=self.rows(root/'final_results/differential/edgeR/test/dda/edgeR_CTRL_vs_B_modifications_by_gene.csv')
  self.assertGreater(len(rows),0)
  self.assertTrue(all(r['normalization']=='TMMwsp' and 'small-count' in r['model'] for r in rows))
  self.assertNotIn('P119',[r['FeatureID'] for r in rows])
  report=self.rows(root/'final_results/analysis_summary.csv')
  sparse=next(r for r in report if r['FeatureType']=='modification_gene')
  self.assertGreater(float(sparse['zero_fraction']),.5)
  self.assertEqual(sparse['min_count'],'1')
  zs=self.rows(root/'zscore/test/dda/modification_gene_zscore.csv')
  for row in zs:
   values=[float(v) for k,v in row.items() if k!='FeatureID'];self.assertAlmostEqual(sum(values),0,places=8)
  self.assertTrue(all(float(v)==0 for k,v in next(r for r in zs if r['FeatureID']=='P119').items() if k!='FeatureID'))

 def test_optional_filter_keeps_low_counts_by_default(self):
  for flag in [None,False,True]:
   root=self.run_model(['CTRL']*4+['B']*4,modifications=True,low_counts=True,filter_counts=flag)
   for category in ['proteins','genes','peptides','modifications_by_gene','modifications_by_protein']:
    report={row['FeatureID']:row for row in self.rows(root/f'raw/models/{category}/counts_test_filter.csv')}
    self.assertEqual(report['P119']['reason'],'all_zero')
    self.assertEqual(report['P118']['tested'],'FALSE' if flag else 'TRUE')
    self.assertEqual(report['P118']['reason'],'low_count' if flag else 'retained')
    for method in ['edgeR','limma']:
     result=self.rows(root/f'final_results/differential/{method}/test/dda/{method}_CTRL_vs_B_{category}.csv')
     self.assertEqual('P118' in {row['FeatureID'] for row in result},not bool(flag))
     self.assertNotIn('P119',{row['FeatureID'] for row in result})

 def test_optional_gene_switch_leaves_protein_statistics_unchanged(self):
  enabled=self.run_model(['CTRL']*3+['B']*3,methods=['limma'],gene_counts=True)
  disabled=self.run_model(['CTRL']*3+['B']*3,methods=['limma'],gene_counts=False)
  protein='final_results/differential/limma/test/dda/limma_CTRL_vs_B_proteins.csv'
  self.assertEqual((enabled/protein).read_bytes(),(disabled/protein).read_bytes())
  gene='final_results/differential/limma/test/dda/limma_CTRL_vs_B_genes.csv'
  self.assertTrue((enabled/gene).is_file())
  self.assertFalse((disabled/gene).exists())
