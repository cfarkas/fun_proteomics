import unittest, shutil
from pathlib import Path
import pandas as pd
import test_differential
from pipeline_lib.result_reports import organize

class ReportTests(unittest.TestCase):
 run_model=test_differential.DifferentialTests.run_model
 def test_empty_feature_matrices_produce_valid_empty_reports(self):
  root=self.run_model(['CTRL']*3+['B']*3,feature_count=0,methods=['edgeR','limma'])
  out=root/'final_results';shutil.copy2(root/'job.json',out/'analysis_settings.json');organize(out)
  index=pd.read_csv(out/'differential/results_index.csv')
  self.assertTrue((index.tested_features==0).all())
  self.assertTrue((index.reported_features==0).all())
  table=pd.read_csv(out/'differential/edgeR/pairwise/test/dda/peptides/edgeR_CTRL_vs_B.csv')
  self.assertEqual(len(table),0)
  self.assertTrue({'FeatureID','FDR','test_status'}<=set(table.columns))
 def test_insufficient_evidence_tables_remain_linked_and_not_counted_as_tested(self):
  root=self.run_model(['CTRL']*3+['B']*3,zero_library=True,methods=['edgeR','limma'])
  out=root/'final_results';shutil.copy2(root/'job.json',out/'analysis_settings.json')
  organize(out)
  index=pd.read_csv(out/'differential/results_index.csv')
  self.assertTrue((index.tested_features==0).all())
  self.assertTrue((index.reported_features==120).all())
  self.assertTrue((index.significant_features==0).all())
  table=pd.read_csv(out/'differential/edgeR/pairwise/test/dda/peptides/edgeR_CTRL_vs_B.csv')
  self.assertTrue(table.FDR.isna().all())
  self.assertTrue(table.test_status.str.contains('insufficient evidence').all())
  self.assertIn('Low coverage is not a run failure',(out/'START_HERE.html').read_text())
 def test_reports_preserve_statistics_and_explain_conditions(self):
  root=self.run_model(['CTRL']*3+['B']*3+['C']*3,methods=['edgeR','edgeRun','limma'],modifications=True,labels={'CTRL':'Normal tissue'})
  out=root/'final_results'
  shutil.copy2(root/'job.json',out/'analysis_settings.json')
  before=pd.read_csv(out/'differential/edgeR/test/dda/edgeR_Normal_tissue_vs_B_genes.csv').set_index('FeatureID')
  self.assertEqual(organize(out),40)
  self.assertFalse(list((out/'differential').glob('*/significant_results.csv')))
  after=pd.read_csv(out/'differential/edgeR/pairwise/test/dda/genes/edgeR_Normal_tissue_vs_B.csv').set_index('FeatureID')
  pd.testing.assert_series_equal(before.adjusted_p_value,after.FDR,check_names=False)
  pd.testing.assert_series_equal(before.PValue,after.PValue)
  self.assertAlmostEqual(after.loc['P0','mean_count_B'],1012)
  self.assertEqual(after.loc['P0','detected_samples_B'],3)
  self.assertEqual(after.loc['P0','total_samples_B'],3)
  for method in ['edgeR','limma']:
   filename='edgeR_GLM_omnibus.csv' if method=='edgeR' else 'limma_omnibus.csv'
   frame=pd.read_csv(out/f'differential/{method}/omnibus/test/dda/genes/{filename}')
   self.assertFalse(any(c.lower().startswith('logfc') for c in frame))
   self.assertIn('mean_count_Normal_tissue',frame)
   self.assertIn('mean_normalized_CPM_C',frame)
  for path in (out/'differential').rglob('*.csv'):
   if path.name in ['results_index.csv','condition_columns.csv']:continue
   cols=pd.read_csv(path,nrows=0).columns
   self.assertIn('FDR',cols)
   self.assertFalse(set(cols)&{'adjusted_p_value','adj.P.Val','P.Value','qvalue','q_value'})
  self.assertTrue((out/'START_HERE.html').is_file())
  self.assertIn('Overall condition test',(out/'START_HERE.html').read_text())
