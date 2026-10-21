from pathlib import Path
import sys,tempfile,unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from sample_detection import detect
class DetectionTests(unittest.TestCase):
 def raw(self,root,name):
  p=root/name;p.mkdir(parents=True)
  for f in ['analysis.tdf','analysis.tdf_bin']:(p/f).touch()
  return p
 def test_discovery_inference_and_single_directory(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);a=self.raw(root,'cells/DrugA_rep2.d');self.raw(root,'cells/DrugA_rep1.d');(root/'bad.d').mkdir()
   result=detect(root,recursive=True);self.assertEqual(len(result['samples']),2);self.assertEqual(len(result['skipped']),1)
   self.assertEqual([(s['condition'],s['replicate'],s['cohort']) for s in result['samples']],[('DrugA',1,'cell'),('DrugA',2,'cell')])
   self.assertEqual(len(detect(a)['samples']),1);self.assertEqual(detect(root,recursive=False)['samples'],[])
 def test_preserve_user_assignments(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);p=self.raw(root,'C1_Control.d');row={'sample_id':'Mine','condition':'Vehicle','replicate':5,'cohort':'cell','raw_path':str(p)}
   self.assertEqual(detect(root,[row])['samples'],[row])
 def test_instrument_suffix_and_unknown_cohort(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);self.raw(root,'DDA_PASEF_EVOSEP_40SPDWZ_CF_TN3_Dual treatment_S1-B3_1_26089.d')
   row=detect(root)['samples'][0];self.assertEqual((row['condition'],row['replicate'],row['cohort']),('Dual treatment',3,'unassigned'))
 def test_unique_ids_and_incomplete_folder(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);self.raw(root,'first/Drug_rep1.d');self.raw(root,'second/Drug_rep1.d')
   rows=detect(root,recursive=True)['samples'];self.assertEqual(len({s['sample_id'] for s in rows}),2)
 def test_selected_folder_only_and_empty_folder_clears_existing(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);direct=self.raw(root,'Drug_rep1.d');self.raw(root,'results/copied/Drug_rep2.d')
   (root/'mock').mkdir()
   rows=detect(root)['samples'];self.assertEqual([s['raw_path'] for s in rows],[str(direct)])
   self.assertEqual(detect(root/'mock',existing=rows)['samples'],[])
   self.assertEqual(len(detect(root,recursive=True)['samples']),2)
 def test_blank_or_relative_path_never_scans_working_directory(self):
  for value in ['', '  ', '.', 'mock', None]:
   with self.subTest(value=value),self.assertRaises(ValueError):detect(value)
 def test_other_ms_formats_are_reported_without_becoming_bruker_samples(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);self.raw(root,'cells_rep1.d')
   for name in ['Thermo.RAW','sample.mzML','sample.mzXML','sample.wiff','sample.wiff2']:(root/name).touch()
   (root/'Waters.raw').mkdir();self.raw(root,'Waters.raw/hidden_rep1.d')
   result=detect(root,recursive=True)
   self.assertEqual(len(result['samples']),1)
   self.assertEqual(len(result['unsupported']),6)
   self.assertEqual(result['skipped'],[])
if __name__=='__main__':unittest.main()
