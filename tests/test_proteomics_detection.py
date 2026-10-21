import sqlite3,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from proteomics_detection import inspect,detect_proteomics,confirm_config
class ProteomicsDetectionTests(unittest.TestCase):
 def make(self,root,name,types):
  p=root/name;p.mkdir();(p/'analysis.tdf_bin').touch()
  with sqlite3.connect(p/'analysis.tdf') as db:
   db.execute('CREATE TABLE Frames (MsMsType INTEGER)');db.executemany('INSERT INTO Frames VALUES (?)',[(t,) for t in types])
  return {'sample_id':name,'raw_path':str(p)}
 def test_metadata_not_filename(self):
  with tempfile.TemporaryDirectory() as tmp:
   r=Path(tmp);a=self.make(r,'misleading_DIA.d',[0,8,8]);b=self.make(r,'misleading_DDA.d',[0,9]);c=self.make(r,'ms1.d',[0])
   self.assertEqual(inspect(a['raw_path'])['acquisition'],'DDA-PASEF');self.assertEqual(inspect(b['raw_path'])['acquisition'],'DIA-PASEF');self.assertEqual(inspect(c['raw_path'])['acquisition'],'MS1-only')
   self.assertEqual(detect_proteomics([a,b])['detected'],'Mixed')
 def test_confirmation_and_supported_framework(self):
  with tempfile.TemporaryDirectory() as tmp:
   a=self.make(Path(tmp),'raw.d',[0,8]);c={'samples':[a]}
   with self.assertRaisesRegex(ValueError,'Confirm'):confirm_config(c)
   c['proteomics']={'acquisition':'DDA-PASEF','analysis':'global-label-free','confirmed':True};confirm_config(c)
   c['proteomics']['analysis']='glycoproteomics'
   with self.assertRaisesRegex(ValueError,'not implemented'):confirm_config(c)
 def test_corrupt_file_stays_unknown(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp)/'bad.d';p.mkdir();(p/'analysis.tdf_bin').touch();(p/'analysis.tdf').write_text('not sqlite')
   self.assertEqual(inspect(p)['acquisition'],'Unknown')
