import os,sys,tempfile,unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import native_folders as native
class NativeFolderTests(unittest.TestCase):
 def invoke(self,initial,target='input',code=0,output=None):
  with patch.dict(os.environ,{'DISPLAY':':1'}),patch.object(native.shutil,'which',side_effect=lambda n:'/usr/bin/zenity' if n=='zenity' else None),patch.object(native.subprocess,'run',return_value=SimpleNamespace(returncode=code,stdout=output or initial)) as run:
   result=native.choose_folder(initial,target)
   return result,run.call_args.args[0]
 def test_folder_with_spaces_and_create_location(self):
  with tempfile.TemporaryDirectory(prefix='folder with spaces ') as folder:
   result,command=self.invoke(folder,'output')
   self.assertEqual(result,{'cancelled':False,'path':folder})
   self.assertIn('--directory',command)
   self.assertIn('--filename='+folder+'/',command)
 def test_cancel_does_not_change_path(self):
  result,_=self.invoke('/tmp',code=1)
  self.assertEqual(result,{'cancelled':True})
 def test_nonexistent_initial_opens_existing_parent(self):
  with tempfile.TemporaryDirectory() as folder:
   result,command=self.invoke(folder+'/not/created','output',output=folder)
   self.assertIn('--filename='+folder+'/',command)
 def test_headless_has_clear_error(self):
  with patch.dict(os.environ,{},clear=True),self.assertRaisesRegex(ValueError,'No desktop'):
   native.choose_folder('/tmp')
 def test_failure_releases_lock(self):
  with self.assertRaisesRegex(ValueError,'Could not open'):self.invoke('/tmp',code=2)
  self.assertFalse(native._DIALOG_LOCK.locked())
 def test_concurrent_dialog_rejected(self):
  native._DIALOG_LOCK.acquire()
  try:
   with self.assertRaisesRegex(ValueError,'already open'):self.invoke('/tmp')
  finally:native._DIALOG_LOCK.release()

 def test_tools_folder_has_own_title(self):
  result,command=self.invoke('/tmp','tools')
  self.assertEqual(result['path'],'/tmp')
  self.assertIn('--title=Choose tools folder',command)

 def test_native_file_selection_preserves_spaces_and_checks_file_type(self):
  with tempfile.TemporaryDirectory(prefix='native files ') as folder:
   path=Path(folder)/'protein matrix.csv';path.write_text('sample,value\nA,1\n')
   with patch.dict(os.environ,{'DISPLAY':':1'}),patch.object(native.shutil,'which',side_effect=lambda n:'/usr/bin/zenity' if n=='zenity' else None),patch.object(native.subprocess,'run',return_value=SimpleNamespace(returncode=0,stdout=str(path))) as run:
    self.assertEqual(native.choose_file(str(path)),{'cancelled':False,'path':str(path)})
    self.assertNotIn('--directory',run.call_args.args[0])
    self.assertIn('--filename='+str(path),run.call_args.args[0])
    run.return_value.stdout=folder
    with self.assertRaisesRegex(ValueError,'wrong type'):native.choose_file(folder)
   self.assertFalse(native._DIALOG_LOCK.locked())

 def test_file_dialog_supports_kdialog(self):
  with tempfile.NamedTemporaryFile() as file:
   with patch.dict(os.environ,{'DISPLAY':':1'}),patch.object(native.shutil,'which',side_effect=lambda n:'/usr/bin/kdialog' if n=='kdialog' else None),patch.object(native.subprocess,'run',return_value=SimpleNamespace(returncode=0,stdout=file.name)) as run:
    self.assertEqual(native.choose_file(file.name)['path'],file.name)
    self.assertIn('--getopenfilename',run.call_args.args[0])
