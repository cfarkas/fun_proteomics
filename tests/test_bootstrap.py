import json,os,shutil,subprocess,sys,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
class BootstrapTests(unittest.TestCase):
 def fixture(self,ready):
  tmp=tempfile.TemporaryDirectory();self.addCleanup(tmp.cleanup);p=Path(tmp.name)
  (p/'scripts').mkdir();(p/'runtime/env/bin').mkdir(parents=True)
  shutil.copy2(ROOT/'scripts/bootstrap_runtime.sh',p/'scripts/bootstrap_runtime.sh')
  (p/'runtime/env/bin/python').symlink_to(sys.executable)
  (p/'setup_environment.sh').write_text('#!/bin/bash\nexit 0\n')
  installer=p/'scripts/install_differential_tools.sh'
  installer.write_text('#!/bin/bash\ntouch "'+str(p/'ready')+'"\necho installed >> "'+str(p/'installs')+'"\n')
  r=p/'Rscript';r.write_text('#!/bin/bash\ntest -f "'+str(p/'ready')+'"\n');r.chmod(0o755)
  if ready:(p/'ready').touch()
  return p,r
 def launch(self,p,r):
  return subprocess.run(['bash',str(p/'scripts/bootstrap_runtime.sh'),str(p/'runtime'),str(r),str(p/'metadata.json'),'true'],capture_output=True,text=True)
 def test_reuses_working_r(self):
  p,r=self.fixture(True);result=self.launch(p,r)
  self.assertEqual(result.returncode,0,result.stderr)
  self.assertEqual(json.loads((p/'metadata.json').read_text())['rscript'],str(r))
  self.assertFalse((p/'installs').exists())
 def test_installs_missing_packages_once(self):
  p,r=self.fixture(False)
  for _ in range(2):
   result=self.launch(p,r);self.assertEqual(result.returncode,0,result.stderr)
  self.assertEqual((p/'installs').read_text().splitlines(),['installed'])
 def test_invalid_explicit_r_fails(self):
  p,r=self.fixture(True);result=self.launch(p,p/'missing-Rscript')
  self.assertNotEqual(result.returncode,0)
  self.assertIn('Rscript not found',result.stderr)
  self.assertFalse((p/'metadata.json').exists())
