import json,sys,tempfile,time,unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from ui_runner import Runner

class RunLogTests(unittest.TestCase):
 def check_logs(self,relative):
  with tempfile.TemporaryDirectory() as t:
   root=Path(t);configs=root/'studies';configs.mkdir()
   output=configs/'chosen results' if relative else root/'chosen results'
   config=configs/'study.json';config.write_text(json.dumps({'outdir':'chosen results' if relative else str(output)}))
   runner=Runner(root,command_prefix=['/bin/echo'])
   with patch.object(runner,'external_run',return_value=None):runner.start(config)
   for _ in range(100):
    if runner.status()['state']!='running':break
    time.sleep(.01)
   self.assertEqual(runner.status()['state'],'completed')
   self.assertEqual(runner.log_path.parent,output/'logs')
   self.assertTrue(runner.log_path.is_file())
   for flag in ['-log','-with-trace','-with-report','-with-timeline']:
    self.assertEqual(Path(runner.command[runner.command.index(flag)+1]).parent,output/'logs')
   self.assertLess(runner.command.index('-log'),runner.command.index('run'))
   self.assertFalse((root/'logs').exists())
 def test_absolute_output(self):self.check_logs(False)
 def test_output_relative_to_study(self):self.check_logs(True)
