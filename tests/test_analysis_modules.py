import json,sys,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import analysis_modules as m
ROOT=Path(__file__).resolve().parents[1]
class ModuleTests(unittest.TestCase):
 def test_install_integrity_and_tampering(self):
  config=ROOT/'config.local.json'
  if not config.exists(): self.skipTest('Local registered archive configuration unavailable')
  tools=json.loads(config.read_text())['tools']
  if not Path(tools.get('fragpipe_zip','')).is_file(): self.skipTest('Registered FragPipe archive unavailable')
  with tempfile.TemporaryDirectory() as d,patch.object(m,'TARGET',Path(d)/'module/fragpipe.workflow'):
   self.assertFalse(m.status()['installed'])
   result=m.install(tools);self.assertTrue(result['installed'])
   text=m.TARGET.read_text();self.assertIn('ptmprophet.run-ptmprophet=true',text);self.assertIn('STY',text);self.assertIn('workflow.input.data-type.im-ms=true',text)
   m.TARGET.write_text(text+'# changed');self.assertFalse(m.status()['installed'])
 def test_missing_archive(self):
  with patch('tool_manager.ensure_archive',side_effect=ValueError('Download unavailable')),self.assertRaisesRegex(ValueError,'Download unavailable'):m.install({})
