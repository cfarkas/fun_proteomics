import hashlib,io,json,sys,tempfile,unittest,zipfile
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import tool_manager as tm
import analysis_modules as am
class ToolTests(unittest.TestCase):
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup);self.root=Path(self.tmp.name)
  self.enter=tm.ENTRIES['MSFragger'];self.archive=self.root/'MSFragger-4.4.1.zip'
  with zipfile.ZipFile(self.archive,'w') as z:z.writestr(self.enter,'verified content')
  specs=json.loads(json.dumps(tm.SPECS));specs['MSFragger']['archive_sha256']=hashlib.sha256(self.archive.read_bytes()).hexdigest()
  for obj,name,value in [(tm,'ROOT',self.root),(tm,'SPECS',specs),(tm,'roots',lambda extra='':[self.root])]:
   p=patch.object(obj,name,value);p.start();self.addCleanup(p.stop)
 def test_discovery_checks_archive_hash(self):
  self.assertEqual(tm.discover()['tools']['msfragger_zip'],str(self.archive))
  self.archive.write_bytes(b'wrong version')
  self.assertNotIn('msfragger_zip',tm.discover()['tools'])
 def test_force_reinstall_replaces_only_managed_installation(self):
  result=tm.install('MSFragger');entry=Path(result['path'])/self.enter
  entry.write_text('damaged')
  tm.install('MSFragger',force=True)
  self.assertEqual(entry.read_text(),'verified content')
 def test_registration_required_without_valid_archive(self):
  self.archive.unlink()
  with patch.object(tm.urllib.request,'urlopen') as fetch,self.assertRaisesRegex(ValueError,'register at'):
   tm.install('MSFragger',force=True)
  fetch.assert_not_called()
 def test_bad_download_does_not_replace_cached_archive(self):
  cache=self.root/'.runtime/downloads';cache.mkdir(parents=True)
  existing=cache/Path(tm.SPECS['FragPipe']['url']).name;existing.write_bytes(b'previous')
  with patch.object(tm.urllib.request,'urlopen',return_value=io.BytesIO(b'bad download')),self.assertRaisesRegex(ValueError,'checksum'):
   tm.ensure_archive('FragPipe',force=True)
  self.assertEqual(existing.read_bytes(),b'previous')
 def test_module_finds_archive_without_database_step(self):
  archive=self.root/'FragPipe.zip'
  with zipfile.ZipFile(archive,'w') as z:z.writestr('fragpipe-24.0/workflows/LFQ-phospho.workflow','ptmprophet.run-ptmprophet=true\nptmprophet.cmdline=STY\n')
  assets=self.root/'assets';assets.mkdir();(assets/'tool_versions.json').write_text(json.dumps({'FragPipe':{'version':'24.0','archive_sha256':hashlib.sha256(archive.read_bytes()).hexdigest()}}))
  with patch.object(am,'ROOT',self.root),patch.object(am,'TARGET',self.root/'modules/fragpipe.workflow'),patch.object(tm,'ensure_archive',return_value=archive) as ensure:
   self.assertTrue(am.install({},force=True)['installed'])
   ensure.assert_called_once_with('FragPipe',{},force=True)
