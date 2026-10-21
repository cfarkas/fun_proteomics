import sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from ui_runner import Runner
class ResumeTests(unittest.TestCase):
 def test_history_selection(self):
  with tempfile.TemporaryDirectory() as t:
   root=Path(t);(root/'.nextflow/cache/session1').mkdir(parents=True);(root/'config.local.json').write_text('{}');h=root/'.nextflow/history'
   analysis='date\t1m\told_analysis\tERR\trev\tsession1\tnextflow run main.nf -w custom_work\n'
   h.write_text(analysis+'date\t1s\tcheck\tOK\trev\tsession2\tnextflow run main.nf --validate_only\n')
   r=Runner(root);c=r.resume_candidate();self.assertTrue(c['available']);self.assertEqual(c['session'],'session1');self.assertEqual(c['workdir'],'custom_work')
   h.write_text(analysis+'date\t1s\tcheck\tERR\trev\tsession2\tnextflow run main.nf --dry_run=true\n');self.assertEqual(r.resume_session(),'session1')
   h.write_text(analysis.replace('\tERR\t','\tOK\t'));self.assertFalse(r.resume_candidate()['available'])
   h.write_text(analysis);(root/'config.local.json').unlink();self.assertFalse(r.resume_candidate()['available'])
if __name__=='__main__':unittest.main()
