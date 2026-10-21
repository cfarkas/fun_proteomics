import contextlib,io,json,os,sys,tempfile,time,unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from ui_log_follow import LogMirror, follow
from ui_run_state import process_identity, write_json
class LogTests(unittest.TestCase):
 def test_existing_appended_and_next_run_logs(self):
  with tempfile.TemporaryDirectory() as tmp:
   out=io.StringIO();m=LogMirror(out);p=Path(tmp)/'run1.log';p.write_text('starting\n')
   m.update({'log_path':str(p)});m.update({'log_path':str(p)})
   with p.open('a') as h:h.write('progress\n')
   m.update({'log_path':str(p)})
   p=Path(tmp)/'run2.log';p.write_text('next run\n');m.update({'log_path':str(p)})
   self.assertEqual(out.getvalue(),'starting\nprogress\nnext run\n')

class FollowTests(unittest.TestCase):
 def test_temporary_disconnect_reconnects_without_duplicate_lines(self):
  with tempfile.TemporaryDirectory() as tmp:
   log=Path(tmp)/'run.log';log.write_text('running\n')
   status=json.dumps({'state':'running','log_path':str(log)})
   output=io.StringIO()
   replies=[io.StringIO(status),OSError('refused'),io.StringIO(status),KeyboardInterrupt()]
   with patch('ui_log_follow.urllib.request.urlopen',side_effect=replies),patch('ui_log_follow.time.sleep'),contextlib.redirect_stdout(output):
    self.assertEqual(follow('http://localhost'),0)
   self.assertIn('[INTERFACE RECONNECTED]',output.getvalue())
   self.assertEqual(output.getvalue().count('running\n'),1)

 def test_permanent_disconnect_without_saved_status_returns_failure(self):
  output=io.StringIO()
  with patch('ui_log_follow.urllib.request.urlopen',side_effect=OSError('refused')),contextlib.redirect_stdout(output):
   self.assertEqual(follow('http://localhost',retry_seconds=0),1)
  self.assertIn('completion could not be verified',output.getvalue())

 def test_follows_durable_log_to_success_or_failure_while_offline(self):
  for exit_code in (0,7):
   with self.subTest(exit_code=exit_code),tempfile.TemporaryDirectory() as tmp:
    root=Path(tmp);log=root/'run.log';log.write_text('started\n');state=root/'run.json'
    record={'state':'running','started':time.time(),'supervisor_pid':os.getpid(),
            'supervisor_identity':process_identity(os.getpid()),'log_path':str(log),
            'state_path':str(state),'stop_path':str(root/'stop'),'exit_code':None}
    write_json(state,record)
    def disconnected(*args,**kwargs):
     with log.open('a') as h:h.write('final lines\n')
     record.update(state='completed' if exit_code==0 else 'failed',exit_code=exit_code)
     write_json(state,record)
     raise OSError('refused')
    calls=iter([lambda:io.StringIO(json.dumps(record)),lambda:(_ for _ in ()).throw(OSError('refused')),disconnected])
    output=io.StringIO()
    with patch('ui_log_follow.urllib.request.urlopen',side_effect=lambda *a,**kw:next(calls)()),patch('ui_log_follow.time.sleep'),contextlib.redirect_stdout(output):
     self.assertEqual(follow('http://localhost'),0 if exit_code==0 else 1)
    self.assertEqual(output.getvalue().count('final lines'),1)
    self.assertIn(f'| exit {exit_code}',output.getvalue())
