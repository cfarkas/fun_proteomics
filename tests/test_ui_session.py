import http.client,json,sys,tempfile,threading,unittest
from pathlib import Path
from http.server import ThreadingHTTPServer
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import config_server

class SessionTests(unittest.TestCase):
 def setUp(self):
  self.server=ThreadingHTTPServer(('127.0.0.1',0),config_server.Handler)
  self.thread=threading.Thread(target=self.server.serve_forever,daemon=True);self.thread.start()
 def tearDown(self):
  self.server.shutdown();self.server.server_close();self.thread.join()
 def request(self,method,path,headers=None,data=None):
  c=http.client.HTTPConnection('127.0.0.1',self.server.server_port)
  c.request(method,path,body=json.dumps(data or {}) if method=='POST' else None,headers=headers or {})
  r=c.getresponse();status=r.status;data=json.loads(r.read());c.close();return status,data
 def test_changed_session_rejects_before_action_and_can_refresh(self):
  with patch.object(config_server,'TOKEN','new-token'):
   status,data=self.request('POST','/api/not-an-action',{'X-Config-Token':'old-token'})
   self.assertEqual(status,403);self.assertEqual(data['code'],'TOKEN_STALE')
   status,data=self.request('GET','/api/session');self.assertEqual(status,200)
   status,_=self.request('POST','/api/not-an-action',{'X-Config-Token':data['token']})
   self.assertEqual(status,404)
 def test_foreign_origin_and_host_cannot_refresh(self):
  for headers in [{'Origin':'https://example.org'},{'Host':'example.org'},{'Sec-Fetch-Site':'cross-site'}]:
   status,_=self.request('GET','/api/session',headers);self.assertEqual(status,403)
 def test_saved_command_respects_continue_checkbox(self):
  with tempfile.TemporaryDirectory() as folder,patch.object(config_server,'ROOT',Path(folder)),patch.object(config_server,'check',side_effect=lambda cfg:cfg):
   for resume in (False,True,None):
    data={'config':{'outdir':'results'}}
    if resume is not None:data['resume']=resume
    status,result=self.request('POST','/api/save',{'X-Config-Token':config_server.TOKEN},data)
    self.assertEqual(status,200,result)
    self.assertEqual('--continue' in result['command'],resume is not False)
    self.assertTrue(Path(result['path']).is_file())
