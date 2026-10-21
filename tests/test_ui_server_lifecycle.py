"""The setup server must outlive Ctrl+C in its launcher without starting analysis."""
import json
import os
from pathlib import Path
import signal
import socket
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from ui_server_lifecycle import server_info, start_server

SERVER = '''\
import argparse,json
from http.server import BaseHTTPRequestHandler,ThreadingHTTPServer
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument('--port',type=int)
a,_=p.parse_known_args()
class Handler(BaseHTTPRequestHandler):
    def log_message(self,*args):pass
    def do_GET(self):
        data=json.dumps({'root':str(Path.cwd()),'state':'idle'}).encode()
        self.send_response(200);self.send_header('Content-Length',str(len(data)))
        self.end_headers();self.wfile.write(data)
ThreadingHTTPServer(('127.0.0.1',a.port),Handler).serve_forever()
'''


def free_port():
    with socket.socket() as sock:
        sock.bind(('127.0.0.1', 0))
        return sock.getsockname()[1]


class ServerLifecycleTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        (self.root / 'scripts').mkdir()
        (self.root / 'scripts/config_server.py').write_text(SERVER)
        self.config = self.root / 'study.json'
        self.config.write_text('{}')
        self.port = free_port()

    def tearDown(self):
        self.temp.cleanup()

    def test_launcher_interrupt_keeps_server_and_reopening_reuses_it(self):
        code = '''
import json,sys,time
sys.path.insert(0,sys.argv[1])
from ui_server_lifecycle import start_server
url,log,process=start_server(sys.argv[2],sys.argv[3],int(sys.argv[4]))
print(json.dumps({'pid':process.pid,'url':url}),flush=True)
try:
    while True:time.sleep(1)
except KeyboardInterrupt:pass
'''
        parent = subprocess.Popen([sys.executable, '-c', code, str(ROOT / 'scripts'),
                                   str(self.root), str(self.config), str(self.port)],
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                  text=True, start_new_session=True)
        child_pid = None
        try:
            # Bound the readiness wait so a failed launch cannot hang the suite.
            import select
            self.assertTrue(select.select([parent.stdout], [], [], 15)[0])
            ready = json.loads(parent.stdout.readline())
            child_pid = ready['pid']
            self.assertEqual(os.getsid(child_pid), child_pid)
            self.assertNotEqual(os.getpgid(parent.pid), os.getpgid(child_pid))
            os.killpg(parent.pid, signal.SIGINT)
            parent.wait(timeout=5)
            self.assertEqual(parent.returncode, 0)
            self.assertEqual(server_info(ready['url'], self.root)['state'], 'idle')
            url, _, duplicate = start_server(self.root, self.config, self.port)
            self.assertEqual(url, ready['url'])
            self.assertIsNone(duplicate)
        finally:
            if parent.poll() is None:
                parent.kill();parent.wait(timeout=5)
            if child_pid:
                os.kill(child_pid, signal.SIGTERM)
            parent.stdout.close();parent.stderr.close()

    def test_port_from_another_project_is_rejected_without_replacing_server(self):
        url, _, process = start_server(self.root, self.config, self.port)
        try:
            with self.assertRaisesRegex(ValueError, 'another project'):
                start_server(self.root / 'other', self.config, self.port)
            self.assertEqual(server_info(url, self.root)['state'], 'idle')
        finally:
            process.terminate();process.wait(timeout=5)

    def test_startup_failure_reports_persistent_diagnostic_log(self):
        (self.root / 'scripts/config_server.py').write_text("print('fixture failure',flush=True)\nraise SystemExit(7)\n")
        with self.assertRaisesRegex(RuntimeError, 'exited with code 7'):
            start_server(self.root, self.config, self.port)
        self.assertIn('fixture failure', (self.root / f'.runtime/ui_server_{self.port}.log').read_text())
