"""Guided setup must recognize raw data and survive its launcher exiting."""
import json
import os
from pathlib import Path
import signal
import sqlite3
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
if not os.environ.get('FUN_PROTEOMICS_TEST_INSTALLED'):
    sys.path.insert(0, str(ROOT / 'src'))
from fun_proteomics import cli
from fun_proteomics.setup_web import Application
from fun_proteomics.setup_session import start, reachable


class GuidedSessionTests(unittest.TestCase):
    def test_raw_working_folder_is_detected_without_overriding_explicit_inputs(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            for name in ('C1', 'C2', 'BLANK'):
                raw = root / (name + '.d')
                raw.mkdir()
                (raw / 'analysis.tdf_bin').write_bytes(b'fixture')
                with sqlite3.connect(raw / 'analysis.tdf') as db:
                    db.execute('CREATE TABLE Frames (MsMsType INTEGER)')
                    db.executemany('INSERT INTO Frames VALUES (?)', [(0,), (9,)])
            with patch('pathlib.Path.cwd', return_value=root):
                app = Application(cli.parser().parse_args(['setup']), ROOT)
                self.assertEqual(app.initial['mode'], 'dia')
                self.assertEqual(app.initial['raw_root'], str(root))
                for flags in (['--mode', 'matrix'], ['--matrix', '/example/table.csv']):
                    app = Application(cli.parser().parse_args(['setup', *flags]), ROOT)
                    self.assertEqual(app.initial['mode'], 'matrix')
                    self.assertNotIn('raw_root', app.initial)
            # Explicit raw-root works even when launched elsewhere.
            app = Application(cli.parser().parse_args(['setup', '--raw-root', str(root)]), ROOT)
            self.assertEqual(app.initial['mode'], 'dia')
            with sqlite3.connect(root / 'C2.d/analysis.tdf') as db:
                db.execute('UPDATE Frames SET MsMsType=8 WHERE MsMsType=9')
            app = Application(cli.parser().parse_args(['setup', '--raw-root', str(root)]), ROOT)
            self.assertEqual(app.initial['mode'], 'matrix')

    def test_server_survives_launcher_hangup_and_reopening_reuses_session(self):
        with tempfile.TemporaryDirectory(prefix='guided-session-') as temp:
            env = os.environ.copy()
            if not os.environ.get('FUN_PROTEOMICS_TEST_INSTALLED'):
                env['PYTHONPATH'] = str(ROOT / 'src')
            else:
                env.pop('PYTHONPATH', None)
            code = """import json,sys,time
from pathlib import Path
from fun_proteomics import cli
from fun_proteomics.setup_session import start
print(json.dumps(start(cli.parser().parse_args(['setup']), Path(sys.argv[1]))), flush=True)
time.sleep(120)
"""
            parent = subprocess.Popen([sys.executable, '-c', code, str(ROOT)], cwd=temp, env=env,
                                      stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, start_new_session=True)
            record = None
            try:
                import select
                self.assertTrue(select.select([parent.stdout], [], [], 30)[0])
                line = parent.stdout.readline()
                self.assertTrue(line, parent.stderr.read() if parent.poll() is not None else '')
                record = json.loads(line)
                self.assertEqual(os.getsid(record['pid']), record['pid'])
                os.killpg(parent.pid, signal.SIGHUP)
                parent.wait(timeout=5)
                self.assertTrue(reachable(record))
                with patch('pathlib.Path.cwd', return_value=Path(temp)):
                    reopened = start(cli.parser().parse_args(['setup']), ROOT)
                self.assertEqual(reopened['url'], record['url'])
                self.assertEqual(reopened['pid'], record['pid'])
                for _ in range(3):
                    request = Request(record['url'].split('/#')[0] + '/api/browse',
                                      json.dumps({'path': temp}).encode(),
                                      headers={'X-Session-Token': record['token'], 'Content-Type': 'application/json'})
                    with urlopen(request, timeout=5) as response:
                        self.assertEqual(json.load(response)['path'], temp)
            finally:
                if parent.poll() is None:
                    parent.terminate(); parent.wait(timeout=5)
                parent.stdout.close(); parent.stderr.close()
                if record:
                    os.kill(record['pid'], signal.SIGTERM)


if __name__ == '__main__':
    unittest.main()
