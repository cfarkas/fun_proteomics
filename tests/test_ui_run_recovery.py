"""Exercise real subprocesses: losing the UI must not lose logs or run status."""
import contextlib
import io
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import threading
import time
import unittest
from unittest.mock import patch

SCRIPTS = Path(__file__).resolve().parents[1] / 'scripts'
sys.path.insert(0, str(SCRIPTS))
from ui_log_follow import mirror_runner
from ui_runner import Runner
from ui_run_state import process_identity, read_state, write_json


class RunRecoveryTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.config = self.root / 'study.json'
        self.config.write_text(json.dumps({'outdir': str(self.root / 'results')}))
        self.fake = self.root / 'analysis.py'
        self.runner = Runner(self.root, command_prefix=[sys.executable, str(self.fake)])
        self.addCleanup(self.stop_if_active)

    def until(self, predicate, timeout=10):
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            result = predicate()
            if result:
                return result
            time.sleep(0.02)
        self.fail('Timed out waiting for synthetic analysis')

    def stop_if_active(self):
        if self.runner.status()['state'] in ('running', 'stopping'):
            self.runner.stop()
            self.until(lambda: self.runner.status()['state'] not in ('running', 'stopping'))

    def start(self):
        with patch.object(self.runner, 'external_run', return_value=None):
            return self.runner.start(self.config)

    def write_waiting_analysis(self):
        self.fake.write_text("""from pathlib import Path
import time
print('analysis started', flush=True)
Path('started').touch()
deadline = time.monotonic() + 20
while not Path('release').exists() and time.monotonic() < deadline:
    time.sleep(0.02)
print('final output after interface loss', flush=True)
""")

    def test_killed_interface_keeps_logging_and_restores_completion(self):
        self.write_waiting_analysis()
        owner_code = """import sys, threading
sys.path.insert(0, sys.argv[1])
from ui_runner import Runner
from ui_log_follow import mirror_runner
from pathlib import Path
runner = Runner(sys.argv[2], command_prefix=[sys.executable, sys.argv[3]])
runner.start(sys.argv[4])
Path(sys.argv[2], 'owner_ready').touch()
mirror_runner(runner, threading.Event())
"""
        owner = subprocess.Popen([sys.executable, '-c', owner_code, str(SCRIPTS), str(self.root),
                                  str(self.fake), str(self.config)], stdout=subprocess.DEVNULL,
                                 stderr=subprocess.DEVNULL, start_new_session=True)
        try:
            self.until(lambda: (self.root / 'owner_ready').exists() and (self.root / 'started').exists())
            owner.kill()
            owner.wait(timeout=5)
            restored = Runner(self.root)
            self.assertEqual(restored.status()['state'], 'running')
            with self.assertRaisesRegex(ValueError, 'already active'):
                restored.start(self.config)
            (self.root / 'release').touch()
            self.until(lambda: restored.status()['state'] == 'completed')
            final = Runner(self.root).status()
            self.assertEqual(final['exit_code'], 0)
            self.assertEqual(final['config'], str(self.config))
            log = Path(final['log_path']).read_text()
            self.assertIn('analysis started', log)
            self.assertIn('final output after interface loss', log)
            self.assertIn('[RUN COMPLETED] exit 0', log)
        finally:
            if owner.poll() is None:
                owner.kill()
                owner.wait(timeout=5)

    def test_closed_terminal_does_not_affect_analysis(self):
        class ClosedTerminal(io.StringIO):
            def write(self, text):
                raise BrokenPipeError('terminal closed')
        self.write_waiting_analysis()
        self.start()
        self.until(lambda: (self.root / 'started').exists())
        with contextlib.redirect_stdout(ClosedTerminal()):
            mirror_runner(self.runner, threading.Event())
        (self.root / 'release').touch()
        self.until(lambda: self.runner.status()['state'] == 'completed')
        self.assertIn('final output after interface loss', self.runner.log_path.read_text())

    def test_failure_exit_code_is_restored(self):
        self.fake.write_text("import sys\nprint('analysis error', file=sys.stderr)\nsys.exit(7)\n")
        self.start()
        self.until(lambda: self.runner.status()['state'] == 'failed')
        state = Runner(self.root).status()
        self.assertEqual(state['exit_code'], 7)
        self.assertIn('analysis error', Path(state['log_path']).read_text())
        self.assertNotIn('[RUN COMPLETED]', Path(state['log_path']).read_text())

    def test_restored_stop_cleans_children_that_ignore_term(self):
        child_code = """import signal, time
from pathlib import Path
signal.signal(signal.SIGTERM, signal.SIG_IGN)
Path('child_ready').touch()
time.sleep(30)
"""
        self.fake.write_text(f"""import subprocess, sys, time
from pathlib import Path
child = subprocess.Popen([sys.executable, '-c', {child_code!r}])
Path('child_pid').write_text(str(child.pid))
time.sleep(30)
""")
        self.start()
        self.until(lambda: (self.root / 'child_ready').exists())
        child_pid = int((self.root / 'child_pid').read_text())
        self.assertIsNotNone(process_identity(child_pid))
        restored = Runner(self.root)
        self.assertEqual(restored.stop()['state'], 'stopping')
        self.assertEqual(restored.stop()['state'], 'stopping')
        with self.assertRaisesRegex(ValueError, 'already active'):
            restored.start(self.config)
        self.until(lambda: restored.status()['state'] == 'stopped')
        self.until(lambda: process_identity(child_pid) is None)
        self.assertEqual(Runner(self.root).status()['state'], 'stopped')
        self.assertNotEqual(restored.status()['exit_code'], 0)

    def test_missing_process_and_reused_pid_never_imply_success(self):
        state = self.root / 'lost.json'
        write_json(state, {'state': 'running', 'started': time.time() - 30,
                           'pid': os.getpid(), 'pid_identity': 'different-start-time',
                           'supervisor_pid': 2147483647, 'supervisor_identity': 'gone',
                           'stop_path': str(self.root / 'stop'), 'exit_code': None})
        actual = read_state(state)
        self.assertEqual(actual['state'], 'interrupted')
        self.assertIsNone(actual['exit_code'])


if __name__ == '__main__':
    unittest.main()
