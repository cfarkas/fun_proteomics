"""Read realistic NUL-terminated /proc command lines without touching live runs."""
import json
from pathlib import Path
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from ui_runner import Runner


class RunDetectionTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name) / 'project'
        self.root.mkdir()
        self.entry = Path(tmp.name) / '218628'
        self.entry.mkdir()
        (self.entry / 'cwd').symlink_to(self.root, target_is_directory=True)
        self.runner = Runner(self.root, command_prefix=['/bin/echo'])

    def command(self, arguments):
        (self.entry / 'cmdline').write_bytes(('\0'.join(arguments) + '\0').encode())

    def detected(self):
        with patch.object(Path, 'iterdir', return_value=iter([self.entry])):
            return self.runner.external_run()

    def test_interactive_launcher_is_ignored_with_trailing_flag(self):
        for prefix in [['nextflow'], ['java', '-jar', '/runtime/nextflow.jar']]:
            for flags in [['--interactively'], ['--interactively', 'true'], ['--interactively=true'],
                          ['--interactively', '-ansi-log', 'false']]:
                with self.subTest(prefix=prefix, flags=flags):
                    self.command(prefix + ['run', 'main.nf'] + flags)
                    self.assertIsNone(self.detected())

    def test_real_analysis_is_still_detected(self):
        for flags in [[], ['--interactively', 'false'], ['--interactively=false']]:
            with self.subTest(flags=flags):
                self.command(['java', '-jar', '/runtime/nextflow.jar', 'run', 'main.nf'] + flags)
                self.assertEqual(self.detected(), 218628)

    def test_absolute_workflow_matches_even_from_another_directory(self):
        (self.entry / 'cwd').unlink()
        (self.entry / 'cwd').symlink_to(self.root.parent, target_is_directory=True)
        self.command(['nextflow', 'run', str(self.root / 'main.nf')])
        self.assertEqual(self.detected(), 218628)
        self.command(['nextflow', 'run', 'another/main.nf', '--note', str(self.root / 'main.nf')])
        self.assertIsNone(self.detected())

    def test_completed_run_can_start_again_with_interactive_launcher_present(self):
        config = self.root / 'study.json'
        config.write_text(json.dumps({'outdir': 'results'}))
        self.command(['nextflow', 'run', 'main.nf', '--interactively'])
        for _ in range(2):
            with patch.object(Path, 'iterdir', return_value=iter([self.entry])):
                self.runner.start(config)
            deadline = time.monotonic() + 5
            while self.runner.status()['state'] == 'running' and time.monotonic() < deadline:
                time.sleep(0.01)
            self.assertEqual(self.runner.status()['state'], 'completed')
            state = self.runner.stop()
            self.assertEqual(state['state'], 'completed')
            self.assertEqual(state['exit_code'], 0)
            self.assertIn('already completed', state['message'])
            self.assertFalse(Path(self.runner.record['stop_path']).exists())

    def test_stop_without_an_active_run_is_harmless(self):
        self.assertEqual(self.runner.stop()['state'], 'idle')
        self.assertFalse(self.runner.index.exists())
