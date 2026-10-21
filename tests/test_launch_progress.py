"""The live log reader must expose progress before tasks finish and skip old runs."""
import contextlib
import importlib.util
import io
import os
from pathlib import Path
import tempfile
import time
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location('stream_progress', Path(__file__).resolve().parents[1] / 'scripts/stream_progress.py')
progress = importlib.util.module_from_spec(spec)
spec.loader.exec_module(progress)


class ProgressTests(unittest.TestCase):
    def test_heartbeat_reports_total_silence_and_last_operation(self):
        with tempfile.TemporaryDirectory() as tmp:
            task = Path(tmp) / 'aa' / 'task'
            task.mkdir(parents=True)
            (task / '.command.begin').touch()
            (task / '.command.run').write_text("### name: 'FINAL_RESULTS'\n")
            log = task / '.command.out'
            log.write_text('[DIFFERENTIAL] peptides edgeRun A_vs_B\n')
            monitor = progress.Monitor(Path(tmp), time.time() - 1)
            output = io.StringIO()
            with patch.object(progress.time, 'monotonic', return_value=0) as timer, contextlib.redirect_stdout(output):
                monitor.poll()
                timer.return_value = 31
                monitor.poll()
                timer.return_value = 62
                monitor.poll()
                self.assertIn('no task log output for 1m 2s', output.getvalue())
                self.assertIn('Last output: [DIFFERENTIAL] peptides edgeRun A_vs_B', output.getvalue())
                log.write_text('[DIFFERENTIAL] peptides edgeRun A_vs_B\nnext operation\n')
                timer.return_value = 63
                monitor.poll()
                timer.return_value = 94
                monitor.poll()
                self.assertIn('for 31s. Last output: next operation', output.getvalue())

    def test_streams_before_exit_without_duplicate_output_and_drains_final_log(self):
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            task = work / 'aa' / 'test'
            (task / 'run').mkdir(parents=True)
            (task / '.command.begin').touch()
            (task / '.command.run').write_text("### name: 'FRAGPIPE (cell DDA+)'\n")
            log = task / 'run/pipeline.log'
            log.write_text('IonQuant running\n')
            monitor = progress.Monitor(work, time.time() - 1)
            output = io.StringIO()
            with contextlib.redirect_stdout(output):
                monitor.poll()
                self.assertIn('IonQuant running', output.getvalue())
                self.assertNotIn('[END]', output.getvalue())
                monitor.poll()
                self.assertEqual(output.getvalue().count('IonQuant running'), 1)
                with log.open('a') as handle:
                    handle.write('Quantification finished\n')
                (task / '.exitcode').write_text('0\n')
                monitor.poll()
            self.assertIn('Quantification finished', output.getvalue())
            self.assertIn('[END] FRAGPIPE (cell DDA+) | exit=0', output.getvalue())

    def test_skips_previous_run_logs(self):
        with tempfile.TemporaryDirectory() as tmp:
            task = Path(tmp) / 'aa' / 'old'
            task.mkdir(parents=True)
            marker = task / '.command.begin'
            marker.touch()
            os.utime(marker, (1, 1))
            output = io.StringIO()
            with contextlib.redirect_stdout(output):
                progress.Monitor(Path(tmp), time.time()).poll()
            self.assertEqual(output.getvalue(), '')

    def test_parent_lifecycle(self):
        self.assertTrue(progress.alive(os.getpid()))
        self.assertFalse(progress.alive(2147483647))

class VerbosityTests(unittest.TestCase):
    def test_levels_filter_logs_but_keep_warnings(self):
        outputs = {}
        for level in range(1, 5):
            with tempfile.TemporaryDirectory() as tmp:
                work = Path(tmp)
                task = work / 'aa' / 'test'
                (task / 'run').mkdir(parents=True)
                (task / '.command.begin').touch()
                (task / '.command.run').write_text("### name: 'FRAGPIPE (cell DDA+)'\n")
                (task / '.command.sh').write_text('python search_command.py\n')
                log = task / 'run/pipeline.log'
                log.write_text('Starting IonQuant\nintensity coefficient: 0.23\nWARN diagnostic\n')
                monitor = progress.Monitor(work, time.time() - 1, level)
                output = io.StringIO()
                with contextlib.redirect_stdout(output):
                    monitor.poll()
                outputs[level] = output.getvalue()
        self.assertNotIn('[START]', outputs[1])
        self.assertIn('[START]', outputs[2])
        self.assertNotIn('Starting IonQuant', outputs[2])
        self.assertIn('Starting IonQuant', outputs[3])
        self.assertNotIn('intensity coefficient', outputs[3])
        self.assertIn('intensity coefficient', outputs[4])
        self.assertIn('[COMMAND]', outputs[4])
        for output in outputs.values():
            self.assertIn('WARN diagnostic', output)

    def test_warning_split_across_reads_is_not_lost(self):
        with tempfile.TemporaryDirectory() as tmp:
            task = Path(tmp) / 'aa' / 'test'
            (task / 'run').mkdir(parents=True)
            (task / '.command.begin').touch()
            log = task / 'run/pipeline.log'
            log.write_text('WA')
            monitor = progress.Monitor(Path(tmp), time.time() - 1, 1)
            output = io.StringIO()
            with contextlib.redirect_stdout(output):
                monitor.poll()
                with log.open('a') as handle:
                    handle.write('RN split warning\n')
                monitor.poll()
            self.assertIn('WARN split warning', output.getvalue())


if __name__ == '__main__':
    unittest.main()
