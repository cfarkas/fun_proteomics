"""Completed DIA results must not launch a second search on a repeated click."""
import contextlib
import io
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
if not os.environ.get('FUN_PROTEOMICS_TEST_INSTALLED'):
    sys.path.insert(0, str(ROOT / 'src'))
from fun_proteomics import cli
from fun_proteomics.continuation import dia_input_snapshot, completed_dia_inputs_match
from fun_proteomics.project import digest


class CompletedDiaTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.raw = self.root / 'raw samples' / 'patient one.d'
        self.raw.mkdir(parents=True)
        (self.raw / 'analysis.tdf').write_bytes(b'metadata')
        (self.raw / 'analysis.tdf_bin').write_bytes(b'spectra')
        self.sheet = self.root / 'samples.tsv'
        self.sheet.write_text('sample_id\tsample_type\traw_path\tlibrary_reference\nS1\tstudy\tpatient one.d\ttrue\n')
        self.fasta = self.root / 'database.fasta'
        self.fasta.write_text('>P1\nPEPTIDE\n')
        self.method = self.root / 'method.toml'; self.method.write_text('test = 1\n')
        self.lock = self.root / 'requirements.txt'; self.lock.write_text('example==1\n')
        self.config = self.root / 'native.json'
        self.config.write_text(json.dumps(dict(samplesheet=str(self.sheet), raw_root=str(self.raw.parent),
            canonical_fasta=str(self.fasta), method_config=str(self.method), requirements_lock=str(self.lock))))
        self.output = self.root / 'results'
        final = self.output / 'final_results'; final.mkdir(parents=True)
        (final / 'matrix.csv').write_text('protein,S1\nP1,10\n')
        (final / 'FINAL_RESULTS_COMPLETE').write_text('complete\n')
        (final / 'OUTPUT_SHA256SUMS').write_text(''.join(digest(p) + '  ' + p.name + '\n' for p in sorted(final.iterdir())))
        controller = self.output / '.fun_proteomics'
        launch = controller / 'nextflow'; (launch / '.nextflow/cache/session').mkdir(parents=True)
        self.state = dict(engine='nextflow', status='completed', project=str(ROOT), outdir=str(self.output),
            config=str(self.config), config_sha256=digest(self.config), launch_dir=str(launch), session='session',
            command=['nextflow', 'run', str(ROOT / 'main.nf'), '--dia_config', str(self.config)],
            dia_inputs=dia_input_snapshot(self.config, ROOT))
        self.state_path = controller / 'run.json'
        self.state_path.write_text(json.dumps(self.state))

    def test_completed_dia_and_unrelated_folder_changes_do_not_start_nextflow(self):
        args = cli.parser().parse_args(['--continue', str(self.output)])
        before = self.state_path.read_bytes()
        with patch.object(cli, 'execute', side_effect=AssertionError('must not start a search')):
            for extra in (False, True):
                if extra:
                    (self.raw.parent / 'unrelated app folder').mkdir()
                    (self.raw.parent / 'launcher.sh').write_text('example')
                with contextlib.redirect_stdout(io.StringIO()) as output:
                    self.assertEqual(cli.continue_run(args, ROOT), 0)
                self.assertIn('No analysis needs repeating', output.getvalue())
        self.assertEqual(before, self.state_path.read_bytes())

    def test_changed_selected_inputs_or_damaged_outputs_resume_instead_of_skip(self):
        cases = [self.raw / 'analysis.tdf_bin', self.raw / 'analysis.tdf', self.fasta, self.sheet,
                 self.method, self.config, self.output / 'final_results/matrix.csv']
        args = cli.parser().parse_args(['--continue', str(self.output), '--print-command'])
        for path in cases:
            with self.subTest(path=path.name):
                original = path.read_bytes(); stat = path.stat()
                # Same-size raw rewrites still invalidate identity even if mtime is restored.
                if path.name == 'analysis.tdf_bin':
                    path.write_bytes(b'changed'); os.utime(path, ns=(stat.st_atime_ns, stat.st_mtime_ns))
                elif path == self.config:
                    data = json.loads(original); data['cpus'] = 2; path.write_text(json.dumps(data))
                elif path == self.sheet:
                    path.write_text(original.decode().replace('S1', 'S2'))
                else:
                    path.write_bytes(original + b'changed')
                with patch.object(cli, 'nextflow_binary', return_value='nextflow'), patch.object(cli, 'execute', return_value=0) as execute:
                    self.assertEqual(cli.continue_run(args, ROOT), 0)
                    self.assertEqual(execute.call_args.args[0][-2:], ['-resume', 'session'])
                path.write_bytes(original)
                self.state['dia_inputs'] = dia_input_snapshot(self.config, ROOT)
                self.state_path.write_text(json.dumps(self.state))

    def test_missing_baseline_and_interrupted_runs_cannot_skip(self):
        for change in ({'dia_inputs': None}, {'status': 'failed'}):
            self.assertFalse(completed_dia_inputs_match(dict(self.state, **change)))
        self.assertFalse(completed_dia_inputs_match(self.state, hash_raw=True))
        hashed = dict(self.state, dia_inputs=dia_input_snapshot(self.config, ROOT, True))
        self.assertTrue(completed_dia_inputs_match(hashed, hash_raw=True))

    def test_new_native_run_records_inputs_before_starting(self):
        target = self.root / 'new results'
        with patch.object(cli, 'resolve_project', return_value=ROOT), patch.object(cli, 'nextflow_binary', return_value='nextflow'), patch.object(cli, 'execute', return_value=0) as execute:
            self.assertEqual(cli.main(['run', '--dia-config', str(self.config), '--outdir', str(target)]), 0)
            self.assertEqual(execute.call_args.args[4]['dia_inputs'], self.state['dia_inputs'])


if __name__ == '__main__':
    unittest.main()
