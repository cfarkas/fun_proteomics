"""Browser DIA conversion keeps study identity and never launches DDA settings."""
import copy
import csv
import hashlib
import http.client
import json
from pathlib import Path
import sqlite3
import sys
import tempfile
import threading
import time
import unittest
import zipfile
from http.server import ThreadingHTTPServer
from unittest.mock import Mock, patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import config_server
from dia_browser import browser_inputs, prepare_browser_study, validate_browser_study
from ui_runner import Runner
from ui_run_state import write_json


class _DiaFixture(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="dia browser's ")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / 'modules').mkdir()
        (self.root / 'modules/diapasef').symlink_to(ROOT / 'modules/diapasef', target_is_directory=True)
        (self.root / 'assets').symlink_to(ROOT / 'assets', target_is_directory=True)
        fasta = self.root / 'dda-search.fasta'
        fasta.write_text('>sp|P1|P1 Target\nMPEPTIDER\n>rev_sp|P1|P1 Target\nMEDITPEPR\n'
                         '>contam_KERATIN Contaminant\nMPEPTIDEK\n>rev_contam_KERATIN Contaminant\nMEDITPEPK\n')
        self.cfg = {'outdir': str(self.root / 'results'), 'modes': ['DDA+', 'DDA'],
                    'samples': [], 'tools': {'old_fragpipe': 'retained'},
                    'database': {'fasta': str(fasta), 'sha256': hashlib.sha256(fasta.read_bytes()).hexdigest(), 'contains_decoys': True},
                    'resources': {'cpus': 1, 'memory_gib': 4},
                    'differential': {'enabled': True, 'methods': ['edgeR']},
                    'optional_analyses': {'drugs': {'applied': True}},
                    'proteomics': {'acquisition': 'DIA-PASEF', 'analysis': 'global-label-free', 'confirmed': True},
                    'dia': {'confirmed': True, 'python': '/usr/bin/python3'}}
        for index in range(3):
            self.add_sample(index)
        self.path = self.root / 'study.json'
        self.save()

    def add_sample(self, index):
        raw = self.root / 'raw' / f'run {index:03}.d'
        raw.mkdir(parents=True)
        with sqlite3.connect(raw / 'analysis.tdf') as db:
            db.execute('CREATE TABLE Frames (MsMsType INTEGER)')
            db.executemany('INSERT INTO Frames VALUES (?)', [(0,), (9,)])
        (raw / 'analysis.tdf_bin').write_bytes(b'fixture only')
        self.cfg['samples'].append({'sample_id': f'Original_{index:03}', 'cohort': 'Serum',
                                    'condition': 'Drug A / vehicle', 'replicate': index+1, 'raw_path': str(raw)})

    def save(self):
        self.path.write_text(json.dumps(self.cfg))

class DiaBrowserTests(_DiaFixture):
    def test_preserves_study_and_names_strips_only_paired_decoys_and_reuses_inputs(self):
        before = self.path.read_bytes()
        native = prepare_browser_study(self.path, self.root)
        self.assertEqual(prepare_browser_study(self.path, self.root), native)
        self.assertEqual(self.path.read_bytes(), before)
        cfg = json.loads(native.read_text())
        self.assertNotIn('tools', cfg)
        self.assertEqual(cfg['differential'], {'enabled': False})
        fasta = Path(cfg['canonical_fasta']).read_text()
        self.assertIn('>contam_KERATIN', fasta)
        self.assertNotIn('>rev_', fasta)
        with Path(cfg['samplesheet']).open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        for original, row in zip(self.cfg['samples'], rows):
            self.assertEqual(row['sample_id'], original['sample_id'])
            self.assertEqual(row['condition'], original['condition'])
            self.assertEqual(row['sample_type'], original['cohort'])
        self.assertEqual(cfg['expected_samples'], 3)
        self.assertEqual(validate_browser_study(self.cfg, self.root)['library_references']['Serum'], [s['sample_id'] for s in self.cfg['samples']])
        self.assertFalse((self.root / 'results').exists())

    def test_library_fallback_is_bounded_deterministic_and_override_is_honored(self):
        for index in range(3, 25):self.add_sample(index)
        self.cfg['samples'].reverse()
        _, rows, _, _ = browser_inputs(self.cfg, self.root)
        self.assertEqual(sorted(r['sample_id'] for r in rows if r['library_reference']=='true'), [f'Original_{i:03}' for i in range(20)])
        self.cfg['samples'][0]['library_reference'] = True
        _, rows, _, _ = browser_inputs(self.cfg, self.root)
        self.assertEqual([r['sample_id'] for r in rows if r['library_reference']=='true'], ['Original_024'])

    def test_dia_acknowledgment_and_global_analysis_required(self):
        self.cfg['dia']['confirmed'] = False
        with self.assertRaisesRegex(ValueError, 'Confirm the DIA method'):browser_inputs(self.cfg, self.root)
        self.cfg['dia']['confirmed'] = True
        self.cfg['proteomics']['analysis'] = 'phosphoproteomics'
        with self.assertRaisesRegex(ValueError, 'global label-free'):browser_inputs(self.cfg, self.root)

    def test_fasta_checksum_and_unrecognized_decoys_rejected(self):
        self.cfg['database']['sha256'] = '0'*64
        with self.assertRaisesRegex(ValueError, 'SHA256'):browser_inputs(self.cfg, self.root)
        fasta = Path(self.cfg['database']['fasta'])
        fasta.write_text('>decoy_P1\nMPEPTIDER\n')
        self.cfg['database']['sha256'] = hashlib.sha256(fasta.read_bytes()).hexdigest()
        with self.assertRaisesRegex(ValueError, 'decoy metadata'):browser_inputs(self.cfg, self.root)

    def test_existing_contaminant_is_not_duplicated_by_archive(self):
        archive = self.root / 'contaminants.zip'
        with zipfile.ZipFile(archive, 'w') as handle:
            handle.writestr('crap-gpmdb.fas', '>KERATIN Contaminant\nMPEPTIDEK\n')
        self.cfg['database']['contaminants_zip'] = str(archive)
        _, _, targets, _ = browser_inputs(self.cfg, self.root)
        self.assertEqual([header for header, _ in targets].count('contam_KERATIN Contaminant'), 1)
        self.assertEqual(len(targets), 2)

    def test_invalid_python_executable_reports_configuration_error(self):
        self.cfg['dia']['python'] = '/bin/false'
        with self.assertRaisesRegex(ValueError, 'selected DIA Python executable'):
            validate_browser_study(self.cfg, self.root)

    def test_native_frame_validation_is_strict_even_with_manual_override(self):
        raw = Path(self.cfg['samples'][0]['raw_path'])
        with sqlite3.connect(raw / 'analysis.tdf') as db:db.execute('UPDATE Frames SET MsMsType=8 WHERE MsMsType=9')
        self.cfg['proteomics']['manual_override'] = True
        with self.assertRaisesRegex(ValueError, 'expected DIA-PASEF'):validate_browser_study(self.cfg, self.root)

    def test_relative_raw_paths_use_explicit_raw_root(self):
        self.cfg['dia']['raw_root'] = str(self.root / 'raw')
        for sample in self.cfg['samples']:sample['raw_path'] = Path(sample['raw_path']).name
        cfg, rows, _, _ = browser_inputs(self.cfg, self.root)
        self.assertEqual(cfg['raw_root'], str(self.root / 'raw'))
        self.assertTrue(all(Path(row['raw_path']).is_absolute() for row in rows))

    def test_backend_dia_validation_does_not_require_dda_tools_or_R(self):
        before = copy.deepcopy(self.cfg)
        with patch.object(config_server, 'ROOT', self.root), patch.object(config_server, 'capabilities', return_value={'available_cpus': 4, 'max_memory_gib': 32}), patch.object(config_server, 'differential_capabilities', side_effect=AssertionError('DIA must not check DDA R tools')):
            self.assertEqual(config_server.check(self.cfg), before)
        self.assertFalse((self.root / 'ui_configs').exists())

    def test_backend_resolves_relative_paths_before_metadata_check_and_save(self):
        self.cfg['dia']['raw_root'] = 'raw'
        self.cfg['database']['fasta'] = 'dda-search.fasta'
        self.cfg['outdir'] = 'results'
        for sample in self.cfg['samples']:sample['raw_path'] = Path(sample['raw_path']).name
        with patch.object(config_server, 'ROOT', self.root), patch.object(config_server, 'capabilities', return_value={'available_cpus': 4, 'max_memory_gib': 32}):
            checked = config_server.check(self.cfg)
        self.assertTrue(all(Path(sample['raw_path']).is_absolute() for sample in checked['samples']))
        self.assertEqual(checked['database']['fasta'], str(self.root / 'dda-search.fasta'))
        self.assertEqual(checked['outdir'], str(self.root / 'results'))

    def test_runner_uses_native_dia_flag_and_resume_retains_browser_config(self):
        native = prepare_browser_study(self.path, self.root)
        runner = Runner(self.root, command_prefix=['/bin/echo'])
        with patch.object(runner, 'external_run', return_value=None):runner.start(self.path, dia_config=native, resume=True, resume_session='fixture-session')
        for _ in range(200):
            if runner.status()['state']!='running':break
            time.sleep(.01)
        state = runner.status()
        self.assertEqual(state['state'], 'completed')
        self.assertEqual(state['config'], str(self.path))
        self.assertIn('--dia_config', state['command'])
        self.assertNotIn('--config', state['command'])
        self.assertNotIn('--continue', state['command'])
        self.assertEqual(state['command'][state['command'].index('-resume')+1], 'fixture-session')
        with self.assertRaisesRegex(ValueError, 'unavailable for DIA'):runner.start(self.path, report_only=True)

    def test_managed_dia_resume_selects_its_own_native_session(self):
        native = prepare_browser_study(self.path, self.root)
        runner = Runner(self.root)
        session='dia-session'
        (self.root / '.nextflow/cache' / session).mkdir(parents=True)
        (self.root / '.nextflow/history').write_text('date\t1m\tDIA\tERR\trev\t'+session+'\tnextflow run main.nf --dia_config '+__import__('shlex').quote(str(native))+'\n')
        state = self.root / 'state.json'
        write_json(state, {'workflow': 'dia', 'kind': 'analysis', 'state': 'failed', 'analysis_config': str(native), 'config': str(self.path)})
        write_json(runner.index, {'state_path': str(state)})
        candidate = runner.resume_candidate()
        self.assertTrue(candidate['available'])
        self.assertEqual(candidate['session'], session)
        self.assertEqual(candidate['dia_config'], str(native))


class DiaBrowserApiTests(_DiaFixture):
    def setUp(self):
        super().setUp()
        self.server = ThreadingHTTPServer(('127.0.0.1', 0), config_server.Handler)
        self.server.config = self.path
        self.server.runner = Mock()
        self.thread = threading.Thread(target=self.server.serve_forever, daemon=True)
        self.thread.start()
        self.root_patch = patch.object(config_server, 'ROOT', self.root)
        self.root_patch.start()
        self.addCleanup(self.root_patch.stop)
        self.addCleanup(self.close_server)

    def close_server(self):
        self.server.shutdown();self.server.server_close();self.thread.join()

    def request(self, method, path, data=None):
        client = http.client.HTTPConnection('127.0.0.1', self.server.server_port)
        client.request(method, path, body=json.dumps(data) if data is not None else None,
                       headers={'X-Config-Token': config_server.TOKEN})
        response = client.getresponse();result = response.status, json.loads(response.read());client.close()
        return result

    def test_draft_accepts_incomplete_config_keeps_prior_file_and_roundtrips_ui(self):
        before = self.path.read_bytes()
        draft = {'samples': [{'sample_id': 'Unsaved_name'}]}
        status, result = self.request('POST', '/api/draft', {'config': draft, 'ui_state': {'currentStep': 2, 'fields': {'input': 'unfinished'}}})
        self.assertEqual(status, 200, result)
        self.assertEqual(self.path.read_bytes(), before)
        status, loaded = self.request('GET', '/api/config')
        self.assertEqual(status, 200)
        self.assertEqual(loaded['config'], draft)
        self.assertEqual(loaded['ui_state']['currentStep'], 2)
        self.server.runner.start.assert_not_called()

    def test_save_and_run_preserve_complete_ui_field_state(self):
        before = self.path.read_bytes()
        ui_state = {'step': 7, 'fields': {
            'analysis-phospho': {'value': 'phosphoproteomics', 'checked': True},
            'drug-study-type': {'value': 'synergy', 'checked': False},
            'dia-python': {'value': '/usr/bin/python3', 'checked': False}}}
        self.server.runner.start.return_value = {'state': 'running', 'command': ['fixture-only']}
        with patch.object(config_server, 'check', side_effect=lambda cfg:cfg):
            for endpoint in ('/api/save', '/api/run'):
                with self.subTest(endpoint=endpoint):
                    status, result = self.request('POST', endpoint, {'config': self.cfg, 'ui_state': ui_state})
                    self.assertEqual(status, 200, result)
                    saved = Path(result['path'])
                    self.assertEqual(json.loads(saved.with_suffix('.ui.json').read_text()), ui_state)
                    status, loaded = self.request('GET', '/api/config')
                    self.assertEqual(status, 200)
                    self.assertEqual(loaded['ui_state'], ui_state)
                    self.assertEqual(json.loads(saved.read_text()), self.cfg)
                    self.assertEqual(loaded['config']['samples'], self.cfg['samples'])
        self.assertEqual(self.path.read_bytes(), before)
        self.assertEqual(self.server.runner.start.call_count, 1)

    def test_invalid_ui_state_is_rejected_before_saving_or_running(self):
        for endpoint in ('/api/save', '/api/run'):
            with self.subTest(endpoint=endpoint):
                status, result = self.request('POST', endpoint, {'config': self.cfg, 'ui_state': []})
                self.assertEqual(status, 400, result)
                self.assertIn('UI state must be an object', result['error'])
        self.assertFalse((self.root / 'ui_configs').exists())
        self.server.runner.start.assert_not_called()

    def test_saved_dia_command_and_library_preview(self):
        with patch.object(config_server, 'check', side_effect=lambda cfg:cfg):
            status, result = self.request('POST', '/api/save', {'config': self.cfg, 'resume': True})
        self.assertEqual(status, 200, result)
        self.assertTrue(result['command'].startswith('nextflow run main.nf --dia_config '))
        self.assertIn(' -resume', result['command'])
        self.assertNotIn('--continue', result['command'])
        self.assertEqual(json.loads(Path(result['path']).read_text()), self.cfg)
        self.assertEqual(len(result['dia']['library_references']['Serum']), 3)
        self.server.runner.start.assert_not_called()


if __name__ == '__main__':
    unittest.main()
