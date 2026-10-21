"""DIA/DDA continuation separation, directory renames and publication checks."""
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from continue_run import select_session
from ui_runner import Runner
from ui_server_lifecycle import server_info
import validate_repo


class DiaIntegrationTests(unittest.TestCase):
    def test_browser_resume_ignores_dia_history(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for session in ('dda_session', 'dia_session'):
                (root / '.nextflow/cache' / session).mkdir(parents=True)
            (root / 'config.local.json').write_text('{}')
            (root / '.nextflow/history').write_text(
                'date\t1m\tdda\tERR\trev\tdda_session\tnextflow run main.nf\n'
                'date\t1m\tdia\tERR\trev\tdia_session\tnextflow run main.nf --dia_config dia.json\n')
            self.assertEqual(Runner(root).resume_candidate()['session'], 'dda_session')

    def test_dda_continue_ignores_dia_history(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            dda = '17ae0893-c3fe-4e41-9181-5ea3feefbc1b'
            dia = '52fa9dc5-c813-46e0-94b4-73a6c51ba57c'
            for session in (dda, dia):
                (root / '.nextflow/cache' / session).mkdir(parents=True)
            config = root / 'config.local.json'
            config.write_text('{"outdir":"results"}')
            (root / '.nextflow/history').write_text(
                f'date\t1m\tdda\tOK\trev\t{dda}\tnextflow run main.nf\n'
                f'date\t1m\tdia\tERR\trev\t{dia}\tnextflow run main.nf --dia_config=dia.json\n')
            self.assertEqual(select_session(root, root, config)['session'], dda)

    def test_setup_server_accepts_renamed_directory_alias(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / 'new-name'
            root.mkdir()
            alias = root.parent / 'old-name'
            alias.symlink_to(root, target_is_directory=True)
            response = io.BytesIO(json.dumps({'root': str(alias)}).encode())
            with patch('ui_server_lifecycle.urllib.request.urlopen', return_value=response):
                self.assertEqual(server_info('http://localhost:1234', root)['root'], str(alias))

    def test_setup_server_rejects_different_project(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            response = io.BytesIO(json.dumps({'root': str(root / 'different')}).encode())
            with patch('ui_server_lifecycle.urllib.request.urlopen', return_value=response):
                with self.assertRaisesRegex(ValueError, 'another project'):
                    server_info('http://localhost:1234', root)

    def test_publication_only_includes_source_module(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for relative in ('modules/diapasef/main.nf', 'modules/installed-tool/private.json',
                             'results/measurements.tsv', 'final_results_dated/measurements.tsv',
                             'config.local.json', 'README.md'):
                path = root / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text('fixture')
            (root / 'modules/diapasef/linked.py').symlink_to(root / 'config.local.json')
            with patch.object(validate_repo, 'ROOT', root):
                files = {str(p.relative_to(root)) for p in validate_repo.source_files(root)}
            self.assertEqual(files, {'README.md', 'modules/diapasef/main.nf'})



if __name__ == '__main__':
    unittest.main()
