"""Folder navigation and discovery must work before an output is selected."""
import json
import os
from pathlib import Path
import sys
import tempfile
import threading
import unittest
from urllib.error import HTTPError
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
if not os.environ.get('FUN_PROTEOMICS_TEST_INSTALLED'):
    sys.path.insert(0, str(ROOT / 'src'))
from fun_proteomics import cli, study
from fun_proteomics.setup_web import Application, make_server


class FolderSetupTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='setup-folders-')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.app = Application(cli.parser().parse_args(['setup']), ROOT)
        self.server = make_server(self.app)
        thread = threading.Thread(target=self.server.serve_forever, daemon=True)
        thread.start()
        self.addCleanup(thread.join)
        self.addCleanup(self.server.server_close)
        self.addCleanup(self.server.shutdown)

    def post(self, route, data, token=None):
        url = f'http://127.0.0.1:{self.server.server_port}/api/{route}'
        headers = {'Content-Type': 'application/json',
                   'X-Session-Token': self.app.token if token is None else token}
        with urlopen(Request(url, json.dumps(data).encode(), headers), timeout=15) as response:
            return json.load(response)

    def test_search_and_folder_filter_precede_the_display_limit(self):
        for index in range(505):
            (self.root / f'folder-{index:03}').mkdir()
        (self.root / 'Z-last folder').mkdir()
        (self.root / 'Z-last file.csv').touch()
        listing = self.post('browse', {'path': str(self.root), 'kind': 'directory'})
        self.assertTrue(listing['truncated'])
        self.assertEqual(listing['total'], 506)
        result = self.post('browse', {'path': str(self.root), 'kind': 'directory', 'query': 'z-LAST'})
        self.assertEqual([v['name'] for v in result['entries']], ['Z-last folder'])
        self.assertFalse(result['truncated'])
        self.assertEqual(result['breadcrumbs'][-1]['path'], str(self.root))
        self.assertTrue(any(p['name'] == 'Home' for p in result['locations']))
        result = self.post('browse', {'path': str(self.root / 'new' / 'project')})
        self.assertEqual(result['path'], str(self.root))

    def test_new_folder_requires_session_and_cannot_overwrite_or_escape_parent(self):
        with self.assertRaises(HTTPError) as error:
            self.post('mkdir', {'path': str(self.root), 'name': 'unauthorized'}, token='')
        self.assertEqual(error.exception.code, 403)
        self.assertFalse((self.root / 'unauthorized').exists())
        data = {'path': str(self.root), 'name': 'My analysis'}
        result = self.post('mkdir', data)
        self.assertEqual(result['path'], str(self.root / 'My analysis'))
        sentinel = self.root / 'My analysis' / 'keep.txt'
        sentinel.write_text('preserve')
        for name in ['My analysis', '', '.', '..', '../escape', 'nested/child', 'nested\\child']:
            with self.subTest(name=name), self.assertRaises(HTTPError) as error:
                self.post('mkdir', dict(data, name=name))
            self.assertEqual(error.exception.code, 400)
        self.assertEqual(sentinel.read_text(), 'preserve')

    def test_matrix_discovery_does_not_require_or_create_a_project(self):
        config = {'mode': 'matrix', 'matrix': str(ROOT / 'examples/quantification_matrix.csv')}
        result = self.post('discover', {'config': config})
        self.assertEqual(len(result['board']['rows']), 6)
        self.assertIsNone(self.app.project)
        with self.assertRaisesRegex(ValueError, 'project folder'):
            study.normalize(config)
        self.assertEqual(list(self.root.iterdir()), [])

    def test_raw_table_resolves_relative_paths_without_a_project(self):
        table = self.root / 'samples.csv'
        table.write_text('sample_id,condition,raw_path\nA,Control,A.d\nB,Treatment,B.d\n')
        result = self.post('discover', {'config': {'mode': 'dia', 'samples': str(table)}})
        self.assertEqual([r['raw_path'] for r in result['board']['rows']],
                         [str(self.root / 'A.d'), str(self.root / 'B.d')])


if __name__ == '__main__':
    unittest.main()
