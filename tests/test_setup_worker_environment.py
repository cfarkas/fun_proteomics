"""Browser workers must preserve Python's standard-library import precedence."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
if os.environ.get('FUN_PROTEOMICS_TEST_INSTALLED') != '1':
    sys.path.insert(0, str(ROOT/'src'))
from fun_proteomics import setup_web


class WorkerEnvironmentTests(unittest.TestCase):
    def run_python(self, env, code, cwd):
        result = subprocess.run([sys.executable, '-c', code], env=env, cwd=cwd, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout+'\n'+result.stderr)
        return result.stdout.strip()

    def test_installed_worker_does_not_shadow_stdlib_with_argparse_backport(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            installed = root/'site-packages'
            package = installed/'fun_proteomics'
            package.mkdir(parents=True)
            (installed/'argparse.py').write_text("raise RuntimeError('backport shadowed the standard library')\n")
            env = dict(os.environ)
            env.pop('PYTHONPATH', None)
            with patch.object(setup_web, '__file__', str(package/'setup_web.py')), \
                 patch.object(setup_web.site, 'getsitepackages', return_value=[str(installed)]), \
                 patch('fun_proteomics.cli.environment', return_value=env):
                child = setup_web.worker_environment(ROOT, object())
            self.run_python(child, 'import argparse; p=argparse.ArgumentParser(); p.add_subparsers(); p.parse_args([])', root)

    def test_source_worker_retains_source_import_and_existing_pythonpath(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root/'src'
            package = source/'fun_proteomics'
            package.mkdir(parents=True)
            (package/'__init__.py').write_text("ORIGIN='source-checkout'\n")
            extra = root/'additional_modules'
            extra.mkdir()
            (extra/'worker_extra.py').write_text("VALUE='existing-path'\n")
            env = dict(os.environ, PYTHONPATH=str(extra))
            with patch.object(setup_web, '__file__', str(package/'setup_web.py')), \
                 patch('fun_proteomics.cli.environment', return_value=env):
                child = setup_web.worker_environment(ROOT, object())
            output = self.run_python(child, 'import fun_proteomics, worker_extra; print(fun_proteomics.ORIGIN, worker_extra.VALUE)', root)
            self.assertEqual(output, 'source-checkout existing-path')


if __name__ == '__main__':
    unittest.main()
