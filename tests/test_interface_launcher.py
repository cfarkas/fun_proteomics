"""The interactive launcher must never claim the analysis workflow completed."""
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = ROOT / '.runtime/env/bin/nextflow'


@unittest.skipUnless(NEXTFLOW.is_file(), 'Project Nextflow runtime unavailable')
class InterfaceLauncherTests(unittest.TestCase):
    def test_interface_exit_is_distinct_from_workflow_completion(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            scripts = root / 'scripts'
            scripts.mkdir()
            shutil.copy2(ROOT / 'main.nf', root / 'main.nf')
            shutil.copytree(ROOT / 'modules/diapasef', root / 'modules/diapasef')
            (root / 'config.local.json').write_text('{}')
            runtime = root / 'runtime'
            (runtime / 'env/bin').mkdir(parents=True)
            (runtime / 'env/bin/python').symlink_to(sys.executable)
            (root / 'nextflow.config').write_text(f"params.runtime_dir = '{runtime}'\n")
            # Exercise real Nextflow control flow without installing dependencies,
            # opening a browser, or launching any scientific tasks.
            (scripts / 'bootstrap_runtime.sh').write_text('printf \'{"rscript":"/bin/true"}\\n\' > "$3"\n')
            (scripts / 'system_capabilities.py').write_text('print(\'{"available_cpus":1,"max_memory_gib":8}\')\n')
            env = os.environ.copy()
            env.update(JAVA_HOME=str(ROOT / '.runtime/env'), NXF_HOME=str(ROOT / '.runtime/nextflow'),
                       NXF_VER='26.04.2', NXF_DISABLE_CHECK_LATEST='true')
            for code, label in [(0, 'CLOSED'), (3, 'FAILED')]:
                with self.subTest(code=code):
                    (scripts / 'config_server.py').write_text(f'raise SystemExit({code})\n')
                    result = subprocess.run([str(NEXTFLOW), 'run', 'main.nf', '--interactively', '-ansi-log', 'false'],
                                            cwd=root, env=env, capture_output=True, text=True, timeout=60)
                    output = result.stdout + result.stderr
                    self.assertEqual(result.returncode == 0, code == 0, output)
                    self.assertIn(f'[INTERFACE {label}]', output)
                    self.assertNotIn('[WORKFLOW COMPLETE]', output)
                    self.assertNotIn('Submitted process', output)
