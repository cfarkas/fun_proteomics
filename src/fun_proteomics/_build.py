"""Bundle only reviewed workflow source, never local tools or study data."""
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
from setuptools.command.build_py import build_py
from setuptools.command.sdist import sdist

ROOT = Path(__file__).resolve().parents[2]

def reviewed_files():
    spec = importlib.util.spec_from_file_location('fun_source_validation', ROOT/'scripts/validate_repo.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.main()
    return list(module.source_files(ROOT))


class BuildPy(build_py):
    def run(self):
        super().run()
        target = Path(self.build_lib)/'fun_proteomics/_workflow'
        if target.exists():
            shutil.rmtree(target)
        manifest = {}
        for path in reviewed_files():
            relative = path.relative_to(ROOT)
            if relative.parts[0] in {'src', 'tests', '.github'} or path.name in {'SOURCE_SHA256SUMS', 'fun_build.py', 'pyproject.toml', 'MANIFEST.in'}:
                continue
            destination = target/relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, destination)
            manifest[relative.as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
        (target/'WORKFLOW_SHA256.json').write_text(json.dumps(manifest, sort_keys=True, indent=2)+'\n')


class SourceDistribution(sdist):
    def make_release_tree(self, base_dir, files):
        allowed = {path.relative_to(ROOT).as_posix() for path in reviewed_files()}
        selected = [name for name in files if name in allowed or name == 'PKG-INFO' or name.startswith('src/fun_proteomics.egg-info/')]
        super().make_release_tree(base_dir, selected)
