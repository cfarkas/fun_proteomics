"""Resolve a checkout or materialize the packaged workflow in a writable cache."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024*1024), b''):
            value.update(block)
    return value.hexdigest()


def verify_files(root, manifest):
    root = Path(root).resolve()
    for name, expected in manifest.items():
        path = root/name
        if Path(name).is_absolute() or '..' in Path(name).parts or not path.resolve().is_relative_to(root):
            raise ValueError('Unsafe workflow manifest path: '+name)
        if not path.is_file() or digest(path) != expected:
            raise ValueError('Workflow source checksum mismatch: '+str(path))


def resolve_project(project=None, cache_dir=None):
    if project:
        root = Path(project).expanduser().resolve()
        if not (root/'main.nf').is_file() or not (root/'scripts/reanalyze_dia.py').is_file():
            raise ValueError('Not a fun_proteomics workflow checkout: '+str(root))
        return root
    checkout = Path(__file__).resolve().parents[2]
    if (checkout/'main.nf').is_file() and (checkout/'pyproject.toml').is_file():
        return checkout
    payload = Path(__file__).with_name('_workflow')
    manifest_path = payload/'WORKFLOW_SHA256.json'
    if not manifest_path.is_file():
        raise ValueError('Packaged workflow is missing; reinstall fun_proteomics')
    manifest = json.loads(manifest_path.read_text())
    verify_files(payload, manifest)
    from . import __version__
    cache = Path(cache_dir or os.environ.get('FUN_PROTEOMICS_CACHE') or Path(os.environ.get('XDG_CACHE_HOME', Path.home()/'.cache'))/'fun_proteomics').expanduser().resolve()
    root = cache/(__version__+'-'+digest(manifest_path)[:16])
    if root.is_dir():
        verify_files(root, manifest)
        return root
    cache.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix='.preparing-', dir=cache))
    try:
        shutil.copytree(payload, temporary, dirs_exist_ok=True)
        verify_files(temporary, manifest)
        try:
            temporary.rename(root)
        except FileExistsError:
            verify_files(root, manifest)
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)
    return root
