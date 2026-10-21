#!/usr/bin/env python3
"""Verify wheel/sdist membership and bytes against the reviewed source files."""
import argparse
from email.parser import BytesParser
import hashlib
import json
from pathlib import Path
import tarfile
import zipfile
from validate_repo import ROOT, main as validate, source_files


def sha(data):
    return hashlib.sha256(data).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory',nargs='?',default=str(ROOT/'dist'))
    args=parser.parse_args()
    validate()
    reviewed={p.relative_to(ROOT).as_posix():p.read_bytes() for p in source_files(ROOT)}
    expected_manifest=''.join(sha(data)+'  '+name+'\n' for name,data in sorted(reviewed.items()) if name!='SOURCE_SHA256SUMS')
    assert reviewed['SOURCE_SHA256SUMS'].decode()==expected_manifest, 'Regenerate source checksums with scripts/package_source.py'
    wheels=list(Path(args.directory).glob('*.whl'));sdists=list(Path(args.directory).glob('*.tar.gz'))
    assert len(wheels)==len(sdists)==1, 'Validate exactly one wheel and one source distribution'
    with tarfile.open(sdists[0]) as archive:
        entries={}
        for item in archive.getmembers():
            assert not item.issym() and not item.islnk(), item.name
            name=item.name.partition('/')[2]
            assert not Path(name).is_absolute() and '..' not in Path(name).parts,item.name
            if item.isfile():
                assert name not in entries,name
                entries[name]=archive.extractfile(item).read()
        for name,data in reviewed.items():
            assert entries.get(name)==data, 'Missing or changed source: '+name
        assert all(name in reviewed or name in {'PKG-INFO','setup.cfg'} or name.startswith('src/fun_proteomics.egg-info/') for name in entries), 'Unexpected source distribution file'
        assert entries.get('setup.cfg') == b'[egg_info]\ntag_build = \ntag_date = 0\n\n', 'Unexpected generated setuptools configuration'
    with zipfile.ZipFile(wheels[0]) as archive:
        names=archive.namelist();assert len(names)==len(set(names)), 'Duplicate wheel members'
        prefix='fun_proteomics/_workflow/'
        workflow=json.loads(archive.read(prefix+'WORKFLOW_SHA256.json'))
        expected={name:data for name,data in reviewed.items() if Path(name).parts[0] not in {'src','tests','.github'} and name not in {'SOURCE_SHA256SUMS','pyproject.toml','MANIFEST.in'}}
        assert set(workflow)==set(expected), 'Workflow inventory differs from reviewed source'
        for name,data in expected.items():
            assert archive.read(prefix+name)==data and workflow[name]==sha(data),name
        allowed={prefix+name for name in expected}|{prefix+'WORKFLOW_SHA256.json'}
        for name,data in reviewed.items():
            if name.startswith('src/'):
                member=name.removeprefix('src/')
                assert archive.read(member)==data,name
                allowed.add(member)
        metadata_files=[name for name in names if '.dist-info/' in name]
        assert set(names)==allowed|set(metadata_files), 'Unexpected wheel payload'
        metadata=BytesParser().parsebytes(archive.read(next(name for name in metadata_files if name.endswith('/METADATA'))))
        assert metadata['Name']=='fun_proteomics' and metadata['License-Expression']=='MIT'
        assert metadata['Author']=='cfarkas'
        assert any(name.endswith('/licenses/LICENSE') for name in names)
    print(f'PASS: wheel and sdist contain only reviewed source; {len(expected)} workflow files; MIT metadata; exact source checksums')


if __name__=='__main__':main()
