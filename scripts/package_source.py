#!/usr/bin/env python3
"""Package reviewed source files without following private input/tool links."""
import hashlib
import zipfile
from pathlib import Path
from validate_repo import ROOT, main as validate, source_files

validate()
files = sorted(p for p in source_files(ROOT) if p.name != 'SOURCE_SHA256SUMS')
manifest = ROOT / 'SOURCE_SHA256SUMS'
manifest.write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest() + '  ' +
                            str(p.relative_to(ROOT)) + '\n' for p in files))
archive = ROOT.parent / (ROOT.name + '-source.zip')
with zipfile.ZipFile(archive, 'w', compression=zipfile.ZIP_DEFLATED) as output:
    for path in files + [manifest]:
        output.write(path, ROOT.name + '/' + str(path.relative_to(ROOT)))
print(f'Packaged {len(files) + 1} source files: {archive}')
checksum = hashlib.sha256(archive.read_bytes()).hexdigest()
archive.with_suffix(archive.suffix + '.sha256').write_text(checksum + '  ' + archive.name + '\n')
print('Archive SHA256: ' + checksum)
