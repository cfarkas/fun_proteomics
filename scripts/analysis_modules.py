"""Install the supported phosphoproteomics adapter from the pinned FragPipe archive."""
import hashlib
import json
from pathlib import Path
import zipfile
import os
import secrets
ROOT=Path(__file__).resolve().parents[1]
TARGET=ROOT/'modules/phosphoproteomics/fragpipe.workflow'
def digest(data): return hashlib.sha256(data).hexdigest()
def status():
    marker=TARGET.with_suffix('.json')
    try:
        metadata=json.loads(marker.read_text())
        if digest(TARGET.read_bytes())!=metadata['sha256']: return {'installed':False}
        return {'installed':True,**metadata}
    except (OSError,ValueError,KeyError): return {'installed':False}
def install(tools, force=False):
    from tool_manager import ensure_archive, LOCK
    with LOCK:
        archive=ensure_archive('FragPipe',tools,force=force)
    spec=json.loads((ROOT/'assets/tool_versions.json').read_text())['FragPipe']
    with archive.open('rb') as f:
        h=hashlib.file_digest(f,'sha256').hexdigest()
    if h!=spec['archive_sha256']: raise ValueError('FragPipe archive checksum does not match the supported version.')
    with zipfile.ZipFile(archive) as z:
        text=z.read('fragpipe-24.0/workflows/LFQ-phospho.workflow').decode()
    properties=dict(line.split('=',1) for line in text.splitlines() if line and not line.startswith('#') and '=' in line)
    if properties.get('ptmprophet.run-ptmprophet')!='true' or 'STY' not in properties.get('ptmprophet.cmdline',''): raise ValueError('Archive lacks the expected phosphorylation localization workflow.')
    properties.update({'workflow.input.data-type.im-ms':'true','workflow.input.data-type.regular-ms':'false'})
    payload=('\n'.join(k+'='+v for k,v in properties.items())+'\n').encode()
    TARGET.parent.mkdir(parents=True,exist_ok=True)
    temporary=TARGET.with_name('.'+secrets.token_hex(8)+'.tmp');temporary.write_bytes(payload);os.replace(temporary,TARGET)
    metadata={'sha256':digest(payload),'name':'DDA-PASEF LFQ-phospho','fragpipe_version':spec['version']}
    temporary=TARGET.with_name('.'+secrets.token_hex(8)+'.json');temporary.write_text(json.dumps(metadata));os.replace(temporary,TARGET.with_suffix('.json'))
    return status()
