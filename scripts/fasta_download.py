"""Source-agnostic FASTA acquisition and database preparation (metadata stays internal)."""
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import time
import urllib.parse
import urllib.request
import zipfile

CRAP_URL = 'https://raw.githubusercontent.com/Nesvilab/philosopher/master/lib/dat/cRAPS.zip'
CRAP_SHA256 = '629a91aa5d72730e6a0ffcb9a092380ae72ab39ab3635702835f9217573018af'


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def checked(path, expected=None):
    path = Path(path)
    if not path.is_file():
        raise ValueError(f'Missing database input: {path}')
    observed = sha256(path)
    if expected and observed != expected:
        raise ValueError(f'FASTA/contaminant SHA256 mismatch: expected {expected}, got {observed}')
    return observed


def normalize(db, base):
    db = dict(db)
    db.setdefault('source', 'local' if db.get('fasta') else 'url' if db.get('url') else 'uniprot')
    db.setdefault('contains_decoys', 'auto')
    db.setdefault('contaminants', 'custom' if db.get('contaminants_zip') or db.get('contaminants_fasta') else 'auto')
    db.setdefault('reviewed', True)
    db.setdefault('include_isoforms', False)
    db.setdefault('offline', False)
    db.setdefault('refresh', False)
    db.setdefault('timeout', 120)
    db.setdefault('retries', 3)
    db.setdefault('cache_dir', '.database_cache')
    for key in ['fasta', 'contaminants_zip', 'contaminants_fasta', 'cache_dir']:
        if db.get(key):
            path = Path(db[key]).expanduser()
            db[key] = str((path if path.is_absolute() else Path(base) / path).resolve())
    if db['source'] == 'local':
        if not db.get('fasta') or db.get('url') or db.get('taxon') or db.get('proteome'):
            raise ValueError('Local database needs fasta only; do not combine local, URL and UniProt selectors')
        db['sha256'] = checked(db['fasta'], db.get('sha256'))
    elif db['source'] == 'url':
        if not db.get('url') or db.get('fasta') or db.get('taxon') or db.get('proteome'):
            raise ValueError('URL database needs url only; do not combine database selectors')
        require_http(db['url'])
    elif db['source'] == 'uniprot':
        if db.get('fasta') or db.get('url') or bool(db.get('taxon')) == bool(db.get('proteome')):
            raise ValueError('UniProt download requires exactly one of --taxon or --proteome, not a local FASTA/URL')
    else:
        raise ValueError('database.source must be local, url or uniprot')
    if db.get('contaminants_zip') and db.get('contaminants_fasta'):
        raise ValueError('Choose contaminants_fasta OR contaminants_zip')
    if db['contaminants'] == 'custom':
        path = db.get('contaminants_fasta') or db.get('contaminants_zip')
        if not path:
            raise ValueError('Custom contaminants need --contaminants_fasta or --contaminants_zip')
        db['contaminants_sha256'] = checked(path, db.get('contaminants_sha256'))
    elif db.get('contaminants_zip') or db.get('contaminants_fasta'):
        raise ValueError('Custom contaminant files cannot be combined with contaminants auto/none')
    if db['offline'] and db['refresh']:
        raise ValueError('Cannot refresh database downloads in offline mode')
    return db


def require_http(url):
    parsed = urllib.parse.urlparse(url)
    if parsed.scheme not in ('http', 'https') or not parsed.netloc or parsed.username or parsed.password:
        raise ValueError('Database URL must be an HTTP(S) URL without embedded credentials; use --fasta for local files')


def uniprot_url(db):
    selector = f'organism_id:{db["taxon"]}' if db.get('taxon') else f'proteome:{db["proteome"]}'
    query = selector + (' AND reviewed:true' if db.get('reviewed', True) else '')
    return 'https://rest.uniprot.org/uniprotkb/stream?' + urllib.parse.urlencode({
        'query': query, 'format': 'fasta', 'includeIsoform': str(db.get('include_isoforms', False)).lower(), 'compressed': 'false'})


def download(url, db, expected=None, validator=None):
    require_http(url)
    cache = Path(db['cache_dir'])
    cache.mkdir(parents=True, exist_ok=True)
    key = hashlib.sha256(url.encode()).hexdigest()
    data, receipt = cache / (key + '.data'), cache / (key + '.json')
    with (cache / (key + '.lock')).open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if data.exists() and receipt.exists() and not db.get('refresh'):
            info = json.loads(receipt.read_text())
            if info.get('url') != url:
                raise ValueError('Database cache URL mismatch')
            checked(data, info['sha256'])
            checked(data, expected)
            if validator:
                validator(data)
            print(f'[DATABASE] Reusing verified download {data}', flush=True)
            return data, info
        if db.get('offline'):
            raise ValueError('Database is not in the verified cache and --database_offline is enabled')
        partial = cache / (key + f'.{os.getpid()}.partial')
        print(f'[DATABASE] Downloading {url}', flush=True)
        try:
            for attempt in range(db.get('retries', 3) + 1):
                try:
                    req = urllib.request.Request(url, headers={'User-Agent': 'dda-pasef-pipeline/0.1', 'Accept-Encoding': 'identity'})
                    with urllib.request.urlopen(req, timeout=db.get('timeout', 120)) as response, partial.open('wb') as output:
                        shutil.copyfileobj(response, output)
                        headers = dict(response.headers.items())
                    break
                except (OSError, TimeoutError):
                    if attempt == db.get('retries', 3):
                        raise
                    time.sleep(min(2 ** attempt, 8))
            if partial.stat().st_size == 0:
                raise ValueError('Downloaded database is empty; check the organism/proteome/query')
            observed = checked(partial, expected)
            if validator:
                validator(partial)
            info = {'url': url, 'sha256': observed, 'headers': headers, 'downloaded_at_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}
            temp_receipt = cache / (key + f'.{os.getpid()}.receipt.partial')
            temp_receipt.write_text(json.dumps(info, indent=2) + '\n')
            partial.replace(data)
            temp_receipt.replace(receipt)
            return data, info
        finally:
            partial.unlink(missing_ok=True)


def records(text):
    header, sequence = None, []
    for line in text.splitlines():
        if line.startswith('>'):
            if header is not None:
                yield header, ''.join(sequence)
            header, sequence = line[1:].strip(), []
        elif line.strip():
            if header is None:
                raise ValueError('Not a FASTA: sequence/text before first header')
            sequence.append(line.strip())
    if header is not None:
        yield header, ''.join(sequence)


def read_records(path, zipped=False):
    if zipped:
        with zipfile.ZipFile(path) as archive:
            text = archive.read('crap-gpmdb.fas').decode('utf-8')
    else:
        data = Path(path).read_bytes()
        text = (gzip.decompress(data) if data.startswith(b'\x1f\x8b') else data).decode('utf-8')
    result = list(records(text))
    if not result or any(not h or not seq or not re.fullmatch('[A-Z]+', seq) for h, seq in result):
        raise ValueError('FASTA must contain nonempty headers and uppercase amino-acid sequences')
    names = [h.split()[0] for h, _ in result]
    if len(names) != len(set(names)):
        raise ValueError('Duplicate FASTA identifiers')
    return result


def accession(header):
    fields = header.split()[0].split('|')
    return fields[1] if len(fields) > 1 else fields[0]


def decoy(sequence):
    segments = re.findall(r'.*?[KR]|.+$', sequence)
    return ''.join(s[0] + s[1:-1][::-1] + s[-1] if len(s) > 1 else s for s in segments)


def build(db):
    source_info = {}
    if db['source'] == 'local':
        source = Path(db['fasta'])
        checked(source, db.get('sha256'))
    else:
        url = uniprot_url(db) if db['source'] == 'uniprot' else db['url']
        source, source_info = download(url, db, db.get('sha256'), read_records)
    proteins = read_records(source)
    has_decoys = any(h.startswith('rev_') for h, _ in proteins)
    expected_decoys = db.get('contains_decoys', 'auto')
    if expected_decoys != 'auto' and expected_decoys != has_decoys:
        raise ValueError('FASTA contents disagree with --fasta_contains_decoys (required prefix: rev_)')
    targets = [(h, s) for h, s in proteins if not h.startswith('rev_')]
    if has_decoys:
        if {h.split()[0] for h, _ in targets} != {h.split()[0][4:] for h, _ in proteins if h.startswith('rev_')}:
            raise ValueError('Prebuilt FASTA must have exactly one rev_ decoy per target identifier')
    if not targets:
        raise ValueError('FASTA contains no target proteins')
    additions, contaminant_info = [], {}
    mode = db.get('contaminants', 'auto')
    if mode == 'custom' or (mode == 'auto' and not any(h.startswith('contam_') for h, _ in targets)):
        if mode == 'auto':
            contamination, contaminant_info = download(CRAP_URL, db, CRAP_SHA256, lambda path: read_records(path, zipped=True))
            contaminants = read_records(contamination, zipped=True)
        else:
            contamination = db.get('contaminants_fasta') or db['contaminants_zip']
            checked(contamination, db.get('contaminants_sha256'))
            contaminants = read_records(contamination, zipped=bool(db.get('contaminants_zip')))
        seen = {accession(h) for h, _ in targets}
        for header, seq in contaminants:
            if header.startswith('rev_'):
                raise ValueError('Contaminant inputs must not contain decoys')
            if accession(header) not in seen:
                additions.append((header if header.startswith('contam_') else 'contam_' + header, seq))
                seen.add(accession(header))
    if has_decoys:
        proteins += additions + [('rev_' + h, decoy(s)) for h, s in additions]
    else:
        targets += additions
        proteins = targets + [('rev_' + h, decoy(s)) for h, s in targets]
    names = [h.split()[0] for h, _ in proteins]
    if len(names) != len(set(names)):
        raise ValueError('Duplicate identifiers after contaminant/decoy construction')
    if has_decoys and not additions and not source.read_bytes().startswith(b'\x1f\x8b'):
        shutil.copyfile(source, 'search.fasta')
    else:
        Path('search.fasta').write_text(''.join(f'>{h}\n{s}\n' for h, s in proteins))
    info = {**db, 'source_sha256': sha256(source), 'download': source_info,
            'contaminant_download': contaminant_info, 'search_fasta_sha256': sha256('search.fasta'),
            'entries': len(proteins), 'decoy_method': 'Existing paired rev_ decoys retained; missing decoys use KR segment internal reversal'}
    Path('database_provenance.json').write_text(json.dumps(info, indent=2) + '\n')
    print(f'[DATABASE] Ready: {len(proteins)} entries; SHA256 {info["search_fasta_sha256"]}', flush=True)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Prepare a UniProt FASTA for the study editor')
    parser.add_argument('--taxon', type=int, required=True)
    parser.add_argument('--cache', required=True)
    parser.add_argument('--unreviewed', action='store_true')
    args = parser.parse_args()
    if args.taxon <= 0:
        parser.error('Taxonomy ID must be positive')
    build(normalize({'source': 'uniprot', 'taxon': args.taxon,
                     'reviewed': not args.unreviewed, 'cache_dir': args.cache}, Path.cwd()))
