#!/usr/bin/env python3
"""Create tiny synthetic native-input fixtures; never real spectra or tool binaries."""
import hashlib
import json
from pathlib import Path
import sqlite3
import sys


def make_fixture(out):
    out = Path(out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    samples = []
    for sid, cohort in [('Control_A', 'cell'), ('Drug_B', 'cell'), ('Serum_A', 'serum')]:
        raw = out / 'raw with spaces' / (sid + '.d')
        raw.mkdir(parents=True, exist_ok=True)
        with sqlite3.connect(raw / 'analysis.tdf') as db:
            db.executescript('DROP TABLE IF EXISTS Frames; DROP TABLE IF EXISTS Precursors; CREATE TABLE Frames (MsMsType INTEGER); INSERT INTO Frames VALUES (0), (8); CREATE TABLE Precursors (Id INTEGER); INSERT INTO Precursors VALUES (1);')
        (raw / 'analysis.tdf_bin').write_bytes(b'SYNTHETIC STUB ONLY')
        samples.append({'sample_id': sid, 'cohort': cohort, 'condition': 'Control' if sid != 'Drug_B' else 'Treatment', 'replicate': 1, 'raw_path': str(raw)})
    fasta = out / 'test.fasta'
    fasta.write_text('>sp|TEST1|TEST_HUMAN OS=Homo sapiens\nMPEPTIDER\n>contam_sp|TEST2|CONTAM\nMPEPTIDEK\n>rev_sp|TEST1|TEST_HUMAN OS=Homo sapiens\nMEDITPEPR\n>rev_contam_sp|TEST2|CONTAM\nMEDITPEPK\n')
    tools = {}
    for key in ['msfragger_zip', 'ionquant_zip', 'diatracer_zip']:
        archive = out / (key + '.placeholder')
        archive.write_text('STUB — NOT A TOOL ARCHIVE\n')
        tools[key] = str(archive)
    cfg = {'outdir': str(out / 'results'), 'modes': ['DDA+', 'DDA'], 'resources': {'cpus': 1, 'memory_gib': 4},
           'samples': samples, 'database': {'fasta': str(fasta), 'sha256': hashlib.sha256(fasta.read_bytes()).hexdigest(), 'contains_decoys': True}, 'tools': tools}
    config = out / 'test.local.json'
    config.write_text(json.dumps(cfg, indent=2)+'\n')
    return config


if __name__ == '__main__':
    print(make_fixture(sys.argv[1]))
