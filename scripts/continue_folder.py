#!/usr/bin/env python3
"""Bridge the direct Nextflow folder-continuation interface to native Python."""
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from fun_proteomics.cli import main
if __name__ == '__main__':
    raise SystemExit(main(['--project',str(ROOT),*sys.argv[1:]]))
