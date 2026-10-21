#!/usr/bin/env python3
"""Open the same guided setup from a source checkout or direct Nextflow."""
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from fun_proteomics.cli import main
if __name__ == '__main__':
    raise SystemExit(main(['--workflow-dir',str(ROOT),'setup',*sys.argv[1:]]))
