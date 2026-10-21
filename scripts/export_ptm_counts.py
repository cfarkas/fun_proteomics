#!/usr/bin/env python3
"""Export gene-level modified peptide counts from a completed cohort/mode search."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'pipeline_lib'))
from ptm_counts import main
if __name__ == '__main__': main()
