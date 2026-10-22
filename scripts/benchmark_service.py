#!/usr/bin/env python3
"""Benchmark matched service and pipeline quantities from any completed run."""
import sys
from analyze_quantification import main
if __name__ == '__main__':
    if not any(v in sys.argv for v in ('-h','--help')) and not any(v=='--service-matrix' or v.startswith('--service-matrix=') for v in sys.argv):
        raise SystemExit('benchmark requires --service-matrix, --matrix, --samples and --outdir')
    main()
