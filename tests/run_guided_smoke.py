#!/usr/bin/env python3
"""Exercise installed setup, result links, continuation and failed integrity checks."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[1]


def run(*args,success=True):
    command=[sys.executable,'-m','fun_proteomics',*map(str,args)]
    result=subprocess.run(command,env=os.environ.copy(),text=True,capture_output=True)
    if (result.returncode==0)!=success:raise AssertionError(result.stdout+'\n'+result.stderr)
    return result.stdout+result.stderr


def main():
    with tempfile.TemporaryDirectory(prefix='proteomics-guided-installed-') as temp:
        base=Path(temp);project=base/'study with spaces';inputs=base/'inputs';inputs.mkdir()
        for name in ('quantification_matrix.csv','quantification_samples.csv','service_quantification.csv'):
            shutil.copy2(ROOT/'examples'/name,inputs/name)
        run('setup','--non-interactive','--project',project,'--mode','matrix',
            '--matrix',inputs/'quantification_matrix.csv','--samples',inputs/'quantification_samples.csv',
            '--service-matrix',inputs/'service_quantification.csv','--independent-samples','--run')
        config=project/'config/run.json';final=project/'results/final_results'
        assert len(list((final/'differential_expression').glob('*/*/all_results.csv')))==6
        assert (project/'MATRICES.html').is_file() and (project/'START_HERE.html').is_file()
        baseline=(final/'OUTPUT_SHA256SUMS').read_bytes()
        assert 'No analysis needs repeating' in run('--continue',project)
        assert (final/'OUTPUT_SHA256SUMS').read_bytes()==baseline
        target=final/'matrices/protein_intensities.csv';target.write_text('damaged delivery\n')
        run('run','--config',config)
        info=json.loads(run('status','--project',project,'--json'))
        recovered=Path(info['results']).parent
        assert 'continuations' in str(recovered)
        assert target.read_text()=='damaged delivery\n'
        assert (recovered/'matrices/protein_intensities.tsv').is_file()
        assert recovered.as_uri() in (project/'MATRICES.html').read_text()
        assert 'No analysis needs repeating' in run('--continue',project)
        samples=inputs/'quantification_samples.csv';samples.write_text(samples.read_text().replace('Control','Altered'))
        assert 'Recorded input changed' in run('--continue',project,success=False)
        cfg=json.loads(config.read_text());cfg['normalization']='median';config.write_text(json.dumps(cfg))
        assert 'settings changed' in run('--continue',project,success=False)
        print('PASS guided setup: six protein/gene contrasts, CSV/TSV, no-op reuse, damaged-output recovery, updated shortcuts and changed-input/settings refusal')


if __name__=='__main__':main()
