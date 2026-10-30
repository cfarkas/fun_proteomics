#!/usr/bin/env python3
"""Record the browser quick start using the public DEP UbiLength dataset.

Requires the installed package, its R runtime, Playwright and ffmpeg. Input
CSV files come from the documented public example archive; no private study
files or raw spectra are used. The clip records an actual matrix analysis.
"""
import argparse
import csv
from contextlib import suppress
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess

from record_walkthroughs import Clip, drag_sample, launch


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--python', required=True)
    parser.add_argument('--dataset', type=Path, required=True)
    parser.add_argument('--work-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--browser-path')
    args = parser.parse_args()
    args.work_dir = args.work_dir.resolve()
    args.output = args.output.resolve()
    args.dataset = args.dataset.resolve()
    provenance = json.loads((args.dataset/'provenance.json').read_text())
    assert provenance['dataset'] == 'DEP UbiLength', 'Use the documented public example'
    with (args.dataset/'samples.csv').open(newline='') as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 12 and {r['condition'] for r in rows} == {'Ctrl', 'Ubi1', 'Ubi4', 'Ubi6'}
    args.work_dir.mkdir(parents=True, exist_ok=False)
    args.output.mkdir(parents=True, exist_ok=False)
    inputs = args.work_dir/'ubiquitin-lfq-example'
    inputs.mkdir()
    for name in ('proteins.csv', 'samples.csv'):
        shutil.copy2(args.dataset/name, inputs/name)
    env = os.environ.copy()
    env.pop('PYTHONPATH', None)
    env.update(FUN_PROTEOMICS_CACHE=str(args.work_dir/'workflow-cache'), PYTHONUNBUFFERED='1')
    version = subprocess.check_output([args.python, '-m', 'fun_proteomics', '--version'], env=env, text=True).strip()
    process, url = launch(args, args.work_dir, env)
    context = None
    from playwright.sync_api import sync_playwright
    try:
        with sync_playwright() as p:
            browser = p.chromium.launch(headless=True, **({'executable_path': args.browser_path} if args.browser_path else {}))
            context = browser.new_context(viewport={'width': 1400, 'height': 1040}, record_video_dir=str(args.output/'capture'), record_video_size={'width': 1400, 'height': 1040})
            page = context.new_page()
            errors = []
            page.on('pageerror', lambda error: errors.append(str(error)))
            page.goto(url)
            page.wait_for_function('document.getElementById("project").value.length>0')
            clip = Clip(page, 'ubiquitin-quickstart', 'Interactive quick start · real published LFQ data', args.output)
            page.locator('[data-step="2"]').click();page.locator('#project').fill(str(args.work_dir/'ubiquitin-study'));page.locator('[data-step="1"]').click()
            clip.say('Choose a new output folder. This walkthrough uses the downloaded 12-sample ubiquitin-affinity dataset.', 6)
            for field, filename, caption in [
                ('matrix', 'proteins.csv', 'Browse to the extracted example folder and select proteins.csv. This is the measured protein quantity table.'),
                ('samples', 'samples.csv', 'Select samples.csv to prefill the four conditions and their experimental replicates.')
            ]:
                page.locator('[data-browse='+field+']').click()
                page.locator('#browse-path').fill(str(inputs))
                page.locator('#browse-go').click()
                page.locator('#entries').get_by_role('button', name=filename, exact=True).wait_for()
                clip.say(caption, 5)
                page.locator('#entries').get_by_role('button', name=filename, exact=True).click()
            page.locator('#discover').click()
            page.wait_for_function('document.querySelectorAll(".sample").length===12')
            clip.say('Discover samples. Review three cards in each condition: Ctrl, Ubi1, Ubi4 and Ubi6.', 6)
            drag_sample(page, 'Ubi4_1', '')
            clip.say('Drag a card to Unassigned to exclude it. Here we will restore it so all 12 samples remain included.', 5)
            drag_sample(page, 'Ubi4_1', 'Ubi4')
            clip.say('Drag the card back into its condition. Each card also has a condition dropdown.', 5)
            page.locator('#next').click()
            page.locator('#independent_samples').check()
            clip.say('Confirm independent experimental replicates for this tutorial. The source design supplies no donor pairing.', 5)
            advanced = page.locator('details').filter(has=page.locator('#normalization'))
            advanced.locator('summary').click()
            page.locator('#normalization').select_option('none')
            page.locator('#normalization').scroll_into_view_if_needed()
            clip.say('Keep all condition pairs and no additional normalization for these LFQ quantities. Missing values remain missing.', 6)
            page.locator('#check').click()
            page.wait_for_function('document.getElementById("preview").textContent.includes("6 condition pairs")', timeout=60000)
            preview = page.locator('#preview').inner_text()
            assert '2902' in preview and '2696' in preview, preview
            page.locator('#preview').scroll_into_view_if_needed()
            clip.say('Check inputs: 12 samples, six condition pairs, 2,902 protein-group rows and 2,696 unambiguous gene rows.', 7)
            page.locator('#save').click()
            page.wait_for_function('!document.getElementById("run").disabled', timeout=60000)
            clip.say('Save configuration and check, then click Run / continue analysis. The settings and sample assignments stay in your output folder.', 6)
            page.locator('#run').click()
            clip.say('The real matrix analysis is running. It creates six protein and six gene comparisons; no raw spectrum search is needed.', 5)
            page.wait_for_function('document.getElementById("state").textContent==="completed"', timeout=240000)
            clip.say('Analysis complete. Open the result dashboard or go directly to the CSV and TSV matrices.', 6)
            base_url = url.split('/#')[0]
            page.goto(base_url+'/results/MATRICES.html')
            clip.say('The matrix catalogue provides protein and gene quantities in both CSV and TSV. The result dashboard links every comparison and QC plot.', 7)
            assert 'CSV' in page.locator('body').inner_text() and 'TSV' in page.locator('body').inner_text()
            page.goto(url)
            page.wait_for_function('!document.getElementById("run").disabled')
            page.locator('#run').click()
            page.wait_for_function('document.getElementById("log").textContent.includes("No analysis needs repeating")', timeout=60000)
            clip.say('Run / continue verifies and reuses the completed results. Keep this project folder to return to the same analysis.', 6)
            final = args.work_dir/'ubiquitin-study/results/final_results'
            assert len(list((final/'differential_expression').glob('*/*/all_results.csv'))) == 12
            assert not errors, errors
            recorded = clip.finish(context)
            context = None
            browser.close()
        metadata = {
            'version': version,
            'input_kind': 'real public DEP UbiLength quantification; no raw spectra',
            'dataset_doi': '10.1016/j.molcel.2017.01.004',
            'dataset_license': 'Artistic-2.0',
            'input_sha256': {name: hashlib.sha256((inputs/name).read_bytes()).hexdigest() for name in ('proteins.csv', 'samples.csv')},
            'recorder_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'clips': [recorded],
            'checked_outputs': {'samples': 12, 'conditions': 4, 'protein_rows': 2902, 'gene_rows': 2696, 'comparisons': 12, 'completed_continuation': True}
        }
        (args.output/'recording-public-quickstart.json').write_text(json.dumps(metadata, indent=2)+'\n')
        print('PASS: recorded the real public-data browser analysis and verified completed continuation', flush=True)
    finally:
        if context:
            with suppress(Exception):
                context.close()
        if process.poll() is None:
            process.send_signal(signal.SIGINT)
            try:
                process.wait(timeout=15)
            except subprocess.TimeoutExpired:
                process.terminate()
                process.wait(timeout=5)
        process.stdout.close()


if __name__ == '__main__':
    main()
