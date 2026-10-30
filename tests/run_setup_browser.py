#!/usr/bin/env python3
"""Set up and reopen an analysis using folder/file buttons, with no typed paths."""
import argparse
import json
import os
from pathlib import Path
import signal
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from record_walkthroughs import fixture, launch


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--python', default=sys.executable)
    parser.add_argument('--browser-path')
    parser.add_argument('--source', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    from playwright.sync_api import sync_playwright
    with tempfile.TemporaryDirectory(prefix='proteomics-simple-setup-') as temp:
        root = Path(temp)
        fixture(root)
        raw = root / 'raw samples'
        raw.mkdir()
        for name in ['C1', 'C2', 'T1', 'T2', 'BLANK']:
            folder = raw / (name + '.d')
            folder.mkdir()
            (folder / 'analysis.tdf').write_bytes(b'fixture')
            (folder / 'analysis.tdf_bin').write_bytes(b'fixture')
        env = os.environ.copy()
        if args.source:
            env['PYTHONPATH'] = str(ROOT / 'src')
        else:
            env.pop('PYTHONPATH', None)
        env.update(FUN_PROTEOMICS_CACHE=str(root / 'cache'), PYTHONUNBUFFERED='1')
        process, url = launch(args, root, env, project=False)
        try:
            with sync_playwright() as playwright:
                browser = playwright.chromium.launch(headless=True, **({'executable_path': args.browser_path} if args.browser_path else {}))
                page = browser.new_page(viewport={'width': 1440, 'height': 1000})
                errors = []
                page.on('pageerror', lambda error: errors.append(str(error)))
                page.goto(url)
                page.wait_for_function("lastFolder.length > 0")
                assert page.locator('#discover').is_visible()
                assert not page.locator('#project').is_visible()
                assert not page.locator('#samples').is_visible()
                assert page.locator('#discover').bounding_box()['y'] < page.locator('#mode').bounding_box()['y']
                page.locator('#mode').select_option('dia')
                if args.output:
                    args.output.parent.mkdir(parents=True, exist_ok=True)
                    page.screenshot(path=str(args.output), full_page=True)
                # The first Discover click opens a chooser when no data is selected.
                page.locator('#discover').click()
                page.get_by_role('button', name='Open folder raw samples', exact=True).click()
                page.locator('#choose-folder').click()
                page.wait_for_function("document.querySelectorAll('.sample').length === 4")
                assert page.locator('#project').input_value() == ''
                assert not list(raw.glob('config'))
                # Change data type and select a matrix by navigating the file picker.
                page.locator('#mode').select_option('matrix')
                page.locator('[data-browse=matrix]').click()
                page.get_by_role('button', name='Working folder', exact=True).click()
                page.get_by_role('button', name='Open folder inputs', exact=True).click()
                page.locator('#browse-search').fill('PROTEINS')
                page.get_by_role('button', name='Choose file proteins.csv', exact=True).click()
                page.wait_for_function("document.querySelectorAll('.sample').length === 6")
                for sample in ['C1', 'C2', 'A1', 'A2', 'B1', 'B2']:
                    page.get_by_label('Condition for ' + sample, exact=True).select_option('Control' if sample.startswith('C') else 'Treatment')
                page.locator('#next').click()
                assert not page.locator('#service_matrix').is_visible()
                assert not page.locator('#subject_column').is_visible()
                page.locator('[data-browse=project]').click()
                page.get_by_role('button', name='Working folder', exact=True).click()
                page.locator('#new-folder').click()
                page.locator('#new-folder-name').fill('My analysis')
                page.locator('#create-folder').click()
                page.wait_for_function("!document.getElementById('picker').open")
                assert page.locator('#project').input_value() == str(root / 'My analysis')
                page.locator('#independent_samples').check()
                page.locator('#save').click()
                page.wait_for_function("!document.getElementById('run').disabled")
                config = json.loads((root / 'My analysis/config/run.json').read_text())
                assert len(config['assignment_board']['rows']) == 6
                # Opening a saved project also uses the folder picker.
                page.locator('[data-step="1"]').click()
                page.locator('#load').click()
                page.locator('#choose-folder').click()
                page.wait_for_function("document.querySelector('[data-panel=\"3\"]').hidden === false")
                assert page.locator('#project-name').inner_text() == str(root / 'My analysis')
                page.locator('[data-step="1"]').click()
                page.set_viewport_size({'width': 390, 'height': 844})
                assert page.evaluate('document.documentElement.scrollWidth <= window.innerWidth')
                assert not errors, errors
                browser.close()
                print('PASS: discovery first, raw/blank scan, file search, automatic discovery, create folder, save, reopen and mobile layout without typing a path')
        finally:
            process.send_signal(signal.SIGINT)
            process.wait(timeout=15)
            process.stdout.close()


if __name__ == '__main__':
    main()
