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
    parser.add_argument('--native', action='store_true', help='Test the system-dialog integration with a scripted Zenity executable')
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
        # Exercise browser navigation without opening a real desktop dialog.
        env.pop('DISPLAY', None)
        env.pop('WAYLAND_DISPLAY', None)
        if args.source:
            env['PYTHONPATH'] = str(ROOT / 'src')
        else:
            env.pop('PYTHONPATH', None)
        env.update(FUN_PROTEOMICS_CACHE=str(root / 'cache'), PYTHONUNBUFFERED='1')
        if args.native:
            # Keep the real browser/API/subprocess path; script only the OS dialog.
            dialog_bin = root / 'bin'
            dialog_bin.mkdir()
            response = root / 'dialog-response.json'
            calls = root / 'dialog-calls.jsonl'
            executable = dialog_bin / 'zenity'
            executable.write_text('#!' + sys.executable + '\n' + r"""import json, os, sys
from pathlib import Path
with Path(os.environ['TEST_DIALOG_CALLS']).open('a') as handle:
    handle.write(json.dumps(sys.argv[1:]) + '\n')
result = json.loads(Path(os.environ['TEST_DIALOG_RESPONSE']).read_text())
print(result.get('path', ''))
raise SystemExit(result.get('code', 0))
""")
            executable.chmod(0o755)
            env.update(DISPLAY=':test', PATH=str(dialog_bin) + os.pathsep + env['PATH'],
                       TEST_DIALOG_CALLS=str(calls), TEST_DIALOG_RESPONSE=str(response))
            response.write_text(json.dumps({'path': str(raw)}))
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
                page.locator('#input-options > summary').click()
                assert page.get_by_label('Example sample table for a matrix').inner_text().splitlines() == [
                    'sample_id,condition', 'C1,Control', 'C2,Control', 'T1,Treatment', 'T2,Treatment']
                for mode in ['dda', 'baf', 'dia']:
                    page.locator('#mode').select_option(mode)
                    example = page.get_by_label('Example sample table for raw data')
                    assert example.is_visible()
                    assert example.inner_text().splitlines()[0] == 'sample_id,condition,raw_path'
                    assert not page.get_by_label('Example sample table for a matrix').is_visible()
                if args.output:
                    args.output.parent.mkdir(parents=True, exist_ok=True)
                    page.screenshot(path=str(args.output), full_page=True)
                page.locator('#input-options > summary').click()
                # The first Discover click opens a chooser when no data is selected.
                page.locator('#discover').click()
                if not args.native:
                    page.get_by_role('button', name='Open folder raw samples', exact=True).click()
                    page.locator('#choose-folder').click()
                page.wait_for_function("document.querySelectorAll('.sample').length === 4")
                assert page.locator('#project').input_value() == ''
                assert not list(raw.glob('config'))
                # A disconnected server must report recovery guidance and keep data.
                page.route('**/api/choose-folder', lambda route: route.abort('connectionrefused'))
                page.locator('[data-browse=raw_root]').click()
                page.wait_for_function('!operationBusy')
                assert 'Run fun_proteomics setup again' in page.locator('#message').inner_text()
                assert page.locator('#raw_root').input_value() == str(raw)
                assert page.locator('.sample').count() == 4
                page.unroute('**/api/choose-folder')
                if args.native:
                    assert not page.locator('#picker').is_visible()
                    command = json.loads(calls.read_text().splitlines()[-1])
                    assert '--directory' in command
                    assert '--filename=' + str(root) + os.sep in command
                    # Cancel must leave both the path and discovered board intact.
                    response.write_text(json.dumps({'code': 1}))
                    page.locator('[data-browse=raw_root]').click()
                    page.wait_for_function('!operationBusy')
                    assert page.locator('#raw_root').input_value() == str(raw)
                    assert page.locator('.sample').count() == 4
                    assert not page.locator('#picker').is_visible()
                    # Errors are visible and unlock the UI without clearing data.
                    response.write_text(json.dumps({'code': 2}))
                    page.locator('[data-browse=raw_root]').click()
                    page.wait_for_function('!operationBusy')
                    assert 'Could not open' in page.locator('#message').inner_text()
                    assert page.locator('#raw_root').input_value() == str(raw)
                    assert page.locator('.sample').count() == 4
                    assert not page.locator('#picker').is_visible()
                    response.write_text(json.dumps({'path': str(raw)}))
                    page.locator('[data-browse=raw_root]').click()
                    page.wait_for_function('!operationBusy')
                    assert 'Detected 4 samples' in page.locator('#message').inner_text()
                    # Existing analysis runs must also use the native folder dialog.
                    page.locator('#mode').select_option('continue')
                    page.locator('[data-browse=source_run]').click()
                    page.wait_for_function('!operationBusy')
                    assert page.locator('#source_run').input_value() == str(raw)
                    assert not page.locator('#picker').is_visible()
                # Change data type and select a matrix by navigating the file picker.
                page.locator('#mode').select_option('matrix')
                if args.native:
                    response.write_text(json.dumps({'path': str(root / 'inputs/proteins.csv')}))
                page.locator('[data-browse=matrix]').click()
                if args.native:
                    page.wait_for_function('!operationBusy')
                    command = json.loads(calls.read_text().splitlines()[-1])
                    assert '--directory' not in command
                    assert '--title=Choose a protein quantification matrix' in command
                    assert not page.locator('#picker').is_visible()
                else:
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
                if args.native:
                    (root / 'My analysis').mkdir()
                    response.write_text(json.dumps({'path': str(root / 'My analysis')}))
                page.locator('[data-browse=project]').click()
                if args.native:
                    page.wait_for_function('!operationBusy')
                    command = json.loads(calls.read_text().splitlines()[-1])
                    assert '--title=Choose or create output folder' in command
                else:
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
                if not args.native:
                    page.locator('#choose-folder').click()
                page.wait_for_function("document.querySelector('[data-panel=\"3\"]').hidden === false")
                assert page.locator('#project-name').inner_text() == str(root / 'My analysis')
                page.locator('[data-step="1"]').click()
                page.set_viewport_size({'width': 390, 'height': 844})
                assert page.evaluate('document.documentElement.scrollWidth <= window.innerWidth')
                assert not errors, errors
                browser.close()
                print('PASS: ' + ('native dialog selection/cancel/error, ' if args.native else 'headless folder fallback, ') + 'CSV examples, raw/blank scan, file search, automatic discovery, output folder, save, reopen and mobile layout')
        finally:
            process.send_signal(signal.SIGINT)
            process.wait(timeout=15)
            process.stdout.close()


if __name__ == '__main__':
    main()
