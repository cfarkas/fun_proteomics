#!/usr/bin/env python3
"""Exercise real pointer dragging, dropdowns, persistence and all-pairs execution.

Requires Playwright and Chromium. Uses only synthetic inputs and an installed
package (or the source package when --source is supplied).
"""
import argparse
import json
import os
from pathlib import Path
import signal
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from record_walkthroughs import fixture, launch, drag_sample as drag


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--python',default=sys.executable);parser.add_argument('--browser-path');parser.add_argument('--source',action='store_true');parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    from playwright.sync_api import sync_playwright
    with tempfile.TemporaryDirectory(prefix='proteomics-browser-') as temp:
        root=Path(temp);inputs=fixture(root);env=os.environ.copy()
        if args.source:env['PYTHONPATH']=str(ROOT/'src')
        else:env.pop('PYTHONPATH',None)
        env.update(FUN_PROTEOMICS_CACHE=str(root/'cache'),PYTHONUNBUFFERED='1')
        process,url=launch(args,root,env)
        try:
            with sync_playwright() as p:
                browser=p.chromium.launch(headless=True,**({'executable_path':args.browser_path} if args.browser_path else {}))
                page=browser.new_page(viewport={'width':1440,'height':1150});errors=[];page.on('pageerror',lambda error:errors.append(str(error)))
                page.goto(url);page.wait_for_function('document.getElementById("project").value.length>0')
                page.locator('#matrix').fill(str(inputs/'proteins.csv'));page.locator('#discover').click();page.wait_for_selector('.sample')
                assert page.locator('.sample').count()==6
                assert '0 / 6 assigned' in page.locator('#sample-count').inner_text()
                name=page.get_by_role('textbox',name='Rename condition Treatment',exact=True);name.fill('Disease');name.press('Tab')
                page.locator('#new-condition').fill('Recovery');page.locator('#add-condition').click()
                for sample,condition in [('C1','Control'),('C2','Control'),('A1','Disease'),('A2','Disease'),('B1','Recovery'),('B2','Recovery')]:drag(page,sample,condition)
                drag(page,'C1','');assert '5 / 6 assigned' in page.locator('#sample-count').inner_text();drag(page,'C1','Control')
                page.get_by_label('Condition for B2',exact=True).select_option('');assert '5 / 6 assigned' in page.locator('#sample-count').inner_text()
                page.get_by_label('Condition for B2',exact=True).select_option('Recovery')
                page.get_by_label('Sample name for C1',exact=True).fill('RenamedControl')
                page.locator('#sample-filter').fill('A1');assert page.locator('.sample:visible').count()==1;page.locator('#sample-filter').fill('')
                if args.output:args.output.parent.mkdir(parents=True,exist_ok=True);page.screenshot(path=str(args.output),full_page=True)
                page.locator('#next').click();page.locator('#independent_samples').check();page.locator('#service_matrix').fill(str(inputs/'service.csv'))
                page.locator('#check').click();page.wait_for_function('document.getElementById("preview").textContent.includes("3 condition pairs")')
                page.locator('#save').click();page.wait_for_function('!document.getElementById("run").disabled')
                cfg=json.loads((root/'study/config/run.json').read_text())
                assert [r['condition'] for r in cfg['assignment_board']['rows']]==['Control','Control','Disease','Disease','Recovery','Recovery']
                page.reload();page.wait_for_function('!document.getElementById("run").disabled')
                page.locator('[data-step="1"]').click();assert page.get_by_label('Sample name for C1',exact=True).input_value()=='RenamedControl'
                assert '6 / 6 assigned' in page.locator('#sample-count').inner_text()
                page.locator('[data-step="3"]').click();page.locator('#run').click()
                page.wait_for_function('document.getElementById("state").textContent==="completed"',timeout=120000)
                import pandas as pd
                samples=pd.read_csv(root/'study/config/samples.csv');assert samples.sample_id.tolist()[0]=='RenamedControl'
                assert samples.condition.tolist()==['Control','Control','Disease','Disease','Recovery','Recovery']
                assert len(list((root/'study/results/final_results/differential_expression').glob('*/*/all_results.csv')))==6
                a=pd.read_csv(root/'study/config/sample_assignments.csv');b=pd.read_csv(root/'study/config/sample_assignments.tsv',sep='\t');pd.testing.assert_frame_equal(a,b)
                page.goto(url.split('/#')[0]+'/results/MATRICES.html');assert 'CSV' in page.locator('body').inner_text()
                page.goto(url);page.wait_for_function('!document.getElementById("run").disabled');page.locator('#run').click()
                page.wait_for_function('document.getElementById("log").textContent.includes("No analysis needs repeating")',timeout=60000)
                # A new design in a new project must exclude unassigned cards from
                # actual matrices and the provider benchmark, not just the preview.
                page.locator('[data-step="1"]').click()
                page.locator('#project').fill(str(root/'subset-study'))
                for sample in ('B1','B2'):page.get_by_label('Condition for '+sample,exact=True).select_option('')
                page.locator('#next').click();page.locator('#save').click()
                page.wait_for_function('document.getElementById("state").textContent==="saved"')
                page.locator('#run').click();page.wait_for_function('document.getElementById("state").textContent==="completed"',timeout=120000)
                subset=root/'subset-study/results/final_results'
                quantities=pd.read_csv(subset/'matrices/protein_intensities.csv')
                assert 'B1' not in quantities and 'B2' not in quantities and 'RenamedControl' in quantities
                assert len(list((subset/'differential_expression').glob('*/*/all_results.csv')))==2
                benchmark=pd.read_csv(subset/'benchmark/detection_summary.csv')
                assert set(benchmark.sample_id)=={'RenamedControl','C2','A1','A2','__union__'}
                assert not errors,errors
                browser.close()
                print('PASS: pointer drag, unassign, dropdown, rename, reload, six protein/gene contrasts, provider benchmark, CSV/TSV, verified continuation and an excluded-sample rerun')
        finally:
            process.send_signal(signal.SIGINT);process.wait(timeout=15);process.stdout.close()


if __name__=='__main__':main()
