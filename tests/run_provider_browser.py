#!/usr/bin/env python3
"""Exercise Excel discovery, paired donor cards and the output preview in Chromium."""
import argparse
import json
import os
from pathlib import Path
import signal
import sys
import tempfile
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from record_walkthroughs import fixture, launch


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--python',default=sys.executable);parser.add_argument('--browser-path');parser.add_argument('--source',action='store_true');parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    from playwright.sync_api import sync_playwright
    with tempfile.TemporaryDirectory(prefix='proteomics-provider-browser-') as temp:
        root=Path(temp);inputs=fixture(root);env=os.environ.copy()
        if args.source: env['PYTHONPATH']=str(ROOT/'src')
        else: env.pop('PYTHONPATH',None)
        env.update(FUN_PROTEOMICS_CACHE=str(root/'cache'),PYTHONUNBUFFERED='1')
        data=pd.read_csv(inputs/'proteins.csv').rename(columns={'Protein.Group':'Protein IDs','Genes':'Gene names'})
        names=['C1','C2','A1','A2','B1','B2']
        data=data.rename(columns={name:'LFQ intensity '+name for name in names})
        book=inputs/'service.xlsx'
        with pd.ExcelWriter(book) as writer:
            pd.DataFrame({'Title':['Service report']}).to_excel(writer,sheet_name='Cover',index=False)
            data.to_excel(writer,sheet_name='LFQ proteins',index=False,startrow=2)
        process,url=launch(args,root,env)
        try:
            with sync_playwright() as p:
                browser=p.chromium.launch(headless=True,**({'executable_path':args.browser_path} if args.browser_path else {}))
                page=browser.new_page(viewport={'width':1440,'height':1150});errors=[];page.on('pageerror',lambda error:errors.append(str(error)))
                page.goto(url);page.wait_for_function('document.getElementById("project").value.length>0')
                assert str(root/'study/results/final_results/matrices') in page.locator('#output-preview').inner_text()
                page.locator('#matrix').fill(str(book));page.locator('#discover').click();page.wait_for_selector('.sample')
                assert page.locator('.sample').count()==6
                for name in names: page.get_by_label('Condition for '+name,exact=True).select_option('Control' if name.startswith('C') else 'Treatment')
                page.locator('#next').click();page.locator('#paired_samples').check()
                page.locator('[data-step="1"]').click()
                for name in names: page.get_by_label('Donor for '+name,exact=True).fill('D'+name[-1])
                # Two treatment measurements from one donor must be rejected.
                page.locator('#next').click();page.locator('#check').click()
                page.wait_for_function('document.getElementById("message").classList.contains("error")')
                assert 'replicate' in page.locator('#message').inner_text() or 'technical' in page.locator('#message').inner_text()
                page.locator('[data-step="1"]').click()
                page.locator('#new-condition').fill('Recovery');page.locator('#add-condition').click()
                for name in ['B1','B2']:page.get_by_label('Condition for '+name,exact=True).select_option('Recovery')
                if args.output: page.screenshot(path=str(args.output),full_page=True)
                page.locator('#next').click();page.locator('#save').click()
                page.wait_for_function('!document.getElementById("run").disabled')
                config=json.loads((root/'study/config/run.json').read_text())
                assert config['subject_column']=='subject' and not config['independent_samples']
                assert [r['subject'] for r in config['assignment_board']['rows']]==['D1','D2','D1','D2','D1','D2']
                page.locator('#run').click();page.wait_for_function('document.getElementById("state").textContent==="completed"',timeout=120000)
                final=root/'study/results/final_results';index=pd.read_csv(final/'results_index.csv')
                assert len(index)==6 and set(index.method)=={'paired limma-trend'}
                audit=json.loads((final/'quality_control/import_audit.json').read_text())
                assert audit['import']['worksheet']=='LFQ proteins' and audit['import']['header_row']==3
                assert (final/'quality_control/input_sample_mapping.tsv').is_file()
                page.reload();page.wait_for_function('!document.getElementById("run").disabled');page.locator('#run').click()
                page.wait_for_function('document.getElementById("log").textContent.includes("No analysis needs repeating")',timeout=60000)
                assert not errors,errors
                browser.close()
                print('PASS: Excel cover/header discovery, output preview, donor-card edits, technical-replicate rejection, paired contrasts, import audit, CSV/TSV and verified continuation')
        finally:
            process.send_signal(signal.SIGINT);process.wait(timeout=15);process.stdout.close()


if __name__=='__main__': main()
