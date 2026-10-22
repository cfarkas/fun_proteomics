#!/usr/bin/env python3
"""Use the real browser to group BAF technical injections and edit paired donors."""
import argparse
import os
from pathlib import Path
import signal
import sys
import tempfile
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from record_walkthroughs import launch


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--python',default=sys.executable);parser.add_argument('--browser-path');parser.add_argument('--source',action='store_true')
    args=parser.parse_args()
    from playwright.sync_api import sync_playwright
    with tempfile.TemporaryDirectory(prefix='proteomics-baf-browser-') as directory:
        root=Path(directory);raw=root/'raw';raw.mkdir()
        for name in ['D1_A_R1','D1_A_R2','D1_B_R1','D1_B_R2']:
            folder=raw/(name+'.d');folder.mkdir();(folder/'analysis.baf').write_bytes(b'synthetic fixture, not searchable')
        env=os.environ.copy()
        if args.source:env['PYTHONPATH']=str(ROOT/'src')
        else:env.pop('PYTHONPATH',None)
        env['FUN_PROTEOMICS_CACHE']=str(root/'cache')
        process,url=launch(args,root,env)
        try:
            with sync_playwright() as p:
                browser=p.chromium.launch(headless=True,**({'executable_path':args.browser_path} if args.browser_path else {}))
                page=browser.new_page(viewport={'width':1440,'height':1150});errors=[];page.on('pageerror',lambda error:errors.append(str(error)))
                page.goto(url);page.wait_for_function('document.getElementById("project").value.length>0')
                page.locator('input[name=mode][value=baf]').check();page.locator('#raw_root').fill(str(raw));page.locator('#discover').click()
                page.wait_for_selector('.sample');assert page.locator('.sample').count()==4
                for name in ['D1_A_R1','D1_A_R2','D1_B_R1','D1_B_R2']:
                    page.get_by_label('Condition for '+name,exact=True).select_option('Control' if '_A_' in name else 'Treatment')
                    page.get_by_label('Biological sample for '+name,exact=True).fill(name.rsplit('_',1)[0])
                page.locator('#next').click();page.locator('#paired_samples').check()
                assert page.locator('#service_matrix').is_visible()
                assert page.locator('#python').is_visible()
                page.locator('[data-step="1"]').click()
                for name in ['D1_A_R1','D1_A_R2','D1_B_R1','D1_B_R2']:page.get_by_label('Donor for '+name,exact=True).fill('D1')
                cfg=page.evaluate('collect()');assert cfg['subject_column']=='subject' and not cfg['independent_samples']
                rows=cfg['assignment_board']['rows'];assert len({r['biological_sample'] for r in rows})==2
                assert {r['subject'] for r in rows}=={'D1'}
                assert not errors,errors
                browser.close()
                print('PASS: Raw BAF discovery, condition assignment, technical-injection grouping, donor fields and service/output settings')
        finally:
            process.send_signal(signal.SIGINT);process.wait(timeout=15);process.stdout.close()


if __name__=='__main__':main()
