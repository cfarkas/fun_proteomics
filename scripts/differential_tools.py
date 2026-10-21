"""Probe optional R differential-analysis dependencies without installing globally."""
import json,os,shutil,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def capabilities():
    executable=os.environ.get('PROTEOMICS_RSCRIPT') or shutil.which('Rscript')
    if not executable:return {'available':False,'methods':{},'message':'Rscript not found'}
    env=os.environ.copy();env['R_LIBS']=str(ROOT/'.runtime/R-library')
    code='pk<-c("edgeR","edgeRun","limma","jsonlite","statmod");for(p in pk)cat(p,if(requireNamespace(p,quietly=TRUE))as.character(packageVersion(p)) else "missing",sep="=",fill=TRUE)'
    result=subprocess.run([executable,'-e',code],capture_output=True,text=True,env=env,timeout=30)
    versions=dict(line.split('=',1) for line in result.stdout.splitlines() if '=' in line)
    common=all(versions.get(p,'missing')!='missing' for p in ['edgeR','limma','jsonlite','statmod'])
    common=common and tuple(int(v) for v in versions['edgeR'].split('.')) >= (4,4,2) and tuple(int(v) for v in versions['limma'].split('.')) >= (3,62,2)
    return {'available':result.returncode==0 and common,'rscript':executable,'methods':{p:common and versions.get(p,'missing')!='missing' for p in ['edgeR','edgeRun','limma']},'versions':versions}
