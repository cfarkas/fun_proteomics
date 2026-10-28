"""Readable dependency checks; installations use the existing pinned bootstrap."""
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

from . import __version__


def system_info(path):
    target=Path(path).expanduser().resolve()
    while not target.exists(): target=target.parent
    memory=0
    try:
        for line in Path('/proc/meminfo').read_text().splitlines():
            if line.startswith('MemAvailable:'): memory=int(line.split()[1])*1024
    except OSError: pass
    cpus=len(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity') else os.cpu_count()
    return dict(cpus=cpus,memory_gib=round(memory/2**30,1),disk_free_gib=round(shutil.disk_usage(target).free/2**30,1),path=str(target))


def readiness(project,args,mode='matrix',statistics=True,fdr_method='BH'):
    from .cli import environment, nextflow_binary
    env=environment(project,args)
    checks={name:importlib.util.find_spec(name) is not None for name in ('numpy','pandas','matplotlib','pyarrow','openpyxl','jsonschema')}
    missing=[name for name,found in checks.items() if not found]
    details=dict(version=__version__,python=sys.executable,workflow=str(project),python_dependencies=checks)
    if statistics:
        r=env.get('PROTEOMICS_RSCRIPT')
        if not r or not shutil.which(r): missing.append('Rscript')
        else:
            packages=['limma','jsonlite'] if mode!='dda' else ['edgeR','limma','jsonlite','statmod','data.table']
            if fdr_method!='BH': packages.append('fdrtool')
            expression='p<-c('+','.join('"'+p+'"' for p in packages)+');bad<-p[!vapply(p,requireNamespace,logical(1),quietly=TRUE)];'
            if mode=='dda':
                expression+='for(n in c("edgeR","limma")){if(!n %in% bad && packageVersion(n)<if(n=="edgeR")"4.4.2" else "3.62.2")bad<-c(bad,paste0(n," (outdated)"))};'
            expression+='cat(paste(bad,collapse=","));quit(status=as.integer(length(bad)>0))'
            try:
                result=subprocess.run([r,'--vanilla','-e',expression],env=env,text=True,capture_output=True,timeout=45)
                if result.returncode: missing.append('R packages: '+(result.stdout.strip() or 'Rscript check failed'))
            except (OSError,subprocess.SubprocessError) as error: missing.append('R check: '+str(error))
        details['Rscript']=r
    if mode=='baf':
        reader=project/'.runtime/baf-python/bin/python'
        if not reader.is_file(): missing.append('BAF reader (install --baf or setup --python PATH)')
        if not shutil.which('java',path=env.get('PATH')): missing.append('Java')
    if mode in ('dia','dda'):
        try: details['nextflow']=nextflow_binary(project,args,env)
        except ValueError: missing.append('Nextflow')
        java=Path(env.get('JAVA_HOME',''))/'bin/java'
        if not java.is_file() and not shutil.which('java',path=env.get('PATH')): missing.append('Java')
        if mode=='dia':
            try:
                from .study import imports
                imports(project)
                from dia_browser import python312
                local=project/'.runtime/dia-python/bin/python'
                if local.is_file():
                    probe=subprocess.check_output([str(local),'-c','import sys; print("%d.%d" % sys.version_info[:2])'],text=True,timeout=10).strip()
                    if probe!='3.12': raise ValueError('Local DIA Python is not 3.12')
                    details['dia_python']=str(local)
                else: details['dia_python']=python312()
            except ValueError: missing.append('Python 3.12 for DIA')
    details['missing']=missing
    details['ready']=not missing
    return details


def doctor(project,args):
    info=readiness(project,args,args.mode)
    if args.json: print(json.dumps(info,indent=2))
    else:
        print('fun_proteomics '+__version__+' | '+args.mode+' analysis')
        for key in ('python','workflow','Rscript','nextflow','dia_python'):
            if info.get(key): print(key+': '+str(info[key]))
        if info['missing']:
            print('Missing: '+', '.join(info['missing']))
            print('Install tools: fun_proteomics install --conda')
            print('Existing R only: fun_proteomics install --r-only')
        else: print('Analysis dependencies are ready. Raw search archives and inputs are checked during setup.')
    return 0 if info['ready'] else 2


def install(project,args):
    from .cli import environment
    if getattr(args,'baf',False):
        if args.conda or args.r_only: raise ValueError('Run --baf separately from --conda/--r-only')
        target=project/'.runtime/baf-python'
        subprocess.run([sys.executable,'-m','venv',str(target)],check=True)
        SDK_COMMIT='b9fdcca8e78693245af304a8561f5abf8e18e95c'
        return subprocess.run([str(target/'bin/python'),'-m','pip','install','numpy>=1.26,<3','pandas>=2.2,<3',
            'psims==1.3.5','git+https://github.com/gtluu/pyBaf2Sql@'+SDK_COMMIT]).returncode
    if args.conda and args.r_only: raise ValueError('Choose --conda or --r-only')
    if not args.conda and not args.r_only: raise ValueError('Use install --conda (full runtime) or install --r-only (existing R)')
    env=environment(project,args)
    if args.r_only:
        r=env.get('PROTEOMICS_RSCRIPT')
        if not r or not shutil.which(r): raise ValueError('Rscript is missing. Use install --conda to install R locally.')
        command=['bash',str(project/'scripts/install_differential_tools.sh'),r]
    else:
        command=['bash',str(project/'scripts/bootstrap_runtime.sh'),str(project/'.runtime'),
                 args.rscript or 'auto',str(project/'.runtime/runtime-settings.json'),'true']
    folder=project/'.runtime';folder.mkdir(parents=True,exist_ok=True)
    logfile=folder/'install.log'
    print('Installing analysis tools. Log: '+str(logfile),flush=True)
    with logfile.open('a') as log:
        process=subprocess.Popen(command,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        try:
            for line in process.stdout:
                print(line,end='',flush=True);log.write(line);log.flush()
            code=process.wait();process.stdout.close()
        except KeyboardInterrupt:
            process.terminate();process.wait();raise
    if code==0 and args.conda:
        from .study import imports
        imports(project)
        from dia_browser import python312
        try: python312()
        except ValueError:
            local=project/'.runtime/dia-python/bin/python'
            if not local.is_file():
                print('Installing Python 3.12 for native DIA searches.',flush=True)
                command=[str(project/'.runtime/bin/micromamba'),'--root-prefix',str(project/'.runtime/mamba'),'create','--yes','--prefix',str(local.parent.parent),'--channel','conda-forge','--strict-channel-priority','python=3.12']
                with logfile.open('a') as log:
                    process=subprocess.Popen(command,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
                    try:
                        with process.stdout:
                            for line in process.stdout: print(line,end='',flush=True);log.write(line);log.flush()
                        code=process.wait()
                    except KeyboardInterrupt: process.terminate();process.wait();raise
    if code==0: print('Installed. Run fun_proteomics doctor --mode matrix (or dia / dda).')
    return code
