"""Small, explicit terminal wizard; every answer is saved before execution."""
from pathlib import Path
import shlex
from . import study


def ask(label,default=None,required=False):
    while True:
        value=input(label+(f' [{default}]' if default not in (None,'') else '')+': ').strip()
        value=value or default
        if value or not required: return value or ''
        print('Please enter a value.')


def yes(label,default=False):
    while True:
        value=ask(label+' (yes/no)','yes' if default else 'no').lower()
        if value in ('yes','y'): return True
        if value in ('no','n'): return False
        print('Enter yes or no.')


def from_args(args):
    cfg={field:getattr(args,field,None) for field in study.PATH_FIELDS if getattr(args,field,None)}
    cfg.update(project=getattr(args,'study_project',None),mode=getattr(args,'mode',None) or 'matrix',
               independent_samples=getattr(args,'study_independent',False) or False,
               differential=not getattr(args,'no_differential',False),quantity=getattr(args,'quantity','auto'),
               normalization=getattr(args,'normalization','none'),cpus=getattr(args,'cpus',8),memory_gib=getattr(args,'memory_gib',32))
    for field in ('sheet','service_sheet','subject_column'):
        if getattr(args,field,None): cfg[field]=getattr(args,field)
    return cfg


def setup(args,workflow):
    interactive=not args.non_interactive
    print('\nfun_proteomics | Set up an analysis\n')
    project=args.study_project or (ask('Output/project folder',required=True) if interactive else None)
    if not project: raise ValueError('Supply --project for non-interactive setup')
    path=Path(project).expanduser().resolve()/'config/run.json'
    if path.is_file():
        cfg=study.load(path)
        supplied=from_args(args)
        for key in study.PATH_FIELDS:
            if supplied.get(key) and str(Path(supplied[key]).expanduser().resolve())!=cfg.get(key):
                raise ValueError('Saved projects reuse their inputs; choose a new project for changed settings')
        if args.mode and args.mode!=cfg['mode']: raise ValueError('The saved project has a different mode; choose a new project')
        print('Loaded saved settings: '+str(path))
        study.display(study.check(cfg,workflow))
    else:
        cfg=from_args(args);cfg['project']=project
        if interactive and not args.mode:
            print('1. Quantification matrix (DIA-NN, FragPipe or service table)')
            print('2. DIA-PASEF raw acquisitions')
            print('3. DDA-PASEF raw acquisitions')
            print('4. Continue an existing run folder')
            value=ask('Analysis','1')
            cfg['mode']={'1':'matrix','2':'dia','3':'dda','4':'continue'}.get(value,value)
            if cfg['mode'] not in study.MODES: raise ValueError('Choose 1, 2, 3 or 4')
        mode=cfg['mode']
        fields=['source_run'] if mode=='continue' else ['study_config'] if cfg.get('study_config') else ['matrix','samples'] if mode=='matrix' else ['raw_root','samples','fasta']
        labels={'source_run':'Existing run folder','matrix':'Protein matrix (CSV, TSV or XLSX)',
                'samples':'Sample table (sample_id, condition; raw_path for acquisitions)',
                'fasta':'FASTA (DIA targets; DDA search database with contaminants and rev_ decoys)',
                'raw_root':'Raw-data folder containing .d acquisitions'}
        for field in fields:
            if not cfg.get(field) and interactive:
                cfg[field]=ask(labels.get(field,field),required=field not in ('raw_root',))
        if mode!='continue' and not cfg.get('study_config'):
            if interactive:
                if mode=='matrix' and not cfg.get('fasta'):
                    cfg['fasta']=ask('Matching FASTA for missing gene symbols (optional)')
                if mode!='dda' and not cfg.get('service_matrix'):
                    cfg['service_matrix']=ask('Service matrix for a detection benchmark (optional)')
                if not args.no_differential:
                    cfg['differential']=yes('Compare every condition against every other condition?',True)
                if cfg['differential'] and not cfg['independent_samples'] and not cfg.get('subject_column'):
                    cfg['independent_samples']=yes('Does each sample come from a different patient or biological replicate?')
                    if not cfg['independent_samples'] and mode=='matrix':
                        cfg['subject_column']=ask('Donor/subject column for paired comparisons', 'subject', required=True)
                if mode=='dda':
                    for key in ('msfragger_zip','ionquant_zip','diatracer_zip'):
                        if not cfg.get(key): cfg[key]=ask('Registered '+key+' archive',required=True)
                if mode in ('dia','dda'):
                    cfg['cpus']=int(ask('CPU threads',str(cfg['cpus'])))
                    cfg['memory_gib']=int(ask('Memory budget (GiB)',str(cfg['memory_gib'])))
        path,report=study.save(cfg,workflow)
        study.display(report)
        print('Saved: '+str(path))
        print('Output folder: '+str(path.parent.parent))
        print('Matrices and comparisons: '+str(path.parent.parent/'results/final_results'))
    if args.run:
        action='run'
    elif interactive:
        action=ask('Next action (save/run)','save').lower()
        if action not in ('save','run'): raise ValueError('Choose save or run')
    else: action='save'
    if action=='run': return study.run_study(path,args,workflow)
    print('Check: fun_proteomics check --config '+shlex.quote(str(path)))
    print('Run:   fun_proteomics run --config '+shlex.quote(str(path)))
    print('Later: fun_proteomics --continue '+shlex.quote(str(path.parent.parent)))
    return 0
