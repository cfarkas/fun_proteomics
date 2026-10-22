"""One saved project shared by browser setup, terminal setup and scripted runs."""
import argparse
from collections import Counter
from contextlib import contextmanager
import copy
import csv
import fcntl
import html
import json
import os
from pathlib import Path
import shlex
import shutil
import sqlite3
import sys
import tempfile

from .continuation import read_json, write_json, run_root
from .project import digest

FORMAT = 1
MODES = ('matrix', 'dia', 'dda', 'continue')
PATH_FIELDS = ('matrix', 'samples', 'fasta', 'service_matrix', 'service_samples', 'source_run',
               'study_config', 'raw_root', 'msfragger_zip', 'ionquant_zip', 'diatracer_zip')
FIELDS = set(PATH_FIELDS) | {'fun_proteomics_project', 'project', 'mode', 'independent_samples',
         'differential', 'quantity', 'normalization', 'fdr', 'min_replicates', 'cpus', 'memory_gib',
         'python', 'sample_rows', 'assignment_board', 'contains_decoys', 'dda_plus'}


def add_commands(sub):
    p = sub.add_parser('setup', help='Open guided browser setup; add --terminal for SSH')
    p.add_argument('--project', dest='study_project', help='Folder containing config/ and results/')
    p.add_argument('--terminal', action='store_true', help='Answer questions in this terminal')
    p.add_argument('--non-interactive', action='store_true', help='Use explicit flags without prompting')
    p.add_argument('--run', action='store_true', help='Check and run after saving; resume if already saved')
    p.add_argument('--mode', choices=MODES)
    for field in PATH_FIELDS:
        p.add_argument('--'+field.replace('_','-'))
    p.add_argument('--independent-samples', dest='study_independent', action='store_true', default=None,
                   help='Confirm each sample is a different biological replicate')
    p.add_argument('--no-differential', action='store_true', help='Produce matrices and QC without statistical comparisons')
    p.add_argument('--quantity', choices=['auto','value','maxlfq','intensity','spectral_count'], default='auto')
    p.add_argument('--normalization', choices=['none','median'], default='none')
    p.add_argument('--cpus', type=int, default=min(8, os.cpu_count() or 1))
    p.add_argument('--memory-gib', type=int, default=32)
    p.add_argument('--port', type=int, default=0, help='Local browser port; default selects an available port')
    p.add_argument('--no-browser', action='store_true')
    c = sub.add_parser('check', help='Validate saved inputs and explain the planned comparisons')
    c.add_argument('--config', required=True)
    c.add_argument('--json', action='store_true')
    s = sub.add_parser('system', help='Show CPU, RAM and available disk before installation')
    s.add_argument('--path', default='.')
    s.add_argument('--json', action='store_true')
    s = sub.add_parser('status', help='Find the saved configuration, progress and result links')
    s.add_argument('--project', dest='study_project', required=True)
    s.add_argument('--json', action='store_true')
    d = sub.add_parser('demo', help='Copy small synthetic inputs, then use the normal setup workflow')
    d.add_argument('--project', dest='study_project', required=True)


def is_study(path):
    try:
        value = json.loads(Path(path).read_text())
        return isinstance(value, dict) and 'fun_proteomics_project' in value
    except (ValueError, OSError):
        return False


def normalize(data):
    if not isinstance(data, dict) or set(data)-FIELDS:
        raise ValueError('Unknown setup settings: '+str(sorted(set(data)-FIELDS) if isinstance(data,dict) else 'expected an object'))
    cfg = dict(fun_proteomics_project=FORMAT, mode='matrix', independent_samples=False, differential=True,
               quantity='auto', normalization='none', fdr=.05, min_replicates=2, cpus=min(8,os.cpu_count() or 1), memory_gib=32)
    cfg.update(data)
    if cfg['fun_proteomics_project'] != FORMAT:
        raise ValueError('Unsupported saved project format')
    if cfg['mode'] not in MODES:
        raise ValueError('Choose matrix, dia, dda or continue')
    for field in ('independent_samples','differential'):
        if not isinstance(cfg[field],bool): raise ValueError(field+' must be true or false')
    for field in ('project', *PATH_FIELDS):
        if cfg.get(field):
            value = str(cfg[field]).strip()
            if any(c in value for c in '\n\r\t'): raise ValueError('Paths cannot contain tabs or line breaks')
            cfg[field] = str(Path(value).expanduser().resolve())
    if not cfg.get('project'): raise ValueError('Choose a project folder')
    if cfg['quantity'] not in ('auto','value','maxlfq','intensity','spectral_count'):
        raise ValueError('Unsupported quantity type')
    return cfg


def load(path):
    path=Path(path).expanduser().resolve()
    cfg=normalize(read_json(path))
    if path != Path(cfg['project'])/'config/run.json':
        raise ValueError('Use the saved PROJECT/config/run.json at its recorded project location')
    return cfg


def imports(workflow):
    # The existing analysis implementation is also bundled in installed wheels.
    for folder in (workflow/'scripts', workflow/'pipeline_lib'):
        if str(folder) not in sys.path: sys.path.insert(0,str(folder))


def samples_for(cfg, workflow):
    imports(workflow)
    from quantification import read_table
    import pandas as pd
    if 'assignment_board' in cfg:
        from .assignments import rows_for
        frame=pd.DataFrame(rows_for(cfg)).fillna('').astype(str)
    elif cfg.get('samples'):
        frame=read_table(cfg['samples'])
    elif cfg.get('sample_rows'):
        frame=pd.DataFrame(cfg['sample_rows']).fillna('').astype(str)
    else:
        raise ValueError('Select a sample table, or scan acquisitions and assign their conditions')
    if 'sample_id' not in frame or frame.empty or frame.sample_id.eq('').any() or not frame.sample_id.is_unique:
        raise ValueError('The sample table needs nonempty, unique sample_id values')
    if cfg['independent_samples']:
        if 'replicate' in frame: frame['source_replicate']=frame['replicate']
        # The user has explicitly confirmed independence; retain the original labels.
        frame['replicate']=[str(i+1) for i in range(len(frame))]
    if cfg['mode'] in ('dia','dda'):
        if 'condition' not in frame or frame.condition.eq('').any(): raise ValueError('Assign a condition to every raw sample')
        if 'replicate' not in frame: frame['replicate']=[str(i+1) for i in range(len(frame))]
        if 'raw_path' not in frame: raise ValueError('Raw analysis sample tables need raw_path')
        base=Path(cfg.get('raw_root') or Path(cfg.get('samples',cfg['project'])).parent)
        frame['raw_path']=[str((base/Path(p).expanduser()).resolve()) for p in frame.raw_path]
        if len(set(frame.raw_path))!=len(frame): raise ValueError('A raw acquisition appears more than once')
        if 'cohort' not in frame: frame['cohort']=frame['sample_type'] if 'sample_type' in frame else 'study'
    return frame


def raw_study(cfg, samples, workflow):
    imports(workflow)
    if cfg['mode']=='dda':
        samples=samples[[c for c in samples if c in {'sample_id','cohort','condition','replicate','raw_path','library_reference','service_sample'}]].copy();samples['replicate']=samples['replicate'].astype(int)
    if 'library_reference' in samples:
        samples=samples.copy()
        if not samples.library_reference.str.lower().isin(['true','false']).all(): raise ValueError('library_reference must be true or false')
        samples['library_reference']=samples.library_reference.str.lower().eq('true')
    base={'outdir':str(Path(cfg['project'])/'results'),
          'resources':{'cpus':cfg['cpus'],'memory_gib':cfg['memory_gib']},
          'samples':samples.to_dict('records')}
    if cfg['mode']=='dia':
        base.update(proteomics={'acquisition':'DIA-PASEF','analysis':'global-label-free','confirmed':True},
             dia={'confirmed':True, 'canonical_fasta':cfg.get('fasta'),
                  'install_cache':str(workflow/'.runtime/dia-cache'),
                  'differential':dict(enabled=cfg['differential'],comparisons='all_pairwise',gene_level=True,
                       fdr=cfg['fdr'],min_replicates=cfg['min_replicates'],normalization=cfg['normalization'])})
        if cfg.get('raw_root'): base['dia']['raw_root']=cfg['raw_root']
        local=workflow/'.runtime/dia-python/bin/python'
        if cfg.get('python'): base['dia']['python']=cfg['python']
        elif local.is_file(): base['dia']['python']=str(local)
        if cfg.get('service_matrix'):
            groups=set(samples.cohort)
            if len(groups)!=1: raise ValueError('Use one cohort for a single service benchmark, or configure a native DIA JSON')
            base['dia']['service_matrices']={next(iter(groups)):cfg['service_matrix']}
        return base
    if cfg.get('service_matrix'):
        raise ValueError('For DDA, import the resulting quantity matrix with its service table in a separate matrix project')
    base.update(modes=['DDA','DDA+'] if cfg.get('dda_plus') else ['DDA'], organism=None,
                database={'fasta':cfg.get('fasta'),'sha256':digest(cfg['fasta']),'contains_decoys':True},
                tools={key:cfg[key] for key in ('msfragger_zip','ionquant_zip','diatracer_zip') if cfg.get(key)},
                differential={'enabled':cfg['differential'],'methods':['edgeR','limma'],'comparisons':'all_pairwise',
                              'fdr':cfg['fdr'],'workers':min(4,cfg['cpus'])},
                optional_analyses={'gene_counts':True})
    return base


def check(cfg, workflow):
    cfg=normalize(cfg); imports(workflow)
    from quantification_config import differential_options, validate_design_rows
    project=Path(cfg['project']); mode=cfg['mode']
    if project == workflow or project.is_relative_to(workflow):
        raise ValueError('Choose a study folder outside the software checkout/cache')
    if cfg.get('raw_root'):
        raw=Path(cfg['raw_root'])
        if project.is_relative_to(raw) or raw.is_relative_to(project):
            raise ValueError('Keep the raw-data folder and project in separate directory trees')
    if mode=='continue':
        if not cfg.get('source_run'): raise ValueError('Select an existing run folder')
        source=run_root(cfg['source_run'])
        if source==project or (source/'config/run.json').is_file():
            raise ValueError('For a guided project, open that project directly; select a legacy run/results folder here')
        if not source.is_dir(): raise ValueError('Run folder does not exist')
        from .continuation import dia_settings, latest_legacy_nextflow
        if not ((source/'.fun_proteomics/run.json').is_file() or (source/'final_results/FINAL_RESULTS_COMPLETE').is_file()
                or dia_settings(source) or latest_legacy_nextflow(source)):
            raise ValueError('No saved analysis or reusable search detected in this folder')
        return {'mode':mode,'source_run':str(source),'summary':'Continue the saved design and verify reusable outputs before executing.'}
    if cfg.get('study_config') and 'assignment_board' in cfg:
        raise ValueError('A native configuration supplies its own assignments; clear it to use the sample board')
    if cfg.get('study_config'):
        if not Path(cfg['study_config']).is_file(): raise ValueError('Native configuration file does not exist')
        native=read_json(cfg['study_config'])
        if mode=='matrix': raise ValueError('Native configurations are for DIA or DDA')
        if mode=='dia':
            from diapasef_config import normalize_config
            with tempfile.TemporaryDirectory(prefix='proteomics-check-') as tmp:
                # Only an absent software-cache directory is substituted for validation.
                original=Path(cfg['study_config'])
                # Native paths remain relative to their original JSON, never the project.
                for key in ('samplesheet','raw_root','canonical_fasta','contaminant_fasta','method_config','requirements_lock','install_cache'):
                    if native.get(key): native[key]=str((original.parent/Path(native[key]).expanduser()).resolve())
                native['service_matrices']={k:str((original.parent/Path(v).expanduser()).resolve()) for k,v in native.get('service_matrices',{}).items()}
                if not Path(native.get('install_cache','')).is_dir(): native['install_cache']=tmp
                temp=Path(tmp)/'dia.json';write_json(temp,native)
                audit=normalize_config(temp,workflow/'modules/diapasef',outdir=str(project/'results'))
            return {'mode':mode,'samples':sum(v['samples'] for v in audit['groups'].values()),'summary':'Use the supplied native DIA configuration and its saved comparison settings.'}
        from config_server import absolute_config, check as legacy_check
        native=absolute_config(Path(cfg['study_config']));native['outdir']=str(project/'results');legacy_check(native)
        return {'mode':mode,'samples':len(native['samples']),'summary':'Use the supplied native DDA configuration and its saved comparison settings.'}
    if cfg['differential'] and not cfg['independent_samples']:
        raise ValueError('Confirm that each sample is an independent biological replicate. Repeated patients need a tailored model; choose matrices only to omit testing.')
    frame=samples_for(cfg,workflow)
    options=differential_options(dict(enabled=cfg['differential'],comparisons='all_pairwise',
                  fdr=cfg['fdr'],min_replicates=cfg['min_replicates'],normalization=cfg['normalization']))
    for key in ('sample_type','cohort'):
        if mode=='matrix' and key in frame and frame[key].nunique()>1:
            raise ValueError('Use a separate matrix project for each '+key)
    groups=frame.groupby('cohort') if mode in ('dia','dda') else [('study',frame)]
    pairs=0; conditions={}
    for name, group in groups:
        validate_design_rows(group.to_dict('records'),options)
        counts=dict(Counter(group.condition)) if 'condition' in group else {}
        conditions[name]=counts
        pairs+=len(counts)*(len(counts)-1)//2 if cfg['differential'] else 0
    report=dict(mode=mode,samples=len(frame),conditions=conditions,comparisons=pairs,
                summary=f'{len(frame)} independent samples; {pairs} condition pairs at protein and available gene levels.' if cfg['independent_samples'] else f'{len(frame)} samples; matrices and quality control only.',warnings=[])
    if 'assignment_board' in cfg:
        report['unassigned']=len(cfg['assignment_board']['rows'])-len(frame)
        report['summary']+=f" {report['unassigned']} unassigned samples excluded."
    if mode=='matrix':
        from quantification import load_matrix, annotate_groups, gene_matrix
        if not cfg.get('matrix'): raise ValueError('Select the protein quantification matrix')
        matrix, annotations, mapping=load_matrix(cfg['matrix'],frame,cfg['quantity'])
        if matrix.empty: raise ValueError('No target protein groups remain in the matrix')
        if cfg.get('fasta'): annotations=annotate_groups(annotations.reset_index(),matrix.index,cfg['fasta'])
        if cfg['differential'] and (cfg['quantity']=='spectral_count' or mapping.quantity.eq('spectral_count').any()):
            raise ValueError('Intensity comparisons require intensities. Use counts for delivered DDA spectral counts, or disable differential testing.')
        report.update(protein_groups=len(matrix),genes=len(gene_matrix(matrix,annotations)))
        if not report['genes']: report['warnings'].append('No unambiguous gene symbols found. Supply the matching FASTA with GN= annotations for gene reports.')
        if cfg.get('service_matrix'):
            from quantification import read_table, benchmark_detection
            from .assignments import service_samples_for
            service_samples=service_samples_for(cfg,frame,workflow)
            service,_,_=load_matrix(cfg['service_matrix'],service_samples)
            report['benchmark']=benchmark_detection(matrix,service)[0]
            # Preview stays small and JSON serializable.
            report['benchmark']=json.loads(report['benchmark'].to_json(orient='records'))
        return report
    if not cfg.get('fasta') or not Path(cfg['fasta']).is_file(): raise ValueError('Select a matching FASTA database')
    native=raw_study(cfg,frame,workflow)
    if mode=='dia':
        from dia_browser import validate_browser_study
        report['raw']=validate_browser_study(native,workflow)
    else:
        import jsonschema
        try: jsonschema.validate(native,read_json(workflow/'assets/config.schema.json'))
        except jsonschema.ValidationError as error: raise ValueError('DDA settings: '+error.message) from error
        for key in ('msfragger_zip','ionquant_zip','diatracer_zip'):
            if not cfg.get(key) or not Path(cfg[key]).is_file(): raise ValueError('Select the registered '+key+' archive')
        from fasta_download import read_records
        records=read_records(cfg['fasta'])
        targets={h.split()[0] for h,_ in records if not h.startswith('rev_')}
        decoys={h.split()[0][4:] for h,_ in records if h.startswith('rev_')}
        if not targets or targets!=decoys or not any(h.startswith('contam_') for h in targets):
            raise ValueError('DDA setup expects a target/decoy search FASTA containing rev_ decoys; use the advanced raw editor to prepare one')
        for path in frame.raw_path:
            raw=Path(path)
            if raw.suffix.lower()!='.d' or any(not (raw/n).is_file() or not (raw/n).stat().st_size for n in ('analysis.tdf','analysis.tdf_bin')):
                raise ValueError('Incomplete Bruker acquisition: '+str(raw))
            with sqlite3.connect((raw/'analysis.tdf').as_uri()+'?mode=ro',uri=True) as db:
                if {v[0] for v in db.execute('SELECT DISTINCT MsMsType FROM Frames')}!={0,8}:
                    raise ValueError('Expected DDA-PASEF acquisition: '+str(raw))
    return report


@contextmanager
def locked(project):
    parent=Path(project)/'.fun_proteomics';parent.mkdir(parents=True,exist_ok=True)
    with (parent/'project.lock').open('a') as handle:
        try: fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError as error: raise ValueError('This project is already running or being saved') from error
        yield


def save(cfg, workflow):
    cfg=normalize(cfg);report=check(cfg,workflow)
    project=Path(cfg['project']);path=project/'config/run.json'
    if not path.is_file() and project.exists() and any(project.iterdir()):
        raise ValueError('Choose a new or empty project folder; to reopen a project, select one with config/run.json')
    with locked(project):
        if path.is_file():
            if load(path)!=cfg:
                # Drafts can be corrected, but a launched design remains immutable.
                if ((project/'.fun_proteomics/launch.json').exists()
                        or (project/'.fun_proteomics/browser-launch.json').exists()
                        or (project/'results').exists()):
                    raise ValueError('This project has already been launched. Use a new project for changed inputs or design.')
                from datetime import datetime, timezone
                backup=path.parent/'backups'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
                backup.mkdir(parents=True)
                for name in ('run.json','sample_assignments.csv','sample_assignments.tsv'):
                    if (path.parent/name).is_file(): shutil.copy2(path.parent/name,backup/name)
                write_json(path,cfg)
                if 'assignment_board' not in cfg:
                    for name in ('sample_assignments.csv','sample_assignments.tsv'):
                        (path.parent/name).unlink(missing_ok=True)
        else: write_json(path,cfg)
        if 'assignment_board' in cfg:
            import pandas as pd
            table=pd.DataFrame(cfg['assignment_board']['rows'])
            table['included']=table.condition.ne('')
            table=table[['included',*[c for c in table if c!='included']]]
            for suffix,sep in [('csv',','),('tsv','\t')]:
                table.to_csv(path.parent/('sample_assignments.'+suffix),sep=sep,index=False)
    return path,report


def result_location(cfg):
    root=run_root(cfg['source_run']) if cfg['mode']=='continue' else Path(cfg['project'])/'results'
    state=read_json(root/'.fun_proteomics/run.json') if (root/'.fun_proteomics/run.json').is_file() else {}
    return Path(state.get('outdir',root)),state


def status(project):
    config=Path(project).expanduser().resolve()/'config/run.json'
    cfg=load(config);root,state=result_location(cfg);final=root/'final_results'
    log=Path(cfg['project'])/'logs/analysis.log'
    if not log.is_file(): log=root/'.fun_proteomics/analysis.log'
    state_name=state.get('status','completed' if (final/'FINAL_RESULTS_COMPLETE').is_file() else 'saved')
    if state_name=='running':
        from .cli import process_identity
        if not state.get('process_identity') or process_identity(state.get('pid'))!=state['process_identity']:
            state_name='interrupted'
    return dict(project=cfg['project'],config=str(config),mode=cfg['mode'],state=state_name,
                results=str(final/'START_HERE.html') if (final/'START_HERE.html').is_file() else None,
                matrices=str(final/'MATRICES.html') if (final/'MATRICES.html').is_file() else None,
                log=str(log))


def runtime_flags(args):
    result=[]
    for key,flag in [('project','--workflow-dir'),('cache_dir','--cache-dir'),('rscript','--rscript'),('r_library','--r-library'),('nextflow','--nextflow')]:
        if getattr(args,key,None): result += [flag,str(getattr(args,key))]
    if getattr(args,'print_command',False): result.append('--print-command')
    if getattr(args,'hash_raw',False): result.append('--hash-raw')
    return result


def display(report):
    print(report['summary'])
    if 'protein_groups' in report: print(f"{report['protein_groups']} protein groups; {report['genes']} unambiguous genes.")
    for warning in report.get('warnings',[]): print('Note: '+warning)


def run_study(path,args,workflow,check_only=False):
    from . import cli
    cfg=load(path);report=check(cfg,workflow);display(report)
    if check_only: return 0
    project=Path(cfg['project']);config=project/'config/run.json';launch=project/'.fun_proteomics/launch.json'
    flags=runtime_flags(args)
    root,state=result_location(cfg)
    if cfg['mode']!='continue' and not state and not getattr(args,'print_command',False):
        from .requirements import readiness
        ready=readiness(workflow,args,cfg['mode'],cfg['differential'] or cfg['mode']=='dda')
        if ready['missing']:
            raise ValueError('Missing tools: '+', '.join(ready['missing'])+'. Run fun_proteomics install --conda, then doctor --mode '+cfg['mode'])
    with locked(project):
        if launch.is_file():
            saved=read_json(launch)
            if saved['config_sha256']!=digest(config): raise ValueError('Saved analysis settings changed; use a new project')
            for item in saved['inputs']:
                if not Path(item['path']).is_file() or digest(item['path'])!=item['sha256']:
                    raise ValueError('Recorded input changed: '+item['path']+'; use a new project')
        elif not getattr(args,'print_command',False):
            paths=[cfg[f] for f in PATH_FIELDS if cfg.get(f) and Path(cfg[f]).is_file()]
            write_json(launch,dict(config_sha256=digest(config),inputs=cli.recorded_files(paths)))
        if state or cfg['mode']=='continue':
            command=[*flags,'--continue',str(run_root(cfg['source_run']) if cfg['mode']=='continue' else project/'results')]
        elif cfg['mode']=='matrix':
            samples=samples_for(cfg,workflow);sheet=project/'config/samples.csv'
            if not getattr(args,'print_command',False): samples.to_csv(sheet,index=False)
            command=[*flags,'quantify','--matrix',cfg['matrix'],'--samples',str(sheet),'--outdir',str(root),
                     '--quantity',cfg['quantity'],'--normalization',cfg['normalization'],'--fdr',str(cfg['fdr']),
                     '--min-replicates',str(cfg['min_replicates'])]
            if cfg['differential']: command.append('--differential')
            for field in ('fasta','service_matrix'):
                if cfg.get(field): command += ['--'+field.replace('_','-'),cfg[field]]
            if cfg.get('service_samples'):
                from .assignments import service_samples_for
                provider=project/'config/service_samples.csv'
                if not getattr(args,'print_command',False): service_samples_for(cfg,samples,workflow).to_csv(provider,index=False)
                command += ['--service-samples',str(provider)]
        else:
            native=Path(cfg['study_config']) if cfg.get('study_config') else project/'config/native.json'
            if not cfg.get('study_config') and not getattr(args,'print_command',False):
                frame=samples_for(cfg,workflow);data=raw_study(cfg,frame,workflow)
                if cfg['mode']=='dia':
                    from dia_browser import browser_inputs, _write_inputs
                    data,rows,targets,provenance=browser_inputs(data,workflow)
                    Path(data['install_cache']).mkdir(parents=True,exist_ok=True)
                    with tempfile.TemporaryDirectory(prefix='proteomics-native-') as directory:
                        temporary=Path(directory)/native.name
                        _write_inputs(data,rows,temporary,targets,provenance,logical_path=native)
                        for item in Path(directory).iterdir():
                            dest=native.parent/item.name
                            if dest.exists() and dest.read_bytes()!=item.read_bytes(): raise ValueError('Generated native inputs changed: '+str(dest))
                            if not dest.exists(): shutil.copy2(item,dest)
                else: write_json(native,data)
            command=[*flags,'run','--dia-config' if cfg['mode']=='dia' else '--config',str(native),'--outdir',str(root)]
        print('Configuration: '+str(config),flush=True)
        if getattr(args,'print_command',False):
            print(shlex.join(['fun_proteomics',*command]));return 0
        code=cli.main(command)
        if code==0 and not getattr(args,'print_command',False):
            info=status(project)
            if info['results']:
                for name,key in [('START_HERE.html','results'),('MATRICES.html','matrices')]:
                    if not info[key]: continue
                    destination=Path(info[key]).as_uri()
                    (project/name).write_text('<!doctype html><meta charset="utf-8"><title>Proteomics results</title><meta http-equiv="refresh" content="0;url='+html.escape(destination,quote=True)+'"><a href="'+html.escape(destination,quote=True)+'">Open '+name+'</a>')
                print('Open results: '+str(project/'START_HERE.html'))
                print('Open matrices (CSV + TSV): '+str(project/'MATRICES.html'))
        return code


def scan(folder):
    root=Path(folder).expanduser().resolve()
    if not root.is_dir(): raise ValueError('Select an existing raw-data folder')
    rows=[]
    for path in sorted(root.iterdir()):
        if path.is_dir() and path.suffix.lower()=='.d' and 'BLANK' not in path.name.upper():
            rows.append(dict(sample_id=path.stem,condition='',raw_path=str(path),cohort='study'))
    if not rows: raise ValueError('No biological .d acquisitions found directly inside this folder')
    if len(rows)>2000: raise ValueError('Select a smaller folder (at most 2000 samples)')
    return rows


def dispatch(args,workflow):
    if args.command=='system':
        from .requirements import system_info
        info=system_info(args.path)
        print(json.dumps(info,indent=2) if args.json else f"{info['cpus']} available CPUs; {info['memory_gib']} GiB RAM; {info['disk_free_gib']} GiB free at {info['path']}.")
        return 0
    if args.command=='status':
        info=status(args.study_project)
        if args.json: print(json.dumps(info,indent=2))
        else:
            for key in ('state','config','results','matrices','log'):
                if info.get(key): print(key.capitalize()+': '+info[key])
        return 0
    if args.command=='check':
        report=check(load(args.config),workflow)
        print(json.dumps(report,indent=2) if args.json else 'Inputs and design checked.')
        if not args.json: display(report)
        return 0
    if args.command=='demo':
        target=Path(args.study_project).expanduser().resolve()
        # Keep demo inputs separate from the empty study directory expected by setup.
        inputs=target.parent/(target.name+'-inputs')
        if inputs.exists(): raise ValueError('Demo input folder already exists: '+str(inputs))
        inputs.mkdir(parents=True)
        for old,new in [('quantification_matrix.csv','proteins.csv'),('quantification_samples.csv','samples.csv'),('service_quantification.csv','service.csv')]:
            shutil.copy2(workflow/'examples'/old,inputs/new)
        print('Synthetic example inputs: '+str(inputs))
        print('Start: fun_proteomics setup --project '+shlex.quote(str(target)))
        print('Choose Quantification matrix and select proteins.csv, samples.csv and optional service.csv.')
        return 0
    if args.terminal or args.non_interactive:
        from .setup_terminal import setup
        return setup(args,workflow)
    from .setup_web import serve
    return serve(args,workflow)
