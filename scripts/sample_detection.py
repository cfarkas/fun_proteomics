"""Discover Bruker acquisitions and suggest editable sample annotations."""
import os
from pathlib import Path
import re
import time


MS_FORMATS={'.raw':'Thermo / Waters RAW', '.mzml':'mzML', '.mzxml':'mzXML', '.wiff':'SCIEX WIFF', '.wiff2':'SCIEX WIFF2'}

def ms_format(path):
    return MS_FORMATS.get(Path(path).suffix.lower())

def infer(path):
    name=path.stem
    # Remove acquisition plate position, injection index and acquisition ID.
    stem=re.sub(r'_S\d+-[A-Z]\d+(?:_\d+)*$','',name,flags=re.I)
    if '_CF_' in stem: stem=stem.split('_CF_',1)[1]
    else: stem=re.sub(r'^(?:(?:DDA|PASEF|EVOSEP|\d+SPDWZ)[_-])+','',stem,flags=re.I)
    code=re.match(r'^([A-Za-z]+)(\d+)(?:_(.*))?$',stem)
    repmatch=re.search(r'(?:^|[_ -])(?:rep(?:licate)?|r)[_ -]?(\d+)(?=$|[_ -])',stem,re.I)
    trailing=re.search(r'[_ -](\d+)$',stem)
    notes=[]
    if code:
        replicate=int(code[2]);condition=code[3] or code[1]
    elif repmatch:
        replicate=int(repmatch[1]);condition=(stem[:repmatch.start()]+' '+stem[repmatch.end():]).strip('_ -')
    elif trailing:
        replicate=int(trailing[1]);condition=stem[:trailing.start()]
    else:
        replicate=None;condition=stem;notes.append('Replicate not encoded in name; provisional sequence number assigned')
    condition=re.sub(r'[_]+',' ',condition).strip() or 'Unassigned'
    if replicate is not None and replicate<1:replicate=None;notes.append('Replicate zero needs review')
    # Explicit tissue tokens only: never assume an unfamiliar sample is a cell.
    tokens=set(re.findall(r'[a-z]+',path.stem.lower()+' '+path.parent.name.lower()))
    cohorts=set()
    if tokens & {'cell','cells','lysate','lysates'}:cohorts.add('cell')
    if 'serum' in tokens:cohorts.add('serum')
    if 'plasma' in tokens:cohorts.add('plasma')
    cohort=next(iter(cohorts)) if len(cohorts)==1 else 'unassigned'
    if cohort=='unassigned':notes.append('Cohort cannot be inferred reliably; set it before running')
    return condition,replicate,cohort,notes


def detect(parent, existing=None, recursive=False):
    if not isinstance(parent,(str,Path)) or not str(parent).strip():
        raise ValueError('Choose an acquisition folder first')
    parent=Path(str(parent).strip()).expanduser()
    if not parent.is_absolute():raise ValueError('Choose an absolute server folder path')
    parent=parent.resolve()
    if not parent.is_dir():raise ValueError('Choose an existing server folder')
    matches=[];skipped=[];unsupported=[];started=time.monotonic();visited=0
    def consider(p):
        if all((p/n).is_file() for n in ('analysis.tdf','analysis.tdf_bin')):matches.append(p.resolve())
        else:skipped.append({'path':str(p),'reason':'Missing analysis.tdf or analysis.tdf_bin'})
    if parent.suffix.lower()=='.d':consider(parent)
    else:
        def walk_error(error):raise ValueError(f'Cannot scan directory: {error}')
        for base,dirs,files in os.walk(parent,followlinks=False,onerror=walk_error):
            visited+=1
            if visited>10000 or time.monotonic()-started>10:raise ValueError('Folder is too broad; choose a more specific raw-data directory')
            for name in sorted(files):
                p=Path(base)/name
                if ms_format(p):unsupported.append({'path':str(p.resolve()),'format':ms_format(p)})
            remaining=[]
            for name in sorted(dirs):
                p=Path(base)/name
                if name.lower().endswith('.d'):consider(p)
                elif ms_format(p):unsupported.append({'path':str(p.resolve()),'format':ms_format(p)})
                elif recursive and not name.startswith('.') and not p.is_symlink():remaining.append(name)
            dirs[:]=remaining
            if len(matches)+len(skipped)+len(unsupported)>2000:raise ValueError('More than 2000 acquisitions; choose a smaller folder')
    known={str(Path(s['raw_path']).expanduser().resolve()):s for s in existing or []}
    samples=[];suggestions=[];ids=set();used={}
    for path in sorted(set(matches)):
        old=known.get(str(path))
        if old:
            row={k:old[k] for k in ['sample_id','condition','replicate','cohort','raw_path']};row['raw_path']=str(path)
            notes=['Kept your existing assignments'];source='existing'
        else:
            condition,replicate,cohort,notes=infer(path);key=(cohort,condition)
            if replicate is None:
                replicate=1
                while replicate in used.get(key,set()):replicate+=1
            stem=re.sub(r'[^A-Za-z0-9_]+','_',condition).strip('_') or 'Sample'
            if not stem[0].isalpha():stem='S_'+stem
            row={'sample_id':stem+'_'+str(replicate),'condition':condition,'replicate':replicate,'cohort':cohort,'raw_path':str(path)};source='filename'
        original=row['sample_id'];suffix=2
        while row['sample_id'] in ids:
            row['sample_id']=original+'_'+str(suffix);suffix+=1
        ids.add(row['sample_id'])
        key=(row['cohort'],row['condition']);used.setdefault(key,set()).add(row['replicate'])
        samples.append(row);suggestions.append({'sample_id':row['sample_id'],'source':source,'notes':notes})
    return {'samples':samples,'suggestions':suggestions,'skipped':skipped,'unsupported':unsupported,'parent':str(parent)}
