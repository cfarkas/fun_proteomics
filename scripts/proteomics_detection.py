"""Identify acquisition methods from native Bruker frame metadata, read-only."""
from pathlib import Path
import sqlite3

FRAMEWORKS={
 'DDA-PASEF': {'name':'FragPipe DDA-PASEF label-free proteomics','engine':'MSFragger → MSBooster/Percolator → ProteinProphet → IonQuant','available':True},
 'DIA-PASEF': {'name':'DIA-NN DIA-PASEF quantification','engine':'DIA-NN library prediction, cohort-specific empirical libraries and protein quantification through main.nf.','available':True},
 'MS1-only': {'name':'MS1-only acquisition','engine':'No tandem-MS identification workflow is configured for these files.','available':False},
 'Mixed': {'name':'Mixed acquisition methods','engine':'Split acquisitions by method before analysis.','available':False},
 'Unknown': {'name':'Unrecognized acquisition','engine':'Inspect the acquisition metadata before selecting an analysis workflow.','available':False},
}

def inspect(path):
    path=Path(path).expanduser().resolve()
    result={'path':str(path),'acquisition':'Unknown','frame_types':{},'instrument':'','method':'','error':''}
    if path.suffix.lower()!='.d' or not (path/'analysis.tdf_bin').is_file():
        result['error']='Not a complete native Bruker .d acquisition';return result
    try:
        with sqlite3.connect((path/'analysis.tdf').as_uri()+'?mode=ro',uri=True,timeout=3) as db:
            db.execute('PRAGMA query_only=ON')
            result['frame_types']={int(k):int(v) for k,v in db.execute('SELECT MsMsType, COUNT(*) FROM Frames GROUP BY MsMsType')}
            types=set(result['frame_types'])
            if types=={0,8}:result['acquisition']='DDA-PASEF'
            elif types=={0,9}:result['acquisition']='DIA-PASEF'
            elif types=={0}:result['acquisition']='MS1-only'
            elif 8 in types and 9 in types:result['acquisition']='Mixed'
            tables={r[0] for r in db.execute("SELECT name FROM sqlite_master WHERE type='table'")}
            if 'GlobalMetadata' in tables:
                meta=dict(db.execute("SELECT Key,Value FROM GlobalMetadata WHERE Key IN ('InstrumentName','MethodName')"))
                result.update(instrument=meta.get('InstrumentName',''),method=meta.get('MethodName',''))
    except (sqlite3.Error,OSError,ValueError) as e:result['error']=str(e)
    return result

def detect_proteomics(samples):
    if not samples:raise ValueError('Select samples before detecting the acquisition method')
    rows=[dict(inspect(s['raw_path']),sample_id=s['sample_id']) for s in samples]
    methods={r['acquisition'] for r in rows}
    detected=next(iter(methods)) if len(methods)==1 else 'Mixed'
    return {'detected':detected,'samples':rows,'framework':FRAMEWORKS[detected],
            'analysis_note':'Acquisition metadata identifies DDA/DIA. It does not establish enrichment, labeling, or the biological analysis goal; choose that explicitly.'}

def confirm_config(cfg):
    choice=cfg.get('proteomics',{})
    if choice.get('confirmed') is not True:raise ValueError('Confirm the acquisition method and analysis type in the Proteomics step.')
    detection=detect_proteomics(cfg['samples'])
    if choice.get('acquisition')!=detection['detected'] and choice.get('manual_override') is not True:raise ValueError('Selected acquisition type does not match the raw metadata. Detect again and review your samples.')
    if choice.get('acquisition')=='DIA-PASEF':
        if detection['detected']!='DIA-PASEF':
            raise ValueError('DIA-PASEF requires native DIA raw metadata for every sample; a manual override cannot change acquisition data.')
        if choice.get('analysis')!='global-label-free':
            raise ValueError('DIA-PASEF currently supports global-label-free quantification; the selected optional search is not implemented.')
        return detection
    if choice.get('acquisition')!='DDA-PASEF' or choice.get('analysis') not in ('global-label-free','phosphoproteomics','immunopeptidomics'):
        raise ValueError('The selected workflow adapter is not implemented; choose a supported DDA-PASEF analysis.')
    if choice.get('analysis')=='phosphoproteomics':
        from analysis_modules import status
        module=status()
        if not module['installed'] or choice.get('module_sha256')!=module['sha256']:
            raise ValueError('Install the phosphoproteomics module and confirm it again.')
    return detection
