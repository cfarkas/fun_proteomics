"""Native BAF conversion; vendor libraries stay in the user's separate runtime."""
import base64
from collections import Counter
import importlib.metadata
from pathlib import Path
import sqlite3
import subprocess
import xml.etree.ElementTree as ET
import zlib
import numpy as np
from pipeline import digest, read, save

CONVERTER_VERSION = 2
SDK_COMMIT = 'b9fdcca8e78693245af304a8561f5abf8e18e95c'


def source_files(raw):
    raw = Path(raw)
    if raw.suffix.lower() != '.d' or not (raw/'analysis.baf').is_file():
        raise ValueError('Expected a Bruker .d acquisition containing analysis.baf: '+str(raw))
    if (raw/'analysis.tdf').exists(): raise ValueError('Mixed BAF/TDF acquisition is ambiguous: '+str(raw))
    return sorted([*raw.glob('analysis.baf*'), *raw.glob('*.ami')])


def validate_mzml(path):
    """Read the entire file, decode every binary array and check MS1/MS2 metadata."""
    counts = Counter(); ids = set(); declared = None; previous = -1.
    for event, node in ET.iterparse(path, events=('start', 'end')):
        tag = node.tag.rsplit('}', 1)[-1]
        if event == 'start':
            if tag == 'spectrumList': declared = int(node.attrib['count'])
            continue
        if tag != 'spectrum': continue
        cv = {p.get('accession'): p for p in node.iter() if p.tag.endswith('}cvParam')}
        level = int(cv['MS:1000511'].get('value'))
        if level not in (1, 2) or 'MS:1000127' not in cv: raise ValueError('Expected centroided MS1/MS2 spectra')
        if int(node.attrib['defaultArrayLength'])==0: raise ValueError('Empty centroid scans are incompatible with the search reader')
        sid = node.attrib['id']
        if sid in ids: raise ValueError('Duplicate spectrum ID')
        ids.add(sid)
        rt = float(cv['MS:1000016'].get('value'))
        if cv['MS:1000016'].get('unitAccession') == 'UO:0000031': rt *= 60
        if not np.isfinite(rt) or rt < previous: raise ValueError('Invalid retention time order')
        previous = rt
        if level == 2:
            if 'MS:1000744' not in cv or float(cv['MS:1000744'].get('value')) <= 0:
                raise ValueError('MS2 spectrum missing precursor m/z')
            precursor = next(p for p in node.iter() if p.tag.endswith('}precursor'))
            if precursor.get('spectrumRef') not in ids: raise ValueError('Missing parent spectrum')
        arrays = {}
        for item in node.iter():
            if not item.tag.endswith('}binaryDataArray'): continue
            terms = {p.get('accession') for p in item if p.tag.endswith('}cvParam')}
            kind = 'mz' if 'MS:1000514' in terms else 'intensity' if 'MS:1000515' in terms else None
            if not kind: continue
            payload = next(p.text or '' for p in item if p.tag.endswith('}binary'))
            decoded = base64.b64decode(payload, validate=True)
            if 'MS:1000574' in terms: decoded = zlib.decompress(decoded)
            values = np.frombuffer(decoded, dtype='<f8' if 'MS:1000523' in terms else '<f4')
            if len(values) != int(node.attrib['defaultArrayLength']) or not np.isfinite(values).all():
                raise ValueError('Corrupt mzML binary array')
            if (values < 0).any(): raise ValueError('Negative m/z or intensity')
            arrays[kind] = values
        if set(arrays) != {'mz', 'intensity'}: raise ValueError('Missing spectral arrays')
        if np.any(np.diff(arrays['mz']) < 0): raise ValueError('Unsorted spectrum')
        counts[level] += 1
        node.clear()
    if declared != sum(counts.values()) or not counts[1] or not counts[2]:
        raise ValueError('Truncated mzML or missing MS1/MS2 spectra')
    return {'ms1_spectra': counts[1], 'ms2_spectra': counts[2], 'last_rt_seconds': previous}


def nonempty_spectra(spectra, keys, length):
    """Exclude genuinely empty scans, without inventing peaks or parent references."""
    retained=[]; omitted=[]
    for row in spectra:
        if row['LineMzId'] is None or row['LineIntensityId'] is None:
            raise ValueError('Vendor centroid arrays missing; cannot substitute profile data')
        lengths=[length(int(row[key])) for key in ('LineMzId','LineIntensityId')]
        if lengths[0]!=lengths[1]: raise ValueError('BAF m/z and intensity array lengths differ')
        if not lengths[0]:
            omitted.append(dict(spectrum_id=row['Id'],ms_level=int(keys[row['AcquisitionKey']]['MsLevel'])+1,reason='both centroid arrays empty'))
        else: retained.append(row)
    empty_ids={r['spectrum_id'] for r in omitted}
    if any(r['Parent'] in empty_ids for r in retained):
        raise ValueError('A nonempty MS2 spectrum refers to an empty MS1 scan; cannot safely omit that parent')
    return retained, omitted


def convert(raw, output):
    """Convert an isolated acquisition copy, with atomic completion and full hashes."""
    from pyBaf2Sql.init_baf2sql import init_baf2sql_api
    from pyBaf2Sql.baf import open_storage, read_double, get_num_elements, get_sqlite_cache_filename_v2
    from pyBaf2Sql.error import throw_last_baf2sql_error
    from psims.mzml import MzMLWriter
    raw, output = Path(raw).resolve(), Path(output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    api = init_baf2sql_api()
    api.baf2sql_set_num_threads(1)
    inputs = {p.name: digest(p) for p in source_files(raw)}
    identity = dict(converter=CONVERTER_VERSION, inputs=inputs, sdk_sha256=digest(api._name),
                    psims=importlib.metadata.version('psims'), pyBaf2Sql=importlib.metadata.version('pyBaf2Sql'),
                    centroid=True, raw_calibration=False)
    marker = output.with_suffix('.conversion.json')
    if marker.is_file():
        old = read(marker)
        if old['identity'] != identity: raise ValueError('BAF input/converter changed; use a new project')
        if output.is_file() and digest(output) == old['mzml_sha256']:
            print('Verified converted acquisition: '+output.name, flush=True)
            return old
    staging = output.parent/(output.stem+'.source.d')
    staging.mkdir(exist_ok=True)
    for p in source_files(raw):
        dest = staging/p.name
        if not dest.is_file() or digest(dest) != inputs[p.name]:
            subprocess.run(['cp', '--reflink=auto', str(p), str(dest)], check=True)
        if digest(dest) != inputs[p.name]: raise ValueError('BAF staging integrity failure')
    handle = open_storage(api, str(staging), False)
    if not handle: throw_last_baf2sql_error(api)
    try:
        get_sqlite_cache_filename_v2(api, str(staging), False)
        with sqlite3.connect((staging/'analysis.sqlite').as_uri()+'?mode=ro', uri=True) as db:
            db.row_factory = sqlite3.Row
            keys = {r['Id']: dict(r) for r in db.execute('SELECT * FROM AcquisitionKeys')}
            spectra = [dict(r) for r in db.execute('SELECT * FROM Spectra ORDER BY Id')]
            variables = {(r['Spectrum'], r['Variable']): r['Value'] for r in db.execute('SELECT * FROM Variables WHERE Variable IN (5,6,7,8)')}
            steps = {}
            for r in db.execute('SELECT * FROM Steps'): steps.setdefault(r['TargetSpectrum'], []).append(dict(r))
        if any(k['ScanMode'] not in (0, 2) for k in keys.values()):
            raise ValueError('This route supports conventional DDA BAF only, not bbCID/DIA')
        spectra, omitted = nonempty_spectra(spectra, keys, lambda key:get_num_elements(api,handle,key))
        temporary = output.with_suffix('.partial.mzML')
        counts = Counter()
        with MzMLWriter(str(temporary)) as writer:
            writer.controlled_vocabularies()
            writer.file_description(['MS1 spectrum', 'MSn spectrum', 'centroid spectrum'],
                [dict(id='BAF', name='analysis.baf', location=raw.as_uri(), params=['Bruker BAF format'])])
            writer.software_list([dict(id='fun_proteomics', version=str(CONVERTER_VERSION), params=['custom unreleased software tool'])])
            writer.instrument_configuration_list([writer.InstrumentConfiguration(id='instrument', component_list=[
                writer.Source(1, ['electrospray ionization']), writer.Analyzer(2, ['quadrupole']),
                writer.Analyzer(3, ['time-of-flight']), writer.Detector(4, ['microchannel plate detector'])])])
            writer.data_processing_list([writer.DataProcessing([writer.ProcessingMethod(
                order=0, software_reference='fun_proteomics', params=['Conversion to mzML'])], id='conversion')])
            with writer.run(id=output.stem, instrument_configuration='instrument'):
                with writer.spectrum_list(count=len(spectra), data_processing_method='conversion'):
                    for row in spectra:
                        sid = row['Id']; key = keys[row['AcquisitionKey']]; level = int(key['MsLevel'])+1
                        if level not in (1, 2): raise ValueError('Unsupported MS level')
                        if row['LineMzId'] is None or row['LineIntensityId'] is None:
                            raise ValueError('Vendor centroid arrays missing; cannot substitute profile data')
                        mz = read_double(api, handle, int(row['LineMzId']))
                        intensity = read_double(api, handle, int(row['LineIntensityId']))
                        precursor = None
                        if level == 2:
                            candidates = steps.get(sid, [])
                            if len(candidates) != 1: raise ValueError('Ambiguous/missing DDA precursor')
                            target = float(variables[(sid, 7)]); width = float(variables[(sid, 8)])
                            precursor = dict(mz=float(candidates[0]['Mass']),
                                spectrum_reference='scan='+str(int(row['Parent'])),
                                activation=['collision-induced dissociation', {'collision energy':float(variables[(sid,5)])}],
                                isolation_window_args=dict(target=target, lower=width/2, upper=width/2))
                            charge = int(variables.get((sid, 6), 0))
                            if charge: precursor['charge'] = charge
                        writer.write_spectrum(mz, intensity, id='scan='+str(sid), centroided=True,
                            polarity='+' if key['Polarity']==0 else '-', scan_start_time=float(row['Rt'])/60,
                            params=[{'ms level':level}, {'total ion current':float(intensity.sum())}],
                            precursor_information=precursor, encoding={'m/z array':np.float64,'intensity array':np.float64}, compression='zlib')
                        counts[level] += 1
        audit = validate_mzml(temporary)
        if (audit['ms1_spectra'], audit['ms2_spectra']) != (counts[1], counts[2]):
            raise ValueError('Converted spectrum count does not match BAF')
        temporary.replace(output)
        result = dict(identity=identity, source=str(raw), mzml_sha256=digest(output), omitted_empty_spectra=omitted, source_spectra=len(spectra)+len(omitted), **audit)
        save(marker, result)
        print('Converted '+output.name+': '+str(audit), flush=True)
        return result
    finally:
        api.baf2sql_array_close_storage(handle)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw', required=True)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    convert(args.raw, args.out)
