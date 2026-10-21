"""Reuse complete FragPipe searches only when acquisition and search inputs match."""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace

from pipeline import digest, read, save, search_properties, selected_samples

# FragPipe rewrites these installation/working paths when saving its workflow.
# Actual executable versions and archive checksums are checked separately.
IGNORED_PROPERTIES = {'database.db-path', 'workflow.threads', 'workflow.ram', 'workflow.description',
                      'workflow.saved-with-ver', 'workdir', 'fragpipe-config.tools-folder',
                      'fragpipe-config.bin-diann', 'fragpipe-config.bin-python'}


def normalized_properties(properties):
    return {k: v for k, v in properties.items() if k not in IGNORED_PROPERTIES}


def identity(cfg, cohort, mode, properties, database_sha256, versions):
    samples = selected_samples(cfg, cohort)
    audit = {row['sample_id']: row for row in cfg.get('raw_audit', [])}
    if any(s['sample_id'] not in audit for s in samples):
        raise ValueError('Cannot reuse a search without the original acquisition audit.')
    return {'version': 1, 'cohort': cohort, 'mode': mode,
            'analysis': cfg.get('proteomics', {}).get('analysis', 'global-label-free'),
            'samples': [{'sample_id': s['sample_id'], 'raw_path': s['raw_path'],
                         'audit': audit[s['sample_id']]} for s in samples],
            'database_sha256': database_sha256, 'properties': normalized_properties(properties),
            'tools': {k: v for k, v in versions.items() if k != 'runtime'}}


def required_outputs(root, samples):
    paths = [root/'combined_protein.tsv', root/'combined_peptide.tsv']
    paths.extend(root/s['sample_id']/'psm.tsv' for s in samples)
    if not (root/'SEARCH_COMPLETE').is_file() or (root/'STUB_NOT_SCIENTIFIC_RESULTS').exists():
        raise ValueError('Search is incomplete or synthetic.')
    if any(not p.is_file() or p.stat().st_size == 0 for p in paths):
        raise ValueError('Completed search is missing required protein, peptide or sample PSM tables.')
    return paths


def write_checkpoint(root, cfg, cohort, mode, properties, database_sha256, versions):
    root = Path(root)
    signature = identity(cfg, cohort, mode, properties, database_sha256, versions)
    paths = required_outputs(root, selected_samples(cfg, cohort))
    save(root/'search_checkpoint.json', {'identity': signature,
        'outputs': {str(p.relative_to(root)): {'size': p.stat().st_size, 'sha256': digest(p)} for p in paths}})


def candidate_matches(root, wanted):
    root = Path(root).resolve()
    required = required_outputs(root, wanted['samples'])
    checkpoint = root/'search_checkpoint.json'
    if checkpoint.is_file():
        record = read(checkpoint)
        if record['identity'] != wanted:
            return False
        outputs = record['outputs']
        return all(str(p.relative_to(root)) in outputs and
                   p.stat().st_size == outputs[str(p.relative_to(root))]['size'] and
                   digest(p) == outputs[str(p.relative_to(root))]['sha256'] for p in required)
    # Previous releases kept complete provenance and the validated plan in work/.
    # They are sufficient to migrate a successful search without searching again.
    provenance = read(root/'run_provenance.json')
    if provenance.get('stub') or provenance.get('dry_run'):
        return False
    previous = read(root.parent/'plan.json')
    properties = dict(line.split('=', 1) for line in (root/'fragpipe.workflow').read_text().splitlines()
                      if line and not line.startswith('#') and '=' in line)
    observed = identity(previous, provenance['cohort'], provenance['mode'], properties,
                        provenance['database_sha256'], read(root/'software_versions.json'))
    return observed == wanted


def find_completed(cfg, cohort, mode, properties, database_sha256, versions, workdir):
    wanted = identity(cfg, cohort, mode, properties, database_sha256, versions)
    slug = 'ddaplus' if mode == 'DDA+' else 'dda'
    published = Path(cfg['outdir'])/'analysis_details'/cohort/slug/'run'
    candidates = [published]
    work = Path(workdir)
    if work.is_dir():
        candidates += sorted((p.parent for p in work.glob('*/*/run/run_provenance.json')),
                             key=lambda p: (p/'SEARCH_COMPLETE').stat().st_mtime_ns if (p/'SEARCH_COMPLETE').exists() else 0,
                             reverse=True)
    for root in candidates:
        try:
            if candidate_matches(root, wanted):
                # Reporting needs this metadata; legacy published copies lack it.
                if not (root/'run_provenance.json').is_file():
                    continue
                return root.resolve()
        except (OSError, ValueError, KeyError, TypeError):
            continue
    return None


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('plan', 'cohort', 'mode', 'tools', 'database', 'template', 'workdir'):
        parser.add_argument('--'+name, required=True)
    parser.add_argument('--reuse', action='store_true')
    args = parser.parse_args()
    cfg = read(args.plan)
    properties = search_properties(SimpleNamespace(**vars(args), cpus=1, memory=4), cfg)
    found = find_completed(cfg, args.cohort, args.mode, properties, digest(args.database),
                           read(Path(args.tools)/'software_versions.json'), args.workdir) if args.reuse else None
    save('search_reuse.json', {'run': str(found) if found else None})
    print(f'[CONTINUE] {args.cohort} {args.mode}: '+('reuse completed search and quantification: '+str(found)
          if found else 'search required; no verified matching completed search' if args.reuse else 'new search'))
