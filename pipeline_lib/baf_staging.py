"""Optional fast storage for disposable BAF search caches."""
from contextlib import contextmanager
from pathlib import Path
import shutil
import tempfile
from pipeline import read, save


@contextmanager
def search_staging(searchdir, converted, acquisitions, scratch_base=None):
    """Keep completed search tables in the project while relocating large caches."""
    if not scratch_base:
        yield False
        return
    searchdir = Path(searchdir)
    base = Path(scratch_base).expanduser().resolve()
    if any(base.is_relative_to(Path(source).resolve()) for source in acquisitions):
        raise ValueError('Search scratch must be outside every source acquisition')
    base.mkdir(parents=True, exist_ok=True)
    # mzBIN and IonQuant indexes need space in addition to the linked input data.
    estimate = int(sum(Path(p).stat().st_size for p in converted) * 1.4) + 2 * 1024**3
    if shutil.disk_usage(base).free < estimate:
        raise ValueError(f'Insufficient search scratch space; estimated requirement {estimate/1024**3:.1f} GiB')
    scratch = Path(tempfile.mkdtemp(prefix='fun-proteomics-baf-', dir=base))
    raw = scratch/'raw'; raw.mkdir()
    staging = searchdir/'raw'
    if staging.exists() or staging.is_symlink():
        shutil.rmtree(scratch)
        raise ValueError('Search staging directory already exists')
    staging.symlink_to(raw, target_is_directory=True)
    record = dict(scratch=str(scratch), estimated_bytes=estimate, status='running',
        persistent_converted=[str(Path(p).resolve()) for p in converted],
        purpose='Temporary mzBIN and IonQuant indexes; completed search tables remain in the project')
    metadata = searchdir/'run/scratch_staging.json'; save(metadata, record)
    print('Using local search caches: '+str(scratch), flush=True)
    try:
        yield True
    except BaseException:
        record['status']='interrupted; caches retained for inspection'
        save(metadata, record)
        raise
    else:
        # Conversion files and all tables needed for verified continuation persist.
        markers = {p.name:read(p) for p in raw.glob('*.copy_complete.json')}
        staging.unlink(); staging.mkdir()
        for name, marker in markers.items():
            save(staging/name, marker)
            source=Path(marker['source'])
            (staging/source.name).symlink_to(source)
        (staging/'README.txt').write_text('Converted spectra are linked from ../../mzML. Temporary mzBIN and IonQuant caches were removed after the search completed. The verified search tables, logs, workflow and checkpoint remain in ../run. --continue reuses the completed search without these temporary caches.\n')
        shutil.rmtree(scratch)
        record['status']='completed; disposable caches removed'
        save(metadata, record)
