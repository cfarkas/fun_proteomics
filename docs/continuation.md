# Continue and recover

From any working directory:

```bash
fun_proteomics --continue /data/my-study
```

You can also reopen the browser with `setup --project /data/my-study` and select **Run / continue analysis**, or repeat `run --config /data/my-study/config/run.json`.

## What is reused

For guided projects, the launcher checks saved settings and recorded input hashes before entering the existing continuation engine. Native deliveries are verified against `OUTPUT_SHA256SUMS` and their input provenance. Valid completed outputs are reused without recomputing the analysis. Missing or altered native reports are rebuilt into `continuations/<timestamp>/`, preserving previous attempts. The report and matrix shortcuts are updated after successful recovery.

Native BAF checks source components, calibration inputs and converted mzML hashes before reusing spectra. Complete search checkpoints verify settings, database, tools and required reports. Damaged conversions are regenerated; incomplete searches and invalid reports are preserved under `previous_attempts/` before rebuilding. Source-data or analysis-setting changes require a new project. See [raw BAF](raw_baf.md).

Nextflow projects reuse the recorded session and work cache. Matching tasks are cached; unfinished tasks run. The work directory and controller cache are needed for this form of reuse. A damaged or missing cache requires a fresh raw run; it is not treated as a completed spectrum search.

For completed DIA searches, reuse verifies task completion, report checksums, configuration, sample assignments and readable Parquet content. Raw acquisition metadata are checked; changed components with recorded hashes are hashed again. Add `--hash-raw` to request a full raw-file hash check. A newly created raw hash baseline cannot establish historical identity with a run that never recorded one.

## Earlier versions

```bash
fun_proteomics --continue /data/earlier-run
fun_proteomics --continue /data/earlier-run/final_results
```

Setup's **Existing run folder** mode creates a small project that points to an earlier run. It reuses that run's recorded design. It does not recopy acquisitions or overwrite your original reporting attempts.

## Changed inputs and moved folders

For changed samples, quantities or design, create a new project. A saved project cannot silently switch to new inputs. The original files and paths must remain accessible; the launcher can run from any working directory, but it does not automatically relocate a study between filesystems. It also refuses a second active launch for the same project.

Stopping preserves the project. Resolve the error shown in the log and run the same continuation command. For missing tools, see [installation](installation.md); for missing cache or changed-input errors, see [troubleshooting](troubleshooting.md).
