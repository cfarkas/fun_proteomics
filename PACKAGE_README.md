# fun_proteomics

Portable Python and Nextflow launchers for DDA-PASEF and DIA-PASEF proteomics. MIT licensed, maintained by [cfarkas](https://github.com/cfarkas).

```bash
python -m pip install fun_proteomics
fun_proteomics --help
```

Python 3.11 or newer is required. Linux is the supported execution platform. Native DIA spectrum searches use their separately pinned Python 3.12 environment. The Python package includes the workflow source, analysis/reporting code and examples; raw acquisitions, study results and registered search binaries are excluded.

## Continue a folder from anywhere

```bash
fun_proteomics --continue /data/my_study
# Equivalent native Python entry point:
python -m fun_proteomics --continue /data/my_study
```

The folder may be a run root or its `final_results/` directory. The launcher detects its saved state and reusable DIA searches. Completed native deliveries are checked against their output checksums and recorded inputs; valid results are reused. An incomplete or damaged native delivery is rebuilt from its verified inputs or completed searches into `continuations/<timestamp>/`, preserving earlier attempts. Nextflow-managed runs use their recorded session and work cache. Changed raw acquisitions fail reuse when a recorded checksum baseline exists. `--hash-raw` requests a full raw-file hash check; otherwise recorded sizes/timestamps are checked and changed files are hashed. A legacy run without historical raw hashes receives a new baseline, not a claim of retrospective raw identity.

```bash
fun_proteomics --continue /data/old_run --outdir /data/new_delivery \
  --service-matrix tissue=/data/service.tsv --independent-samples --hash-raw
```

Use `--independent-samples` only when each row represents a different biological sample. This preserves old replicate labels and uses unique sample IDs for modeling. The native DIA reporting path runs all unordered condition pairs at protein and unambiguous gene level, using limma-trend on log2 positive intensities without imputation. CSV and TSV matrices, results, service detection gains/losses, quality control and provenance have separate folders.

## Start a raw-data workflow

```bash
fun_proteomics run --config /data/dda.json --outdir /data/results
fun_proteomics run --dia-config /data/dia.json --outdir /data/results_dia
fun_proteomics ui
```

These commands invoke Nextflow. Java and Nextflow are required; registered spectrum-search programs remain separate installations with their own licenses. Use `fun_proteomics doctor` to inspect dependencies, or `--project /path/to/checkout` to reuse an existing checkout and its runtime. Packaged workflow files are verified and extracted to a writable user cache; `--cache-dir` changes that location.

Direct Nextflow remains supported:

```bash
nextflow run "$(fun_proteomics workflow-path)/main.nf" --continue /data/my_study
```

Native matrix analysis is available without Nextflow:

```bash
fun_proteomics quantify --matrix /data/proteins.tsv --samples /data/samples.csv \
  --outdir /data/quantification --differential
```

The installed `quantify` command writes its organized delivery to `OUTDIR/final_results/` and stores its continuation state in `OUTDIR/.fun_proteomics/`. Use the same `--continue OUTDIR` command for matrix analyses.

Differential analysis also requires R with `limma` and `jsonlite`; DDA count modeling requires the corresponding Bioconductor packages. Configure existing installations using `--rscript` and `--r-library` before the subcommand. R and registered search engines are not installed by pip.

```bash
fun_proteomics --rscript /usr/bin/Rscript --r-library /data/R-library \
  --continue /data/my_study
```

Read the full workflow documentation at [github.com/cfarkas/fun_proteomics](https://github.com/cfarkas/fun_proteomics). More detected protein groups do not alone demonstrate better accuracy: service searches may use different databases, grouping and confidence settings.
