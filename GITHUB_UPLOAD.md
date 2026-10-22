# GitHub preparation

The public repository is **cfarkas/fun_proteomics**. The source package supports native
Bruker DDA-PASEF through FragPipe and DIA-PASEF through DIA-NN, using the same
`main.nf` entry point. `fun_proteomics setup --project /data/my-study` opens the
drag-and-drop sample setup. The direct Nextflow equivalent is
`nextflow run main.nf --interactive --project /data/my-study`. Direct DIA runs use `nextflow run main.nf --dia_config study.json`;
the browser generates that native DIA configuration while preserving the editable study.

## Validate and package

From the repository directory, after installing the runtime described in README:

```bash
.runtime/env/bin/python scripts/validate_repo.py
bash tests/run_all.sh
bash tests/run_smoke.sh
python3.12 tests/run_diapasef_smoke.py
python3 scripts/package_source.py
```

Packaging creates `../pasef-proteomics-source.zip`, its `.sha256` checksum and `SOURCE_SHA256SUMS`. It
uses an explicit source allowlist and excludes symlinks, local configurations,
results, runtime environments, input data and installed search engines. Only
the source-controlled DIA module is included from `modules/`.

The GitHub Actions workflow runs source checks, unit tests and DDA/DIA synthetic
execution/resume checks. Synthetic checks produce no scientific measurements.
Real searches require user-supplied acquisitions, FASTA and separately obtained
search engines.

## Upload the source repository

Use the public GitHub repository `cfarkas/fun_proteomics`. If Git has already been initialized and a source commit prepared,
start with the remote command below. Otherwise initialize a source-only extracted
archive with `git init -b main`, review `git status --short`, add the source files
and commit them.

For a clean source-only checkout:

```bash
git remote add origin git@github.com:cfarkas/fun_proteomics.git
git push -u origin main
```

No remote repository is created and nothing is uploaded by the local preparation
or packaging commands. The archive also supports GitHub's file-upload interface
after extraction.

## What stays local

Keep native acquisitions, study/sample identifiers, FASTA databases, reference
matrices, downloaded or registered tools, licenses/registration credentials,
local configuration JSON, browser drafts and UI state, generated `ui_configs/dia/`
inputs, logs, cached work and results outside Git. Use the
templates under `examples/` for public configuration examples.

Original source and documentation use the [MIT license](LICENSE). Separately installed tools retain their own licenses; see [LICENSE_NOTES.md](LICENSE_NOTES.md).

## Clean installation location

Extract the source ZIP into the parent directory where you want the repository.
For example, `unzip pasef-proteomics-source.zip -d /path/to/projects` creates
`/path/to/projects/pasef-proteomics`. Use a new destination so existing files are preserved.
The package is small because it excludes raw data, outputs, work caches and runtimes.
Choose an output directory outside the source repository for each real study.

Version 0.4.0 adds a top-level matrix catalogue, annotated DIA protein/gene intensities,
all-pairs intensity comparisons, post-search count reanalysis and service detection
benchmarks. The separate DIA CI job exercises actual reporting/statistical code on
synthetic precursor tables; spectrum searches in CI remain explicitly stubbed.

## Python distribution

The public package name is `fun_proteomics` (PyPI normalizes it to `fun-proteomics`). The repository is `cfarkas/fun_proteomics`; the license is MIT. Build and inspect both distributions before publishing:

```bash
python -m build
python -m twine check dist/*
```

The wheel contains a checked workflow source bundle. Both source distribution and wheel exclude local acquisitions, results, tools and credentials. Tests install the wheel outside the checkout and exercise native Python reporting and folder continuation. Publishing to PyPI requires an authorized account or a configured trusted publisher. A release workflow is provided for the `pypi` GitHub environment.
