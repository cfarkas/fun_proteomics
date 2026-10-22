# Installation

Linux is the supported execution platform. Python 3.11 or newer runs the launcher; native DIA spectrum searches use the separately pinned Python 3.12 environment.

```bash
python3 -m venv proteomics-env
source proteomics-env/bin/activate
python -m pip install fun_proteomics
fun_proteomics system --path /data/my-study
fun_proteomics install --conda
fun_proteomics doctor
```

Keep the environment activated when running the command. `python -m fun_proteomics` is equivalent. The package includes the scientific workflow and examples. `install --conda` downloads the pinned micromamba runtime, Python, Java, Nextflow and R reporting packages into the writable workflow cache. It does not modify your system package manager. Full progress is shown and saved in `.runtime/install.log` inside that workflow cache. Installing source R packages requires a C/C++ toolchain and the system libraries listed in the [advanced installation reference](advanced_reference.md).

`system` reports available CPU, RAM and disk. Matrix analysis needs much less memory and space than a raw spectrum search. Plan raw searches according to acquisition size, sample number and database size; the displayed hardware capacity is not a prediction of runtime.

Run a mode-specific dependency check:

```bash
fun_proteomics doctor --mode matrix
fun_proteomics doctor --mode dia
fun_proteomics doctor --mode dda
# For conventional Bruker DDA BAF:
fun_proteomics install --baf
fun_proteomics doctor --mode baf
```

`doctor --json` gives machine-readable output and returns a nonzero exit code for missing dependencies. Raw input checks also validate selected acquisitions, reference files and registered archives.

## Reuse an installation

With an existing system R, install only the analysis packages:

```bash
fun_proteomics install --r-only
```

Or point to installed tools and a package library (global flags precede the subcommand):

```bash
fun_proteomics --rscript /usr/bin/Rscript --r-library /data/R-library \
  setup --project /data/my-study
```

Use `--workflow-dir /path/to/pasef-proteomics` to reuse a checkout and its `.runtime/`. The older global `--project` spelling remains an alias for this software directory. The `setup --project` option always means the study directory.

The installed workflow is verified and extracted under `~/.cache/fun_proteomics/`. `--cache-dir /data/software-cache` changes this location. A new version has its own verified source cache; use an existing compatible checkout/runtime explicitly if desired.

## Raw search tools

DDA-PASEF and native BAF use the registered versions in [tool_versions.json](../assets/tool_versions.json): MSFragger 4.4.1, IonQuant 1.11.20 and diaTracer 2.2.1. Obtain these archives under their respective licenses and select them during setup. They are not distributed by this MIT package. The complete pinned DIA environment and method are in [DIA method](../modules/diapasef/method.toml) and [dependency lock](../modules/diapasef/requirements-lock.txt).

BAF additionally uses a separately installed Bruker Baf2Sql reader through pinned pyBaf2Sql and psims. `install --baf` creates its own runtime; it does not bundle or relicense the vendor library. See [raw BAF](raw_baf.md) for supported acquisitions, technical-injection grouping and MS1 LFQ.

## Source installation

```bash
git clone https://github.com/cfarkas/fun_proteomics.git
python -m pip install -e ./fun_proteomics
fun_proteomics setup --project /data/my-study
```

Keep the checkout for an editable installation. Direct Nextflow remains available with the same pinned workflow and its [advanced configuration](advanced_reference.md).
