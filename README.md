# fun_proteomics

Analyze Bruker DDA-PASEF and DIA-PASEF data, or import a MaxQuant Excel workbook / protein quantification table. Get protein and gene matrices, all-versus-all comparisons, quality control, and detection benchmarks against a service provider. Every matrix is delivered as **CSV and TSV**.

```bash
pip install fun_proteomics
fun_proteomics setup --project /data/my-study
```

**`--project` is your output folder.** `--outdir` is an equivalent setup option. The browser previews the full destination before saving: `PROJECT/results/final_results/matrices/` holds protein and gene CSV/TSV files.

Already installed? `pip install --upgrade fun_proteomics` updates the launcher.
[PyPI package](https://pypi.org/project/fun-proteomics/) · [Watch drag-and-drop setup](https://cfarkas.github.io/fun_proteomics/#videos)


## Start here

| What you want to do | Guide |
| --- | --- |
| Try a small example | [Quick start](docs/quick_start.md) |
| Watch interactive setup | [Walkthrough videos](https://cfarkas.github.io/fun_proteomics/) · [transcripts](docs/videos.md) |
| Import MaxQuant or Bruker service LFQ, with paired donors | [Provider workbooks](docs/provider_workbooks.md) |
| Analyze your own data | [Set up a project](docs/setup.md) |
| Install analysis tools | [Installation](docs/installation.md) |
| Find matrices or understand comparisons | [Outputs](docs/outputs.md) |
| Continue an earlier analysis | [Continue and recover](docs/continuation.md) |
| Resolve an error | [Troubleshooting](docs/troubleshooting.md) |

## Install

Linux and Python 3.11 or newer:

```bash
python3 -m venv proteomics-env
source proteomics-env/bin/activate
python -m pip install fun_proteomics
fun_proteomics system --path /data/my-study
fun_proteomics install --conda
fun_proteomics doctor
```

The Python package installs the launcher and reporting code. `install --conda` prepares the pinned analysis runtime and R packages in a user-writable cache; it does not use sudo. Registered DDA search archives remain separate. If you already have R and the analysis tools, see [reuse an installation](docs/installation.md#reuse-an-installation).

## Set up your analysis

```bash
fun_proteomics setup --project /data/my-study
```

Keep the terminal open. Setup opens a local browser: choose your input, click **Discover samples**, and drag sample cards into condition columns, just as in OncoTracer. Add or rename conditions, move cards between them, or return them to **Unassigned** to exclude them. Each card also has a condition dropdown. Then click **Save configuration and check** and **Run / continue analysis**. Results and progress appear in the same page. Use the Browse buttons to select files on the computer running the analysis.

For terminal questions, including over SSH:

```bash
fun_proteomics setup --terminal --project /data/my-study --run
```

Both interfaces save `my-study/config/run.json`. Browser assignments are also exported as `config/sample_assignments.csv` and `.tsv`, including exclusions. To review first, omit `--run` in terminal setup. Then:

```bash
fun_proteomics check --config /data/my-study/config/run.json
fun_proteomics run --config /data/my-study/config/run.json
```

`check` validates inputs and explains the planned comparisons. `run` starts the saved analysis or continues it. Choose independent biological replicates or paired matrix comparisons with explicit donor IDs. Four conditions produce six comparisons at each feature level.

## Find results and continue

Open **`my-study/START_HERE.html`** for results or **`my-study/MATRICES.html`** for CSV/TSV downloads. Detailed deliverables live in `my-study/results/final_results/`.

```bash
fun_proteomics --continue /data/my-study
```

Completed native results are verified and reused. Incomplete reports are rebuilt from verified inputs; raw workflows resume their saved Nextflow session. Existing runs from earlier releases also accept `--continue /path/to/run`. [Integrity checks and recovery details](docs/continuation.md).

Native Python and direct Nextflow remain supported. See [advanced commands](docs/advanced_reference.md) for raw JSON configurations, the advanced raw editor, phosphoproteomics, DDA counts, and specialist settings.

MIT licensed. Maintained by [cfarkas](https://github.com/cfarkas). The setup/check/run experience follows [OncoTracer](https://github.com/cfarkas/oncotracer).
