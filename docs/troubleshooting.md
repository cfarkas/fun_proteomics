# Troubleshooting

| Message or symptom | Action |
| --- | --- |
| `fun_proteomics: command not found` | Activate the Python environment used for installation, or run `python -m fun_proteomics`. |
| Missing R, limma or jsonlite | Run `fun_proteomics install --conda`, or `install --r-only` with an existing R installation. Then run `doctor`. |
| Missing Nextflow or Java | Run `install --conda`; raw acquisitions need these tools. Matrix imports use native Python. |
| Browser says invalid session | Reopen the complete URL printed by the launcher, including the code after `#`. |
| Browser runs on a remote server | Use the SSH forwarding instructions in [setup](setup.md), or add `--terminal`. |
| Sample matched zero or multiple columns | Make `sample_id` match the quantity column. Use a separate `service_samples` mapping with `service_sample` when provider names differ. |
| No unambiguous genes | Supply the matching FASTA with `GN=` annotations. Review excluded and unmapped gene assignments; gene names are never guessed. |
| Independent samples not confirmed | Confirm only if each sample is a different patient or biological replicate. For paired matrix data, select Paired comparisons and supply donor IDs. Otherwise omit comparisons. |
| An input or saved setting changed | Create a new project for the new analysis. Restore the original input only if you intend to continue the original analysis. |
| Output folder is not empty | Use `--continue` for an existing run, or select a new project folder. |
| BAF / Impact QTOF acquisitions detected | Choose Raw BAF in setup and run install --baf. Alternatively import the delivered LFQ workbook. MGF alone cannot regenerate MS1 LFQ. See [BAF guide](raw_baf.md). |
| Multiple protein quantity worksheets | Select the worksheet by name in setup or with `--sheet`. |
| A raw acquisition has the wrong frame type | Check the acquisition choice and the raw files. DDA requires frame types 0/8; DIA requires 0/9. |
| Registered tool archive fails validation | Supply the pinned registered archive; see [installation](installation.md). |
| Missing Nextflow session/cache | Restore the original work/cache or run a fresh project. Incomplete searches are not reused as completed results. |
| A report disappeared or was edited | `--continue` checks integrity and rebuilds native reporting into a new generation when its inputs remain valid. |

Find progress and result locations:

```bash
fun_proteomics status --project /data/my-study
fun_proteomics doctor --json
```

`check --config PROJECT/config/run.json` checks input structure and analysis design. `doctor` separately checks installed dependencies. Raw workflow validation and the actual search perform additional engine-specific checks. A successful input preview does not mean a spectrum search has already completed.

## File-format checks

For a sample-name mismatch, these two files must agree:

**`samples.csv`**:

```csv
sample_id,condition
C1,Control
C2,Control
```

**`proteins.csv`**:

```csv
Protein.Group,Genes,C1,C2
P04637,TP53,120000,135000
```

These two rows illustrate matching only; differential testing needs replicated conditions. See [all input previews](input_files.md) for service aliases, paired donors, FASTA and saved configurations.

**`config/run.json`** — excerpt of the setup-generated file checked above:

```json
{
  "project": "/data/my-study",
  "mode": "matrix",
  "matrix": "/data/proteins.csv",
  "samples": "/data/samples.csv"
}
```
