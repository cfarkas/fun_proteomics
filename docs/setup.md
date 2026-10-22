# Set up your own data

```bash
fun_proteomics setup --project /data/my-study
```

Setup opens in your browser on a local, automatically selected port. Keep the launcher terminal open while configuring. Copy the entire URL, including its session code, when opening another browser. Analyses launched with **Run / continue** keep running if you close the setup server; reopen the same project to view progress.

## 1. Choose data

Choose one input type:

| Input | What to select |
| --- | --- |
| Quantification matrix | A CSV, TSV or XLSX protein table. Discover sample columns, or import an optional sample table. This also accepts a service table as the main matrix. |
| DIA-PASEF raw data | A folder of `.d` acquisitions, a sample table or scanned sample assignments, and a target-only FASTA. |
| DDA-PASEF raw data | A folder of `.d` acquisitions, sample assignments, a validated target/contaminant/decoy FASTA, and registered MSFragger, IonQuant and diaTracer archives. |
| Existing run folder | A previous run root or its `final_results/` directory. Recorded settings and verified intermediates are reused. |

The guided raw workflow covers global label-free analysis. For phosphoproteomics, immunopeptidomics, DDA+, database preparation or specialist options, use `fun_proteomics ui --advanced`; see the [advanced reference](advanced_reference.md). An existing native DIA/DDA JSON can also be selected under Advanced settings.

The project folder holds `config/`, `results/`, logs and launch state. Keep raw acquisitions in a separate directory tree. Save new analyses to a new or empty project. **Load saved project** reopens earlier settings.

## 2. Review samples and comparisons

Click **Discover samples** after selecting a matrix or raw-data folder. The assignment board follows OncoTracer's current interaction:

1. Detected acquisitions or matrix samples appear as cards in **Unassigned**. Imported sample tables prefill their conditions.
2. Drag cards into **Control**, **Treatment**, or a condition you add. Edit the column heading to rename a condition. Use a card's dropdown as a keyboard/touch alternative.
3. Drag cards between columns to correct assignments, or back to **Unassigned** to exclude them. Removing a condition returns its samples to Unassigned. **Clear assignments** returns every card there.
4. Edit sample names if needed; their original acquisition paths or matrix columns stay linked. The filter helps find samples in larger studies.
5. Review the assigned/excluded counts and the number of unique comparisons before continuing.

Only assigned samples enter matrices, comparisons and service benchmarks. Conditions are never inferred from filenames. For plain matrices, discovery lists numeric candidate columns; review them and leave annotation columns unassigned. Use a sample table for explicit mappings, including all-missing columns. FragPipe quantity suffixes are recognized automatically.

An optional matrix sample table needs:

```csv
sample_id,condition
C1,Control
C2,Control
T1,Treatment
T2,Treatment
```

Sample IDs must match matrix columns. A protein table needs a recognized identifier such as `Protein.Group`; optional `Genes` and `Protein.Names` are retained. Use one cohort or tissue type per matrix analysis. A separate service mapping can supply `sample_id,service_sample` when provider names differ.

Raw sample tables additionally need `raw_path`; paths can be absolute or relative to the selected raw-data folder. Alternatively, **Scan .d folders** discovers acquisition cards for the same drag-and-drop board. No condition is inferred from a filename. Blank acquisitions are excluded from the scan. You can rename sample IDs; DDA identifiers must start with a letter and contain only letters, digits and underscores.

All-versus-all testing is the default. Three conditions produce three unique comparisons; four produce six. Confirm that each sample is from a different patient or independent biological replicate. Original replicate labels are retained in the normalized matrix sample sheet; unique replicate numbers are used for the independent design. Repeated measures require a tailored statistical model: turn comparisons off if you only need matrices and QC.

DIA and imported intensity tables use log2 positive quantities with limma-trend, no imputation, and no extra normalization by default. Median centering is optional for suitable unnormalized input. DDA uses count models. Native JSON imports retain their explicitly configured comparisons and design; they do not inherit the guided defaults.

Supply the matching FASTA with `GN=` annotations when gene symbols are absent. Ambiguous gene assignments remain visible in the audit and are excluded from gene rollups. A service benchmark reports detection gains and losses; it is separate from differential abundance.

## 3. Save and run

**Check inputs** previews sample numbers, groups, comparisons and gene availability without saving. **Save configuration and check** saves `PROJECT/config/run.json` and the full assignment audit as `config/sample_assignments.csv` and `.tsv`. The audit includes Unassigned samples marked `included=false`; the analysis uses only assigned rows. Draft settings can be corrected before the first launch, with earlier versions retained in `config/backups/`. After launch, use a new project for a changed design. **Run / continue analysis** launches the same native Python or Nextflow code used from the terminal. Missing tools produce an actionable error in the progress log.

**Stop analysis** requests a stop and preserves inputs and outputs. Use Run / continue after resolving the interruption. Changed saved inputs are rejected rather than silently reused; use a new project for a changed analysis.

## Terminal and automation

```bash
fun_proteomics setup --terminal --project /data/my-study --run
fun_proteomics check --config /data/my-study/config/run.json
fun_proteomics run --config /data/my-study/config/run.json
fun_proteomics --continue /data/my-study
```

`fun_proteomics --interactive --project /data/my-study` is a browser shortcut. `python -m fun_proteomics` accepts the same arguments. `--non-interactive` accepts setup answers as flags; see [the complete example](quick_start.md).

From a source checkout with Python dependencies installed:

```bash
nextflow run main.nf --interactive --project /data/my-study
nextflow run main.nf --continue /data/my-study
```

For SSH browser access, launch with `setup --no-browser --port 8765 --project /data/my-study`, then forward that local port with `ssh -L 8765:127.0.0.1:8765 user@server`. Open the printed URL on your computer. Browse paths refer to the server. `--terminal` needs no forwarding.
