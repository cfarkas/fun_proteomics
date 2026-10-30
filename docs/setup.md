# Set up your own data

Choose browser setup, terminal questions or explicit command-line options. All three methods save the same project configuration.

## A) Interactive settings

Watch the browser walkthrough: select files and an output folder, assign conditions, run, find matrices, and continue.

https://github.com/user-attachments/assets/bcd925bb-c021-4ccb-92fd-8ac3f4603049

Start in your working directory:

```bash
cd /path/to/dir
```

Invoke interactive setup:

```bash
fun_proteomics setup
```

Start with **Discover samples** at the top. Choose your data type, then use **Choose folder** for raw `.d` data or **Choose file** for a matrix. Selecting the input discovers samples automatically. Assign their conditions, then choose the output folder in **Review & save**. No full paths need to be typed. The video uses the [real quick-start download](quick_start.md); follow the same steps with your own files.

Setup opens in your browser on a local, automatically selected port. The setup server runs independently of the launcher terminal. Repeating the same setup command reopens its existing session. Use `--no-browser` to run a foreground server without opening a browser; in that mode, keep the terminal open. Copy the entire URL, including its session code, when opening another browser. Analyses launched with **Run / continue** keep running if you close the setup server; reopen the same project to view progress.

### 1. Discover and assign samples

Choose one input type:

| Input | What to select |
| --- | --- |
| Quantification matrix | A CSV, TSV or XLSX protein table. Discover sample columns, or import an optional sample table. This also accepts a service table as the main matrix. |
| DIA-PASEF raw data | A folder of `.d` acquisitions, a sample table or scanned sample assignments, and a target-only FASTA. |
| DDA-PASEF raw data | A folder of `.d` acquisitions, sample assignments, a validated target/contaminant/decoy FASTA, and registered MSFragger, IonQuant and diaTracer archives. |
| Raw BAF | Conventional DDA `.d` folders containing `analysis.baf`, a search FASTA, registered archives and biological sample assignments. [BAF guide](raw_baf.md). |
| Existing run folder | A previous run root or its `final_results/` directory. Recorded settings and verified intermediates are reused. |

The guided raw workflow covers global label-free analysis. For phosphoproteomics, immunopeptidomics, DDA+, database preparation or specialist options, use `fun_proteomics ui --advanced`; see the [advanced reference](advanced_reference.md). An existing native DIA/DDA JSON can also be selected under Optional analysis settings.

**Choose folder** and **Choose file** open system dialogs on the computer running the analysis: Zenity on Ubuntu, or KDialog when Zenity is unavailable. This applies to raw-data folders, output folders, existing runs and **Open saved project**. Starting setup from a raw-data folder detects its acquisition type from metadata and discovers its samples. Selecting raw data also discovers samples automatically; cancelling preserves the current selection. Rediscovering the same input keeps your edited sample names and condition assignments. Use the dialog’s new-folder control to create an output folder.

New DIA runs record the configuration, reference files and selected acquisitions. Completed results with matching inputs and valid output checksums reopen without starting another search. Unrelated files beside the acquisitions do not trigger a new run; changed inputs or damaged outputs retain normal Nextflow continuation.

Without a desktop session or either dialog tool, setup falls back to its browser picker. It provides **Working folder**, **Home**, mounted-drive shortcuts, clickable breadcrumbs, **Up one folder**, search and **New folder**. Search filters before the display limit is applied. Direct path entry remains optional.

The top **Discover samples** button opens the picker if no input has been selected, or refreshes discovery for the selected input. An output folder is only required when reviewing and saving. Expand **Optional: import sample groups or choose an Excel sheet** to supply a sample table or worksheet name. A CSV code box below **Sample table** shows the expected columns for the selected data type, including `raw_path` for raw acquisitions. Relative raw paths are resolved from the selected raw-data folder, or from the CSV’s folder when no raw-data folder is selected.

Assign samples to conditions using the interactive board:

1. Detected acquisitions or matrix samples appear as cards in **Unassigned**. Imported sample tables prefill their conditions.
2. Drag cards into **Control**, **Treatment**, or a condition you add. Edit the column heading to rename a condition. Use a card's dropdown as a keyboard/touch alternative.
3. Drag cards between columns to correct assignments, or back to **Unassigned** to exclude them. Removing a condition returns its samples to Unassigned. **Clear assignments** returns every card there.
4. Edit sample names if needed; their original acquisition paths or matrix columns stay linked. The filter helps find samples in larger studies.
5. Review the assigned/excluded counts and the number of unique comparisons before continuing.

Only assigned samples enter matrices, comparisons and service benchmarks. Conditions are never inferred from filenames. For plain matrices, discovery lists numeric candidate columns; review them and leave annotation columns unassigned. Use a sample table for explicit mappings, including all-missing columns. FragPipe suffixes and MaxQuant `LFQ intensity SAMPLE` / `Intensity SAMPLE` prefixes are recognized automatically. Excel cover pages are skipped; use **Excel worksheet** or `--sheet` when a workbook contains more than one eligible quantity sheet. [Provider workbook details](provider_workbooks.md).

An optional **`samples.csv`** needs:

```csv
sample_id,condition
C1,Control
C2,Control
T1,Treatment
T2,Treatment
```

**`proteins.csv`** — illustrative quantities; sample names match the table above:

```csv
Protein.Group,Genes,C1,C2,T1,T2
P04637,TP53,120000,135000,260000,245000
P60709,ACTB,480000,505000,490000,510000
P02768,ALB,90000,85000,,100000
```

Sample IDs must match matrix columns. A protein table needs a recognized identifier such as `Protein.Group`; optional `Genes` and `Protein.Names` are retained. Use one cohort or tissue type per matrix analysis. A separate service mapping can supply `sample_id,service_sample` when provider names differ.

**`service_samples.csv`** — optional service aliases for these sample IDs:

```csv
sample_id,service_sample
C1,Lab_C1
C2,Lab_C2
T1,Lab_T1
T2,Lab_T2
```

**`service.csv`** — illustrative quantities matching the aliases above:

```csv
Protein.Group,Genes,Lab_C1,Lab_C2,Lab_T1,Lab_T2
P04637,TP53,118000,132000,255000,240000
P60709,ACTB,475000,500000,485000,505000
```

Raw sample tables additionally need `raw_path`; paths can be absolute or relative to the selected raw-data folder. Alternatively, **Choose folder** discovers acquisition cards automatically for the same assignment board. No condition is inferred from a filename. Blank acquisitions are excluded from the scan. You can rename sample IDs; DDA identifiers must start with a letter and contain only letters, digits and underscores.

**`samples.csv`** — raw PASEF sample-table version:

```csv
sample_id,condition,replicate,raw_path
C1,Control,1,C1.d
C2,Control,2,C2.d
T1,Treatment,1,T1.d
T2,Treatment,2,T2.d
```

For BAF biological-sample grouping, use the [acquisition table in the complete BAF chapter](raw_baf.md#1-install-and-gather-inputs).

All-versus-all testing is the default. Three conditions produce three unique comparisons; four produce six. Confirm that each sample is from a different patient or independent biological replicate. Original replicate labels are retained in the normalized matrix sample sheet; unique replicate numbers are used for the independent design. For repeated donors in matrix imports or native BAF, select **Paired comparisons** and specify the donor column, or enter donor IDs on the sample cards. The model uses within-donor log2 differences; at least two complete pairs per feature are required. Aggregate technical injections before matrix import. In native BAF, assign technical injections the same biological sample ID; the workflow groups them before statistics. Other complex repeated-measures designs still need a tailored model.

DIA and imported intensity tables use log2 positive quantities with limma-trend, no imputation, and no extra normalization by default. Median centering is optional for suitable unnormalized input. BH remains the default FDR method. Optional analysis settings offer `fdrtool_pvalue` on unchanged limma P values and `fdrtool_normal` as an empirical-null sensitivity analysis; [the FDR guide](fdr_methods.md) explains their distinct P-value columns and diagnostics. DDA uses count models. Native JSON imports retain their explicitly configured comparisons and design; they do not inherit the guided defaults.

Supply the matching FASTA with `GN=` annotations when gene symbols are absent. Ambiguous gene assignments remain visible in the audit and are excluded from gene rollups. A service benchmark reports detection gains and losses; it is separate from differential abundance.

**`targets.fasta`** — gene-annotation format illustration; use your full matching database:

```fasta
>sp|DEMO001|EXAMPLE_A GN=GENE_A
MPEPTIDEK
>sp|DEMO002|EXAMPLE_B GN=GENE_B
MPEPTIDER
```

### 2. Review and save

Use **Next: review & save** to choose where results will be saved. Click **Choose folder**, browse to the desired parent, then use **New folder** for a new project. Select any required FASTA or registered archives with their **Choose file** buttons. Confirm independent samples or paired donors when comparing conditions. Service benchmarks, compute budgets and other optional settings are in **Optional analysis settings**.

**Result locations** shows where matrices and reports will be written. `setup --outdir /data/my-study` remains an alias for `setup --project /data/my-study`. Keep raw acquisitions separate from the project folder. **Save and continue** checks inputs, saves settings and opens the run step. **Open saved project** at the top reopens earlier settings through the same picker.

### 3. Run and find results

**Check inputs** previews sample numbers, groups, comparisons and gene availability without saving. **Save and continue** saves `PROJECT/config/run.json` and the full assignment audit as `config/sample_assignments.csv` and `.tsv`. The audit includes Unassigned samples marked `included=false`; the analysis uses only assigned rows. Draft settings can be corrected before the first launch, with earlier versions retained in `config/backups/`. After launch, use a new project for a changed design. **Run / continue analysis** launches the same native Python or Nextflow code used from the terminal. Missing tools produce an actionable error in the progress log.

**Stop analysis** requests a stop and preserves inputs and outputs. Use Run / continue after resolving the interruption. Changed saved inputs are rejected rather than silently reused; use a new project for a changed analysis.

**`config/sample_assignments.csv`** — selected audit columns; setup also writes TSV:

```csv
sample_id,condition,included
C1,Control,True
C2,Control,True
T1,Treatment,True
Unused,,False
```

**`config/run.json`** — generated excerpt; use setup to write the complete configuration:

```json
{
  "project": "/data/my-study",
  "mode": "matrix",
  "matrix": "/data/proteins.csv",
  "samples": "/data/samples.csv"
}
```

## B) Terminal questions

```bash
fun_proteomics setup --terminal --project /data/my-study --run
```

Answer the input, condition, design and output-folder questions. Omit `--run` to save for review. Terminal setup works over SSH without browser forwarding.

## C) Full command-line options

```bash
fun_proteomics setup --non-interactive --project /data/my-study \
  --mode matrix --matrix /data/proteins.csv --samples /data/samples.csv \
  --independent-samples --normalization none --cpus 8 --memory-gib 32 --run
```

For paired matrix/BAF data, replace `--independent-samples` with `--subject-column donor`. [All options and raw-data commands](command_line.md) · [Real-data quick start](quick_start.md).

Check and run a saved configuration, or continue its results:

```bash
fun_proteomics check --config /data/my-study/config/run.json
fun_proteomics run --config /data/my-study/config/run.json
fun_proteomics --continue /data/my-study
```

`python -m fun_proteomics` accepts the same arguments. From a source checkout, `nextflow run main.nf --interactive --project /data/my-study` opens the same setup interface; see the [advanced reference](advanced_reference.md).

For SSH browser access, launch with `setup --no-browser --port 8765 --project /data/my-study`, then forward the port with `ssh -L 8765:127.0.0.1:8765 user@server`. Open the printed URL on your computer. Browse paths refer to the server.
