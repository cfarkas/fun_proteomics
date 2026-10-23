# Raw BAF: one complete workflow

Use this chapter for conventional Bruker Impact/QTOF DDA acquisitions containing `analysis.baf`. The native Python workflow converts MS1/MS2 spectra, runs FragPipe/MSFragger, quantifies with IonQuant, and exports protein/gene matrices and all condition pairs. For PASEF, use the [DDA/DIA commands](command_line.md#raw-dia-pasef).

## 1. Install and gather inputs

```bash
pip install --upgrade fun_proteomics
fun_proteomics install --conda
fun_proteomics install --baf
```

You need your complete `.d` folders, a matching search FASTA with `GN=` annotations, contaminants and `rev_` decoys, and registered MSFragger 4.4.1, IonQuant 1.11.20 and diaTracer 2.2.1 ZIP archives. Select these archives during setup; they retain their own licenses.

Keep each acquisition folder intact; the files inside it are binary instrument data:

```text
raw/
  D1_CTL_R1.d/
    analysis.baf
    ...other original acquisition components...
  D1_CTL_R2.d/
  D1_Treatment_R1.d/
```

**`acquisitions.csv`** — optional sample assignments; illustrative rows:

```csv
sample_id,raw_path,condition,biological_sample,donor
D1_CTL_R1,D1_CTL_R1.d,Control,D1_CTL,D1
D1_CTL_R2,D1_CTL_R2.d,Control,D1_CTL,D1
D1_Treatment_R1,D1_Treatment_R1.d,Treatment,D1_Treatment,D1
D2_CTL_R1,D2_CTL_R1.d,Control,D2_CTL,D2
D2_Treatment_R1,D2_Treatment_R1.d,Treatment,D2_Treatment,D2
```

`raw_path` is relative to the raw folder or absolute. Technical injections share `biological_sample`. Paired samples share `donor` across conditions. Technical injections do not increase the biological replicate count. Omit `donor` for an independent design.

<details>
<summary>Search FASTA format</summary>

**`search.fasta`** — format illustration only:

```fasta
>sp|DEMO001|EXAMPLE_A GN=GENE_A
MPEPTIDEK
>rev_sp|DEMO001|EXAMPLE_A GN=GENE_A
KEDITPEPM
>contam_DEMO
MPEPTIDER
>rev_contam_DEMO
REDITPEPM
```

Use a complete, validated organism database with matching decoys. These tiny demonstration sequences are not an analysis database. The advanced database setup can prepare a suitable search FASTA.

</details>

## 2. Configure in the browser

The embedded recording shows condition assignment, technical-injection grouping and donor pairing. Its raw files are synthetic placeholders, so it demonstrates setup only.

https://github.com/user-attachments/assets/a422425b-a9ad-4b8b-929f-b2bb9255a20e

Open setup with `fun_proteomics setup --mode baf`, then:

1. Choose the **output/project folder** and raw-data folder.
2. Click **Discover samples** and drag acquisitions into conditions, or import the sample table above.
3. Set biological sample IDs. Choose independent replicates or paired donors, and select the FASTA and registered tools.
4. Optionally select service quantification, then **Save configuration and check** and **Run / continue analysis**.

The optional service table needs one column per **biological sample**, after the service has combined any technical injections:

**`service.csv`** — illustrative biological-sample LFQ quantities:

```csv
Protein.Group,Genes,D1_CTL,D1_Treatment,D2_CTL,D2_Treatment
P04637,TP53,120000,260000,135000,245000
P60709,ACTB,480000,490000,505000,510000
```

For Excel sheets or different service sample names, use the [service benchmark guide](benchmarking.md), which includes workbook and mapping previews.

## 3. Terminal or scripted launch

For terminal questions, including SSH:

```bash
fun_proteomics setup --terminal --mode baf --project /data/my-baf-study --run
```

<details>
<summary>Full command-line alternative</summary>

For a paired launch using the files above:

```bash
fun_proteomics setup --non-interactive --mode baf --project /data/my-baf-study \
  --raw-root /data/raw --samples /data/acquisitions.csv --fasta /data/search.fasta \
  --msfragger-zip /data/tools/MSFragger-4.4.1.zip \
  --ionquant-zip /data/tools/IonQuant-1.11.20.zip \
  --diatracer-zip /data/tools/diatracer-2.2.1.zip \
  --subject-column donor --service-matrix /data/service.csv \
  --cpus 16 --memory-gib 64 --run
```

`--project` is the output folder. Use `--independent-samples` instead of `--subject-column donor` for an independent design. Omit `--service-matrix` when no service table is available.

</details>

## 4. Open results and continue

| Open inside your project | What you get |
| --- | --- |
| `START_HERE.html` | All comparisons, quality control and service benchmark. |
| `MATRICES.html` | Protein/gene LFQ and separate protein/peptide/gene spectral counts, each in CSV and TSV. |
| `results/final_results/RAW_EXECUTION.html` | Search settings, conversion totals, detection counts and execution evidence. |

Four conditions produce six protein and six gene comparisons when gene annotations are available. The models use LFQ abundance. See [output table previews](outputs.md#table-previews).

```bash
fun_proteomics --continue /data/my-baf-study
```

Keep the entire project. Continuation verifies source components, conversions and complete search reports; valid work is reused. Damaged conversions are regenerated, and incomplete attempts are preserved before rebuilding. Changed inputs or analysis settings require a new project. Original acquisition folders are read from isolated copies.

<details>
<summary>Search defaults, supported data and optional fast scratch</summary>

The search uses trypsin, up to two missed cleavages, peptide length at least seven, fixed carbamidomethyl C, variable oxidation M and protein N-terminal acetylation. Initial precursor/fragment tolerances are ±20 ppm / 0.05 Da; MSFragger calibrates and optimizes them. Protein and PSM filtering use 1% FDR. IonQuant MBR uses 1% ion-transfer FDR; MSBooster is disabled. Submitted settings and optimized values are retained in the search provenance and log. LFQ can include match-between-runs transfers; direct MS2 spectral counts are reported separately.

MaxLFQ needs at least three biological samples; smaller studies export ordinary MS1 intensity and label it accordingly. Technical injections do not increase that count; see the [IonQuant changelog](https://github.com/Nesvilab/IonQuant/blob/master/CHANGELOG.md). Only conventional centroid DDA BAF is supported; bbCID/DIA BAF, profile-only data and MGF-only inputs cannot supply this MS1 LFQ workflow. Empty scans are audited, and corrupt or unsupported spectra fail validation.

For slow shared storage, keep disposable caches on a fast disk:

```bash
FUN_PROTEOMICS_SCRATCH=/fast/scratch fun_proteomics --continue /data/my-baf-study
```

Converted spectra and complete search outputs stay in the project. Successful runs remove disposable scratch caches. Use `--python /path/to/python` only to select an existing converter runtime instead of the one installed above.

</details>
