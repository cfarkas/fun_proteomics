# MaxQuant and Bruker service LFQ

Choose your output folder explicitly, then select the provider workbook:

```bash
pip install --upgrade fun_proteomics
fun_proteomics setup --outdir /data/ev-analysis \
  --mode matrix --matrix /data/service-report.xlsx
```

The browser finds the protein quantity worksheet, reads MaxQuant `LFQ intensity SAMPLE` columns, and shows the samples on the condition board. Drag them into your conditions. Cover pages and statistical-result sheets are skipped. Use **Excel worksheet** or `--sheet 'Identified proteins'` to select a specific sheet. CSV/TSV `proteinGroups.txt` reports with `Protein IDs` and `Gene names` also work.

`--outdir` and `--project` mean the same thing in **setup**: the complete study folder. Configuration, logs and results are contained there. Open `ev-analysis/START_HERE.html` or `ev-analysis/MATRICES.html`; protein and gene CSV/TSV files are in `ev-analysis/results/final_results/matrices/`.

## Repeated donors and technical injections

Use one biological quantity per donor and condition. If a provider has already combined triplicate LC-MS injections into one LFQ column, use that column once. Three injections from one donor are not three independent donors. Prepare an optional sample table:

```csv
sample_id,condition,donor
D1_CTL,Control,D1
D1_Treatment,Treatment,D1
D2_CTL,Control,D2
D2_Treatment,Treatment,D2
```

In the browser select **Paired comparisons**, specify `donor` as the subject column, and review the IDs on the sample cards. For automation:

```bash
fun_proteomics setup --non-interactive --outdir /data/ev-analysis \
  --mode matrix --matrix /data/service-report.xlsx \
  --sheet 'Identified proteins' --samples /data/samples.csv \
  --subject-column donor --run
fun_proteomics --continue /data/ev-analysis
```

Every unordered condition pair is tested at protein and gene levels: four conditions give six protein comparisons and six gene comparisons when gene annotations are available. Each paired contrast fits an intercept to within-donor log2 differences with limma-trend. This removes the donor effect, as in a paired linear model. The variance trend uses mean abundance across complete pairs. Each feature needs at least two observed donor pairs; missing values remain missing. BH FDR is separate for each contrast and feature level. No provider-imputed values or precomputed p-values are imported as raw intensities. See the [limma guide, paired samples](https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf).

For independent studies, use `--independent-samples` instead of `--subject-column`. Raw PASEF searches currently retain the independent design; import aggregated quantities for paired downstream analysis.

## Import audit and genes

`quality_control/import_audit.json` records the worksheet and header row, provider flags and any removed empty duplicate rows. MaxQuant reverse, contaminant and site-only rows are excluded. Multiple nonempty rows for the same canonical group are rejected; they are never silently summed. An entirely zero/missing duplicate row in a FASTA-header service report can be removed without changing measured quantities, and its source row is recorded.

Accession and `GN=` fields are read from UniProt FASTA headers when explicit identifier/gene columns are absent. The original headers remain in the protein matrix. Clipped or incomplete protein-group headers are marked `partial_header`; reported genes remain visible, but those groups are excluded from gene sums and exact-group detection benchmarks. A FASTA can supplement annotation, but cannot reconstruct group members missing from a clipped report. Request the complete `proteinGroups.txt` export when complete gene/group coverage is required.

LFQ values must be linear, finite and nonnegative. Log-transformed or imputed result tables need their original pre-transformation quantities; the importer does not undo an unknown transformation. Zero indicates no positive quantity for detection and becomes missing in the log2 statistical model.

## Impact / QTOF BAF acquisitions

A `.d` extension does not identify the acquisition type. `analysis.baf` files from Bruker Impact/QTOF instruments differ from timsTOF `analysis.tdf` PASEF data. Setup detects this distinction and directs BAF studies to the supported provider-matrix route. **Raw BAF re-search is not implemented.** MGF exports can support identification searches, but lack the MS1 data needed for LFQ quantification; see the [FragPipe input guide](https://fragpipe.nesvilab.org/docs/tutorial_fragpipe).

This route analyzes delivered quantities and does not claim new identifications or better detection than the provider. Detection gains/losses require a separate, compatible quantification matrix with explicit sample matching.
