# Execution validation — 22 October 2026

The EV_PPC raw BAF study completed conversion, database search, protein inference, MS1 quantification, protein/gene exports, paired all-versus-all analysis and service benchmarking with fun_proteomics 0.9.0. This page reports aggregate execution evidence. Acquisitions, sample-level quantities, local file paths and service files are not distributed in this repository.

## Full EV_PPC execution

- Instrument/data: Bruker Impact II nanoLC-QTOF, conventional DDA/autoMSMS BAF. This is not PASEF.
- Design: 48 acquisitions, three technical injections per biological sample, 16 biological samples, four donors and four conditions (CTL, 2DG, Oligo, TI). All four donors were retained. Technical injections share a FragPipe experiment; paired comparisons use donor IDs.
- Conversion: 470,977 retained MS1 and 1,386,238 retained MS2 spectra. 286 genuinely empty centroid scans were omitted and audited by original spectrum ID. Every binary array and precursor reference was validated. Original acquisition components were read from isolated copies.
- Direct MS2 protein groups after filtering: 1,285. Positive MS1 LFQ protein groups: 872. Positive unambiguous gene rollups: 856. These are distinct measurements; LFQ can include transferred quantification without direct MS2 in that sample.
- Every condition pair was tested at protein and gene levels: six contrasts per level, 12 result tables in total. Full tables include nonsignificant and untestable rows.
- Final search/report execution completed at 2026-10-22T21:19:33.069766+00:00. The search was launched from [source commit a9a6de9](https://github.com/cfarkas/fun_proteomics/commit/a9a6de9477b7ecf483d0390cc6c377edf0494254).
- Delivery integrity: the complete output SHA-256 inventory passed verification; 79 CSV/TSV pairs contain matching tables. The full-study `--continue` invocation verified and reused all 48 conversions. The completed full-study search passed the same checkpoint predicate used for continuation. Real pilot continuations reused completed searches and rebuilt reports without changing their search outputs. Configuration, database, tool versions, conversion audits, executed search workflow and logs are retained in the local run.

The workflow used FragPipe 24.0, MSFragger 4.4.1, IonQuant 1.11.20, Percolator and ProteinProphet. The native converter used pinned pyBaf2Sql 0.2.2 / Bruker Baf2Sql and psims 1.3.5. The database contained 20,431 human target entries, 49 contaminant targets and 20,480 matching decoys; its SHA-256 was `c8e0f4d450356987746be60ecb2cdd9c26b18e7805fe1a3a080edef72f6bb8b5`.

The closed search used trypsin, up to two missed cleavages, minimum peptide length seven, fixed carbamidomethyl C, variable oxidation M and protein N-terminal acetylation. Initial precursor tolerance was ±20 ppm and fragment tolerance 0.05 Da; MSFragger calibration and parameter optimization selected a final fragment tolerance of 25 ppm and 200 retained peaks per spectrum, minimum intensity ratio 0, no intensity transform, and precursor-peak removal; the precursor search window remained ±20 ppm. The complete optimized settings are recorded in the search log. Protein and PSM identification filtering used 1% FDR. IonQuant used MaxLFQ and 1% ion-transfer FDR for match-between-runs. MSBooster was disabled and regular LC–MS mode was selected.

Paired limma-trend tests used within-donor log2 differences, at least two complete donor pairs per feature, no imputation and no additional median normalization. Gene LFQ is an unambiguous protein-group rollup; supporting gene spectral counts use unambiguous peptide-to-gene assignments. BH correction is separate for each feature level and contrast. Positive logFC means comparison minus reference.

| Feature | Reference | Comparison | Tested | BH FDR < 0.05 |
| --- | --- | --- | ---: | ---: |
| gene | 2DG | CTL | 295 | 0 |
| gene | 2DG | Oligo | 247 | 0 |
| gene | 2DG | TI | 237 | 0 |
| gene | CTL | Oligo | 286 | 0 |
| gene | CTL | TI | 273 | 0 |
| gene | Oligo | TI | 234 | 0 |
| protein | 2DG | CTL | 303 | 0 |
| protein | 2DG | Oligo | 254 | 0 |
| protein | 2DG | TI | 244 | 0 |
| protein | CTL | Oligo | 292 | 0 |
| protein | CTL | TI | 283 | 0 |
| protein | Oligo | TI | 242 | 0 |

## Service benchmarks across available runs

| Study | Search | Matched samples | Pipeline positive groups | Service positive groups | Net change |
| --- | --- | ---: | ---: | ---: | ---: |
| Cell-line study | DDA | 16 | 6,083 | 6,216 | -133 (-2.14%) |
| Cell-line study | DDA+ | 16 | 7,117 | 6,216 | +901 (+14.49%) |
| Independent-patient study | DIA | 49 | 3,406 | 2,974 | +432 (+14.53%) |
| EV_PPC | BAF DDA | 16 | 872 | 445 | +427 (+95.96%) |

These are unions of positive reported groups over exactly matched samples. Zero/NA quantities do not count as detections. The cell-line DDA and DDA+ searches remain separate analyses of the same 16 samples. The 49-patient DIA study retains its independent-sample design. Existing searches were reused for those benchmarks.

For EV_PPC, the service matrix had 445 positive reported groups, including 163 groups reconstructed from truncated headers. Incomplete groups remain in the reported-count comparison but are excluded from exact identity matching. Among complete groups, the pipeline had 872, the service had 282, 149 were shared, 723 were pipeline-only and 133 were service-only. Supporting accession overlap is provided separately and is not an independent-protein count.

The provider described MaxQuant 1.6.17.0, a February 2021 human UniProt database, one missed cleavage, minimum peptide length five, variable N deamidation/M oxidation, and 2.5% protein identification FDR. Its original statistical analysis excluded D4 and imputed missing values. This execution retains D4 and does not impute, so it does not reproduce that statistical design. Differences in database, modifications, protein inference, FDR and transfer rules limit accuracy claims: additional coverage alone does not establish superiority. One cell-line service sample alias was resolved against its original acquisition manifest; the exact mapping is retained in the local benchmark audit.

## Validation and recovery

The [release CI run](https://github.com/cfarkas/fun_proteomics/actions/runs/39846508759) passed all three jobs. It ran 265 source tests with five environment-dependent skips; four skipped DIA tests were exercised in the separate DIA job, while the registered-archive check requires private tool archives. The local source suite and the added scratch-storage and small-study regressions also passed. Additional guided BAF continuation, converted-file staging/integrity, installed-package, real-browser, count/DIA reporting and damaged-output recovery checks also passed. Real BAF tests covered ordinary MS1 intensity for two biological samples, MaxLFQ with technical-injection grouping for three biological samples, the full 16-sample analysis, converted-file corruption recovery and report recovery without repeating completed searches. One- and two-sample studies are labeled as ordinary intensity because IonQuant disables MaxLFQ below three biological experiments. The two-sample pilot produced matching matrices with standard staging and local scratch (relative tolerance 1e-8); scratch caches were removed after successful search, and completed-search reuse still passed. The two raw-data walkthroughs are actual UI recordings using explicitly synthetic placeholders; they claim no spectrum search.

The public [real-data quick start](quick_start.md) also completed: 12 measured ubiquitin-affinity samples, 2,902 protein-group rows, 2,696 unambiguous gene rows and six comparisons at each feature level. CSV/TSV tables matched, delivery checksums passed, and `--continue` verified the completed results without repeating the analysis. The downloadable archive includes the original pinned DEP data and a reproducible base-R exporter; its upstream data license is separate from the MIT pipeline code.

## Run and find results

```bash
pip install --upgrade fun_proteomics
fun_proteomics install --conda
fun_proteomics install --baf
fun_proteomics setup --mode baf --project /data/my-baf-study
fun_proteomics --continue /data/my-baf-study
```

Select the output folder, raw acquisitions, biological sample/donor assignments, matching search FASTA, registered tools and optional service matrix in setup. See the [BAF guide](raw_baf.md), [benchmarking guide](benchmarking.md) and [walkthroughs](videos.md).

Open `PROJECT/results/final_results/RAW_EXECUTION.html` for the complete local execution report, `START_HERE.html` for results, and `MATRICES.html` for CSV/TSV downloads. The matrix directory includes protein/gene LFQ abundance and separately labeled protein/peptide/gene spectral counts. Service summaries, overlap lists and all differential result tables also have matching CSV/TSV versions.
