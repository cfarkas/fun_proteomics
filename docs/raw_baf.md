# Raw BAF re-search

Choose **Raw BAF · Impact/QTOF** in setup for conventional Bruker DDA acquisitions containing `analysis.baf`. This route converts vendor centroid MS1 and MS2 spectra to indexed mzML, searches with MSFragger/FragPipe, and quantifies MS1 abundance with IonQuant MaxLFQ. It runs through native Python. PASEF continues through Nextflow.

```bash
pip install --upgrade fun_proteomics
fun_proteomics install --conda
fun_proteomics install --baf
fun_proteomics setup --mode baf --project /data/my-baf-study
```

`--project` is the output folder. Select the raw folder, discover acquisitions, and drag each card into its condition. Every acquisition has a unique analysis ID. Give technical injections from the same biological sample the **same biological sample ID**; they share one FragPipe experiment. For donors measured across conditions, select paired comparisons and enter each donor ID. Technical injections never increase the biological replicate count.

Select a human or other appropriate search FASTA containing `GN=` gene annotations, contaminants and matching `rev_` decoys, plus the registered MSFragger, IonQuant and diaTracer archives listed in [installation](installation.md). An existing prepared FASTA can be reused. Registered tools and the Bruker reader are separate dependencies; the repository's MIT license does not replace their licenses.

Supply the service LFQ workbook, CSV or TSV in the same setup screen to benchmark detections automatically. The service table must contain one quantity per biological sample. A sample table can provide `service_sample` aliases; an optional separate service mapping is also supported. Each technical injection of a biological sample must agree on condition, donor and service alias.

A scripted launch uses the same saved configuration:

```bash
fun_proteomics run --config /data/my-baf-study/config/run.json
fun_proteomics --continue /data/my-baf-study
```

For SSH, add `--terminal` to setup. Sample tables use `sample_id,raw_path,condition,biological_sample,donor,service_sample`; `donor` and `service_sample` are optional for independent samples with matching names. Specify `--subject-column donor` for pairing. `--python /path/to/python` selects an existing runtime with pyBaf2Sql and psims; otherwise the runtime installed by `install --baf` is used.

[Watch the paired BAF setup walkthrough](https://cfarkas.github.io/fun_proteomics/#baf-video). Its synthetic files demonstrate setup only; they are never searched.

## Outputs and reuse

- `results/final_results/START_HERE.html`: result dashboard.
- `results/final_results/MATRICES.html`: protein and gene matrices, with CSV/TSV downloads.
- `results/final_results/matrices/*_spectral_counts.csv` and `.tsv`: direct protein, peptide and unambiguous peptide-to-gene spectral counts. These are supporting evidence; all-pairs models use LFQ intensities.
- `results/final_results/differential_expression/`: every unordered condition pair, independently for proteins and available unambiguous genes.
- `results/final_results/benchmark/`: matched service detection summaries, gains/losses, sample mapping and group completeness audit, in CSV and TSV.
- `results/inputs/`: acquisition-to-biological-sample mapping, source provenance and MS1/MS2 conversion counts.
- `results/intermediates/mzML/`: converted spectra, per-acquisition logs, hashes and isolated vendor input copies.
- `results/intermediates/search/run/`: complete FragPipe search and LFQ outputs, logs, workflow and manifest.

Continuation hashes the BAF components, calibration inputs, converter library, mzML, search inputs and required search reports. Successful conversions and complete matching searches are reused. A damaged conversion is regenerated; an incomplete search is preserved under `previous_attempts/` and rerun. Changed source data or settings require a new project. Original acquisition directories are never opened for writing. Empty centroid scans are omitted with their original IDs and MS levels recorded; no peaks are invented. A nonempty MS2 scan referring to an empty MS1 parent is rejected. Conversion fails on unsupported scan modes, missing centroid arrays, missing precursors or corrupt/truncated spectra.

The default closed search uses trypsin, up to two missed cleavages, peptide length at least seven, initial precursor tolerance ±20 ppm and fragment tolerance 0.05 Da. MSFragger calibration and parameter optimization select final tolerances/peak filtering from the data; inspect the delivered fragger.params and pipeline.log for the values actually used. Carbamidomethyl C is fixed; oxidation M and protein N-terminal acetylation are variable. Protein and PSM filtering use 1% FDR. MSBooster is disabled for this non-PASEF preset. IonQuant MBR uses its 1% ion transfer FDR; the exact workflow is delivered with each search. These settings are recorded rather than adjusted to maximize a provider benchmark.

MS1 LFQ includes transferred quantifications when MBR is enabled. It is not a count of independently identified proteins. Provider databases, search settings, protein inference, FDR and transfer rules can differ. See [service benchmarking](benchmarking.md).

Only conventional DDA BAF scan modes are supported. bbCID/DIA BAF and profile-only BAF fail explicitly. MGF exports lack the MS1 data needed for this LFQ route. See the [FragPipe input-format guidance](https://fragpipe.nesvilab.org/docs/tutorial_fragpipe.html), [pyBaf2Sql reader](https://github.com/gtluu/pyBaf2Sql), and [TIMSCONVERT methods](https://doi.org/10.1093/bioinformatics/btac419).
