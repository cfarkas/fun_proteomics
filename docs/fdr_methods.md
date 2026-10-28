# FDR methods for intensity comparisons

DIA-PASEF and imported quantification matrices support three explicit methods. **BH is the default.** The log2 transformation, missing-value handling, normalization, observation filter and limma model are the same for all three choices. Each contrast and feature level is a separate testing family.

| Setting | Input | Meaning of `FDR` | Meaning of `FDR_PValue` |
| --- | --- | --- | --- |
| `BH` | Original limma P values | Benjamini–Hochberg adjusted P | Original limma P |
| `fdrtool_pvalue` | Original limma P values | fdrtool tail-area q | Original limma P |
| `fdrtool_normal` | Moderated t statistics | fdrtool tail-area q under an estimated normal null | Recalibrated normal-null P |

`PValue` **always** contains the original limma P value. `BH_FDR` always contains its BH adjustment. `FDR_method` identifies the selected method, and `significant` applies the configured threshold to `FDR`. `moderated_t` and `df_total` make the original test auditable. `local_fdr` is a density-based estimate, distinct from the tail-area q used for significance.

The P-value-input method estimates the null proportion and P-value distribution while retaining limma's P values. The empirical-normal method instead estimates the scale of a zero-centred normal null from the t statistics, using the `fndr` cutoff. Its recalibrated P value is `2 * pnorm(-abs(t) / estimated_sd)`. Consequently, its q value can be much smaller than the **original** limma P value. This does not mean that BH reduced that original P value: the null calibration changed.

Treat empirical-normal results as a sensitivity analysis. Moderated t statistics have their own t-based null and may have different degrees of freedom across features. Review the fitted null and warnings; do not choose a method because it returns more significant features. A smaller q value alone does not validate the empirical-null assumption. No combined FDR claim is made across contrasts, tissues or feature levels.

## Browser and terminal

In `fun_proteomics setup`, select **Advanced settings → Intensity FDR method** for DIA or matrix analysis. The advanced DIA editor also saves this choice. The terminal/scripted setup accepts the same setting:

```bash
fun_proteomics setup --terminal --non-interactive \
  --mode matrix --project ./fdrtool-analysis \
  --matrix /path/to/protein_intensities.csv \
  --samples /path/to/samples.csv \
  --independent-samples --normalization none \
  --fdr-method fdrtool_pvalue --run
```

For a direct import:

```bash
fun_proteomics quantify \
  --matrix /path/to/protein_intensities.csv \
  --samples /path/to/samples.csv --outdir ./fdrtool-results \
  --differential --fdr-method fdrtool_pvalue
```

Use `--fdr-method fdrtool_normal` only when you explicitly want empirical-normal recalibration. For paired matrix data, use `--subject-column donor` and omit the independent-sample confirmation in guided setup. Changing the FDR method requires a new output/project; an existing run retains its saved settings.

The native DIA JSON setting is:

```json
{
  "differential": {
    "enabled": true,
    "comparisons": "all_pairwise",
    "fdr": 0.05,
    "fdr_method": "fdrtool_pvalue",
    "min_replicates": 2,
    "normalization": "none",
    "gene_level": true
  }
}
```

These methods apply to intensity comparisons. DDA count models and other workflow modes retain their existing methods.

## Installation and outputs

Install the R dependency with `fun_proteomics install --r-only` (or `install --conda` for the full runtime). With a custom R library, install `fdrtool` into that library as well. A missing dependency or failed empirical-null fit stops the requested analysis; BH is never silently substituted.

Alongside `all_results.csv/.tsv` and `significant_results.csv/.tsv`, each fdrtool contrast contains:

- `fdrtool_parameters.csv/.tsv`: input type, cutoff, tested features, package version, fitted null parameters and their standard errors.
- `fdrtool_diagnostics.csv/.tsv`: the input statistic, original limma P, calibrated P, q and local fdr for every tested feature.
- `fdrtool_diagnostics.pdf`: distribution and FDR diagnostics from the same fitted object.
- `fdrtool_warnings.txt`: all warnings emitted during fitting, or an explicit statement that none were emitted.

Untestable features remain in the full table with blank inferential values. The input matrix and sample design remain available; no missing values are imputed.

Method references: [fdrtool documentation](https://strimmerlab.github.io/software/fdrtool/), [package manual](https://cran.r-project.org/web/packages/fdrtool/fdrtool.pdf), and [R's BH adjustment documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html).
