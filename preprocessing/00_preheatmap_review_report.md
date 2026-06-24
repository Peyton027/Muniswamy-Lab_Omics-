# RNA-seq limma-voom pre-heatmap review

## Scope completed

This package contains the requested **pre-heatmap** workflow outputs only. No final all-target or category heatmap has been generated.

All 35 samples were retained and ordered as **Normal → Obese → Diabetic → HFrEF**. `D8S` was treated as Diabetic.

## Step 1 — input validation

- Raw matrix: **58,735 genes**, `gene_id` and `gene_symbol`, plus **35** sample-count columns.
- Metadata: **35 unique sample IDs**, exactly matching the count-matrix columns and their order.
- Cohorts: 9 Normal, 9 Obese, 7 Diabetic, and 10 HFrEF.
- Count values: all finite, nonnegative, and integer-like.
- Raw gene IDs: unique. Raw gene symbols include repeated symbols across different gene IDs, so all target mapping was performed primarily by gene ID.
- Target annotation: **167 unique target gene IDs and symbols**; all map uniquely to the raw count matrix.
- `TMEM94`: exactly one annotation row and exactly one raw-matrix row; required flag is `Yes`.

## Step 3 — low-expression filtering

The cohort-aware filter used a minimum of 10 counts per million threshold equivalent, requiring expression in at least the smallest cohort size (7 samples).

| Metric | Result |
|---|---:|
| Genes before filtering | 58,735 |
| Genes retained | 19,546 |
| Genes removed | 39,189 |
| Target genes retained | 152 |
| Target genes filtered for low expression | 15 |

Target genes filtered for low expression: `CACNA1S`, `CACNA1I`, `SLC8A3`, `PPP3R2`, `KCND2`, `KCNE2`, `DPP10`, `KCNN1`, `NIPAL4`, `CLDN16`, `CLDN19`, `CLCN1`, `SLC4A1`, `SLC12A5`, and `SLC12A1`.

## Step 4 — TMM normalization and QC

- TMM reference sample: `D3`.
- TMM normalization factors ranged from **0.855 to 1.177**.
- Raw, filtered, and TMM-effective library sizes are included in the TMM table.
- Sample-correlation review flags are **review-only**. No sample was removed.
- The lowest median correlations were observed for `O10`, `HF7`, and `O2`; the robust-MAD review rule also flagged `HF1`, `D8S`, `D7`, `N7`, and `D5` for visual review. These are not technical-exclusion decisions.

## Step 5 — voom-scale expression

The filtered genes were transformed to normalized voom-scale log2-CPM values. The mean-variance plot and the complete precision-weight matrix are included.

## Step 6 — formal differential-expression status

The archived files marked `METHOD_MATCHED_PROVISIONAL` were generated in a Python implementation because this execution environment did not contain the official R `edgeR` and `limma` packages. They must **not** be cited as final official limma results. The included R script is the exact edgeR/limma implementation to execute in a suitable R environment for final cohort-level claims.

## Step 7 — individual-level Normal-centered matrix

The pre-heatmap matrix is based on voom-scale log2-CPM, not raw counts and not cohort logFC:

```
centered expression = individual voom log2-CPM − mean voom log2-CPM of the 9 Normal samples for that gene
```

- Matrix dimensions: **152 retained target genes × 35 individual samples**.
- The column order is the required Normal → Obese → Diabetic → HFrEF order.
- The mean centered value across the nine Normal samples is essentially zero for every gene (maximum numerical rounding residual: ~2.37 × 10⁻¹⁵).
- `TMEM94` is retained. Its Normal-group mean voom expression is **5.4124 log2-CPM**; the separate TMEM94 file contains every individual's voom and centered value.

## Review gate before heatmaps

Please review the final Normal-centered matrix, the filtering/mapping audit, TMEM94 file, TMM table, and QC plots. Once approved, the next step is the all-target heatmap and the Calcium, Potassium, Magnesium, Sodium, and Chloride panels using a zero-centered diverging scale.
