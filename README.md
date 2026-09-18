# ctDNA_STAD

Code for the gastric cancer ctDNA Figure 1 analyses (panels B onward).

## Run

Place the following input files in `data/`:

- `gastric_all_samples_clinical_unknown_filled.txt`
- `临床信息.xlsx`
- `sample_VAF_statistics.txt`

Then run:

```bash
Rscript Figure1_analysis.R
```

Figures are written to `output/`. Intermediate analysis tables are written to
`data/` and `output/vaf_stage_subtype/`.
