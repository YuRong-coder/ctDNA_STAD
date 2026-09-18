# ctDNA_STAD

Code for the gastric cancer ctDNA Figure 1 and Figure 2 analyses.

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

## Figure 2

Place the following input files in `data/`:

- `Gastric_all_sample_multianno.txt`
- `oncoplot_matrix_CFDNA_Gastric.txt`
- `oncoplot_matrix_TCGA_STAD.txt`
- `patient_info.txt`

Then run:

```bash
Rscript Figure2_analysis.R
```

Panels A-D and their intermediate tables are written to `output/Figure2/`.
