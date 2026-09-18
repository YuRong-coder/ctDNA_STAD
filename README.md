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

## Figure 3

Place the following input files in `data/`:

- `gastric_maf_with_clinical.rds`
- `oncoplot_matrix_TCGA_STAD.txt`
- `gastric_TME_gene_literature_support_intersection_optimized_60.csv`

Then run:

```bash
Rscript Figure3_analysis.R
```

Panels A-D are written to separate subdirectories under `output/Figure3/`.

## Figure 4

Place the following TCGA-STAD/Xena input files in `data/`:

- `TCGA-STAD.somaticmutation_wxs.tsv`
- `TCGA-STAD.clinical.gz`
- `TCGA-STAD.survival.gz`
- `oncoplot_matrix_TCGA_STAD.txt`

Then run:

```bash
Rscript Figure4_analysis.R
```

Panel A uses the typed 20-gene waterfall workflow. Panels B-D use the Xena
Stage, pT and pN dot-plot workflows. Outputs are written to separate
subdirectories under `output/Figure4/`.

## Figure 5

Place the following TCGA-STAD input files in `data/`:

- `TCGA-STAD.star_tpm.tsv`
- `TCGA-STAD.survival.gz`
- `TCGA-STAD.clinical.gz`

Then run:

```bash
Rscript Figure5_analysis.R
```

The script produces the univariable Cox forest plot (A), multivariable Cox
forest plot (B), LASSO risk-group Kaplan-Meier curve (C), and paired
tumor-versus-normal expression analysis (D) under `output/Figure5/`.

## Figure 6

Copy the existing GEO directory trees into:

- `data/Figure6/GEO/`
- `data/Figure6/GEO2/`

Then run:

```bash
Rscript Figure6_analysis.R
```

Independent GEO-cohort survival validation results, audit tables, and
intermediate data are written under `output/Figure6/`.

## Supplementary figures

Run the supplementary analyses independently:

```bash
Rscript FigureS2_analysis.R
Rscript FigureS3_analysis.R
Rscript FigureS4_analysis.R
```

- Figure S2 reads `gastric_all_samples_clinical_unknown_filled.txt` and
  `all_exonic_VAF.txt` from `data/`.
- Figure S3 reads `TCGA-STAD.somaticmutation_wxs.tsv.gz` and
  `TCGA-STAD.clinical.gz` from `data/`.
- Figure S4 reads the Figure 5 merged expression-survival table generated under
  `output/Figure5/AC_prognostic_model/`.

Results are written independently to `output/FigureS2/`, `output/FigureS3/`,
and `output/FigureS4/`.

Figure S5 is a combined three-panel LASSO/Cox workflow:

```bash
Rscript FigureS5_analysis.R
```

It generates the LASSO cross-validation curve (A), selected-gene coefficient
plot (B), and univariable Cox forest plot (C) under `output/FigureS5/`.
