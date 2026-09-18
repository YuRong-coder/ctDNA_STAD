############################################################
## Figure S2: Clinicopathological and cfDNA-derived features
## Run with: Rscript FigureS2_analysis.R
############################################################

required_packages <- c("data.table", "dplyr", "ggplot2", "cowplot")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) stop("缺少R包：", paste(missing_packages, collapse = ", "))

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
})
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  normalizePath(getwd())
}

script_dir <- get_script_dir()
clinical_file <- file.path(script_dir, "data", "gastric_all_samples_clinical_unknown_filled.txt")
vaf_file <- file.path(script_dir, "data", "all_exonic_VAF.txt")
out_dir <- file.path(script_dir, "output", "FigureS2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(clinical_file)) stop("Input file not found: ", clinical_file)
if (!file.exists(vaf_file)) stop("Input file not found: ", vaf_file)
setwd(out_dir)
clinical_raw <- fread(clinical_file, sep = "\t", header = TRUE, data.table = FALSE)
vaf_raw <- fread(vaf_file, sep = "\t", header = TRUE, data.table = FALSE,
                 select = c("Sample", "VAF"))

vaf_summary <- vaf_raw %>%
  mutate(VAF = as.numeric(VAF)) %>%
  filter(!is.na(VAF)) %>%
  group_by(Sample) %>%
  summarise(Mutation_Count = n(), Avg_VAF = mean(VAF), .groups = "drop")

df <- clinical_raw %>%
  distinct(Tumor_Sample_Barcode, .keep_all = TRUE) %>%
  left_join(vaf_summary, by = c("Tumor_Sample_Barcode" = "Sample")) %>%
  mutate(
    pT_table = case_when(
      grepl("^T1", pT, ignore.case = TRUE) ~ "T1",
      grepl("^T2", pT, ignore.case = TRUE) ~ "T2",
      grepl("^T3", pT, ignore.case = TRUE) ~ "T3",
      grepl("^T4", pT, ignore.case = TRUE) ~ "T4",
      TRUE ~ NA_character_
    ),
    pT_table = factor(pT_table, levels = c("T1", "T2", "T3", "T4")),
    cN_table = case_when(
      pN == "N0" ~ "Negative",
      grepl("^N[1-3]", pN, ignore.case = TRUE) ~ "Positive",
      TRUE ~ NA_character_
    ),
    cN_table = factor(cN_table, levels = c("Negative", "Positive")),
    Stage_table = case_when(
      grepl("^IV", Stage, ignore.case = TRUE) ~ "IV",
      grepl("^III", Stage, ignore.case = TRUE) ~ "III",
      grepl("^II", Stage, ignore.case = TRUE) ~ "II",
      grepl("^I", Stage, ignore.case = TRUE) ~ "I",
      TRUE ~ NA_character_
    ),
    Stage_table = factor(Stage_table, levels = c("I", "II", "III", "IV")),
    Lauren_table = ifelse(Lauren %in% c("Intestinal", "Diffuse", "Mixed"), Lauren, NA_character_),
    Lauren_table = factor(Lauren_table, levels = c("Intestinal", "Diffuse", "Mixed")),
    Differentiation_table = case_when(
      grepl("Highly", Differentiation, ignore.case = TRUE) ~ "Well/Moderate",
      Differentiation == "moderately" ~ "Well/Moderate",
      grepl("Pooly|poorly", Differentiation, ignore.case = TRUE) ~ "Poor/Moderate-poor",
      TRUE ~ NA_character_
    ),
    Differentiation_table = factor(Differentiation_table,
                                   levels = c("Well/Moderate", "Poor/Moderate-poor")),
    log_mut_count = log2(Mutation_Count + 1)
  )

clinical_cols <- c(
  "T1" = "#F3B56B", "T2" = "#E84A3C", "T3" = "#8B1E2D", "T4" = "#5B1020",
  "Negative" = "#67A9CF", "Positive" = "#D63B37"
)
stage_cols <- c("I" = "#F8C991", "II" = "#F4B36A", "III" = "#E6423A", "IV" = "#8B1E2D")
lauren_cols <- c("Intestinal" = "#4DBBD5", "Diffuse" = "#E64B35", "Mixed" = "#8491B4")
diff_cols <- c("Well/Moderate" = "#91D1C2", "Poor/Moderate-poor" = "#DC0000")

fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) "<0.001" else formatC(p, format = "g", digits = 2)
}

bar_count_plot <- function(data, var, title, fill_values) {
  x <- data %>% filter(!is.na(.data[[var]])) %>% count(.data[[var]], name = "n")
  names(x)[1] <- "group"
  ggplot(x, aes(group, n, fill = group)) +
    geom_col(width = 0.68, color = "white", linewidth = 0.4) +
    geom_text(aes(label = n), vjust = -0.25, size = 3.1) +
    scale_fill_manual(values = fill_values, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
    labs(title = title, x = NULL, y = "Number") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 11),
          axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")
}

box_plot_pairwise <- function(data, xvar, yvar, title, ylab, fill_values, ref_group) {
  dd <- data %>% filter(!is.na(.data[[xvar]]), !is.na(.data[[yvar]]))
  dd[[xvar]] <- factor(dd[[xvar]], levels = names(fill_values))
  p <- ggplot(dd, aes(.data[[xvar]], .data[[yvar]], fill = .data[[xvar]])) +
    geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.72,
                 color = "black", linewidth = 0.45) +
    geom_jitter(aes(color = .data[[xvar]]), width = 0.17, size = 1.15,
                alpha = 0.65, show.legend = FALSE) +
    scale_fill_manual(values = fill_values, drop = FALSE) +
    scale_color_manual(values = fill_values, drop = FALSE) +
    labs(title = title, x = NULL, y = ylab) +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 11),
          axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")

  lev <- levels(dd[[xvar]])
  compare_groups <- setdiff(lev, ref_group)
  y_min <- min(dd[[yvar]], na.rm = TRUE)
  y_max <- max(dd[[yvar]], na.rm = TRUE)
  span <- max(y_max - y_min, 0.1)
  step <- span * 0.13
  base_y <- y_max + step * 0.35
  tick <- step * 0.12
  for (i in seq_along(compare_groups)) {
    g <- compare_groups[i]
    ref_vals <- dd[dd[[xvar]] == ref_group, yvar, drop = TRUE]
    cmp_vals <- dd[dd[[xvar]] == g, yvar, drop = TRUE]
    pv <- tryCatch(wilcox.test(ref_vals, cmp_vals, exact = FALSE)$p.value,
                   error = function(e) NA_real_)
    x1 <- match(ref_group, lev); x2 <- match(g, lev); yy <- base_y + (i - 1) * step
    p <- p +
      annotate("segment", x = x1, xend = x2, y = yy, yend = yy, linewidth = 0.4) +
      annotate("segment", x = x1, xend = x1, y = yy, yend = yy - tick, linewidth = 0.4) +
      annotate("segment", x = x2, xend = x2, y = yy, yend = yy - tick, linewidth = 0.4) +
      annotate("text", x = (x1 + x2) / 2, y = yy + tick,
               label = paste0("p=", fmt_p(pv)), size = 2.6)
  }
  p + coord_cartesian(ylim = c(y_min, base_y + length(compare_groups) * step + step * 0.35), clip = "off")
}

pA <- bar_count_plot(df, "pT_table", "Pathological tumor category", clinical_cols)
pB <- bar_count_plot(df, "cN_table", "Nodal status", clinical_cols)
pC <- box_plot_pairwise(df, "Stage_table", "log_mut_count", "Mutation count by stage",
                        "log2(mutation count + 1)", stage_cols, "I")
pD <- box_plot_pairwise(df, "Lauren_table", "Avg_VAF", "Mean VAF by Lauren subtype",
                        "Mean VAF", lauren_cols, "Intestinal")
pE <- box_plot_pairwise(df, "Differentiation_table", "Avg_VAF", "Mean VAF by differentiation",
                        "Mean VAF", diff_cols, "Well/Moderate")

fig <- plot_grid(
  pA, pB, pC, pD, pE, NULL,
  labels = c("A", "B", "C", "D", "E", ""),
  label_fontface = "bold", label_size = 15,
  ncol = 2, align = "hv"
)

ggsave("FigureS2_gastric_clinical_cfDNA_supplement.pdf", fig,
       width = 8.5, height = 11.5, device = cairo_pdf)
ggsave("FigureS2_gastric_clinical_cfDNA_supplement.png", fig,
       width = 8.5, height = 11.5, dpi = 600, bg = "white")

plots <- list(A = pA, B = pB, C = pC, D = pD, E = pE)
plot_names <- c(
  A = "FigureS2A_pathological_tumor_category",
  B = "FigureS2B_nodal_status",
  C = "FigureS2C_mutation_count_by_stage",
  D = "FigureS2D_VAF_by_Lauren_subtype",
  E = "FigureS2E_VAF_by_differentiation"
)
for (nm in names(plots)) {
  ggsave(paste0(plot_names[[nm]], ".pdf"), plots[[nm]], width = 4.5, height = 4.0, device = cairo_pdf)
  ggsave(paste0(plot_names[[nm]], ".png"), plots[[nm]], width = 4.5, height = 4.0, dpi = 300, bg = "white")
}

write.csv(df, "FigureS2_gastric_analysis_data.csv", row.names = FALSE)
write.csv(vaf_summary, "FigureS2_sample_mutation_VAF_summary.csv", row.names = FALSE)

writeLines(c(
  paste0("Clinical records: ", nrow(clinical_raw)),
  paste0("Unique clinical samples: ", nrow(df)),
  paste0("Samples with exonic mutation/VAF data: ", sum(!is.na(df$Mutation_Count))),
  "Mutation count is the number of rows per sample in all_exonic_VAF.txt.",
  "Average VAF is the arithmetic mean of exonic variant VAFs per sample.",
  "Pairwise comparisons use two-sided Wilcoxon rank-sum tests.",
  "The healthy-control pre-library concentration panel in the breast reference was not reproduced because the supplied gastric input has no corresponding concentration variable."
), "FigureS2_run_notes.txt")

message("Figure S2 outputs written to: ", out_dir)
