############################################################
## Figure S4: Individual survival analysis of seven risk genes
## Run with: Rscript FigureS4_analysis.R
############################################################

required_packages <- c("survival", "survminer", "ggplot2", "gridExtra")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) stop("缺少R包：", paste(missing_packages, collapse = ", "))

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(ggplot2)
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
input_file <- file.path(
  script_dir, "output", "Figure5", "AC_prognostic_model",
  "04_expression_survival_clinical_merged.csv"
)
output_dir <- file.path(script_dir, "output", "FigureS4")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_file)) stop("Input file not found: ", input_file)
setwd(output_dir)
genes <- c("COL15A1", "PLXNC1", "CDH11", "CCL16", "CLDN11", "ERBB4", "CDH6")
dat <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)

required_columns <- c("patient_id", "OS_time", "OS_event", genes)
missing_columns <- setdiff(required_columns, colnames(dat))
if (length(missing_columns) > 0) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

dat$OS_time <- as.numeric(dat$OS_time)
dat$OS_event <- as.numeric(dat$OS_event)
dat <- dat[!is.na(dat$OS_time) & dat$OS_time > 0 & dat$OS_event %in% c(0, 1), ]

stats_list <- list()
plot_list <- list()
group_data_list <- list()

for (gene in genes) {
  gene_dat <- dat[, c("patient_id", "OS_time", "OS_event", gene)]
  colnames(gene_dat)[4] <- "Expression"
  gene_dat <- gene_dat[complete.cases(gene_dat), ]
  cutoff <- median(gene_dat$Expression, na.rm = TRUE)
  gene_dat$Group <- factor(
    ifelse(gene_dat$Expression > cutoff, "High", "Low"),
    levels = c("Low", "High")
  )
  gene_dat$Gene <- gene

  km_fit <- survfit(Surv(OS_time, OS_event) ~ Group, data = gene_dat)
  surv_diff <- survdiff(Surv(OS_time, OS_event) ~ Group, data = gene_dat)
  logrank_p <- 1 - pchisq(surv_diff$chisq, df = length(surv_diff$n) - 1)

  group_cox <- summary(coxph(Surv(OS_time, OS_event) ~ Group, data = gene_dat))
  continuous_cox <- summary(coxph(Surv(OS_time, OS_event) ~ Expression, data = gene_dat))

  stats_list[[gene]] <- data.frame(
    Gene = gene,
    N = nrow(gene_dat),
    Events = sum(gene_dat$OS_event),
    Median_cutoff = cutoff,
    Low_N = sum(gene_dat$Group == "Low"),
    High_N = sum(gene_dat$Group == "High"),
    High_vs_Low_HR = unname(group_cox$coefficients[1, "exp(coef)"]),
    Group_lower95 = unname(group_cox$conf.int[1, "lower .95"]),
    Group_upper95 = unname(group_cox$conf.int[1, "upper .95"]),
    Group_Cox_P = unname(group_cox$coefficients[1, "Pr(>|z|)"]),
    Logrank_P = logrank_p,
    Continuous_HR = unname(continuous_cox$coefficients[1, "exp(coef)"]),
    Continuous_lower95 = unname(continuous_cox$conf.int[1, "lower .95"]),
    Continuous_upper95 = unname(continuous_cox$conf.int[1, "upper .95"]),
    Continuous_Cox_P = unname(continuous_cox$coefficients[1, "Pr(>|z|)"]),
    stringsAsFactors = FALSE
  )

  hr_label <- sprintf(
    "HR = %.2f (95%% CI %.2f-%.2f)",
    group_cox$coefficients[1, "exp(coef)"],
    group_cox$conf.int[1, "lower .95"],
    group_cox$conf.int[1, "upper .95"]
  )

  gene_dat$OS_months <- gene_dat$OS_time / 30.4375
  km_fit_plot <- survfit(Surv(OS_months, OS_event) ~ Group, data = gene_dat)
  xmax <- ceiling(max(gene_dat$OS_months) / 25) * 25
  p_label <- if (logrank_p < 0.0001) "p < 0.0001" else sprintf("p = %.3g", logrank_p)

  km <- ggsurvplot(
    km_fit_plot,
    data = gene_dat,
    risk.table = TRUE,
    risk.table.height = 0.24,
    risk.table.title = "Number at risk",
    risk.table.y.text.col = TRUE,
    risk.table.y.text = TRUE,
    pval = FALSE,
    pval.method = FALSE,
    conf.int = FALSE,
    censor = TRUE,
    palette = c("#4DBBD5", "#E64B35"),
    legend = "none",
    xlim = c(0, xmax),
    break.time.by = 25,
    xlab = "",
    ylab = "OS probability",
    title = gene,
    ggtheme = theme_classic(base_size = 12),
    tables.theme = theme_cleantable(base_size = 10)
  )
  km$plot <- km$plot +
    annotate("text", x = xmax * 0.055, y = 0.16,
             label = p_label, hjust = 0, size = 3.8) +
    annotate("text", x = xmax * 0.055, y = 0.08,
             label = hr_label, hjust = 0, size = 3.5) +
    annotate("text", x = xmax * 0.67, y = 0.94,
             label = "Low expression", color = "#4DBBD5", hjust = 0, size = 4.0) +
    annotate("text", x = xmax * 0.67, y = 0.87,
             label = "High expression", color = "#E64B35", hjust = 0, size = 4.0) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.border = element_rect(color = "#333333", fill = NA, linewidth = 0.8),
      axis.line = element_blank(),
      plot.margin = margin(6, 8, 2, 8)
    )
  km$table <- km$table +
    labs(x = "OS time (months)", y = NULL) +
    theme_classic(base_size = 11) +
    theme(
      panel.border = element_rect(color = "#333333", fill = NA, linewidth = 0.8),
      axis.line = element_blank(),
      plot.title = element_text(hjust = 0, size = 11),
      legend.position = "none",
      plot.margin = margin(2, 8, 6, 8)
    )

  full_plot <- gridExtra::arrangeGrob(
    km$plot, km$table, ncol = 1, heights = c(3.7, 1.15)
  )

  pdf_file <- file.path(output_dir, paste0("TCGA_STAD_", gene, "_OS_KM.pdf"))
  png_file <- file.path(output_dir, paste0("TCGA_STAD_", gene, "_OS_KM.png"))
  ggsave(pdf_file, full_plot, width = 6.4, height = 6.3, device = cairo_pdf)
  ggsave(png_file, full_plot, width = 6.4, height = 6.3, dpi = 300, bg = "white")

  # The main KM panel is retained separately for the combined multi-panel figure.
  plot_list[[gene]] <- full_plot
  group_data_list[[gene]] <- gene_dat
}

stats_df <- do.call(rbind, stats_list)
rownames(stats_df) <- NULL
stats_df$FDR_Logrank <- p.adjust(stats_df$Logrank_P, method = "BH")
stats_df$FDR_Group_Cox <- p.adjust(stats_df$Group_Cox_P, method = "BH")
stats_df$FDR_Continuous_Cox <- p.adjust(stats_df$Continuous_Cox_P, method = "BH")

write.csv(
  stats_df,
  file.path(output_dir, "TCGA_STAD_7genes_individual_OS_statistics.csv"),
  row.names = FALSE
)
write.csv(
  do.call(rbind, group_data_list),
  file.path(output_dir, "TCGA_STAD_7genes_individual_OS_group_data.csv"),
  row.names = FALSE
)

# Combined seven-panel KM figure, retaining the framed risk tables.
if (!requireNamespace("gridExtra", quietly = TRUE)) {
  stop("Package 'gridExtra' is required to create the combined figure.")
}
combined <- gridExtra::arrangeGrob(grobs = unname(plot_list), ncol = 2, nrow = 4)

ggsave(
  file.path(output_dir, "TCGA_STAD_7genes_individual_OS_KM_combined.pdf"),
  combined,
  width = 12.8,
  height = 24,
  device = cairo_pdf
)
ggsave(
  file.path(output_dir, "TCGA_STAD_7genes_individual_OS_KM_combined.png"),
  combined,
  width = 12.8,
  height = 24,
  dpi = 300,
  bg = "white"
)

writeLines(
  c(
    "TCGA-STAD seven-gene individual overall-survival analysis",
    paste0("Input: ", input_file),
    paste0("Eligible patients: ", nrow(dat)),
    paste0("Deaths: ", sum(dat$OS_event)),
    "Grouping: expression above the median = High; expression at or below the median = Low.",
    "Group HR compares High with Low expression.",
    "P values are reported for log-rank, dichotomized Cox, and continuous-expression Cox analyses.",
    "FDR values were calculated across the seven genes using the Benjamini-Hochberg method."
  ),
  file.path(output_dir, "TCGA_STAD_7genes_individual_OS_run_notes.txt")
)

print(stats_df)
message("All outputs were written to: ", output_dir)
