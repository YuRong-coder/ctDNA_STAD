# ================================================================
# 胃癌平均 VAF 密度图与各分期 Lauren 分子亚型堆叠图
# ================================================================
# 本脚本读取原始 txt 文件，生成两张 PDF，并保留绘图中间数据。
# 原始文件不会被修改。
#
# 输入：
#   sample_VAF_statistics.txt
#   gastric_all_samples_clinical_unknown_filled.txt
#
# 输出目录：output_vaf_stage_subtype_R
#   Figure2_mean_VAF_density_R.pdf
#   Figure3_stage_molecular_subtype_stacked_R.pdf
#   mean_VAF_plot_data.csv
#   stage_lauren_patient_data.csv
#   stage_lauren_summary.csv
# ================================================================
# 清空Global环境所有变量、数据框、列表
rm(list = ls())
# ---------------------- 1. 加载绘图包 ---------------------------
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("缺少 ggplot2，请先运行 install.packages('ggplot2')。")
}
library(ggplot2)

# ---------------------- 2. 【修改点】手动指定文件路径，注释自动获取脚本目录函数 -------------------------
# 注释原有自动获取目录代码
# get_script_dir <- function() {
#   args <- commandArgs(trailingOnly = FALSE)
#   file_arg <- grep("^--file=", args, value = TRUE)
#   if (length(file_arg) == 1L) {
#     return(dirname(normalizePath(sub("^--file=", "", file_arg))))
#   }
#   normalizePath(getwd())
# }
# base_dir <- get_script_dir()

# 你的文件存放根目录（Windows路径统一用正斜杠 /）
base_dir <- "data"
# 输出文件夹放在同目录下
out_dir <- file.path(base_dir, "output_vaf_stage_subtype_R3")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 输入文件完整路径
vaf_file <- file.path(base_dir, "sample_VAF_statistics.txt")
clinical_file <- file.path(base_dir, "gastric_all_samples_clinical_unknown_filled.txt")

# 校验文件是否存在并打印路径，方便排查
cat("VAF文件检索路径：", vaf_file, "\n")
cat("临床文件检索路径：", clinical_file, "\n")
stopifnot(file.exists(vaf_file), file.exists(clinical_file))

# ---------------------- 3. 平均 VAF 数据 ------------------------
# Avg_VAF 已在原始统计表中按样本计算；转换为数值后剔除缺失值。
vaf_raw <- read.delim(
  vaf_file, header = TRUE, sep = "\t", check.names = FALSE,
  stringsAsFactors = FALSE, fileEncoding = "UTF-8"
)
stopifnot(all(c("Sample", "Avg_VAF") %in% names(vaf_raw)))

vaf_data <- data.frame(
  Sample = as.character(vaf_raw$Sample),
  Avg_VAF = suppressWarnings(as.numeric(vaf_raw$Avg_VAF)),
  stringsAsFactors = FALSE
)
vaf_data <- vaf_data[
  nzchar(vaf_data$Sample) & is.finite(vaf_data$Avg_VAF), , drop = FALSE
]

# 保存真正进入密度图的样本级中间数据。
write.csv(
  vaf_data, file.path(out_dir, "mean_VAF_plot_data.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

vaf_n <- nrow(vaf_data)
vaf_mean <- mean(vaf_data$Avg_VAF)
vaf_median <- median(vaf_data$Avg_VAF)

# ---------------------- 4. 平均 VAF 密度图 ----------------------
# 蓝色曲线为高斯核密度估计；红色虚线为均值，灰色点线为中位数。
p_vaf <- ggplot(vaf_data, aes(x = Avg_VAF)) +
  geom_density(
    fill = "#80AFC8", color = "#1D6BA8",
    alpha = 0.45, linewidth = 1.15, adjust = 1
  ) +
  geom_vline(
    xintercept = vaf_mean, color = "#CA4B48",
    linewidth = 0.75, linetype = "dashed"
  ) +
  geom_vline(
    xintercept = vaf_median, color = "#4A4A4A",
    linewidth = 0.75, linetype = "dotted"
  ) +
  coord_cartesian(xlim = c(0.38, 0.69), expand = FALSE) +
  labs(
    title = "Mean VAF distribution",
    subtitle = sprintf(
      "n = %d   Mean = %.3f   Median = %.3f",
      vaf_n, vaf_mean, vaf_median
    ),
    x = "Average variant allele frequency (VAF)",
    y = "Density"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(color = "#555555", margin = margin(b = 8)),
    panel.grid.major.y = element_line(color = "#E1E5EB", linetype = "dashed"),
    panel.grid.minor = element_blank()
  )

# 使用标准 pdf() 设备，避免 cairo 依赖造成兼容性问题。
pdf(
  file.path(out_dir, "Figure2_mean_VAF_density_R.pdf"),
  width = 7.0, height = 5.2, onefile = FALSE,
  family = "Helvetica", useDingbats = FALSE
)
print(p_vaf)
dev.off()

# ---------------------- 5. 临床分期与 Lauren 分型 --------------
clinical <- read.delim(
  clinical_file, header = TRUE, sep = "\t", check.names = FALSE,
  stringsAsFactors = FALSE, fileEncoding = "UTF-8"
)
needed <- c("Tumor_Sample_Barcode", "Patient_ID", "Stage", "Lauren")
stopifnot(all(needed %in% names(clinical)))

# 将亚分期合并为 I、II、III、IV 四组。
# Ia/Ib -> I；IIa/IIb -> II；IIIa/IIIb/IIIc -> III。
stage_raw <- trimws(as.character(clinical$Stage))
stage_group <- rep(NA_character_, length(stage_raw))
stage_group[grepl("^I[ab]?$", stage_raw, ignore.case = TRUE)] <- "I"
stage_group[grepl("^II[ab]?$", stage_raw, ignore.case = TRUE)] <- "II"
stage_group[grepl("^III[abc]?$", stage_raw, ignore.case = TRUE)] <- "III"
stage_group[toupper(stage_raw) == "IV"] <- "IV"

patient_data <- data.frame(
  Sample = as.character(clinical$Tumor_Sample_Barcode),
  Patient_ID = as.character(clinical$Patient_ID),
  Stage_raw = stage_raw,
  Stage_group = stage_group,
  Lauren = trimws(as.character(clinical$Lauren)),
  stringsAsFactors = FALSE
)

# 仅纳入分期明确且 Lauren 属于三种标准类型的患者。
valid_lauren <- c("Intestinal", "Diffuse", "Mixed")
patient_data <- patient_data[
  !is.na(patient_data$Stage_group) & patient_data$Lauren %in% valid_lauren,
  , drop = FALSE
]
patient_data$Stage_group <- factor(
  patient_data$Stage_group, levels = c("I", "II", "III", "IV")
)
patient_data$Lauren <- factor(patient_data$Lauren, levels = valid_lauren)

# 保存进入堆叠图的患者级中间数据。
write.csv(
  patient_data, file.path(out_dir, "stage_lauren_patient_data.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

# ---------------------- 6. 计算分期内人数与比例 ----------------
# expand.grid 确保 4 个分期 × 3 个亚型的组合全部存在。
summary_data <- expand.grid(
  Stage_group = factor(c("I", "II", "III", "IV"),
                       levels = c("I", "II", "III", "IV")),
  Lauren = factor(valid_lauren, levels = valid_lauren),
  KEEP.OUT.ATTRS = FALSE
)

summary_data$count <- mapply(
  function(s, l) sum(
    patient_data$Stage_group == s & patient_data$Lauren == l,
    na.rm = TRUE
  ),
  as.character(summary_data$Stage_group),
  as.character(summary_data$Lauren)
)

stage_totals <- tapply(
  summary_data$count, summary_data$Stage_group, sum
)
summary_data$stage_total <- as.integer(
  stage_totals[as.character(summary_data$Stage_group)]
)
summary_data$proportion <- ifelse(
  summary_data$stage_total > 0,
  summary_data$count / summary_data$stage_total,
  0
)
summary_data$label <- ifelse(
  summary_data$count > 0,
  sprintf("%d\n(%.0f%%)", summary_data$count, 100 * summary_data$proportion),
  ""
)

# 按分期、亚型顺序排序并保存汇总中间数据。
summary_data <- summary_data[
  order(summary_data$Stage_group, summary_data$Lauren), , drop = FALSE
]
write.csv(
  summary_data, file.path(out_dir, "stage_lauren_summary.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

# ---------------------- 7. 100% 堆叠柱状图 ---------------------
subtype_colors <- c(
  Intestinal = "#6BAED6",
  Diffuse = "#FDBF6F",
  Mixed = "#a6a6a6"
)
axis_labels <- paste0(
  names(stage_totals), "\n", "n=", as.integer(stage_totals)
)
names(axis_labels) <- names(stage_totals)

p_stage <- ggplot(
  summary_data,
  aes(x = Stage_group, y = proportion, fill = Lauren)
) +
  geom_col(width = 0.62, color = "#333333", linewidth = 0.45) +
  geom_text(
    aes(label = label), position = position_stack(vjust = 0.5),
    color = "white", fontface = "bold", size = 4.0, lineheight = 0.95
  ) +
  scale_fill_manual(values = subtype_colors, drop = FALSE) +
  scale_x_discrete(labels = axis_labels, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, 1), breaks = seq(0, 1, 0.25),
    labels = function(x) paste0(round(100 * x), "%"),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "Molecular subtype composition by stage",
    subtitle = "Lauren classification; Unknown values excluded",
    x = "Pathological stage", y = "Proportion within stage", fill = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(color = "#555555", margin = margin(b = 8)),
    axis.text.x = element_text(face = "bold", lineheight = 1.2),
    panel.grid.major.y = element_line(color = "#E1E5EB"),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.justification = "right"
  )

pdf(
  file.path(out_dir, "Figure3_stage_molecular_subtype_stacked_R.pdf"),
  width = 7.2, height = 5.4, onefile = FALSE,
  family = "Helvetica", useDingbats = FALSE
)
print(p_stage)
dev.off()

# ---------------------- 8. Lauren 分型的平均 VAF 密度图 --------
# 将样本平均 VAF 与临床表按样本名精确匹配。
# 仅保留 Lauren 分型明确的 Intestinal、Diffuse、Mixed 三组。
vaf_lauren_data <- merge(
  vaf_data,
  data.frame(
    Sample = as.character(clinical$Tumor_Sample_Barcode),
    Lauren = trimws(as.character(clinical$Lauren)),
    stringsAsFactors = FALSE
  ),
  by = "Sample",
  all = FALSE
)
vaf_lauren_data <- vaf_lauren_data[
  vaf_lauren_data$Lauren %in% valid_lauren, , drop = FALSE
]
vaf_lauren_data$Lauren <- factor(
  vaf_lauren_data$Lauren,
  levels = valid_lauren
)

# 保存进入 Lauren 分组密度图的样本级数据。
write.csv(
  vaf_lauren_data,
  file.path(out_dir, "mean_VAF_by_Lauren_plot_data.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# 计算各 Lauren 分型的样本量、均值、中位数和标准差。
lauren_summary <- do.call(
  rbind,
  lapply(valid_lauren, function(type) {
    x <- vaf_lauren_data$Avg_VAF[vaf_lauren_data$Lauren == type]
    data.frame(
      Lauren = type,
      n = length(x),
      Mean_VAF = mean(x),
      Median_VAF = median(x),
      SD_VAF = sd(x),
      stringsAsFactors = FALSE
    )
  })
)
write.csv(
  lauren_summary,
  file.path(out_dir, "mean_VAF_by_Lauren_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# 在图例中同时标注各组样本量。
lauren_legend_labels <- setNames(
  paste0(lauren_summary$Lauren, " (n=", lauren_summary$n, ")"),
  lauren_summary$Lauren
)

# 三条密度曲线使用相同带宽规则和相同横轴，便于比较分布位置和形状。
# 半透明填充可显示不同 Lauren 分型之间的重叠区域。
p_vaf_lauren <- ggplot(
  vaf_lauren_data,
  aes(x = Avg_VAF, color = Lauren, fill = Lauren)
) +
  geom_density(alpha = 0.20, linewidth = 1.05, adjust = 1) +
  scale_color_manual(
    values = subtype_colors,
    labels = lauren_legend_labels,
    drop = FALSE
  ) +
  scale_fill_manual(
    values = subtype_colors,
    labels = lauren_legend_labels,
    drop = FALSE
  ) +
  coord_cartesian(xlim = c(0.38, 0.69), expand = FALSE) +
  labs(
    title = "Mean VAF distribution by Lauren subtype",
    subtitle = "Sample-level mean VAF; Unknown Lauren values excluded",
    x = "Average variant allele frequency (VAF)",
    y = "Density",
    color = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(color = "#555555", margin = margin(b = 8)),
    panel.grid.major.y = element_line(color = "#E1E5EB", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.justification = "left"
  )

pdf(
  file.path(out_dir, "Figure2_mean_VAF_density_by_Lauren_R.pdf"),
  width = 7.2, height = 5.4, onefile = FALSE,
  family = "Helvetica", useDingbats = FALSE
)
print(p_vaf_lauren)
dev.off()

# 同步输出 PNG 预览，用于快速检查 PDF 的版式与配色。
png(
  file.path(out_dir, "Figure2_mean_VAF_density_by_Lauren_R.png"),
  width = 7.2, height = 5.4, units = "in", res = 300,
  bg = "white", type = "windows"
)
print(p_vaf_lauren)
dev.off()

# ---------------------- 9. 输出运行摘要 -------------------------
message("R 绘图完成。")
message(sprintf(
  "平均 VAF：n=%d，mean=%.4f，median=%.4f",
  vaf_n, vaf_mean, vaf_median
))
message("分期/Lauren 有效患者数：", nrow(patient_data))
message(
  "Lauren 分型 VAF 图样本量：",
  paste0(lauren_summary$Lauren, "=", lauren_summary$n, collapse = "；")
)
message("输出目录：", out_dir)
