# ================================================================
# 胃癌患者原发肿瘤部位分布图（患者水平）
# ================================================================
# 输入文件：
#   1. gastric_all_samples_clinical_unknown_filled.txt
#      用于确定本研究纳入的患者ID；同一患者可能有多个样本。
#   2. 临床信息.xlsx
#      使用“patients with gastric cancer”工作表中的Tumor location(site)。
#
# 输出文件（均保存在Figure1文件夹）：
#   Figure_J_tumor_location.pdf
#   Figure_J_tumor_location.png
#   tumor_location_patient_data.csv
#   tumor_location_summary.csv
#
# 统计原则：
#   - 以Patient_ID为单位去重，避免重复采样导致重复计数。
#   - 缺失或Unknown部位不纳入百分比的分母。
#   - E和EGJ合并为EGJ；跨部位记录归为Overlapping。
# ================================================================

rm(list = ls())

# ---------------------- 1. 加载依赖包 -----------------------------
required_packages <- c("ggplot2", "readxl")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "缺少R包：", paste(missing_packages, collapse = ", "),
    "。请先运行 install.packages(c(",
    paste(sprintf("'%s'", missing_packages), collapse = ", "), "))."
  )
}
library(ggplot2)

# ---------------------- 2. 文件路径 --------------------------------
data_dir <- "data"
out_dir <- "output"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cohort_file <- file.path(
  data_dir, "gastric_all_samples_clinical_unknown_filled.txt"
)
clinical_xlsx <- file.path(data_dir, "临床信息.xlsx")
stopifnot(file.exists(cohort_file), file.exists(clinical_xlsx))

# ---------------------- 3. 确定研究队列（患者水平） -----------------
cohort <- read.delim(
  cohort_file, header = TRUE, sep = "\t", check.names = FALSE,
  stringsAsFactors = FALSE, fileEncoding = "UTF-8"
)
stopifnot("Patient_ID" %in% names(cohort))

cohort_ids <- unique(trimws(as.character(cohort$Patient_ID)))
cohort_ids <- cohort_ids[
  nzchar(cohort_ids) & !is.na(cohort_ids) & cohort_ids != "Unknown"
]

# ---------------------- 4. 读取并匹配肿瘤部位 -----------------------
needed <- c("Patient ID", "Tumor location(site)")

# 不固定工作表名称：不同版本的Excel可能命名为Sheet1，或使用疾病名称。
# 自动选择同时包含患者ID和肿瘤部位字段的工作表，避免因改名而报错。
sheet_names <- readxl::excel_sheets(clinical_xlsx)
matched_sheet <- NA_character_
for (sheet_name in sheet_names) {
  header_names <- names(readxl::read_excel(
    clinical_xlsx, sheet = sheet_name, n_max = 0
  ))
  if (all(needed %in% header_names)) {
    matched_sheet <- sheet_name
    break
  }
}
if (is.na(matched_sheet)) {
  stop(
    "未找到同时包含以下字段的工作表：",
    paste(needed, collapse = ", "),
    "。当前工作表：", paste(sheet_names, collapse = ", ")
  )
}
message("读取临床信息工作表：", matched_sheet)
clinical <- readxl::read_excel(clinical_xlsx, sheet = matched_sheet)

location_data <- data.frame(
  Patient_ID = trimws(as.character(clinical[["Patient ID"]])),
  Location_raw = trimws(as.character(clinical[["Tumor location(site)"]])),
  stringsAsFactors = FALSE
)
location_data <- location_data[
  location_data$Patient_ID %in% cohort_ids, , drop = FALSE
]

# 如果Excel中同一患者出现多行，仅保留第一条；同时检查部位记录是否冲突。
conflict_check <- tapply(
  location_data$Location_raw, location_data$Patient_ID,
  function(x) length(unique(x[!is.na(x) & nzchar(x)]))
)
if (any(conflict_check > 1L, na.rm = TRUE)) {
  warning("部分患者存在多个不同的肿瘤部位记录，请检查患者级数据。")
}
location_data <- location_data[!duplicated(location_data$Patient_ID), ]

# ---------------------- 5. 统一部位名称 -----------------------------
location_data$Tumor_location <- NA_character_
location_data$Tumor_location[location_data$Location_raw == "L"] <- "Lower"
location_data$Tumor_location[location_data$Location_raw == "M"] <- "Middle"
location_data$Tumor_location[location_data$Location_raw == "U"] <- "Upper"
location_data$Tumor_location[
  location_data$Location_raw %in% c("E", "EGJ")
] <- "EGJ"
location_data$Tumor_location[
  location_data$Location_raw %in% c("M/L", "U/L", "U/M", "U/M/L")
] <- "Overlapping"

valid_levels <- c("Lower", "Middle", "Upper", "EGJ", "Overlapping")
plot_data <- location_data[
  location_data$Tumor_location %in% valid_levels, , drop = FALSE
]
plot_data$Tumor_location <- factor(
  plot_data$Tumor_location, levels = valid_levels
)

# 保存真正进入统计的患者级数据，便于复核。
write.csv(
  plot_data, file.path(out_dir, "tumor_location_patient_data.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

# ---------------------- 6. 汇总人数和比例 ---------------------------
location_count <- table(plot_data$Tumor_location)
summary_data <- data.frame(
  Tumor_location = factor(names(location_count), levels = valid_levels),
  Count = as.integer(location_count),
  stringsAsFactors = FALSE
)
summary_data$Percentage <- 100 * summary_data$Count / sum(summary_data$Count)
summary_data$Label <- sprintf(
  "%d (%.1f%%)", summary_data$Count, summary_data$Percentage
)

write.csv(
  summary_data, file.path(out_dir, "tumor_location_summary.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

# ---------------------- 7. 绘制横向柱状图 ---------------------------
# 使用Figure 1中已有的蓝色体系；横向柱状图便于显示较长的分类名称。
p_location <- ggplot(
  summary_data,
  aes(x = Tumor_location, y = Percentage, fill = Tumor_location)
) +
  geom_col(width = 0.66, color = "white", linewidth = 0.45) +
  geom_text(
    aes(label = Label), hjust = -0.12,
    size = 3.6, fontface = "bold", color = "#333333"
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = c(
    Lower = "#4AA9C0",
    Middle = "#77BDD0",
    Upper = "#A2D1DB",
    EGJ = "#E6B84F",
    Overlapping = "#B9B9B9"
  )) +
  scale_y_continuous(
    limits = c(0, max(summary_data$Percentage) * 1.30),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "Tumor location",
    x = NULL,
    y = "Percentage"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold", size = 13, hjust = 0.5,
      margin = margin(b = 7)
    ),
    axis.text = element_text(color = "black"),
    axis.title.y = element_blank(),
    legend.position = "none",
    plot.margin = margin(6, 20, 6, 6)
  )

# ---------------------- 8. 输出矢量PDF和高分辨率PNG -----------------
ggsave(
  file.path(out_dir, "Figure_J_tumor_location.pdf"),
  plot = p_location, width = 4.2, height = 3.3,
  units = "in", device = "pdf", useDingbats = FALSE
)
ggsave(
  file.path(out_dir, "Figure_J_tumor_location.png"),
  plot = p_location, width = 4.2, height = 3.3,
  units = "in", dpi = 600, bg = "white"
)

# ---------------------- 9. 打印运行摘要 -----------------------------
message("Tumor location绘图完成。")
message("研究队列患者数：", length(cohort_ids))
message("匹配到Excel的患者数：", nrow(location_data))
message("具有有效部位信息的患者数：", nrow(plot_data))
message("缺失/未识别部位数：", nrow(location_data) - nrow(plot_data))
message("输出目录：", out_dir)
