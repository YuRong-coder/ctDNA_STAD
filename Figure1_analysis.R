############################################################
## Figure 1 unified analysis pipeline (panels B onward)
##
## Put all input files in ./data:
##   - gastric_all_samples_clinical_unknown_filled.txt
##   - 临床信息.xlsx
##   - sample_VAF_statistics.txt
##
## Run with: Rscript Figure1_analysis.R
## Results are written under ./output and intermediate tables under ./data.
############################################################

required_packages <- c("dplyr", "ggplot2", "patchwork", "readxl")
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

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readxl)
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
data_dir <- file.path(script_dir, "data")
root_output_dir <- file.path(script_dir, "output")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(root_output_dir, recursive = TRUE, showWarnings = FALSE)

############################################################
## Part 1: Clinical characteristics
############################################################
source_file <- file.path(data_dir, "gastric_all_samples_clinical_unknown_filled.txt")
raw_clinical_file <- file.path(data_dir, "临床信息.xlsx")
out_dir <- root_output_dir
if (!file.exists(source_file)) stop("Input file not found: ", source_file)
if (!file.exists(raw_clinical_file)) stop("Raw clinical Excel file not found: ", raw_clinical_file)
clinical_raw <- read.delim(
  source_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "UTF-8"
)

# KI67不再使用TXT中预先写死的High/Low，而是读取原始临床Excel后动态计算中位数。
# 原始值存在0.7、75%+、50%-75%等混合格式，统一转换为百分数：
#   0到1之间的小数乘以100；区间取两端均值；“75%+”等取报告的数值75。
ki67_to_percent <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "-", "/", "NA", "N/A", "Unknown", "unknown")] <- NA_character_
  number_list <- regmatches(x, gregexpr("[0-9]+(?:\\.[0-9]+)?", x, perl = TRUE))
  value <- vapply(number_list, function(z) {
    if (length(z) == 0) return(NA_real_)
    mean(as.numeric(z))
  }, numeric(1))
  value[!is.na(value) & value <= 1] <- value[!is.na(value) & value <= 1] * 100
  value
}

ki67_source <- read_excel(raw_clinical_file, sheet = 1)
ki67_required <- c("Patient ID", "KI67")
ki67_missing <- setdiff(ki67_required, names(ki67_source))
if (length(ki67_missing) > 0) {
  stop("Missing KI67 source columns: ", paste(ki67_missing, collapse = ", "))
}
ki67_source <- ki67_source %>%
  transmute(
    Patient_ID = trimws(as.character(.data[["Patient ID"]])),
    Ki67_raw = trimws(as.character(.data[["KI67"]])),
    Ki67_percent = ki67_to_percent(.data[["KI67"]])
  )
# 原始Excel中同一患者可能有空值行和有效复测行；优先保留非缺失值。
# 若同一患者存在两个互相矛盾的有效数值，则停止运行，避免静默选取。
ki67_conflict <- ki67_source %>%
  filter(!is.na(Patient_ID), Patient_ID != "", !is.na(Ki67_percent)) %>%
  group_by(Patient_ID) %>%
  summarise(n_value = n_distinct(Ki67_percent), .groups = "drop") %>%
  filter(n_value > 1)
if (nrow(ki67_conflict) > 0) {
  stop("Conflicting valid KI67 values for Patient ID: ", paste(ki67_conflict$Patient_ID, collapse = ", "))
}
ki67_source <- ki67_source %>%
  filter(!is.na(Patient_ID), Patient_ID != "") %>%
  arrange(is.na(Ki67_percent)) %>%
  distinct(Patient_ID, .keep_all = TRUE)

ki67_cutoff <- median(ki67_source$Ki67_percent, na.rm = TRUE)
if (!is.finite(ki67_cutoff)) stop("No valid numeric KI67 values were found in the source Excel.")
message("KI67 median cutoff: ", ki67_cutoff, "%")

# 按Patient_ID把原始KI67回填到当前绘图样本。未匹配或原始值缺失者保留为NA，
# 后续KI67百分比的分母仍只包含具有有效KI67的样本。
clinical_raw <- clinical_raw %>%
  mutate(Patient_ID = trimws(as.character(Patient_ID))) %>%
  select(-Ki67) %>%
  left_join(ki67_source, by = "Patient_ID")

# 绘图和分组所必需的原始字段。
# 在分析开始前统一检查，可以避免后续mutate或绘图时出现难以定位的列名错误。
required_cols <- c(
  "Tumor_Sample_Barcode", "Patient_ID", "Age", "Gender", "M",
  "Stage", "Lauren", "Differentiation", "HER2", "Ki67_raw", "Ki67_percent"
)
missing_cols <- setdiff(required_cols, names(clinical_raw))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# 将不同形式的无效值统一转换为真正的NA。
# trimws先清除字段两端的空格，防止“ Unknown”等值未被识别。
# “0”在Lauren等文本字段中代表无资料；实际M0、HER2_0不会受影响，因为它们不是单独的“0”。
clean_missing <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "/", "NA", "N/A", "Unknown", "unknown", "0", "-")] <- NA_character_
  x
}

# 在不改动clinical_raw的前提下创建绘图用清洗表clinical_clean。
# 同时保留Stage_raw、HER2_raw等原始标准化字段，方便追溯归并前后的变化。
clinical_clean <- clinical_raw %>%
  mutate(
    # 将Ia、Ib等亚分期归并到I，将IIa、IIb归并到II，III和IV同理。
    # 这样与参考乳腺癌Figure 1的I–IV主分期展示方式一致。
    Stage_raw = clean_missing(Stage),
    Stage_group = case_when(
      grepl("^I($|a|b)", Stage_raw, ignore.case = TRUE) ~ "I",
      grepl("^II($|a|b)", Stage_raw, ignore.case = TRUE) ~ "II",
      grepl("^III($|a|b|c)", Stage_raw, ignore.case = TRUE) ~ "III",
      grepl("^IV($|a|b)", Stage_raw, ignore.case = TRUE) ~ "IV",
      TRUE ~ NA_character_
    ),

    # Lauren分型直接保留Intestinal、Diffuse、Mixed三类；Unknown和0作为缺失值排除。
    Lauren_group = clean_missing(Lauren),

    # 分化程度原始表存在Pooly拼写错误以及“moderately and poorly”等复合描述。
    # 归并优先级为：含poorly/pooly -> Poorly；其次含moderately -> Moderately；
    # 含highly/well -> Well。无法识别的非空值保留原文字，便于发现新的类别。
    Differentiation_raw = clean_missing(Differentiation),
    Differentiation_group = case_when(
      is.na(Differentiation_raw) ~ NA_character_,
      grepl("pooly|poorly", Differentiation_raw, ignore.case = TRUE) ~ "Poorly differentiated",
      grepl("moderately", Differentiation_raw, ignore.case = TRUE) ~ "Moderately differentiated",
      grepl("highly|well", Differentiation_raw, ignore.case = TRUE) ~ "Well differentiated",
      TRUE ~ Differentiation_raw
    ),

    # HER2免疫组化结果按临床常用方式归并：
    # 0/1+为阴性，2+为临界，3+为阳性。
    HER2_raw = clean_missing(HER2),
    HER2_group = case_when(
      HER2_raw %in% c("HER2_0", "HER2_1") ~ "Negative (0/1+)",
      HER2_raw == "HER2_2" ~ "Equivocal (2+)",
      HER2_raw == "HER2_3" ~ "Positive (3+)",
      TRUE ~ NA_character_
    ),

    # 以原始临床数据KI67中位数动态分组：低于中位数为Low，中位数及以上为High。
    Ki67_cutoff = ki67_cutoff,
    Ki67_group = case_when(
      is.na(Ki67_percent) ~ NA_character_,
      Ki67_percent < ki67_cutoff ~ "Low",
      Ki67_percent >= ki67_cutoff ~ "High"
    ),
    Age_group = clean_missing(Age),
    Gender_group = clean_missing(Gender),

    # M分期只保留M0和M1，分别表示无远处转移和存在远处转移。
    M_raw = clean_missing(M),
    M_group = case_when(
      M_raw == "M0" ~ "M0",
      M_raw == "M1" ~ "M1",
      TRUE ~ NA_character_
    )
  )

# 每个子图的配置集中保存在plot_specs中，包括：
#   var    ：clinical_clean中的绘图字段名；
#   title  ：子图标题；
#   levels ：类别顺序，控制图例、扇区和柱子的排列；
#   colors ：类别与颜色的一一对应关系，确保重复运行时颜色不发生变化。
plot_specs <- list(
  list(var = "Stage_group", title = "Stage", levels = c("I", "II", "III", "IV"),
       colors = c("I" = "#FDD49E", "II" = "#FDBF6F", "III" = "#EF3B2C", "IV" = "#99000D")),
  list(var = "Lauren_group", title = "Lauren classification",
       levels = c("Intestinal", "Diffuse", "Mixed"),
       colors = c("Intestinal" = "#4DBBD5", "Diffuse" = "#E64B35", "Mixed" = "#E7B64F")),
  list(var = "Differentiation_group", title = "Differentiation",
       levels = c("Well differentiated", "Moderately differentiated", "Poorly differentiated"),
       colors = c("Well differentiated" = "#91D1C2", "Moderately differentiated" = "#F0B75E", "Poorly differentiated" = "#DC0000")),
  list(var = "HER2_group", title = "HER2",
       levels = c("Negative (0/1+)", "Equivocal (2+)", "Positive (3+)"),
       colors = c("Negative (0/1+)" = "#6BAED6", "Equivocal (2+)" = "#F0B75E", "Positive (3+)" = "#DE2D26")),
  list(var = "Ki67_group", title = "Ki67", levels = c("Low", "High"),
       colors = c("Low" = "#6BAED6", "High" = "#DE2D26")),
  list(var = "Age_group", title = "Age", levels = c("<60", ">=60"),
       colors = c("<60" = "#91D1C2", ">=60" = "#DC0000")),
  list(var = "Gender_group", title = "Sex", levels = c("Male", "Female"),
       colors = c("Male" = "#4DBBD5", "Female" = "#E64B35")),
  list(var = "M_group", title = "Distant metastasis", levels = c("M0", "M1"),
       colors = c("M0" = "#6BAED6", "M1" = "#DE2D26"))
)

# 对单个临床指标进行计数并计算百分比。
# factor(levels=...)用于固定类别顺序；.drop=TRUE仅保留数据中实际出现的类别。
# valid_n是该指标排除缺失值后的有效样本数，也是百分比计算的分母。
summarize_variable <- function(data, spec) {
  values <- data[[spec$var]]
  values <- factor(values, levels = spec$levels)
  tibble(group = values) %>%
    filter(!is.na(group)) %>%
    count(group, name = "n", .drop = TRUE) %>%
    mutate(
      variable = spec$var,
      title = spec$title,
      percent = n / sum(n) * 100,
      valid_n = sum(n),
      legend_label = paste0(as.character(group), " (n=", n, ", ", sprintf("%.1f", percent), "%)")
    ) %>%
    select(variable, title, group, n, percent, valid_n, legend_label)
}

# 依次汇总8个指标，并同时生成：
#   summary_list：按指标拆分的列表，供单独绘图使用；
#   figure1_summary：合并后的长表，输出为CSV供人工核查或后续制表。
summary_list <- lapply(plot_specs, function(x) summarize_variable(clinical_clean, x))
names(summary_list) <- vapply(plot_specs, `[[`, character(1), "var")
figure1_summary <- bind_rows(summary_list)

# 通用饼图函数。
# 输入为某个指标的汇总表和对应配置，输出为一个ggplot对象。
# 饼图内部不放文字，样本量和百分比统一显示在右侧图例中，避免小扇区标签重叠。
plot_pie <- function(summary_df, spec) {
  if (nrow(summary_df) == 0) stop("No valid observations for ", spec$var)
  summary_df$group <- factor(summary_df$group, levels = spec$levels)
  label_map <- setNames(summary_df$legend_label, as.character(summary_df$group))
  ggplot(summary_df, aes(x = "", y = n, fill = group)) +
    # 先绘制堆叠柱，再通过极坐标转换成饼图；白色边线用于区分相邻扇区。
    geom_col(width = 1, color = "white", linewidth = 0.55) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = spec$colors, labels = label_map, drop = TRUE) +
    labs(title = spec$title, fill = NULL) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      legend.position = "right",
      legend.text = element_text(size = 7.8, color = "black"),
      legend.key.size = grid::unit(0.38, "cm"),
      plot.margin = margin(4, 4, 4, 4)
    )
}

# 通用竖向柱状图函数，用于分化程度、HER2和Ki67。
# 横轴为临床类别，纵轴为百分比；柱顶显示“n=样本数（百分比）”。
plot_bar <- function(summary_df, spec) {
  if (nrow(summary_df) == 0) stop("No valid observations for ", spec$var)
  summary_df$group <- factor(summary_df$group, levels = spec$levels)

  # 生成柱顶标签，例如“n=49 (77.8%)”。
  summary_df$bar_label <- paste0("n=", summary_df$n, " (", sprintf("%.1f", summary_df$percent), "%)")

  # 分化程度类别名称较长，因此横轴显示为Well、Moderately、Poorly。
  # HER2标签使用换行符，将检测结果和评分分成两行，防止相邻标签重叠。
  # Ki67等未单独指定的指标保持原始类别名称。
  axis_labels <- switch(
    spec$var,
    "Differentiation_group" = c(
      "Well differentiated" = "Well",
      "Moderately differentiated" = "Moderately",
      "Poorly differentiated" = "Poorly"
    ),
    "HER2_group" = c(
      "Negative (0/1+)" = "Negative\n(0/1+)",
      "Equivocal (2+)" = "Equivocal\n(2+)",
      "Positive (3+)" = "Positive\n(3+)"
    ),
    setNames(spec$levels, spec$levels)
  )
  ggplot(summary_df, aes(x = group, y = percent, fill = group)) +
    # 柱高使用百分比而非原始样本数，使不同有效样本量的指标仍可直观比较组成比例。
    geom_col(width = 0.68, color = "white", linewidth = 0.4) +
    # vjust为负数时文字位于柱顶上方；黑色标签在浅色和深色柱子上均清晰可见。
    geom_text(aes(label = bar_label), vjust = -0.35, size = 2.45, color = "black") +
    scale_fill_manual(values = spec$colors, drop = TRUE) +
    scale_x_discrete(labels = axis_labels) +
    scale_y_continuous(
      # 纵轴上限额外增加22%，为柱顶的n和百分比标签预留空间。
      limits = c(0, max(summary_df$percent) * 1.22),
      breaks = scales::pretty_breaks(n = 4),
      labels = function(x) paste0(x, "%"),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(title = spec$title, x = NULL, y = "Percentage") +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      axis.text.x = element_text(size = 7.6, color = "black", lineheight = 0.9),
      axis.text.y = element_text(size = 8, color = "black"),
      axis.title.y = element_text(size = 9),
      legend.position = "none",
      plot.margin = margin(8, 12, 6, 6)
    )
}

# 指定需要使用柱状图的三个变量；其余变量自动使用饼图。
# 如需增加柱状图，只需把对应的var字段名加入此向量。
bar_variables <- c("Differentiation_group", "HER2_group", "Ki67_group")

# Map按相同位置同时遍历plot_specs和summary_list。
# 根据变量是否出现在bar_variables中，自动选择柱状图或饼图函数。
plots <- Map(
  function(spec, summary_df) {
    if (spec$var %in% bar_variables) plot_bar(summary_df, spec) else plot_pie(summary_df, spec)
  },
  plot_specs,
  summary_list
)

# 将8张子图固定排列为2行×4列，并自动添加A–H面板标签。
# wrap_plots显式指定ncol和nrow，避免“+”和“/”运算符优先级导致布局改变。
figure1_gastric <-
  patchwork::wrap_plots(plots, ncol = 4, nrow = 2, byrow = TRUE) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 15))

# 保存中间数据：
#   clinical_clean：逐样本清洗及分组结果；
#   figure1_summary：逐指标计数、有效样本数和百分比；
#   RData：一次性保存原始表、清洗表、汇总表及绘图配置；
#   RDS：保存最终patchwork绘图对象，便于以后在R中直接载入修改。
write.csv(clinical_clean, file.path(data_dir, "Figure1_gastric_clinical_clean.csv"), row.names = FALSE, na = "")
write.csv(figure1_summary, file.path(data_dir, "Figure1_gastric_summary.csv"), row.names = FALSE, na = "")
save(clinical_raw, clinical_clean, figure1_summary, plot_specs,
     file = file.path(data_dir, "Figure1_gastric_data.RData"))
saveRDS(figure1_gastric, file.path(data_dir, "Figure1_gastric_plot.rds"))

# 同时输出三种图片格式：
#   PNG：300 dpi位图，适合预览、Word和常规投稿；
#   PDF/SVG：矢量图，放大不失真，适合排版软件进一步编辑。
ggsave(file.path(out_dir, "Figure1_gastric.png"), figure1_gastric,
       width = 16, height = 8.5, units = "in", dpi = 300, bg = "white")
ggsave(file.path(out_dir, "Figure1_gastric.pdf"), figure1_gastric,
       width = 16, height = 8.5, units = "in", device = cairo_pdf, bg = "white")
ggsave(file.path(out_dir, "Figure1_gastric.svg"), figure1_gastric,
       width = 16, height = 8.5, units = "in", bg = "white")

message("Completed Figure 1 gastric cancer plots.")
message("Input rows: ", nrow(clinical_raw))
message("Output directory: ", out_dir)

############################################################
## Part 2: Primary tumor location
############################################################
data_dir <- file.path(script_dir, "data")
out_dir <- root_output_dir
cohort_file <- file.path(data_dir, "gastric_all_samples_clinical_unknown_filled.txt")
clinical_xlsx <- file.path(data_dir, "临床信息.xlsx")
if (!file.exists(cohort_file)) stop("Input file not found: ", cohort_file)
if (!file.exists(clinical_xlsx)) stop("Clinical Excel file not found: ", clinical_xlsx)
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

############################################################
## Part 3: Mean VAF and stage/Lauren subtype plots
############################################################
base_dir <- file.path(script_dir, "data")
out_dir <- file.path(root_output_dir, "vaf_stage_subtype")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
vaf_file <- file.path(base_dir, "sample_VAF_statistics.txt")
clinical_file <- file.path(base_dir, "gastric_all_samples_clinical_unknown_filled.txt")
cat("VAF文件检索路径：", vaf_file, "\n")
cat("临床文件检索路径：", clinical_file, "\n")
if (!file.exists(vaf_file)) stop("VAF input file not found: ", vaf_file)
if (!file.exists(clinical_file)) stop("Clinical input file not found: ", clinical_file)
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
