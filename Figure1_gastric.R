############################################################
## Figure 1: Gastric cancer clinical characteristics
## Stage / Lauren / Differentiation / HER2 / Ki67 / Age / Sex / M
## Input: ../gastric_txt/gastric_all_samples_clinical_unknown_filled.txt
##        ../gastric_txt/临床信息.xlsx（KI67原始值）
## Outputs: data/ (intermediate data) and output/ (figures)
##
## 图形类型：
##   1. Stage、Lauren、Age、Sex、M：饼图
##   2. Differentiation、HER2、Ki67：竖向百分比柱状图
##
## 统计口径：
##   - 原始表中的全部样本都会保存在清洗后数据中；
##   - 每个临床指标分别排除该指标的 Unknown、空字符串等无效值；
##   - 百分比的分母是该指标的有效样本数，而不是原始数据总行数；
##   - 图中 n 表示该类别的样本数，括号内为该类别在有效样本中的比例。
############################################################

# 加载数据整理、绘图和拼图所需的R包。
# suppressPackageStartupMessages用于隐藏载入包时的普通提示，便于查看真正的报错和运行结果。
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readxl)
})

get_script_dir <- function() {
  # Rscript运行时，commandArgs中通常包含“--file=脚本路径”。
  # 通过该参数定位脚本自身所在目录，使脚本不依赖当前工作目录，可以从任意位置运行。
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  normalizePath(getwd())
}

script_dir <- get_script_dir()

# 输入文件位于“代码”目录的同级目录“gastric_txt”中。
# 使用相对脚本位置构造路径，移动整个“胃癌”文件夹后仍可正常运行。
source_file <- file.path(
  dirname(script_dir), "gastric_txt",
  "gastric_all_samples_clinical_unknown_filled.txt"
)
raw_clinical_file <- file.path(dirname(script_dir), "gastric_txt", "临床信息.xlsx")
data_dir <- file.path(script_dir, "data")
out_dir <- file.path(script_dir, "output")

# recursive=TRUE表示父目录不存在时一并创建；showWarnings=FALSE避免目录已存在时产生警告。
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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
