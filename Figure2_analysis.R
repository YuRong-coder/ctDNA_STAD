############################################################
## Figure 2 unified analysis pipeline (panels A-D)
##
## Put all input files in ./data:
##   - Gastric_all_sample_multianno.txt
##   - oncoplot_matrix_CFDNA_Gastric.txt
##   - oncoplot_matrix_TCGA_STAD.txt
##   - patient_info.txt
##
## Run with: Rscript Figure2_analysis.R
## Outputs are written to ./output/Figure2.
############################################################

required_packages <- c(
  "data.table", "dplyr", "maftools", "tidyr", "stringr", "ggplot2",
  "GenomicRanges", "GenomeInfoDb", "IRanges", "BSgenome",
  "BSgenome.Hsapiens.UCSC.hg38", "MutationalPatterns"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "缺少R包：", paste(missing_packages, collapse = ", "),
    "。请先安装这些 CRAN/Bioconductor 依赖后再运行。"
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(maftools)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(IRanges)
  library(BSgenome)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(MutationalPatterns)
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
input_dir <- file.path(script_dir, "data")
figure2_output_dir <- file.path(script_dir, "output", "Figure2")
dir.create(figure2_output_dir, recursive = TRUE, showWarnings = FALSE)

required_input_files <- c(
  "Gastric_all_sample_multianno.txt",
  "oncoplot_matrix_CFDNA_Gastric.txt",
  "oncoplot_matrix_TCGA_STAD.txt",
  "patient_info.txt"
)
missing_input_files <- required_input_files[
  !file.exists(file.path(input_dir, required_input_files))
]
if (length(missing_input_files) > 0L) {
  stop("data/ 中缺少输入文件：", paste(missing_input_files, collapse = ", "))
}

############################################################
## Panel A: Mutation waterfall / oncoplot
############################################################
run_figure2a <- function(input_dir, output_dir) {
  outdir <- output_dir
## ===================== 1. 读取annovar原始突变注释文件 =====================
# 定义原始突变文件完整路径
infile <- file.path(input_dir, "Gastric_all_sample_multianno.txt")

# fread快速读取制表符分隔txt，data.table=FALSE输出普通数据框，适配dplyr语法
mut <- fread(
  infile,
  sep = "\t",
  header = TRUE,
  data.table = FALSE
)

# 打印前6行，快速查看数据结构
head(mut)
# 打印全部列名，校验Chr、Gene.refGene、ExonicFunc.refGene等关键字段是否存在
colnames(mut)

## ===================== 2. 自定义函数：转换突变分类Variant_Classification =====================
# maftools要求固定标准突变分类名称，annovar原生注释为描述文本，需要一一映射
# 入参：func=Func.refGene(基因区域注释)、exonic_func=ExonicFunc.refGene(外显子突变功能注释)
# 出参：maftools标准突变分类字符串
get_variant_classification <- function(func, exonic_func) {

  # 统一转为字符型，避免因子匹配失败
  func <- as.character(func)
  exonic_func <- as.character(exonic_func)

  # 分优先级匹配，将annovar注释翻译成maftools规范标签
  case_when(
    # 外显子功能突变（临床关注的非同义突变）
    exonic_func == "nonsynonymous SNV" ~ "Missense_Mutation",       # 错义突变
    exonic_func == "synonymous SNV" ~ "Silent",                     # 同义沉默突变
    exonic_func == "stopgain" ~ "Nonsense_Mutation",                # 无义突变，提前产生终止密码子
    exonic_func == "stoploss" ~ "Nonstop_Mutation",                 # 终止密码子丢失
    exonic_func == "frameshift deletion" ~ "Frame_Shift_Del",       # 移码缺失
    exonic_func == "frameshift insertion" ~ "Frame_Shift_Ins",      # 移码插入
    exonic_func == "nonframeshift deletion" ~ "In_Frame_Del",       # 框内缺失，阅读框不变
    exonic_func == "nonframeshift insertion" ~ "In_Frame_Ins",      # 框内插入，阅读框不变

    # 非编码区/剪接区域突变
    func == "splicing" ~ "Splice_Site",                             # 剪接位点突变
    func == "exonic" ~ "Missense_Mutation",                          # 普通外显子变异默认归为错义
    func == "intronic" ~ "Intron",                                  # 内含子突变
    func == "intergenic" ~ "IGR",                                   # 基因间区
    func == "upstream" ~ "5'Flank",                                 # 基因5'上游区域
    func == "downstream" ~ "3'Flank",                               # 基因3'下游区域
    func == "UTR3" ~ "3'UTR",                                       # 3'非翻译区
    func == "UTR5" ~ "5'UTR",                                       # 5'非翻译区
    grepl("ncRNA", func) ~ "RNA",                                   # 非编码RNA突变

    # 未匹配到的全部默认归类为基因间区
    TRUE ~ "IGR"
  )
}

## ===================== 3. 自定义函数：生成必填字段Variant_Type =====================
# maftools必填，区分单碱基替换、插入、缺失、多碱基替换等突变类型
# 判断逻辑：根据参考碱基Ref、突变碱基Alt的字符长度区分
get_variant_type <- function(ref, alt) {

  # 计算参考、突变碱基的字符长度
  ref_len <- nchar(ref)
  alt_len <- nchar(alt)

  case_when(
    ref_len == 1 & alt_len == 1 ~ "SNP",    # 单碱基替换
    ref_len < alt_len ~ "INS",              # 插入突变（突变序列更长）
    ref_len > alt_len ~ "DEL",              # 缺失突变（参考序列更长）
    ref_len == 2 & alt_len == 2 ~ "DNP",    # 双碱基替换
    ref_len == 3 & alt_len == 3 ~ "TNP",    # 三碱基替换
    ref_len > 3 & ref_len == alt_len ~ "ONP",# 大于3bp多碱基替换
    TRUE ~ "SNP"                            # 异常情况兜底归为SNP
  )
}

## ===================== 4. 标准化处理，生成完整标准MAF表格 =====================
# 基于原始mut表衍生所有maftools强制要求字段，输出规范MAF数据框
maf_df <- mut %>%
  mutate(
    # 处理基因名：多基因以;分隔只保留第一个；空值/点填充为Unknown
    Hugo_Symbol = ifelse(
      is.na(Gene.refGene) | Gene.refGene == "." | Gene.refGene == "",
      "Unknown",
      sapply(strsplit(Gene.refGene, ";"), `[`, 1)
    ),

    # 染色体标准化，去除chr前缀（chr1 → 1，适配hg38标准）
    Chromosome = gsub("^chr", "", Chr),
    # 突变起始坐标，直接映射原始Start列
    Start_Position = Start,
    # 突变终止坐标，直接映射原始End列
    End_Position = End,
    # 标准参考碱基列
    Reference_Allele = Ref,
    # Tumor_Seq_Allele1固定为参考碱基
    Tumor_Seq_Allele1 = Ref,
    # Tumor_Seq_Allele2固定为突变碱基Alt
    Tumor_Seq_Allele2 = Alt,

    # 调用自定义函数生成标准突变分类
    Variant_Classification = get_variant_classification(
      Func.refGene,
      ExonicFunc.refGene
    ),

    # 调用自定义函数生成突变类型SNP/INS/DEL
    Variant_Type = get_variant_type(Ref, Alt),

    # 样本ID列原生命名符合标准，直接保留
    Tumor_Sample_Barcode = Tumor_Sample_Barcode,
    # 参考基因组版本固定hg38
    NCBI_Build = "hg38",
    # 自定义队列标识
    Center = "Gastric"
  ) %>%
  # 筛选输出列：MAF标准必填字段 + 原始annovar注释字段用于溯源
  dplyr::select(
    Hugo_Symbol,
    Chromosome,
    Start_Position,
    End_Position,
    Reference_Allele,
    Tumor_Seq_Allele1,
    Tumor_Seq_Allele2,
    Variant_Classification,
    Variant_Type,
    Tumor_Sample_Barcode,
    NCBI_Build,
    Center,
    Func.refGene,
    Gene.refGene,
    GeneDetail.refGene,
    ExonicFunc.refGene,
    AAChange.refGene
  )

## ===================== 5. 导出标准.maf文本文件 =====================
# 定义输出MAF文件保存路径
outfile <- file.path(outdir, "Gastric_mutations.maf")

# 写出制表符分隔MAF文件，无引号、无行号，通用标准格式
write.table(
  maf_df,
  file = outfile,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

## ===================== 6. maftools读取MAF，生成分析对象 =====================
## ===================== 6. 读取并整理胃癌临床信息 =====================
# maftools临床注释要求第一列为Tumor_Sample_Barcode，且样本ID需要与MAF一致
clinical_file <- file.path(input_dir, "patient_info.txt")

clinical_raw <- fread(
  clinical_file,
  sep = "\t",
  header = TRUE,
  data.table = FALSE,
  check.names = FALSE
)

clean_value <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "-", "NA", "N/A", "na", "null", "NULL")] <- NA
  x
}

clinical_df <- clinical_raw %>%
  transmute(
    Tumor_Sample_Barcode = clean_value(research),
    Patient_ID = clean_value(`Patient ID`),
    Age = case_when(
      suppressWarnings(as.numeric(age)) < 60 ~ "<60",
      suppressWarnings(as.numeric(age)) >= 60 ~ ">=60",
      TRUE ~ NA_character_
    ),
    Gender = recode(clean_value(gender), "m" = "Male", "f" = "Female", .default = clean_value(gender)),
    pT = ifelse(is.na(clean_value(pT)), NA, paste0("T", clean_value(pT))),
    pN = ifelse(is.na(clean_value(pN)), NA, paste0("N", clean_value(pN))),
    M = ifelse(is.na(clean_value(M)), NA, paste0("M", clean_value(M))),
    Stage = clean_value(Stage),
    Metastasis = clean_value(`Metastasis type`),
    Lauren = recode(
      clean_value(`lauren type`),
      "1" = "Intestinal",
      "2" = "Diffuse",
      "3" = "Mixed",
      .default = clean_value(`lauren type`)
    ),
    Differentiation = clean_value(`Differentiation degree`),
    HER2 = ifelse(is.na(clean_value(HER2)), NA, paste0("HER2_", clean_value(HER2))),
    Ki67 = case_when(
      is.na(clean_value(KI67)) ~ NA_character_,
      grepl("^[0-9]+", clean_value(KI67)) &
        suppressWarnings(as.numeric(sub("^([0-9]+).*", "\\1", clean_value(KI67)))) >= 70 ~ "High",
      grepl("^[0-9]+", clean_value(KI67)) ~ "Low",
      TRUE ~ clean_value(KI67)
    )
  ) %>%
  filter(!is.na(Tumor_Sample_Barcode)) %>%
  distinct(Tumor_Sample_Barcode, .keep_all = TRUE)

maf_samples <- unique(as.character(maf_df$Tumor_Sample_Barcode))
clinical_samples <- unique(as.character(clinical_df$Tumor_Sample_Barcode))

# 本图的样本总体以 MAF 中的 127 个样本为准；临床表中不属于该总体的
# 记录不加入绘图。缺少临床信息的 MAF 样本稍后补为 Unknown。
clinical_df <- clinical_df %>%
  filter(Tumor_Sample_Barcode %in% maf_samples)

matched_sample_count <- nrow(clinical_df)
missing_samples <- setdiff(maf_samples, clinical_df$Tumor_Sample_Barcode)

if (length(missing_samples) > 0) {
  unknown_clinical_df <- data.frame(
    Tumor_Sample_Barcode = missing_samples,
    Patient_ID = "Unknown",
    Age = "Unknown",
    Gender = "Unknown",
    pT = "Unknown",
    pN = "Unknown",
    M = "Unknown",
    Stage = "Unknown",
    Metastasis = "Unknown",
    Lauren = "Unknown",
    Differentiation = "Unknown",
    HER2 = "Unknown",
    Ki67 = "Unknown",
    stringsAsFactors = FALSE
  )

  clinical_df <- bind_rows(clinical_df, unknown_clinical_df)
}

clinical_df <- clinical_df %>%
  mutate(
    across(
      everything(),
      ~ ifelse(is.na(.) | trimws(as.character(.)) == "", "Unknown", as.character(.))
    )
  )

all_samples <- clinical_df$Tumor_Sample_Barcode

matched_clinical_file <- file.path(outdir, "gastric_matched_clinical_info.txt")
write.table(
  clinical_df,
  file = matched_clinical_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nMAF样本数:", length(maf_samples), "\n")
cat("最终展示样本数:", length(all_samples), "\n")
cat("匹配到临床信息的样本数:", matched_sample_count, "\n")
cat("设为Unknown的MAF样本数:", length(missing_samples), "\n")

# 加载标准MAF文件，verbose=TRUE打印运行日志方便排错
gastric_maf <- read.maf(
  maf = outfile,
  clinicalData = clinical_df,
  verbose = TRUE
)

gastric_maf_rds_file <- file.path(outdir, "gastric_maf_with_clinical.rds")
saveRDS(gastric_maf, gastric_maf_rds_file)

gastric_maf_clinical_file <- file.path(outdir, "gastric_maf_clinical_data_in_object.txt")
write.table(
  gastric_maf@clinical.data,
  file = gastric_maf_clinical_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# 打印MAF对象基础统计信息（样本总数、突变总数、各类型突变占比等）
gastric_maf

## ===================== 7. 绘制突变瀑布图oncoplot =====================
# 查看MAF对象内临床信息和各临床分组，便于确认读入是否正常
head(gastric_maf@clinical.data)

for (f in c("Age", "Gender", "pT", "pN", "M", "Stage", "Metastasis", "Lauren", "Differentiation", "HER2", "Ki67")) {
  cat("\n临床变量:", f, "\n")
  print(table(gastric_maf@clinical.data[[f]], useNA = "ifany"))
}

variant_colors <- c(
  Missense_Mutation = "#43A2CA",
  Nonsense_Mutation = "#F8D042",
  Frame_Shift_Del = "#377EB8",
  Frame_Shift_Ins = "#C2E2C2",
  In_Frame_Del = "#F39B7F",
  In_Frame_Ins = "#8491B4",
  Splice_Site = "#B0D295",
  Translation_Start_Site = "#F0A0A0",
  Nonstop_Mutation = "#7E6148",
  Multi_Hit = "#B4DCEC"
)

ann_colors <- list(
  Age = c("<60" = "#2A9D8F", ">=60" = "#E9C46A"),
  Gender = c("Male" = "#2EBCB2", "Female" = "#E76F51"),
  pT = c("T0" = "#F5FFFA", "T1" = "#D9F0D3", "T1a" = "#C7E9C0", "T1b" = "#A1D99B", "T2" = "#74C476", "T3" = "#31A354", "T4" = "#006D2C", "T4a" = "#005A32", "T4b" = "#00441B"),
  pN = c("N0" = "#D8DEF7", "N1" = "#BCBDDC", "N2" = "#807DBA", "N3" = "#54278F", "N3a" = "#3F007D", "N3b" = "#2D004B"),
  M = c("M0" = "#8CC6ED", "M1" = "#F0A0A0"),
  Stage = c("0" = "#F5F5DC", "Ia" = "#FFF3B0", "Ib" = "#FAD643", "IIa" = "#F6BD60", "IIb" = "#F4A261", "IIIa" = "#E76F51", "IIIb" = "#D1495B", "IIIc" = "#9D0208", "IV" = "#6A040F"),
  Metastasis = c("CY" = "#FFB000", "LYM" = "#6A994E", "HEP" = "#A7C957", "PER" = "#EF476F", "PUL" = "#8338EC", "LYM/CY" = "#F77F00", "LYM/H" = "#BC6C25", "LYM/P" = "#D62828"),
  Lauren = c("Intestinal" = "#66C2A5", "Diffuse" = "#FC8D62", "Mixed" = "#8DA0CB"),
  HER2 = c("HER2_0" = "#F2E8CF", "HER2_1" = "#F6BD60", "HER2_2" = "#F28482", "HER2_3" = "#C1121F"),
  Ki67 = c("Low" = "#90BE6D", "High" = "#D45659")
)

ann_colors <- lapply(ann_colors, function(x) c(x, "Unknown" = "#ECF0F1"))

clinicalFeatures <- c("Age", "Gender", "pT", "pN", "M", "Stage", "Metastasis", "Lauren", "HER2", "Ki67")

clinicalFeatures_use <- clinicalFeatures[
  sapply(clinicalFeatures, function(x) {
    vals <- gastric_maf@clinical.data[[x]]
    vals <- vals[!is.na(vals) & vals != ""]
    length(unique(vals)) > 0
  })
]

ann_colors_use <- ann_colors[clinicalFeatures_use]

# Interleave samples without clinical information among clinically annotated
# samples instead of placing all Unknown samples together at one end.
known_samples <- clinical_df$Tumor_Sample_Barcode[
  clinical_df$Patient_ID != "Unknown"
]
unknown_samples <- clinical_df$Tumor_Sample_Barcode[
  clinical_df$Patient_ID == "Unknown"
]

# Shuffle both groups reproducibly before interleaving. This keeps Unknown
# samples dispersed even when oncoplot displays only a mutation-based subset.
set.seed(20260827)
known_samples <- sample(known_samples, length(known_samples))
unknown_samples <- sample(unknown_samples, length(unknown_samples))

interleave_samples <- function(known, unknown) {
  if (length(unknown) == 0) return(known)
  insert_after <- floor(seq(0, length(known), length.out = length(unknown) + 2))[-c(1, length(unknown) + 2)]
  result <- character(0)
  unknown_index <- 1L
  for (i in 0:length(known)) {
    if (i > 0) result <- c(result, known[i])
    while (unknown_index <= length(unknown) && insert_after[unknown_index] == i) {
      result <- c(result, unknown[unknown_index])
      unknown_index <- unknown_index + 1L
    }
  }
  if (unknown_index <= length(unknown)) {
    result <- c(result, unknown[unknown_index:length(unknown)])
  }
  result
}

sample_order_interleaved <- interleave_samples(known_samples, unknown_samples)

# 输出样本审计表，便于核对每个展示样本是否有突变及临床记录。
sample_display_audit <- data.frame(
  Tumor_Sample_Barcode = sample_order_interleaved,
  Has_mutation = sample_order_interleaved %in% maf_samples,
  Has_clinical_record = sample_order_interleaved %in% clinical_samples,
  stringsAsFactors = FALSE
)
write.table(
  sample_display_audit,
  file = file.path(outdir, "all_samples_display_audit.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

plot_file <- file.path(
  outdir,
  paste0("gastric_maf_clinical_oncoplot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
)

pdf(
  plot_file,
  width = 18,
  height = 12,
  useDingbats = FALSE
)

# Use one size for the main plot, clinical annotations, and all legends.
uniform_legend_size <- 0.7

# 展示队列突变分布，可视化各样本高频突变基因
oncoplot(
  maf = gastric_maf,
  top = 30,
  clinicalFeatures = clinicalFeatures_use,
  colors = variant_colors,
  annotationColor = ann_colors_use,
  sortByAnnotation = FALSE,
  sortByMutation = FALSE,
  sampleOrder = sample_order_interleaved,
  # 必须设为 FALSE，保留在前 30 个基因中没有突变的 MAF 样本。
  removeNonMutated = FALSE,
  showTumorSampleBarcodes = FALSE,
  drawRowBar = TRUE,
  drawColBar = TRUE,
  legend_height = 10,
  fontSize = uniform_legend_size,
  legendFontSize = uniform_legend_size,
  annotationFontSize = uniform_legend_size,
  titleText = "Clinical annotated gastric cancer mutation landscape"
)

dev.off()
}

############################################################
## Panel B: Stage-specific top-five mutation dot plot
############################################################
run_figure2b <- function(input_dir, output_dir) {
  matrix_file <- file.path(input_dir, "oncoplot_matrix_CFDNA_Gastric.txt")
  multianno_file <- file.path(input_dir, "Gastric_all_sample_multianno.txt")
  clinical_file <- file.path(input_dir, "patient_info.txt")
  outdir <- output_dir
# ====================== 路径打印 + 容错文件检测 ======================
cat("==== 文件路径校验 ====\n")
cat("突变矩阵：", matrix_file, "\n")
cat("Annovar注释：", multianno_file, "\n")
cat("临床信息表：", clinical_file, "\n\n")

check_file <- function(path){
  if(!file.exists(path)){
    stop(paste0("错误：未找到文件，请检查路径\n", path))
  }
}
check_file(matrix_file)
check_file(multianno_file)
check_file(clinical_file)

# ===================== 清洗工具函数 =====================
clean_blank <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "-", "/", "NA", "N/A", "na", "n/a")] <- NA
  x
}
norm_name <- function(x) tolower(gsub("[^A-Za-z0-9]", "", as.character(x)))
stage_clean <- function(x) {
  x <- clean_blank(x)
  x <- toupper(x)
  x <- gsub("^([IV]+)[A-C]?$", "\\1", x)
  x
}

variant_class_from_annovar <- function(func, exonic_func) {
  func <- as.character(func)
  exonic_func <- as.character(exonic_func)
  dplyr::case_when(
    grepl("splicing", func, ignore.case = TRUE) ~ "Splice_Site",
    grepl("stopgain", exonic_func, ignore.case = TRUE) ~ "Nonsense_Mutation",
    grepl("stoploss", exonic_func, ignore.case = TRUE) ~ "Nonstop_Mutation",
    grepl("frameshift deletion", exonic_func, ignore.case = TRUE) ~ "Frame_Shift_Del",
    grepl("frameshift insertion", exonic_func, ignore.case = TRUE) ~ "Frame_Shift_Ins",
    grepl("nonframeshift deletion", exonic_func, ignore.case = TRUE) ~ "In_Frame_Del",
    grepl("nonframeshift insertion", exonic_func, ignore.case = TRUE) ~ "In_Frame_Ins",
    grepl("nonsynonymous SNV", exonic_func, ignore.case = TRUE) ~ "Missense_Mutation",
    grepl("synonymous SNV", exonic_func, ignore.case = TRUE) ~ "Silent",
    TRUE ~ "Missense_Mutation"
  )
}
variant_type_from_ref_alt <- function(ref, alt) {
  ref <- as.character(ref)
  alt <- as.character(alt)
  dplyr::case_when(
    nchar(ref) == 1 & nchar(alt) == 1 ~ "SNP",
    nchar(ref) > nchar(alt) ~ "DEL",
    nchar(ref) < nchar(alt) ~ "INS",
    TRUE ~ "ONP"
  )
}
first_gene <- function(x) {
  x <- as.character(x)
  x <- sub(";.*$", "", x)
  x[x %in% c("", ".")] <- NA
  toupper(x)
}

# ===================== 0. 读取突变矩阵+Annovar注释，生成标准MAF =====================
mut_matrix <- read.delim(
  matrix_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  row.names = 1
)
mut_matrix[] <- lapply(mut_matrix, function(x) as.numeric(as.character(x)))
mut_matrix[is.na(mut_matrix)] <- 0
mut_matrix <- as.matrix(mut_matrix)
sample_ids <- rownames(mut_matrix)
gene_ids <- colnames(mut_matrix)

multianno_raw <- read.delim(
  multianno_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# 构建MAF长表
mut_long <- multianno_raw %>%
  dplyr::transmute(
    Hugo_Symbol = first_gene(.data[["Gene.refGene"]]),
    Chromosome = as.character(.data[["Chr"]]),
    Start_Position = suppressWarnings(as.integer(.data[["Start"]])),
    End_Position = suppressWarnings(as.integer(.data[["End"]])),
    Variant_Classification = variant_class_from_annovar(.data[["Func.refGene"]], .data[["ExonicFunc.refGene"]]),
    Variant_Type = variant_type_from_ref_alt(.data[["Ref"]], .data[["Alt"]]),
    Reference_Allele = toupper(as.character(.data[["Ref"]])),
    Tumor_Seq_Allele1 = toupper(as.character(.data[["Ref"]])),
    Tumor_Seq_Allele2 = toupper(as.character(.data[["Alt"]])),
    Tumor_Sample_Barcode = as.character(.data[["Tumor_Sample_Barcode"]]),
    Data_Source = as.character(.data[["Data_Source"]])
  ) %>%
  dplyr::filter(
    Tumor_Sample_Barcode %in% sample_ids,
    !is.na(Hugo_Symbol),
    Hugo_Symbol != "",
    Hugo_Symbol != "NONE",
    Hugo_Symbol %in% toupper(gene_ids),
    Variant_Classification != "Silent",
    !is.na(Start_Position),
    !is.na(End_Position),
    Reference_Allele %in% c("A", "T", "C", "G"),
    Tumor_Seq_Allele2 %in% c("A", "T", "C", "G"),
    Reference_Allele != Tumor_Seq_Allele2
  ) %>%
  dplyr::distinct(
    Hugo_Symbol,
    Chromosome,
    Start_Position,
    End_Position,
    Reference_Allele,
    Tumor_Seq_Allele2,
    Tumor_Sample_Barcode,
    .keep_all = TRUE
  ) %>%
  dplyr::select(-Data_Source)

if (nrow(mut_long) == 0) stop("过滤后无有效突变数据，请核对样本ID匹配情况！")

# ===================== 1. 临床数据清洗匹配（修复norm笔误为norm_name） =====================
clinical_raw <- read.delim(
  clinical_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
names(clinical_raw) <- trimws(names(clinical_raw))
clinical_raw <- clinical_raw[, !is.na(names(clinical_raw)) & names(clinical_raw) != "", drop = FALSE]
names(clinical_raw) <- make.unique(names(clinical_raw))

# 【修复处：norm → norm_name】
clinical_raw <- clinical_raw %>%
  dplyr::mutate(
    research = clean_blank(.data[["research"]]),
    Patient_ID = clean_blank(.data[["Patient ID"]]),
    sample_key = norm_name(research)
  )
sample_map <- data.frame(
  Tumor_Sample_Barcode = sample_ids,
  sample_key = norm_name(sample_ids),
  stringsAsFactors = FALSE
)
clinical_joined <- sample_map %>% dplyr::left_join(clinical_raw, by = "sample_key")

clinical_data <- clinical_joined %>%
  dplyr::transmute(
    Tumor_Sample_Barcode = Tumor_Sample_Barcode,
    cohort = "Gastric_cfDNA",
    Age = dplyr::case_when(
      suppressWarnings(as.numeric(age)) < 60 ~ "<60",
      suppressWarnings(as.numeric(age)) >= 60 ~ ">=60",
      TRUE ~ NA_character_
    ),
    Gender = dplyr::case_when(
      tolower(clean_blank(gender)) %in% c("m", "male") ~ "Male",
      tolower(clean_blank(gender)) %in% c("f", "female") ~ "Female",
      TRUE ~ NA_character_
    ),
    Stage = stage_clean(.data[["Stage"]])
  )

# ===================== 2. 构建maf核心对象 =====================
maf_sub <- maftools::read.maf(
  maf = mut_long,
  clinicalData = clinical_data,
  verbose = FALSE
)

# ===================== 【Figure2B 气泡图绘图逻辑】 =====================
mut_data <- maf_sub@data %>% as.data.frame()
clinical_data <- maf_sub@clinical.data %>% as.data.frame()

sample_col <- "Tumor_Sample_Barcode"
gene_col <- "Hugo_Symbol"
subtype_col <- "Stage"
exome_size_mb <- 38

# 1. 计算每个样本总突变数、TMB
sample_burden_df <- mut_data %>%
  dplyr::filter(
    !is.na(.data[[sample_col]]),
    !is.na(.data[[gene_col]]),
    .data[[gene_col]] != ""
  ) %>%
  dplyr::group_by(sample_id = .data[[sample_col]]) %>%
  dplyr::summarise(
    mutation_count = dplyr::n(),
    mutated_genes = dplyr::n_distinct(.data[[gene_col]]),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    estimated_TMB = mutation_count / exome_size_mb,
    mutation_burden = mutation_count
  )

# 2. 提取分期信息，过滤无分期样本
clinical_stage_df <- clinical_data %>%
  dplyr::select(
    sample_id = all_of(sample_col),
    stage_group = all_of(subtype_col)
  ) %>%
  dplyr::mutate(
    sample_id = as.character(sample_id),
    stage_group = as.character(stage_group),
    stage_group = ifelse(stage_group %in% c("", "/", "NA", "N/A"), NA, stage_group)
  ) %>%
  dplyr::filter(!is.na(stage_group))

# 3. 突变+分期+TMB合并
mut_stage_df <- mut_data %>%
  dplyr::select(
    sample_id = all_of(sample_col),
    Gene = all_of(gene_col)
  ) %>%
  dplyr::mutate(
    sample_id = as.character(sample_id),
    Gene = toupper(as.character(Gene))
  ) %>%
  dplyr::filter(
    !is.na(sample_id),
    !is.na(Gene),
    Gene != ""
  ) %>%
  dplyr::left_join(clinical_stage_df, by = "sample_id") %>%
  dplyr::left_join(sample_burden_df, by = "sample_id") %>%
  dplyr::filter(!is.na(stage_group))

# 4. 各分期总样本数
stage_n_df <- clinical_stage_df %>%
  dplyr::count(stage_group, name = "total_samples")

# 5. 统计各分期基因突变频率，每组取TOP5
top5_gene_stage_df <- mut_stage_df %>%
  dplyr::distinct(
    sample_id,
    stage_group,
    Gene,
    estimated_TMB,
    mutation_burden
  ) %>%
  dplyr::group_by(stage_group, Gene) %>%
  dplyr::summarise(
    mutated_samples = dplyr::n_distinct(sample_id),
    median_estimated_TMB = median(estimated_TMB, na.rm = TRUE),
    mean_estimated_TMB = mean(estimated_TMB, na.rm = TRUE),
    median_mutation_burden = median(mutation_burden, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(stage_n_df, by = "stage_group") %>%
  dplyr::mutate(
    mutation_frequency = mutated_samples / total_samples * 100
  ) %>%
  dplyr::filter(mutated_samples >= 1) %>%
  dplyr::group_by(stage_group) %>%
  dplyr::arrange(
    desc(mutation_frequency),
    desc(mutated_samples),
    desc(median_estimated_TMB),
    .by_group = TRUE
  ) %>%
  dplyr::slice_head(n = 5) %>%
  dplyr::ungroup()

# 输出统计表格
write.csv(
  top5_gene_stage_df,
  file.path(outdir, "Figure2B_gastric_stage_top5_dotplot_data.csv"),
  row.names = FALSE
)

# 6. 绘图坐标轴顺序设定
stage_order <- c("I", "II", "III", "IV")
stage_order <- stage_order[stage_order %in% unique(as.character(top5_gene_stage_df$stage_group))]
stage_order <- c(stage_order, setdiff(unique(as.character(top5_gene_stage_df$stage_group)), stage_order))

gene_order_by_stage <- top5_gene_stage_df %>%
  dplyr::mutate(
    stage_group = as.character(stage_group),
    Gene = as.character(Gene),
    stage_group = factor(stage_group, levels = stage_order)
  ) %>%
  dplyr::arrange(
    stage_group,
    desc(mutation_frequency),
    desc(mutated_samples),
    desc(median_estimated_TMB)
  ) %>%
  dplyr::pull(Gene) %>%
  unique()

top5_gene_stage_df_plot <- top5_gene_stage_df %>%
  dplyr::mutate(
    stage_group = factor(as.character(stage_group), levels = stage_order),
    Gene = factor(as.character(Gene), levels = rev(gene_order_by_stage))
  )

# 7. ggplot气泡图绘制
p_top5_gastric_dot <- ggplot(
  top5_gene_stage_df_plot,
  aes(
    x = stage_group,
    y = Gene
  )
) +
  geom_point(
    aes(
      size = mutation_frequency,
      color = median_estimated_TMB
    ),
    alpha = 0.9
  ) +
  scale_size_continuous(
    range = c(3, 10),
    name = "Mutation frequency (%)"
  ) +
  scale_color_gradient(
    low = "lightblue",
    high = "#D73027",
    name = "Median estimated TMB"
  ) +
  labs(
    x = "Gastric cancer stage",
    y = NULL,
    title = "Top 5 mutated genes across gastric cancer stages",
    subtitle = "Dot size = mutation frequency; Dot color = median estimated TMB"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey30"),
    axis.text.x = element_text(angle = 35, hjust = 1, color = "black", face = "bold", size = 11),
    axis.text.y = element_text(color = "black", face = "italic", size = 10),
    axis.line = element_line(linewidth = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    plot.margin = margin(10, 15, 10, 10)
  )

# 8. 导出PDF+高清PNG
ggsave(
  filename = file.path(outdir, "Figure2B_gastric_top5_by_stage_dotplot.pdf"),
  plot = p_top5_gastric_dot,
  width = 7,
  height = 5,
  useDingbats = FALSE
)
ggsave(
  filename = file.path(outdir, "Figure2B_gastric_top5_by_stage_dotplot.png"),
  plot = p_top5_gastric_dot,
  width = 7,
  height = 5,
  dpi = 300
)

cat("Figure2B 绘图完成！输出文件夹：", outdir, "\n")
print(p_top5_gastric_dot)

# ----------------------新增代码开始----------------------
# 统计绘制气泡图所使用的有效样本
total_plot_samples <- nrow(clinical_stage_df)
stage_sample_count <- clinical_stage_df %>% count(stage_group, name = "样本数量")
cat("\n================= 气泡图样本统计结果 =================")
cat("\n气泡图总有效样本量：", total_plot_samples, "\n")
cat("各胃癌分期样本分布：\n")
print(stage_sample_count, row.names = FALSE)
# ----------------------新增代码结束----------------------
}

############################################################
## Panel C: MAF summary
############################################################
run_figure2c <- function(input_dir, output_dir) {
  matrix_file <- file.path(input_dir, "oncoplot_matrix_CFDNA_Gastric.txt")
  multianno_file <- file.path(input_dir, "Gastric_all_sample_multianno.txt")
  clinical_file <- file.path(input_dir, "patient_info.txt")
  outdir <- output_dir
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
# 校验文件存在
stopifnot(file.exists(matrix_file))
stopifnot(file.exists(multianno_file))
stopifnot(file.exists(clinical_file))

# ==================== 读取突变矩阵，获取样本/基因列表 ====================
mut_matrix <- read.delim(
  matrix_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  row.names = 1
)
mut_matrix[] <- lapply(mut_matrix, function(x) as.numeric(as.character(x)))
mut_matrix[is.na(mut_matrix)] <- 0
mut_matrix <- as.matrix(mut_matrix)

sample_ids <- rownames(mut_matrix)
gene_ids <- colnames(mut_matrix)

cat("样本数量：", length(sample_ids), "\n基因数量：", length(gene_ids), "\n")

# ==================== 解析annovar注释，转为maftools标准MAF格式 ====================
multianno_raw <- read.delim(
  multianno_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# 突变类型转换函数
variant_class_from_annovar <- function(func, exonic_func) {
  func <- as.character(func)
  exonic_func <- as.character(exonic_func)
  dplyr::case_when(
    grepl("splicing", func, ignore.case = TRUE) ~ "Splice_Site",
    grepl("stopgain", exonic_func, ignore.case = TRUE) ~ "Nonsense_Mutation",
    grepl("stoploss", exonic_func, ignore.case = TRUE) ~ "Nonstop_Mutation",
    grepl("frameshift deletion", exonic_func, ignore.case = TRUE) ~ "Frame_Shift_Del",
    grepl("frameshift insertion", exonic_func, ignore.case = TRUE) ~ "Frame_Shift_Ins",
    grepl("nonframeshift deletion", exonic_func, ignore.case = TRUE) ~ "In_Frame_Del",
    grepl("nonframeshift insertion", exonic_func, ignore.case = TRUE) ~ "In_Frame_Ins",
    grepl("nonsynonymous SNV", exonic_func, ignore.case = TRUE) ~ "Missense_Mutation",
    grepl("synonymous SNV", exonic_func, ignore.case = TRUE) ~ "Silent",
    TRUE ~ "Missense_Mutation"
  )
}

variant_type_from_ref_alt <- function(ref, alt) {
  ref <- as.character(ref)
  alt <- as.character(alt)
  dplyr::case_when(
    nchar(ref) == 1 & nchar(alt) == 1 ~ "SNP",
    nchar(ref) > nchar(alt) ~ "DEL",
    nchar(ref) < nchar(alt) ~ "INS",
    TRUE ~ "ONP"
  )
}

first_gene <- function(x) {
  x <- as.character(x)
  x <- sub(";.*$", "", x)
  x[x %in% c("", ".")] <- NA
  toupper(x)
}

# 清洗过滤突变数据，生成标准MAF输入
mut_long <- multianno_raw %>%
  dplyr::transmute(
    Hugo_Symbol = first_gene(.data[["Gene.refGene"]]),
    Chromosome = as.character(.data[["Chr"]]),
    Start_Position = suppressWarnings(as.integer(.data[["Start"]])),
    End_Position = suppressWarnings(as.integer(.data[["End"]])),
    Variant_Classification = variant_class_from_annovar(.data[["Func.refGene"]], .data[["ExonicFunc.refGene"]]),
    Variant_Type = variant_type_from_ref_alt(.data[["Ref"]], .data[["Alt"]]),
    Reference_Allele = toupper(as.character(.data[["Ref"]])),
    Tumor_Seq_Allele1 = toupper(as.character(.data[["Ref"]])),
    Tumor_Seq_Allele2 = toupper(as.character(.data[["Alt"]])),
    Tumor_Sample_Barcode = as.character(.data[["Tumor_Sample_Barcode"]]),
    Data_Source = as.character(.data[["Data_Source"]])
  ) %>%
  dplyr::filter(
    Tumor_Sample_Barcode %in% sample_ids,
    !is.na(Hugo_Symbol),
    Hugo_Symbol != "",
    Hugo_Symbol != "NONE",
    Hugo_Symbol %in% toupper(gene_ids),
    Variant_Classification != "Silent",
    !is.na(Start_Position),
    !is.na(End_Position),
    Reference_Allele %in% c("A", "T", "C", "G"),
    Tumor_Seq_Allele2 %in% c("A", "T", "C", "G"),
    Reference_Allele != Tumor_Seq_Allele2
  ) %>%
  dplyr::distinct(
    Hugo_Symbol,
    Chromosome,
    Start_Position,
    End_Position,
    Reference_Allele,
    Tumor_Seq_Allele2,
    Tumor_Sample_Barcode,
    .keep_all = TRUE
  ) %>%
  dplyr::select(
    Hugo_Symbol,
    Chromosome,
    Start_Position,
    End_Position,
    Variant_Classification,
    Variant_Type,
    Reference_Allele,
    Tumor_Seq_Allele1,
    Tumor_Seq_Allele2,
    Tumor_Sample_Barcode,
    Data_Source
  )

if (nrow(mut_long) == 0) {
  stop("匹配样本后无有效突变数据，请核对注释文件与矩阵样本名！")
}

# 输出中间MAF表格
write.table(
  mut_long,
  file.path(outdir, "gastric_cfDNA_maf_from_multianno.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ==================== 临床数据清洗与样本匹配 ====================
clinical_raw <- read.delim(
  clinical_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
names(clinical_raw) <- trimws(names(clinical_raw))
clinical_raw <- clinical_raw[, !is.na(names(clinical_raw)) & names(clinical_raw) != "", drop = FALSE]
names(clinical_raw) <- make.unique(names(clinical_raw))

# 清洗工具函数
clean_blank <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "-", "/", "NA", "N/A", "na", "n/a")] <- NA
  x
}
norm_name <- function(x) {
  tolower(gsub("[^A-Za-z0-9]", "", as.character(x)))
}
stage_clean <- function(x) {
  x <- clean_blank(x)
  x <- toupper(x)
  x <- gsub("^([IV]+)[A-C]?$", "\\1", x)
  x
}
pt_clean <- function(x) {
  x <- clean_blank(x)
  x <- toupper(x)
  x <- sub("^([0-9]+).*", "T\\1", x)
  x <- ifelse(is.na(x), NA_character_, ifelse(grepl("^T", x), x, paste0("T", x)))
  x
}
pn_clean <- function(x) {
  x <- clean_blank(x)
  x <- toupper(x)
  x <- sub("^([0-9]+).*", "N\\1", x)
  x <- ifelse(is.na(x), NA_character_, ifelse(grepl("^N", x), x, paste0("N", x)))
}
pm_clean <- function(x) {
  x <- clean_blank(x)
  x <- toupper(x)
  x <- ifelse(is.na(x), NA_character_, ifelse(grepl("^M", x), x, paste0("M", x)))
  x
}
ihc_clean <- function(x) {
  x <- clean_blank(x)
  x <- ifelse(is.na(x), NA, x)
  x <- ifelse(x %in% c("0", "1", "1+", "2", "2+", "3", "3+"), x, x)
  x
}
ki67_group <- function(x) {
  x <- clean_blank(x)
  pct <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", x)))
  dplyr::case_when(
    is.na(pct) ~ NA_character_,
    pct < 50 ~ "<50%",
    pct >= 50 ~ ">=50%"
  )
}

# 临床数据预处理
clinical_raw <- clinical_raw %>%
  dplyr::mutate(
    research = clean_blank(.data[["research"]]),
    Patient_ID = clean_blank(.data[["Patient ID"]]),
    sample_key = norm_name(research)
  )

sample_map <- data.frame(
  Tumor_Sample_Barcode = sample_ids,
  sample_key = norm_name(sample_ids),
  stringsAsFactors = FALSE
)

clinical_joined <- sample_map %>%
  dplyr::left_join(clinical_raw, by = "sample_key")

clinical_data <- clinical_joined %>%
  dplyr::transmute(
    Tumor_Sample_Barcode = Tumor_Sample_Barcode,
    cohort = "Gastric_cfDNA",
    Age = dplyr::case_when(
      suppressWarnings(as.numeric(age)) < 60 ~ "<60",
      suppressWarnings(as.numeric(age)) >= 60 ~ ">=60",
      TRUE ~ NA_character_
    ),
    Gender = dplyr::case_when(
      tolower(clean_blank(gender)) %in% c("m", "male") ~ "Male",
      tolower(clean_blank(gender)) %in% c("f", "female") ~ "Female",
      TRUE ~ NA_character_
    ),
    pT = pt_clean(.data[["pT"]]),
    pN = pn_clean(.data[["pN"]]),
    M = pm_clean(.data[["M"]]),
    Stage = stage_clean(.data[["Stage"]]),
    Histology = clean_blank(.data[["histology"]]),
    Differentiation = clean_blank(.data[["Differentiation degree"]]),
    Lauren_type = dplyr::case_when(
      clean_blank(.data[["lauren type"]]) == "1" ~ "Intestinal",
      clean_blank(.data[["lauren type"]]) == "2" ~ "Diffuse",
      clean_blank(.data[["lauren type"]]) == "3" ~ "Mixed",
      clean_blank(.data[["lauren type"]]) == "0" ~ NA_character_,
      TRUE ~ clean_blank(.data[["lauren type"]])
    ),
    Tumor_site = clean_blank(.data[["Tumor location(site)"]]),
    HER2 = ihc_clean(.data[["HER2"]]),
    KI67 = ki67_group(.data[["KI67"]]),
    MLH1 = ihc_clean(.data[["MLH1"]]),
    MSH2 = ihc_clean(.data[["MSH2"]]),
    MSH6 = ihc_clean(.data[["MSH6"]]),
    PMS2 = ihc_clean(.data[["PMS2"]])
  )

# 输出匹配临床表
write.csv(
  clinical_data,
  file.path(outdir, "gastric_cfDNA_clinical_matched.csv"),
  row.names = FALSE
)

# ==================== 构建MAF对象 ====================
maf_sub <- maftools::read.maf(
  maf = mut_long %>% dplyr::select(-Data_Source),
  clinicalData = clinical_data,
  verbose = FALSE
)

# ==================== Figure2C：plotmafSummary 突变汇总图 ====================
pdf(
  file.path(outdir, "Figure2C_gastric_maf_summary.pdf"),
  width = 8,
  height = 5,
  useDingbats = FALSE
)

plotmafSummary(
  maf = maf_sub,
  rmOutlier = TRUE,
  addStat = "median",
  dashboard = TRUE,
  titvRaw = FALSE
)

dev.off()

cat("Figure2C 绘图完成，输出文件夹：", normalizePath(outdir), "\n")
}

############################################################
## Panel D: cfDNA 96-trinucleotide mutation profile
############################################################
run_figure2d <- function(input_dir, output_dir) {
  multianno_file <- file.path(input_dir, "Gastric_all_sample_multianno.txt")
  cfdna_matrix_file <- file.path(input_dir, "oncoplot_matrix_CFDNA_Gastric.txt")
  tcga_matrix_file <- file.path(input_dir, "oncoplot_matrix_TCGA_STAD.txt")
  clinical_file <- file.path(input_dir, "patient_info.txt")
  outdir <- output_dir
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
read_table_auto <- function(file_path, row_names = NULL) {
  tryCatch(
    read.delim(
      file_path,
      header = TRUE,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      row.names = row_names,
      fileEncoding = "UTF-8-BOM"
    ),
    error = function(e) {
      read.delim(
        file_path,
        header = TRUE,
        check.names = FALSE,
        stringsAsFactors = FALSE,
        row.names = row_names
      )
    }
  )
}

variant_type_from_ref_alt <- function(ref, alt) {
  ref <- as.character(ref)
  alt <- as.character(alt)
  dplyr::case_when(
    nchar(ref) == 1 & nchar(alt) == 1 ~ "SNP",
    nchar(ref) > nchar(alt) ~ "DEL",
    nchar(ref) < nchar(alt) ~ "INS",
    TRUE ~ "ONP"
  )
}

cfdna_matrix <- read_table_auto(cfdna_matrix_file, row_names = 1)
cfdna_matrix[] <- lapply(cfdna_matrix, function(x) as.numeric(as.character(x)))
cfdna_matrix[is.na(cfdna_matrix)] <- 0
cfdna_matrix <- as.matrix(cfdna_matrix)

tcga_matrix <- read_table_auto(tcga_matrix_file, row_names = 1)

sample_ids <- trimws(as.character(rownames(cfdna_matrix)))

cat("Gastric cfDNA samples: ", length(sample_ids), "\n", sep = "")
cat("TCGA-STAD matrix samples loaded for reference: ", nrow(tcga_matrix), "\n", sep = "")

multianno_raw <- read_table_auto(multianno_file)

cfDNA_snv <- multianno_raw %>%
  dplyr::transmute(
    Chromosome = as.character(.data[["Chr"]]),
    Start_Position = suppressWarnings(as.integer(.data[["Start"]])),
    End_Position = suppressWarnings(as.integer(.data[["End"]])),
    Reference_Allele = toupper(as.character(.data[["Ref"]])),
    Tumor_Seq_Allele2 = toupper(as.character(.data[["Alt"]])),
    Tumor_Sample_Barcode = as.character(.data[["Tumor_Sample_Barcode"]]),
    Variant_Type = variant_type_from_ref_alt(.data[["Ref"]], .data[["Alt"]])
  ) %>%
  dplyr::filter(
    Tumor_Sample_Barcode %in% sample_ids,
    Variant_Type == "SNP",
    !is.na(Chromosome),
    !is.na(Start_Position),
    !is.na(End_Position),
    Reference_Allele %in% c("A", "T", "C", "G"),
    Tumor_Seq_Allele2 %in% c("A", "T", "C", "G"),
    Reference_Allele != Tumor_Seq_Allele2
  ) %>%
  dplyr::mutate(
    Chromosome = ifelse(grepl("^chr", Chromosome), Chromosome, paste0("chr", Chromosome)),
    sample_id_use = Tumor_Sample_Barcode
  ) %>%
  dplyr::distinct(
    Chromosome,
    Start_Position,
    End_Position,
    Reference_Allele,
    Tumor_Seq_Allele2,
    sample_id_use,
    .keep_all = TRUE
  )

main_chr <- paste0("chr", c(1:22, "X", "Y"))
cfDNA_snv <- cfDNA_snv %>%
  dplyr::filter(Chromosome %in% main_chr)

if (nrow(cfDNA_snv) == 0) {
  stop("没有可用于 96 突变谱的胃癌 cfDNA SNP。")
}

write.table(
  cfDNA_snv,
  file.path(outdir, "gastric_cfDNA_snv_for_96_profile.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("SNVs used for 96 profile: ", nrow(cfDNA_snv), "\n", sep = "")
cat("Samples used for 96 profile: ", length(unique(cfDNA_snv$sample_id_use)), "\n", sep = "")

vcfs_list <- split(cfDNA_snv, cfDNA_snv$sample_id_use)

vcfs_gr <- lapply(vcfs_list, function(df) {
  gr <- GenomicRanges::GRanges(
    seqnames = df$Chromosome,
    ranges = IRanges::IRanges(
      start = df$Start_Position,
      end = df$Start_Position
    )
  )
  GenomicRanges::mcols(gr)$REF <- df$Reference_Allele
  GenomicRanges::mcols(gr)$ALT <- df$Tumor_Seq_Allele2
  GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"
  GenomeInfoDb::genome(gr) <- "hg38"
  gr
})

ref_genome <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
mut_mat <- MutationalPatterns::mut_matrix(
  vcf_list = vcfs_gr,
  ref_genome = ref_genome
)

gastric_profile <- rowSums(mut_mat)
mut_mat_with_total <- cbind(mut_mat, Gastric_cfDNA = gastric_profile)

write.csv(
  mut_mat_with_total,
  file.path(outdir, "gastric_cfDNA_96_mutation_matrix.csv"),
  row.names = TRUE,
  fileEncoding = "UTF-8"
)

profile_plot <- MutationalPatterns::plot_96_profile(
  mut_mat_with_total[, "Gastric_cfDNA", drop = FALSE],
  condensed = TRUE
) +
  ggplot2::ggtitle("Gastric cfDNA 96 trinucleotide mutation profile")

pdf(
  file.path(outdir, "Figure2D_gastric_96_profile.pdf"),
  width = 10,
  height = 5,
  useDingbats = FALSE
)

print(profile_plot)

dev.off()

png(
  file.path(outdir, "Figure2D_gastric_96_profile.png"),
  width = 3000,
  height = 1500,
  res = 300
)

print(profile_plot)

dev.off()

cat("Figure2D gastric 96 profile saved: ", file.path(outdir, "Figure2D_gastric_96_profile.pdf"), "\n", sep = "")
}

run_figure2a(input_dir, figure2_output_dir)
run_figure2b(input_dir, figure2_output_dir)
run_figure2c(input_dir, figure2_output_dir)
run_figure2d(input_dir, figure2_output_dir)

message("Figure 2 panels A-D completed. Output directory: ", figure2_output_dir)
