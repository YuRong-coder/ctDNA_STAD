############################################################
## Figure 4 unified TCGA-STAD analysis pipeline (panels A-D)
##
## Panel A: typed 20-gene waterfall from the first source script.
## Panels B-D: Stage, pT and pN dot plots from the Xena source script.
##
## Put all input files in ./data:
##   - TCGA-STAD.somaticmutation_wxs.tsv
##   - TCGA-STAD.clinical.gz
##   - TCGA-STAD.survival.gz
##   - oncoplot_matrix_TCGA_STAD.txt
##
## Run with: Rscript Figure4_analysis.R
## Outputs are written under ./output/Figure4.
############################################################

required_packages <- c("dplyr", "maftools", "tidyr", "ggplot2", "scales")
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
  library(dplyr)
  library(maftools)
  library(tidyr)
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
input_dir <- file.path(script_dir, "data")
figure4_output_dir <- file.path(script_dir, "output", "Figure4")
dir.create(figure4_output_dir, recursive = TRUE, showWarnings = FALSE)

required_input_files <- c(
  "TCGA-STAD.somaticmutation_wxs.tsv",
  "TCGA-STAD.clinical.gz",
  "TCGA-STAD.survival.gz",
  "oncoplot_matrix_TCGA_STAD.txt"
)
missing_input_files <- required_input_files[
  !file.exists(file.path(input_dir, required_input_files))
]
if (length(missing_input_files) > 0L) {
  stop("data/ 中缺少输入文件：", paste(missing_input_files, collapse = ", "))
}

############################################################
## Panel A: Typed 20-gene waterfall (first source script)
############################################################
run_figure4a <- function(input_dir, output_root) {
  out_dir <- file.path(output_root, "A_typed_20gene_waterfall")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  mutation_file <- file.path(input_dir, "TCGA-STAD.somaticmutation_wxs.tsv")
  clinical_file <- file.path(input_dir, "TCGA-STAD.clinical.gz")
  survival_file <- file.path(input_dir, "TCGA-STAD.survival.gz")
genes <- c(
  "ARID1B", "ATM", "BRCA2", "C1S",
  "CCL16", "CD72", "CDH11", "CDH6",
  "CLDN11", "COL12A1", "COL15A1", "COL23A1",
  "EP300", "EPHB2", "EPHB6", "ERBB4",
  "HSPG2", "IL10RB", "IL16", "ITGA4"
)

mut_col <- c(

  Missense_Mutation = "#43A2CA",
  Nonsense_Mutation = "#F8D042",
  Frame_Shift_Del = "#377EB8",
  Frame_Shift_Ins = "#C2E2C2",
  In_Frame_Del = "#F39B7F",
  In_Frame_Ins = "#8491B4",
  Splice_Site = "#B0D295",
  Translation_Start_Site = "#F0A0A0",
  Nonstop_Mutation = "#7E6148"

  # Missense_Mutation = "#377EB8", Nonsense_Mutation = "#E64B35",
  # Frame_Shift_Del = "#4DBBD5", Frame_Shift_Ins = "#00A087",
  # In_Frame_Del = "#F39B7F", In_Frame_Ins = "#8491B4",
  # Splice_Site = "#91D1C2", Translation_Start_Site = "#DC0000",
  # Nonstop_Mutation = "#7E6148"
)

mutation_raw <- read.delim(mutation_file, check.names = FALSE, stringsAsFactors = FALSE)
clinical_raw <- read.delim(clinical_file, check.names = FALSE, stringsAsFactors = FALSE)
survival_raw <- read.delim(survival_file, check.names = FALSE, stringsAsFactors = FALSE)

classify_effect <- function(effect, ref, alt) {
  case_when(
    grepl("stop_gained", effect) ~ "Nonsense_Mutation",
    grepl("stop_lost", effect) ~ "Nonstop_Mutation",
    grepl("start_lost", effect) ~ "Translation_Start_Site",
    grepl("splice_acceptor|splice_donor", effect) ~ "Splice_Site",
    grepl("frameshift", effect) & nchar(alt) > nchar(ref) ~ "Frame_Shift_Ins",
    grepl("frameshift", effect) ~ "Frame_Shift_Del",
    grepl("inframe_insertion", effect) ~ "In_Frame_Ins",
    grepl("inframe_deletion", effect) ~ "In_Frame_Del",
    grepl("missense", effect) ~ "Missense_Mutation",
    TRUE ~ NA_character_
  )
}

maf_selected <- mutation_raw %>%
  filter(gene %in% genes) %>%
  mutate(
    Variant_Classification = classify_effect(effect, ref, alt),
    Variant_Type = case_when(nchar(ref) == 1 & nchar(alt) == 1 ~ "SNP",
                             nchar(ref) > nchar(alt) ~ "DEL", TRUE ~ "INS")
  ) %>%
  filter(!is.na(Variant_Classification)) %>%
  transmute(
    Hugo_Symbol = gene, Chromosome = sub("^chr", "", chrom),
    Start_Position = as.integer(start), End_Position = as.integer(end),
    Reference_Allele = ref, Tumor_Seq_Allele1 = ref, Tumor_Seq_Allele2 = alt,
    Variant_Classification, Variant_Type,
    Tumor_Sample_Barcode = sample,
    Full_Tumor_Barcode = Tumor_Sample_Barcode,
    Amino_Acid_Change, Original_effect = effect, callers, dna_vaf
  )

# Keep only samples with at least one qualifying mutation in the selected genes.
mutated_samples <- unique(maf_selected$Tumor_Sample_Barcode)
if (length(mutated_samples) == 0) stop("No qualifying mutations found in the 20 target genes.")

clean_tnm <- function(x, letter) {
  x <- toupper(trimws(as.character(x)))
  ans <- sub(paste0(".*(", letter, "[0-4X][A-D]?).*"), "\\1", x)
  ifelse(is.na(x) | x == "" | !grepl(paste0(letter, "[0-4X]"), x), "Unknown", ans)
}
clean_stage <- function(x) {
  case_when(
    grepl("Stage IV", x, ignore.case = TRUE) ~ "IV",
    grepl("Stage III", x, ignore.case = TRUE) ~ "III",
    grepl("Stage II", x, ignore.case = TRUE) ~ "II",
    grepl("Stage I", x, ignore.case = TRUE) ~ "I",
    TRUE ~ "Unknown"
  )
}

clinical_clean <- clinical_raw %>%
  transmute(
    Tumor_Sample_Barcode = sample,
    Patient = substr(sample, 1, 12),
    Gender = case_when(tolower(gender.demographic) == "female" ~ "Female",
                       tolower(gender.demographic) == "male" ~ "Male", TRUE ~ "Unknown"),
    Age = case_when(is.na(as.numeric(age_at_index.demographic)) ~ "Unknown",
                    as.numeric(age_at_index.demographic) < 60 ~ "<60", TRUE ~ ">=60"),
    pT = clean_tnm(ajcc_pathologic_t.diagnoses, "T"),
    pN = clean_tnm(ajcc_pathologic_n.diagnoses, "N"),
    pM = clean_tnm(ajcc_pathologic_m.diagnoses, "M"),
    Stage = clean_stage(ajcc_pathologic_stage.diagnoses),
    Sample_Type = ifelse(is.na(sample_type.samples) | sample_type.samples == "",
                         "Unknown", sample_type.samples)
  ) %>% distinct(Tumor_Sample_Barcode, .keep_all = TRUE)

survival_clean <- survival_raw %>%
  transmute(Tumor_Sample_Barcode = sample,
            OS_time_days = suppressWarnings(as.numeric(OS.time)),
            OS_event = suppressWarnings(as.integer(OS)),
            OS_status = case_when(OS_event == 1 ~ "Dead", OS_event == 0 ~ "Alive", TRUE ~ "Unknown")) %>%
  distinct(Tumor_Sample_Barcode, .keep_all = TRUE)

clinical_selected <- data.frame(Tumor_Sample_Barcode = mutated_samples) %>%
  left_join(clinical_clean, by = "Tumor_Sample_Barcode") %>%
  left_join(survival_clean, by = "Tumor_Sample_Barcode") %>%
  mutate(across(c(Gender, Age, pT, pN, pM, Stage, Sample_Type, OS_status),
                ~ifelse(is.na(.x) | .x == "", "Unknown", .x)))

write.table(maf_selected, file.path(out_dir, "01_TCGA_STAD_20genes_typed_mutations.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(clinical_selected, file.path(out_dir, "02_mutated_samples_clinical_survival_TNM.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(as.data.frame(table(maf_selected$Hugo_Symbol, maf_selected$Variant_Classification)),
            file.path(out_dir, "03_gene_by_mutation_type_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(sample = mutated_samples), file.path(out_dir, "04_retained_mutated_samples.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

maf_obj <- read.maf(maf = maf_selected, clinicalData = clinical_selected, verbose = FALSE)

ann_colors <- list(

  Age = c("<60" = "#B2D990", ">=60" = "#F0A0A0", "Unknown" = "#D9D9D9"),
  Gender = c("Male" = "#82D4D7", "Female" = "#DC7579", "Unknown" = "#D9D9D9"),
  pT = c("T0" = "#F5FFFA", "T1" = "#D9F0D3", "T1A" = "#C7E9C0", "T1B" = "#CDE5CD", "T2" = "#C2E2C2", "T3" = "#B9DFB9", "T4" = "#AFD1BF", "T4A" = "#90BFCF", "T4B" = "#88AFD8", "TX" = "#AAAAAA", "Unknown" = "#D9D9D9"),
  pN = c("N0" = "#D8DEF7", "N1" = "#BCBDDC", "N2" = "#807DBA", "N3" = "#9898DC", "N3A" = "#897CD3", "N3B" = "#5E379D", "NX" = "#AAAAAA", "Unknown" = "#D9D9D9"),
  pM = c("M0" = "#82B6E8", "M1" = "#EEAAAA", "MX" = "#AAAAAA", "Unknown" = "#D9D9D9"),
  Stage = c("I" = "#FFF3B0", "II" = "#FBE3D2", "III" = "#F8D0B0", "IV" = "#F2B382", "Unknown" = "#D9D9D9"),
  #
  #
  # Gender = c(Female="#E64B35", Male="#4DBBD5", Unknown="#D9D9D9"),
  # Age = c("<60"="#F4A582", ">=60"="#92C5DE", Unknown="#D9D9D9"),
  # pT = c(T0="#F2E2F2", T1="#D8B7D8", T1A="#D8B7D8", T1B="#CFA3D1", T2="#B987C0", T3="#8B4B9C", T4="#5A236E", T4A="#5A236E", T4B="#40184F", TX="#AAAAAA", Unknown="#D9D9D9"),
  # pN = c(N0="#D9D9D9", N1="#B6A6D8", N2="#8C6BB1", N3="#5E3C99", N3A="#5E3C99", N3B="#42246F", NX="#AAAAAA", Unknown="#D9D9D9"),
  # pM = c(M0="#BDBDBD", M1="#E64B35", MX="#AAAAAA", Unknown="#D9D9D9"),
  # Stage = c(I="#DDECC9", II="#A8DDB5", III="#43A2CA", IV="#0868AC", Unknown="#D9D9D9"),

  OS_status = c(Alive="#8CB26C", Dead="#EA9E58", Unknown="#ECF0F1")
)

features <- c("Gender", "Age", "pT", "pN", "pM", "Stage", "OS_status")

# Guarantee that every value present in the selected clinical data has a color.
for (feature in features) {
  observed <- sort(unique(as.character(clinical_selected[[feature]])))
  observed <- observed[!is.na(observed) & observed != ""]
  missing_levels <- setdiff(observed, names(ann_colors[[feature]]))
  if (length(missing_levels) > 0) {
    extra_colors <- grDevices::hcl.colors(length(missing_levels), "Set 3")
    names(extra_colors) <- missing_levels
    ann_colors[[feature]] <- c(ann_colors[[feature]], extra_colors)
  }
}

plot_fun <- function() oncoplot(
  maf = maf_obj, genes = genes, colors = mut_col,
  removeNonMutated = TRUE, clinicalFeatures = features,
  annotationColor = ann_colors, sortByAnnotation = FALSE,
  showTumorSampleBarcodes = FALSE, drawRowBar = TRUE, drawColBar = TRUE,
  fontSize = 0.7,
  titleText = paste0("TCGA-STAD: 20-gene mutation landscape (n = ", length(mutated_samples), ")")
)

grDevices::cairo_pdf(
  file.path(out_dir, "TCGA_STAD_20genes_typed_TNM_waterfall_compatible.pdf"),
  width = 15, height = 10, onefile = TRUE
)
plot_fun(); dev.off()
png(file.path(out_dir, "TCGA_STAD_20genes_typed_TNM_waterfall.png"), width = 4500, height = 3000, res = 300)
plot_fun(); dev.off()
tiff(file.path(out_dir, "TCGA_STAD_20genes_typed_TNM_waterfall.tiff"), width = 4500, height = 3000,
     res = 300, compression = "lzw")
plot_fun(); dev.off()

download_page <- "https://xenabrowser.net/datapages/?dataset=TCGA-STAD.somaticmutation_wxs.tsv&host=https%3A%2F%2Ftcga.xenahubs.net"
writeLines(c(
  paste("Mutation input:", mutation_file), paste("Clinical input:", clinical_file),
  paste("Survival input:", survival_file), paste("Retained mutated samples:", length(mutated_samples)),
  paste("Selected typed mutations:", nrow(maf_selected)),
  paste("Original data download page:", download_page),
  "Only samples mutated in at least one of the 20 target genes are retained.",
  "Synonymous, intronic and other effects outside mut_col are excluded."
), file.path(out_dir, "README_inputs_outputs_download_URL.txt"))
capture.output(sessionInfo(), file = file.path(out_dir, "05_sessionInfo.txt"))
message("Completed: ", out_dir)
}

############################################################
## Panels B-D: Stage, pT and pN dot plots (second source script)
############################################################
run_figure4bcd <- function(input_dir, output_root) {
  out_dir <- file.path(output_root, "BCD_Xena_dotplots")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- list(
    clinical = file.path(input_dir, "TCGA-STAD.clinical.gz"),
    survival = file.path(input_dir, "TCGA-STAD.survival.gz"),
    mutation_matrix = file.path(input_dir, "oncoplot_matrix_TCGA_STAD.txt")
  )
# 定义本次分析候选驱动基因集，toupper统一大写去重
candidate_genes <- unique(toupper(c(
  "ARID1B", "ATM", "BRCA2", "C1S",
  "CCL16", "CD72", "CDH11", "CDH6",
  "CLDN11", "COL12A1", "COL15A1", "COL23A1",
  "EP300", "EPHB2", "EPHB6", "ERBB4",
  "HSPG2", "IL10RB", "IL16", "ITGA4"
)))

# ===================== 2. 自定义通用清洗函数 =====================
## 读取制表符分隔tsv/gz文件，修复空列名、重命名重复列
read_tsv <- function(path) {
  x <- read.delim(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  # 识别空列名，填充临时sample_id编号
  empty_names <- is.na(names(x)) | names(x) == ""
  names(x)[empty_names] <- paste0("sample_id", seq_len(sum(empty_names)))
  # 处理重复列名，自动加后缀区分
  names(x) <- make.unique(names(x), sep = "_")
  x
}

## 将空字符串、未知标记统一替换为标准NA缺失值
na_if_blank <- function(x) {
  x <- as.character(x)
  x[x %in% c("", ".", "--", "NA", "NaN", "not reported", "Not Reported", "unknown", "Unknown")] <- NA_character_
  x
}

## AJCC病理分期标准化：I/II/III/IV四分类
clean_stage <- function(x) {
  x <- na_if_blank(x)
  case_when(
    grepl("Stage I[A-C]*$", x, ignore.case = TRUE) ~ "Stage I",
    grepl("Stage II[A-C]*$", x, ignore.case = TRUE) ~ "Stage II",
    grepl("Stage III[A-C]*$", x, ignore.case = TRUE) ~ "Stage III",
    grepl("Stage IV", x, ignore.case = TRUE) ~ "Stage IV",
    TRUE ~ NA_character_
  )
}

## T原发肿瘤分期标准化：T0/is/T1/T2/T3/T4
clean_t <- function(x) {
  x <- na_if_blank(x)
  case_when(
    grepl("Tis|T0", x, ignore.case = TRUE) ~ "T0/is",
    grepl("T1", x, ignore.case = TRUE) ~ "T1",
    grepl("T2", x, ignore.case = TRUE) ~ "T2",
    grepl("T3", x, ignore.case = TRUE) ~ "T3",
    grepl("T4", x, ignore.case = TRUE) ~ "T4",
    TRUE ~ NA_character_
  )
}

## N淋巴结转移分期标准化：N0/N1/N2/N3
clean_n <- function(x) {
  x <- na_if_blank(x)
  case_when(
    grepl("N0", x, ignore.case = TRUE) ~ "N0",
    grepl("N1", x, ignore.case = TRUE) ~ "N1",
    grepl("N2", x, ignore.case = TRUE) ~ "N2",
    grepl("N3", x, ignore.case = TRUE) ~ "N3",
    TRUE ~ NA_character_
  )
}

## M远处转移分期标准化：M0无转移/M1远处转移
clean_m <- function(x) {
  x <- na_if_blank(x)
  case_when(
    grepl("M0", x, ignore.case = TRUE) ~ "M0",
    grepl("M1", x, ignore.case = TRUE) ~ "M1",
    TRUE ~ NA_character_
  )
}

## 输出清洗后表格到指定输出目录，tsv无引号、无行号
save_table <- function(x, filename) {
  write.table(x, file.path(out_dir, filename), sep = "\t", quote = FALSE, row.names = FALSE)
}

# ===================== 3. 临床、生存数据读取与清洗 =====================
# 读取原始临床表、生存表
clinical_raw <- read_tsv(paths$clinical)
survival_raw <- read_tsv(paths$survival)

# 临床信息精简清洗：提取关键临床变量、标准化分期、年龄分组
clinical_clean <- clinical_raw %>%
  transmute(
    sample_id = substr(sample, 1, 16), # 截取TCGA标准16位样本ID用于匹配突变矩阵
    patient_id = substr(submitter_id, 1, 12), # 12位患者ID
    Gender = na_if_blank(tolower(gender.demographic)), # 性别统一小写+缺失清洗
    Age = suppressWarnings(as.numeric(age_at_index.demographic)), # 年龄转为数值，抑制转换警告
    Age_Group = case_when(is.na(Age) ~ NA_character_, Age < 60 ~ "<60", TRUE ~ ">=60"), # 年龄二分类
    Stage = clean_stage(ajcc_pathologic_stage.diagnoses), # 标准化总病理分期
    pT = clean_t(ajcc_pathologic_t.diagnoses), # 标准化T分期
    pN = clean_n(ajcc_pathologic_n.diagnoses), # 标准化N分期
    pM = clean_m(ajcc_pathologic_m.diagnoses), # 标准化M分期
    Grade = na_if_blank(tumor_grade.diagnoses), # 肿瘤分化等级
    Sample_Type = na_if_blank(sample_type.samples), # 样本类型：原发瘤/转移瘤/正常组织
    Vital_Status = na_if_blank(vital_status.demographic) # 生存状态：存活/死亡
  ) %>%
  filter(!is.na(sample_id), sample_id != "") %>% # 过滤无样本ID无效行
  distinct(sample_id, .keep_all = TRUE) # 去重，保留第一条完整记录

# 生存数据清洗：提取OS总生存时间、死亡事件状态
survival_clean <- survival_raw %>%
  transmute(
    sample_id = substr(sample, 1, 16), # 统一样本ID格式用于合并
    patient_id = `_PATIENT`,
    OS_time = suppressWarnings(as.numeric(OS.time)), # 生存随访天数
    OS_event = suppressWarnings(as.integer(OS)), # 生存事件：1=死亡，0=截尾存活
    OS_status = ifelse(OS_event == 1, "Dead", "Alive/Censored") # 文字化生存标签
  ) %>%
  filter(!is.na(sample_id), sample_id != "") %>%
  distinct(sample_id, .keep_all = TRUE)

# 临床表 + 生存表左连接整合，一份完整临床-生存注释表
clinical_anno <- clinical_clean %>%
  left_join(survival_clean %>% dplyr::select(sample_id, OS_time, OS_event, OS_status), by = "sample_id")

# ===================== 4. 突变矩阵处理 + TMB肿瘤突变负荷计算 =====================
mutation_matrix <- read_tsv(paths$mutation_matrix)
names(mutation_matrix)[1] <- "sample_id" # 首列强制命名为sample_id，统一匹配键
mutation_matrix <- mutation_matrix %>%
  filter(!is.na(sample_id), sample_id != "") %>%
  # 所有基因列转为整数0(野生型)/1(突变型)，缺失填充为0
  mutate(across(-sample_id, ~ {
    y <- suppressWarnings(as.integer(.x))
    y[is.na(y)] <- 0L
    y
  })) %>%
  distinct(sample_id, .keep_all = TRUE)

# 提取所有基因列名（排除样本ID列）
gene_cols <- setdiff(names(mutation_matrix), "sample_id")
# 逐样本统计突变基因总数，估算TMB（除以38Mb外显子覆盖长度）
sample_burden <- mutation_matrix %>%
  transmute(
    sample_id,
    mutation_count = rowSums(across(all_of(gene_cols)) > 0, na.rm = TRUE), # 单个样本突变基因数量
    estimated_TMB = mutation_count / 38 # 估算肿瘤突变负荷TMB，单位：突变/Mb
  )

# 候选基因矩阵转为长格式，仅保留突变样本（Mutated=1），用于气泡图统计
matrix_long_candidate <- mutation_matrix %>%
  dplyr::select(sample_id, any_of(intersect(candidate_genes, gene_cols))) %>%
  pivot_longer(-sample_id, names_to = "Gene", values_to = "Mutated") %>%
  filter(Mutated > 0)

# ===================== 5. 输出全部中间清洗表格（溯源备用） =====================
save_table(clinical_raw, "intermediate_01_TCGA_STAD_clinical_raw.txt")
save_table(survival_raw, "intermediate_02_TCGA_STAD_survival_raw.txt")
save_table(clinical_anno, "intermediate_03_TCGA_STAD_clinical_survival_clean.txt")
save_table(sample_burden, "intermediate_04_TCGA_STAD_sample_TMB_from_matrix.txt")
save_table(matrix_long_candidate, "intermediate_05_TCGA_STAD_candidate_mutation_long.txt")

# ===================== 7. 函数：气泡图数据统计（分临床亚组计算突变频率、中位TMB） =====================
# group_var：分组变量名(Stage/pT/pN/Sample_Type)；top_n：每组取突变频率前N基因
make_group_dotplot_data <- function(group_var, top_n = 5) {
  # 临床分组样本过滤，剔除缺失分组样本
  group_df <- clinical_anno %>%
    dplyr::select(sample_id, Group = all_of(group_var)) %>%
    filter(!is.na(Group), Group != "") %>%
    distinct(sample_id, .keep_all = TRUE)

  # 统计每个分组总样本量
  group_n <- group_df %>% count(Group, name = "total_samples")
  # 筛选交集：同时存在于突变矩阵和候选集的基因
  genes_use <- intersect(candidate_genes, gene_cols)

  # 突变矩阵转长格式，仅保留突变样本，关联分组、TMB数据
  mutation_matrix %>%
    dplyr::select(sample_id, all_of(genes_use)) %>%
    pivot_longer(-sample_id, names_to = "Gene", values_to = "Mutated") %>%
    filter(Mutated > 0) %>%
    inner_join(group_df, by = "sample_id") %>%
    left_join(sample_burden, by = "sample_id") %>%
    # 按分组+基因聚合统计指标
    group_by(Group, Gene) %>%
    summarise(
      mutated_samples = n_distinct(sample_id), # 该分组携带该基因突变的样本数
      median_estimated_TMB = median(estimated_TMB, na.rm = TRUE), # 突变样本中位TMB（气泡颜色映射）
      mean_estimated_TMB = mean(estimated_TMB, na.rm = TRUE),
      median_mutation_burden = median(mutation_count, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(group_n, by = "Group") %>%
    mutate(mutation_frequency = mutated_samples / total_samples * 100) %>% # 突变频率%（气泡大小映射）
    # 分组内按突变频率降序，取top5高频基因
    group_by(Group) %>%
    arrange(desc(mutation_frequency), desc(mutated_samples), desc(median_estimated_TMB), .by_group = TRUE) %>%
    slice_head(n = top_n) %>%
    ungroup()
}

# ===================== 8. 函数：绘制分组TMB突变气泡图（Figure4B/C/D/E共用绘图模板） =====================
# dot_df：make_group_dotplot_data输出统计表格；group_order：X轴分组自定义顺序；title：子图标题
plot_group_dotplot <- function(dot_df, group_order, title) {
  # 过滤掉数据中不存在的分组，保证绘图顺序正确
  group_order <- group_order[group_order %in% unique(dot_df$Group)]
  group_order <- c(group_order, setdiff(unique(dot_df$Group), group_order))
  # 基因绘图排序
  gene_order <- dot_df %>%
    mutate(Group = factor(Group, levels = group_order)) %>%
    arrange(Group, desc(mutation_frequency), desc(mutated_samples), desc(median_estimated_TMB)) %>%
    pull(Gene) %>%
    unique()

  # 气泡图主绘图代码
  dot_df %>%
    mutate(
      Group = factor(Group, levels = group_order),
      Gene = factor(Gene, levels = rev(gene_order))
    ) %>%
    ggplot(aes(Group, Gene)) +
    geom_point(aes(size = mutation_frequency, color = median_estimated_TMB), alpha = 0.9) +
    scale_size_continuous(range = c(3, 10), name = "Mutation frequency (%)") + # 气泡大小范围
    scale_color_gradient(low = "lightblue", high = "#D73027", name = "Median estimated TMB") + # 蓝-红TMB渐变色
    labs(
      x = NULL,
      y = NULL,
      title = title,
      subtitle = "Dot color represents median estimated TMB of samples carrying each gene mutation"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      plot.subtitle = element_text(hjust = 0.5, size = 8, color = "grey30"),
      axis.text.x = element_text(angle = 35, hjust = 1, color = "black", face = "bold"),
      axis.text.y = element_text(color = "black", face = "italic"),
      legend.position = "right"
    )
}

# ===================== 9. 绘图输出通用函数：同时导出PDF矢量图+300dpi PNG位图 =====================
save_plot <- function(plot, filename, width, height) {
  ggsave(file.path(out_dir, paste0(filename, ".pdf")), plot = plot, width = width, height = height)
  ggsave(file.path(out_dir, paste0(filename, ".png")), plot = plot, width = width, height = height, dpi = 300)
}

# ===================== 11. 统计 Figure4B/C/D 气泡图数据 =====================
dot_stage <- make_group_dotplot_data("Stage")      # 按总病理Stage分期分组
dot_pt <- make_group_dotplot_data("pT")            # 按T浸润深度分组
dot_pn <- make_group_dotplot_data("pN")            # 按淋巴结N转移分组

# 保存3套气泡图统计中间表
save_table(dot_stage, "intermediate_06_Figure4B_TCGA_STAD_stage_dotplot_data.txt")
save_table(dot_pt, "intermediate_07_Figure4C_TCGA_STAD_pT_dotplot_data.txt")
save_table(dot_pn, "intermediate_08_Figure4D_TCGA_STAD_pN_dotplot_data.txt")

# ===================== 12. 分别绘制 Figure4B/C/D =====================
p_stage <- plot_group_dotplot(dot_stage, c("Stage I", "Stage II", "Stage III", "Stage IV"),
                              "Top 5 mutated candidate genes across TCGA-STAD stages")
p_pt <- plot_group_dotplot(dot_pt, c("T0/is", "T1", "T2", "T3", "T4"),
                           "Top 5 mutated candidate genes across TCGA-STAD pT stages")
p_pn <- plot_group_dotplot(dot_pn, c("N0", "N1", "N2", "N3"),
                           "Top 5 mutated candidate genes across TCGA-STAD pN stages")

# 单张子图单独导出
save_plot(p_stage, "Figure4B_TCGA_STAD_stage_TMB_dotplot", 6, 7)
save_plot(p_pt, "Figure4C_TCGA_STAD_pT_TMB_dotplot", 6, 7)
save_plot(p_pn, "Figure4D_TCGA_STAD_pN_TMB_dotplot", 6, 7)
  message("TCGA-STAD Figure4 panels B-D finished. Output: ", out_dir)
}

run_figure4a(input_dir, figure4_output_dir)
run_figure4bcd(input_dir, figure4_output_dir)

message("Figure 4 panels A-D completed. Output directory: ", figure4_output_dir)
