############################################################
## Figure 3 unified analysis pipeline (panels A-D)
##
## Put all input files in ./data:
##   - gastric_maf_with_clinical.rds
##   - oncoplot_matrix_TCGA_STAD.txt
##   - gastric_TME_gene_literature_support_intersection_optimized_60.csv
##
## Run with: Rscript Figure3_analysis.R
## Outputs are written under ./output/Figure3.
############################################################

required_packages <- c(
  "maftools", "data.table", "dplyr", "ggplot2", "ggVennDiagram",
  "tibble", "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi"
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
  library(maftools)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggVennDiagram)
  library(tibble)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
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
figure3_output_dir <- file.path(script_dir, "output", "Figure3")
dir.create(figure3_output_dir, recursive = TRUE, showWarnings = FALSE)

required_input_files <- c(
  "gastric_maf_with_clinical.rds",
  "oncoplot_matrix_TCGA_STAD.txt",
  "gastric_TME_gene_literature_support_intersection_optimized_60.csv"
)
missing_input_files <- required_input_files[
  !file.exists(file.path(input_dir, required_input_files))
]
if (length(missing_input_files) > 0L) {
  stop("data/ 中缺少输入文件：", paste(missing_input_files, collapse = ", "))
}

############################################################
## Panel A: Gene-selection Venn diagram
############################################################
run_figure3a <- function(input_dir, output_root) {
  outdir <- file.path(output_root, "A_gene_selection_venn")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  own_maf_rds_file <- file.path(input_dir, "gastric_maf_with_clinical.rds")
  tcga_matrix_file <- file.path(input_dir, "oncoplot_matrix_TCGA_STAD.txt")
  literature_file <- file.path(input_dir, "gastric_TME_gene_literature_support_intersection_optimized_60.csv")
# ===================== 3. 自定义通用工具函数 =====================
## 函数：标准化基因名（去空格、全部大写，避免匹配失败）
norm_gene <- function(x) {
  toupper(trimws(as.character(x)))
}

## 函数：读取TCGA二元突变矩阵（第一列为样本名，剩余列为基因）
read_binary_matrix <- function(path) {
  mat <- fread(path, data.table = FALSE, check.names = FALSE)
  row_ids <- mat[[1]]          # 提取样本行名
  mat <- mat[, -1, drop = FALSE] # 删除样本名列，只保留0/1突变矩阵
  rownames(mat) <- row_ids
  # 全部转为数值，缺失值填充0（无突变）
  mat[] <- lapply(mat, function(x) as.numeric(as.character(x)))
  mat[is.na(mat)] <- 0
  as.matrix(mat)
}

# ===================== 4. 提取三组独立基因集 =====================
## 4.1 提取本队列所有非同义突变基因
maf_sub <- readRDS(own_maf_rds_file) # 读取保存好的MAF对象

own_genes <- as.data.frame(maf_sub@data) %>%
  transmute(Gene = norm_gene(Hugo_Symbol)) %>% # 标准化基因名
  filter(!is.na(Gene), Gene != "") %>%         # 过滤空基因名
  distinct(Gene) %>%                           # 去重，每个基因只保留一次
  arrange(Gene) %>%                            # 字母排序
  pull(Gene)                                   # 转为字符向量存储

## 4.2 提取TCGA-STAD中发生突变的全部基因
tcga_matrix <- read_binary_matrix(tcga_matrix_file)
colnames(tcga_matrix) <- norm_gene(colnames(tcga_matrix)) # 标准化列名（基因）
# 统计每列总和>0：至少1例样本突变的基因
tcga_genes <- colnames(tcga_matrix)[colSums(tcga_matrix > 0, na.rm = TRUE) > 0] %>%
  unique() %>%
  sort()

## 4.3 提取文献TME相关基因
literature_genes <- fread(literature_file, data.table = FALSE, check.names = FALSE) %>%
  transmute(Gene = norm_gene(Gene)) %>%
  filter(!is.na(Gene), Gene != "") %>%
  distinct(Gene) %>%
  arrange(Gene) %>%
  pull(Gene)

# 组装三组基因列表，用于韦恩图输入
sets <- list(
  `Own non-synonymous` = own_genes,
  `TCGA-STAD` = tcga_genes,
  `Literature TME` = literature_genes
)

# ===================== 5. 计算各类交集/独有基因 =====================
# 三者共同交集（最终候选基因）
venn_intersection <- Reduce(intersect, sets)

# 两组交集
own_tcga <- intersect(own_genes, tcga_genes)
own_literature <- intersect(own_genes, literature_genes)
tcga_literature <- intersect(tcga_genes, literature_genes)

# 拆分韦恩图7个独立分区：仅A、仅B、仅C、AB独有、AC独有、BC独有、三者共有
region_lists <- list(
  own_only = setdiff(own_genes, union(tcga_genes, literature_genes)),
  tcga_only = setdiff(tcga_genes, union(own_genes, literature_genes)),
  literature_only = setdiff(literature_genes, union(own_genes, tcga_genes)),
  own_tcga_only = setdiff(own_tcga, literature_genes),
  own_literature_only = setdiff(own_literature, tcga_genes),
  tcga_literature_only = setdiff(tcga_literature, own_genes),
  three_way_intersection = venn_intersection
)

# ===================== 6. 导出全部基因csv文件 =====================
## 6.1 三组原始全部基因
write.csv(data.frame(Gene = own_genes), file.path(outdir, "set_own_nonsyn_genes.csv"), row.names = FALSE)
write.csv(data.frame(Gene = tcga_genes), file.path(outdir, "set_TCGA_STAD_mutated_genes.csv"), row.names = FALSE)
write.csv(data.frame(Gene = literature_genes), file.path(outdir, "set_literature_TME_genes.csv"), row.names = FALSE)

## 6.2 韦恩图7个分区各自基因清单
for (nm in names(region_lists)) {
  write.csv(
    data.frame(Gene = sort(region_lists[[nm]])),
    file.path(outdir, paste0("venn_region_", nm, ".csv")),
    row.names = FALSE
  )
}

## 6.3 汇总计数表格（每组基因数量）
summary_df <- data.frame(
  region = c(
    "Own non-synonymous",
    "TCGA-STAD",
    "Literature TME",
    "Own_and_TCGA",
    "Own_and_Literature",
    "TCGA_and_Literature",
    "Own_and_TCGA_and_Literature"
  ),
  count = c(
    length(own_genes),
    length(tcga_genes),
    length(literature_genes),
    length(own_tcga),
    length(own_literature),
    length(tcga_literature),
    length(venn_intersection)
  )
)
write.csv(summary_df, file.path(outdir, "venn_set_count_summary.csv"), row.names = FALSE)

# ===================== 7. 绘制三色韦恩图（低饱和粉/绿/蓝） =====================
# 低饱和科研三色：粉、绿、蓝
venn_pal <- c(
  "#D39AA2",  # 雾粉 Own non-synonymous
  "#8FB7A8",  # 鼠尾草绿 TCGA-STAD
  "#8096BC"   # 雾霾蓝 Literature TME
)

venn_plot <- ggVennDiagram(
  sets,
  set_color = venn_pal,
  label_alpha = 0,
  label = "count",
  edge_lty = 1,
  edge_size = 0.6
) +
  scale_fill_gradientn(colors = venn_pal) +
  labs(
    title = paste0("Selection using ", length(literature_genes), " literature-supported gastric TME genes"),
    subtitle = "Own non-synonymous mutations ∩ TCGA-STAD mutations ∩ literature-supported TME genes"
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey30"),
    legend.position = "none"
  )

# 导出矢量PDF（论文投稿）
# Override the default count-gradient plot with seven explicitly colored
# exclusive regions. Pairwise colors are 50:50 RGB mixtures; the central
# three-way intersection is white.
venn_data <- process_data(Venn(sets))
region_edge <- venn_regionedge(venn_data)
region_label <- venn_regionlabel(venn_data)
set_edge <- venn_setedge(venn_data)
set_label <- venn_setlabel(venn_data)

region_pal <- c(
  "1" = "#D39AA2",
  "2" = "#8FB7A8",
  "3" = "#8096BC",
  "1/2" = "#B1A9A5",
  "1/3" = "#AA98AF",
  "2/3" = "#88A7B2",
  "1/2/3" = "#FFFFFF"
)

venn_plot <- ggplot() +
  geom_polygon(
    data = region_edge,
    aes(x = X, y = Y, group = id, fill = id),
    color = NA
  ) +
  geom_path(
    data = set_edge,
    aes(x = X, y = Y, group = id, color = id),
    linewidth = 0.8
  ) +
  geom_text(
    data = region_label,
    aes(x = X, y = Y, label = count),
    size = 5,
    color = "black"
  ) +
  geom_text(
    data = set_label,
    aes(x = X, y = Y, label = name, color = id),
    size = 5,
    fontface = "bold"
  ) +
  scale_fill_manual(values = region_pal) +
  scale_color_manual(values = c(
    "1" = "#D39AA2",
    "2" = "#8FB7A8",
    "3" = "#8096BC"
  )) +
  coord_equal(clip = "off") +
  labs(
    title = paste0("Selection using ", length(literature_genes), " literature-supported gastric TME genes"),
    subtitle = "Own non-synonymous mutations / TCGA-STAD mutations / literature-supported TME genes"
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey30"),
    legend.position = "none",
    plot.margin = margin(20, 55, 25, 55)
  )

ggsave(
  file.path(outdir, "Figure3_gene_selection_venn_TME60_7region.pdf"),
  venn_plot,
  width = 7,
  height = 6,
  useDingbats = FALSE
)

# 导出300DPI高清PNG（PPT使用）
ggsave(
  file.path(outdir, "Figure3_gene_selection_venn_TME60_7region.png"),
  venn_plot,
  width = 7,
  height = 6,
  dpi = 300
)

# ===================== 8. 生成运行日志txt =====================
# 记录输入路径、各组基因数量、三者共有基因名称
writeLines(c(
  "Gastric Figure3 gene selection Venn run notes",
  paste0("Own maftools RDS: ", own_maf_rds_file),
  paste0("TCGA mutation matrix: ", tcga_matrix_file),
  paste0("Literature support file: ", literature_file),
  paste0("Own non-synonymous genes: ", length(own_genes)),
  paste0("TCGA-STAD mutated genes: ", length(tcga_genes)),
  paste0("Literature TME genes: ", length(literature_genes)),
  paste0("Three-way intersection genes count: ", length(venn_intersection)),
  paste0("Three-way intersection list: ", paste(sort(venn_intersection), collapse = ", "))
), file.path(outdir, "Figure3_gene_selection_venn_TME60_run_notes.txt"))

# ===================== 9. 控制台输出汇总信息 =====================
cat("Done.\n")
cat("所有输出文件存放目录：", outdir, "\n")
# 打印各组基因数量统计表
print(summary_df)
# 打印最终筛选得到的全部候选TME基因
cat("三组共有候选TME基因列表：", paste(sort(venn_intersection), collapse = ", "), "\n")
}

############################################################
## Panel B: Selected-gene oncoplot
############################################################
run_figure3b <- function(input_dir, output_root) {
# =======================================================================
# 用户自定义参数：通常只需要修改本区域
# =======================================================================
# 1. 按希望在图中显示的顺序填写基因名；可使用大小写混合，代码会自动转为大写。
# 2. 可以增加或删除任意基因，例如：custom_genes <- c("TP53", "ARID1A", "PIK3CA")。
# 3. 重复基因会自动去重，并保留第一次出现的位置。
custom_genes <- c(

  "ARID1B","ATM","BRCA2","C1S",
  "CCL16","CD72","CDH11","CDH6",
  "CLDN11","COL12A1","COL15A1","COL23A1",
  "EP300","EPHB2","EPHB6","ERBB4",
  "HSPG2","IL10RB","IL16","ITGA4"
)

#   "TTN", "MUC16", "DNAH5", "DST", "PLEC",
#   "SACS", "XIRP2", "LRRK2", "UBR5", "CUBN",
#   "DNAH7", "USH2A", "AHNAK2", "ADGRV1", "DNAH9"
# )

# 自定义图标题和输出文件前缀；请仅使用适合文件名的英文、数字或下划线。
custom_plot_title <- "Custom non-synonymous mutated genes in gastric cfDNA"
output_prefix <- "Figure3B_unknown_randomly_mixed_v3"

  outdir <- file.path(output_root, "B_selected_gene_oncoplot")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  intermediate_dir <- file.path(outdir, "intermediate_data")
  dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)
  own_maf_rds_file <- file.path(input_dir, "gastric_maf_with_clinical.rds")
  tcga_matrix_file <- file.path(input_dir, "oncoplot_matrix_TCGA_STAD.txt")
norm_gene <- function(x) {
  toupper(trimws(as.character(x)))
}

clean_unknown <- function(x) {
  x <- trimws(as.character(x))
  unknown_tokens <- c("", "-", "/", ".", "na", "n/a", "unknown", "unkown")
  x[is.na(x) | tolower(x) %in% unknown_tokens] <- "Unknown"
  x
}

read_binary_matrix <- function(path) {
  mat <- fread(path, data.table = FALSE, check.names = FALSE)
  row_ids <- mat[[1]]
  mat <- mat[, -1, drop = FALSE]
  rownames(mat) <- row_ids
  mat[] <- lapply(mat, function(x) as.numeric(as.character(x)))
  mat[is.na(mat)] <- 0
  as.matrix(mat)
}

freq_from_long <- function(mut_df, all_samples, cohort_name) {
  mut_df %>%
    distinct(Tumor_Sample_Barcode, Gene) %>%
    count(Gene, name = "mutated_samples") %>%
    mutate(
      total_samples = length(unique(all_samples)),
      mutation_frequency = mutated_samples / total_samples * 100,
      cohort = cohort_name
    )
}

freq_from_matrix <- function(mat, cohort_name) {
  tibble::tibble(
    Gene = norm_gene(colnames(mat)),
    mutated_samples = colSums(mat > 0, na.rm = TRUE),
    total_samples = nrow(mat),
    mutation_frequency = mutated_samples / total_samples * 100,
    cohort = cohort_name
  ) %>%
    filter(mutated_samples > 0)
}

maf_sub <- readRDS(own_maf_rds_file)

own_mut <- as.data.frame(maf_sub@data) %>%
  mutate(
    Tumor_Sample_Barcode = as.character(Tumor_Sample_Barcode),
    Gene = norm_gene(Hugo_Symbol)
  ) %>%
  filter(!is.na(Gene), Gene != "")

own_clinical <- as.data.frame(maf_sub@clinical.data) %>%
  mutate(across(everything(), clean_unknown)) %>%
  distinct(Tumor_Sample_Barcode, .keep_all = TRUE)

# Use the cleaned annotations in the plot as well, so every spelling/blank form
# of unknown is displayed consistently as "Unknown".
maf_sub@clinical.data <- data.table::as.data.table(own_clinical)

tcga_matrix <- read_binary_matrix(tcga_matrix_file)
colnames(tcga_matrix) <- norm_gene(colnames(tcga_matrix))

# 标准化自定义基因列表：去除空值、转为大写、去重并生成绘图顺序 Rank。
selected_gene_df <- tibble::tibble(Gene = norm_gene(custom_genes)) %>%
  filter(!is.na(Gene), Gene != "") %>%
  distinct(Gene, .keep_all = TRUE) %>%
  mutate(Rank = row_number()) %>%
  select(Rank, Gene)

if (nrow(selected_gene_df) == 0) {
  stop("custom_genes 为空；请至少填写一个用于作图的基因。")
}

own_all_samples <- unique(own_clinical$Tumor_Sample_Barcode)
own_mutated_genes <- sort(unique(own_mut$Gene))
tcga_mutated_genes <- sort(colnames(tcga_matrix)[colSums(tcga_matrix > 0, na.rm = TRUE) > 0])
own_freq <- freq_from_long(own_mut, own_all_samples, "Own cfDNA")
tcga_freq <- freq_from_matrix(tcga_matrix, "TCGA-STAD")

# 使用用户自定义基因作为唯一作图基因集，顺序与 custom_genes 中的顺序一致。
selected_genes <- selected_gene_df$Gene

# 检查所选基因是否能在本地 MAF 和 TCGA 矩阵中找到。
missing_in_own <- setdiff(selected_genes, own_mutated_genes)
missing_in_tcga <- setdiff(selected_genes, tcga_mutated_genes)
if (length(missing_in_own) > 0) {
  warning("以下基因未在本地 MAF 突变数据中找到：", paste(missing_in_own, collapse = ", "))
}
if (length(missing_in_tcga) > 0) {
  warning("以下基因未在 TCGA 突变矩阵中找到：", paste(missing_in_tcga, collapse = ", "))
}

intersection_stats <- full_join(
  own_freq %>% filter(Gene %in% selected_genes) %>%
    select(Gene, own_mutated_samples = mutated_samples, own_total_samples = total_samples, own_mutation_frequency = mutation_frequency),
  tcga_freq %>% filter(Gene %in% selected_genes) %>%
    select(Gene, tcga_mutated_samples = mutated_samples, tcga_total_samples = total_samples, tcga_mutation_frequency = mutation_frequency),
  by = "Gene"
) %>%
  right_join(selected_gene_df, by = "Gene") %>%
  mutate(
    across(where(is.numeric), ~ ifelse(is.na(.), 0, .)),
    combined_frequency = own_mutation_frequency + tcga_mutation_frequency
  ) %>%
  arrange(Rank)

# 仅将确实存在于本地 MAF 中的所选基因交给 maftools，顺序仍按 Rank 保持。
genes_oncoplot <- selected_genes[selected_genes %in% own_mutated_genes]
if (length(genes_oncoplot) == 0) {
  stop("custom_genes 中没有任何基因存在于本地 MAF 突变数据，无法绘制 oncoplot。")
}

# 保存自定义基因对应的本地非同义突变明细，便于逐条复核瀑布图中的方格。
selected_mutation_records <- own_mut %>%
  filter(Gene %in% selected_genes) %>%
  arrange(match(Gene, selected_genes), Tumor_Sample_Barcode)

# 保存所有关键中间数据。编号反映数据处理顺序。
write.csv(own_clinical, file.path(intermediate_dir, "01_own_clinical_unknown_filled.csv"), row.names = FALSE)
write.csv(own_freq, file.path(intermediate_dir, "02_own_nonsyn_gene_mutation_frequency.csv"), row.names = FALSE)
write.csv(tcga_freq, file.path(intermediate_dir, "03_TCGA_STAD_gene_mutation_frequency.csv"), row.names = FALSE)
write.csv(selected_gene_df, file.path(intermediate_dir, "04_custom_gene_list_input.csv"), row.names = FALSE)
write.csv(data.frame(Gene = missing_in_own), file.path(intermediate_dir, "05_selected_genes_missing_in_own_MAF.csv"), row.names = FALSE)
write.csv(data.frame(Gene = missing_in_tcga), file.path(intermediate_dir, "06_selected_genes_missing_in_TCGA.csv"), row.names = FALSE)
write.csv(intersection_stats, file.path(intermediate_dir, "07_custom_gene_statistics.csv"), row.names = FALSE)
write.csv(selected_mutation_records, file.path(intermediate_dir, "08_custom_gene_mutation_records.csv"), row.names = FALSE)
write.csv(data.frame(Rank = seq_along(genes_oncoplot), Gene = genes_oncoplot),
          file.path(outdir, "genes_used_in_Figure3B_oncoplot.csv"), row.names = FALSE)

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

# Unknown is always neutral gray in every clinical annotation track.
ann_colors <- lapply(
  ann_colors,
  function(x) c(x, "Unknown" = "#ECF0F1")
)
#
# ann_colors <- list(
#   Age = c("<60" = "#A6CEE3", ">=60" = "#1F78B4", "Unknown" = "#D9D9D9"),
#   Gender = c("Male" = "#4DBBD5", "Female" = "#E64B35", "Unknown" = "#D9D9D9"),
#   pT = c("T0" = "#F7F7F7", "T1" = "#D9F0D3", "T1a" = "#C7E9C0", "T1b" = "#A1D99B", "T2" = "#74C476", "T3" = "#31A354", "T4" = "#006D2C", "T4a" = "#005A32", "T4b" = "#00441B", "Unknown" = "#D9D9D9"),
#   pN = c("N0" = "#F0F0F0", "N1" = "#BCBDDC", "N2" = "#807DBA", "N3" = "#54278F", "N3a" = "#3F007D", "N3b" = "#2D004B", "Unknown" = "#D9D9D9"),
#   M = c("M0" = "#BDBDBD", "M1" = "#E41A1C", "Unknown" = "#D9D9D9"),
#   Stage = c("Ia" = "#D9F0A3", "Ib" = "#ADDD8E", "IIa" = "#78C679", "IIb" = "#41AB5D", "IIIa" = "#238443", "IIIb" = "#006837", "IIIc" = "#004529", "IV" = "#08519C", "Unknown" = "#D9D9D9"),
#   Metastasis = c("CY" = "#FDB462", "LYM" = "#80B1D3", "LYM/CY" = "#B3DE69", "LYM/H" = "#FB8072", "LYM/P" = "#BEBADA", "Unknown" = "#D9D9D9"),
#   Lauren = c("Intestinal" = "#66C2A5", "Diffuse" = "#FC8D62", "Mixed" = "#8DA0CB", "0" = "#BDBDBD", "Unknown" = "#D9D9D9"),
#   HER2 = c("HER2_0" = "#F0F0F0", "HER2_1" = "#BDD7E7", "HER2_2" = "#6BAED6", "HER2_3" = "#2171B5", "Unknown" = "#D9D9D9"),
#   Ki67 = c("Low" = "#9ECAE1", "High" = "#DE2D26", "Unknown" = "#D9D9D9")
# )

clinical_features <- c("Age", "Gender", "pT", "pN", "M", "Stage", "Metastasis", "Lauren", "HER2", "Ki67")
clinical_features_use <- clinical_features[
  sapply(clinical_features, function(x) {
    vals <- maf_sub@clinical.data[[x]]
    vals <- vals[!is.na(vals) & vals != ""]
    length(unique(vals)) > 0
  })
]

# Build a deterministic sample order that spreads samples containing any
# Unknown clinical annotation evenly among samples with complete annotations.
# The normalized positions make this work even when the two groups differ in size.
clinical_for_order <- own_clinical %>%
  filter(Tumor_Sample_Barcode %in% own_all_samples) %>%
  mutate(
    unknown_count = apply(
      select(., all_of(clinical_features_use)),
      1,
      function(x) sum(tolower(trimws(as.character(x))) %in% c("unknown", "unkown"))
    ),
    clinical_group = case_when(
      unknown_count == 0 ~ "Complete",
      unknown_count == length(clinical_features_use) ~ "All_unknown",
      TRUE ~ "Partial_unknown"
    )
  )

# Randomly shuffle all samples together. A fixed seed makes the result
# reproducible while retaining a genuinely randomized sample arrangement.
set.seed(20260829)
sample_order_df <- clinical_for_order %>%
  slice_sample(prop = 1) %>%
  mutate(random_position = row_number())

mixed_sample_order <- sample_order_df$Tumor_Sample_Barcode
write.csv(
  select(sample_order_df, Tumor_Sample_Barcode, clinical_group,
         unknown_count, random_position),
  file.path(intermediate_dir, "09_randomly_mixed_sample_order.csv"),
  row.names = FALSE
)

ann_colors_use <- ann_colors[clinical_features_use]
for (feature in names(ann_colors_use)) {
  values <- sort(unique(as.character(maf_sub@clinical.data[[feature]])))
  values <- values[!is.na(values) & values != ""]
  missing_values <- setdiff(values, names(ann_colors_use[[feature]]))
  if (length(missing_values) > 0) {
    extra <- grDevices::hcl.colors(length(missing_values), "Set 3")
    names(extra) <- missing_values
    ann_colors_use[[feature]] <- c(ann_colors_use[[feature]], extra)
  }
}

# 根据实际作图基因数自动调整图高；基因越多，图形高度越大。
plot_height <- max(8, 6 + length(genes_oncoplot) * 0.25)
pdf(file.path(outdir, paste0(output_prefix, "_oncoplot.pdf")),
    width = 14, height = plot_height, useDingbats = FALSE)
oncoplot(
  maf = maf_sub,
  genes = genes_oncoplot,
  sampleOrder = mixed_sample_order,
  clinicalFeatures = clinical_features_use,
  colors = variant_colors,
  annotationColor = ann_colors_use,
  sortByAnnotation = FALSE,
  sortByMutation = FALSE,
  showTumorSampleBarcodes = FALSE,
  drawRowBar = TRUE,
  drawColBar = TRUE,
  legend_height = 9,
  fontSize = 0.75,
  legendFontSize = 0.6,
  annotationFontSize = 0.7,
  removeNonMutated = FALSE,
  titleText = custom_plot_title
)
dev.off()

writeLines(c(
  "Gastric Figure3B custom-gene non-synonymous oncoplot run notes",
  paste0("Own maftools RDS: ", own_maf_rds_file),
  paste0("TCGA mutation matrix: ", tcga_matrix_file),
  paste0("User-specified genes: ", paste(selected_genes, collapse = ", ")),
  paste0("Own clinical samples retained: ", length(own_all_samples)),
  paste0("Own mutated samples in maftools data: ", length(unique(own_mut$Tumor_Sample_Barcode))),
  paste0("Own non-synonymous mutated genes: ", length(own_mutated_genes)),
  paste0("TCGA mutated genes: ", length(tcga_mutated_genes)),
  paste0("Selected genes in input file: ", length(selected_genes)),
  paste0("Selected genes missing in own MAF: ", paste(missing_in_own, collapse = ", ")),
  paste0("Selected genes missing in TCGA: ", paste(missing_in_tcga, collapse = ", ")),
  paste0("Genes plotted: ", paste(genes_oncoplot, collapse = ", "))
), file.path(outdir, paste0(output_prefix, "_run_notes.txt")))

cat("Done.\n")
cat("Output directory:", outdir, "\n")
cat("User-specified genes:", length(selected_genes), "\n")
cat("Genes plotted:", paste(genes_oncoplot, collapse = ", "), "\n")

# 将本次实际运行的 R 脚本复制到输出目录，确保结果、参数和代码一起归档。
args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
if (length(script_arg) > 0) {
  current_script <- normalizePath(sub("^--file=", "", script_arg[1]),
                                  winslash = "/", mustWork = TRUE)
  archived_script <- normalizePath(
    file.path(outdir, paste0(output_prefix, "_code.R")),
    winslash = "/", mustWork = FALSE
  )
  if (!identical(current_script, archived_script)) {
    file.copy(current_script, archived_script, overwrite = TRUE)
  }
}
}

############################################################
## Panel C: Lauren-subtype top-five mutation dot plot
############################################################
run_figure3c <- function(input_dir, output_root) {
  outdir <- file.path(output_root, "C_lauren_top5")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  own_maf_rds_file <- file.path(input_dir, "gastric_maf_with_clinical.rds")
norm_gene <- function(x) {
  toupper(trimws(as.character(x)))
}

clean_unknown <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "-", "/", ".", "NA", "N/A", "na", "n/a")] <- "Unknown"
  x[is.na(x)] <- "Unknown"
  x
}

maf_sub <- readRDS(own_maf_rds_file)

target_genes <- c(
  "ARID1B", "ATM", "BRCA2", "C1S",
  "CCL16", "CD72", "CDH11", "CDH6",
  "CLDN11", "COL12A1", "COL15A1", "COL23A1",
  "EP300", "EPHB2", "EPHB6", "ERBB4",
  "HSPG2", "IL10RB", "IL16", "ITGA4"
) %>% norm_gene()

mut_data <- as.data.frame(maf_sub@data) %>%
  transmute(
    sample_id = as.character(Tumor_Sample_Barcode),
    Gene = norm_gene(Hugo_Symbol)
  ) %>%
  filter(!is.na(sample_id), !is.na(Gene), Gene %in% target_genes)

clinical_data <- as.data.frame(maf_sub@clinical.data) %>%
  transmute(
    sample_id = as.character(Tumor_Sample_Barcode),
    Lauren = clean_unknown(Lauren)
  ) %>%
  mutate(
    Lauren = case_when(
      Lauren %in% c("Intestinal", "Diffuse", "Mixed") ~ Lauren,
      TRUE ~ "Unknown"
    )
  ) %>%
  distinct(sample_id, .keep_all = TRUE)

sample_burden_df <- as.data.frame(maf_sub@data) %>%
  transmute(
    sample_id = as.character(Tumor_Sample_Barcode),
    Gene = norm_gene(Hugo_Symbol)
  ) %>%
  filter(!is.na(sample_id), !is.na(Gene), Gene != "") %>%
  group_by(sample_id) %>%
  summarise(
    mutation_count = n(),
    mutated_genes = n_distinct(Gene),
    .groups = "drop"
  ) %>%
  mutate(
    estimated_TMB = mutation_count / 38,
    mutation_burden = mutation_count
  )

lauren_samples <- clinical_data %>%
  filter(Lauren %in% c("Intestinal", "Diffuse", "Mixed"))

lauren_n_df <- lauren_samples %>%
  count(Lauren, name = "total_samples")

top5_lauren_df <- mut_data %>%
  left_join(lauren_samples, by = "sample_id") %>%
  left_join(sample_burden_df, by = "sample_id") %>%
  filter(Lauren %in% c("Intestinal", "Diffuse", "Mixed")) %>%
  distinct(sample_id, Lauren, Gene, estimated_TMB, mutation_burden) %>%
  group_by(Lauren, Gene) %>%
  summarise(
    mutated_samples = n_distinct(sample_id),
    median_estimated_TMB = median(estimated_TMB, na.rm = TRUE),
    mean_estimated_TMB = mean(estimated_TMB, na.rm = TRUE),
    median_mutation_burden = median(mutation_burden, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(lauren_n_df, by = "Lauren") %>%
  mutate(
    mutation_frequency = mutated_samples / total_samples * 100
  ) %>%
  filter(mutated_samples >= 1) %>%
  group_by(Lauren) %>%
  arrange(
    desc(mutation_frequency),
    desc(mutated_samples),
    desc(median_estimated_TMB),
    Gene,
    .by_group = TRUE
  ) %>%
  slice_head(n = 5) %>%
  ungroup()

lauren_order <- c("Intestinal", "Diffuse", "Mixed")
gene_order <- top5_lauren_df %>%
  mutate(Lauren = factor(Lauren, levels = lauren_order)) %>%
  arrange(
    Lauren,
    desc(mutation_frequency),
    desc(mutated_samples),
    desc(median_estimated_TMB),
    Gene
  ) %>%
  pull(Gene) %>%
  unique()

plot_df <- top5_lauren_df %>%
  mutate(
    Lauren = factor(Lauren, levels = lauren_order),
    Gene = factor(Gene, levels = rev(gene_order))
  )

write.csv(data.frame(Gene = target_genes), file.path(outdir, "input_20_target_genes.csv"), row.names = FALSE)
write.csv(clinical_data, file.path(outdir, "own_clinical_lauren_unknown_filled.csv"), row.names = FALSE)
write.csv(lauren_n_df, file.path(outdir, "lauren_subtype_sample_counts.csv"), row.names = FALSE)
write.csv(sample_burden_df, file.path(outdir, "sample_mutation_burden_estimated_TMB.csv"), row.names = FALSE)
write.csv(top5_lauren_df, file.path(outdir, "Figure3C_lauren_top5_target20_dotplot_data.csv"), row.names = FALSE)

p_lauren_top5 <- ggplot(plot_df, aes(x = Lauren, y = Gene)) +
  geom_point(
    aes(size = mutation_frequency, color = median_estimated_TMB),
    alpha = 0.9
  ) +
  scale_size_continuous(
    range = c(3, 11),
    name = "Mutation frequency (%)"
  ) +
  scale_color_gradient(
    low = "#9ECAE1",
    high = "#D73027",
    name = "Median estimated TMB"
  ) +
  labs(
    x = "Lauren subtype",
    y = NULL,
    title = "Top 5 mutated target genes across Lauren subtypes",
    subtitle = "Restricted to the 20 specified target genes"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey30"),
    axis.text.x = element_text(color = "black", face = "bold"),
    axis.text.y = element_text(color = "black", face = "italic"),
    legend.position = "right"
  )

ggsave(
  file.path(outdir, "Figure3C_lauren_top5_target20_dotplot.pdf"),
  p_lauren_top5,
  width = 7.5,
  height = 5.5,
  useDingbats = FALSE
)

ggsave(
  file.path(outdir, "Figure3C_lauren_top5_target20_dotplot.png"),
  p_lauren_top5,
  width = 7.5,
  height = 5.5,
  dpi = 300
)

writeLines(c(
  "Gastric Figure3C Lauren Top5 analysis of 20 target genes",
  paste0("Own maftools RDS: ", own_maf_rds_file),
  paste0("Target genes: ", paste(target_genes, collapse = ", ")),
  paste0("Lauren subtypes used: ", paste(lauren_order, collapse = ", ")),
  paste0("Lauren sample counts: ", paste(lauren_n_df$Lauren, lauren_n_df$total_samples, sep = "=", collapse = "; ")),
  paste0("Plotted rows: ", nrow(top5_lauren_df))
), file.path(outdir, "Figure3C_lauren_top5_target20_run_notes.txt"))

cat("Done.\n")
cat("Output directory:", outdir, "\n")
cat("Lauren sample counts:\n")
print(lauren_n_df)
cat("Plotted rows:", nrow(top5_lauren_df), "\n")
}

############################################################
## Panel D: GO-BP enrichment
############################################################
run_figure3d <- function(output_root) {
  outdir <- file.path(output_root, "D_GO_BP")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
norm_gene <- function(x) {
  toupper(trimws(as.character(x)))
}

ratio_to_numeric <- function(x) {
  if (is.numeric(x)) {
    return(x)
  }
  vapply(strsplit(as.character(x), "/"), function(z) {
    if (length(z) == 2) {
      as.numeric(z[1]) / as.numeric(z[2])
    } else {
      as.numeric(z[1])
    }
  }, numeric(1))
}

genes <- c(
  "ARID1B", "ATM", "BRCA2", "C1S",
  "CCL16", "CD72", "CDH11", "CDH6",
  "CLDN11", "COL12A1", "COL15A1", "COL23A1",
  "EP300", "EPHB2", "EPHB6", "ERBB4",
  "HSPG2", "IL10RB", "IL16", "ITGA4"
) %>% norm_gene()

gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = genes,
  keytype = "SYMBOL",
  columns = c("SYMBOL", "ENTREZID", "GENENAME")
) %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(SYMBOL, ENTREZID, .keep_all = TRUE)

write.csv(data.frame(Gene = genes), file.path(outdir, "input_20_target_genes.csv"), row.names = FALSE)
write.csv(gene_map, file.path(outdir, "input_20_gene_symbol_to_entrez_mapping.csv"), row.names = FALSE)

if (nrow(gene_map) < 3) {
  stop("Fewer than three input genes could be mapped to Entrez IDs; GO enrichment is not reliable.")
}

ego_bp <- enrichGO(
  gene = unique(gene_map$ENTREZID),
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  readable = TRUE
)

go_all <- as.data.frame(ego_bp)

write.csv(go_all, file.path(outdir, "GO_BP_enrichment_all_results.csv"), row.names = FALSE)

if (nrow(go_all) == 0) {
  writeLines(c(
    "No GO-BP enrichment terms were returned by clusterProfiler for the 20-gene input.",
    paste0("Input genes: ", paste(genes, collapse = ", ")),
    paste0("Mapped genes: ", paste(gene_map$SYMBOL, collapse = ", "))
  ), file.path(outdir, "GO_BP_no_terms_returned.txt"))
  stop("No GO-BP enrichment terms returned. Intermediate files were saved.")
}

plot_df <- go_all %>%
  mutate(
    GeneRatio_num = ratio_to_numeric(GeneRatio),
    neg_log10_FDR = -log10(p.adjust)
  ) %>%
  arrange(p.adjust, desc(Count), desc(GeneRatio_num)) %>%
  slice_head(n = min(14, nrow(go_all))) %>%
  arrange(GeneRatio_num)

write.csv(plot_df, file.path(outdir, "GO_BP_enrichment_top_terms_for_plot.csv"), row.names = FALSE)

panel_d <- ggplot(
  plot_df,
  aes(
    x = GeneRatio_num,
    y = factor(Description, levels = Description)
  )
) +
  geom_point(
    aes(size = Count, color = neg_log10_FDR),
    alpha = 0.95
  ) +
  scale_color_gradient(
    low = "#2C5A99",
    high = "#9D2F2F",
    name = expression(-log[10](FDR))
  ) +
  scale_size_continuous(name = "Gene count", range = c(3, 9)) +
  labs(
    x = "Gene ratio",
    y = NULL,
    title = "GO-BP enrichment of 20 target genes"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    plot.margin = margin(5, 10, 5, 5)
  )

ggsave(
  file.path(outdir, "Figure3D_GO_BP_target20_dotplot.pdf"),
  panel_d,
  width = 7.5,
  height = 5.5,
  useDingbats = FALSE
)

ggsave(
  file.path(outdir, "Figure3D_GO_BP_target20_dotplot.png"),
  panel_d,
  width = 7.5,
  height = 5.5,
  dpi = 300
)

writeLines(c(
  "Gastric Figure3D GO-BP analysis of 20 target genes",
  paste0("Input genes: ", paste(genes, collapse = ", ")),
  paste0("Mapped Entrez IDs: ", nrow(gene_map)),
  paste0("GO-BP terms returned: ", nrow(go_all)),
  paste0("GO-BP terms plotted: ", nrow(plot_df))
), file.path(outdir, "Figure3D_GO_target20_run_notes.txt"))

cat("Done.\n")
cat("Output directory:", outdir, "\n")
cat("Input genes:", paste(genes, collapse = ", "), "\n")
cat("GO-BP terms returned:", nrow(go_all), "\n")
cat("GO-BP terms plotted:", nrow(plot_df), "\n")
}

run_figure3a(input_dir, figure3_output_dir)
run_figure3b(input_dir, figure3_output_dir)
run_figure3c(input_dir, figure3_output_dir)
run_figure3d(figure3_output_dir)

message("Figure 3 panels A-D completed. Output directory: ", figure3_output_dir)
