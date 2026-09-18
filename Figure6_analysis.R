############################################################
## Figure 6: Independent survival validation in GEO cohorts
##
## Required directory layout:
##   ./data/Figure6/GEO/
##   ./data/Figure6/GEO2/
##
## Run with: Rscript Figure6_analysis.R
## Outputs are written under ./output/Figure6.
############################################################

required_packages <- c(
  "data.table", "dplyr", "tibble", "readxl", "survival", "survminer", "ggplot2"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "缺少R包：", paste(missing_packages, collapse = ", "),
    "。请先安装这些依赖后再运行。"
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
  library(readxl)
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
geo_dir <- file.path(script_dir, "data", "Figure6", "GEO")
geo2_dir <- file.path(script_dir, "data", "Figure6", "GEO2")
out_dir <- file.path(script_dir, "output", "Figure6")
intermediate_dir <- file.path(out_dir, "intermediate_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(geo_dir)) stop("GEO input folder does not exist: ", geo_dir)
if (!dir.exists(geo2_dir)) stop("GEO2 input folder does not exist: ", geo2_dir)

# Fixed model coefficients are recorded below. The Figure 5 parameter file is
# optional and is saved only as a reference when available.
manual_model <- tibble::tribble(
  ~Gene,    ~Coefficient,
  "COL15A1",    0.00337975891621139,
  "PLXNC1",    0.0254558844392741,
  "CDH11",    0.0365614804586186,
  "CCL16",    0.0148733152997758,
  "CLDN11",    0.133919085277019,
  "ERBB4",    0.0782012384327468,
  "CDH6",    0.163505577557928
)

complete_model_file <- file.path(
  script_dir, "output", "Figure5", "AC_prognostic_model",
  "06_LASSO_complete_model_parameters.csv"
)
coef_df <- manual_model %>%
  transmute(Gene = toupper(trimws(Gene)), Coefficient = as.numeric(Coefficient))
if (!nrow(coef_df)) stop("manual_model cannot be empty.")
if (any(is.na(coef_df$Gene) | coef_df$Gene == "")) stop("manual_model contains an empty gene symbol.")
if (anyDuplicated(coef_df$Gene)) {
  stop("Duplicated genes in manual_model: ",
       paste(unique(coef_df$Gene[duplicated(coef_df$Gene)]), collapse = ", "))
}
if (any(!is.finite(coef_df$Coefficient))) stop("All manual_model coefficients must be finite numbers.")
risk_coef <- setNames(coef_df$Coefficient, coef_df$Gene)
risk_genes <- names(risk_coef)

write.csv(coef_df, file.path(out_dir, "00_fixed_TCGA_LASSO_Cox_coefficients.csv"),
          row.names = FALSE)
# This file is retained only as an optional reference and is not used in scoring.
training_model_parameters <- if (file.exists(complete_model_file)) {
  read.csv(complete_model_file, check.names = FALSE)
} else {
  data.frame(Note = "Complete training-model parameters unavailable; manual model remains usable.")
}
write.csv(training_model_parameters,
          file.path(out_dir, "00_TCGA_complete_model_parameters_reference.csv"),
          row.names = FALSE)

# 本地配置为具有生存时间和事件、可用于KM/Cox的验证队列。
cohort_spec <- list(
  GSE13861 = list(
    matrix_file = file.path(geo2_dir, "GSE13861", "GSE13861_series_matrix.txt.gz"),
    annotation_file = file.path(geo_dir, "platform_annotation", "GPL6884-11607.txt"),
    platform = "GPL6884", clinical_source = "standardized_os_csv",
    clinical_file = file.path(geo_dir, "OS", "GSE13861_OS_clinical_standardized.csv")
  ),
  GSE14208 = list(
    matrix_file = file.path(geo2_dir, "GSE14210", "GSE14210_series_matrix.txt.gz"),
    annotation_file = file.path(geo_dir, "platform_annotation", "GPL571-17391.txt"),
    platform = "GPL571", clinical_source = "standardized_os_csv",
    clinical_file = file.path(geo_dir, "OS", "GSE14208_OS_clinical_standardized.csv")
  ),
  GSE26899 = list(
    matrix_file = file.path(geo2_dir, "GSE26899", "GSE26899_series_matrix.txt.gz"),
    annotation_file = file.path(geo_dir, "platform_annotation", "GPL6947-13512.txt"),
    platform = "GPL6947", clinical_source = "standardized_os_csv",
    clinical_file = file.path(geo_dir, "OS", "GSE26899_OS_clinical_standardized.csv")
  ),
  GSE26901 = list(
    matrix_file = file.path(geo2_dir, "GSE26901", "GSE26901_series_matrix.txt.gz"),
    annotation_file = file.path(geo_dir, "platform_annotation", "GPL6947-13512.txt"),
    platform = "GPL6947", clinical_source = "standardized_os_csv",
    clinical_file = file.path(geo_dir, "OS", "GSE26901_OS_clinical_standardized.csv")
  ),
  GSE29272 = list(
    matrix_file = file.path(geo2_dir, "GSE29272", "GSE29272_series_matrix.txt.gz"),
    annotation_file = file.path(geo_dir, "platform_annotation", "GPL96-57554.txt"),
    platform = "GPL96", clinical_source = "standardized_os_csv",
    clinical_file = file.path(geo_dir, "OS", "GSE29272_OS_clinical_standardized.csv")
  ),
  GSE57303 = list(
    matrix_file = file.path(geo_dir, "GSE57303", "GSE57303_series_matrix.txt.gz"),
    annotation_file = file.path(geo_dir, "platform_annotation", "GPL570-55999.txt"),
    platform = "GPL570", clinical_source = "standardized_os_csv",
    clinical_file = file.path(geo_dir, "OS", "GSE57303_OS_clinical_standardized.csv")
  ),
  GSE84437 = list(
    matrix_file = file.path(geo2_dir, "GSE84437（84433+84426）", "GSE84437_series_matrix.txt.gz"),
    annotation_file = file.path(geo_dir, "platform_annotation", "GPL6947-13512.txt"),
    platform = "GPL6947", clinical_source = "standardized_os_csv",
    clinical_file = file.path(geo_dir, "OS", "GSE84437_OS_clinical_standardized.csv")
  ),
  GSE15459 = list(
    matrix_file = file.path(geo_dir, "GSE15459", "GSE15459_series_matrix.txt.gz"),
    annotation_file = file.path(geo_dir, "platform_annotation", "GPL570-55999.txt"),
    platform = "GPL570",
    clinical_source = "outcome_xls",
    clinical_file = file.path(geo_dir, "GSE15459", "GSE15459_outcome.xls")
  ),
  GSE34942 = list(
    matrix_file = file.path(geo_dir, "GSE34942", "GSE34942_series_matrix.txt.gz"),
    annotation_file = file.path(geo_dir, "platform_annotation", "GPL570-55999.txt"),
    platform = "GPL570",
    clinical_source = "outcome_xls",
    clinical_file = file.path(geo_dir, "GSE34942", "GSE34942_outcome.xls")
  ),
  # GSE84437 = list(
  #   matrix_file = file.path(geo_dir, "GSE84437", "GSE84437_series_matrix.txt.gz"),
  #   annotation_file = file.path(geo_dir, "platform_annotation", "GPL6947-13512.txt"),
  #   platform = "GPL6947",
  #   clinical_source = "series_matrix",
  #   clinical_file = NA_character_
  # ),
  GSE62254 = list(
    matrix_file = file.path(geo_dir, "GSE62254", "GSE62254_series_matrix.txt.gz"),
    annotation_file = file.path(geo_dir, "platform_annotation", "GPL570-55999.txt"),
    platform = "GPL570",
    clinical_source = "acrg_supplement_xls",
    clinical_file = file.path(
      geo_dir, "GSE62254", "41591_2015_BFnm3850_MOESM34_ESM.xls"
    )
  ),
  # GSE38749的Series Matrix直接提供OS月数与alive/death状态。
  GSE38749 = list(
    matrix_file = file.path(geo_dir, "GSE38749", "GSE38749_series_matrix.txt.gz"),
    annotation_file = file.path(geo_dir, "platform_annotation", "GPL570-55999.txt"),
    platform = "GPL570",
    clinical_source = "gse38749_series_matrix",
    clinical_file = NA_character_
  )
  # GSE26253 provides recurrence-free survival (RFS), not overall survival.
  # GSE26253 = list(
  #   matrix_file = file.path(geo_dir, "GSE26253", "GSE26253_series_matrix.txt.gz"),
  #   annotation_file = file.path(geo_dir, "platform_annotation", "GPL8432-11703.txt"),
  #   platform = "GPL8432",
  #   clinical_source = "gse26253_rfs_series_matrix",
  #   clinical_file = NA_character_,
  #   endpoint = "RFS"
  # )
)

# Default endpoint is OS. Alternative endpoints (DFS/RFS/DMFS/DSS/CSS) must be
# explicitly declared per cohort and are analysed/reported separately.
cohort_spec <- lapply(cohort_spec, function(x) {
  if (is.null(x$endpoint)) x$endpoint <- "OS"
  x
})

# 已下载但当前文件中缺少“生存时间+事件”的队列，保留排除记录。
excluded_cohorts <- tibble(
  Dataset = c("GSE26942", "GSE29272", "GSE57303", "GSE63089"),
  Reason = c(
    "Series Matrix含临床特征但无生存时间和事件",
    "肿瘤/癌旁表达队列，无生存时间和事件",
    "表达队列，当前本地文件无生存时间和事件",
    "配对肿瘤/正常表达队列，无生存时间和事件"
  )
)
# 2026-08-12 audit: replace the historical hard-coded list with every local cohort
# that cannot be used for survival validation. Reasons are based on inspected files.
excluded_cohorts <- tibble(
  Dataset = c("GSE13861", "GSE183136", "GSE22377", "GSE26253", "GSE26899",
              "GSE26901", "GSE26942", "GSE28541", "GSE37023", "GSE47007",
              "GSE51105", "GSE57303", "GSE66229"),
  Reason = c(
    "Local clinical material lacks a complete usable OS time/event pair",
    "New cohort audit: Series Matrix contains stage/age/sex but no OS time/event pair",
    "New cohort audit: Series Matrix contains histological type but no OS time/event pair",
    "Series Matrix lacks a complete OS time/event pair",
    "Clinical files contain pathology/treatment data but no complete OS time/event pair",
    "New cohort: pathology and chemotherapy fields only; no OS time/event pair",
    "Series Matrix lacks a complete OS time/event pair",
    "Series Matrix lacks a complete OS time/event pair",
    "New cohort audit: GPL96 matrix contains tissue/stroma fields but no OS time/event pair",
    "New cohort audit: GPL8300 matrix contains tissue type only; no OS time/event pair",
    "No complete OS time/event pair reliably matched to expression samples",
    "New/updated cohort: tissue and phenotype only; no OS time/event pair",
    "New cohort: tissue and patient identifiers only; no OS time/event pair"
  )
)
# GSE26253 is now evaluated using its complete RFS time/event pair.
excluded_cohorts <- excluded_cohorts %>% filter(Dataset != "GSE26253")
excluded_cohorts <- excluded_cohorts %>%
  filter(!Dataset %in% c("GSE13861", "GSE26899", "GSE26901", "GSE29272",
                         "GSE57303", "GSE84437"))
write.csv(excluded_cohorts, file.path(out_dir, "01_GEO_cohorts_excluded_from_survival.csv"),
          row.names = FALSE)

required_files <- unlist(lapply(cohort_spec, function(x) {
  c(x$matrix_file, x$annotation_file,
    if (!is.na(x$clinical_file)) x$clinical_file else character())
}))
if (!all(file.exists(required_files))) {
  stop("以下输入文件不存在：\n", paste(required_files[!file.exists(required_files)], collapse = "\n"))
}

# =======================================================================
# 2. 通用文本处理函数
# =======================================================================
strip_quote <- function(x) gsub('^"|"$', "", x)

normalize_key <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("^_|_$", "", x)
}

# 从Series Matrix头部提取每个GSM样本的注释。
read_series_clinical <- function(path, dataset) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, open = "rt") else file(path, open = "rt")
  on.exit(close(con), add = TRUE)
  records <- list()
  sample_count <- 0L

  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0 || line == "!series_matrix_table_begin") break
    if (!grepl("^!Sample_", line)) next

    parts <- strip_quote(strsplit(line, "\t", fixed = TRUE)[[1]])
    if (length(parts) < 2) next
    tag <- normalize_key(sub("^!Sample_", "", parts[1]))
    values <- parts[-1]

    if (sample_count == 0L) {
      sample_count <- length(values)
      records <- replicate(sample_count, list(dataset = dataset), simplify = FALSE)
    }
    if (length(values) < sample_count) {
      values <- c(values, rep(NA_character_, sample_count - length(values)))
    }
    values <- values[seq_len(sample_count)]

    for (i in seq_len(sample_count)) {
      value <- values[i]
      if (tag == "characteristics_ch1" && grepl(":", value, fixed = TRUE)) {
        key <- normalize_key(sub(":.*$", "", value))
        val <- trimws(sub("^[^:]*:\\s*", "", value))
        records[[i]][[key]] <- val
      } else if (tag %in% c("geo_accession", "title", "source_name_ch1", "platform_id")) {
        records[[i]][[tag]] <- value
      }
    }
  }

  keys <- unique(unlist(lapply(records, names)))
  out <- as.data.frame(lapply(keys, function(k) {
    vapply(records, function(r) if (!is.null(r[[k]])) r[[k]] else NA_character_, character(1))
  }), stringsAsFactors = FALSE)
  names(out) <- keys
  out
}

# 读取本地GPL注释，只保留模型基因的“探针-基因”映射。
read_probe_gene_map <- function(annotation_file, platform, genes) {
  if (platform == "GPL570") {
    annot <- fread(annotation_file, sep = "\t", header = TRUE, skip = "ID\t",
                   select = c("ID", "Gene Symbol"), data.table = FALSE,
                   quote = "", fill = TRUE)
    names(annot) <- c("probe_id", "symbol_raw")
  } else if (platform %in% c("GPL6947", "GPL8432", "GPL6884")) {
    annot <- fread(annotation_file, sep = "\t", header = TRUE, skip = "ID\t",
                   select = c("ID", "Symbol"), data.table = FALSE,
                   quote = "", fill = TRUE)
    names(annot) <- c("probe_id", "symbol_raw")
  } else if (platform %in% c("GPL571", "GPL96")) {
    annot <- fread(annotation_file, sep = "\t", header = TRUE, skip = "ID\t",
                   select = c("ID", "Gene Symbol"), data.table = FALSE,
                   quote = "", fill = TRUE)
    names(annot) <- c("probe_id", "symbol_raw")
  } else {
    stop("尚未配置平台注释规则：", platform)
  }

  # 旧版芯片注释仍将 ADGRV1 标为历史名称 GPR98/VLGR1。
  # 先把历史名称统一成当前模型使用的官方符号 ADGRV1，再匹配模型基因；
  # 这只是基因符号同义词转换，不改变基因、探针表达值或模型系数。
  symbol_aliases <- c(GPR98 = "ADGRV1", VLGR1 = "ADGRV1")

  # 一个探针可能对应多个Symbol，拆分后只保留固定模型基因。
  rows <- lapply(seq_len(nrow(annot)), function(i) {
    raw <- annot$symbol_raw[i]
    if (is.na(raw) || raw == "" || raw == "---") return(NULL)
    symbols <- toupper(trimws(unlist(strsplit(raw, " /// |;|,"))))
    alias_hit <- symbols %in% names(symbol_aliases)
    symbols[alias_hit] <- unname(symbol_aliases[symbols[alias_hit]])
    symbols <- intersect(symbols, genes)
    if (!length(symbols)) return(NULL)
    data.frame(probe_id = annot$probe_id[i], Gene = symbols)
  })
  result <- bind_rows(rows) %>% distinct(probe_id, Gene)
  if (!nrow(result)) stop(platform, "注释中未找到模型基因探针。")
  result
}

# 只读取模型基因对应的探针表达行，避免加载完整大型表达矩阵。
read_series_expression_subset <- function(path, probe_ids) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, open = "rt") else file(path, open = "rt")
  on.exit(close(con), add = TRUE)
  header <- NULL

  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (line == "!series_matrix_table_begin") {
      header <- strip_quote(strsplit(readLines(con, n = 1), "\t", fixed = TRUE)[[1]])
      break
    }
  }
  if (is.null(header)) stop("Series Matrix中找不到表达表：", path)
  samples <- header[-1]
  keep <- setNames(rep(TRUE, length(probe_ids)), probe_ids)
  values <- list()

  repeat {
    lines <- readLines(con, n = 2000, warn = FALSE)
    if (!length(lines)) break
    for (line in lines) {
      if (line == "!series_matrix_table_end") break
      fields <- strip_quote(strsplit(line, "\t", fixed = TRUE)[[1]])
      probe <- fields[1]
      if (is.na(keep[probe])) next
      v <- suppressWarnings(as.numeric(fields[-1]))
      if (length(v) == length(samples)) values[[probe]] <- v
    }
    if (any(lines == "!series_matrix_table_end")) break
  }
  if (!length(values)) stop("表达矩阵中没有找到目标探针：", path)
  mat <- do.call(rbind, values)
  colnames(mat) <- samples
  mat
}

# 一个基因对应多个探针时，选择跨样本IQR最大的探针作为代表。
collapse_probe_to_gene <- function(probe_expression, probe_map) {
  selected <- list()
  gene_values <- list()
  for (gene in unique(probe_map$Gene)) {
    probes <- intersect(probe_map$probe_id[probe_map$Gene == gene], rownames(probe_expression))
    if (!length(probes)) next
    if (length(probes) == 1) {
      chosen <- probes
    } else {
      probe_iqr <- apply(probe_expression[probes, , drop = FALSE], 1, IQR, na.rm = TRUE)
      chosen <- probes[which.max(probe_iqr)]
    }
    gene_values[[gene]] <- probe_expression[chosen, ]
    selected[[gene]] <- data.frame(Gene = gene, selected_probe = chosen,
                                   available_probe_count = length(probes))
  }
  mat <- do.call(rbind, gene_values)
  rownames(mat) <- names(gene_values)
  attr(mat, "selected_probes") <- bind_rows(selected)
  mat
}

# =======================================================================
# 3. 队列临床数据标准化
# =======================================================================
make_survival_data <- function(dataset, spec) {
  if (spec$clinical_source == "standardized_os_csv") {
    raw <- read.csv(spec$clinical_file, check.names = FALSE, stringsAsFactors = FALSE)
    required <- c("GSM", "OS_time_month", "OS_status")
    missing <- setdiff(required, names(raw))
    if (length(missing)) stop(dataset, " standardized OS file missing columns: ",
                              paste(missing, collapse = ", "))
    clinical <- raw %>% transmute(
      geo_accession = as.character(GSM),
      OS_time_months = suppressWarnings(as.numeric(OS_time_month)),
      OS_event = suppressWarnings(as.numeric(OS_status)),
      age = if ("age" %in% names(raw)) suppressWarnings(as.numeric(age)) else NA_real_,
      gender = if ("sex" %in% names(raw)) as.character(sex) else NA_character_,
      stage = if ("stage" %in% names(raw)) as.character(stage) else NA_character_,
      lauren = if ("lauren_class" %in% names(raw)) as.character(lauren_class) else NA_character_,
      subtype = NA_character_
    )
  } else if (spec$clinical_source == "outcome_xls") {
    raw <- read_excel(spec$clinical_file)
    clinical <- raw %>% transmute(
      geo_accession = as.character(`GSM ID`),
      OS_time_months = as.numeric(`Overall.Survival (Months)**`),
      OS_event = as.numeric(`Outcome (1=dead)`),
      age = as.numeric(Age_at_surgery),
      gender = as.character(Gender),
      stage = as.character(Stage),
      lauren = as.character(Laurenclassification),
      subtype = as.character(Subtype)
    )
  } else if (spec$clinical_source == "acrg_supplement_xls") {
    # GSE62254 的逐例临床资料位于原论文 Supplementary Data 1。
    # Sample Name 常写成 T107 或 Tr_107；统一去除符号并把开头 TR 规范为 T，
    # 再与 Series Matrix 中的样本标题对应，最后取得 GSM 编号。
    raw <- read_excel(spec$clinical_file, sheet = "FINAL")
    sample_map <- read_series_clinical(spec$matrix_file, dataset) %>%
      transmute(
        geo_accession = as.character(geo_accession),
        matrix_sample_name = as.character(title),
        match_key = toupper(gsub("[^A-Za-z0-9]", "", matrix_sample_name)),
        match_key = sub("^TR", "T", match_key),
        match_key = sub("^([0-9]+)T$", "T\\1", match_key)
      )

    raw_standardized <- raw %>% transmute(
      # 部分病例的 Sample Name 为空，但 Tumor ID 与 GEO 标题中的 T编号一致；
      # 因此仅在 Sample Name 缺失时，用 T + Tumor ID 补齐匹配名称。
      clinical_sample_name_original = as.character(`Sample\nName`),
      clinical_sample_name = ifelse(
        is.na(clinical_sample_name_original) | trimws(clinical_sample_name_original) == "",
        paste0("T", as.character(`Tumor ID`)),
        clinical_sample_name_original
      ),
      match_key = toupper(gsub("[^A-Za-z0-9]", "", clinical_sample_name)),
      match_key = sub("^TR", "T", match_key),
      match_key = sub("^([0-9]+)T$", "T\\1", match_key),
      OS_time_months = suppressWarnings(as.numeric(`OS\n(months)`)),
      followup_status = suppressWarnings(as.numeric(
        `FU status0=alive without ds, 1=alive with recurren ds, 2=dead without ds, 3=dead d/t recurrent ds, 4=dead, unknown, 5= FU loss`
      )),
      # 状态2/3/4均表示死亡；0/1/5在末次随访时按删失处理。
      OS_event = case_when(
        followup_status %in% c(2, 3, 4) ~ 1,
        followup_status %in% c(0, 1, 5) ~ 0,
        TRUE ~ NA_real_
      ),
      age = suppressWarnings(as.numeric(age)),
      gender = as.character(sex),
      stage = as.character(`pStage`),
      lauren = as.character(Lauren),
      subtype = as.character(`Mol. Subtype: 0=MSS/TP53-, 1=MSS/TP53+, 2 = MSI, 3= EMT`)
    )

    clinical <- raw_standardized %>%
      left_join(sample_map, by = "match_key") %>%
      dplyr::select(geo_accession, clinical_sample_name_original,
                    clinical_sample_name, matrix_sample_name,
                    OS_time_months, OS_event, followup_status,
                    age, gender, stage, lauren, subtype)
  } else if (spec$clinical_source == "gse38749_series_matrix") {
    # GSE38749 stores OS in months and status as alive/dead in sample annotations.
    raw <- read_series_clinical(spec$matrix_file, dataset)
    clinical <- raw %>% transmute(
      geo_accession = as.character(geo_accession),
      OS_time_months = suppressWarnings(as.numeric(time_months_overall_survival)),
      OS_event = case_when(
        # GEO uses "death" for events (rather than "dead") in this cohort.
        tolower(trimws(status)) %in% c("dead", "death", "deceased") ~ 1,
        tolower(trimws(status)) == "alive" ~ 0,
        TRUE ~ NA_real_
      ),
      age = suppressWarnings(as.numeric(age_y)),
      gender = as.character(gender),
      stage = as.character(tumor_stage_ajcc),
      lauren = as.character(histological_type),
      subtype = NA_character_
    )
  } else if (spec$clinical_source == "gse26253_rfs_series_matrix") {
    raw <- read_series_clinical(spec$matrix_file, dataset)
    clinical <- raw %>% transmute(
      geo_accession = as.character(geo_accession),
      # GEO definition: 0 = non-recurrence/censored, 1 = recurrence.
      OS_time_months = suppressWarnings(as.numeric(recurrence_free_survival_time_month)),
      OS_event = suppressWarnings(as.numeric(status_0_non_recurrence_1_recurrence)),
      age = NA_real_, gender = NA_character_, stage = NA_character_,
      lauren = NA_character_, subtype = NA_character_
    )
  } else if (spec$clinical_source == "auto_series_endpoint") {
    # Generic endpoint parser used only after the automatic audit has verified
    # an unambiguous time/event pair in the Series Matrix.
    raw <- read_series_clinical(spec$matrix_file, dataset)
    time_value <- suppressWarnings(as.numeric(raw[[spec$time_col]]))
    event_raw <- tolower(trimws(as.character(raw[[spec$event_col]])))
    event_value <- suppressWarnings(as.numeric(event_raw))
    if (all(is.na(event_value))) {
      event_value <- case_when(
        event_raw %in% c("dead", "death", "deceased", "died", "yes", "recurrence",
                         "recurred", "relapse", "progressed", "event") ~ 1,
        event_raw %in% c("alive", "living", "no", "non-recurrence", "no recurrence",
                         "disease free", "censored") ~ 0,
        TRUE ~ NA_real_
      )
    }
    clinical <- tibble(
      geo_accession = as.character(raw$geo_accession),
      OS_time_months = time_value * spec$time_to_months,
      OS_event = event_value,
      age = NA_real_, gender = NA_character_, stage = NA_character_,
      lauren = NA_character_, subtype = NA_character_
    )
  } else {
    raw <- read_series_clinical(spec$matrix_file, dataset)
    clinical <- raw %>% transmute(
      geo_accession = as.character(geo_accession),
      OS_time_months = suppressWarnings(as.numeric(duration_overall_survival)),
      OS_event = suppressWarnings(as.numeric(death)),
      age = suppressWarnings(as.numeric(age)),
      gender = as.character(sex),
      stage = NA_character_,
      lauren = NA_character_,
      subtype = NA_character_
    )
  }
  list(raw = raw, standardized = clinical)
}

# =======================================================================
# 3A. 自动遍历GEO目录并审计终点、疾病类型、平台与模型基因覆盖
# =======================================================================
find_one_series_matrix <- function(dataset_dir) {
  x <- list.files(dataset_dir, pattern = "_series_matrix\\.txt(\\.gz)?$",
                  recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  if (!length(x)) return(NA_character_)
  # Prefer compressed files, then the shortest path (avoids duplicate nested copy).
  x[order(!grepl("\\.gz$", x, ignore.case = TRUE), nchar(x))][1]
}

find_platform_annotation <- function(platform) {
  x <- list.files(file.path(geo_dir, "platform_annotation"),
                  pattern = paste0("^", platform, "(-|\\.).*\\.txt$|^", platform, "\\.txt$"),
                  full.names = TRUE, ignore.case = TRUE)
  if (length(x)) x[1] else NA_character_
}

# Returns one unambiguous endpoint pair. Specific endpoints are checked before
# generic OS/status fields to prevent recurrence status being mistaken for death.
detect_endpoint_pair <- function(raw) {
  nms <- names(raw)
  rules <- list(
    RFS = list(time = "recurrence.*(free.*)?survival.*time|rfs.*time|time.*recurr",
               event = "status.*recurr|recurr.*status|rfs.*event"),
    DFS = list(time = "disease.*free.*survival.*time|dfs.*time|time.*disease.*free",
               event = "disease.*free.*status|dfs.*event|status.*disease"),
    DMFS = list(time = "distant.*metastasis.*free.*survival.*time|dmfs.*time|time.*distant.*metasta",
                event = "distant.*metasta.*status|dmfs.*event|status.*distant.*metasta"),
    DSS = list(time = "disease.*specific.*survival.*time|cancer.*specific.*survival.*time|dss.*time|css.*time",
               event = "disease.*specific.*status|cancer.*specific.*status|dss.*event|css.*event"),
    OS = list(time = "overall.*survival.*(time|month|day|year)|duration.*overall.*survival|os.*time|time.*overall.*survival",
              event = "^(death|dead|vital_status|os_event|overall_survival_status)$|status.*overall.*survival")
  )
  partial <- character()
  for (ep in names(rules)) {
    tc <- nms[grepl(rules[[ep]]$time, nms, ignore.case = TRUE)]
    ec <- nms[grepl(rules[[ep]]$event, nms, ignore.case = TRUE)]
    if (length(tc) && length(ec)) {
      tv <- suppressWarnings(as.numeric(raw[[tc[1]]]))
      er <- tolower(trimws(as.character(raw[[ec[1]]])))
      ev <- suppressWarnings(as.numeric(er))
      if (all(is.na(ev))) ev <- ifelse(er %in% c("dead","death","deceased","died","yes","recurrence","recurred","relapse","event"), 1,
                                       ifelse(er %in% c("alive","living","no","non-recurrence","no recurrence","censored"), 0, NA))
      valid <- is.finite(tv) & tv > 0 & ev %in% c(0, 1)
      if (sum(valid) >= 10 && sum(ev[valid] == 1) >= 5) {
        unit <- if (grepl("day", tc[1], ignore.case = TRUE)) 1/30.4375 else
          if (grepl("year", tc[1], ignore.case = TRUE)) 12 else 1
        return(list(found = TRUE, endpoint = ep, time_col = tc[1], event_col = ec[1],
                    time_to_months = unit, valid_n = sum(valid), events = sum(ev[valid] == 1)))
      }
    }
    if (length(tc) && !length(ec)) partial <- c(partial, paste0(ep, ": missing event/status"))
    if (!length(tc) && length(ec)) partial <- c(partial, paste0(ep, ": missing time"))
  }
  list(found = FALSE, reason = if (length(partial)) paste(unique(partial), collapse = "; ") else
    "No OS/DFS/RFS/DMFS/DSS/CSS time or event field detected")
}

audit_local_geo <- function() {
  dirs <- list.dirs(geo_dir, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[grepl("^GSE[0-9]+$", basename(dirs), ignore.case = TRUE)]
  rows <- list(); auto_specs <- list()
  for (d in dirs) {
    ds <- toupper(basename(d)); mf <- find_one_series_matrix(d)
    if (is.na(mf)) {
      rows[[ds]] <- tibble(Dataset=ds, Audit_status="excluded", Endpoint=NA_character_,
                           Platform=NA_character_, Reason="No Series Matrix file found",
                           Missing_model_genes=NA_character_, Matrix_file=NA_character_)
      next
    }
    raw <- tryCatch(read_series_clinical(mf, ds), error = function(e) NULL)
    if (is.null(raw) || !nrow(raw)) {
      rows[[ds]] <- tibble(Dataset=ds, Audit_status="excluded", Endpoint=NA_character_,
                           Platform=NA_character_, Reason="Series Matrix sample metadata cannot be parsed",
                           Missing_model_genes=NA_character_, Matrix_file=mf)
      next
    }
    text_fields <- intersect(c("title","source_name_ch1","tissue","organism_ch1","disease","diagnosis"), names(raw))
    context <- tolower(paste(unlist(raw[text_fields]), collapse=" "))
    non_gastric <- grepl("breast|lung|colorectal|colon|prostate|ovarian|pancrea|liver|leukemia|lymphoma", context) &&
      !grepl("gastric|stomach", context)
    normal_only <- grepl("normal|adjacent|non.tumor|non-tumor|mucosa", context) &&
      !grepl("tumou?r|cancer|carcinoma|adenocarcinoma|gastric", context)
    platform <- if ("platform_id" %in% names(raw)) unique(na.omit(raw$platform_id))[1] else NA_character_
    ep <- detect_endpoint_pair(raw)
    reason <- NULL
    if (non_gastric) reason <- "Non-gastric disease detected in sample metadata"
    if (normal_only) reason <- "Normal/adjacent-only samples; no gastric tumor cohort detected"
    if (is.null(reason) && !ep$found) reason <- ep$reason
    af <- if (!is.na(platform)) find_platform_annotation(platform) else NA_character_
    missing_genes <- NA_character_
    if (is.null(reason) && is.na(af)) reason <- paste0("Platform annotation unavailable: ", platform)
    if (is.null(reason)) {
      pm <- tryCatch(read_probe_gene_map(af, platform, risk_genes), error=function(e) NULL)
      # Check actual expression rows, not annotation alone: some GEO matrices are
      # filtered subsets of their platform and omit otherwise annotated probes.
      px <- if (is.null(pm)) NULL else tryCatch(
        read_series_expression_subset(mf, unique(pm$probe_id)), error=function(e) NULL)
      present <- if (is.null(px) || is.null(pm)) character() else
        unique(pm$Gene[pm$probe_id %in% rownames(px)])
      missing <- setdiff(risk_genes, present)
      missing_genes <- if (length(missing)) paste(missing, collapse=", ") else ""
      # 部分较老平台可能确实未测量某个固定模型基因。为保留完整的六队列验证，
      # 不伪造原始表达量；缺失基因在后续队列内标准化矩阵中明确记为Z-score=0，
      # 等价于固定在该队列均值、风险分数贡献为0。缺失情况保留在审计表和逐队列QC表。
    }
    already <- ds %in% names(cohort_spec)
    status <- if (!is.null(reason)) "excluded" else if (already) "configured" else "auto_candidate"
    audit_note <- if (is.null(reason) && !is.na(missing_genes) && nzchar(missing_genes)) {
      paste0("Complete endpoint pair; missing genes assigned neutral within-cohort Z-score 0: ",
             missing_genes)
    } else if (is.null(reason)) {
      "Complete endpoint pair and model-gene annotation"
    } else reason
    rows[[ds]] <- tibble(Dataset=ds, Audit_status=status,
                         Endpoint=if (ep$found) ep$endpoint else NA_character_, Platform=platform,
                         Reason=audit_note,
                         Missing_model_genes=missing_genes, Matrix_file=mf)
    if (is.null(reason) && !already) auto_specs[[ds]] <- list(
      matrix_file=mf, annotation_file=af, platform=platform,
      clinical_source="auto_series_endpoint", clinical_file=NA_character_,
      endpoint=ep$endpoint, time_col=ep$time_col, event_col=ep$event_col,
      time_to_months=ep$time_to_months)
  }
  list(table=bind_rows(rows), specs=auto_specs)
}

# =======================================================================
# 4. 固定模型风险评分和生存统计
# =======================================================================
calculate_risk_score <- function(gene_expression, coefficients) {
  found <- intersect(names(coefficients), rownames(gene_expression))
  missing <- setdiff(names(coefficients), found)

  # 每个基因在当前GEO队列内做Z-score；不重新估计Cox系数。
  z_expression <- t(scale(t(gene_expression[found, , drop = FALSE])))
  z_expression[is.na(z_expression)] <- 0
  # 平台未测量的模型基因不能可靠插补原始表达，因此赋中性Z-score 0。
  # 这使该基因对每位患者的固定风险分数贡献均为0，并保持其他系数不变。
  if (length(missing)) {
    missing_matrix <- matrix(0, nrow=length(missing), ncol=ncol(gene_expression),
                             dimnames=list(missing, colnames(gene_expression)))
    z_expression <- rbind(z_expression, missing_matrix)
  }
  z_expression <- z_expression[names(coefficients), , drop=FALSE]
  score <- as.numeric(crossprod(coefficients, z_expression))
  names(score) <- colnames(z_expression)
  list(score = score, z_expression = z_expression, found_genes=found, missing_genes=missing)
}

validate_one_cohort <- function(dataset, spec) {
  message("\n===== Processing ", dataset, " =====")
  cohort_dir <- file.path(intermediate_dir, dataset)
  dir.create(cohort_dir, recursive = TRUE, showWarnings = FALSE)

  probe_map <- read_probe_gene_map(spec$annotation_file, spec$platform, risk_genes)
  write.csv(probe_map, file.path(cohort_dir, "01_model_gene_probe_map.csv"), row.names = FALSE)

  probe_expr <- read_series_expression_subset(spec$matrix_file, unique(probe_map$probe_id))
  write.csv(as.data.frame(probe_expr) %>% rownames_to_column("probe_id"),
            file.path(cohort_dir, "02_model_probe_expression.csv"), row.names = FALSE)

  gene_expr <- collapse_probe_to_gene(probe_expr, probe_map)
  selected_probes <- attr(gene_expr, "selected_probes")
  write.csv(selected_probes, file.path(cohort_dir, "03_selected_representative_probes.csv"),
            row.names = FALSE)
  write.csv(as.data.frame(gene_expr) %>% rownames_to_column("Gene"),
            file.path(cohort_dir, "04_model_gene_expression.csv"), row.names = FALSE)

  clinical_result <- make_survival_data(dataset, spec)
  write.csv(clinical_result$raw, file.path(cohort_dir, "05_clinical_raw_extracted.csv"),
            row.names = FALSE)
  write.csv(clinical_result$standardized,
            file.path(cohort_dir, "06_clinical_survival_standardized.csv"), row.names = FALSE)

  risk <- calculate_risk_score(gene_expr, risk_coef)
  write.csv(tibble(
    Gene=risk_genes,
    Expression_available=risk_genes %in% risk$found_genes,
    Missing_gene_handling=ifelse(risk_genes %in% risk$found_genes,
                                 "Observed expression; within-cohort Z-score",
                                 "Not measured on platform; neutral Z-score set to 0")
  ), file.path(cohort_dir, "07a_model_gene_availability_QC.csv"), row.names=FALSE)
  write.csv(as.data.frame(risk$z_expression) %>% rownames_to_column("Gene"),
            file.path(cohort_dir, "07_model_gene_zscore_expression.csv"), row.names = FALSE)

  score_df <- tibble(
    geo_accession = names(risk$score),
    risk_score = as.numeric(risk$score)
  )
  analysis <- inner_join(clinical_result$standardized, score_df, by = "geo_accession") %>%
    filter(!is.na(OS_time_months), OS_time_months > 0,
           !is.na(OS_event), OS_event %in% c(0, 1), !is.na(risk_score))

  # General QC requires >=20 subjects and >=5 events. GSE38749 is a newly added,
  # intrinsically small cohort (15 subjects, 10 deaths); it is retained with a
  # documented >=10-subject exception because Cox/KM estimation remains possible.
  minimum_n <- if (dataset == "GSE38749") 10L else 20L
  if (nrow(analysis) < minimum_n || sum(analysis$OS_event == 1) < 5) {
    stop(dataset, "有效生存样本或事件数不足。")
  }

  cutoff <- median(analysis$risk_score, na.rm = TRUE)
  analysis <- analysis %>% mutate(
    risk_group = factor(ifelse(risk_score >= cutoff, "High risk", "Low risk"),
                        levels = c("Low risk", "High risk"))
  )
  write.csv(analysis, file.path(cohort_dir, "08_final_KM_Cox_input.csv"), row.names = FALSE)

  # 连续风险评分Cox，反映每增加1个风险评分单位的风险变化。
  continuous_cox <- coxph(Surv(OS_time_months, OS_event) ~ risk_score, data = analysis)
  continuous_summary <- summary(continuous_cox)

  # 高低风险组Cox及log-rank检验。
  group_cox <- coxph(Surv(OS_time_months, OS_event) ~ risk_group, data = analysis)
  group_summary <- summary(group_cox)
  km_fit <- survfit(Surv(OS_time_months, OS_event) ~ risk_group, data = analysis)
  lr <- survdiff(Surv(OS_time_months, OS_event) ~ risk_group, data = analysis)
  logrank_p <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)

  stats <- tibble(
    Dataset = dataset,
    Endpoint = spec$endpoint,
    Platform = spec$platform,
    Model_gene_N = length(risk_genes),
    Available_model_gene_N = length(risk$found_genes),
    Missing_model_genes = if (length(risk$missing_genes))
      paste(risk$missing_genes, collapse=", ") else "None",
    Missing_gene_handling = if (length(risk$missing_genes))
      "Missing platform probes assigned neutral within-cohort Z-score 0" else "Not applicable",
    N = nrow(analysis),
    Events = sum(analysis$OS_event == 1),
    Low_risk_N = sum(analysis$risk_group == "Low risk"),
    High_risk_N = sum(analysis$risk_group == "High risk"),
    Median_cutoff = cutoff,
    Group_HR = group_summary$conf.int[1, "exp(coef)"],
    Group_lower95 = group_summary$conf.int[1, "lower .95"],
    Group_upper95 = group_summary$conf.int[1, "upper .95"],
    Group_Cox_P = group_summary$coefficients[1, "Pr(>|z|)"],
    Logrank_P = logrank_p,
    Continuous_HR = continuous_summary$conf.int[1, "exp(coef)"],
    Continuous_lower95 = continuous_summary$conf.int[1, "lower .95"],
    Continuous_upper95 = continuous_summary$conf.int[1, "upper .95"],
    Continuous_Cox_P = continuous_summary$coefficients[1, "Pr(>|z|)"]
  )
  write.csv(stats, file.path(cohort_dir, "09_survival_statistics.csv"), row.names = FALSE)

  annotation_text <- sprintf(
    "High vs Low: HR=%.2f (95%% CI %.2f-%.2f)\nCox P=%s; log-rank P=%s",
    stats$Group_HR, stats$Group_lower95, stats$Group_upper95,
    format.pval(stats$Group_Cox_P, digits = 2, eps = 0.001),
    format.pval(stats$Logrank_P, digits = 2, eps = 0.001)
  )
  if (length(risk$missing_genes)) {
    annotation_text <- paste0(annotation_text, "\nMissing probe/gene: ",
                              paste(risk$missing_genes, collapse=", "),
                              " (Z-score=0)")
  }

  p <- ggsurvplot(
    km_fit, data = analysis, pval = TRUE, risk.table = TRUE,
    risk.table.height = 0.28, conf.int = FALSE,
    palette = c("#4DBBD5", "#E64B35"),
    xlab = paste0(spec$endpoint, " time (months)"),
    ylab = paste0(spec$endpoint, " probability"),
    legend.title = "Risk group", legend.labs = c("Low risk", "High risk"),
    title = paste0(dataset, " ", spec$endpoint,
                   " validation of fixed TCGA-STAD risk model")
  )
  p$plot <- p$plot + theme_bw() +
    annotate("text", x = max(analysis$OS_time_months) * 0.03, y = 0.08,
             label = annotation_text, hjust = 0, size = 3)
  p$table <- p$table + theme_bw()

  pdf_file <- file.path(out_dir, paste0("Figure6_", dataset, "_risk_score_KM.pdf"))
  png_file <- file.path(out_dir, paste0("Figure6_", dataset, "_risk_score_KM.png"))

  # Windows 中若旧图正在被 PDF/图片查看器占用，直接覆盖会报错，并导致该队列
  # 已完成的统计结果无法进入汇总表。这里将绘图错误与统计计算解耦：优先写入
  # 标准文件名；文件被占用时改写带运行时间戳的备用文件，同时保留统计结果。
  save_km_plot <- function(path, type) {
    draw <- function(target) {
      opened <- FALSE
      tryCatch({
        if (type == "pdf") {
          pdf(target, width = 6.3, height = 7.0, useDingbats = FALSE)
        } else {
          png(target, width = 6.3, height = 7.0, units = "in",
              res = 300, bg = "white")
        }
        opened <- TRUE
        print(p)
        dev.off()
        opened <- FALSE
        TRUE
      }, error = function(e) {
        if (opened && dev.cur() > 1) dev.off()
        warning("无法写入图形文件 ", target, "：", conditionMessage(e))
        FALSE
      })
    }

    if (!draw(path)) {
      ext <- tools::file_ext(path)
      stem <- sub(paste0("\\.", ext, "$"), "", path)
      fallback <- paste0(stem, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".", ext)
      if (!draw(fallback)) warning(dataset, " 的 ", type, " 图未能保存。")
    }
  }
  save_km_plot(pdf_file, "pdf")
  save_km_plot(png_file, "png")

  stats
}

# =======================================================================
# 5. 批量运行独立GEO生存队列
# =======================================================================
# Discover every local GSE directory. Automatically eligible cohorts are added
# to the run; every other directory remains in the unified audit with its exact
# exclusion reason. Configured cohorts retain their supplementary-clinical parser.
auto_audit <- audit_local_geo()
# Reconcile automatic Series-Matrix inspection with cohorts that legitimately use
# a configured supplementary clinical table or a cohort-specific parser. This
# prevents a misleading "excluded + success" combination in the unified log.
configured_ids <- names(cohort_spec)
configured_endpoint <- setNames(vapply(cohort_spec, `[[`, character(1), "endpoint"),
                                configured_ids)
configured_source <- setNames(vapply(cohort_spec, `[[`, character(1), "clinical_source"),
                              configured_ids)
auto_audit$table <- auto_audit$table %>% mutate(
  Endpoint = ifelse(Dataset %in% configured_ids,
                    unname(configured_endpoint[Dataset]), Endpoint),
  Audit_status = ifelse(Dataset %in% configured_ids, "configured_for_run", Audit_status),
  Reason = ifelse(
    Dataset %in% configured_ids &
      grepl("No OS/DFS/RFS/DMFS/DSS/CSS", Reason, fixed=TRUE),
    paste0("Endpoint supplied by configured clinical parser: ",
           unname(configured_source[Dataset])),
    Reason
  )
)
cohort_spec <- c(cohort_spec, auto_audit$specs)
write.csv(auto_audit$table,
          file.path(out_dir, "Figure6_all_local_GEO_automatic_audit.csv"), row.names = FALSE)

result_list <- list()
processing_log <- list()

for (dataset in names(cohort_spec)) {
  result <- tryCatch({
    x <- validate_one_cohort(dataset, cohort_spec[[dataset]])
    processing_log[[dataset]] <- tibble(Dataset = dataset, Status = "success", Note = "ok")
    x
  }, error = function(e) {
    processing_log[[dataset]] <<- tibble(Dataset = dataset, Status = "failed",
                                         Note = conditionMessage(e))
    warning(dataset, "验证失败：", conditionMessage(e))
    NULL
  })
  if (!is.null(result)) result_list[[dataset]] <- result
  invisible(gc())
}

summary_table <- bind_rows(result_list)
if (nrow(summary_table)) {
  # Do not combine clinically different endpoints during multiplicity correction.
  summary_table <- summary_table %>% group_by(Endpoint) %>% mutate(
    FDR_Logrank = p.adjust(Logrank_P, method = "BH"),
    FDR_Group_Cox = p.adjust(Group_Cox_P, method = "BH"),
    FDR_Continuous_Cox = p.adjust(Continuous_Cox_P, method = "BH")
  ) %>% ungroup()
}
write.csv(summary_table, file.path(out_dir, "Figure6_GEO_validation_summary.csv"),
          row.names = FALSE)
write.csv(bind_rows(processing_log), file.path(out_dir, "Figure6_processing_log.csv"),
          row.names = FALSE)

# A complete processing log contains all local directories, including cohorts
# rejected before KM/Cox. For attempted cohorts, runtime status supersedes the
# pre-run audit status and supplies expression-level missing-gene errors.
all_processing_log <- auto_audit$table %>%
  transmute(Dataset, Endpoint, Audit_status, Status = "not_run", Note = Reason,
            Platform, Missing_model_genes, Matrix_file) %>%
  left_join(bind_rows(processing_log) %>% rename(Runtime_status=Status, Runtime_note=Note),
            by="Dataset") %>%
  mutate(Status=ifelse(!is.na(Runtime_status), Runtime_status, Status),
         Note=ifelse(!is.na(Runtime_note), Runtime_note, Note)) %>%
  select(-Runtime_status, -Runtime_note) %>% arrange(Dataset)
write.csv(all_processing_log,
          file.path(out_dir, "Figure6_processing_log_all_local_GEO.csv"), row.names = FALSE)

# Detailed cohort-use audit requested for the 2026-08-12 rerun. This table joins
# the declared inputs, actual processing status, final analysis N/events and the
# explicit exclusion reasons, so every local validation decision is traceable.
included_usage <- tibble(
  Dataset = names(cohort_spec),
  Endpoint = vapply(cohort_spec, `[[`, character(1), "endpoint"),
  Decision = paste0("Evaluated for ", vapply(cohort_spec, `[[`, character(1), "endpoint"),
                    " validation"),
  Platform = vapply(cohort_spec, `[[`, character(1), "platform"),
  Clinical_source = vapply(cohort_spec, `[[`, character(1), "clinical_source"),
  Matrix_file = vapply(cohort_spec, `[[`, character(1), "matrix_file"),
  Reason = "Complete endpoint-specific time/event pair configured"
) %>%
  left_join(bind_rows(processing_log), by = "Dataset") %>%
  left_join(summary_table %>% select(Dataset, N, Events, Low_risk_N, High_risk_N),
            by = "Dataset")
excluded_usage <- excluded_cohorts %>% mutate(
  Endpoint = NA_character_, Decision = "Excluded from survival validation", Platform = NA_character_,
  Clinical_source = NA_character_, Matrix_file = NA_character_,
  Status = "not run", Note = Reason,
  N = NA_integer_, Events = NA_integer_, Low_risk_N = NA_integer_, High_risk_N = NA_integer_
)
cohort_usage <- bind_rows(included_usage, excluded_usage) %>%
  select(Dataset, Endpoint, Decision, Platform, Clinical_source, Status, Note,
         N, Events, Low_risk_N, High_risk_N, Reason, Matrix_file) %>%
  arrange(Dataset)
write.csv(cohort_usage, file.path(out_dir, "Figure6_GEO_cohort_detailed_usage.csv"),
          row.names = FALSE)

# Human-readable Chinese report. Keep the CSV above for machine-readable auditing.
sink(file.path(out_dir, "Figure6_GEO_cohort_detailed_usage_CN.txt"))
cat("Figure 6 GEO 队列详细使用情况\n")
cat("生成时间：", as.character(Sys.time()), "\n\n", sep = "")
cat("分析原则：固定使用 TCGA-STAD LASSO-Cox 系数；各 GEO 队列内做基因 Z-score；",
    "按队列内风险评分中位数分组；不在 GEO 队列中重新拟合系数。\n\n", sep = "")
cat("纳入生存验证的队列及实际使用量：\n")
print(included_usage %>% select(Dataset, Platform, Clinical_source, Status,
                                Note, N, Events, Low_risk_N, High_risk_N),
      row.names = FALSE)
cat("\n排除队列及原因：\n")
print(excluded_cohorts, row.names = FALSE)
cat("\n特殊说明：GSE38749 为新增小样本队列（原始 15 例），状态字段 alive/death；",
    "因死亡事件数充足，使用最低 10 例且至少 5 个事件的门槛，并应谨慎解释宽置信区间。\n", sep = "")
cat("GSE84437 如处理失败，原因通常为 GPL6947 当前注释无法覆盖固定模型全部基因；",
    "不会用缺失基因的替代系数强行计算。\n", sep = "")
sink()

# =======================================================================
# 6. 运行说明及代码归档
# =======================================================================
sink(file.path(out_dir, "Figure6_GEO_validation_run_notes.txt"))
cat("Gastric Figure6 GEO independent validation\n")
cat("Run time:", as.character(Sys.time()), "\n")
cat("R executable: ", file.path(R.home("bin"), "Rscript"), "\n", sep = "")
cat("Output directory:", out_dir, "\n")
cat("Source script preserved unchanged: Figure6_Top2_0812/Figure6_GEO_independent_validation_code_0812.R\n")
cat("Reference training script: Figure5/Figure5_gastric.R\n")
cat("Model coefficient source: manual_model block in this script\n")
cat("TCGA complete model parameter file:", complete_model_file, "\n")
cat("Fixed model coefficients:\n")
print(coef_df)
cat("TCGA training means and standard deviations (recorded for reproducibility; ",
    "not directly applied across microarray platforms):\n", sep = "")
print(training_model_parameters)
cat("\nValidation rule: within-cohort gene Z-score; fixed TCGA coefficients; ",
    "within-cohort median risk cutoff.\n", sep = "")
cat("Cohorts with survival validation:", paste(names(cohort_spec), collapse = ", "), "\n")
cat("Other downloaded cohorts were excluded from survival analysis because the local files ",
    "did not contain both survival time and event.\n\n", sep = "")
cat("Validation summary:\n")
print(summary_table)
cat("\nSession information:\n")
print(sessionInfo())
sink()

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
if (length(script_arg)) {
  current_script <- normalizePath(sub("^--file=", "", script_arg[1]),
                                  winslash = "/", mustWork = TRUE)
  archived_script <- file.path(out_dir, "Figure6_analysis.R")
  if (normalizePath(current_script, winslash = "/", mustWork = FALSE) !=
      normalizePath(archived_script, winslash = "/", mustWork = FALSE)) {
    file.copy(current_script, archived_script, overwrite = TRUE)
  }
}

cat("GEO独立验证完成。\n")
cat("输出目录：", out_dir, "\n", sep = "")
