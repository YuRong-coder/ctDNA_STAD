############################################################
## Figure S5: LASSO-Cox model selection and coefficients
##
## Panel A: LASSO cross-validation curve.
## Panel B: Coefficients of the seven selected genes.
## Panel C: Univariable Cox forest plot for clinical predictors and risk group.
##
## Required files in ./data:
##   - TCGA-STAD.star_tpm.tsv
##   - TCGA-STAD.survival.gz
##   - TCGA-STAD.clinical.gz
##
## Run with: Rscript FigureS5_analysis.R
## Outputs are written under ./output/FigureS5.
############################################################

required_packages <- c(
  "data.table", "dplyr", "tibble", "stringr", "survival", "glmnet",
  "ggplot2", "org.Hs.eg.db", "AnnotationDbi"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) stop("缺少R包：", paste(missing_packages, collapse = ", "))

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(survival)
  library(glmnet)
  library(ggplot2)
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

set.seed(20260717)
script_dir <- get_script_dir()
input_dir <- file.path(script_dir, "data")
out_dir <- file.path(script_dir, "output", "FigureS5")
output_dir <- out_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expr_file <- file.path(input_dir, "TCGA-STAD.star_tpm.tsv")
surv_file <- file.path(input_dir, "TCGA-STAD.survival.gz")
clin_file <- file.path(input_dir, "TCGA-STAD.clinical.gz")
for (f in c(expr_file, surv_file, clin_file)) {
  if (!file.exists(f)) stop("Input file not found: ", f)
}

############################################################
## Panel A: LASSO cross-validation
############################################################
candidate_genes <- c(    "ARID1B","ATM","BRCA2","C1S",
                         "CCL16","CD72","CDH11","CDH6",
                         "CLDN11","COL12A1","COL15A1","COL23A1",
                         "EP300","EPHB2","EPHB6","ERBB4",
                         "HSPG2","IL10RB","IL16","ITGA4",
                         "ITGA6","ITGAL","KDR","LGALS3BP",
                         "LRP1","MMP24","MMP25","MMP28",
                         "MUC16","MUC17","MUC21","MUC4",
                         "MUC5B","MUC6","PIK3C2A","PLXNC1",
                         "POLE","TET1","TET3","TNFRSF21")




write.csv(data.frame(Gene=candidate_genes),file.path(out_dir,"00_candidate_14_genes.csv"),row.names=FALSE)

# 2. 表达矩阵：仅保留样本类型01的原发肿瘤 ------------------------------
message("Reading expression matrix ...")
expr_raw <- fread(expr_file,data.table=FALSE,check.names=FALSE)
colnames(expr_raw)[1] <- "Ensembl_ID_version"
expr_raw$Ensembl_ID <- sub("\\..*$","",expr_raw$Ensembl_ID_version)
tcga_cols <- grep("^TCGA-",colnames(expr_raw),value=TRUE)
tumor_cols <- tcga_cols[substr(tcga_cols,14,15)=="01"]
if(!length(tumor_cols)) stop("未识别到TCGA原发肿瘤样本（01）。")

# 仅转换14个候选基因的Symbol，减少注释歧义和运行负担。
ens <- AnnotationDbi::mapIds(org.Hs.eg.db,keys=candidate_genes,column="ENSEMBL",
                             keytype="SYMBOL",multiVals="first")
gene_map <- data.frame(Gene=names(ens),Ensembl_ID=unname(ens))
write.csv(gene_map,file.path(out_dir,"01_candidate_gene_Ensembl_map.csv"),row.names=FALSE)
expr_gene <- expr_raw %>% inner_join(gene_map,by="Ensembl_ID") %>% dplyr::select(Gene,all_of(tumor_cols))
matched_genes <- sort(unique(expr_gene$Gene)); missing_genes <- setdiff(candidate_genes,matched_genes)
write.csv(data.frame(Gene=matched_genes),file.path(out_dir,"02_expression_matched_genes.csv"),row.names=FALSE)
write.csv(data.frame(Gene=missing_genes),file.path(out_dir,"02_expression_missing_genes.csv"),row.names=FALSE)
if(length(matched_genes)<2) stop("表达矩阵中匹配的候选基因少于2个。")
expr_gene[,tumor_cols] <- lapply(expr_gene[,tumor_cols,drop=FALSE],as.numeric)
expr_gene <- expr_gene %>% group_by(Gene) %>%
  summarise(across(all_of(tumor_cols),~mean(.x,na.rm=TRUE)),.groups="drop")
expr_t <- as.data.frame(t(column_to_rownames(expr_gene,"Gene"))) %>% rownames_to_column("sample_id") %>%
  mutate(patient_id=substr(sample_id,1,12)) %>% dplyr::select(patient_id,all_of(matched_genes)) %>%
  group_by(patient_id) %>% summarise(across(all_of(matched_genes),~mean(as.numeric(.x),na.rm=TRUE)),.groups="drop")
# Xena star_fpkm矩阵已经是对数尺度，因此不再重复log2转换。
write.csv(expr_t,file.path(out_dir,"03_candidate_expression_primary_tumor.csv"),row.names=FALSE)
rm(expr_raw,expr_gene); invisible(gc())

# 3. 合并OS、临床资料和表达数据 --------------------------------------
surv_raw <- fread(surv_file,data.table=FALSE)
surv_df <- surv_raw %>% transmute(patient_id=substr(as.character(sample),1,12),
  OS_time=as.numeric(OS.time),OS_event=as.numeric(OS)) %>%
  filter(!is.na(OS_time),OS_time>0,OS_event %in% c(0,1)) %>% distinct(patient_id,.keep_all=TRUE)
clean_na <- function(x){x<-trimws(as.character(x));x[tolower(x)%in%c("","na","not reported","not available","unknown","--")]<-NA;x}
clin_raw <- fread(clin_file,data.table=FALSE,check.names=FALSE)
clin_df <- clin_raw %>% mutate(patient_id=substr(as.character(sample),1,12),
  age=suppressWarnings(as.numeric(clean_na(age_at_index.demographic))),
  age_group=case_when(age<60~"<60",age>=60~">=60",TRUE~NA_character_),
  gender=str_to_title(clean_na(gender.demographic)),
  stage_raw=str_to_upper(clean_na(ajcc_pathologic_stage.diagnoses)),
  stage=case_when(str_detect(stage_raw,"STAGE IV")~"Stage IV",str_detect(stage_raw,"STAGE III")~"Stage III",
                  str_detect(stage_raw,"STAGE II")~"Stage II",str_detect(stage_raw,"STAGE I")~"Stage I",TRUE~NA_character_)) %>%
  dplyr::select(patient_id,age,age_group,gender,stage) %>% distinct(patient_id,.keep_all=TRUE)
analysis_df <- inner_join(surv_df,expr_t,by="patient_id") %>% left_join(clin_df,by="patient_id")
genes_use <- intersect(candidate_genes,colnames(analysis_df))
write.csv(analysis_df,file.path(out_dir,"04_expression_survival_clinical_merged.csv"),row.names=FALSE)

# 计算单因素 Cox 结果，仅用于确定进入 LASSO 的候选基因；不绘制 Figure 5A。
univ_res <- bind_rows(lapply(genes_use,function(g){
  z<-analysis_df[,c("OS_time","OS_event",g)];colnames(z)[3]<-"expression";z<-z[complete.cases(z),]
  fit<-tryCatch(coxph(Surv(OS_time,OS_event)~expression,data=z),error=function(e)NULL)
  if(is.null(fit)) return(data.frame(Gene=g,HR=NA,lower95=NA,upper95=NA,pvalue=NA,n=nrow(z),events=sum(z$OS_event)))
  s<-summary(fit);data.frame(Gene=g,HR=s$conf.int[1,"exp(coef)"],lower95=s$conf.int[1,"lower .95"],
    upper95=s$conf.int[1,"upper .95"],pvalue=s$coefficients[1,"Pr(>|z|)"],n=nrow(z),events=sum(z$OS_event))
})) %>% mutate(FDR=p.adjust(pvalue,"BH"),Direction=case_when(pvalue<.1&HR>1~"Risk",pvalue<.1&HR<1~"Protective",TRUE~"NS")) %>% arrange(HR)
write.csv(univ_res,file.path(out_dir,"05_univariate_Cox_14_genes.csv"),row.names=FALSE)
# 5. LASSO-Cox：先取单因素P<0.1；不足2个时用全部候选基因 ----------
lasso_genes <- univ_res%>%filter(!is.na(pvalue),pvalue<.1)%>%pull(Gene)
lasso_rule <- "Univariate Cox P < 0.1"
if(length(lasso_genes)<2){lasso_genes<-genes_use;lasso_rule<-"P<0.1 genes <2; all candidates used"}
# 对LASSO候选基因使用TCGA训练队列的均值和标准差进行Z-score标准化。
# scale()返回的scaled:center和scaled:scale必须保存，才能完整复现训练模型。
x<-scale(as.matrix(analysis_df[,lasso_genes,drop=FALSE]));y<-Surv(analysis_df$OS_time,analysis_df$OS_event)

# 保存TCGA训练集的标准化参数。Mean为训练集均值，SD为训练集标准差。
# 注意：这些参数适用于与TCGA相同表达处理流程的数据；不同芯片平台不能直接套用。
training_scaling_df<-data.frame(
  Gene=colnames(x),
  Training_mean=as.numeric(attr(x,"scaled:center")),
  Training_SD=as.numeric(attr(x,"scaled:scale")),
  stringsAsFactors=FALSE
)
write.csv(training_scaling_df,
          file.path(out_dir,"06_LASSO_training_expression_scaling.csv"),
          row.names=FALSE)
cvfit<-cv.glmnet(x,y,family="cox",alpha=1,nfolds=10,type.measure="deviance",standardize=FALSE)
cm<-as.matrix(coef(cvfit,s="lambda.min"));risk_coef<-cm[cm[,1]!=0,1]
if(!length(risk_coef)){cm<-as.matrix(coef(cvfit,s="lambda.1se"));risk_coef<-cm[cm[,1]!=0,1]}
if(!length(risk_coef)) stop("LASSO未选择任何基因。")
coef_df<-data.frame(Gene=names(risk_coef),Coefficient=as.numeric(risk_coef))
write.csv(coef_df,file.path(out_dir,"06_LASSO_coefficients.csv"),row.names=FALSE)

# 将最终非零系数与对应的TCGA训练均值/标准差合并保存为完整模型参数表。
# 该文件便于以后对同一表达体系的独立样本直接复现风险评分。
final_model_parameters<-coef_df%>%
  left_join(training_scaling_df,by="Gene")%>%
  mutate(Score_formula=paste0("(",Gene," - ",signif(Training_mean,8),") / ",
                              signif(Training_SD,8)," * ",signif(Coefficient,8)))
write.csv(final_model_parameters,
          file.path(out_dir,"06_LASSO_complete_model_parameters.csv"),
          row.names=FALSE)
pdf(file.path(out_dir,"FigureS5A_LASSO_cross_validation.pdf"),6,5);plot(cvfit);dev.off()

# Figure S5C 使用与主 Figure 5 相同的 LASSO 风险评分数据。
analysis_df$risk_score<-as.vector(x[,names(risk_coef),drop=FALSE]%*%risk_coef)
cutoff<-median(analysis_df$risk_score,na.rm=TRUE)
analysis_df$risk_group<-factor(ifelse(analysis_df$risk_score>=cutoff,"High risk","Low risk"),levels=c("Low risk","High risk"))
write.csv(analysis_df,file.path(out_dir,"07_LASSO_risk_score_data.csv"),row.names=FALSE)

############################################################
## Panel B: Selected-gene coefficients
############################################################
genes <- c("COL15A1", "PLXNC1", "CDH11", "CCL16", "CLDN11", "ERBB4", "CDH6")
coefficients <- c(
  0.00337975891621139,
  0.0254558844392741,
  0.0365614804586186,
  0.0148733152997758,
  0.133919085277019,
  0.0782012384327468,
  0.163505577557928
)

plot_lasso_coefficients <- function() {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  par(
    mar = c(6.2, 7.2, 3.9, 1.2),
    mgp = c(2.7, 0.75, 0),
    tcl = -0.25,
    family = "sans",
    xaxs = "i"
  )

  # Reverse the vectors so the first gene is displayed at the top.
  plot_genes <- rev(genes)
  plot_coef <- rev(coefficients)
  bar_colors <- ifelse(plot_coef >= 0, "#D65243", "#4DABC0")
  x_limit <- c(min(-0.01, min(plot_coef) * 1.08), max(plot_coef) * 1.16)

  barplot(
    height = plot_coef,
    names.arg = plot_genes,
    horiz = TRUE,
    col = bar_colors,
    border = "black",
    lwd = 0.9,
    space = 0.42,
    xlim = x_limit,
    axes = FALSE,
    cex.names = 1.05,
    las = 1
  )

  axis(
    side = 1,
    at = seq(0, 0.16, by = 0.04),
    labels = sprintf("%.2f", seq(0, 0.16, by = 0.04)),
    cex.axis = 0.95,
    lwd = 1.1,
    lwd.ticks = 1.1
  )
  abline(v = 0, col = "black", lwd = 1.1)
  box(bty = "l", lwd = 1.1)

  title(
    main = "LASSO Cox selected genes",
    cex.main = 1.25,
    font.main = 2,
    line = 1.0
  )
  mtext("LASSO coefficient", side = 1, line = 3.35, cex = 1.05)
  mtext("B", side = 3, adj = -0.19, line = 1.15, cex = 1.45, font = 2)
}



pdf(
  file = file.path(output_dir, "FigureS5B_LASSO_Cox_selected_genes.pdf"),
  width = 8.0,
  height = 5.6,
  useDingbats = FALSE
)
plot_lasso_coefficients()
dev.off()

png(
  filename = file.path(output_dir, "FigureS5B_LASSO_Cox_selected_genes.png"),
  width = 2400,
  height = 1680,
  res = 300,
  type = "windows"
)
plot_lasso_coefficients()
dev.off()

write.csv(
  data.frame(Gene = genes, Coefficient = coefficients),
  file = file.path(output_dir, "FigureS5B_LASSO_Cox_coefficients.csv"),
  row.names = FALSE,
  quote = TRUE
)

message("Figure files were written to: ", output_dir)

############################################################
## Panel C: Univariable Cox forest plot
############################################################
input_file <- file.path(out_dir, "07_LASSO_risk_score_data.csv")
output_dir <- out_dir
raw <- read.csv(input_file, check.names = FALSE)
required <- c("patient_id", "OS_time", "OS_event", "risk_group",
              "age_group", "gender", "stage")
missing_columns <- setdiff(required, names(raw))
if (length(missing_columns) > 0) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

clean_text <- function(x) {
  x <- trimws(as.character(x))
  x[x == "" | tolower(x) %in% c("na", "n/a", "unknown", "not reported")] <- NA
  x
}

dat <- raw %>%
  transmute(
    patient_id = as.character(patient_id),
    OS_time = suppressWarnings(as.numeric(OS_time)),
    OS_event = suppressWarnings(as.numeric(OS_event)),
    risk_group = clean_text(risk_group),
    age_group = clean_text(age_group),
    gender = clean_text(gender),
    stage = clean_text(stage)
  ) %>%
  filter(!is.na(OS_time), OS_time > 0, OS_event %in% c(0, 1)) %>%
  mutate(
    risk_group = case_when(
      grepl("high", risk_group, ignore.case = TRUE) ~ "High risk",
      grepl("low", risk_group, ignore.case = TRUE) ~ "Low risk",
      TRUE ~ NA_character_
    ),
    age_group = case_when(
      age_group %in% c(">=60", "≥60") ~ ">=60",
      age_group == "<60" ~ "<60",
      TRUE ~ NA_character_
    ),
    gender = case_when(
      grepl("^male$", gender, ignore.case = TRUE) ~ "Male",
      grepl("^female$", gender, ignore.case = TRUE) ~ "Female",
      TRUE ~ NA_character_
    ),
    stage = case_when(
      grepl("stage[[:space:]]*iv", stage, ignore.case = TRUE) ~ "Stage IV",
      grepl("stage[[:space:]]*iii", stage, ignore.case = TRUE) ~ "Stage III",
      grepl("stage[[:space:]]*ii", stage, ignore.case = TRUE) ~ "Stage II",
      grepl("stage[[:space:]]*i", stage, ignore.case = TRUE) ~ "Stage I",
      TRUE ~ NA_character_
    ),
    risk_group = factor(risk_group, levels = c("Low risk", "High risk")),
    age_group = factor(age_group, levels = c("<60", ">=60")),
    gender = factor(gender, levels = c("Female", "Male")),
    stage = factor(stage, levels = c("Stage I", "Stage II", "Stage III", "Stage IV"))
  )

fit_one <- function(variable) {
  z <- dat[, c("OS_time", "OS_event", variable)]
  z <- z[complete.cases(z), , drop = FALSE]
  if (nrow(z) == 0 || length(unique(z[[variable]])) < 2) {
    stop("Insufficient usable data for variable: ", variable)
  }
  fit <- coxph(as.formula(paste0("Surv(OS_time, OS_event) ~ ", variable)), data = z)
  s <- summary(fit)
  data.frame(
    variable = variable,
    term = rownames(s$coefficients),
    HR = s$conf.int[, "exp(coef)"],
    lower95 = s$conf.int[, "lower .95"],
    upper95 = s$conf.int[, "upper .95"],
    pvalue = s$coefficients[, "Pr(>|z|)"],
    n = s$n,
    events = s$nevent,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

res <- bind_rows(lapply(c("risk_group", "age_group", "gender", "stage"), fit_one)) %>%
  mutate(
    label = case_when(
      variable == "risk_group" ~ "Risk group (High vs Low)",
      variable == "age_group" ~ "Age group (>=60 vs <60)",
      variable == "gender" ~ "Gender (Male vs Female)",
      term == "stageStage II" ~ "Stage II vs Stage I",
      term == "stageStage III" ~ "Stage III vs Stage I",
      term == "stageStage IV" ~ "Stage IV vs Stage I",
      TRUE ~ term
    ),
    significance = case_when(
      pvalue < 0.05 & HR > 1 ~ "Risk",
      pvalue < 0.05 & HR < 1 ~ "Protective",
      TRUE ~ "NS"
    ),
    estimate_label = sprintf("%.2f (%.2f-%.2f)", HR, lower95, upper95),
    p_label = ifelse(pvalue < 0.001, "<0.001", sprintf("%.3f", pvalue))
  )

label_order <- c(
  "Risk group (High vs Low)", "Age group (>=60 vs <60)",
  "Gender (Male vs Female)", "Stage II vs Stage I",
  "Stage III vs Stage I", "Stage IV vs Stage I"
)
res$label <- factor(res$label, levels = rev(label_order))

write.csv(dat, file.path(output_dir, "FigureS5C_TCGA_STAD_analysis_data.csv"), row.names = FALSE)
write.csv(res, file.path(output_dir, "FigureS5C_TCGA_STAD_univariate_Cox_results.csv"), row.names = FALSE)

x_min <- min(c(res$lower95, 1), na.rm = TRUE)
x_max <- max(c(res$upper95, 1), na.rm = TRUE)
label_x <- exp(log(res$upper95) + 0.035 * (log(x_max) - log(x_min)))
plot_x_max <- max(label_x * 1.65, na.rm = TRUE)

p <- ggplot(res, aes(x = HR, y = label, color = significance)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.55, color = "#8C8C8C") +
  geom_errorbar(aes(xmin = lower95, xmax = upper95), orientation = "y",
                width = 0.18, linewidth = 0.8) +
  geom_point(size = 3.2) +
  geom_text(aes(x = label_x, label = estimate_label), hjust = 0,
            color = "#303030", size = 4.0, show.legend = FALSE) +
  scale_x_log10(limits = c(x_min * 0.85, plot_x_max)) +
  scale_color_manual(values = c(Risk = "#E64B35", Protective = "#4DBBD5", NS = "#7A7A7A"),
                     breaks = c("Risk", "Protective", "NS"), drop = TRUE) +
  labs(title = "Univariate Cox regression", x = "Hazard ratio (log scale)", y = NULL, color = NULL) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "plain"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 12),
    legend.position = c(0.12, 0.18),
    legend.background = element_blank(),
    legend.key = element_blank(),
    plot.margin = margin(12, 24, 8, 8)
  )

ggsave(file.path(output_dir, "FigureS5C_TCGA_STAD_univariate_Cox_forest.pdf"),
       p, width = 10.2, height = 5.8, device = cairo_pdf)
ggsave(file.path(output_dir, "FigureS5C_TCGA_STAD_univariate_Cox_forest.png"),
       p, width = 10.2, height = 5.8, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "FigureS5C_TCGA_STAD_univariate_Cox_forest.tiff"),
       p, width = 10.2, height = 5.8, dpi = 300, compression = "lzw", bg = "white")

notes <- c(
  "TCGA-STAD Figure S5C: univariate Cox regression forest plot",
  paste("Generated:", Sys.time()),
  paste("Input:", input_file),
  paste("Eligible survival samples:", nrow(dat)),
  paste("OS events:", sum(dat$OS_event)),
  "Separate univariate Cox models were fitted for risk group, age group, gender, and stage.",
  "Reference levels: Low risk, age <60, Female, and Stage I.",
  "Red: P < 0.05 and HR > 1; blue: P < 0.05 and HR < 1; gray: P >= 0.05.",
  "Risk group was inherited from the upstream TCGA-STAD LASSO-Cox analysis and median cutoff.",
  "",
  capture.output(sessionInfo())
)
writeLines(notes, file.path(output_dir, "FigureS5C_run_notes.txt"))

message("Completed: ", output_dir)
