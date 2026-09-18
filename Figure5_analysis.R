############################################################
## Figure 5 unified TCGA-STAD analysis pipeline (panels A-D)
##
## Panel A: univariable Cox forest plot.
## Panel B: multivariable Cox forest plot (former Figure S6).
## Panel C: LASSO risk-group Kaplan-Meier curve.
## Panel D: paired tumor-versus-normal expression of final risk genes.
##
## Put all input files in ./data:
##   - TCGA-STAD.star_tpm.tsv
##   - TCGA-STAD.survival.gz
##   - TCGA-STAD.clinical.gz
##
## Run with: Rscript Figure5_analysis.R
## Outputs are written under ./output/Figure5.
############################################################

required_packages <- c(
  "data.table", "dplyr", "tibble", "stringr", "survival", "survminer",
  "glmnet", "ggplot2", "org.Hs.eg.db", "AnnotationDbi", "tidyr"
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
  library(tibble)
  library(stringr)
  library(survival)
  library(survminer)
  library(glmnet)
  library(ggplot2)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(tidyr)
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
figure5_output_dir <- file.path(script_dir, "output", "Figure5")
dir.create(figure5_output_dir, recursive = TRUE, showWarnings = FALSE)

required_input_files <- c(
  "TCGA-STAD.star_tpm.tsv",
  "TCGA-STAD.survival.gz",
  "TCGA-STAD.clinical.gz"
)
missing_input_files <- required_input_files[
  !file.exists(file.path(input_dir, required_input_files))
]
if (length(missing_input_files) > 0L) {
  stop("data/ 中缺少输入文件：", paste(missing_input_files, collapse = ", "))
}

############################################################
## Panels A and C: Prognostic model and LASSO risk KM
############################################################
run_figure5ac <- function(input_dir, output_root) {
  out_dir <- file.path(output_root, "AC_prognostic_model")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  expr_file <- file.path(input_dir, "TCGA-STAD.star_tpm.tsv")
  surv_file <- file.path(input_dir, "TCGA-STAD.survival.gz")
  clin_file <- file.path(input_dir, "TCGA-STAD.clinical.gz")
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

# 4. Figure5A：14基因单因素Cox森林图 ---------------------------------
univ_res <- bind_rows(lapply(genes_use,function(g){
  z<-analysis_df[,c("OS_time","OS_event",g)];colnames(z)[3]<-"expression";z<-z[complete.cases(z),]
  fit<-tryCatch(coxph(Surv(OS_time,OS_event)~expression,data=z),error=function(e)NULL)
  if(is.null(fit)) return(data.frame(Gene=g,HR=NA,lower95=NA,upper95=NA,pvalue=NA,n=nrow(z),events=sum(z$OS_event)))
  s<-summary(fit);data.frame(Gene=g,HR=s$conf.int[1,"exp(coef)"],lower95=s$conf.int[1,"lower .95"],
    upper95=s$conf.int[1,"upper .95"],pvalue=s$coefficients[1,"Pr(>|z|)"],n=nrow(z),events=sum(z$OS_event))
})) %>% mutate(FDR=p.adjust(pvalue,"BH"),Direction=case_when(pvalue<.1&HR>1~"Risk",pvalue<.1&HR<1~"Protective",TRUE~"NS")) %>% arrange(HR)
write.csv(univ_res,file.path(out_dir,"05_univariate_Cox_14_genes.csv"),row.names=FALSE)
pAdf<-univ_res%>%filter(is.finite(HR),is.finite(lower95),is.finite(upper95))%>%mutate(Gene=factor(Gene,levels=Gene))
pA<-ggplot(pAdf,aes(HR,Gene,color=Direction))+geom_vline(xintercept=1,linetype=2,color="grey45")+
  geom_errorbar(aes(xmin=lower95,xmax=upper95),orientation="y",width=.2)+geom_point(size=2.7)+scale_x_log10()+
  scale_color_manual(values=c(Risk="#E64B35",Protective="#4DBBD5",NS="#7A7A7A"))+
  labs(x="Hazard ratio per expression unit (log scale)",y=NULL,color=NULL,title="Univariate Cox regression",
       subtitle="TCGA-STAD primary tumors")+theme_classic(base_size=11)+
  theme(plot.title=element_text(hjust=.5,face="bold"),plot.subtitle=element_text(hjust=.5))
ggsave(file.path(out_dir,"Figure5A_univariate_Cox.pdf"),pA,width=6.5,height=5.5)
ggsave(file.path(out_dir,"Figure5A_univariate_Cox.png"),pA,width=6.5,height=5.5,dpi=300)

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
analysis_df$risk_score<-as.vector(x[,names(risk_coef),drop=FALSE]%*%risk_coef)
cutoff<-median(analysis_df$risk_score,na.rm=TRUE)
analysis_df$risk_group<-factor(ifelse(analysis_df$risk_score>=cutoff,"High risk","Low risk"),levels=c("Low risk","High risk"))
write.csv(analysis_df,file.path(out_dir,"07_LASSO_risk_score_data.csv"),row.names=FALSE)

# 6. Figure5C：LASSO风险组KM曲线 --------------------------------------
fitC<-survfit(Surv(OS_time,OS_event)~risk_group,data=analysis_df)
coxC<-summary(coxph(Surv(OS_time,OS_event)~risk_group,data=analysis_df))
hrC<-sprintf("High vs Low: HR=%.2f (95%% CI %.2f-%.2f), P=%s",coxC$conf.int[1,"exp(coef)"],
  coxC$conf.int[1,"lower .95"],coxC$conf.int[1,"upper .95"],
  ifelse(coxC$coefficients[1,"Pr(>|z|)"]<.001,
         formatC(coxC$coefficients[1,"Pr(>|z|)"],format="e",digits=2),
         sprintf("%.3f",coxC$coefficients[1,"Pr(>|z|)"])))
pC<-ggsurvplot(fitC,data=analysis_df,pval=TRUE,risk.table=TRUE,palette=c("#4DBBD5","#E64B35"),
  xlab="Time (days)",ylab="Overall survival probability",legend.title="",legend.labs=c("Low risk","High risk"),title="LASSO risk score")
pC$plot<-pC$plot+theme_bw()+annotate("text",x=max(analysis_df$OS_time)*.03,y=.08,label=hrC,hjust=0,size=3)
pdf(file.path(out_dir,"Figure5C_LASSO_risk_KM.pdf"),6,6.5);print(pC);dev.off()
}

############################################################
## Panel B: Multivariable Cox forest plot
############################################################
run_figure5b <- function(output_root) {
  input_file <- file.path(output_root, "AC_prognostic_model", "07_LASSO_risk_score_data.csv")
  output_dir <- file.path(output_root, "B_multivariable_Cox")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
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
  ) %>%
  filter(complete.cases(OS_time, OS_event, risk_group, age_group, gender, stage))

if (nrow(dat) == 0) stop("No complete cases available for multivariable Cox regression")

fit <- coxph(
  Surv(OS_time, OS_event) ~ risk_group + age_group + gender + stage,
  data = dat,
  ties = "efron"
)
s <- summary(fit)

res <- data.frame(
  term = rownames(s$coefficients),
  HR = s$conf.int[, "exp(coef)"],
  lower95 = s$conf.int[, "lower .95"],
  upper95 = s$conf.int[, "upper .95"],
  pvalue = s$coefficients[, "Pr(>|z|)"],
  n = s$n,
  events = s$nevent,
  stringsAsFactors = FALSE,
  row.names = NULL
) %>%
  mutate(
    label = case_when(
      term == "risk_groupHigh risk" ~ "Risk group (High vs Low)",
      term == "age_group>=60" ~ "Age group (>=60 vs <60)",
      term == "genderMale" ~ "Gender (Male vs Female)",
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

write.csv(dat, file.path(output_dir, "Figure5B_TCGA_STAD_multivariate_analysis_data.csv"),
          row.names = FALSE)
write.csv(res, file.path(output_dir, "Figure5B_TCGA_STAD_multivariate_Cox_results.csv"),
          row.names = FALSE)

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
  labs(title = "Multivariable Cox regression", x = "Hazard ratio (log scale)",
       y = NULL, color = NULL) +
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

ggsave(file.path(output_dir, "Figure5B_TCGA_STAD_multivariate_Cox_forest.pdf"),
       p, width = 10.2, height = 5.8, device = cairo_pdf)
ggsave(file.path(output_dir, "Figure5B_TCGA_STAD_multivariate_Cox_forest.png"),
       p, width = 10.2, height = 5.8, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Figure5B_TCGA_STAD_multivariate_Cox_forest.tiff"),
       p, width = 10.2, height = 5.8, dpi = 300, compression = "lzw", bg = "white")

ph_test <- cox.zph(fit)
capture.output(ph_test,
               file = file.path(output_dir, "Figure5B_TCGA_STAD_multivariate_PH_test.txt"))

notes <- c(
  "TCGA-STAD Figure 5B: multivariable Cox regression forest plot",
  paste("Generated:", Sys.time()),
  paste("Input:", input_file),
  paste("Complete-case samples:", nrow(dat)),
  paste("OS events:", sum(dat$OS_event)),
  "One multivariable Cox model included risk group, age group, gender, and stage simultaneously.",
  "Reference levels: Low risk, age <60, Female, and Stage I.",
  "Red: P < 0.05 and HR > 1; blue: P < 0.05 and HR < 1; gray: P >= 0.05.",
  "A Schoenfeld-residual proportional-hazards test was also saved.",
  "",
  capture.output(summary(fit)),
  "",
  capture.output(sessionInfo())
)
writeLines(notes, file.path(output_dir, "Figure5B_multivariate_run_notes.txt"))

message("Completed: ", output_dir)
}

############################################################
## Panel D: Paired tumor-versus-normal risk-gene expression
############################################################
run_figure5d <- function(input_dir, output_root) {
  model_dir <- file.path(output_root, "AC_prognostic_model")
  expr_file <- file.path(input_dir, "TCGA-STAD.star_tpm.tsv")
  map_file <- file.path(model_dir, "01_candidate_gene_Ensembl_map.csv")
  model_file <- file.path(model_dir, "06_LASSO_complete_model_parameters.csv")
  out_dir <- file.path(output_root, "D_tumor_vs_normal")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
risk_genes <- read.csv(model_file, check.names = FALSE)$Gene
gene_map <- read.csv(map_file, check.names = FALSE) %>%
  filter(Gene %in% risk_genes) %>%
  distinct(Gene, .keep_all = TRUE)
if (!setequal(risk_genes, gene_map$Gene)) stop("Incomplete Ensembl mapping for the risk genes.")

message("Reading TCGA-STAD bulk expression matrix ...")
expr <- fread(expr_file, data.table = FALSE, check.names = FALSE)
names(expr)[1] <- "Ensembl_ID_version"
expr$Ensembl_ID <- sub("\\..*$", "", expr$Ensembl_ID_version)
expr7 <- expr %>%
  inner_join(gene_map, by = "Ensembl_ID") %>%
  select(Gene, starts_with("TCGA-"))
rm(expr); invisible(gc())

sample_cols <- setdiff(names(expr7), "Gene")
expr7[sample_cols] <- lapply(expr7[sample_cols], as.numeric)
expr7 <- expr7 %>%
  group_by(Gene) %>%
  summarise(across(all_of(sample_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
if (!setequal(risk_genes, expr7$Gene)) stop("Not all seven genes were found in expression data.")

long_df <- expr7 %>%
  pivot_longer(-Gene, names_to = "sample_id", values_to = "expression") %>%
  mutate(patient_id = substr(sample_id, 1, 12),
         sample_code = substr(sample_id, 14, 15),
         Tissue = recode(sample_code,
                         `01` = "Primary tumor",
                         `11` = "Adjacent normal",
                         .default = NA_character_)) %>%
  filter(sample_code %in% c("01", "11"), is.finite(expression)) %>%
  group_by(patient_id, sample_code, Tissue, Gene) %>%
  summarise(expression = mean(expression), .groups = "drop") %>%
  mutate(Tissue = factor(Tissue, levels = c("Adjacent normal", "Primary tumor")),
         Gene = factor(Gene, levels = risk_genes))

sample_counts <- long_df %>%
  distinct(patient_id, sample_code, Tissue) %>%
  count(sample_code, Tissue, name = "patient_n")
write.csv(sample_counts, file.path(out_dir, "00_sample_counts.csv"), row.names = FALSE)

paired_wide <- long_df %>%
  select(patient_id, Gene, Tissue, expression) %>%
  pivot_wider(names_from = Tissue, values_from = expression) %>%
  filter(!is.na(`Adjacent normal`), !is.na(`Primary tumor`)) %>%
  mutate(paired_difference = `Primary tumor` - `Adjacent normal`)
write.csv(paired_wide, file.path(out_dir, "01_paired_expression_data.csv"), row.names = FALSE)

paired_stats <- bind_rows(lapply(risk_genes, function(g) {
  d <- paired_wide %>% filter(Gene == g)
  wt <- wilcox.test(d$`Primary tumor`, d$`Adjacent normal`,
                    paired = TRUE, exact = FALSE, conf.int = FALSE)
  data.frame(
    Gene = g,
    Paired_N = nrow(d),
    Normal_median = median(d$`Adjacent normal`, na.rm = TRUE),
    Tumor_median = median(d$`Primary tumor`, na.rm = TRUE),
    Median_paired_difference = median(d$paired_difference, na.rm = TRUE),
    Mean_paired_difference = mean(d$paired_difference, na.rm = TRUE),
    Wilcoxon_V = unname(wt$statistic),
    P_value = wt$p.value
  )
})) %>%
  mutate(Direction = case_when(
           Median_paired_difference > 0 ~ "Up in tumor",
           Median_paired_difference < 0 ~ "Down in tumor",
           TRUE ~ "No median change"
         )) %>%
  arrange(P_value, desc(abs(Median_paired_difference)))
write.csv(paired_stats, file.path(out_dir, "02_primary_paired_differential_statistics.csv"),
          row.names = FALSE)

unpaired_stats <- bind_rows(lapply(risk_genes, function(g) {
  d <- long_df %>% filter(Gene == g)
  x <- d$expression[d$sample_code == "01"]
  y <- d$expression[d$sample_code == "11"]
  wt <- wilcox.test(x, y, paired = FALSE, exact = FALSE)
  data.frame(Gene = g, Tumor_N = length(x), Normal_N = length(y),
             Tumor_median = median(x), Normal_median = median(y),
             Median_difference = median(x) - median(y),
             Wilcoxon_W = unname(wt$statistic), P_value = wt$p.value)
})) %>%
  mutate(Direction = ifelse(Median_difference > 0, "Up in tumor", "Down in tumor")) %>%
  arrange(P_value)
write.csv(unpaired_stats, file.path(out_dir, "03_sensitivity_all_samples_unpaired_statistics.csv"),
          row.names = FALSE)

paired_long <- paired_wide %>%
  select(patient_id, Gene, `Adjacent normal`, `Primary tumor`) %>%
  pivot_longer(c(`Adjacent normal`, `Primary tumor`),
               names_to = "Tissue", values_to = "expression") %>%
  mutate(Tissue = factor(Tissue, levels = c("Adjacent normal", "Primary tumor")),
         Gene = factor(Gene, levels = risk_genes))

annotation_df <- paired_long %>%
  group_by(Gene) %>%
  summarise(y = max(expression, na.rm = TRUE) +
              max(diff(range(expression, na.rm = TRUE)) * 0.10, 0.15), .groups = "drop") %>%
  left_join(paired_stats %>% select(Gene, P_value), by = "Gene") %>%
  mutate(label = ifelse(P_value < 0.001,
                        paste0("P=", format(P_value, scientific = TRUE, digits = 2)),
                        paste0("P=", sprintf("%.3f", P_value))),
         Gene = factor(Gene, levels = risk_genes))

p1 <- ggplot(paired_long, aes(Tissue, expression)) +
  geom_line(aes(group = patient_id), color = "grey60", alpha = 0.38, linewidth = 0.35) +
  geom_boxplot(aes(fill = Tissue), width = 0.48, outlier.shape = NA,
               alpha = 0.72, linewidth = 0.5) +
  geom_point(aes(color = Tissue), position = position_jitter(width = 0.055),
             size = 1.25, alpha = 0.78) +
  geom_text(data = annotation_df, aes(x = 1.5, y = y, label = label),
            inherit.aes = FALSE, size = 3.2, fontface = "bold") +
  facet_wrap(~ Gene, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("Adjacent normal" = "#4DBBD5", "Primary tumor" = "#E64B35")) +
  scale_color_manual(values = c("Adjacent normal" = "#4DBBD5", "Primary tumor" = "#E64B35")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  labs(x = NULL, y = "Expression: log2(TPM + 0.001)",
       title = "TCGA-STAD paired tumor vs adjacent-normal expression of seven risk genes",
       subtitle = paste0("Patient-matched analysis; n = ",
                         length(unique(paired_wide$patient_id)), " pairs; paired Wilcoxon test")) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 25, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

ggsave(file.path(out_dir, "Figure5D_paired_7risk_genes_tumor_vs_normal.pdf"),
       p1, width = 12, height = 7.5)
ggsave(file.path(out_dir, "Figure5D_paired_7risk_genes_tumor_vs_normal.png"),
       p1, width = 12, height = 7.5, dpi = 300)
ggsave(file.path(out_dir, "Figure5D_paired_7risk_genes_tumor_vs_normal.tiff"),
       p1, width = 12, height = 7.5, dpi = 300, compression = "lzw")

effect_df <- paired_stats %>%
  mutate(Gene = factor(Gene, levels = rev(risk_genes)))
p2 <- ggplot(effect_df, aes(Median_paired_difference, Gene, fill = Direction)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
  geom_col(width = 0.65) +
  geom_text(aes(label = paste0("P=", format(P_value, digits = 2, scientific = TRUE))),
            hjust = ifelse(effect_df$Median_paired_difference >= 0, -0.08, 1.08), size = 3.2) +
  scale_fill_manual(values = c("Up in tumor" = "#E64B35", "Down in tumor" = "#4DBBD5")) +
  scale_x_continuous(expand = expansion(mult = c(0.20, 0.25))) +
  labs(x = "Median paired difference: tumor - adjacent normal", y = NULL,
       title = "Direction and magnitude of paired expression differences", fill = NULL) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "top")
ggsave(file.path(out_dir, "Figure5D_paired_7risk_genes_effect_size.pdf"), p2, width = 8, height = 5)
ggsave(file.path(out_dir, "Figure5D_paired_7risk_genes_effect_size.png"), p2,
       width = 8, height = 5, dpi = 300)

capture.output(
  list(R_version = R.version.string,
       risk_genes = risk_genes,
       sample_counts = sample_counts,
       paired_patient_n = length(unique(paired_wide$patient_id)),
       primary_analysis = "Paired Wilcoxon signed-rank test; raw P values are reported.",
       sensitivity_analysis = "All available samples, unpaired Wilcoxon rank-sum test.",
       expression_scale = "UCSC Xena TCGA-STAD STAR TPM matrix: log2(TPM + 0.001).",
       paired_results = paired_stats,
       session_info = sessionInfo()),
  file = file.path(out_dir, "TCGA_STAD_7risk_genes_tumor_vs_normal_run_notes.txt")
)

message("Completed. Output directory: ", out_dir)
print(sample_counts)
print(paired_stats)
}

run_figure5ac(input_dir, figure5_output_dir)
run_figure5b(figure5_output_dir)
run_figure5d(input_dir, figure5_output_dir)

message("Figure 5 panels A-D completed. Output directory: ", figure5_output_dir)
