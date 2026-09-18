############################################################
## Figure S3: TCGA-STAD mutation landscape and spectrum
## Run with: Rscript FigureS3_analysis.R
############################################################

required_packages <- c(
  "dplyr", "ggplot2", "maftools", "patchwork",
  "BSgenome.Hsapiens.UCSC.hg38", "MutationalPatterns", "png"
)
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) stop("缺少R包：", paste(missing_packages, collapse = ", "))

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(maftools)
  library(patchwork)
  library(BSgenome.Hsapiens.UCSC.hg38)
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
mutation_file <- file.path(script_dir, "data", "TCGA-STAD.somaticmutation_wxs.tsv.gz")
clinical_file <- file.path(script_dir, "data", "TCGA-STAD.clinical.gz")
output_dir <- file.path(script_dir, "output", "FigureS3")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(mutation_file)) stop("Input file not found: ", mutation_file)
if (!file.exists(clinical_file)) stop("Input file not found: ", clinical_file)
setwd(output_dir)
mut <- read.delim(gzfile(mutation_file), check.names = FALSE, stringsAsFactors = FALSE)
clinical_df <- read.delim(gzfile(clinical_file), check.names = FALSE, stringsAsFactors = FALSE)

# Convert Xena mutation effects into standard MAF classes. Frameshift direction
# is determined from allele lengths rather than assigning all events as insertions.
maf_df <- mut %>%
  mutate(
    effect_lower = tolower(effect),
    ref_len = nchar(ref),
    alt_len = nchar(alt),
    Variant_Classification = case_when(
      grepl("frameshift_variant", effect_lower) & alt_len < ref_len ~ "Frame_Shift_Del",
      grepl("frameshift_variant", effect_lower) & alt_len > ref_len ~ "Frame_Shift_Ins",
      grepl("frameshift_variant", effect_lower) ~ "Frame_Shift_Del",
      grepl("stop_gained", effect_lower) ~ "Nonsense_Mutation",
      grepl("stop_lost", effect_lower) ~ "Nonstop_Mutation",
      grepl("start_lost", effect_lower) ~ "Translation_Start_Site",
      grepl("splice_acceptor|splice_donor|splice_region", effect_lower) ~ "Splice_Site",
      grepl("inframe_deletion", effect_lower) ~ "In_Frame_Del",
      grepl("inframe_insertion", effect_lower) ~ "In_Frame_Ins",
      grepl("missense_variant|protein_altering_variant", effect_lower) ~ "Missense_Mutation",
      grepl("synonymous_variant|stop_retained_variant", effect_lower) ~ "Silent",
      TRUE ~ NA_character_
    ),
    Variant_Type = case_when(
      ref_len == 1 & alt_len == 1 ~ "SNP",
      ref_len > alt_len ~ "DEL",
      ref_len < alt_len ~ "INS",
      TRUE ~ "ONP"
    )
  ) %>%
  transmute(
    Hugo_Symbol = gene,
    Chromosome = sub("^chr", "", chrom, ignore.case = TRUE),
    Start_Position = as.numeric(start),
    End_Position = as.numeric(end),
    Reference_Allele = ref,
    Tumor_Seq_Allele1 = ref,
    Tumor_Seq_Allele2 = alt,
    Tumor_Sample_Barcode = substr(sample, 1, 16),
    Variant_Classification,
    Variant_Type
  ) %>%
  filter(!is.na(Variant_Classification), !is.na(Hugo_Symbol), Hugo_Symbol != "") %>%
  distinct()

clinical_anno <- clinical_df %>%
  transmute(
    Tumor_Sample_Barcode = substr(sample, 1, 16),
    Gender = case_when(
      tolower(gender.demographic) == "male" ~ "Male",
      tolower(gender.demographic) == "female" ~ "Female",
      TRUE ~ NA_character_
    ),
    Age_value = suppressWarnings(as.numeric(age_at_index.demographic)),
    Stage_raw = ajcc_pathologic_stage.diagnoses,
    T_raw = ajcc_pathologic_t.diagnoses,
    N_raw = ajcc_pathologic_n.diagnoses,
    M_raw = ajcc_pathologic_m.diagnoses,
    Sample_Type = sample_type.samples
  ) %>%
  mutate(
    Age = case_when(Age_value < 60 ~ "<60", Age_value >= 60 ~ ">=60", TRUE ~ NA_character_),
    pT = case_when(
      grepl("Tis|T0", T_raw, ignore.case = TRUE) ~ "T0/is",
      grepl("T1", T_raw, ignore.case = TRUE) ~ "T1",
      grepl("T2", T_raw, ignore.case = TRUE) ~ "T2",
      grepl("T3", T_raw, ignore.case = TRUE) ~ "T3",
      grepl("T4", T_raw, ignore.case = TRUE) ~ "T4",
      TRUE ~ NA_character_
    ),
    pN = case_when(
      grepl("N0", N_raw, ignore.case = TRUE) ~ "N0",
      grepl("N1", N_raw, ignore.case = TRUE) ~ "N1",
      grepl("N2", N_raw, ignore.case = TRUE) ~ "N2",
      grepl("N3", N_raw, ignore.case = TRUE) ~ "N3",
      TRUE ~ NA_character_
    ),
    pM = case_when(
      grepl("M0", M_raw, ignore.case = TRUE) ~ "M0",
      grepl("M1", M_raw, ignore.case = TRUE) ~ "M1",
      TRUE ~ NA_character_
    ),
    Stage = case_when(
      grepl("Stage IV", Stage_raw, ignore.case = TRUE) ~ "IV",
      grepl("Stage III", Stage_raw, ignore.case = TRUE) ~ "III",
      grepl("Stage II", Stage_raw, ignore.case = TRUE) ~ "II",
      grepl("Stage I", Stage_raw, ignore.case = TRUE) ~ "I",
      TRUE ~ NA_character_
    )
  ) %>%
  select(Tumor_Sample_Barcode, Age, Gender, pT, pN, pM, Stage, Sample_Type) %>%
  distinct(Tumor_Sample_Barcode, .keep_all = TRUE)

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
)

ann_colors <- list(
  Age = c("<60" = "#B2D990", ">=60" = "#F0A0A0", "Unknown" = "#D9D9D9"),
  Gender = c("Male" = "#82D4D7", "Female" = "#DC7579", "Unknown" = "#D9D9D9"),
  pT = c("T0/is" = "#F7F7F7", "T1" = "#D9F0D3", "T2" = "#C2E2C2", "T3" = "#B9DFB9", "T4" = "#AFD1BF", "Unknown" = "#D9D9D9"),
  pN = c("N0" = "#F0F0F0", "N1" = "#BCBDDC", "N2" = "#807DBA", "N3" = "#9898DC", "N3A" = "#897CD3", "N3B" = "#5E379D", "NX" = "#AAAAAA", "Unknown" = "#D9D9D9"),
  pM = c("M0" = "#82B6E8", "M1" = "#EEAAAA", "MX" = "#AAAAAA", "Unknown" = "#D9D9D9"),
  Stage = c("I" = "#FFF3B0", "II" = "#FBE3D2", "III" = "#F8D0B0", "IV" = "#F2B382", "Unknown" = "#D9D9D9"),

  # Age = c("<60" = "#E7A6B0", ">=60" = "#B94747"),
  # Gender = c("female" = "#E64B35", "male" = "#4DBBD5", "not reported" = "#BDBDBD"),
  # pT = c("T0/is" = "#F2E2F2", "T1" = "#D8B7D8", "T2" = "#B987C0", "T3" = "#8B4B9C", "T4" = "#5A236E"),
  # pN = c("N0" = "#D9D9D9", "N1" = "#B6A6D8", "N2" = "#8C6BB1", "N3" = "#5E3C99"),
  # pM = c("M0" = "#BDBDBD", "M1" = "#E64B35"),
  # Stage = c("I" = "#DDECC9", "II" = "#A8DDB5", "III" = "#43A2CA", "IV" = "#0868AC"),
  Sample_Type = c("Primary Tumor" = "#4DBBD5", "Solid Tissue Normal" = "#7A7A7A", "Metastatic" = "#E64B35")
)

maf <- read.maf(maf = maf_df, clinicalData = clinical_anno, verbose = FALSE)
clinical_features <- c("Age", "Gender", "pT", "pN", "pM", "Stage", "Sample_Type")
clinical_features <- clinical_features[vapply(clinical_features, function(x) {
  vals <- maf@clinical.data[[x]]
  length(unique(vals[!is.na(vals) & vals != ""])) > 0
}, logical(1))]
ann_colors_use <- ann_colors[clinical_features]

# Write device output to a temporary file first. This prevents an interrupted or
# failed plot from leaving behind a zero-byte/corrupt file with the final name.
save_oncoplot <- function(filename, device = c("pdf", "png")) {
  device <- match.arg(device)
  extension <- paste0(".", device)
  tmp <- tempfile(pattern = "FigureS3A_", tmpdir = output_dir, fileext = extension)
  device_open <- FALSE
  completed <- FALSE

  on.exit({
    if (device_open && grDevices::dev.cur() > 1L) {
      try(grDevices::dev.off(), silent = TRUE)
    }
    if (!completed && file.exists(tmp)) unlink(tmp)
  }, add = TRUE)

  if (device == "pdf") {
    grDevices::cairo_pdf(tmp, width = 13, height = 9)
  } else {
    grDevices::png(tmp, width = 3900, height = 2700, res = 300)
  }
  device_open <- TRUE

  oncoplot(
    maf = maf, top = 40, colors = mut_col, removeNonMutated = TRUE,
    clinicalFeatures = clinical_features, annotationColor = ann_colors_use,
    sortByAnnotation = TRUE, showTumorSampleBarcodes = FALSE,
    drawRowBar = TRUE, drawColBar = TRUE, fontSize = 0.65,
    titleText = "Clinical annotated TCGA-STAD mutation landscape"
  )

  grDevices::dev.off()
  device_open <- FALSE
  if (!file.copy(tmp, filename, overwrite = TRUE)) {
    stop("Could not replace output file: ", filename,
         ". Close any PDF viewer or older R process using it, then rerun.")
  }
  unlink(tmp)
  completed <- TRUE
}

# Figure S3A: clinical-annotated TCGA-STAD oncoplot.
save_oncoplot("FigureS3A_TCGA_STAD_clinical_annotated_oncoplot.pdf", "pdf")
save_oncoplot("FigureS3A_TCGA_STAD_clinical_annotated_oncoplot.png", "png")

# Some Windows PDF viewers cannot render the very complex vector output made by
# maftools::oncoplot. Rebuild A as a simple, broadly compatible raster-backed PDF
# from the full-resolution 300 dpi PNG.
png_to_compatible_pdf <- function(png_file, pdf_file) {
  if (!requireNamespace("png", quietly = TRUE)) {
    stop("Package 'png' is required to create the compatible Figure S3A PDF.")
  }
  img <- png::readPNG(png_file)
  tmp <- tempfile(pattern = "FigureS3A_compatible_", tmpdir = output_dir,
                  fileext = ".pdf")
  device_open <- FALSE
  completed <- FALSE
  on.exit({
    if (device_open && grDevices::dev.cur() > 1L) {
      try(grDevices::dev.off(), silent = TRUE)
    }
    if (!completed && file.exists(tmp)) unlink(tmp)
  }, add = TRUE)

  grDevices::cairo_pdf(tmp, width = 13, height = 9, onefile = TRUE)
  device_open <- TRUE
  grid::grid.newpage()
  grid::grid.raster(img, width = grid::unit(1, "npc"),
                    height = grid::unit(1, "npc"), interpolate = TRUE)
  grDevices::dev.off()
  device_open <- FALSE

  if (!file.copy(tmp, pdf_file, overwrite = TRUE)) {
    stop("Could not replace compatible PDF: ", pdf_file)
  }
  unlink(tmp)
  completed <- TRUE
}

png_to_compatible_pdf(
  "FigureS3A_TCGA_STAD_clinical_annotated_oncoplot.png",
  "FigureS3A_TCGA_STAD_clinical_annotated_oncoplot.pdf"
)

# SNVs and hg38 trinucleotide contexts.
snv_df <- maf_df %>%
  filter(
    nchar(Reference_Allele) == 1, nchar(Tumor_Seq_Allele2) == 1,
    Reference_Allele %in% c("A", "C", "G", "T"),
    Tumor_Seq_Allele2 %in% c("A", "C", "G", "T"),
    Reference_Allele != Tumor_Seq_Allele2,
    Chromosome %in% c(as.character(1:22), "X", "Y"),
    Start_Position > 1
  ) %>%
  mutate(chr_use = paste0("chr", Chromosome))

snv_df$tri_context <- as.character(BSgenome::getSeq(
  BSgenome.Hsapiens.UCSC.hg38,
  names = snv_df$chr_use,
  start = snv_df$Start_Position - 1,
  end = snv_df$Start_Position + 1
))
snv_df$tri_context <- toupper(snv_df$tri_context)

comp_base <- function(x) chartr("ACGT", "TGCA", x)
revcomp_tri <- function(x) vapply(x, function(z) {
  paste0(rev(strsplit(comp_base(z), "", fixed = TRUE)[[1]]), collapse = "")
}, character(1))

snv_df <- snv_df %>%
  mutate(
    ref_use = ifelse(Reference_Allele %in% c("C", "T"), Reference_Allele, comp_base(Reference_Allele)),
    alt_use = ifelse(Reference_Allele %in% c("C", "T"), Tumor_Seq_Allele2, comp_base(Tumor_Seq_Allele2)),
    tri_use = ifelse(Reference_Allele %in% c("C", "T"), tri_context, revcomp_tri(tri_context)),
    mut_type = paste0(ref_use, ">", alt_use),
    left_base = substr(tri_use, 1, 1),
    right_base = substr(tri_use, 3, 3),
    context_label = paste0(left_base, "[", mut_type, "]", right_base)
  ) %>%
  filter(
    mut_type %in% c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
    left_base %in% c("A", "C", "G", "T"), right_base %in% c("A", "C", "G", "T")
  )

mut_levels <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
all_context <- expand.grid(
  mut_type = mut_levels, left_base = c("A", "C", "G", "T"),
  right_base = c("A", "C", "G", "T"), stringsAsFactors = FALSE
) %>%
  mutate(context_label = paste0(left_base, "[", mut_type, "]", right_base))

spectrum_df <- snv_df %>%
  count(mut_type, context_label, name = "count") %>%
  right_join(all_context, by = c("mut_type", "context_label")) %>%
  mutate(count = ifelse(is.na(count), 0L, count)) %>%
  select(mut_type, left_base, right_base, context_label, count) %>%
  arrange(factor(mut_type, mut_levels), left_base, right_base) %>%
  mutate(
    relative_contribution = count / sum(count),
    mut_type = factor(mut_type, levels = mut_levels),
    context_label = factor(context_label, levels = context_label)
  )

snv_cols <- c("C>A" = "#2CA9E1", "C>G" = "#000000", "C>T" = "#E64B35",
              "T>A" = "#BDBDBD", "T>C" = "#8BC34A", "T>G" = "#F4A3A3")

pB <- ggplot(spectrum_df, aes(context_label, relative_contribution, fill = mut_type)) +
  geom_col(width = 0.78, color = "black", linewidth = 0.15) +
  facet_grid(. ~ mut_type, scales = "free_x", space = "free_x") +
  # Match MutationalPatterns::plot_96_profile(condensed = TRUE): the facet
  # strip gives the substitution class, so the x axis only shows flanking bases.
  scale_x_discrete(labels = function(x) paste0(substr(x, 1, 1), ".", substr(x, nchar(x), nchar(x)))) +
  scale_fill_manual(values = snv_cols, drop = FALSE) +
  labs(x = NULL, y = "Relative contribution") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none", panel.grid = element_blank(),
        strip.background = element_rect(fill = "grey85", color = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5.5))

variant_df <- maf_df %>%
  filter(Variant_Classification != "Silent") %>%
  count(Variant_Classification, name = "count") %>% arrange(count)
pC <- ggplot(variant_df, aes(count, reorder(Variant_Classification, count), fill = Variant_Classification)) +
  geom_col(width = 0.75, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = mut_col, drop = FALSE) +
  labs(title = "Variant type", x = "Mutation count", y = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"))

snv_type_df <- snv_df %>% count(mut_type, name = "count") %>%
  mutate(mut_type = factor(mut_type, levels = rev(mut_levels)))
pD <- ggplot(snv_type_df, aes(count, mut_type, fill = mut_type)) +
  geom_col(width = 0.75, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = snv_cols, drop = FALSE) +
  labs(title = "SNV type", x = "Mutation count", y = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("FigureS3B_TCGA_STAD_96_channel_spectrum.pdf", pB, width = 13, height = 5, device = cairo_pdf)
ggsave("FigureS3B_TCGA_STAD_96_channel_spectrum.png", pB, width = 13, height = 5, dpi = 300, bg = "white")
ggsave("FigureS3C_TCGA_STAD_variant_type.pdf", pC, width = 6, height = 4.5, device = cairo_pdf)
ggsave("FigureS3C_TCGA_STAD_variant_type.png", pC, width = 6, height = 4.5, dpi = 300, bg = "white")
ggsave("FigureS3D_TCGA_STAD_SNV_type.pdf", pD, width = 5.5, height = 4.5, device = cairo_pdf)
ggsave("FigureS3D_TCGA_STAD_SNV_type.png", pD, width = 5.5, height = 4.5, dpi = 300, bg = "white")

combined_BCD <- pB / (pC | pD) +
  plot_annotation(tag_levels = list(c("B", "C", "D"))) &
  theme(plot.tag = element_text(face = "bold", size = 15))
ggsave("FigureS3BCD_TCGA_STAD_mutation_summaries.pdf", combined_BCD, width = 13, height = 9, device = cairo_pdf)
ggsave("FigureS3BCD_TCGA_STAD_mutation_summaries.png", combined_BCD, width = 13, height = 9, dpi = 300, bg = "white")

write.csv(maf_df, "FigureS3_processed_TCGA_STAD_MAF.csv", row.names = FALSE)
write.csv(clinical_anno, "FigureS3_processed_TCGA_STAD_clinical.csv", row.names = FALSE)
write.csv(spectrum_df, "FigureS3_96_channel_spectrum_data.csv", row.names = FALSE)
write.csv(variant_df, "FigureS3_variant_type_counts.csv", row.names = FALSE)
write.csv(snv_type_df, "FigureS3_SNV_type_counts.csv", row.names = FALSE)

writeLines(c(
  paste0("Raw mutation rows: ", nrow(mut)),
  paste0("MAF-compatible coding rows: ", nrow(maf_df)),
  paste0("Mutated samples in MAF: ", length(unique(maf_df$Tumor_Sample_Barcode))),
  paste0("Context-qualified SNVs: ", nrow(snv_df)),
  "Reference genome for trinucleotide extraction: hg38",
  "Figure S3A displays the 40 most frequently mutated genes."
), "FigureS3_run_notes.txt")

message("Figure S3 outputs written to: ", output_dir)
