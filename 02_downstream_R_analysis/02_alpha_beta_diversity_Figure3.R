# =============================================================================
# Figure 3. Alpha and beta diversity of the milk microbiome by lactation stage
# =============================================================================
#
# Purpose
# -------
# Reproduce the alpha- and beta-diversity analyses used for Figure 3 of the
# manuscript.
#
# Study design
# ------------
# - 42 milk samples
# - 21 Early-lactation samples
# - 21 Late-lactation samples
# - 21 matched Early-Late animal pairs
# - Samples collected from Farm A and Farm B
#
# Analyses
# --------
# Alpha diversity:
# - Observed richness
# - Chao1 richness
# - Shannon diversity
# - Simpson diversity
# - Paired Wilcoxon signed-rank tests
# - Benjamini-Hochberg false-discovery-rate correction
#
# Beta diversity:
# - Bray-Curtis dissimilarity
# - Binary Jaccard dissimilarity
# - Principal coordinates analysis (PCoA)
# - Farm-adjusted PERMANOVA with permutations restricted within animal pairs
# - ANOSIM as a complementary rank-based test
# - PERMDISP as a dispersion diagnostic
#
# Expected input files
# --------------------
# data/processed/feature-table_milk_42samples.tsv
# data/metadata/metadata_milk.csv
#
# Main outputs
# ------------
# results/figure3_alpha_beta/figures/
# results/figure3_alpha_beta/tables/
#
# Reproducibility
# ---------------
# Run this script from the root of the GitHub repository.
# Random seed and permutation count are fixed below.
#
# IMPORTANT:
# Confirm whether the input feature table contains ASVs or species-collapsed
# abundances. The manuscript figure title should use "species level" only if the
# input table has already been collapsed to species-level features.
# =============================================================================

# =============================================================================
# 0. PACKAGES
# =============================================================================

required_packages <- c("tidyverse", "vegan", "patchwork", "stringr")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running this script."
  )
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(patchwork)
  library(stringr)
})


# =============================================================================
# 1. REPOSITORY-RELATIVE FILE PATHS
# =============================================================================

# Run from the root of the repository.
PROJECT_DIR <- "."

FEATURE_FILE <- file.path(
  PROJECT_DIR,
  "data",
  "processed",
  "feature-table_milk_42samples.tsv"
)

METADATA_FILE <- file.path(
  PROJECT_DIR,
  "data",
  "metadata",
  "metadata_milk.csv"
)

OUT_DIR <- file.path(PROJECT_DIR, "results", "figure3_alpha_beta")
OUT_FIG <- file.path(OUT_DIR, "figures")
OUT_TAB <- file.path(OUT_DIR, "tables")

dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(FEATURE_FILE)) {
  stop("Feature table not found: ", FEATURE_FILE)
}

if (!file.exists(METADATA_FILE)) {
  stop("Metadata file not found: ", METADATA_FILE)
}


# =============================================================================
# 2. REPRODUCIBILITY + FIGURE STYLE
# =============================================================================

set.seed(12345)
N_PERM <- 9999

STAGE_COLORS <- c(
  "Early" = "#3274A1",
  "Late"  = "#E0447E"
)

FONT_TITLE        <- 32
FONT_FACET        <- 26
FONT_AXIS_TITLE   <- 28
FONT_AXIS_TEXT    <- 24
FONT_LEGEND_TITLE <- 25
FONT_LEGEND_TEXT  <- 23
FONT_PANEL_TAG    <- 30
FONT_STAT         <- 7.2


# =============================================================================
# 3. READ MILK-ONLY FEATURE TABLE
# =============================================================================
# This TSV was generated from BIOM in R, so skip = 1 is NOT used.
# Rows = features; columns = samples.
# Verify whether these features are ASVs or species-level collapsed taxa.

feature_counts <- read.delim(
  FEATURE_FILE,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = ""
)

feature_id <- as.character(feature_counts[[1]])
count_mat <- as.matrix(feature_counts[, -1, drop = FALSE])
mode(count_mat) <- "numeric"
rownames(count_mat) <- feature_id

# Feature-table QC
if (anyNA(count_mat)) stop("NA values produced when feature counts were converted to numeric.")
if (any(count_mat < 0)) stop("Negative counts detected in feature table.")

integer_check <- all(abs(count_mat - round(count_mat)) < 1e-8)
if (!integer_check) {
  stop("Feature table is not integer-like; richness metrics require raw count data.")
}

cat("\n============================================================\n")
cat("FEATURE-TABLE QC\n")
cat("============================================================\n")
cat("ASVs/features:", nrow(count_mat), "\n")
cat("Samples:", ncol(count_mat), "\n")


# =============================================================================
# 4. READ AND PREPARE METADATA
# =============================================================================

meta <- read.csv(
  METADATA_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Ignore unnamed empty columns at the end of the CSV.
meta <- meta[, nzchar(colnames(meta)), drop = FALSE]

required_meta_cols <- c("SampleID", "Sample_Name", "Lactation", "Farm")
missing_meta_cols <- setdiff(required_meta_cols, colnames(meta))
if (length(missing_meta_cols) > 0) {
  stop("Missing metadata columns: ", paste(missing_meta_cols, collapse = ", "))
}

meta$SampleID    <- str_trim(as.character(meta$SampleID))
meta$Sample_Name <- str_trim(as.character(meta$Sample_Name))
meta$Lactation   <- str_trim(as.character(meta$Lactation))
meta$Farm        <- str_trim(as.character(meta$Farm))


# =============================================================================
# 5. DEFINE PAIR ID, LACTATION STAGE AND FARM LABEL
# =============================================================================
# Prefer a dedicated animal/pair column if present. Otherwise derive PairID
# from Sample_Name using the same suffix-removal rule as the previous analysis.

pair_candidates <- c("AnimalID", "Animal_ID", "CowID", "Cow_ID", "PairID")
available_pair_col <- pair_candidates[pair_candidates %in% colnames(meta)]

if (length(available_pair_col) > 0) {
  pair_column <- available_pair_col[1]
  meta$PairID <- str_trim(as.character(meta[[pair_column]]))
  cat("\nPair ID source:", pair_column, "\n")
} else {
  meta$PairID <- str_remove(meta$Sample_Name, "_(AB|NI)2?$")
  cat("\nPair ID derived from Sample_Name.\n")
}

# Add Farm to avoid accidental pairing of identical IDs across farms.
meta$PairUID <- paste(meta$Farm, meta$PairID, sep = "_")

meta$Stage_full <- case_when(
  tolower(meta$Lactation) == "early" ~ "Early",
  tolower(meta$Lactation) == "late"  ~ "Late",
  TRUE ~ NA_character_
)

meta$Stage_full <- factor(meta$Stage_full, levels = c("Early", "Late"))

meta$Farm_label <- case_when(
  meta$Farm == "A" ~ "Farm A",
  meta$Farm == "B" ~ "Farm B",
  TRUE ~ NA_character_
)

meta$Farm_label <- factor(meta$Farm_label, levels = c("Farm A", "Farm B"))


# =============================================================================
# 6. MATCH FEATURE-TABLE SAMPLES TO METADATA
# =============================================================================

feature_samples <- colnames(count_mat)
matched_samples <- intersect(feature_samples, meta$SampleID)

cat("\n============================================================\n")
cat("MATCHED-SAMPLE QC\n")
cat("============================================================\n")
cat("Feature-table samples:", length(feature_samples), "\n")
cat("Metadata samples:", nrow(meta), "\n")
cat("Matched samples:", length(matched_samples), "\n")

if (length(matched_samples) != 42) {
  stop("Expected exactly 42 matched milk samples, but found ", length(matched_samples), ".")
}

# Put metadata in exactly the same order as feature-table sample columns.
meta_div <- meta[match(feature_samples, meta$SampleID), , drop = FALSE]

if (anyNA(meta_div$SampleID)) {
  stop("One or more feature-table samples could not be matched to metadata.")
}

stopifnot(identical(feature_samples, meta_div$SampleID))


# =============================================================================
# 7. CONFIRM 21 EARLY + 21 LATE SAMPLES
# =============================================================================

cat("\n============================================================\n")
cat("LACTATION-STAGE QC\n")
cat("============================================================\n")
print(table(meta_div$Stage_full, useNA = "ifany"))

n_early <- sum(meta_div$Stage_full == "Early", na.rm = TRUE)
n_late  <- sum(meta_div$Stage_full == "Late",  na.rm = TRUE)

cat("\nEarly samples:", n_early, "\n")
cat("Late samples:", n_late, "\n")

if (n_early != 21 || n_late != 21) {
  stop("Expected 21 Early and 21 Late samples.")
}


# =============================================================================
# 8. CONFIRM 21 COMPLETE EARLY-LATE ANIMAL PAIRS
# =============================================================================
# A complete pair = exactly two samples, one Early and one Late.

pair_qc <- meta_div %>%
  group_by(PairUID, Farm_label) %>%
  summarise(
    N_samples = n_distinct(SampleID),
    Has_Early = any(Stage_full == "Early"),
    Has_Late  = any(Stage_full == "Late"),
    Complete_pair = (N_samples == 2 && Has_Early && Has_Late),
    .groups = "drop"
  )

n_complete_pairs <- sum(pair_qc$Complete_pair)

cat("\n============================================================\n")
cat("ANIMAL-PAIR QC\n")
cat("============================================================\n")
cat("Complete Early-Late pairs:", n_complete_pairs, "\n")
cat("Samples in complete pairs:", n_complete_pairs * 2, "\n")
cat("Samples outside complete pairs:", nrow(meta_div) - n_complete_pairs * 2, "\n")
print(pair_qc)

write.csv(pair_qc, file.path(OUT_TAB, "00_Pair_QC.csv"), row.names = FALSE)

if (n_complete_pairs != 21) {
  stop("Expected exactly 21 complete Early-Late animal pairs.")
}

cat("\nPASS: 42 samples = 21 Early + 21 Late = 21 complete animal pairs.\n")


# =============================================================================
# 9. BUILD SAMPLE-BY-ASV MATRIX
# =============================================================================
# vegan expects samples as rows and ASVs as columns.

count_sample_asv <- t(
  count_mat[, meta_div$SampleID, drop = FALSE]
)

count_sample_asv <- count_sample_asv[
  meta_div$SampleID,
  ,
  drop = FALSE
]

stopifnot(identical(rownames(count_sample_asv), meta_div$SampleID))

# Remove ASVs absent from all 42 milk samples.
count_sample_asv <- count_sample_asv[
  ,
  colSums(count_sample_asv) > 0,
  drop = FALSE
]

sample_library_size <- rowSums(count_sample_asv)
if (any(sample_library_size <= 0)) stop("One or more samples have zero total reads.")

cat("\nASVs retained after milk-only filtering:", ncol(count_sample_asv), "\n")


# =============================================================================
# 10. SAVE LIBRARY-SIZE / SEQUENCING-DEPTH QC
# =============================================================================

library_size_table <- tibble(
  SampleID = rownames(count_sample_asv),
  Total_reads = sample_library_size
) %>%
  left_join(
    meta_div %>% select(SampleID, Farm_label, Stage_full, PairUID),
    by = "SampleID"
  )

write.csv(
  library_size_table,
  file.path(OUT_TAB, "00_Sample_Library_Sizes.csv"),
  row.names = FALSE
)


# =============================================================================
# 11. CALCULATE ALPHA DIVERSITY
# =============================================================================
# Alpha diversity measures diversity within each sample.
# Observed = directly observed ASV richness
# Chao1    = richness estimator sensitive to rare ASVs
# Shannon  = richness + evenness
# Simpson  = diversity weighted more toward dominant taxa

Observed <- vegan::specnumber(count_sample_asv)

richness_estimates <- vegan::estimateR(count_sample_asv)
Chao1 <- richness_estimates["S.chao1", ]

Shannon <- vegan::diversity(count_sample_asv, index = "shannon")
Simpson <- vegan::diversity(count_sample_asv, index = "simpson")


# =============================================================================
# 12. CREATE AND SAVE ALPHA-DIVERSITY VALUES TABLE
# =============================================================================

alpha_df <- meta_div %>%
  select(
    SampleID,
    Sample_Name,
    PairID,
    PairUID,
    Farm,
    Farm_label,
    Stage_full
  ) %>%
  mutate(
    Observed = as.numeric(Observed),
    Chao1    = as.numeric(Chao1),
    Shannon  = as.numeric(Shannon),
    Simpson  = as.numeric(Simpson)
  )

write.csv(
  alpha_df,
  file.path(OUT_TAB, "01_Alpha_Diversity_Values.csv"),
  row.names = FALSE
)


# =============================================================================
# 13. PAIRED ALPHA-DIVERSITY TESTING
# =============================================================================
# Early and Late samples come from the same animals, so a paired Wilcoxon
# signed-rank test is used for each alpha-diversity metric.

complete_pair_ids <- pair_qc %>%
  filter(Complete_pair) %>%
  pull(PairUID)

alpha_complete <- alpha_df %>%
  filter(PairUID %in% complete_pair_ids)

paired_wilcoxon <- function(data, metric) {

  paired <- data %>%
    select(PairUID, Stage_full, all_of(metric)) %>%
    pivot_wider(
      names_from = Stage_full,
      values_from = all_of(metric)
    ) %>%
    filter(!is.na(Early), !is.na(Late))

  if (nrow(paired) < 2) {
    return(
      tibble(
        Metric = metric,
        N_pairs = nrow(paired),
        W = NA_real_,
        P_value = NA_real_
      )
    )
  }

  test <- wilcox.test(
    paired$Late,
    paired$Early,
    paired = TRUE,
    exact = FALSE
  )

  tibble(
    Metric = metric,
    N_pairs = nrow(paired),
    W = unname(test$statistic),
    P_value = test$p.value
  )
}

alpha_metrics <- c("Observed", "Chao1", "Shannon", "Simpson")

alpha_stats <- bind_rows(
  lapply(
    alpha_metrics,
    function(metric) paired_wilcoxon(alpha_complete, metric)
  )
)


# =============================================================================
# 14. APPLY BENJAMINI-HOCHBERG FDR CORRECTION
# =============================================================================
# Four alpha-diversity tests are performed, so BH-FDR is applied to control
# the false-discovery rate across the four alpha metrics.

alpha_stats <- alpha_stats %>%
  mutate(
    FDR = p.adjust(P_value, method = "BH"),
    Significant_FDR_0.05 = ifelse(
      !is.na(FDR) & FDR < 0.05,
      "Yes",
      "No"
    )
  )

cat("\n============================================================\n")
cat("ALPHA-DIVERSITY PAIRED WILCOXON RESULTS\n")
cat("============================================================\n")
print(alpha_stats)

write.csv(
  alpha_stats,
  file.path(OUT_TAB, "02_Alpha_Diversity_Paired_Wilcoxon.csv"),
  row.names = FALSE
)


# =============================================================================
# 15. PREPARE ALPHA PLOTTING DATA
# =============================================================================

alpha_long <- alpha_df %>%
  pivot_longer(
    cols = all_of(alpha_metrics),
    names_to = "Index",
    values_to = "Value"
  ) %>%
  mutate(
    Index = factor(
      Index,
      levels = c("Observed", "Chao1", "Shannon", "Simpson")
    )
  )


# =============================================================================
# 16. ALPHA PLOT FUNCTION
# =============================================================================

plot_alpha_panel <- function(metric, y_label, panel_tag) {

  d <- alpha_long %>% filter(Index == metric)

  ggplot(
    d,
    aes(x = Stage_full, y = Value, fill = Stage_full)
  ) +
    geom_boxplot(
      width = 0.48,
      linewidth = 1.1,
      outlier.shape = 16,
      outlier.size = 3.2
    ) +
    facet_wrap(~ Farm_label, nrow = 1) +
    scale_fill_manual(values = STAGE_COLORS, name = "Lactation") +
    labs(
      tag = panel_tag,
      title = metric,
      x = "Lactation",
      y = y_label
    ) +
    theme_bw(base_size = FONT_AXIS_TEXT) +
    theme(
      plot.tag = element_text(size = FONT_PANEL_TAG, face = "bold.italic"),
      plot.title = element_text(size = FONT_TITLE, face = "bold"),
      strip.text = element_text(size = FONT_FACET, face = "bold"),
      axis.title.x = element_text(size = FONT_AXIS_TITLE, face = "bold"),
      axis.title.y = element_text(size = FONT_AXIS_TITLE, face = "bold"),
      axis.text.x = element_text(size = FONT_AXIS_TEXT, face = "bold"),
      axis.text.y = element_text(size = FONT_AXIS_TEXT, face = "bold"),
      legend.title = element_text(size = FONT_LEGEND_TITLE, face = "bold"),
      legend.text = element_text(size = FONT_LEGEND_TEXT, face = "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(12, 12, 12, 12)
    )
}


# =============================================================================
# 17. CREATE ALPHA PANELS A-D
# =============================================================================

p_observed <- plot_alpha_panel("Observed", "Observed ASVs", "A")
p_chao1    <- plot_alpha_panel("Chao1", "Chao1 richness estimate", "B")
p_shannon  <- plot_alpha_panel("Shannon", "Shannon diversity index", "C")
p_simpson  <- plot_alpha_panel("Simpson", "Simpson diversity index", "D")


# =============================================================================
# 18. COMBINE AND SAVE ALPHA FIGURE
# =============================================================================

p_alpha <- (
  (p_observed | p_shannon) /
    (p_chao1 | p_simpson) +
    plot_layout(guides = "collect")
) &
  theme(legend.position = "right")

print(p_alpha)

ggsave(
  file.path(OUT_FIG, "01_Alpha_Diversity_A-D_Shared_Legend.png"),
  p_alpha,
  width = 22,
  height = 17,
  units = "in",
  dpi = 300,
  limitsize = FALSE
)

ggsave(
  file.path(OUT_FIG, "01_Alpha_Diversity_A-D_Shared_Legend.pdf"),
  p_alpha,
  width = 22,
  height = 17,
  units = "in",
  device = cairo_pdf,
  limitsize = FALSE
)


# =============================================================================
# 19. CALCULATE BETA-DIVERSITY DISTANCES
# =============================================================================
# Beta diversity measures differences in community composition between samples.
# Bray-Curtis = abundance-based dissimilarity.
# Binary Jaccard = presence/absence-based dissimilarity.

# Convert to relative abundance for Bray-Curtis.
beta_rel <- sweep(
  count_sample_asv,
  1,
  rowSums(count_sample_asv),
  "/"
)

bray_dist <- vegan::vegdist(
  beta_rel,
  method = "bray"
)

jaccard_dist <- vegan::vegdist(
  count_sample_asv,
  method = "jaccard",
  binary = TRUE
)


# =============================================================================
# 20. PCoA FUNCTION
# =============================================================================
# PCoA represents the distance matrix in a two-dimensional ordination space.

make_pcoa <- function(dist_object, distance_name) {

  ord <- cmdscale(
    dist_object,
    eig = TRUE,
    k = 2,
    add = TRUE
  )

  positive_eigenvalues <- ord$eig[ord$eig > 0]

  pc1_percent <- (ord$eig[1] / sum(positive_eigenvalues)) * 100
  pc2_percent <- (ord$eig[2] / sum(positive_eigenvalues)) * 100

  coordinates <- as.data.frame(
    ord$points[, 1:2, drop = FALSE]
  )

  colnames(coordinates) <- c("PCoA1", "PCoA2")
  coordinates$SampleID <- rownames(coordinates)

  coordinates <- coordinates %>%
    left_join(
      meta_div %>%
        select(
          SampleID,
          Sample_Name,
          PairID,
          PairUID,
          Farm,
          Farm_label,
          Stage_full
        ),
      by = "SampleID"
    )

  list(
    data = coordinates,
    pc1 = pc1_percent,
    pc2 = pc2_percent,
    name = distance_name
  )
}


# =============================================================================
# 21. RUN PCoA FOR BRAY-CURTIS AND JACCARD
# =============================================================================

bray_pcoa <- make_pcoa(bray_dist, "Bray-Curtis")
jaccard_pcoa <- make_pcoa(jaccard_dist, "Jaccard")

write.csv(
  bray_pcoa$data,
  file.path(OUT_TAB, "03_BrayCurtis_PCoA_Coordinates.csv"),
  row.names = FALSE
)

write.csv(
  jaccard_pcoa$data,
  file.path(OUT_TAB, "04_Jaccard_PCoA_Coordinates.csv"),
  row.names = FALSE
)


# =============================================================================
# 22. FARM-ADJUSTED PAIRED PERMANOVA
# =============================================================================
# Model: distance ~ Farm_label + Stage_full
# Farm is included to account for farm-related community variation.
# Permutations are restricted within PairUID to preserve the paired design.
# by = "margin" gives a separate Stage_full test after accounting for Farm.

set.seed(12345)

bray_permanova <- vegan::adonis2(
  bray_dist ~ Farm_label + Stage_full,
  data = meta_div,
  permutations = N_PERM,
  strata = meta_div$PairUID,
  by = "margin"
)

set.seed(12345)

jaccard_permanova <- vegan::adonis2(
  jaccard_dist ~ Farm_label + Stage_full,
  data = meta_div,
  permutations = N_PERM,
  strata = meta_div$PairUID,
  by = "margin"
)

cat("\n============================================================\n")
cat("BRAY-CURTIS PERMANOVA\n")
cat("============================================================\n")
print(bray_permanova)

cat("\n============================================================\n")
cat("JACCARD PERMANOVA\n")
cat("============================================================\n")
print(jaccard_permanova)


# =============================================================================
# 23. EXTRACT LACTATION-STAGE PERMANOVA RESULTS
# =============================================================================

extract_stage_permanova <- function(result, distance_name) {

  result_df <- as.data.frame(result)

  if (!"Stage_full" %in% rownames(result_df)) {
    stop("Stage_full row not found in ", distance_name, " PERMANOVA.")
  }

  tibble(
    Distance = distance_name,
    F = as.numeric(result_df["Stage_full", "F"]),
    R2 = as.numeric(result_df["Stage_full", "R2"]),
    P_value = as.numeric(result_df["Stage_full", "Pr(>F)"]),
    Permutations = N_PERM
  )
}

permanova_results <- bind_rows(
  extract_stage_permanova(bray_permanova, "Bray-Curtis"),
  extract_stage_permanova(jaccard_permanova, "Jaccard")
)

print(permanova_results)

write.csv(
  permanova_results,
  file.path(OUT_TAB, "05_PERMANOVA_Lactation_Stage.csv"),
  row.names = FALSE
)


# =============================================================================
# 24. ANOSIM
# =============================================================================
# Complementary rank-based test of Early-vs-Late community separation.
# Permutations are also restricted within PairUID.

set.seed(12345)

bray_anosim <- vegan::anosim(
  bray_dist,
  grouping = meta_div$Stage_full,
  permutations = N_PERM,
  strata = meta_div$PairUID
)

set.seed(12345)

jaccard_anosim <- vegan::anosim(
  jaccard_dist,
  grouping = meta_div$Stage_full,
  permutations = N_PERM,
  strata = meta_div$PairUID
)

anosim_results <- tibble(
  Distance = c("Bray-Curtis", "Jaccard"),
  ANOSIM_R = c(
    unname(bray_anosim$statistic),
    unname(jaccard_anosim$statistic)
  ),
  P_value = c(
    bray_anosim$signif,
    jaccard_anosim$signif
  ),
  Permutations = N_PERM
)

print(anosim_results)

write.csv(
  anosim_results,
  file.path(OUT_TAB, "06_ANOSIM_Lactation_Stage.csv"),
  row.names = FALSE
)


# =============================================================================
# 25. PERMDISP DIAGNOSTIC
# =============================================================================
# PERMDISP checks whether Early and Late differ in within-group multivariate
# dispersion, which helps interpret PERMANOVA results.

run_permdisp <- function(dist_object, distance_name) {

  dispersion <- vegan::betadisper(
    dist_object,
    group = meta_div$Stage_full,
    type = "median"
  )

  set.seed(12345)

  perm_test <- vegan::permutest(
    dispersion,
    permutations = N_PERM
  )

  tibble(
    Distance = distance_name,
    F = as.numeric(perm_test$tab[1, "F"]),
    P_value = as.numeric(perm_test$tab[1, "Pr(>F)"]),
    Permutations = N_PERM
  )
}

permdisp_results <- bind_rows(
  run_permdisp(bray_dist, "Bray-Curtis"),
  run_permdisp(jaccard_dist, "Jaccard")
)

print(permdisp_results)

write.csv(
  permdisp_results,
  file.path(OUT_TAB, "07_PERMDISP_Diagnostic.csv"),
  row.names = FALSE
)


# =============================================================================
# 26. CREATE FINAL BETA-STATISTICS TABLE
# =============================================================================

beta_final_statistics <- permanova_results %>%
  select(
    Distance,
    PERMANOVA_F = F,
    PERMANOVA_R2 = R2,
    PERMANOVA_P = P_value
  ) %>%
  left_join(
    anosim_results %>%
      select(
        Distance,
        ANOSIM_R,
        ANOSIM_P = P_value
      ),
    by = "Distance"
  ) %>%
  left_join(
    permdisp_results %>%
      select(
        Distance,
        PERMDISP_F = F,
        PERMDISP_P = P_value
      ),
    by = "Distance"
  )

print(beta_final_statistics)

write.csv(
  beta_final_statistics,
  file.path(OUT_TAB, "08_Beta_Diversity_Final_Statistics.csv"),
  row.names = FALSE
)


# =============================================================================
# 27. BETA-DIVERSITY PCoA PLOT FUNCTION
# =============================================================================

plot_pcoa_panel <- function(pcoa_result, R2_value, p_value, panel_tag) {

  if (is.na(p_value)) {
    p_text <- "p = NA"
  } else if (p_value < 0.001) {
    p_text <- "p < 0.001"
  } else {
    p_text <- paste0(
      "p = ",
      format(round(p_value, 3), nsmall = 3)
    )
  }

  if (is.na(R2_value)) {
    R2_text <- "NA"
  } else {
    R2_text <- format(round(R2_value, 3), nsmall = 3)
  }

  stat_label <- paste0(
    "PERMANOVA\n",
    "R\u00B2 = ",
    R2_text,
    "\n",
    p_text
  )

  ggplot(
    pcoa_result$data,
    aes(
      x = PCoA1,
      y = PCoA2,
      colour = Stage_full,
      fill = Stage_full,
      shape = Stage_full
    )
  ) +
    stat_ellipse(
      aes(group = Stage_full),
      geom = "polygon",
      alpha = 0.17,
      linewidth = 1.3,
      show.legend = FALSE
    ) +
    geom_point(size = 6.2, stroke = 1.1) +
    scale_colour_manual(values = STAGE_COLORS, name = "Lactation") +
    scale_fill_manual(values = STAGE_COLORS, guide = "none") +
    scale_shape_manual(
      values = c("Early" = 17, "Late" = 17),
      name = "Lactation"
    ) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = stat_label,
      hjust = 1.08,
      vjust = 1.20,
      size = FONT_STAT,
      fontface = "bold"
    ) +
    labs(
      tag = panel_tag,
      title = paste0(pcoa_result$name, " PCoA (Milk)"),
      x = paste0("PCoA1 (", round(pcoa_result$pc1, 2), "%)"),
      y = paste0("PCoA2 (", round(pcoa_result$pc2, 2), "%)")
    ) +
    theme_bw(base_size = FONT_AXIS_TEXT) +
    theme(
      plot.tag = element_text(size = FONT_PANEL_TAG, face = "bold.italic"),
      plot.title = element_text(size = FONT_TITLE, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = FONT_AXIS_TITLE, face = "bold"),
      axis.title.y = element_text(size = FONT_AXIS_TITLE, face = "bold"),
      axis.text.x = element_text(size = FONT_AXIS_TEXT, face = "bold"),
      axis.text.y = element_text(size = FONT_AXIS_TEXT, face = "bold"),
      legend.title = element_text(size = FONT_LEGEND_TITLE, face = "bold"),
      legend.text = element_text(size = FONT_LEGEND_TEXT, face = "bold"),
      panel.grid.minor = element_blank(),
      aspect.ratio = 1,
      plot.margin = margin(12, 12, 12, 12)
    )
}


# =============================================================================
# 28. CREATE BETA PANELS E-F
# =============================================================================

bray_R2 <- permanova_results %>%
  filter(Distance == "Bray-Curtis") %>%
  pull(R2)

bray_p <- permanova_results %>%
  filter(Distance == "Bray-Curtis") %>%
  pull(P_value)

jaccard_R2 <- permanova_results %>%
  filter(Distance == "Jaccard") %>%
  pull(R2)

jaccard_p <- permanova_results %>%
  filter(Distance == "Jaccard") %>%
  pull(P_value)

p_bray <- plot_pcoa_panel(
  bray_pcoa,
  bray_R2,
  bray_p,
  "E"
)

p_jaccard <- plot_pcoa_panel(
  jaccard_pcoa,
  jaccard_R2,
  jaccard_p,
  "F"
)


# =============================================================================
# 29. COMBINE AND SAVE BETA-DIVERSITY FIGURE
# =============================================================================

p_beta <- (
  p_bray /
    p_jaccard +
    plot_layout(
      ncol = 1,
      guides = "collect"
    )
) &
  theme(legend.position = "right")

print(p_beta)

ggsave(
  file.path(OUT_FIG, "02_Beta_Diversity_E-F_Vertical_Shared_Legend.png"),
  p_beta,
  width = 14,
  height = 22,
  units = "in",
  dpi = 300,
  limitsize = FALSE
)

ggsave(
  file.path(OUT_FIG, "02_Beta_Diversity_E-F_Vertical_Shared_Legend.pdf"),
  p_beta,
  width = 14,
  height = 22,
  units = "in",
  device = cairo_pdf,
  limitsize = FALSE
)


# =============================================================================
# 30. CREATE COMBINED A-F FIGURE
# =============================================================================
# Layout:
#   A | C | E
#   B | D | F

p_all_six <- (
  (p_observed / p_chao1) |
    (p_shannon / p_simpson) |
    (p_bray / p_jaccard)
) +
  plot_layout(
    guides = "collect",
    widths = c(1, 1, 1)
  )

p_all_six <- p_all_six &
  theme(legend.position = "right")

print(p_all_six)

ggsave(
  file.path(OUT_FIG, "03_Combined_Alpha_Beta_A-F.png"),
  p_all_six,
  width = 30,
  height = 17,
  units = "in",
  dpi = 300,
  limitsize = FALSE
)

ggsave(
  file.path(OUT_FIG, "03_Combined_Alpha_Beta_A-F.pdf"),
  p_all_six,
  width = 30,
  height = 17,
  units = "in",
  device = cairo_pdf,
  limitsize = FALSE
)


# =============================================================================
# 31. PRINT FINAL REPORT-READY RESULTS
# =============================================================================

cat("\n\n============================================================\n")
cat("FINAL ALPHA AND BETA DIVERSITY RESULTS\n")
cat("============================================================\n")

cat("\nSamples analysed:", nrow(meta_div), "\n")
cat("Early-lactation samples:", n_early, "\n")
cat("Late-lactation samples:", n_late, "\n")
cat("Complete Early-Late animal pairs:", n_complete_pairs, "\n")

cat("\n------------------------------------------------------------\n")
cat("ALPHA DIVERSITY: PAIRED WILCOXON TESTS\n")
cat("------------------------------------------------------------\n")
print(alpha_stats)

cat("\n------------------------------------------------------------\n")
cat("BETA DIVERSITY: FINAL STATISTICS\n")
cat("------------------------------------------------------------\n")
print(beta_final_statistics)

cat("\n------------------------------------------------------------\n")
cat("PERMANOVA LACTATION-STAGE RESULTS\n")
cat("------------------------------------------------------------\n")

for (i in seq_len(nrow(permanova_results))) {

  x <- permanova_results[i, ]

  cat(
    x$Distance,
    ": F = ",
    round(x$F, 3),
    ", R2 = ",
    round(x$R2, 4),
    ", p = ",
    format.pval(
      x$P_value,
      digits = 3,
      eps = 0.001
    ),
    "\n",
    sep = ""
  )
}

cat("\nFigures saved to:\n", OUT_FIG, "\n")
cat("\nTables saved to:\n", OUT_TAB, "\n")
cat("============================================================\n")


# =============================================================================
# 32. SAVE R SESSION INFORMATION
# =============================================================================
# Records R version and package versions for reproducibility.

capture.output(
  sessionInfo(),
  file = file.path(OUT_TAB, "09_R_SessionInfo.txt")
)

cat("\nFigure 3 alpha/beta diversity analysis completed successfully.\n")
