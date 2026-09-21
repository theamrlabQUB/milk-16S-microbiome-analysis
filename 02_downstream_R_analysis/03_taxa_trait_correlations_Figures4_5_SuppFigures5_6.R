# =============================================================================
# Figures 4–5 and Supplementary Figures 5–6
# Lactation stage-specific correlations between bacterial taxa and
# milk nutritional and production traits
# =============================================================================
#
# Manuscript figures
# ------------------
# Figure 4
#   Lactation stage-specific correlations between bacterial SPECIES and
#   milk nutritional and production traits.
#   Significance annotations are based on Benjamini-Hochberg FDR-adjusted
#   P values.
#
# Figure 5
#   Lactation stage-specific correlations between bacterial GENERA and
#   milk nutritional and production traits.
#   Significance annotations are based on Benjamini-Hochberg FDR-adjusted
#   P values.
#
# Supplementary Figure 5
#   Species-level correlations BEFORE multiple-testing correction.
#   Significance annotations are based on raw (unadjusted) P values.
#
# Supplementary Figure 6
#   Genus-level correlations BEFORE multiple-testing correction.
#   Significance annotations are based on raw (unadjusted) P values.
#
# Statistical approach
# --------------------
# - Spearman rank correlation (two-sided)
# - Top 50 taxa selected after a 20% prevalence filter
# - Taxa must be detected (>0 relative abundance) in at least 20% of samples
#   in each Farm x Lactation subgroup
# - BH-FDR correction is applied separately within each taxonomic-rank /
#   subgroup analysis across all taxon x trait tests in that subgroup
# - Heatmap colour represents Spearman's rho
#
# Figure layout for each taxonomic rank
# -------------------------------------
# Column 1: Farm A
#   A. Early lactation
#   B. Late lactation
#
# Column 2: Farm B
#   C. Early lactation
#   D. Late lactation
#
# Column 3: Combined farms
#   E. Early lactation
#   F. Late lactation
#
# Data availability
# -----------------
# Data are not stored in this GitHub repository. Update the two input paths
# below after downloading/preparing the ENA-derived/processed study data.
# =============================================================================


# =============================================================================
# 0. REQUIRED PACKAGES
# =============================================================================

required_packages <- c(
  "tidyverse",
  "stringr",
  "tibble",
  "ggplot2",
  "patchwork",
  "scales",
  "grid"
)

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
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(grid)
})


# =============================================================================
# 1. INPUT AND OUTPUT PATHS
# =============================================================================

DATA_DIR <- "path/to/downloaded_or_processed_study_data"

MICROBIOME_FILE <- file.path(DATA_DIR, "species_rel_no30.tsv")
METADATA_FILE   <- file.path(DATA_DIR, "metadata_clean.csv")

OUT_DIR <- file.path("results", "taxa_trait_correlations")
OUT_FIG <- file.path(OUT_DIR, "figures")
OUT_TAB <- file.path(OUT_DIR, "tables")

dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(MICROBIOME_FILE)) {
  stop("Microbiome file not found: ", MICROBIOME_FILE)
}

if (!file.exists(METADATA_FILE)) {
  stop("Metadata file not found: ", METADATA_FILE)
}


# =============================================================================
# 2. ANALYSIS SETTINGS
# =============================================================================

TOP_N <- 50
MIN_PREVALENCE <- 0.20

TRAIT_VARS <- c(
  "Production",
  "Fat",
  "Energy",
  "Lactose",
  "Mineral",
  "Protein"
)

# Significance thresholds
P_CUTS <- c(0.001, 0.01, 0.05)

# Figure typography
OVERALL_TITLE_SIZE  <- 30
PANEL_TITLE_SIZE    <- 22
X_LABEL_SIZE        <- 16
TAXON_LABEL_SIZE    <- 16
LEGEND_TITLE_SIZE   <- 16
LEGEND_TEXT_SIZE    <- 14
PANEL_TAG_SIZE      <- 22
STAR_SIZE           <- 5.2
TILE_LINEWIDTH      <- 0.35


# =============================================================================
# 3. READ INPUT DATA
# =============================================================================

micro_raw <- read.delim(
  MICROBIOME_FILE,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (ncol(micro_raw) < 2) {
  stop("Microbiome table must contain a taxonomy column and sample columns.")
}

colnames(micro_raw)[1] <- "Taxa"

meta <- read.csv(
  METADATA_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_meta_cols <- c(
  "SampleID",
  "Sample_Name",
  "Lactation",
  "Farm",
  TRAIT_VARS
)

missing_meta_cols <- setdiff(required_meta_cols, names(meta))

if (length(missing_meta_cols) > 0) {
  stop(
    "Missing metadata columns: ",
    paste(missing_meta_cols, collapse = ", ")
  )
}

meta <- meta %>%
  mutate(
    SampleID    = str_trim(as.character(SampleID)),
    Sample_Name = str_trim(as.character(Sample_Name)),
    Lactation   = str_trim(as.character(Lactation)),
    Farm        = str_trim(as.character(Farm)),
    across(
      all_of(TRAIT_VARS),
      ~ suppressWarnings(as.numeric(.x))
    ),
    PairID = str_remove(Sample_Name, "_(AB|NI)2?$"),
    Stage_full = case_when(
      tolower(Lactation) == "early" ~ "Early",
      tolower(Lactation) == "late"  ~ "Late",
      TRUE ~ NA_character_
    ),
    Farm_label = case_when(
      Farm == "A" ~ "Farm A",
      Farm == "B" ~ "Farm B",
      TRUE ~ NA_character_
    )
  )


# =============================================================================
# 4. MATCH SAMPLES ACROSS DATASETS
# =============================================================================

micro_samples <- colnames(micro_raw)[-1]

meta_use <- meta %>%
  filter(
    SampleID %in% micro_samples,
    Farm %in% c("A", "B"),
    Stage_full %in% c("Early", "Late")
  )

if (nrow(meta_use) == 0) {
  stop("No metadata samples matched the microbiome table.")
}

if (anyDuplicated(meta_use$SampleID)) {
  stop("Duplicate SampleID values detected in metadata.")
}

cat("\nMatched sample counts:\n")
print(
  meta_use %>%
    count(Farm_label, Stage_full, name = "N_samples")
)
cat("Total matched samples:", nrow(meta_use), "\n")


# =============================================================================
# 5. TAXONOMY EXTRACTION
# =============================================================================
# The first column is expected to contain semicolon-separated taxonomy.
# Greengenes/Greengenes2-style prefixes such as g__ and s__ are removed.

extract_tax_rank <- function(taxonomy, rank = c("Genus", "Species")) {

  rank <- match.arg(rank)

  position <- switch(
    rank,
    Genus   = 6,
    Species = 7
  )

  prefix <- switch(
    rank,
    Genus   = "g__",
    Species = "s__"
  )

  output <- vapply(
    taxonomy,
    function(x) {

      parts <- unlist(strsplit(x, ";", fixed = TRUE))

      if (length(parts) < position) {
        return(NA_character_)
      }

      value <- trimws(parts[position])
      value <- sub(paste0("^", prefix), "", value)
      value <- trimws(value)

      if (
        is.na(value) ||
        value == "" ||
        value == "__" ||
        value == "_"
      ) {
        return(NA_character_)
      }

      gsub("_", " ", value)
    },
    character(1)
  )

  unname(output)
}


# =============================================================================
# 6. BUILD RELATIVE-ABUNDANCE MATRIX
# =============================================================================

sample_ids <- meta_use$SampleID

micro_use <- micro_raw %>%
  select(Taxa, all_of(sample_ids))

micro_mat <- as.matrix(
  micro_use[, -1, drop = FALSE]
)

mode(micro_mat) <- "numeric"
rownames(micro_mat) <- micro_use$Taxa

if (anyNA(micro_mat)) {
  stop("NA values were produced when microbiome abundances were converted to numeric.")
}

if (any(micro_mat < 0)) {
  stop("Negative abundances detected in microbiome table.")
}

sample_totals <- colSums(micro_mat)

if (any(sample_totals <= 0)) {
  stop("One or more matched samples have zero total abundance.")
}

# Re-normalise within each sample to relative abundance (%).
micro_percent <- sweep(
  micro_mat,
  2,
  sample_totals,
  "/"
) * 100


# =============================================================================
# 7. PREVALENCE FILTER AND TOP-50 TAXON SELECTION
# =============================================================================

build_top_taxa <- function(rank, top_n = TOP_N) {

  rank_labels <- extract_tax_rank(
    rownames(micro_percent),
    rank
  )

  rank_df <- as.data.frame(
    micro_percent,
    check.names = FALSE
  ) %>%
    mutate(TaxaName = rank_labels) %>%
    filter(!is.na(TaxaName))

  # Sum features assigned to the same named taxon.
  rank_sum <- rank_df %>%
    group_by(TaxaName) %>%
    summarise(
      across(
        all_of(sample_ids),
        ~ sum(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    )

  rank_mat <- rank_sum %>%
    column_to_rownames("TaxaName") %>%
    as.matrix()

  mode(rank_mat) <- "numeric"

  group_samples <- list(
    FarmA_Early = meta_use %>%
      filter(Farm_label == "Farm A", Stage_full == "Early") %>%
      pull(SampleID),

    FarmA_Late = meta_use %>%
      filter(Farm_label == "Farm A", Stage_full == "Late") %>%
      pull(SampleID),

    FarmB_Early = meta_use %>%
      filter(Farm_label == "Farm B", Stage_full == "Early") %>%
      pull(SampleID),

    FarmB_Late = meta_use %>%
      filter(Farm_label == "Farm B", Stage_full == "Late") %>%
      pull(SampleID)
  )

  prevalence_table <- tibble(
    TaxaName = rownames(rank_mat)
  )

  for (g in names(group_samples)) {

    ids <- intersect(
      group_samples[[g]],
      colnames(rank_mat)
    )

    if (length(ids) == 0) {
      stop("No samples available for subgroup: ", g)
    }

    prevalence_table[[g]] <- rowMeans(
      rank_mat[, ids, drop = FALSE] > 0
    )
  }

  prevalence_table <- prevalence_table %>%
    mutate(
      Keep_robust = if_all(
        all_of(
          c(
            "FarmA_Early",
            "FarmA_Late",
            "FarmB_Early",
            "FarmB_Late"
          )
        ),
        ~ .x >= MIN_PREVALENCE
      )
    )

  robust_taxa <- prevalence_table %>%
    filter(Keep_robust) %>%
    pull(TaxaName)

  cat(
    "\n", rank,
    ": taxa before prevalence filter = ", nrow(rank_mat),
    "; retained = ", length(robust_taxa),
    "\n",
    sep = ""
  )

  if (length(robust_taxa) == 0) {
    stop(
      "No ", rank,
      " taxa passed the prevalence filter."
    )
  }

  robust_mat <- rank_mat[
    robust_taxa,
    ,
    drop = FALSE
  ]

  mean_abundance <- rowMeans(
    robust_mat,
    na.rm = TRUE
  )

  n_use <- min(
    top_n,
    length(mean_abundance)
  )

  top_taxa <- names(
    sort(
      mean_abundance,
      decreasing = TRUE
    )
  )[seq_len(n_use)]

  write.csv(
    prevalence_table %>%
      arrange(desc(Keep_robust), TaxaName),
    file.path(
      OUT_TAB,
      paste0("QC_", rank, "_Prevalence_Filter.csv")
    ),
    row.names = FALSE
  )

  list(
    rank = rank,
    top_taxa = top_taxa,
    matrix = robust_mat[
      top_taxa,
      ,
      drop = FALSE
    ],
    mean_abundance = mean_abundance[top_taxa]
  )
}

species_top50 <- build_top_taxa(
  "Species",
  TOP_N
)

genus_top50 <- build_top_taxa(
  "Genus",
  TOP_N
)

write.csv(
  tibble(
    Rank_order = seq_along(species_top50$top_taxa),
    Species = species_top50$top_taxa,
    Mean_relative_abundance_percent =
      species_top50$mean_abundance
  ),
  file.path(
    OUT_TAB,
    "Top50_Species_List.csv"
  ),
  row.names = FALSE
)

write.csv(
  tibble(
    Rank_order = seq_along(genus_top50$top_taxa),
    Genus = genus_top50$top_taxa,
    Mean_relative_abundance_percent =
      genus_top50$mean_abundance
  ),
  file.path(
    OUT_TAB,
    "Top50_Genera_List.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 8. CREATE ANALYSIS TABLES
# =============================================================================

make_analysis_df <- function(taxa_result) {

  taxa_df <- as.data.frame(
    t(taxa_result$matrix),
    check.names = FALSE
  ) %>%
    rownames_to_column("SampleID")

  meta_use %>%
    select(
      SampleID,
      Farm,
      Farm_label,
      Stage_full,
      all_of(TRAIT_VARS)
    ) %>%
    inner_join(
      taxa_df,
      by = "SampleID"
    )
}

species_analysis_df <- make_analysis_df(
  species_top50
)

genus_analysis_df <- make_analysis_df(
  genus_top50
)


# =============================================================================
# 9. SAFE SPEARMAN CORRELATION
# =============================================================================

safe_spearman <- function(x, y) {

  keep <- complete.cases(x, y)

  x <- x[keep]
  y <- y[keep]

  n_complete <- length(x)

  if (
    n_complete < 4 ||
    length(unique(x)) < 2 ||
    length(unique(y)) < 2
  ) {
    return(
      tibble(
        N = n_complete,
        Rho = NA_real_,
        P_value = NA_real_
      )
    )
  }

  test <- suppressWarnings(
    cor.test(
      x,
      y,
      method = "spearman",
      exact = FALSE,
      alternative = "two.sided"
    )
  )

  tibble(
    N = n_complete,
    Rho = unname(test$estimate),
    P_value = test$p.value
  )
}


# =============================================================================
# 10. CALCULATE CORRELATIONS AND BH-FDR
# =============================================================================
# BH adjustment is performed within each subgroup/rank analysis across all
# Top-50 taxon x trait tests in that subgroup.

calculate_correlations <- function(
    data,
    taxa_names,
    analysis_name,
    rank_name
) {

  results <- vector(
    "list",
    length(taxa_names) * length(TRAIT_VARS)
  )

  k <- 1

  for (taxon in taxa_names) {

    for (trait in TRAIT_VARS) {

      test_result <- safe_spearman(
        data[[taxon]],
        data[[trait]]
      )

      results[[k]] <- tibble(
        Analysis = analysis_name,
        Rank = rank_name,
        TaxaName = taxon,
        Trait = trait,
        N = test_result$N,
        Rho = test_result$Rho,
        P_value = test_result$P_value
      )

      k <- k + 1
    }
  }

  results <- bind_rows(results)

  results$FDR <- NA_real_

  valid <- which(
    !is.na(results$P_value)
  )

  if (length(valid) > 0) {
    results$FDR[valid] <- p.adjust(
      results$P_value[valid],
      method = "BH"
    )
  }

  results %>%
    mutate(
      Raw_significance = case_when(
        is.na(P_value) ~ "",
        P_value < P_CUTS[1] ~ "***",
        P_value < P_CUTS[2] ~ "**",
        P_value < P_CUTS[3] ~ "*",
        TRUE ~ ""
      ),
      FDR_significance = case_when(
        is.na(FDR) ~ "",
        FDR < P_CUTS[1] ~ "***",
        FDR < P_CUTS[2] ~ "**",
        FDR < P_CUTS[3] ~ "*",
        TRUE ~ ""
      ),
      Direction = case_when(
        is.na(Rho) ~ NA_character_,
        Rho > 0 ~ "Positive",
        Rho < 0 ~ "Negative",
        TRUE ~ "Zero"
      )
    )
}


# =============================================================================
# 11. RUN SIX SUBGROUP ANALYSES FOR EACH TAXONOMIC RANK
# =============================================================================

run_rank_analyses <- function(
    analysis_df,
    taxa_names,
    rank_name
) {

  list(
    farmA_early = calculate_correlations(
      analysis_df %>%
        filter(
          Farm_label == "Farm A",
          Stage_full == "Early"
        ),
      taxa_names,
      "Farm A - Early lactation",
      rank_name
    ),

    farmA_late = calculate_correlations(
      analysis_df %>%
        filter(
          Farm_label == "Farm A",
          Stage_full == "Late"
        ),
      taxa_names,
      "Farm A - Late lactation",
      rank_name
    ),

    farmB_early = calculate_correlations(
      analysis_df %>%
        filter(
          Farm_label == "Farm B",
          Stage_full == "Early"
        ),
      taxa_names,
      "Farm B - Early lactation",
      rank_name
    ),

    farmB_late = calculate_correlations(
      analysis_df %>%
        filter(
          Farm_label == "Farm B",
          Stage_full == "Late"
        ),
      taxa_names,
      "Farm B - Late lactation",
      rank_name
    ),

    overall_early = calculate_correlations(
      analysis_df %>%
        filter(
          Stage_full == "Early"
        ),
      taxa_names,
      "Combined - Early lactation",
      rank_name
    ),

    overall_late = calculate_correlations(
      analysis_df %>%
        filter(
          Stage_full == "Late"
        ),
      taxa_names,
      "Combined - Late lactation",
      rank_name
    )
  )
}

species_cor <- run_rank_analyses(
  species_analysis_df,
  species_top50$top_taxa,
  "Species"
)

genus_cor <- run_rank_analyses(
  genus_analysis_df,
  genus_top50$top_taxa,
  "Genus"
)


# =============================================================================
# 12. SAVE COMPLETE CORRELATION TABLES
# =============================================================================
# Each table contains both raw P values and BH-FDR-adjusted P values, allowing
# the main and supplementary figures to be reproduced from the same results.

save_result_tables <- function(
    result_list,
    prefix
) {

  for (nm in names(result_list)) {

    d <- result_list[[nm]]

    write.csv(
      d,
      file.path(
        OUT_TAB,
        paste0(
          prefix,
          "_",
          nm,
          "_Spearman_rawP_and_FDR.csv"
        )
      ),
      row.names = FALSE
    )

    write.csv(
      d %>%
        filter(
          !is.na(P_value),
          P_value < 0.05
        ) %>%
        arrange(
          P_value,
          desc(abs(Rho))
        ),
      file.path(
        OUT_TAB,
        paste0(
          prefix,
          "_",
          nm,
          "_Significant_rawP05.csv"
        )
      ),
      row.names = FALSE
    )

    write.csv(
      d %>%
        filter(
          !is.na(FDR),
          FDR < 0.05
        ) %>%
        arrange(
          FDR,
          desc(abs(Rho))
        ),
      file.path(
        OUT_TAB,
        paste0(
          prefix,
          "_",
          nm,
          "_Significant_FDR05.csv"
        )
      ),
      row.names = FALSE
    )
  }
}

save_result_tables(
  species_cor,
  "Species_Top50"
)

save_result_tables(
  genus_cor,
  "Genera_Top50"
)


# =============================================================================
# 13. COMMON TAXON ORDER
# =============================================================================
# The same row order is used in the raw-P and FDR versions so that the main
# and supplementary figures can be compared directly.

build_common_order <- function(
    result_list,
    taxa_names
) {

  all_results <- bind_rows(
    result_list
  )

  profile <- all_results %>%
    select(
      Analysis,
      TaxaName,
      Trait,
      Rho
    ) %>%
    mutate(
      Key = paste(
        Analysis,
        Trait,
        sep = "__"
      )
    ) %>%
    select(
      TaxaName,
      Key,
      Rho
    ) %>%
    pivot_wider(
      names_from = Key,
      values_from = Rho,
      values_fill = 0
    )

  mat <- profile %>%
    column_to_rownames(
      "TaxaName"
    ) %>%
    as.matrix()

  mat[is.na(mat)] <- 0

  if (nrow(mat) > 1) {
    hc <- hclust(
      dist(mat),
      method = "complete"
    )
    taxon_order <- rownames(mat)[hc$order]
  } else {
    taxon_order <- rownames(mat)
  }

  missing_taxa <- setdiff(
    taxa_names,
    taxon_order
  )

  rev(
    c(
      taxon_order,
      missing_taxa
    )
  )
}

species_order <- build_common_order(
  species_cor,
  species_top50$top_taxa
)

genus_order <- build_common_order(
  genus_cor,
  genus_top50$top_taxa
)


# =============================================================================
# 14. HEATMAP FUNCTION
# =============================================================================

plot_heatmap <- function(
    correlation_data,
    title_text,
    taxa_order,
    rank_name,
    panel_tag,
    significance_mode = c("FDR", "raw")
) {

  significance_mode <- match.arg(
    significance_mode
  )

  d <- correlation_data %>%
    mutate(
      TaxaName = factor(
        TaxaName,
        levels = taxa_order
      ),
      Trait = factor(
        Trait,
        levels = TRAIT_VARS
      ),
      Plot_significance = ifelse(
        significance_mode == "FDR",
        FDR_significance,
        Raw_significance
      )
    )

  ggplot(
    d,
    aes(
      x = Trait,
      y = TaxaName,
      fill = Rho
    )
  ) +
    geom_tile(
      colour = "white",
      linewidth = TILE_LINEWIDTH
    ) +
    geom_text(
      aes(
        label = Plot_significance
      ),
      size = STAR_SIZE,
      fontface = "bold",
      na.rm = TRUE
    ) +
    scale_fill_gradient2(
      low = "#B2182B",
      mid = "white",
      high = "#2166AC",
      midpoint = 0,
      limits = c(-1, 1),
      breaks = c(
        -1,
        -0.5,
        0,
        0.5,
        1
      ),
      name = "Spearman\n\u03C1"
    ) +
    labs(
      tag = panel_tag,
      title = title_text,
      x = NULL,
      y = NULL
    ) +
    theme_minimal(
      base_size = 16
    ) +
    theme(
      plot.tag = element_text(
        size = PANEL_TAG_SIZE,
        face = "bold"
      ),
      plot.title = element_text(
        size = PANEL_TITLE_SIZE,
        face = "bold",
        hjust = 0.5
      ),
      axis.text.x = element_text(
        size = X_LABEL_SIZE,
        face = "bold",
        angle = 45,
        hjust = 1
      ),
      axis.text.y = element_text(
        size = TAXON_LABEL_SIZE,
        face = ifelse(
          rank_name == "Species",
          "italic",
          "plain"
        )
      ),
      legend.title = element_text(
        size = LEGEND_TITLE_SIZE,
        face = "bold"
      ),
      legend.text = element_text(
        size = LEGEND_TEXT_SIZE
      ),
      legend.key.height = unit(
        1.0,
        "cm"
      ),
      legend.key.width = unit(
        0.8,
        "cm"
      ),
      panel.grid = element_blank(),
      plot.margin = margin(
        12,
        12,
        12,
        12
      )
    )
}


# =============================================================================
# 15. SIX-PANEL FIGURE FUNCTION
# =============================================================================

make_six_panel_figure <- function(
    result_list,
    taxa_order,
    rank_name,
    significance_mode = c("FDR", "raw"),
    overall_title,
    file_stem
) {

  significance_mode <- match.arg(
    significance_mode
  )

  p_A <- plot_heatmap(
    result_list$farmA_early,
    "Farm A\nEarly lactation",
    taxa_order,
    rank_name,
    "A",
    significance_mode
  )

  p_B <- plot_heatmap(
    result_list$farmA_late,
    "Farm A\nLate lactation",
    taxa_order,
    rank_name,
    "B",
    significance_mode
  )

  p_C <- plot_heatmap(
    result_list$farmB_early,
    "Farm B\nEarly lactation",
    taxa_order,
    rank_name,
    "C",
    significance_mode
  )

  p_D <- plot_heatmap(
    result_list$farmB_late,
    "Farm B\nLate lactation",
    taxa_order,
    rank_name,
    "D",
    significance_mode
  )

  p_E <- plot_heatmap(
    result_list$overall_early,
    "Combined farms\nEarly lactation",
    taxa_order,
    rank_name,
    "E",
    significance_mode
  )

  p_F <- plot_heatmap(
    result_list$overall_late,
    "Combined farms\nLate lactation",
    taxa_order,
    rank_name,
    "F",
    significance_mode
  )

  combined_plot <- (
    (p_A / p_B) |
      (p_C / p_D) |
      (p_E / p_F)
  ) +
    plot_layout(
      guides = "collect",
      widths = c(1, 1, 1)
    ) +
    plot_annotation(
      title = overall_title,
      theme = theme(
        plot.title = element_text(
          size = OVERALL_TITLE_SIZE,
          face = "bold",
          hjust = 0.5
        )
      )
    ) &
    theme(
      legend.position = "right"
    )

  print(
    combined_plot
  )

  ggsave(
    file.path(
      OUT_FIG,
      paste0(
        file_stem,
        ".png"
      )
    ),
    combined_plot,
    width = 36,
    height = 30,
    units = "in",
    dpi = 300,
    limitsize = FALSE
  )

  ggsave(
    file.path(
      OUT_FIG,
      paste0(
        file_stem,
        ".pdf"
      )
    ),
    combined_plot,
    width = 36,
    height = 30,
    units = "in",
    device = cairo_pdf,
    limitsize = FALSE
  )

  invisible(
    combined_plot
  )
}


# =============================================================================
# 16. FIGURE 4 — SPECIES, AFTER BH-FDR CORRECTION
# =============================================================================

figure4_species <- make_six_panel_figure(
  result_list = species_cor,
  taxa_order = species_order,
  rank_name = "Species",
  significance_mode = "FDR",
  overall_title = paste(
    "Figure 4.",
    "Lactation stage-specific correlations between bacterial species",
    "and milk nutritional and production traits"
  ),
  file_stem = "Figure4_Species_Correlations_FDR"
)


# =============================================================================
# 17. FIGURE 5 — GENERA, AFTER BH-FDR CORRECTION
# =============================================================================

figure5_genera <- make_six_panel_figure(
  result_list = genus_cor,
  taxa_order = genus_order,
  rank_name = "Genus",
  significance_mode = "FDR",
  overall_title = paste(
    "Figure 5.",
    "Lactation stage-specific correlations between bacterial genera",
    "and milk nutritional and production traits"
  ),
  file_stem = "Figure5_Genera_Correlations_FDR"
)


# =============================================================================
# 18. SUPPLEMENTARY FIGURE 5 — SPECIES, BEFORE FDR CORRECTION
# =============================================================================
# Stars represent raw, unadjusted P values:
# * P < 0.05, ** P < 0.01, *** P < 0.001

supp_figure5_species <- make_six_panel_figure(
  result_list = species_cor,
  taxa_order = species_order,
  rank_name = "Species",
  significance_mode = "raw",
  overall_title = paste(
    "Supplementary Figure 5.",
    "Lactation stage-specific correlations between bacterial species",
    "and milk production and nutritional traits"
  ),
  file_stem = "SupplementaryFigure5_Species_Correlations_RawP"
)


# =============================================================================
# 19. SUPPLEMENTARY FIGURE 6 — GENERA, BEFORE FDR CORRECTION
# =============================================================================
# Stars represent raw, unadjusted P values:
# * P < 0.05, ** P < 0.01, *** P < 0.001

supp_figure6_genera <- make_six_panel_figure(
  result_list = genus_cor,
  taxa_order = genus_order,
  rank_name = "Genus",
  significance_mode = "raw",
  overall_title = paste(
    "Supplementary Figure 6.",
    "Lactation stage-specific correlations between bacterial genera",
    "and milk production and nutritional traits"
  ),
  file_stem = "SupplementaryFigure6_Genera_Correlations_RawP"
)


# =============================================================================
# 20. SUMMARY TABLES
# =============================================================================

count_raw_significant <- function(x) {
  sum(
    !is.na(x$P_value) &
      x$P_value < 0.05
  )
}

count_fdr_significant <- function(x) {
  sum(
    !is.na(x$FDR) &
      x$FDR < 0.05
  )
}

summary_table <- tibble(
  Rank = c(
    rep("Species", 6),
    rep("Genus", 6)
  ),
  Analysis = c(
    names(species_cor),
    names(genus_cor)
  ),
  Significant_rawP05 = c(
    vapply(
      species_cor,
      count_raw_significant,
      numeric(1)
    ),
    vapply(
      genus_cor,
      count_raw_significant,
      numeric(1)
    )
  ),
  Significant_FDR05 = c(
    vapply(
      species_cor,
      count_fdr_significant,
      numeric(1)
    ),
    vapply(
      genus_cor,
      count_fdr_significant,
      numeric(1)
    )
  )
)

print(
  summary_table
)

write.csv(
  summary_table,
  file.path(
    OUT_TAB,
    "Summary_Significant_Correlations_RawP_vs_FDR.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 21. SAVE SESSION INFORMATION
# =============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    OUT_TAB,
    "R_SessionInfo.txt"
  )
)


# =============================================================================
# 22. FINAL MESSAGE
# =============================================================================

cat("\n============================================================\n")
cat("CORRELATION ANALYSIS COMPLETED\n")
cat("============================================================\n")
cat("\nMain manuscript figures (BH-FDR significance):\n")
cat("1. Figure4_Species_Correlations_FDR.png / .pdf\n")
cat("2. Figure5_Genera_Correlations_FDR.png / .pdf\n")
cat("\nSupplementary figures (raw P-value significance; before FDR):\n")
cat("3. SupplementaryFigure5_Species_Correlations_RawP.png / .pdf\n")
cat("4. SupplementaryFigure6_Genera_Correlations_RawP.png / .pdf\n")
cat("\nFigures saved to:\n", OUT_FIG, "\n")
cat("\nTables saved to:\n", OUT_TAB, "\n")
cat("============================================================\n")
