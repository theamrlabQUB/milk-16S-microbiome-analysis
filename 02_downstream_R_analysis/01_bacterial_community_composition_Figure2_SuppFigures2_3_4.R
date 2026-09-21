# =============================================================================
# Figure 2 and Supplementary Figures 2–4: milk microbiome composition
# =============================================================================
#
# Purpose
# -------
# Reproduce Figure 2 of the manuscript:
#
# (a) Mean relative abundance of the top 20 species in Early- and
#     Late-lactation milk samples in the combined dataset.
# (b) Species-level composition stratified by farm and lactation stage.
# (c) Phylum-level composition stratified by farm and lactation stage.
# (d) Genus-level composition stratified by farm and lactation stage.
#
# Bars represent mean relative abundance.
#
# Data availability
# -----------------
# Sequence data are deposited in ENA and are not stored in this GitHub
# repository. Update the input paths below after downloading/preparing the
# required files.
#
# Required metadata columns
# -------------------------
# Sample_Id   : unique sample identifier
# SampleType  : sample type; milk samples must be labelled "Milk"
# Farm        : farm identifier (for example, A and B)
# Lactation   : lactation stage; expected values are Early and Late
#
# Required QIIME 2 files
# ----------------------
# - Feature table after mitochondrial/chloroplast filtering
# - Greengenes2 taxonomy assignment
#
# Notes
# -----
# 1. Relative abundances are calculated within each sample.
# 2. Taxa outside the selected top 20 are grouped as "Other", so each
#    stacked bar represents the full community (100%).
# 3. The same top-20 species set is used in panels A and B.
# 4. Top-20 phyla and genera are selected across all milk samples.
# =============================================================================


# =============================================================================
# 0. REQUIRED PACKAGES
# =============================================================================

required_packages <- c(
  "readxl",
  "phyloseq",
  "qiime2R",
  "ggplot2",
  "dplyr",
  "stringr",
  "patchwork"
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
  library(readxl)
  library(phyloseq)
  library(qiime2R)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(patchwork)
})


# =============================================================================
# 1. INPUT AND OUTPUT PATHS
# =============================================================================
# Update DATA_DIR to the location of the downloaded/processed study data.

DATA_DIR <- "path/to/downloaded/study_data"

FEATURE_QZA <- file.path(DATA_DIR, "table_noMitoChloro.qza")
TAXONOMY_QZA <- file.path(DATA_DIR, "taxonomy_gg2_full.qza")

# Use either metadata_all.xlsx or metadata_all.csv in DATA_DIR.
METADATA_XLSX <- file.path(DATA_DIR, "metadata_all.xlsx")
METADATA_CSV  <- file.path(DATA_DIR, "metadata_all.csv")

OUT_DIR <- file.path("results", "figure2_community_composition")
OUT_FIG <- file.path(OUT_DIR, "figures")
OUT_TAB <- file.path(OUT_DIR, "tables")

dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(FEATURE_QZA)) {
  stop("Feature table not found: ", FEATURE_QZA)
}

if (!file.exists(TAXONOMY_QZA)) {
  stop("Taxonomy file not found: ", TAXONOMY_QZA)
}


# =============================================================================
# 2. PLOT SETTINGS
# =============================================================================

TOP_N <- 20

FONT_SIZE_AXIS_TEXT    <- 14
FONT_SIZE_AXIS_TITLE   <- 16
FONT_SIZE_PLOT_TITLE   <- 17
FONT_SIZE_STRIP        <- 15
FONT_SIZE_LEGEND_TEXT  <- 10
FONT_SIZE_LEGEND_TITLE <- 12
FONT_SIZE_PANEL_TAG    <- 18

MASTER_COLORS <- c(
  "#E63946", "#F1A7A1", "#F39C12", "#F5D06D", "#E67E22", "#00C875",
  "#3CB371", "#A8DADC", "#3498DB", "#0000CD", "#00CED1", "#457B9D",
  "#6A5ACD", "#8E44AD", "#FF69B4", "#16A085", "#DDA0DD", "#7B3F00",
  "#2F4F4F", "#FFD700", "#8B0000", "#FF8C00", "#556B2F", "#20B2AA",
  "#1E90FF", "#483D8B", "#9932CC", "#C71585", "#708090", "#228B22"
)

OTHER_COLOR <- "#BDBDBD"

make_top_palette <- function(top_taxa_names) {
  n <- length(top_taxa_names)

  if (n > length(MASTER_COLORS)) {
    stop("Not enough colours available for selected taxa.")
  }

  pal <- MASTER_COLORS[seq_len(n)]
  names(pal) <- top_taxa_names

  c(pal, Other = OTHER_COLOR)
}

publication_theme <- function() {
  theme_bw(base_size = FONT_SIZE_AXIS_TEXT) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = FONT_SIZE_AXIS_TEXT
      ),
      axis.text.y = element_text(size = FONT_SIZE_AXIS_TEXT),
      axis.title = element_text(
        size = FONT_SIZE_AXIS_TITLE,
        face = "bold"
      ),
      plot.title = element_text(
        size = FONT_SIZE_PLOT_TITLE,
        face = "bold",
        hjust = 0.5
      ),
      strip.text = element_text(
        size = FONT_SIZE_STRIP,
        face = "bold"
      ),
      legend.text = element_text(size = FONT_SIZE_LEGEND_TEXT),
      legend.title = element_text(
        size = FONT_SIZE_LEGEND_TITLE,
        face = "bold"
      ),
      plot.tag = element_text(
        size = FONT_SIZE_PANEL_TAG,
        face = "bold"
      ),
      panel.grid.minor = element_blank()
    )
}


# =============================================================================
# 3. LOAD METADATA
# =============================================================================

if (file.exists(METADATA_XLSX)) {
  meta_df <- as.data.frame(read_excel(METADATA_XLSX))
} else if (file.exists(METADATA_CSV)) {
  meta_df <- read.csv(
    METADATA_CSV,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
} else {
  stop(
    "Metadata file not found. Expected metadata_all.xlsx or metadata_all.csv ",
    "inside DATA_DIR."
  )
}

required_meta_cols <- c("Sample_Id", "SampleType", "Farm", "Lactation")
missing_meta_cols <- setdiff(required_meta_cols, colnames(meta_df))

if (length(missing_meta_cols) > 0) {
  stop(
    "Missing required metadata columns: ",
    paste(missing_meta_cols, collapse = ", ")
  )
}

meta_df$Sample_Id  <- str_trim(as.character(meta_df$Sample_Id))
meta_df$SampleType <- str_trim(as.character(meta_df$SampleType))
meta_df$Farm       <- str_trim(as.character(meta_df$Farm))
meta_df$Lactation  <- str_trim(as.character(meta_df$Lactation))

meta_df$Lactation <- case_when(
  tolower(meta_df$Lactation) == "early" ~ "Early",
  tolower(meta_df$Lactation) == "late"  ~ "Late",
  TRUE ~ NA_character_
)

meta_df$Lactation <- factor(
  meta_df$Lactation,
  levels = c("Early", "Late")
)

rownames(meta_df) <- meta_df$Sample_Id


# =============================================================================
# 4. IMPORT QIIME 2 FEATURE TABLE AND TAXONOMY
# =============================================================================

ps <- qza_to_phyloseq(
  features = FEATURE_QZA,
  taxonomy = TAXONOMY_QZA
)


# =============================================================================
# 5. MATCH METADATA AND RETAIN MILK SAMPLES
# =============================================================================

common_samples <- intersect(sample_names(ps), rownames(meta_df))

if (length(common_samples) == 0) {
  stop("No sample IDs are shared between the feature table and metadata.")
}

ps <- prune_samples(common_samples, ps)
meta_matched <- meta_df[common_samples, , drop = FALSE]
sample_data(ps) <- sample_data(meta_matched)

ps_milk <- subset_samples(
  ps,
  SampleType == "Milk" & !is.na(Lactation)
)

ps_milk <- prune_taxa(taxa_sums(ps_milk) > 0, ps_milk)

cat("Milk samples retained:", nsamples(ps_milk), "\n")
print(table(sample_data(ps_milk)$Lactation, useNA = "ifany"))
print(table(sample_data(ps_milk)$Farm, useNA = "ifany"))
print(table(
  Farm = sample_data(ps_milk)$Farm,
  Lactation = sample_data(ps_milk)$Lactation
))

if (nsamples(ps_milk) == 0) {
  stop("No milk samples remained after filtering.")
}


# =============================================================================
# 6. CONVERT COUNTS TO WITHIN-SAMPLE RELATIVE ABUNDANCE
# =============================================================================

if (any(sample_sums(ps_milk) <= 0)) {
  stop("At least one retained milk sample has zero total reads.")
}

ps_rel <- transform_sample_counts(
  ps_milk,
  function(x) x / sum(x)
)


# =============================================================================
# 7. HELPER FUNCTION: PREPARE TOP-N TAXA
# =============================================================================

prepare_top_taxa <- function(ps_object, taxrank, top_n = TOP_N) {

  if (!taxrank %in% rank_names(ps_object)) {
    stop("Taxonomic rank not found in taxonomy table: ", taxrank)
  }

  ps_glom <- tax_glom(
    ps_object,
    taxrank = taxrank,
    NArm = FALSE
  )

  df <- psmelt(ps_glom)

  df[[taxrank]] <- as.character(df[[taxrank]])
  df[[taxrank]][
    is.na(df[[taxrank]]) | df[[taxrank]] == ""
  ] <- "Unclassified"

  top_taxa <- df %>%
    group_by(.data[[taxrank]]) %>%
    summarise(
      MeanAbundance = mean(Abundance),
      .groups = "drop"
    ) %>%
    arrange(desc(MeanAbundance)) %>%
    slice_head(n = top_n) %>%
    pull(.data[[taxrank]])

  df_top <- df %>%
    mutate(
      Taxon = ifelse(
        .data[[taxrank]] %in% top_taxa,
        .data[[taxrank]],
        "Other"
      )
    ) %>%
    group_by(
      Sample,
      Farm,
      Lactation,
      Taxon
    ) %>%
    summarise(
      Abundance = sum(Abundance),
      .groups = "drop"
    )

  df_top$Taxon <- factor(
    df_top$Taxon,
    levels = c(top_taxa, "Other")
  )

  list(
    data = df_top,
    top_taxa = top_taxa,
    palette = make_top_palette(top_taxa)
  )
}


# =============================================================================
# 8. SPECIES LABEL HANDLING
# =============================================================================
# Greengenes2 may store a full binomial in the Species field. The helper below
# avoids creating duplicated labels such as
# "Staphylococcus Staphylococcus equorum".

prepare_top_species <- function(ps_object, top_n = TOP_N) {

  if (!all(c("Genus", "Species") %in% rank_names(ps_object))) {
    stop("Both Genus and Species ranks are required for species-level plots.")
  }

  ps_species <- tax_glom(
    ps_object,
    taxrank = "Species",
    NArm = FALSE
  )

  df <- psmelt(ps_species)

  df$Genus <- as.character(df$Genus)
  df$Species <- as.character(df$Species)

  df$Genus[
    is.na(df$Genus) | df$Genus == ""
  ] <- "Unclassified"

  df$Species[
    is.na(df$Species) | df$Species == ""
  ] <- "Unresolved"

  df$Taxon <- ifelse(
    df$Species == "Unresolved",
    paste(df$Genus, "sp."),
    ifelse(
      startsWith(df$Species, df$Genus),
      df$Species,
      paste(df$Genus, df$Species)
    )
  )

  top_species <- df %>%
    group_by(Taxon) %>%
    summarise(
      MeanAbundance = mean(Abundance),
      .groups = "drop"
    ) %>%
    arrange(desc(MeanAbundance)) %>%
    slice_head(n = top_n) %>%
    pull(Taxon)

  df_top <- df %>%
    mutate(
      Taxon = ifelse(
        Taxon %in% top_species,
        Taxon,
        "Other"
      )
    ) %>%
    group_by(
      Sample,
      Farm,
      Lactation,
      Taxon
    ) %>%
    summarise(
      Abundance = sum(Abundance),
      .groups = "drop"
    )

  df_top$Taxon <- factor(
    df_top$Taxon,
    levels = c(top_species, "Other")
  )

  list(
    data = df_top,
    top_taxa = top_species,
    palette = make_top_palette(top_species)
  )
}


# =============================================================================
# 9. PREPARE SPECIES, PHYLUM AND GENUS DATA
# =============================================================================

species_top <- prepare_top_species(ps_rel, TOP_N)
phylum_top  <- prepare_top_taxa(ps_rel, "Phylum", TOP_N)
genus_top   <- prepare_top_taxa(ps_rel, "Genus", TOP_N)


# =============================================================================
# 10. PANEL A: TOP 20 SPECIES, COMBINED DATASET
# =============================================================================

p_a <- ggplot(
  species_top$data,
  aes(
    x = Lactation,
    y = Abundance * 100,
    fill = Taxon
  )
) +
  geom_bar(
    stat = "summary",
    fun = "mean",
    position = "stack",
    width = 0.75
  ) +
  scale_fill_manual(
    values = species_top$palette,
    name = "Species"
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = c(0, 0)
  ) +
  labs(
    tag = "A",
    title = "Top 20 species",
    x = "Lactation stage",
    y = "Mean relative abundance (%)"
  ) +
  publication_theme()


# =============================================================================
# 11. PANEL B: SPECIES COMPOSITION BY FARM AND LACTATION
# =============================================================================

p_b <- ggplot(
  species_top$data,
  aes(
    x = Lactation,
    y = Abundance * 100,
    fill = Taxon
  )
) +
  geom_bar(
    stat = "summary",
    fun = "mean",
    position = "stack",
    width = 0.75
  ) +
  facet_wrap(~ Farm, nrow = 1) +
  scale_fill_manual(
    values = species_top$palette,
    name = "Species"
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = c(0, 0)
  ) +
  labs(
    tag = "B",
    title = "Species-level composition by farm",
    x = "Lactation stage",
    y = "Mean relative abundance (%)"
  ) +
  publication_theme()


# =============================================================================
# 12. PANEL C: PHYLUM COMPOSITION BY FARM AND LACTATION
# =============================================================================

p_c <- ggplot(
  phylum_top$data,
  aes(
    x = Lactation,
    y = Abundance * 100,
    fill = Taxon
  )
) +
  geom_bar(
    stat = "summary",
    fun = "mean",
    position = "stack",
    width = 0.75
  ) +
  facet_wrap(~ Farm, nrow = 1) +
  scale_fill_manual(
    values = phylum_top$palette,
    name = "Phylum"
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = c(0, 0)
  ) +
  labs(
    tag = "C",
    title = "Phylum-level composition by farm",
    x = "Lactation stage",
    y = "Mean relative abundance (%)"
  ) +
  publication_theme()


# =============================================================================
# 13. PANEL D: GENUS COMPOSITION BY FARM AND LACTATION
# =============================================================================

p_d <- ggplot(
  genus_top$data,
  aes(
    x = Lactation,
    y = Abundance * 100,
    fill = Taxon
  )
) +
  geom_bar(
    stat = "summary",
    fun = "mean",
    position = "stack",
    width = 0.75
  ) +
  facet_wrap(~ Farm, nrow = 1) +
  scale_fill_manual(
    values = genus_top$palette,
    name = "Genus"
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = c(0, 0)
  ) +
  labs(
    tag = "D",
    title = "Genus-level composition by farm",
    x = "Lactation stage",
    y = "Mean relative abundance (%)"
  ) +
  publication_theme()


# =============================================================================
# 14. COMBINE FIGURE 2 PANELS
# =============================================================================
# Separate legends are retained because panels represent different
# taxonomic ranks.

figure2 <- (
  p_a / p_b / p_c / p_d
) +
  plot_layout(
    heights = c(1, 1, 1, 1)
  )

print(figure2)


# =============================================================================
# 15. SAVE FIGURE
# =============================================================================

ggsave(
  filename = file.path(
    OUT_FIG,
    "Figure2_Bacterial_Community_Composition.png"
  ),
  plot = figure2,
  width = 18,
  height = 28,
  units = "in",
  dpi = 300,
  limitsize = FALSE
)

ggsave(
  filename = file.path(
    OUT_FIG,
    "Figure2_Bacterial_Community_Composition.pdf"
  ),
  plot = figure2,
  width = 18,
  height = 28,
  units = "in",
  device = cairo_pdf,
  limitsize = FALSE
)


# =============================================================================
# 16. SAVE TOP-TAXA TABLES
# =============================================================================

write.csv(
  data.frame(
    Rank = "Species",
    Rank_order = seq_along(species_top$top_taxa),
    Taxon = species_top$top_taxa
  ),
  file.path(OUT_TAB, "Figure2_Top20_Species.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    Rank = "Phylum",
    Rank_order = seq_along(phylum_top$top_taxa),
    Taxon = phylum_top$top_taxa
  ),
  file.path(OUT_TAB, "Figure2_Top20_Phyla.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    Rank = "Genus",
    Rank_order = seq_along(genus_top$top_taxa),
    Taxon = genus_top$top_taxa
  ),
  file.path(OUT_TAB, "Figure2_Top20_Genera.csv"),
  row.names = FALSE
)


# =============================================================================
# 17. SAVE SESSION INFORMATION
# =============================================================================

capture.output(
  sessionInfo(),
  file = file.path(OUT_TAB, "Figure2_R_SessionInfo.txt")
)

cat("\nFigure 2 community-composition analysis completed successfully.\n")
cat("Figures saved to:", OUT_FIG, "\n")
cat("Tables saved to:", OUT_TAB, "\n")


# =============================================================================
# SUPPLEMENTARY FIGURE 2
# Identified taxa from raw bovine milk samples (n = 42)
# =============================================================================
#
# Categories are mutually exclusive and represent the deepest taxonomic
# resolution reached for each feature/ASV retained in the milk dataset.
#
# Definitions:
# - Unclassified        : no phylum-level assignment
# - Phylum              : phylum assigned, class unresolved
# - Class               : class assigned, order unresolved
# - Order               : order assigned, family unresolved
# - Family              : family assigned, genus unresolved
# - Genus               : genus assigned, species field absent/blank
# - Species unspecified : species field present but not resolved to a named
#                         species (e.g. "sp.", "uncultured", "unclassified")
# - Species identified  : named species assignment
#
# This reproduces the logic behind Supplementary Figure 2 while also saving
# the category counts as a CSV for transparent reporting.
# =============================================================================

required_tax_ranks <- c(
  "Phylum", "Class", "Order", "Family", "Genus", "Species"
)

missing_tax_ranks <- setdiff(
  required_tax_ranks,
  rank_names(ps_milk)
)

if (length(missing_tax_ranks) > 0) {
  stop(
    "Taxonomy table is missing required ranks: ",
    paste(missing_tax_ranks, collapse = ", ")
  )
}

tax_mat_milk <- as.data.frame(
  tax_table(ps_milk),
  stringsAsFactors = FALSE
)

clean_tax_value <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "__", "_", "NA")] <- NA_character_
  x
}

for (rk in required_tax_ranks) {
  tax_mat_milk[[rk]] <- clean_tax_value(
    tax_mat_milk[[rk]]
  )
}

is_unspecified_species <- function(x) {

  if (is.na(x) || x == "") {
    return(FALSE)
  }

  grepl(
    paste(
      c(
        "(^|\\s)sp\\.?($|\\s)",
        "uncultured",
        "unclassified",
        "unidentified",
        "unknown",
        "bacterium",
        "metagenome",
        "environmental sample"
      ),
      collapse = "|"
    ),
    x,
    ignore.case = TRUE
  )
}

deepest_assignment <- vapply(
  seq_len(nrow(tax_mat_milk)),
  function(i) {

    x <- tax_mat_milk[i, , drop = FALSE]

    if (is.na(x$Phylum)) {
      return("Unclassified")
    }

    if (is.na(x$Class)) {
      return("Phylum")
    }

    if (is.na(x$Order)) {
      return("Class")
    }

    if (is.na(x$Family)) {
      return("Order")
    }

    if (is.na(x$Genus)) {
      return("Family")
    }

    if (is.na(x$Species)) {
      return("Genus")
    }

    if (is_unspecified_species(x$Species)) {
      return("Species unspecified")
    }

    "Species identified"
  },
  character(1)
)

supp2_levels <- c(
  "Species identified",
  "Species unspecified",
  "Genus",
  "Family",
  "Order",
  "Class",
  "Phylum",
  "Unclassified"
)

supp2_counts <- as.data.frame(
  table(
    factor(
      deepest_assignment,
      levels = supp2_levels
    )
  ),
  stringsAsFactors = FALSE
)

colnames(supp2_counts) <- c(
  "Taxonomic_resolution",
  "Number"
)

write.csv(
  supp2_counts,
  file.path(
    OUT_TAB,
    "SupplementaryFigure2_Taxonomic_Resolution_Counts.csv"
  ),
  row.names = FALSE
)

supp2_palette <- setNames(
  MASTER_COLORS[seq_along(supp2_levels)],
  supp2_levels
)

p_supp2 <- ggplot(
  supp2_counts,
  aes(
    x = Taxonomic_resolution,
    y = Number,
    fill = Taxonomic_resolution
  )
) +
  geom_col(
    width = 0.68,
    colour = "black",
    linewidth = 0.35
  ) +
  geom_text(
    aes(label = Number),
    hjust = -0.12,
    size = 4.5,
    fontface = "bold"
  ) +
  coord_flip() +
  scale_fill_manual(
    values = supp2_palette,
    guide = "none"
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(0, 0.10)
    )
  ) +
  labs(
    title = "Supplementary Figure 2. Identified taxa from raw bovine milk samples",
    x = NULL,
    y = "Number"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 16
    ),
    axis.text.y = element_text(
      face = "bold",
      size = 12
    ),
    axis.title.x = element_text(
      face = "bold"
    ),
    panel.grid.minor = element_blank()
  )

print(p_supp2)

ggsave(
  file.path(
    OUT_FIG,
    "SupplementaryFigure2_Identified_Taxa.png"
  ),
  p_supp2,
  width = 9,
  height = 6.5,
  units = "in",
  dpi = 300
)

ggsave(
  file.path(
    OUT_FIG,
    "SupplementaryFigure2_Identified_Taxa.pdf"
  ),
  p_supp2,
  width = 9,
  height = 6.5,
  units = "in",
  device = cairo_pdf
)


# =============================================================================
# HELPER: COW / PAIR IDENTIFIER FOR SUPPLEMENTARY FIGURES 3–4
# =============================================================================
#
# A dedicated cow/animal ID column is preferred. If none is present, the
# previous study naming convention is used to derive the cow ID from
# Sample_Name. The script stops rather than silently inventing a pairing ID.

sample_meta_milk <- as.data.frame(
  sample_data(ps_rel)
)

sample_meta_milk$Sample <- rownames(
  sample_meta_milk
)

pair_candidates <- c(
  "CowID",
  "Cow_ID",
  "AnimalID",
  "Animal_ID",
  "PairID"
)

available_pair_col <- pair_candidates[
  pair_candidates %in% colnames(sample_meta_milk)
]

if (length(available_pair_col) > 0) {

  pair_col <- available_pair_col[1]

  sample_meta_milk$CowID <- str_trim(
    as.character(
      sample_meta_milk[[pair_col]]
    )
  )

} else if (
  "Sample_Name" %in% colnames(sample_meta_milk)
) {

  sample_meta_milk$CowID <- str_remove(
    str_trim(
      as.character(
        sample_meta_milk$Sample_Name
      )
    ),
    "_(AB|NI)2?$"
  )

} else {

  stop(
    paste0(
      "Supplementary Figures 3–4 require a cow/pair identifier. ",
      "Add one of CowID, Cow_ID, AnimalID, Animal_ID or PairID to the ",
      "metadata, or include Sample_Name so CowID can be derived."
    )
  )
}

sample_meta_milk <- sample_meta_milk %>%
  mutate(
    Lactation = factor(
      as.character(Lactation),
      levels = c("Early", "Late")
    ),
    Stage_short = ifelse(
      Lactation == "Early",
      "E",
      "L"
    ),
    Sample_plot_id = paste(
      Farm,
      CowID,
      Lactation,
      sep = "__"
    ),
    Sample_display = paste0(
      CowID,
      "-",
      Stage_short
    )
  ) %>%
  arrange(
    Farm,
    CowID,
    Lactation
  )

sample_plot_levels <- sample_meta_milk$Sample_plot_id
sample_plot_labels <- setNames(
  sample_meta_milk$Sample_display,
  sample_meta_milk$Sample_plot_id
)


# =============================================================================
# SUPPLEMENTARY FIGURE 3
# Relative abundance of the top 20 dominant microbial taxa across individual
# cow milk samples during Early and Late lactation in Farms A and B
# =============================================================================

prepare_individual_sample_data <- function(
    top_result
) {

  top_result$data %>%
    left_join(
      sample_meta_milk %>%
        select(
          Sample,
          CowID,
          Sample_plot_id,
          Sample_display
        ),
      by = "Sample"
    ) %>%
    mutate(
      Sample_plot_id = factor(
        Sample_plot_id,
        levels = sample_plot_levels
      )
    )
}

supp3_phylum_df <- prepare_individual_sample_data(
  phylum_top
)

supp3_genus_df <- prepare_individual_sample_data(
  genus_top
)

supp3_species_df <- prepare_individual_sample_data(
  species_top
)

plot_individual_top20 <- function(
    data,
    palette,
    legend_title,
    title_text,
    panel_tag
) {

  ggplot(
    data,
    aes(
      x = Sample_plot_id,
      y = Abundance * 100,
      fill = Taxon
    )
  ) +
    geom_col(
      width = 0.90
    ) +
    facet_wrap(
      ~ Farm,
      scales = "free_x",
      nrow = 1
    ) +
    scale_fill_manual(
      values = palette,
      name = legend_title
    ) +
    scale_x_discrete(
      labels = sample_plot_labels,
      drop = TRUE
    ) +
    scale_y_continuous(
      limits = c(0, 100),
      expand = c(0, 0)
    ) +
    labs(
      tag = panel_tag,
      title = title_text,
      x = "Cow sample (E = Early, L = Late)",
      y = "Relative abundance (%)"
    ) +
    publication_theme() +
    theme(
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 8
      ),
      legend.text = element_text(
        size = 8
      ),
      legend.title = element_text(
        size = 10,
        face = "bold"
      )
    )
}

p_supp3_a <- plot_individual_top20(
  supp3_phylum_df,
  phylum_top$palette,
  "Phylum",
  "Phylum-level composition",
  "A"
)

p_supp3_b <- plot_individual_top20(
  supp3_genus_df,
  genus_top$palette,
  "Genus",
  "Genus-level composition",
  "B"
)

p_supp3_c <- plot_individual_top20(
  supp3_species_df,
  species_top$palette,
  "Species",
  "Species-level composition",
  "C"
)

supplementary_figure3 <- (
  p_supp3_a /
    p_supp3_b /
    p_supp3_c
)

print(
  supplementary_figure3
)

ggsave(
  file.path(
    OUT_FIG,
    "SupplementaryFigure3_Top20_Dominant_Taxa_Individual_Cows.png"
  ),
  supplementary_figure3,
  width = 22,
  height = 24,
  units = "in",
  dpi = 300,
  limitsize = FALSE
)

ggsave(
  file.path(
    OUT_FIG,
    "SupplementaryFigure3_Top20_Dominant_Taxa_Individual_Cows.pdf"
  ),
  supplementary_figure3,
  width = 22,
  height = 24,
  units = "in",
  device = cairo_pdf,
  limitsize = FALSE
)


# =============================================================================
# SUPPLEMENTARY FIGURE 4
# Distribution of the top 20 rare bacterial taxa across milk samples
# =============================================================================
#
# Operational definition used here:
# A rare species has a mean relative abundance >0 and <= 0.01% across the
# 42 milk samples. The 20 most abundant taxa within this rare set are selected.
#
# This threshold is explicit and can be changed below if the manuscript uses
# a different rare-taxon definition.

RARE_THRESHOLD_PERCENT <- 0.01
RARE_TOP_N <- 20

ps_species_all <- tax_glom(
  ps_rel,
  taxrank = "Species",
  NArm = FALSE
)

rare_species_df <- psmelt(
  ps_species_all
)

rare_species_df$Genus <- as.character(
  rare_species_df$Genus
)

rare_species_df$Species <- as.character(
  rare_species_df$Species
)

rare_species_df$Genus[
  is.na(rare_species_df$Genus) |
    rare_species_df$Genus == ""
] <- "Unclassified"

rare_species_df$Species[
  is.na(rare_species_df$Species) |
    rare_species_df$Species == ""
] <- "Unresolved"

rare_species_df$Taxon <- ifelse(
  rare_species_df$Species == "Unresolved",
  paste(
    rare_species_df$Genus,
    "sp."
  ),
  ifelse(
    startsWith(
      rare_species_df$Species,
      rare_species_df$Genus
    ),
    rare_species_df$Species,
    paste(
      rare_species_df$Genus,
      rare_species_df$Species
    )
  )
)

rare_species_sample <- rare_species_df %>%
  group_by(
    Sample,
    Farm,
    Lactation,
    Taxon
  ) %>%
  summarise(
    Abundance = sum(Abundance),
    .groups = "drop"
  ) %>%
  left_join(
    sample_meta_milk %>%
      select(
        Sample,
        CowID,
        Sample_plot_id,
        Sample_display
      ),
    by = "Sample"
  )

rare_taxa_summary <- rare_species_sample %>%
  group_by(Taxon) %>%
  summarise(
    Mean_relative_abundance_percent =
      mean(Abundance * 100),
    Prevalence =
      mean(Abundance > 0),
    .groups = "drop"
  ) %>%
  filter(
    Mean_relative_abundance_percent > 0,
    Mean_relative_abundance_percent <=
      RARE_THRESHOLD_PERCENT
  ) %>%
  arrange(
    desc(
      Mean_relative_abundance_percent
    )
  )

if (nrow(rare_taxa_summary) == 0) {
  stop(
    paste0(
      "No species met the rare-taxon definition (mean abundance <= ",
      RARE_THRESHOLD_PERCENT,
      "%)."
    )
  )
}

n_rare_use <- min(
  RARE_TOP_N,
  nrow(rare_taxa_summary)
)

top20_rare_species <- rare_taxa_summary %>%
  slice_head(
    n = n_rare_use
  ) %>%
  pull(Taxon)

if (n_rare_use < RARE_TOP_N) {
  warning(
    "Fewer than 20 species met the rare-taxon threshold; plotting ",
    n_rare_use,
    " taxa."
  )
}

write.csv(
  rare_taxa_summary,
  file.path(
    OUT_TAB,
    "SupplementaryFigure4_Rare_Species_QC.csv"
  ),
  row.names = FALSE
)

write.csv(
  rare_taxa_summary %>%
    slice_head(
      n = n_rare_use
    ),
  file.path(
    OUT_TAB,
    "SupplementaryFigure4_Top20_Rare_Species.csv"
  ),
  row.names = FALSE
)

rare_plot_df <- rare_species_sample %>%
  filter(
    Taxon %in% top20_rare_species
  ) %>%
  mutate(
    Taxon = factor(
      Taxon,
      levels = top20_rare_species
    ),
    Sample_plot_id = factor(
      Sample_plot_id,
      levels = sample_plot_levels
    )
  )

rare_palette <- make_top_palette(
  top20_rare_species
)

# Remove the "Other" entry because Supplementary Figure 4 shows only the
# selected rare taxa, rather than the remainder of the community.
rare_palette <- rare_palette[
  names(rare_palette) != "Other"
]

p_supp4_a <- ggplot(
  rare_plot_df,
  aes(
    x = Lactation,
    y = Abundance * 100,
    fill = Taxon
  )
) +
  geom_bar(
    stat = "summary",
    fun = "mean",
    position = "stack",
    width = 0.75
  ) +
  scale_fill_manual(
    values = rare_palette,
    name = "Rare species"
  ) +
  labs(
    tag = "A",
    title = "Rare species by lactation stage",
    x = "Lactation stage",
    y = "Mean relative abundance (%)"
  ) +
  publication_theme()

p_supp4_b <- ggplot(
  rare_plot_df,
  aes(
    x = Farm,
    y = Abundance * 100,
    fill = Taxon
  )
) +
  geom_bar(
    stat = "summary",
    fun = "mean",
    position = "stack",
    width = 0.75
  ) +
  scale_fill_manual(
    values = rare_palette,
    name = "Rare species"
  ) +
  labs(
    tag = "B",
    title = "Rare species by farm",
    x = "Farm",
    y = "Mean relative abundance (%)"
  ) +
  publication_theme()

p_supp4_c <- ggplot(
  rare_plot_df,
  aes(
    x = Sample_plot_id,
    y = Abundance * 100,
    fill = Taxon
  )
) +
  geom_col(
    width = 0.90
  ) +
  facet_wrap(
    ~ Farm,
    scales = "free_x",
    nrow = 1
  ) +
  scale_x_discrete(
    labels = sample_plot_labels,
    drop = TRUE
  ) +
  scale_fill_manual(
    values = rare_palette,
    name = "Rare species"
  ) +
  labs(
    tag = "C",
    title = "Rare species across individual cow milk samples",
    x = "Cow sample (E = Early, L = Late)",
    y = "Relative abundance (%)"
  ) +
  publication_theme() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 8
    ),
    legend.text = element_text(
      size = 8
    )
  )

supplementary_figure4 <- (
  p_supp4_a /
    p_supp4_b /
    p_supp4_c
)

print(
  supplementary_figure4
)

ggsave(
  file.path(
    OUT_FIG,
    "SupplementaryFigure4_Top20_Rare_Species.png"
  ),
  supplementary_figure4,
  width = 20,
  height = 22,
  units = "in",
  dpi = 300,
  limitsize = FALSE
)

ggsave(
  file.path(
    OUT_FIG,
    "SupplementaryFigure4_Top20_Rare_Species.pdf"
  ),
  supplementary_figure4,
  width = 20,
  height = 22,
  units = "in",
  device = cairo_pdf,
  limitsize = FALSE
)


# =============================================================================
# SAVE UPDATED SESSION INFORMATION
# =============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    OUT_TAB,
    "Figure2_and_SuppFigures2_3_4_R_SessionInfo.txt"
  )
)

cat("\n============================================================\n")
cat("FIGURE 2 + SUPPLEMENTARY FIGURES 2–4 COMPLETED\n")
cat("============================================================\n")
cat("Main Figure 2 and Supplementary Figures 2, 3 and 4 were generated.\n")
cat("Figures saved to:\n", OUT_FIG, "\n")
cat("Tables/QC files saved to:\n", OUT_TAB, "\n")
cat("============================================================\n")
