# =============================================================================
# Supplementary Figure 2
# Number of distinct bacterial taxa at each taxonomic rank
# =============================================================================
#
# Purpose
# -------
# Count the ASVs, reads and distinct taxa (phylum to species) in the final
# feature table of the 42 milk samples, i.e. after removal of mitochondrial and
# chloroplast sequences (01_milk_16S_metataxonomic_analysis, step 06).
# A taxon is counted once per rank, however many ASVs are assigned to it.
# ASVs without an assignment at a given rank are not counted as a taxon at that
# rank; their numbers are reported separately.
#
# Input files
# -----------
# data/processed/feature-table_milk_42samples.tsv  ASV counts (first column = ASV ID)
# data/processed/taxonomy.tsv                      Greengenes2 taxonomy
#
# Outputs (results/taxa_per_rank/)
# --------------------------------
# SupplementaryFigure2_taxa_per_rank.tiff / .png   600 dpi
# SupplementaryFigure2_counts.csv                   distinct taxa per rank
# SupplementaryFigure2_summary.txt                  ASVs, reads, unassigned ASVs
#
# Run from the repository root:
#   Rscript 02_downstream_R_analysis/07_taxa_per_rank_SupplementaryFigure2.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

FEATURE_FILE  <- file.path("data", "processed", "feature-table_milk_42samples.tsv")
TAXONOMY_FILE <- file.path("data", "processed", "taxonomy.tsv")

OUT_DIR <- file.path("results", "taxa_per_rank")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

for (f in c(FEATURE_FILE, TAXONOMY_FILE)) {
  if (!file.exists(f)) stop("File not found: ", f)
}


# -----------------------------------------------------------------------------
# 1. Read data
# -----------------------------------------------------------------------------
# The QIIME 2 export starts with the header line "#OTU ID ...";
# comment.char = "" keeps that line as the header.

ft <- read.delim(FEATURE_FILE, check.names = FALSE, row.names = 1,
                 comment.char = "", quote = "")
ft <- as.matrix(ft)
ft <- ft[rowSums(ft) > 0, , drop = FALSE]

tax <- read.delim(TAXONOMY_FILE, check.names = FALSE, quote = "")
names(tax)[1] <- "FeatureID"

ranks <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")

tx <- tax %>%
  filter(FeatureID %in% rownames(ft)) %>%
  separate(Taxon, c("Domain", ranks), sep = ";\\s*", fill = "right") %>%
  mutate(across(all_of(c("Domain", ranks)),
                ~ na_if(trimws(sub("^[a-z]__", "", .x)), "")))

if (nrow(tx) != nrow(ft)) {
  warning(nrow(ft) - nrow(tx), " ASVs in the feature table have no taxonomy entry.")
}


# -----------------------------------------------------------------------------
# 2. Counts
# -----------------------------------------------------------------------------

counts <- tibble(
  Rank   = ranks,
  Number = vapply(ranks, function(r) n_distinct(na.omit(tx[[r]])), integer(1))
)
write.csv(counts, file.path(OUT_DIR, "SupplementaryFigure2_counts.csv"), row.names = FALSE)

n_asv <- nrow(ft)
pct   <- function(x) sprintf("%s (%.1f%%)", format(x, big.mark = ","), 100 * x / n_asv)

summary_lines <- c(
  sprintf("Samples: %d", ncol(ft)),
  sprintf("ASVs: %s", format(n_asv, big.mark = ",")),
  sprintf("Reads: %s", format(sum(ft), big.mark = ",")),
  "",
  "Distinct taxa per rank:",
  sprintf("  %-8s %s", counts$Rank, format(counts$Number, big.mark = ",")),
  "",
  sprintf("ASVs not assigned to a phylum: %s", pct(sum(is.na(tx$Phylum)))),
  sprintf("ASVs not assigned to a genus:  %s", pct(sum(is.na(tx$Genus)))),
  sprintf("ASVs not assigned to a species: %s", pct(sum(is.na(tx$Species))))
)
writeLines(summary_lines, file.path(OUT_DIR, "SupplementaryFigure2_summary.txt"))
cat(summary_lines, sep = "\n")


# -----------------------------------------------------------------------------
# 3. Figure
# -----------------------------------------------------------------------------

plot_data <- counts %>% mutate(Rank = factor(Rank, levels = rev(ranks)))

p <- ggplot(plot_data, aes(x = Number, y = Rank, fill = Rank)) +
  geom_col(colour = "black", width = 0.7) +
  geom_text(aes(label = scales::comma(Number)), hjust = -0.15,
            size = 5, fontface = "bold") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12)),
                     labels = scales::comma) +
  scale_fill_manual(values = c(Phylum = "#E67E22", Class = "#00FA9A",
                               Order = "#3CB371", Family = "#A8DADC",
                               Genus = "#3498DB", Species = "#0000CD")) +
  labs(x = "Number of distinct taxa", y = NULL) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        panel.grid.major.y = element_blank(),
        axis.text = element_text(face = "bold", colour = "black"),
        axis.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "SupplementaryFigure2_taxa_per_rank.tiff"),
       p, width = 8, height = 5.5, dpi = 600, compression = "lzw")
ggsave(file.path(OUT_DIR, "SupplementaryFigure2_taxa_per_rank.png"),
       p, width = 8, height = 5.5, dpi = 600)

cat("\nOutputs saved in:", OUT_DIR, "\n")
