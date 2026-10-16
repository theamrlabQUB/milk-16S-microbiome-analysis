# =============================================================================
# Contaminant screening (decontam) and Lactobacillus helveticus ASVs
# Manuscript: Methods (Bioinformatic analysis), Results (differential
# abundance) and Discussion (limitations)
# =============================================================================
#
# Purpose
# -------
# A. Extract all ASVs assigned to Lactobacillus helveticus, summarise their
#    read counts, prevalence and classifier confidence, and write their
#    sequences to FASTA for BLASTn against the NCBI 16S rRNA RefSeq database.
# B. Identify potential contaminant ASVs with the frequency method of decontam
#    (Davis et al., 2018). Extraction blanks were not sequenced, so the
#    prevalence method cannot be used. The primary input is the post-PCR
#    amplicon concentration; pre-PCR DNA concentration is used as a
#    sensitivity check (threshold = 0.1).
# C. Repeat the lactation-stage PERMANOVA (distance ~ Farm + Stage,
#    9,999 permutations restricted within cow) with and without the flagged
#    ASVs.
#
# Input files
# -----------
# data/processed/sequences.fasta                  representative ASV sequences
# data/processed/taxonomy.tsv                     Greengenes2 taxonomy
# data/processed/feature-table_milk_42samples.tsv ASV counts, 42 milk samples
# data/metadata/metadata_cow.csv                  SampleID, CowID, Farm, Lactation
# DNA and amplicon concentrations of the 42 sequenced extracts are taken from
# Supplementary Table S17 (section 3).
#
# Outputs (results/decontam/)
# ---------------------------
# L_helveticus_ASVs_summary.csv
# L_helveticus_ASVs.fasta
# decontam_results_all_ASVs.csv
# decontam_flagged_ASVs.csv
# PERMANOVA_with_and_without_contaminants.csv
# decontam_summary.txt
#
# Run from the repository root:
#   Rscript 02_downstream_R_analysis/06_decontam_and_L_helveticus_TableS17.R
# =============================================================================


# =============================================================================
# 0. PACKAGES
# =============================================================================
# install.packages(c("tidyverse", "vegan", "BiocManager"))
# BiocManager::install("decontam")

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(decontam)
})


# =============================================================================
# 1. PATHS
# =============================================================================

FASTA_FILE    <- file.path("data", "processed", "sequences.fasta")
TAXONOMY_FILE <- file.path("data", "processed", "taxonomy.tsv")
FEATURE_FILE  <- file.path("data", "processed", "feature-table_milk_42samples.tsv")
METADATA_FILE <- file.path("data", "metadata", "metadata_cow.csv")

OUT_DIR <- file.path("results", "decontam")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

for (f in c(FASTA_FILE, TAXONOMY_FILE, FEATURE_FILE, METADATA_FILE)) {
  if (!file.exists(f)) stop("File not found: ", f)
}

set.seed(12345)


# =============================================================================
# 2. READ DATA
# =============================================================================

# ---- Feature table (ASVs x samples) ----
ft <- read.delim(FEATURE_FILE, check.names = FALSE, stringsAsFactors = FALSE, quote = "")
counts <- as.matrix(ft[, -1, drop = FALSE])
mode(counts) <- "numeric"
rownames(counts) <- as.character(ft[[1]])
counts <- counts[rowSums(counts) > 0, , drop = FALSE]
cat("ASVs:", nrow(counts), "| Samples:", ncol(counts), "\n")

# ---- Taxonomy ----
tax <- read.delim(TAXONOMY_FILE, check.names = FALSE, stringsAsFactors = FALSE, quote = "")
names(tax)[1] <- "FeatureID"
tax <- tax[, c("FeatureID", "Taxon", intersect("Confidence", names(tax)))]

get_rank <- function(taxon, prefix) {
  vapply(taxon, function(x) {
    r <- trimws(unlist(strsplit(x, ";", fixed = TRUE)))
    v <- r[grepl(paste0("^", prefix, "__"), r)]
    if (length(v) == 0) return(NA_character_)
    v <- trimws(sub(paste0("^", prefix, "__"), "", v[1]))
    if (v %in% c("", "_", "__")) NA_character_ else gsub("_", " ", v)
  }, character(1), USE.NAMES = FALSE)
}
tax$Genus   <- get_rank(tax$Taxon, "g")
tax$Species <- get_rank(tax$Taxon, "s")

# ---- Metadata ----
meta <- read.csv(METADATA_FILE, stringsAsFactors = FALSE, check.names = FALSE)
meta <- meta[, !(is.na(names(meta)) | trimws(names(meta)) == ""), drop = FALSE]
names(meta) <- trimws(names(meta))
need <- c("SampleID", "CowID", "Farm", "Lactation")
if (length(setdiff(need, names(meta))) > 0) {
  stop("metadata_cow.csv must contain: ", paste(need, collapse = ", "),
       "\nColumns found: ", paste(names(meta), collapse = ", "))
}
meta <- meta %>%
  mutate(across(all_of(need), ~ trimws(as.character(.x))),
         Stage = ifelse(tolower(Lactation) == "early", "Early", "Late"),
         Key   = paste(Farm, CowID, Stage, sep = "_")) %>%
  filter(SampleID %in% colnames(counts))

counts <- counts[, meta$SampleID, drop = FALSE]
cat("Samples matched to metadata:", ncol(counts), "\n")


# =============================================================================
# 3. DNA CONCENTRATIONS (Supplementary Table S17, 42 sequenced extracts)
# =============================================================================
# DNA_ng_uL     = DNA concentration of the extract before PCR (DeNovix)
# Amplicon_ng_uL = full-length 16S amplicon concentration after PCR

dna <- tribble(
  ~Farm, ~CowID, ~Stage, ~DNA_ng_uL, ~Amplicon_ng_uL,
  "A", "37", "Early", 10.617, 10.8475,
  "A", "37", "Late", 2.836, 8.8875,
  "A", "107", "Early", 5.257, 3.2300,
  "A", "107", "Late", 5.523, 7.1750,
  "A", "142", "Early", 15.582, 4.5125,
  "A", "142", "Late", 3.963, 10.3800,
  "A", "180", "Early", 19.017, 13.9100,
  "A", "180", "Late", 11.621, 4.2625,
  "A", "214", "Early", 8.167, 3.3350,
  "A", "214", "Late", 20.413, 8.5475,
  "A", "282", "Early", 5.813, 3.7350,
  "A", "282", "Late", 8.597, 4.2125,
  "A", "354", "Early", 4.090, 5.3900,
  "A", "354", "Late", 6.366, 3.0275,
  "A", "365", "Early", 3.380, 6.3850,
  "A", "365", "Late", 53.630, 11.1025,
  "B", "142", "Early", 4.570, 3.0825,
  "B", "142", "Late", 6.380, 3.9050,
  "B", "490", "Early", 14.210, 6.9425,
  "B", "490", "Late", 5.776, 3.0975,
  "B", "503", "Early", 6.562, 12.7875,
  "B", "503", "Late", 15.252, 22.9200,
  "B", "520", "Early", 6.382, 17.3950,
  "B", "520", "Late", 3.823, 3.2000,
  "B", "563", "Early", 6.758, 15.6800,
  "B", "563", "Late", 6.484, 3.3525,
  "B", "642", "Early", 4.445, 4.8750,
  "B", "642", "Late", 9.569, 4.7075,
  "B", "759", "Early", 8.794, 21.0125,
  "B", "759", "Late", 37.310, 3.5700,
  "B", "797", "Early", 14.370, 12.8625,
  "B", "797", "Late", 86.840, 9.3500,
  "B", "851", "Early", 42.730, 3.0450,
  "B", "851", "Late", 11.566, 8.3100,
  "B", "973", "Early", 5.712, 3.5550,
  "B", "973", "Late", 3.739, 4.6350,
  "B", "989", "Early", 14.480, 20.5025,
  "B", "989", "Late", 81.010, 7.4975,
  "B", "7358", "Early", 11.439, 13.4700,
  "B", "7358", "Late", 4.027, 8.3175,
  "B", "7403", "Early", 12.082, 20.9350,
  "B", "7403", "Late", 3.963, 6.4900
) %>% mutate(Key = paste(Farm, CowID, Stage, sep = "_"))

meta <- meta %>% left_join(dna %>% select(Key, DNA_ng_uL, Amplicon_ng_uL), by = "Key")

if (anyNA(meta$Amplicon_ng_uL)) {
  print(meta %>% filter(is.na(Amplicon_ng_uL)) %>% select(SampleID, CowID, Farm, Lactation))
  stop("Every sample must match a row of Supplementary Table S17 (CowID, Farm, Lactation).")
}
cat("All 42 samples matched to Table S17 DNA concentrations.\n")


# =============================================================================
# PART A. L. HELVETICUS ASVs FOR BLAST
# =============================================================================

lh <- tax %>% filter(grepl("Lactobacillus helveticus", Species))
lh <- lh %>% filter(FeatureID %in% rownames(counts))

early <- meta$SampleID[meta$Stage == "Early"]
late  <- meta$SampleID[meta$Stage == "Late"]

lh_summary <- lh %>%
  mutate(
    Total_reads       = rowSums(counts[FeatureID, , drop = FALSE]),
    Samples_present   = rowSums(counts[FeatureID, , drop = FALSE] > 0),
    Early_present_n21 = rowSums(counts[FeatureID, early, drop = FALSE] > 0),
    Late_present_n21  = rowSums(counts[FeatureID, late,  drop = FALSE] > 0)
  ) %>%
  arrange(desc(Total_reads))

cat("\nASVs assigned to L. helveticus:", nrow(lh_summary), "\n")
print(lh_summary %>% select(FeatureID, Total_reads, Samples_present,
                            Early_present_n21, Late_present_n21,
                            any_of("Confidence")))
write.csv(lh_summary, file.path(OUT_DIR, "L_helveticus_ASVs_summary.csv"), row.names = FALSE)

# ---- Extract their sequences ----
fa  <- readLines(FASTA_FILE)
hdr <- grep("^>", fa)
ids <- sub("^>(\\S+).*", "\\1", fa[hdr])
seqs <- vapply(seq_along(hdr), function(i) {
  end <- if (i < length(hdr)) hdr[i + 1] - 1 else length(fa)
  paste(fa[(hdr[i] + 1):end], collapse = "")
}, character(1))
names(seqs) <- ids

missing_seq <- setdiff(lh_summary$FeatureID, names(seqs))
if (length(missing_seq) > 0) warning("Sequences not found for: ", paste(missing_seq, collapse = ", "))

keep <- intersect(lh_summary$FeatureID, names(seqs))
fasta_out <- c(rbind(
  paste0(">", keep, "_reads", lh_summary$Total_reads[match(keep, lh_summary$FeatureID)]),
  seqs[keep]
))
writeLines(fasta_out, file.path(OUT_DIR, "L_helveticus_ASVs.fasta"))
cat("Saved L_helveticus_ASVs.fasta (", length(keep), " sequences ) for NCBI BLAST.\n", sep = "")


# =============================================================================
# PART B. DECONTAM (FREQUENCY METHOD)
# =============================================================================
# Primary input: post-PCR amplicon concentration, which reflects bacterial
# template in these host-DNA-rich milk extracts. Pre-PCR DNA concentration is
# run as a sensitivity check. Threshold = 0.1 (decontam default).

seqtab <- t(counts)                       # samples x ASVs

run_decontam <- function(conc, label) {
  res <- isContaminant(seqtab, conc = conc, method = "frequency", threshold = 0.1)
  res$FeatureID <- rownames(res)
  res$Input <- label
  res
}

dc_amp <- run_decontam(meta$Amplicon_ng_uL, "post-PCR amplicon concentration")
dc_dna <- run_decontam(meta$DNA_ng_uL,      "pre-PCR DNA concentration")

total_reads <- sum(counts)

summarise_dc <- function(res) {
  flagged <- res$FeatureID[res$contaminant %in% TRUE]
  list(
    n_flagged   = length(flagged),
    pct_reads   = 100 * sum(counts[flagged, , drop = FALSE]) / total_reads,
    flagged     = flagged,
    lh_flagged  = any(flagged %in% lh_summary$FeatureID)
  )
}
s_amp <- summarise_dc(dc_amp)
s_dna <- summarise_dc(dc_dna)

write.csv(bind_rows(dc_amp, dc_dna), file.path(OUT_DIR, "decontam_results_all_ASVs.csv"),
          row.names = FALSE)

flagged_table <- bind_rows(
  tibble(FeatureID = s_amp$flagged, Input = "post-PCR amplicon"),
  tibble(FeatureID = s_dna$flagged, Input = "pre-PCR DNA")
) %>%
  left_join(tax %>% select(FeatureID, Genus, Species), by = "FeatureID") %>%
  mutate(Total_reads = rowSums(counts[FeatureID, , drop = FALSE]),
         Pct_of_reads = round(100 * Total_reads / total_reads, 4),
         Samples_present = rowSums(counts[FeatureID, , drop = FALSE] > 0)) %>%
  arrange(Input, desc(Total_reads))
write.csv(flagged_table, file.path(OUT_DIR, "decontam_flagged_ASVs.csv"), row.names = FALSE)


# =============================================================================
# PART C. DO THE MAIN COMMUNITY RESULTS CHANGE AFTER REMOVAL?
# =============================================================================
# Same model as Figure 3: distance ~ Farm + Stage, 9,999 permutations
# restricted within cow, marginal tests.

meta$Pair <- factor(paste(meta$Farm, meta$CowID, sep = "_"))
meta$Farm <- factor(meta$Farm)
meta$Stage <- factor(meta$Stage, levels = c("Early", "Late"))

permanova_stage <- function(mat, label) {
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  x <- t(mat)
  rel <- x / rowSums(x)
  out <- list()
  for (d in c("bray", "jaccard")) {
    dist <- if (d == "bray") vegdist(rel, "bray") else vegdist(x, "jaccard", binary = TRUE)
    fit <- adonis2(dist ~ Farm + Stage, data = meta, permutations = 9999,
                   strata = meta$Pair, by = "margin")
    out[[d]] <- tibble(Dataset = label, Distance = d,
                       R2 = fit["Stage", "R2"], F = fit["Stage", "F"], P = fit["Stage", "Pr(>F)"])
  }
  bind_rows(out)
}

perm <- bind_rows(
  permanova_stage(counts, "All ASVs"),
  permanova_stage(counts[setdiff(rownames(counts), s_amp$flagged), , drop = FALSE],
                  "Contaminants removed (amplicon concentration)")
)
print(perm)
write.csv(perm, file.path(OUT_DIR, "PERMANOVA_with_and_without_contaminants.csv"), row.names = FALSE)


# =============================================================================
# SUMMARY
# =============================================================================

summary_lines <- c(
  "DECONTAM (frequency method, threshold 0.1)",
  sprintf("Input post-PCR amplicon concentration: %d ASVs flagged (%.2f%% of reads); L. helveticus flagged: %s",
          s_amp$n_flagged, s_amp$pct_reads, s_amp$lh_flagged),
  sprintf("Input pre-PCR DNA concentration:       %d ASVs flagged (%.2f%% of reads); L. helveticus flagged: %s",
          s_dna$n_flagged, s_dna$pct_reads, s_dna$lh_flagged),
  "",
  "PERMANOVA, lactation stage (Farm + Stage, permutations within cow):",
  capture.output(print(perm)),
  "",
  sprintf("L. helveticus ASVs written for BLAST: %d", length(keep))
)
writeLines(summary_lines, file.path(OUT_DIR, "decontam_summary.txt"))
cat("\n", paste(summary_lines, collapse = "\n"), "\n")
cat("\nAll outputs saved in:", OUT_DIR, "\n")
