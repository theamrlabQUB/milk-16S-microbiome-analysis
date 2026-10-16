# =============================================================================
# Supplementary Table S16 and Supplementary Figure 8
# Agreement between the Lactoscan milk analyser and reference methods
# =============================================================================
#
# Purpose
# -------
# Compare Lactoscan fat and protein with reference methods in the validation
# subset (first-collection early-lactation milk from 50 cows):
#   - fat:     Soxhlet extraction (AOAC 963.15)
#   - protein: Bradford assay (BSA standard)
#
# Unit conversions
# ----------------
# Soxhlet fat is measured on oven-dried milk (% of dry matter, DM):
#   Soxhlet fat (% of milk) = Soxhlet fat (% DM) x total solids (%) / 100
# Bradford protein (mg/mL) is converted to g/100 mL by dividing by 10
# (the input file already holds g/100 mL).
#
# Statistics (overall and per farm)
# ---------------------------------
# - Pearson correlation
# - Bland-Altman analysis (Bland and Altman, 1986):
#     difference = Lactoscan - reference
#     bias = mean difference; 95% limits of agreement = bias +/- 1.96 SD
# - Paired t-test of the bias
# - Proportional bias: regression of the difference on the mean of the two
#   methods
# - Coefficient of variation of each method
#
# Input file
# ----------
# data/validation/validation_lactoscan_soxhlet_bradford.tsv
#   Cow      cow ID
#   Farm     A or B
#   Stage    lactation stage (Early)
#   TS       total solids (%), oven drying
#   SoxDM    Soxhlet fat (% of dry matter)
#   LacFat   Lactoscan fat (%)
#   BradProt Bradford protein (g/100 mL)
#   LacProt  Lactoscan protein (%)
#
# Outputs (results/method_agreement/)
# -----------------------------------
# TableS16A_agreement_summary.csv
# TableS16B_sample_level_data.csv
# SupplementaryFigure8_method_agreement.tiff / .png
#
# Run from the repository root:
#   Rscript 02_downstream_R_analysis/08_method_agreement_TableS16_SupplementaryFigure8.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

INPUT_FILE <- file.path("data", "validation", "validation_lactoscan_soxhlet_bradford.tsv")
OUT_DIR    <- file.path("results", "method_agreement")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(INPUT_FILE)) stop("File not found: ", INPUT_FILE)


# -----------------------------------------------------------------------------
# 1. Data
# -----------------------------------------------------------------------------

val <- read.delim(INPUT_FILE, check.names = FALSE) %>%
  mutate(Farm     = as.character(Farm),
         SoxFat   = SoxDM * TS / 100,
         FatDiff  = LacFat - SoxFat,
         ProtDiff = LacProt - BradProt)

write.csv(val %>% select(Cow, Farm, Stage, TS, SoxDM, SoxFat, LacFat, FatDiff,
                         BradProt, LacProt, ProtDiff),
          file.path(OUT_DIR, "TableS16B_sample_level_data.csv"), row.names = FALSE)


# -----------------------------------------------------------------------------
# 2. Agreement statistics
# -----------------------------------------------------------------------------

fmt_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
cv    <- function(x) 100 * sd(x) / mean(x)

agreement <- function(lac, ref) {
  d   <- lac - ref
  m   <- (lac + ref) / 2
  ct  <- cor.test(lac, ref)
  tt  <- t.test(lac, ref, paired = TRUE)
  pb  <- summary(lm(d ~ m))$coefficients
  bias <- mean(d); sdd <- sd(d)
  tibble(
    n                    = length(d),
    Lactoscan_mean_SD    = sprintf("%.2f ± %.2f", mean(lac), sd(lac)),
    Reference_mean_SD    = sprintf("%.2f ± %.2f", mean(ref), sd(ref)),
    Lactoscan_CV_pct     = round(cv(lac), 1),
    Reference_CV_pct     = round(cv(ref), 1),
    Pearson_r            = round(unname(ct$estimate), 2),
    P_r                  = fmt_p(ct$p.value),
    Bias                 = round(bias, 2),
    SD_differences       = round(sdd, 2),
    LoA_95               = sprintf("%.2f to %.2f", bias - 1.96 * sdd, bias + 1.96 * sdd),
    P_paired_t           = fmt_p(tt$p.value),
    Proportional_bias    = sprintf("%.2f (%s)", pb["m", "Estimate"], fmt_p(pb["m", "Pr(>|t|)"]))
  )
}

groups <- list(All = val, `Farm A` = filter(val, Farm == "A"), `Farm B` = filter(val, Farm == "B"))

summary_tab <- bind_rows(
  imap_dfr(groups, ~ agreement(.x$LacFat,  .x$SoxFat)   %>% mutate(Trait = "Fat (%; Lactoscan vs Soxhlet)", Group = .y)),
  imap_dfr(groups, ~ agreement(.x$LacProt, .x$BradProt) %>% mutate(Trait = "Protein (Lactoscan % vs Bradford g/100 mL)", Group = .y))
) %>% relocate(Trait, Group)

write.csv(summary_tab, file.path(OUT_DIR, "TableS16A_agreement_summary.csv"), row.names = FALSE)
print(summary_tab, width = Inf)


# -----------------------------------------------------------------------------
# 3. Figure: (a) fat scatter, (b) fat Bland-Altman,
#            (c) protein scatter, (d) protein Bland-Altman
# -----------------------------------------------------------------------------

farm_cols   <- c(A = "#1F77B4", B = "#2CA02C")
farm_shapes <- c(A = 16, B = 17)
farm_labels <- c(A = "Farm A", B = "Farm B")

scatter_plot <- function(x, y, xlab, ylab) {
  ggplot(val, aes(.data[[x]], .data[[y]], colour = Farm, shape = Farm)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
    geom_point(size = 2.2, alpha = 0.85) +
    scale_colour_manual(values = farm_cols, labels = farm_labels) +
    scale_shape_manual(values = farm_shapes, labels = farm_labels) +
    labs(x = xlab, y = ylab, colour = NULL, shape = NULL) +
    theme_bw(base_size = 11)
}

ba_plot <- function(lac, ref, ylab) {
  d <- val[[lac]] - val[[ref]]; m <- (val[[lac]] + val[[ref]]) / 2
  bias <- mean(d); lo <- bias - 1.96 * sd(d); hi <- bias + 1.96 * sd(d)
  x_lab <- max(m)
  ggplot(tibble(m, d, Farm = val$Farm), aes(m, d, colour = Farm, shape = Farm)) +
    geom_hline(yintercept = bias, colour = "black") +
    geom_hline(yintercept = c(lo, hi), linetype = "dashed", colour = "black") +
    geom_point(size = 2.2, alpha = 0.85) +
    annotate("text", x = x_lab, y = c(hi, bias, lo), vjust = -0.4, hjust = 1, size = 3,
             label = c(sprintf("+1.96 SD: %.2f", hi), sprintf("mean: %.2f", bias),
                       sprintf("\u22121.96 SD: %.2f", lo))) +
    scale_colour_manual(values = farm_cols, labels = farm_labels) +
    scale_shape_manual(values = farm_shapes, labels = farm_labels) +
    labs(x = "Mean of the two methods", y = ylab, colour = NULL, shape = NULL) +
    theme_bw(base_size = 11)
}

fig <- (scatter_plot("SoxFat", "LacFat", "Soxhlet fat (% of milk)", "Lactoscan fat (%)") |
        ba_plot("LacFat", "SoxFat", "Lactoscan \u2212 Soxhlet fat (percentage points)")) /
       (scatter_plot("BradProt", "LacProt", "Bradford protein (g/100 mL)", "Lactoscan protein (%)") |
        ba_plot("LacProt", "BradProt", "Lactoscan \u2212 Bradford protein")) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom", plot.tag = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "SupplementaryFigure8_method_agreement.tiff"),
       fig, width = 9, height = 8.5, dpi = 600, compression = "lzw")
ggsave(file.path(OUT_DIR, "SupplementaryFigure8_method_agreement.png"),
       fig, width = 9, height = 8.5, dpi = 600)

cat("\nOutputs saved in:", OUT_DIR, "\n")
