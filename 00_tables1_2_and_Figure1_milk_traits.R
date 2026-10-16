# =============================================================================
# Tables 1–2 and Figure 1
# Milk nutritional/production traits and PCA analysis
# =============================================================================
#
# Outputs
# -------
# Table 1:
#   Milk composition by lactation stage, parity group and farm
#   (mean ± SEM; paired lactation-stage tests; parity and farm comparisons)
#
# Table 2:
#   Lactation-stage comparison for Included vs Excluded sample sets
#   (mean ± SEM; paired Early-vs-Late tests using Cow ID)
#
# Figure 1:
#   PCA biplot of milk nutritional and production traits by farm
#
# Figure 1 statistical workflow
# -----------------------------
# 1. Use the six milk traits:
#      fat, protein, lactose, Milk_energy, Milk_Production, Minerals
# 2. Keep complete observations for these traits and Farm.
# 3. Standardise the six traits.
# 4. Calculate Mahalanobis distance in the standardised multivariate space.
# 5. Exclude observations exceeding the chi-square threshold at P < 0.001
#    with df equal to the number of PCA variables.
# 6. Perform PCA using centred and scaled variables.
# 7. Plot PC1 vs PC2, 95% confidence ellipses, and variable loading arrows.
#
# IMPORTANT
# ---------
# Update the two input paths below before running the script.
# =============================================================================


# =============================================================================
# 0. PACKAGES
# =============================================================================

required_packages <- c(
  "readxl",
  "tidyverse",
  "openxlsx",
  "ggplot2",
  "ggrepel"
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
  library(tidyverse)
  library(openxlsx)
  library(ggplot2)
  library(ggrepel)
})


# =============================================================================
# 1. INPUT / OUTPUT PATHS
# =============================================================================

TABLE1_FIG1_FILE <- "path/to/ALL.xlsx"
TABLE2_FILE      <- "path/to/Nutritional data all sample.xlsx"

OUT_DIR <- "results/Table1_Table2_Figure1"
OUT_TAB <- file.path(OUT_DIR, "tables")
OUT_FIG <- file.path(OUT_DIR, "figures")

dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(TABLE1_FIG1_FILE)) {
  stop("Table 1 / Figure 1 input file not found: ", TABLE1_FIG1_FILE)
}

if (!file.exists(TABLE2_FILE)) {
  stop("Table 2 input file not found: ", TABLE2_FILE)
}


# =============================================================================
# 2. SHARED SETTINGS
# =============================================================================

traits <- c(
  "fat",
  "protein",
  "lactose",
  "Milk_energy",
  "Milk_Production",
  "Minerals"
)

trait_labels <- c(
  "fat"             = "Fat (%)",
  "protein"         = "Protein (%)",
  "lactose"         = "Lactose (%)",
  "Milk_energy"     = "Milk energy (kcal/L)",
  "Milk_Production" = "Milk production (L/day)",
  "Minerals"        = "Minerals (%)"
)

pca_labels <- c(
  "fat"             = "fat",
  "protein"         = "protein",
  "lactose"         = "lactose",
  "Milk_energy"     = "milk energy",
  "Milk_Production" = "milk production",
  "Minerals"        = "minerals"
)

format_p <- function(p) {
  dplyr::case_when(
    is.na(p)   ~ "",
    p < 0.0001 ~ "<0.0001",
    TRUE       ~ sprintf("%.4f", p)
  )
}

safe_independent_ttest <- function(x, g) {

  temp <- tibble(x = x, g = g) %>%
    filter(!is.na(x), !is.na(g))

  if (n_distinct(temp$g) != 2) {
    return(NA_real_)
  }

  group_n <- table(temp$g)

  if (any(group_n < 2)) {
    return(NA_real_)
  }

  tryCatch(
    t.test(x ~ g, data = temp)$p.value,
    error = function(e) NA_real_
  )
}


# =============================================================================
# 3. TABLE 1
# Milk composition by lactation stage, parity and farm
# =============================================================================

data_t1 <- read_excel(
  TABLE1_FIG1_FILE
) %>%
  mutate(
    parity = as.numeric(parity),
    `Cow ID` = as.character(`Cow ID`),

    Lactation_days = case_when(
      str_to_lower(trimws(as.character(Lactation_days))) == "early" ~ "Early",
      str_to_lower(trimws(as.character(Lactation_days))) == "late"  ~ "Late",
      TRUE ~ as.character(Lactation_days)
    ),

    Farm = case_when(
      Farm == "Aberdeen" ~ "A",
      Farm == "NI"       ~ "B",
      Farm == "A"        ~ "A",
      Farm == "B"        ~ "B",
      TRUE ~ as.character(Farm)
    ),

    Parity_Group = case_when(
      parity == 1 ~ "P1",
      parity >= 2 ~ "P2+",
      TRUE ~ NA_character_
    ),

    Farm = factor(Farm, levels = c("A", "B")),
    Lactation_days = factor(
      Lactation_days,
      levels = c("Early", "Late")
    ),
    Parity_Group = factor(
      Parity_Group,
      levels = c("P1", "P2+")
    )
  )

missing_t1 <- setdiff(
  c("Cow ID", "Farm", "Lactation_days", "parity", traits),
  names(data_t1)
)

if (length(missing_t1) > 0) {
  stop(
    "Missing Table 1 / Figure 1 columns: ",
    paste(missing_t1, collapse = ", ")
  )
}

df_long_t1 <- data_t1 %>%
  pivot_longer(
    cols = all_of(traits),
    names_to = "Trait_code",
    values_to = "Value"
  ) %>%
  mutate(
    Trait = recode(
      Trait_code,
      !!!trait_labels
    )
  ) %>%
  filter(
    !is.na(Value),
    !is.na(Farm),
    !is.na(Lactation_days),
    !is.na(Parity_Group),
    !is.na(`Cow ID`)
  )

# Mean ± SEM
mean_sem_t1 <- df_long_t1 %>%
  group_by(
    Trait,
    Trait_code,
    Lactation_days,
    Farm,
    Parity_Group
  ) %>%
  summarise(
    N = n(),
    Mean = mean(Value),
    SD = sd(Value),
    SEM = SD / sqrt(N),
    .groups = "drop"
  ) %>%
  mutate(
    Mean_SEM = case_when(
      N <= 1 ~ sprintf("%.2f", Mean),
      TRUE   ~ sprintf("%.2f ± %.2f", Mean, SEM)
    ),
    Group = paste0(
      "Farm",
      Farm,
      "_",
      Parity_Group
    )
  )

mean_sem_wide_t1 <- mean_sem_t1 %>%
  select(
    Trait,
    Trait_code,
    Lactation_days,
    Group,
    Mean_SEM
  ) %>%
  pivot_wider(
    names_from = Group,
    values_from = Mean_SEM
  )

# Lactation-stage effect: paired t-test within Farm × Parity
stage_p_t1 <- df_long_t1 %>%
  group_by(
    Trait,
    Trait_code,
    Farm,
    Parity_Group
  ) %>%
  group_modify(~{

    temp <- .x %>%
      select(
        `Cow ID`,
        Lactation_days,
        Value
      ) %>%
      group_by(
        `Cow ID`,
        Lactation_days
      ) %>%
      summarise(
        Value = mean(Value, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      pivot_wider(
        names_from = Lactation_days,
        values_from = Value
      )

    if (!all(c("Early", "Late") %in% names(temp))) {
      return(
        tibble(
          N_pairs = 0,
          P_stage = NA_real_
        )
      )
    }

    temp <- temp %>%
      filter(
        !is.na(Early),
        !is.na(Late)
      )

    if (nrow(temp) >= 2) {

      test_result <- tryCatch(
        t.test(
          temp$Early,
          temp$Late,
          paired = TRUE
        ),
        error = function(e) NULL
      )

      p_val <- ifelse(
        is.null(test_result),
        NA_real_,
        test_result$p.value
      )

    } else {
      p_val <- NA_real_
    }

    tibble(
      N_pairs = nrow(temp),
      P_stage = p_val
    )
  }) %>%
  ungroup() %>%
  mutate(
    Stage_group = paste0(
      "Stage_",
      Farm,
      "_",
      Parity_Group
    )
  )

stage_p_wide_t1 <- stage_p_t1 %>%
  select(
    Trait,
    Trait_code,
    Stage_group,
    P_stage
  ) %>%
  pivot_wider(
    names_from = Stage_group,
    values_from = P_stage
  )

stage_n_wide_t1 <- stage_p_t1 %>%
  select(
    Trait,
    Trait_code,
    Stage_group,
    N_pairs
  ) %>%
  mutate(
    Stage_group = paste0(
      Stage_group,
      "_N"
    )
  ) %>%
  pivot_wider(
    names_from = Stage_group,
    values_from = N_pairs
  )

# Parity effect: Welch independent t-test within Farm × Lactation
parity_p_t1 <- df_long_t1 %>%
  group_by(
    Trait,
    Trait_code,
    Farm,
    Lactation_days
  ) %>%
  summarise(
    P_parity = safe_independent_ttest(
      Value,
      Parity_Group
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Parity_group = paste0(
      "Parity_Farm",
      Farm
    )
  )

parity_p_wide_t1 <- parity_p_t1 %>%
  select(
    Trait,
    Trait_code,
    Lactation_days,
    Parity_group,
    P_parity
  ) %>%
  pivot_wider(
    names_from = Parity_group,
    values_from = P_parity
  )

# Farm effect: Welch independent t-test within Parity × Lactation
farm_p_t1 <- df_long_t1 %>%
  group_by(
    Trait,
    Trait_code,
    Parity_Group,
    Lactation_days
  ) %>%
  summarise(
    P_farm = safe_independent_ttest(
      Value,
      Farm
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Farm_group = paste0(
      "Farm_",
      Parity_Group
    )
  )

farm_p_wide_t1 <- farm_p_t1 %>%
  select(
    Trait,
    Trait_code,
    Lactation_days,
    Farm_group,
    P_farm
  ) %>%
  pivot_wider(
    names_from = Farm_group,
    values_from = P_farm
  )

table1_combined <- mean_sem_wide_t1 %>%
  left_join(
    stage_p_wide_t1,
    by = c("Trait", "Trait_code")
  ) %>%
  left_join(
    stage_n_wide_t1,
    by = c("Trait", "Trait_code")
  ) %>%
  left_join(
    parity_p_wide_t1,
    by = c(
      "Trait",
      "Trait_code",
      "Lactation_days"
    )
  ) %>%
  left_join(
    farm_p_wide_t1,
    by = c(
      "Trait",
      "Trait_code",
      "Lactation_days"
    )
  )

required_table1_cols <- c(
  "FarmA_P1",
  "FarmA_P2+",
  "FarmB_P1",
  "FarmB_P2+",
  "Stage_A_P1",
  "Stage_A_P2+",
  "Stage_B_P1",
  "Stage_B_P2+",
  "Parity_FarmA",
  "Parity_FarmB",
  "Farm_P1",
  "Farm_P2+"
)

for (nm in required_table1_cols) {
  if (!nm %in% names(table1_combined)) {
    table1_combined[[nm]] <- NA
  }
}

table1_combined <- table1_combined %>%
  mutate(
    Trait = factor(
      Trait,
      levels = unname(trait_labels)
    )
  ) %>%
  arrange(
    Trait,
    Lactation_days
  ) %>%
  mutate(
    `Stage A P1`   = format_p(Stage_A_P1),
    `Stage A P2+`  = format_p(`Stage_A_P2+`),
    `Stage B P1`   = format_p(Stage_B_P1),
    `Stage B P2+`  = format_p(`Stage_B_P2+`),
    `Parity Farm A` = format_p(Parity_FarmA),
    `Parity Farm B` = format_p(Parity_FarmB),
    `Farm effect P1`  = format_p(Farm_P1),
    `Farm effect P2+` = format_p(`Farm_P2+`)
  ) %>%
  mutate(
    `Stage A P1` = if_else(
      Lactation_days == "Early",
      `Stage A P1`,
      ""
    ),
    `Stage A P2+` = if_else(
      Lactation_days == "Early",
      `Stage A P2+`,
      ""
    ),
    `Stage B P1` = if_else(
      Lactation_days == "Early",
      `Stage B P1`,
      ""
    ),
    `Stage B P2+` = if_else(
      Lactation_days == "Early",
      `Stage B P2+`,
      ""
    )
  )

table1 <- table1_combined %>%
  transmute(
    Trait = as.character(Trait),
    `Lactation stage` = as.character(Lactation_days),
    `Farm A, P1` = FarmA_P1,
    `Farm A, P2+` = `FarmA_P2+`,
    `Farm B, P1` = FarmB_P1,
    `Farm B, P2+` = `FarmB_P2+`,
    `Stage effect: A, P1` = `Stage A P1`,
    `Stage effect: A, P2+` = `Stage A P2+`,
    `Stage effect: B, P1` = `Stage B P1`,
    `Stage effect: B, P2+` = `Stage B P2+`,
    `Parity effect: Farm A` = `Parity Farm A`,
    `Parity effect: Farm B` = `Parity Farm B`,
    `Farm effect: P1` = `Farm effect P1`,
    `Farm effect: P2+` = `Farm effect P2+`
  )


# =============================================================================
# 4. TABLE 2
# Included vs Excluded sample sets: Early-vs-Late paired comparison
# =============================================================================

sheet_names <- excel_sheets(
  TABLE2_FILE
)

if (!all(c("Included", "Excluded") %in% sheet_names)) {
  stop(
    "Table 2 workbook must contain sheets named 'Included' and 'Excluded'."
  )
}

analyse_lactation_sheet <- function(sheet_name) {

  dat <- read_excel(
    TABLE2_FILE,
    sheet = sheet_name
  ) %>%
    mutate(
      `Cow ID` = as.character(`Cow ID`),
      Lactation_days = case_when(
        str_to_lower(
          trimws(
            as.character(Lactation_days)
          )
        ) == "early" ~ "Early",

        str_to_lower(
          trimws(
            as.character(Lactation_days)
          )
        ) == "late" ~ "Late",

        TRUE ~ as.character(Lactation_days)
      ),
      Lactation_days = factor(
        Lactation_days,
        levels = c("Early", "Late")
      )
    )

  missing_cols <- setdiff(
    c("Cow ID", "Lactation_days", traits),
    names(dat)
  )

  if (length(missing_cols) > 0) {
    stop(
      "Missing columns in sheet '",
      sheet_name,
      "': ",
      paste(missing_cols, collapse = ", ")
    )
  }

  long <- dat %>%
    pivot_longer(
      cols = all_of(traits),
      names_to = "Trait_code",
      values_to = "Value"
    ) %>%
    mutate(
      Trait = recode(
        Trait_code,
        !!!trait_labels
      )
    ) %>%
    filter(
      !is.na(Value),
      !is.na(`Cow ID`),
      !is.na(Lactation_days)
    )

  mean_sem <- long %>%
    group_by(
      Trait,
      Trait_code,
      Lactation_days
    ) %>%
    summarise(
      N = n(),
      Mean = mean(Value),
      SD = sd(Value),
      SEM = SD / sqrt(N),
      .groups = "drop"
    ) %>%
    mutate(
      Mean_SEM = case_when(
        N <= 1 ~ sprintf("%.2f", Mean),
        TRUE ~ sprintf(
          "%.2f ± %.2f",
          Mean,
          SEM
        )
      )
    )

  mean_sem_wide <- mean_sem %>%
    select(
      Trait,
      Trait_code,
      Lactation_days,
      Mean_SEM
    ) %>%
    pivot_wider(
      names_from = Lactation_days,
      values_from = Mean_SEM
    )

  stage_p <- long %>%
    group_by(
      Trait,
      Trait_code
    ) %>%
    group_modify(~{

      temp <- .x %>%
        select(
          `Cow ID`,
          Lactation_days,
          Value
        ) %>%
        group_by(
          `Cow ID`,
          Lactation_days
        ) %>%
        summarise(
          Value = mean(
            Value,
            na.rm = TRUE
          ),
          .groups = "drop"
        ) %>%
        pivot_wider(
          names_from = Lactation_days,
          values_from = Value
        )

      if (!all(c("Early", "Late") %in% names(temp))) {
        return(
          tibble(
            N_pairs = 0,
            P_value = NA_real_
          )
        )
      }

      temp <- temp %>%
        filter(
          !is.na(Early),
          !is.na(Late)
        )

      if (nrow(temp) >= 2) {

        test_result <- tryCatch(
          t.test(
            temp$Early,
            temp$Late,
            paired = TRUE
          ),
          error = function(e) NULL
        )

        p_val <- ifelse(
          is.null(test_result),
          NA_real_,
          test_result$p.value
        )

      } else {
        p_val <- NA_real_
      }

      tibble(
        N_pairs = nrow(temp),
        P_value = p_val
      )
    }) %>%
    ungroup()

  final_table <- mean_sem_wide %>%
    left_join(
      stage_p,
      by = c(
        "Trait",
        "Trait_code"
      )
    ) %>%
    mutate(
      P = format_p(P_value),
      Trait = factor(
        Trait,
        levels = unname(trait_labels)
      )
    ) %>%
    arrange(Trait) %>%
    transmute(
      Trait = as.character(Trait),
      Early,
      Late,
      `N pairs` = N_pairs,
      P
    )

  list(
    Final_Table = final_table,
    Mean_SEM = mean_sem,
    Paired_Test = stage_p
  )
}

included_results <- analyse_lactation_sheet(
  "Included"
)

excluded_results <- analyse_lactation_sheet(
  "Excluded"
)

table2 <- included_results$Final_Table %>%
  rename(
    `Included Early` = Early,
    `Included Late` = Late,
    `Included N pairs` = `N pairs`,
    `Included P` = P
  ) %>%
  left_join(
    excluded_results$Final_Table %>%
      rename(
        `Excluded Early` = Early,
        `Excluded Late` = Late,
        `Excluded N pairs` = `N pairs`,
        `Excluded P` = P
      ),
    by = "Trait"
  )


# =============================================================================
# 5. FIGURE 1
# PCA biplot of milk traits between farms
# =============================================================================

pca_data <- data_t1 %>%
  select(
    Farm,
    all_of(traits)
  ) %>%
  mutate(
    across(
      all_of(traits),
      ~ suppressWarnings(as.numeric(.x))
    )
  ) %>%
  filter(
    !is.na(Farm)
  ) %>%
  drop_na(
    all_of(traits)
  )

if (nrow(pca_data) <= length(traits) + 1) {
  stop("Too few complete observations for PCA/outlier detection.")
}

# Standardise variables for Mahalanobis distance.
trait_matrix <- as.matrix(
  pca_data[, traits, drop = FALSE]
)

trait_scaled <- scale(
  trait_matrix,
  center = TRUE,
  scale = TRUE
)

cov_mat <- cov(
  trait_scaled
)

if (det(cov_mat) == 0) {
  stop(
    "Covariance matrix is singular; Mahalanobis distance cannot be calculated."
  )
}

mahal_dist <- mahalanobis(
  trait_scaled,
  center = colMeans(trait_scaled),
  cov = cov_mat
)

mahal_cutoff <- qchisq(
  0.999,
  df = length(traits)
)

pca_qc <- pca_data %>%
  mutate(
    Mahalanobis_D2 = mahal_dist,
    Chi_square_cutoff = mahal_cutoff,
    Excluded_as_outlier = Mahalanobis_D2 > Chi_square_cutoff
  )

write.csv(
  pca_qc,
  file.path(
    OUT_TAB,
    "Figure1_Mahalanobis_Outlier_QC.csv"
  ),
  row.names = FALSE
)

pca_clean <- pca_qc %>%
  filter(
    !Excluded_as_outlier
  )

cat("\nFigure 1 PCA sample counts after outlier removal:\n")
print(
  table(
    pca_clean$Farm
  )
)

cat(
  "Outliers removed:",
  sum(pca_qc$Excluded_as_outlier),
  "\n"
)

# PCA is performed on the cleaned raw variables with center=TRUE and scale=TRUE.
pca_fit <- prcomp(
  pca_clean[, traits, drop = FALSE],
  center = TRUE,
  scale. = TRUE
)

variance_explained <- (
  pca_fit$sdev^2 /
    sum(pca_fit$sdev^2)
) * 100

scores <- as.data.frame(
  pca_fit$x[, 1:2, drop = FALSE]
) %>%
  mutate(
    Farm = pca_clean$Farm
  )

loadings <- as.data.frame(
  pca_fit$rotation[, 1:2, drop = FALSE]
) %>%
  rownames_to_column(
    "Trait_code"
  ) %>%
  mutate(
    Trait_label = recode(
      Trait_code,
      !!!pca_labels
    )
  )

# Scale loading vectors to fit naturally within the score plot.
score_radius <- min(
  max(abs(scores$PC1)),
  max(abs(scores$PC2))
)

loading_radius <- max(
  sqrt(
    loadings$PC1^2 +
      loadings$PC2^2
  )
)

arrow_multiplier <- (
  score_radius /
    loading_radius
) * 0.75

loadings <- loadings %>%
  mutate(
    PC1_plot = PC1 * arrow_multiplier,
    PC2_plot = PC2 * arrow_multiplier
  )

farm_shapes <- c(
  "A" = 16,
  "B" = 17
)

farm_labels <- c(
  "A" = "Farm A",
  "B" = "Farm B"
)

p_figure1 <- ggplot(
  scores,
  aes(
    x = PC1,
    y = PC2,
    colour = Farm,
    shape = Farm
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  stat_ellipse(
    aes(group = Farm),
    type = "norm",
    level = 0.95,
    linewidth = 1.0,
    show.legend = FALSE
  ) +
  geom_point(
    size = 3.0,
    alpha = 0.75
  ) +
  geom_segment(
    data = loadings,
    aes(
      x = 0,
      y = 0,
      xend = PC1_plot,
      yend = PC2_plot
    ),
    inherit.aes = FALSE,
    arrow = arrow(
      length = grid::unit(
        0.18,
        "cm"
      )
    ),
    linewidth = 0.6
  ) +
  geom_text_repel(
    data = loadings,
    aes(
      x = PC1_plot,
      y = PC2_plot,
      label = Trait_label
    ),
    inherit.aes = FALSE,
    size = 4.2,
    fontface = "bold",
    min.segment.length = 0,
    box.padding = 0.35
  ) +
  scale_shape_manual(
    values = farm_shapes,
    labels = farm_labels,
    name = NULL
  ) +
  scale_colour_discrete(
    labels = farm_labels,
    name = NULL
  ) +
  labs(
    x = paste0(
      "PC1 (",
      round(
        variance_explained[1],
        1
      ),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(
        variance_explained[2],
        1
      ),
      "%)"
    )
  ) +
  theme_minimal(
    base_size = 14
  ) +
  theme(
    axis.title = element_text(
      face = "bold"
    ),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

print(
  p_figure1
)

ggsave(
  file.path(
    OUT_FIG,
    "Figure1_PCA_Milk_Traits_by_Farm.png"
  ),
  p_figure1,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300
)

ggsave(
  file.path(
    OUT_FIG,
    "Figure1_PCA_Milk_Traits_by_Farm.pdf"
  ),
  p_figure1,
  width = 10,
  height = 7,
  units = "in",
  device = cairo_pdf
)

write.csv(
  scores,
  file.path(
    OUT_TAB,
    "Figure1_PCA_Scores.csv"
  ),
  row.names = FALSE
)

write.csv(
  loadings,
  file.path(
    OUT_TAB,
    "Figure1_PCA_Loadings.csv"
  ),
  row.names = FALSE
)

write.csv(
  tibble(
    Principal_component = paste0(
      "PC",
      seq_along(
        variance_explained
      )
    ),
    Variance_explained_percent =
      variance_explained
  ),
  file.path(
    OUT_TAB,
    "Figure1_PCA_Variance_Explained.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 6. EXPORT TABLES 1–2
# =============================================================================

write.xlsx(
  list(
    "Table1" = table1,
    "Table1_Mean_SEM" = mean_sem_t1,
    "Table1_Stage_Effect" = stage_p_t1,
    "Table1_Parity_Effect" = parity_p_t1,
    "Table1_Farm_Effect" = farm_p_t1,

    "Table2" = table2,
    "Table2_Included_Mean_SEM" =
      included_results$Mean_SEM,
    "Table2_Included_Paired_Test" =
      included_results$Paired_Test,
    "Table2_Excluded_Mean_SEM" =
      excluded_results$Mean_SEM,
    "Table2_Excluded_Paired_Test" =
      excluded_results$Paired_Test
  ),
  file = file.path(
    OUT_TAB,
    "Tables1_2_Milk_Traits.xlsx"
  ),
  overwrite = TRUE
)


# =============================================================================
# 7. SAVE R SESSION INFORMATION
# =============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    OUT_TAB,
    "Tables1_2_Figure1_R_SessionInfo.txt"
  )
)


# =============================================================================
# 8. FINAL REPORT
# =============================================================================

cat("\n============================================================\n")
cat("TABLES 1–2 + FIGURE 1 ANALYSIS COMPLETED\n")
cat("============================================================\n")

cat(
  "\nTable 1 rows:",
  nrow(table1),
  "\n"
)

cat(
  "Table 2 rows:",
  nrow(table2),
  "\n"
)

cat(
  "\nPCA observations before outlier removal:",
  nrow(pca_qc),
  "\n"
)

cat(
  "PCA observations retained:",
  nrow(pca_clean),
  "\n"
)

cat(
  "PCA observations excluded:",
  sum(pca_qc$Excluded_as_outlier),
  "\n"
)

cat(
  "PC1 variance explained:",
  round(variance_explained[1], 1),
  "%\n"
)

cat(
  "PC2 variance explained:",
  round(variance_explained[2], 1),
  "%\n"
)

cat(
  "\nTables saved to:",
  OUT_TAB,
  "\n"
)

cat(
  "Figure saved to:",
  OUT_FIG,
  "\n"
)

cat("============================================================\n")
