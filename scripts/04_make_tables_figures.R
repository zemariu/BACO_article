# =============================================================================
# BACO manuscript — publication tables and figures for the Results section
# =============================================================================
#
# Target journal: Journal of Clinical Epidemiology
#
# This script generates:
#
# MAIN MANUSCRIPT
#   Table 1. Operating characteristics under the null hypothesis
#   Table 2. Detailed empirical comparison of the primary reanalysis and external COVID-19 replication
#   Figure 1. Type I error across critical-event information and denominator strength
#   Figure 2. Size-adjusted power parameterized by |Delta|
#   Figure 3. Original BACO Wald test versus Fisher exact test in the two empirical datasets
#
# SUPPLEMENT
#   Table S1. Simulation design
#   Table S2. Type I error summary across the final null grid
#   Table S3. Size-adjusted power summary across the final Delta grid
#   Table S4. Complete empirical reanalysis
#   Table S5. Comparisons with discordant inferential conclusions
#   Table S6. Comparisons significant by Fisher exact test
#   Table S7. Complete external empirical replication
#   Table S8. External replication sensitivity analyses
#   Figure S1. Sparse-event type I error profile
#   Figure S2. Relationship between BACO and Delta
#
#
# IMPORTANT:
#   - This script does NOT rerun simulations.
#   - It uses the finalized simulation and empirical output CSV files.
#   - Fisher exact is treated as the primary count-based discrepancy test.
#   - Delta is the discrepancy magnitude.
#   - Fieller is used for BACO ratio inference.
#   - Original BACO Wald results are retained only as a historical comparator.
#
# =============================================================================

# -----------------------------------------------------------------------------
# 0. USER CONFIGURATION
# -----------------------------------------------------------------------------

ROOT_DIR <- getwd()

# The helper below searches ROOT_DIR and common result subfolders.
SEARCH_DIRS <- c(
  ROOT_DIR,
  file.path(ROOT_DIR, "results", "empirical"),
  file.path(ROOT_DIR, "results", "simulation_final_delta"),
  file.path(ROOT_DIR, "results", "simulation_zero_event"),
  file.path(ROOT_DIR, "frozen_results", "simulation")
)

OUTPUT_DIR <- file.path(ROOT_DIR, "manuscript_outputs")
MAIN_DIR <- file.path(OUTPUT_DIR, "main_manuscript")
SUPP_DIR <- file.path(OUTPUT_DIR, "supplement")
DATA_DIR <- file.path(OUTPUT_DIR, "figure_table_data")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(MAIN_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUPP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

# Remove obsolete main-table files from earlier script versions so that rerunning
# the script does not leave a stale standalone Table 3 in the output directory.
unlink(
  c(
    file.path(MAIN_DIR, "Table_2_empirical_reanalysis_summary.docx"),
    file.path(MAIN_DIR, "Table_3_external_replication_comparison.docx"),
    file.path(DATA_DIR, "Table_3_underlying_data.csv")
  ),
  force = TRUE
)


# -----------------------------------------------------------------------------
# 1. PACKAGES
# -----------------------------------------------------------------------------

required <- c(
  "readr",
  "dplyr",
  "tidyr",
  "ggplot2",
  "scales",
  "stringr",
  "purrr",
  "flextable",
  "officer",
  "tibble"
)

missing <- required[
  !vapply(required, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing)) {
  stop(
    "Install the following packages first: ",
    paste(missing, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(stringr)
  library(purrr)
  library(flextable)
  library(officer)
  library(tibble)
})


# -----------------------------------------------------------------------------
# 2. HELPERS
# -----------------------------------------------------------------------------

locate_input <- function(filename) {
  candidates <- file.path(SEARCH_DIRS, filename)
  found <- candidates[file.exists(candidates)]

  if (!length(found)) {
    stop(
      "\nRequired input file not found: ",
      filename,
      "\nSearched:\n  ",
      paste(SEARCH_DIRS, collapse = "\n  "),
      "\n"
    )
  }

  normalizePath(found[1], winslash = "/", mustWork = TRUE)
}


pct <- function(n, d, digits = 1) {
  if (is.na(n) || is.na(d) || d == 0) {
    return("")
  }
  sprintf(paste0("%.", digits, "f"), 100 * n / d)
}


n_pct_text <- function(n, d, digits = 1) {
  if (is.na(n) || is.na(d) || d == 0) {
    return("")
  }

  paste0(
    n,
    "/",
    d,
    " (",
    pct(n, d, digits = digits),
    "%)"
  )
}


fmt_num <- function(x, digits = 3) {
  ifelse(
    is.na(x),
    "",
    formatC(x, format = "f", digits = digits)
  )
}


fmt_p <- function(x) {
  ifelse(
    is.na(x),
    "",
    ifelse(
      x < 0.001,
      "<0.001",
      formatC(x, format = "f", digits = 3)
    )
  )
}


method_labels <- c(
  "Fisher" = "Fisher exact",
  "MidP" = "Mid-P",
  "Delta" = expression(Delta ~ "Wald"),
  "Conditional_LR" = "Conditional LR",
  "BACO_Wald" = "Original BACO Wald"
)


method_labels_text <- c(
  "Fisher" = "Fisher exact",
  "MidP" = "Mid-P",
  "Delta" = "Delta Wald",
  "Conditional_LR" = "Conditional LR",
  "BACO_Wald" = "Original BACO Wald"
)


method_order <- c(
  "Fisher",
  "MidP",
  "Delta",
  "Conditional_LR",
  "BACO_Wald"
)


publication_theme <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = base_size + 1,
        margin = margin(b = 6)
      ),
      plot.subtitle = element_text(
        size = base_size - 1,
        margin = margin(b = 8)
      ),
      plot.caption = element_text(
        size = base_size - 2,
        hjust = 0,
        margin = margin(t = 8)
      ),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      legend.box = "vertical"
    )
}


style_ft <- function(ft, font_size = 9) {
  ft <- flextable::theme_booktabs(ft)
  ft <- flextable::font(
    ft,
    fontname = "Arial",
    part = "all"
  )
  ft <- flextable::fontsize(
    ft,
    size = font_size,
    part = "all"
  )
  ft <- flextable::bold(
    ft,
    part = "header"
  )
  ft <- flextable::valign(
    ft,
    valign = "center",
    part = "all"
  )
  ft <- flextable::autofit(ft)

  ft
}


save_ft <- function(ft, path) {
  flextable::save_as_docx(ft, path = path)
  invisible(path)
}


save_plot_both <- function(plot, basename, directory, width, height) {
  pdf_path <- file.path(directory, paste0(basename, ".pdf"))
  png_path <- file.path(directory, paste0(basename, ".png"))

  ggsave(
    filename = pdf_path,
    plot = plot,
    width = width,
    height = height,
    units = "in"
  )

  ggsave(
    filename = png_path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600
  )

  invisible(c(pdf_path, png_path))
}


# -----------------------------------------------------------------------------
# 3. INPUT FILES
# -----------------------------------------------------------------------------

null_file <- locate_input("02_null_calibration_final_sensitivity.csv")
power_file <- locate_input("03_power_final_delta_sensitivity.csv")
empirical_file <- locate_input("empirical_baco_full.csv")
zero_event_file <- locate_input("01_zero_event_test_simulation_results.csv")
replication_file <- locate_input("external_replication_baco_full.csv")

null_dat <- read_csv(null_file, show_col_types = FALSE)
power_dat <- read_csv(power_file, show_col_types = FALSE)
emp <- read_csv(empirical_file, show_col_types = FALSE)
zero_dat <- read_csv(zero_event_file, show_col_types = FALSE)
rep <- read_csv(replication_file, show_col_types = FALSE)

## 3.1 Standardized study identifiers

# The public supplementary workbook already contains a standardized first_author
# field. No unpublished author strings or article-level descriptive fields are
# required to reproduce the manuscript tables and figures.
stopifnot("first_author" %in% names(emp), "first_author" %in% names(rep))


# -----------------------------------------------------------------------------
# 3.2 External replication analysis sets
# -----------------------------------------------------------------------------

rep_primary <- rep %>%
  filter(primary_analysis %in% TRUE)

rep_direct <- rep %>%
  filter(direct_count_sensitivity %in% TRUE)

# The standardized workbook stores the adverse-complement counts for the two
# favorable-composite comparisons, so all 24 rows form the reoriented sensitivity set.
rep_all_reoriented <- rep


# -----------------------------------------------------------------------------
# 4. FREEZE BASIC DATA CHECKS
# -----------------------------------------------------------------------------

stopifnot(
  nrow(emp) == 130,
  dplyr::n_distinct(emp$study_id) == 123,
  sum(emp$aggregate_count_valid, na.rm = TRUE) == 129,
  sum(emp$fisher_available, na.rm = TRUE) == 129,
  sum(emp$baco_estimable, na.rm = TRUE) == 111,
  sum(emp$delta_estimable, na.rm = TRUE) == 114,
  sum(emp$fisher_sig_05, na.rm = TRUE) == 5,
  sum(emp$delta_wald_sig_05, na.rm = TRUE) == 8,
  sum(emp$original_baco_wald_sig_05, na.rm = TRUE) == 27,
  sum(emp$fisher_sig_005, na.rm = TRUE) == 1,
  sum(emp$original_baco_wald_sig_005, na.rm = TRUE) == 11
)


# Freeze the external replication dataset and its prespecified analysis sets.
stopifnot(
  nrow(rep) == 24,
  dplyr::n_distinct(rep$study_id) == 23,
  sum(rep$primary_analysis, na.rm = TRUE) == 22,
  sum(rep$direct_count_sensitivity, na.rm = TRUE) == 20,
  sum(rep$reoriented_for_sensitivity, na.rm = TRUE) == 2,
  sum(rep$data_origin == "reconstructed", na.rm = TRUE) == 2,
  length(intersect(tolower(emp$doi), tolower(rep$doi))) == 0,

  nrow(rep_primary) == 22,
  dplyr::n_distinct(rep_primary$study_id) == 21,
  sum(rep_primary$aggregate_count_valid, na.rm = TRUE) == 22,
  sum(rep_primary$zero_critical_any_arm, na.rm = TRUE) == 0,
  sum(rep_primary$zero_critical_both_arms, na.rm = TRUE) == 0,
  sum(rep_primary$fieller_available, na.rm = TRUE) == 22,
  sum(rep_primary$weak_denominator, na.rm = TRUE) == 21,
  sum(rep_primary$fieller_type == "bounded", na.rm = TRUE) == 1,
  sum(rep_primary$fieller_type == "disjoint_unbounded", na.rm = TRUE) == 5,
  sum(rep_primary$fieller_type == "all_real", na.rm = TRUE) == 16,
  sum(rep_primary$fisher_sig_05, na.rm = TRUE) == 2,
  sum(rep_primary$delta_wald_sig_05, na.rm = TRUE) == 2,
  sum(rep_primary$original_baco_wald_sig_05, na.rm = TRUE) == 4,

  nrow(rep_direct) == 20,
  sum(rep_direct$weak_denominator, na.rm = TRUE) == 19,
  sum(rep_direct$fisher_sig_05, na.rm = TRUE) == 1,
  sum(rep_direct$delta_wald_sig_05, na.rm = TRUE) == 1,
  sum(rep_direct$original_baco_wald_sig_05, na.rm = TRUE) == 3,

  nrow(rep_all_reoriented) == 24,
  sum(rep_all_reoriented$weak_denominator, na.rm = TRUE) == 23,
  sum(rep_all_reoriented$fieller_type == "bounded", na.rm = TRUE) == 1,
  sum(rep_all_reoriented$fieller_type == "disjoint_unbounded", na.rm = TRUE) ==
    5,
  sum(rep_all_reoriented$fieller_type == "all_real", na.rm = TRUE) == 18
)

# =============================================================================
# MAIN MANUSCRIPT TABLES
# =============================================================================

# -----------------------------------------------------------------------------
# TABLE 1
# Operating characteristics under the null hypothesis
# -----------------------------------------------------------------------------

table1_long <- null_dat %>%
  mutate(
    method = factor(method, levels = method_order),
    denominator_class = factor(
      denominator_class,
      levels = c("strong", "weak")
    )
  ) %>%
  group_by(method, denominator_class) %>%
  summarise(
    scenarios = n(),
    mean_type1 = mean(nominal_type1, na.rm = TRUE),
    max_type1 = max(nominal_type1, na.rm = TRUE),
    median_estimable = median(null_estimable_fraction, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    method = as.character(method),
    Method = unname(method_labels_text[method])
  )


table1 <- table1_long %>%
  select(
    Method,
    denominator_class,
    scenarios,
    mean_type1,
    max_type1,
    median_estimable
  ) %>%
  pivot_wider(
    names_from = denominator_class,
    values_from = c(
      scenarios,
      mean_type1,
      max_type1,
      median_estimable
    ),
    names_glue = "{denominator_class}_{.value}"
  ) %>%
  mutate(
    `Strong: mean type I error` = sprintf("%.1f%%", 100 * strong_mean_type1),
    `Strong: maximum type I error` = sprintf("%.1f%%", 100 * strong_max_type1),
    `Strong: median estimable` = sprintf(
      "%.1f%%",
      100 * strong_median_estimable
    ),
    `Weak: mean type I error` = sprintf("%.1f%%", 100 * weak_mean_type1),
    `Weak: maximum type I error` = sprintf("%.1f%%", 100 * weak_max_type1),
    `Weak: median estimable` = sprintf("%.1f%%", 100 * weak_median_estimable)
  ) %>%
  select(
    Method,
    `Strong: mean type I error`,
    `Strong: maximum type I error`,
    `Strong: median estimable`,
    `Weak: mean type I error`,
    `Weak: maximum type I error`,
    `Weak: median estimable`
  )


ft1 <- flextable(table1)

ft1 <- set_header_labels(
  ft1,
  Method = "Method",
  `Strong: mean type I error` = "Mean",
  `Strong: maximum type I error` = "Maximum",
  `Strong: median estimable` = "Estimable",
  `Weak: mean type I error` = "Mean",
  `Weak: maximum type I error` = "Maximum",
  `Weak: median estimable` = "Estimable"
)

ft1 <- add_header_row(
  ft1,
  values = c(
    "",
    "Strong critical-endpoint effect",
    "Weak critical-endpoint effect"
  ),
  colwidths = c(1, 3, 3)
)

ft1 <- style_ft(ft1, font_size = 8.5)

ft1 <- align(
  ft1,
  j = 2:7,
  align = "center",
  part = "all"
)

ft1 <- set_caption(
  ft1,
  caption = paste0(
    "Table 1. Operating characteristics of the inferential procedures ",
    "under the null hypothesis."
  )
)

ft1 <- add_footer_lines(
  ft1,
  values = paste0(
    "Values summarize the final Delta-parameterized null simulation grid. ",
    "Mean and maximum empirical type I error are shown at nominal alpha = 0.05. ",
    "Estimable is the median proportion of simulation replicates for which the ",
    "corresponding statistic could be calculated. Fisher exact is the primary ",
    "count-based test; the original BACO Wald procedure is shown as a historical comparator."
  )
)

save_ft(
  ft1,
  file.path(MAIN_DIR, "Table_1_operating_characteristics.docx")
)

write_csv(
  table1_long,
  file.path(DATA_DIR, "Table_1_underlying_data.csv")
)


# -----------------------------------------------------------------------------
# TABLE 2
# Detailed empirical comparison: primary reanalysis vs external replication
# -----------------------------------------------------------------------------

# The two empirical datasets are presented at the same level of detail in a
# single main-manuscript table. They are not pooled. The primary dataset is the
# 130-comparison empirical reanalysis, whereas the external dataset is the
# independently assembled COVID-19 review. The primary external replication is
# restricted to the 22 comparisons satisfying the nested adverse-endpoint
# structure as originally reported.

# Primary empirical reanalysis -------------------------------------------------

n_all <- nrow(emp)
n_studies <- n_distinct(emp$study_id)
n_published <- sum(emp$source_group == "published_benchmark", na.rm = TRUE)
n_published_studies <- n_distinct(emp$study_id[
  emp$source_group == "published_benchmark"
])
n_additional <- sum(emp$source_group == "additional_restored", na.rm = TRUE)
n_additional_studies <- n_distinct(emp$study_id[
  emp$source_group == "additional_restored"
])

n_valid <- sum(emp$aggregate_count_valid, na.rm = TRUE)
n_valid_studies <- n_distinct(
  emp$study_id[emp$aggregate_count_valid %in% TRUE]
)

n_fisher <- sum(emp$fisher_available, na.rm = TRUE)
n_baco <- sum(emp$baco_estimable, na.rm = TRUE)
n_delta <- sum(emp$delta_estimable, na.rm = TRUE)
n_legacy <- sum(emp$original_baco_wald_available, na.rm = TRUE)

n_zero_any <- sum(emp$zero_critical_any_arm %in% TRUE, na.rm = TRUE)
n_zero_both <- sum(emp$zero_critical_both_arms %in% TRUE, na.rm = TRUE)

n_fieller <- sum(emp$fieller_available, na.rm = TRUE)
n_bounded <- sum(emp$fieller_type == "bounded", na.rm = TRUE)
n_disjoint <- sum(emp$fieller_type == "disjoint_unbounded", na.rm = TRUE)
n_all_real <- sum(emp$fieller_type == "all_real", na.rm = TRUE)
n_weak <- sum(emp$weak_denominator %in% TRUE, na.rm = TRUE)

n_comparable <- sum(
  emp$method_agreement_05 != "not_comparable",
  na.rm = TRUE
)

# External COVID-19 replication -----------------------------------------------

n_rep_all <- nrow(rep)
n_rep_all_studies <- n_distinct(rep$study_id)
n_rep_primary <- nrow(rep_primary)
n_rep_primary_studies <- n_distinct(rep_primary$study_id)
n_rep_valid <- sum(rep_primary$aggregate_count_valid, na.rm = TRUE)
n_rep_valid_studies <- n_distinct(
  rep_primary$study_id[rep_primary$aggregate_count_valid %in% TRUE]
)

n_rep_fisher <- sum(rep_primary$fisher_available, na.rm = TRUE)
n_rep_baco <- sum(rep_primary$baco_estimable, na.rm = TRUE)
n_rep_delta <- sum(rep_primary$delta_estimable, na.rm = TRUE)
n_rep_legacy <- sum(rep_primary$original_baco_wald_available, na.rm = TRUE)

n_rep_zero_any <- sum(
  rep_primary$zero_critical_any_arm %in% TRUE,
  na.rm = TRUE
)
n_rep_zero_both <- sum(
  rep_primary$zero_critical_both_arms %in% TRUE,
  na.rm = TRUE
)

n_rep_fieller <- sum(rep_primary$fieller_available, na.rm = TRUE)
n_rep_bounded <- sum(
  rep_primary$fieller_type == "bounded",
  na.rm = TRUE
)
n_rep_disjoint <- sum(
  rep_primary$fieller_type == "disjoint_unbounded",
  na.rm = TRUE
)
n_rep_all_real <- sum(
  rep_primary$fieller_type == "all_real",
  na.rm = TRUE
)
n_rep_weak <- sum(rep_primary$weak_denominator %in% TRUE, na.rm = TRUE)

n_rep_comparable <- sum(
  rep_primary$method_agreement_05 != "not_comparable",
  na.rm = TRUE
)

# One detailed table for both empirical datasets ------------------------------

table2 <- tribble(
  ~Section                                                                                                  ,
  ~Characteristic                                                                                           ,
  ~`Primary empirical reanalysis`                                                                           ,
  ~`External COVID-19 replication`                                                                          ,

  "Dataset description"                                                                                     ,
  "Studies represented in full dataset"                                                                     ,
  as.character(n_studies)                                                                                   ,
  as.character(n_rep_all_studies)                                                                           ,

  "Dataset description"                                                                                     ,
  "Comparisons identified in full dataset"                                                                  ,
  as.character(n_all)                                                                                       ,
  as.character(n_rep_all)                                                                                   ,

  "Dataset description"                                                                                     ,
  "Studies represented in primary count-based analysis"                                                     ,
  as.character(n_valid_studies)                                                                             ,
  as.character(n_rep_valid_studies)                                                                         ,

  "Dataset description"                                                                                     ,
  "Comparisons included in primary count-based analysis"                                                    ,
  n_pct_text(n_valid, n_all)                                                                                ,
  n_pct_text(n_rep_valid, n_rep_all)                                                                        ,

  "Sparse critical events"                                                                                  ,
  "At least one zero critical-event cell"                                                                   ,
  n_pct_text(n_zero_any, n_valid)                                                                           ,
  n_pct_text(n_rep_zero_any, n_rep_valid)                                                                   ,

  "Sparse critical events"                                                                                  ,
  "Zero critical events in both arms"                                                                       ,
  n_pct_text(n_zero_both, n_valid)                                                                          ,
  n_pct_text(n_rep_zero_both, n_rep_valid)                                                                  ,

  "Availability and estimability"                                                                           ,
  "Fisher exact test available"                                                                             ,
  n_pct_text(n_fisher, n_valid)                                                                             ,
  n_pct_text(n_rep_fisher, n_rep_valid)                                                                     ,

  "Availability and estimability"                                                                           ,
  "BACO point estimate available"                                                                           ,
  n_pct_text(n_baco, n_valid)                                                                               ,
  n_pct_text(n_rep_baco, n_rep_valid)                                                                       ,

  "Availability and estimability"                                                                           ,
  "Delta available"                                                                                         ,
  n_pct_text(n_delta, n_valid)                                                                              ,
  n_pct_text(n_rep_delta, n_rep_valid)                                                                      ,

  "Availability and estimability"                                                                           ,
  "Fieller confidence set available"                                                                        ,
  n_pct_text(n_fieller, n_valid)                                                                            ,
  n_pct_text(n_rep_fieller, n_rep_valid)                                                                    ,

  "Fieller inference"                                                                                       ,
  "Weak denominator"                                                                                        ,
  n_pct_text(n_weak, n_fieller)                                                                             ,
  n_pct_text(n_rep_weak, n_rep_fieller)                                                                     ,

  "Fieller inference"                                                                                       ,
  "Bounded Fieller confidence interval"                                                                     ,
  n_pct_text(n_bounded, n_fieller)                                                                          ,
  n_pct_text(n_rep_bounded, n_rep_fieller)                                                                  ,

  "Fieller inference"                                                                                       ,
  "Disjoint unbounded Fieller set"                                                                          ,
  n_pct_text(n_disjoint, n_fieller)                                                                         ,
  n_pct_text(n_rep_disjoint, n_rep_fieller)                                                                 ,

  "Fieller inference"                                                                                       ,
  "Fieller set covering the entire real line"                                                               ,
  n_pct_text(n_all_real, n_fieller)                                                                         ,
  n_pct_text(n_rep_all_real, n_rep_fieller)                                                                 ,

  "Decisions at alpha = 0.05"                                                                               ,
  "Fisher exact rejections"                                                                                 ,
  n_pct_text(sum(emp$fisher_sig_05, na.rm = TRUE), n_fisher)                                                ,
  n_pct_text(sum(rep_primary$fisher_sig_05, na.rm = TRUE), n_rep_fisher)                                    ,

  "Decisions at alpha = 0.05"                                                                               ,
  "Delta Wald rejections"                                                                                   ,
  n_pct_text(sum(emp$delta_wald_sig_05, na.rm = TRUE), n_delta)                                             ,
  n_pct_text(sum(rep_primary$delta_wald_sig_05, na.rm = TRUE), n_rep_delta)                                 ,

  "Decisions at alpha = 0.05"                                                                               ,
  "Original BACO Wald rejections"                                                                           ,
  n_pct_text(sum(emp$original_baco_wald_sig_05, na.rm = TRUE), n_legacy)                                    ,
  n_pct_text(sum(rep_primary$original_baco_wald_sig_05, na.rm = TRUE), n_rep_legacy)                        ,

  "Wald-Fisher agreement at alpha = 0.05"                                                                   ,
  "Both significant"                                                                                        ,
  n_pct_text(sum(emp$method_agreement_05 == "both_significant", na.rm = TRUE), n_comparable)                ,
  n_pct_text(sum(rep_primary$method_agreement_05 == "both_significant", na.rm = TRUE), n_rep_comparable)    ,

  "Wald-Fisher agreement at alpha = 0.05"                                                                   ,
  "Both nonsignificant"                                                                                     ,
  n_pct_text(sum(emp$method_agreement_05 == "both_nonsignificant", na.rm = TRUE), n_comparable)             ,
  n_pct_text(sum(rep_primary$method_agreement_05 == "both_nonsignificant", na.rm = TRUE), n_rep_comparable) ,

  "Wald-Fisher agreement at alpha = 0.05"                                                                   ,
  "Fisher only"                                                                                             ,
  n_pct_text(sum(emp$method_agreement_05 == "fisher_only", na.rm = TRUE), n_comparable)                     ,
  n_pct_text(sum(rep_primary$method_agreement_05 == "fisher_only", na.rm = TRUE), n_rep_comparable)         ,

  "Wald-Fisher agreement at alpha = 0.05"                                                                   ,
  "Original BACO Wald only"                                                                                 ,
  n_pct_text(sum(emp$method_agreement_05 == "original_wald_only", na.rm = TRUE), n_comparable)              ,
  n_pct_text(sum(rep_primary$method_agreement_05 == "original_wald_only", na.rm = TRUE), n_rep_comparable)  ,

  "Decisions at alpha = 0.005"                                                                              ,
  "Fisher exact rejections"                                                                                 ,
  n_pct_text(sum(emp$fisher_sig_005, na.rm = TRUE), n_fisher)                                               ,
  n_pct_text(sum(rep_primary$fisher_sig_005, na.rm = TRUE), n_rep_fisher)                                   ,

  "Decisions at alpha = 0.005"                                                                              ,
  "Original BACO Wald rejections"                                                                           ,
  n_pct_text(sum(emp$original_baco_wald_sig_005, na.rm = TRUE), n_legacy)                                   ,
  n_pct_text(sum(rep_primary$original_baco_wald_sig_005, na.rm = TRUE), n_rep_legacy)
)

ft2 <- flextable(table2)

ft2 <- merge_v(
  ft2,
  j = "Section"
)

ft2 <- valign(
  ft2,
  j = "Section",
  valign = "top",
  part = "body"
)

ft2 <- style_ft(ft2, font_size = 7.8)
ft2 <- width(ft2, j = "Section", width = 1.45)
ft2 <- width(ft2, j = "Characteristic", width = 3.15)
ft2 <- width(ft2, j = "Primary empirical reanalysis", width = 1.55)
ft2 <- width(ft2, j = "External COVID-19 replication", width = 1.75)

ft2 <- align(
  ft2,
  j = 3:4,
  align = "center",
  part = "all"
)

ft2 <- set_caption(
  ft2,
  caption = paste0(
    "Table 2. Empirical BACO inference in the primary reanalysis and external ",
    "COVID-19 replication."
  )
)

ft2 <- add_footer_lines(
  ft2,
  values = paste0(
    "Values are n/N (%) when a method-specific denominator applies. The primary ",
    "dataset comprised 97 published-benchmark comparisons from 91 studies and 33 ",
    "additional comparisons from 32 studies; 129 of 130 comparisons had valid ",
    "aggregate counts. The external review identified 24 comparisons from 23 studies; ",
    "the primary external replication included 22 comparisons from 21 studies ",
    "satisfying the nested adverse-endpoint structure. The two favorable-composite ",
    "comparisons were excluded from the primary replication and included only after ",
    "explicit adverse reorientation in sensitivity analysis. Dequin and Salama were ",
    "retained using transparently reconstructed aligned counts; the direct-count-only ",
    "sensitivity analysis is reported in Table S8. The datasets were analyzed ",
    "separately and no studies overlapped. Wald-Fisher agreement is restricted to ",
    "comparisons for which both procedures were available."
  )
)

save_ft(
  ft2,
  file.path(
    MAIN_DIR,
    "Table_2_empirical_reanalysis_external_replication.docx"
  )
)

write_csv(
  table2,
  file.path(DATA_DIR, "Table_2_underlying_data.csv")
)


# =============================================================================
# MAIN MANUSCRIPT FIGURES
# =============================================================================

# -----------------------------------------------------------------------------
# FIGURE 1
# Type I error across expected critical-event information
# -----------------------------------------------------------------------------

fig1_data <- null_dat %>%
  filter(
    method %in% c("Fisher", "Delta", "BACO_Wald")
  ) %>%
  mutate(
    Method = factor(
      unname(method_labels_text[method]),
      levels = c(
        "Fisher exact",
        "Delta Wald",
        "Original BACO Wald"
      )
    ),
    `Critical-endpoint effect` = factor(
      denominator_class,
      levels = c("strong", "weak"),
      labels = c(
        "Strong critical-endpoint effect",
        "Weak critical-endpoint effect"
      )
    )
  ) %>%
  group_by(
    `Critical-endpoint effect`,
    Method,
    expected_total_deaths_target
  ) %>%
  summarise(
    mean_type1 = mean(nominal_type1, na.rm = TRUE),
    min_type1 = min(nominal_type1, na.rm = TRUE),
    max_type1 = max(nominal_type1, na.rm = TRUE),
    .groups = "drop"
  )


fig1 <- ggplot(
  fig1_data,
  aes(
    x = expected_total_deaths_target,
    y = mean_type1,
    color = Method,
    linetype = Method,
    shape = Method,
    group = Method
  )
) +
  geom_hline(
    yintercept = 0.05,
    linetype = 3,
    linewidth = 0.55,
    color = "black"
  ) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.4) +
  facet_wrap(
    ~`Critical-endpoint effect`,
    nrow = 1
  ) +
  scale_x_continuous(
    breaks = sort(unique(fig1_data$expected_total_deaths_target))
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  scale_color_manual(
    values = c(
      "Fisher exact" = "#0072B2",
      "Delta Wald" = "#D55E00",
      "Original BACO Wald" = "#009E73"
    )
  ) +
  labs(
    x = "Expected total critical events",
    y = "Empirical type I error",
    color = "Method",
    linetype = "Method",
    shape = "Method"
  ) +
  publication_theme()

save_plot_both(
  fig1,
  "Figure_1_type1_error",
  MAIN_DIR,
  width = 8.2,
  height = 4.8
)

write_csv(
  fig1_data,
  file.path(DATA_DIR, "Figure_1_data.csv")
)

writeLines(
  paste(
    "Fig. 1. Empirical type I error across critical-event information and denominator strength.",
    "Note: The dashed horizontal line indicates the nominal α = 0.05 level.",
    "Points represent average empirical type I error across the final simulation design.",
    "Fisher's exact test remained conservative, whereas the original BACO Wald procedure",
    "showed substantial type I error inflation when the critical-endpoint effect was weak."
  ),
  file.path(
    MAIN_DIR,
    "Fig_1_caption.txt"
  )
)


# -----------------------------------------------------------------------------
# FIGURE 2
# Size-adjusted power parameterized directly by |Delta|
# -----------------------------------------------------------------------------

fig2_data <- power_dat %>%
  filter(method == "Fisher") %>%
  mutate(
    `Critical-endpoint effect` = factor(
      denominator_class,
      levels = c("strong", "weak"),
      labels = c(
        "Strong critical-endpoint effect",
        "Weak critical-endpoint effect"
      )
    ),
    `Expected critical events` = factor(
      expected_total_deaths_target,
      levels = sort(unique(expected_total_deaths_target))
    )
  ) %>%
  group_by(
    `Critical-endpoint effect`,
    `Expected critical events`,
    abs_population_Delta
  ) %>%
  summarise(
    mean_power = mean(size_adjusted_power, na.rm = TRUE),
    .groups = "drop"
  )


fig2 <- ggplot(
  fig2_data,
  aes(
    x = abs_population_Delta,
    y = mean_power,
    color = `Expected critical events`,
    linetype = `Expected critical events`,
    shape = `Expected critical events`,
    group = `Expected critical events`
  )
) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.4) +
  facet_wrap(
    ~`Critical-endpoint effect`,
    nrow = 1
  ) +
  scale_x_continuous(
    breaks = sort(unique(fig2_data$abs_population_Delta))
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.08))
  ) +
  scale_color_manual(
    values = c(
      "2" = "#0072B2",
      "5" = "#D55E00",
      "10" = "#009E73",
      "20" = "#984ea3"
    )
  ) +
  labs(
    x = expression(abs(Delta)),
    y = "Size-adjusted power",
    linetype = "Expected critical events",
    shape = "Expected critical events"
  ) +
  publication_theme()


save_plot_both(
  fig2,
  "Figure_2_size_adjusted_power_Delta",
  MAIN_DIR,
  width = 8.2,
  height = 4.8
)

write_csv(
  fig2_data,
  file.path(DATA_DIR, "Figure_2_data.csv")
)

writeLines(
  paste(
    "Fig. 2. Size-adjusted power parameterized by the log-effect discrepancy.",
    "Note: Fisher's exact test. Alternatives are defined by |Delta| =",
    "|log(RR_composite) - log(RR_critical)|.",
    "For a fixed discrepancy magnitude, power was similar under strong and",
    "weak critical-endpoint effects and increased primarily with critical-event information."
  ),
  file.path(
    MAIN_DIR,
    "Fig_2_caption.txt"
  )
)


# -----------------------------------------------------------------------------
# FIGURE 3
# Original BACO Wald versus Fisher exact in both empirical datasets
# -----------------------------------------------------------------------------

# The main empirical reanalysis and the external COVID-19 replication are shown
# in separate panels using identical axes and decision categories. This allows a
# direct visual assessment of whether the Wald-Fisher discordance pattern
# observed in the primary dataset is reproduced in the independent dataset.

fig3_primary <- emp %>%
  filter(
    fisher_available,
    original_baco_wald_available,
    is.finite(fisher_p),
    is.finite(original_baco_wald_p)
  ) %>%
  transmute(
    comparison_id,
    Dataset = "A. Primary empirical reanalysis",
    fisher_p,
    original_baco_wald_p,
    method_agreement_05
  )

fig3_external <- rep_primary %>%
  filter(
    fisher_available,
    original_baco_wald_available,
    is.finite(fisher_p),
    is.finite(original_baco_wald_p)
  ) %>%
  transmute(
    comparison_id,
    Dataset = "B. External COVID-19 replication",
    fisher_p,
    original_baco_wald_p,
    method_agreement_05
  )

fig3_data <- bind_rows(
  fig3_primary,
  fig3_external
) %>%
  mutate(
    Dataset = factor(
      Dataset,
      levels = c(
        "A. Primary empirical reanalysis",
        "B. External COVID-19 replication"
      )
    ),
    fisher_plot = pmax(fisher_p, 1e-16),
    original_baco_wald_plot = pmax(original_baco_wald_p, 1e-16),
    x = -log10(original_baco_wald_plot),
    y = -log10(fisher_plot),
    Decision = factor(
      method_agreement_05,
      levels = c(
        "both_nonsignificant",
        "both_significant",
        "fisher_only",
        "original_wald_only"
      ),
      labels = c(
        "Both nonsignificant",
        "Both significant",
        "Fisher only",
        "Original BACO Wald only"
      )
    )
  )


cut05 <- -log10(0.05)
cut005 <- -log10(0.005)

# Common coordinates are used for both panels so that the geometry is directly
# comparable. Left-side labels are left-aligned to prevent clipping at the
# panel boundary.

max_x <- max(fig3_data$x, na.rm = TRUE)
max_y <- max(fig3_data$y, na.rm = TRUE)

# Position labels further inside the plotting area.
# hjust = 0 will make left-side labels extend toward the right.
x_left <- cut05 * 0.20
x_right <- cut05 + 0.68 * (max_x - cut05)

y_low <- cut05 * 0.85
y_high <- cut05 + 0.68 * (max_y - cut05)


fig3_annotations <- fig3_data %>%
  count(Dataset, Decision, name = "n", .drop = FALSE) %>%
  tidyr::complete(
    Dataset,
    Decision,
    fill = list(n = 0)
  ) %>%
  mutate(
    x = case_when(
      Decision %in%
        c(
          "Both nonsignificant",
          "Fisher only"
        ) ~ x_left,
      TRUE ~ x_right
    ),

    y = case_when(
      Decision %in%
        c(
          "Both significant",
          "Fisher only"
        ) ~ y_high,
      TRUE ~ y_low
    ),

    # Left quadrants are left-aligned so the text grows into the panel.
    hjust = case_when(
      Decision %in%
        c(
          "Both nonsignificant",
          "Fisher only"
        ) ~ 0,
      TRUE ~ 0.5
    ),

    label = case_when(
      Decision == "Both nonsignificant" ~ paste0(
        "Both nonsignificant\nn = ",
        n
      ),
      Decision == "Both significant" ~ paste0(
        "Both significant\nn = ",
        n
      ),
      Decision == "Fisher only" ~ paste0(
        "Fisher only\nn = ",
        n
      ),
      Decision == "Original BACO Wald only" ~ paste0(
        "Original Wald only\nn = ",
        n
      ),
      TRUE ~ ""
    )
  )


fig3 <- ggplot(
  fig3_data,
  aes(
    x = x,
    y = y,
    color = Decision,
    shape = Decision
  )
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = 3,
    linewidth = 0.45
  ) +
  geom_vline(
    xintercept = cut05,
    linetype = 2,
    linewidth = 0.55
  ) +
  geom_hline(
    yintercept = cut05,
    linetype = 2,
    linewidth = 0.55
  ) +
  geom_vline(
    xintercept = cut005,
    linetype = 3,
    linewidth = 0.45
  ) +
  geom_hline(
    yintercept = cut005,
    linetype = 3,
    linewidth = 0.45
  ) +
  geom_point(
    size = 2.2,
    alpha = 0.78
  ) +
  geom_text(
    data = fig3_annotations,
    aes(
      x = x,
      y = y,
      label = label,
      hjust = hjust
    ),
    inherit.aes = FALSE,
    size = 2.8,
    lineheight = 0.95
  ) +
  facet_wrap(
    ~Dataset,
    nrow = 1,
    scales = "fixed"
  ) +
  scale_color_manual(
    name = "Decision at alpha = 0.05",
    values = c(
      "Both nonsignificant" = "#D55E00",
      "Both significant" = "#984ea3",
      "Original BACO Wald only" = "#0072B2",
      "Fisher only" = "#009E73"
    )
  ) +
  scale_shape_manual(
    name = "Decision at alpha = 0.05",
    values = c(
      "Both nonsignificant" = 16,
      "Both significant" = 17,
      "Original BACO Wald only" = 15,
      "Fisher only" = 18
    )
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.03, 0.05))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(
    x = expression(-log[10]("original BACO Wald p")),
    y = expression(-log[10]("Fisher exact p")),
    color = "Decision at alpha = 0.05",
    shape = "Decision at alpha = 0.05"
  ) +
  publication_theme(base_size = 10) +
  theme(
    strip.text = element_text(
      face = "bold",
      size = 10
    ),
    panel.spacing.x = unit(1.0, "lines")
  )


save_plot_both(
  fig3,
  "Figure_3_Wald_vs_Fisher",
  MAIN_DIR,
  width = 8.2,
  height = 5.2
)

write_csv(
  fig3_data,
  file.path(
    DATA_DIR,
    "Figure_3_data.csv"
  )
)

write_csv(
  fig3_annotations,
  file.path(
    DATA_DIR,
    "Figure_3_quadrant_counts.csv"
  )
)

writeLines(
  paste(
    "Fig. 3. Original BACO Wald test versus Fisher's exact test in the two empirical datasets.",
    "Panel A shows the primary empirical reanalysis and panel B the independent COVID-19 replication.",
    "Dashed horizontal and vertical lines denote p = 0.05, dotted horizontal and vertical lines denote p = 0.005,",
    "and the diagonal dotted line denotes equality of the two p-values.",
    "Quadrant counts refer to decisions at alpha = 0.05. Only comparisons for which both procedures were available are displayed.",
    "The same axes are used in both panels to facilitate direct visual comparison."
  ),
  file.path(
    MAIN_DIR,
    "Fig_3_caption.txt"
  )
)


# =============================================================================
# SUPPLEMENTARY TABLES
# =============================================================================

# -----------------------------------------------------------------------------
# TABLE S1
# Simulation design
# -----------------------------------------------------------------------------

table_s1 <- tribble(
  ~Component                                           , ~Values                                                ,
  "Total sample size"                                  , "500; 2,000; 10,000"                                   ,
  "Allocation ratio"                                   , "1:1; 2:1"                                             ,
  "Composite-event risk in control"                    , "0.05; 0.25; 0.50"                                     ,
  "Expected total critical events"                     , "2; 5; 10; 20"                                         ,
  "Critical-endpoint RR"                               , "0.50; 0.95; 1.05; 2.00"                               ,
  "Denominator class"                                  , "Strong: RRd = 0.50 or 2.00; Weak: RRd = 0.95 or 1.05" ,
  "Delta alternatives"                                 , "-0.50; -0.25; -0.10; 0.10; 0.25; 0.50"                ,
  "Null simulation replicates"                         , "5,000 per scenario"                                   ,
  "Alternative simulation replicates"                  , "3,000 per scenario"                                   ,
  "Nominal alpha"                                      , "0.05"                                                 ,
  "Legacy alpha used in the original BACO application" , "0.005"                                                ,
  "Primary exact test"                                 , "Probability-ordered Fisher exact test"                ,
  "Ratio confidence procedure"                         , "Fieller confidence set"                               ,
  "Sensitivity procedures"                             , "Mid-P; conditional LR; Delta Wald"                    ,
  "Historical comparator"                              , "Original BACO Wald/Taylor procedure"
)


fts1 <- flextable(table_s1)
fts1 <- style_ft(fts1, font_size = 9)
fts1 <- width(fts1, j = "Component", width = 2.7)
fts1 <- width(fts1, j = "Values", width = 4.8)

fts1 <- set_caption(
  fts1,
  caption = "Table S1. Final simulation design and inferential procedures."
)

save_ft(
  fts1,
  file.path(SUPP_DIR, "Table_S1_simulation_design.docx")
)


# -----------------------------------------------------------------------------
# TABLE S2
# Type I error summary across final null grid
# -----------------------------------------------------------------------------

table_s2 <- null_dat %>%
  mutate(
    Method = unname(method_labels_text[method]),
    `Denominator class` = ifelse(
      denominator_class == "strong",
      "Strong",
      "Weak"
    )
  ) %>%
  group_by(
    Method,
    allocation_label,
    `Denominator class`
  ) %>%
  summarise(
    Scenarios = n(),
    `Mean type I error` = mean(nominal_type1, na.rm = TRUE),
    `Maximum type I error` = max(nominal_type1, na.rm = TRUE),
    `Mean calibrated size` = mean(calibrated_size, na.rm = TRUE),
    `Median estimable fraction` = median(
      null_estimable_fraction,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    `Mean type I error` = sprintf("%.2f%%", 100 * `Mean type I error`),
    `Maximum type I error` = sprintf("%.2f%%", 100 * `Maximum type I error`),
    `Mean calibrated size` = sprintf("%.2f%%", 100 * `Mean calibrated size`),
    `Median estimable fraction` = sprintf(
      "%.1f%%",
      100 * `Median estimable fraction`
    )
  ) %>%
  arrange(
    factor(Method, levels = unname(method_labels_text[method_order])),
    allocation_label,
    `Denominator class`
  )


fts2 <- flextable(table_s2)
fts2 <- style_ft(fts2, font_size = 8)
fts2 <- align(
  fts2,
  j = 3:ncol(table_s2),
  align = "center",
  part = "all"
)

fts2 <- set_caption(
  fts2,
  caption = paste0(
    "Table S2. Type I error and estimability across the final ",
    "Delta-parameterized null simulation grid."
  )
)

save_ft(
  fts2,
  file.path(SUPP_DIR, "Table_S2_type1_error_summary.docx")
)

write_csv(
  null_dat,
  file.path(DATA_DIR, "Table_S2_complete_null_grid.csv")
)


# -----------------------------------------------------------------------------
# TABLE S3
# Size-adjusted power summary
# -----------------------------------------------------------------------------

table_s3 <- power_dat %>%
  mutate(
    Method = unname(method_labels_text[method]),
    `Denominator class` = ifelse(
      denominator_class == "strong",
      "Strong",
      "Weak"
    )
  ) %>%
  group_by(
    Method,
    `Denominator class`,
    abs_population_Delta,
    expected_total_deaths_target
  ) %>%
  summarise(
    Scenarios = n(),
    `Mean estimable fraction` = mean(
      estimable_fraction,
      na.rm = TRUE
    ),
    `Mean nominal power` = mean(
      nominal_power,
      na.rm = TRUE
    ),
    `Mean size-adjusted power` = mean(
      size_adjusted_power,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    `|Delta|` = sprintf("%.2f", abs_population_Delta),
    `Expected critical events` = expected_total_deaths_target,
    `Mean estimable fraction` = sprintf(
      "%.1f%%",
      100 * `Mean estimable fraction`
    ),
    `Mean nominal power` = sprintf(
      "%.1f%%",
      100 * `Mean nominal power`
    ),
    `Mean size-adjusted power` = sprintf(
      "%.1f%%",
      100 * `Mean size-adjusted power`
    )
  ) %>%
  select(
    Method,
    `Denominator class`,
    `|Delta|`,
    `Expected critical events`,
    Scenarios,
    `Mean estimable fraction`,
    `Mean nominal power`,
    `Mean size-adjusted power`
  )


fts3 <- flextable(table_s3)
fts3 <- style_ft(fts3, font_size = 7.5)
fts3 <- align(
  fts3,
  j = 2:ncol(table_s3),
  align = "center",
  part = "all"
)

fts3 <- set_caption(
  fts3,
  caption = paste0(
    "Table S3. Statistical power across the final Delta-parameterized ",
    "alternative simulation grid."
  )
)

save_ft(
  fts3,
  file.path(SUPP_DIR, "Table_S3_size_adjusted_power.docx")
)

write_csv(
  power_dat,
  file.path(DATA_DIR, "Table_S3_complete_power_grid.csv")
)


# -----------------------------------------------------------------------------
# TABLE S4
# Complete empirical reanalysis
# -----------------------------------------------------------------------------

table_s4 <- emp %>%
  mutate(
    `First author` = first_author,
    `Critical events` = total_critical_events,
    `RR composite` = fmt_num(rr_composite, 2),
    `RR critical` = fmt_num(rr_critical, 2),
    BACO = fmt_num(baco, 2),
    `Fieller set` = ifelse(is.na(fieller_set), "", fieller_set),
    Delta = fmt_num(delta, 2),
    `RRc/RRd` = fmt_num(effect_ratio_c_over_d, 2),
    QRR = fmt_num(qrr, 3),
    `Fisher p` = fmt_p(fisher_p),
    `Original Wald p` = fmt_p(original_baco_wald_p),
    `Weak denominator` = ifelse(
      weak_denominator %in% TRUE,
      "Yes",
      ifelse(weak_denominator %in% FALSE, "No", "")
    )
  ) %>%
  select(
    `First author`,
    `Critical events`,
    `RR composite`,
    `RR critical`,
    BACO,
    `Fieller set`,
    Delta,
    `RRc/RRd`,
    QRR,
    `Fisher p`,
    `Original Wald p`,
    `Weak denominator`
  ) %>%
  arrange(desc(BACO))

fts4 <- flextable(table_s4)
fts4 <- style_ft(fts4, font_size = 6.5)

fts4 <- set_caption(
  fts4,
  caption = paste0(
    "Table S4. Complete empirical reanalysis of the 130 BACO comparisons."
  )
)

fts4 <- add_footer_lines(
  fts4,
  values = paste0(
    "Frith is retained in study accounting but has no valid aggregate ",
    "composite-event counts and is therefore unavailable for count-based ",
    "BACO, Fieller, Delta, and Fisher calculations."
  )
)

save_ft(
  fts4,
  file.path(SUPP_DIR, "Table_S4_complete_empirical_reanalysis.docx")
)

write_csv(
  emp,
  file.path(DATA_DIR, "Table_S4_complete_empirical_reanalysis.csv")
)


# -----------------------------------------------------------------------------
# TABLE S5
# Discordant conclusions
# -----------------------------------------------------------------------------

table_s5 <- emp %>%
  filter(
    method_agreement_05 %in%
      c(
        "fisher_only",
        "original_wald_only"
      ) |
      method_agreement_005 %in%
        c(
          "fisher_only",
          "original_wald_only"
        )
  ) %>%
  mutate(
    `First author` = first_author,
    `Critical events` = total_critical_events,
    `RR composite` = fmt_num(rr_composite, 3),
    `RR critical` = fmt_num(rr_critical, 3),
    BACO = fmt_num(baco, 3),
    `Fieller type` = ifelse(is.na(fieller_type), "", fieller_type),
    `Fisher p` = fmt_p(fisher_p),
    `Original Wald p` = fmt_p(original_baco_wald_p),
    `Decision at 0.05` = method_agreement_05,
    `Decision at 0.005` = method_agreement_005
  ) %>%
  select(
    `First author`,
    `Critical events`,
    `RR composite`,
    `RR critical`,
    BACO,
    `Fieller type`,
    `Fisher p`,
    `Original Wald p`,
    `Decision at 0.05`,
    `Decision at 0.005`
  ) %>%
  arrange(`Decision at 0.05`, `Fisher p`)


fts5 <- flextable(table_s5)
fts5 <- style_ft(fts5, font_size = 7)

fts5 <- set_caption(
  fts5,
  caption = paste0(
    "Table S5. Empirical comparisons with discordant conclusions between ",
    "the original BACO Wald procedure and Fisher exact testing."
  )
)

save_ft(
  fts5,
  file.path(SUPP_DIR, "Table_S5_method_disagreement.docx")
)

write_csv(
  table_s5,
  file.path(DATA_DIR, "Table_S5_method_disagreement.csv")
)


# -----------------------------------------------------------------------------
# TABLE S6
# Fisher-significant comparisons
# -----------------------------------------------------------------------------

table_s6 <- emp %>%
  filter(fisher_sig_05) %>%
  mutate(
    `First author` = first_author,
    `Critical events` = total_critical_events,
    `RR composite` = fmt_num(rr_composite, 3),
    `RR critical` = fmt_num(rr_critical, 3),
    BACO = fmt_num(baco, 3),
    `Fieller type` = ifelse(is.na(fieller_type), "", fieller_type),
    Delta = fmt_num(delta, 3),
    `RRc/RRd` = fmt_num(effect_ratio_c_over_d, 3),
    QRR = fmt_num(qrr, 3),
    `Fisher p` = fmt_p(fisher_p),
    `Original Wald p` = fmt_p(original_baco_wald_p)
  ) %>%
  select(
    `First author`,
    `Critical events`,
    `RR composite`,
    `RR critical`,
    BACO,
    `Fieller type`,
    Delta,
    `RRc/RRd`,
    QRR,
    `Fisher p`,
    `Original Wald p`
  ) %>%
  arrange(`Fisher p`)


fts6 <- flextable(table_s6)
fts6 <- style_ft(fts6, font_size = 7.5)

fts6 <- set_caption(
  fts6,
  caption = paste0(
    "Table S6. Comparisons with evidence of composite-critical effect ",
    "discrepancy according to Fisher exact testing at alpha = 0.05."
  )
)

save_ft(
  fts6,
  file.path(SUPP_DIR, "Table_S6_Fisher_significant_comparisons.docx")
)

write_csv(
  table_s6,
  file.path(DATA_DIR, "Table_S6_Fisher_significant_comparisons.csv")
)


# -----------------------------------------------------------------------------
# TABLE S7
# Complete external empirical replication
# -----------------------------------------------------------------------------

table_s7 <- rep %>%
  arrange(
    desc(primary_analysis),
    desc(baco)
  ) %>%
  mutate(
    `First author` = first_author,
    `Analysis set` = ifelse(
      primary_analysis,
      "Primary external replication",
      "Reoriented sensitivity only"
    ),
    `Data origin` = ifelse(
      data_origin == "reconstructed",
      "Reconstructed",
      "Direct"
    ),
    `Critical events` = total_critical_events,
    `RR composite` = fmt_num(rr_composite, 3),
    `RR critical` = fmt_num(rr_critical, 3),
    BACO = fmt_num(baco, 3),
    `Fieller type` = recode(
      fieller_type,
      bounded = "Bounded",
      disjoint_unbounded = "Disjoint unbounded",
      all_real = "All real",
      .default = fieller_type
    ),
    Delta = fmt_num(delta, 3),
    `RRc/RRd` = fmt_num(effect_ratio_c_over_d, 2),
    QRR = fmt_num(qrr, 3),
    `Fisher p` = fmt_p(fisher_p),
    `Fieller set` = ifelse(is.na(fieller_set), "", fieller_set),
    `Original Wald p` = fmt_p(original_baco_wald_p),
    `Weak denominator` = ifelse(
      weak_denominator %in% TRUE,
      "Yes",
      ifelse(weak_denominator %in% FALSE, "No", "")
    )
  ) %>%
  select(
    `First author`,
    `Critical events`,
    `RR composite`,
    `RR critical`,
    BACO,
    `Fieller set`,
    Delta,
    `RRc/RRd`,
    QRR,
    `Fisher p`,
    `Original Wald p`,
    `Weak denominator`
  )

fts7 <- flextable(table_s7)
fts7 <- style_ft(fts7, font_size = 6.5)

fts7 <- set_caption(
  fts7,
  caption = paste0(
    "Table S7. Complete external empirical replication of BACO comparisons ",
    "from COVID-19 randomized trials."
  )
)

fts7 <- add_footer_lines(
  fts7,
  values = paste0(
    "The primary external replication includes 22 comparisons with an adverse ",
    "composite endpoint and a nested critical component. Cremer and Temesgen ",
    "used favorable composite endpoints and are shown only after explicit ",
    "reorientation to the adverse complement for sensitivity analysis. Dequin ",
    "and Salama use transparently reconstructed aligned counts and are identified ",
    "as reconstructed. No continuity correction was applied."
  )
)

save_ft(
  fts7,
  file.path(SUPP_DIR, "Table_S7_external_replication_complete.docx")
)

write_csv(
  table_s7,
  file.path(DATA_DIR, "Table_S7_external_replication_complete.csv")
)


# -----------------------------------------------------------------------------
# TABLE S8
# External replication sensitivity analyses
# -----------------------------------------------------------------------------

summarize_replication_set <- function(dat, label) {
  n_set <- nrow(dat)
  n_fieller_set <- sum(dat$fieller_available, na.rm = TRUE)
  n_fisher_set <- sum(dat$fisher_available, na.rm = TRUE)
  n_delta_set <- sum(dat$delta_estimable, na.rm = TRUE)
  n_legacy_set <- sum(dat$original_baco_wald_available, na.rm = TRUE)

  tibble(
    `Analysis set` = label,
    Studies = as.character(n_distinct(dat$study_id)),
    Comparisons = as.character(n_set),
    `Weak denominator` = n_pct_text(
      sum(dat$weak_denominator, na.rm = TRUE),
      n_fieller_set
    ),
    `Bounded Fieller` = n_pct_text(
      sum(dat$fieller_type == "bounded", na.rm = TRUE),
      n_fieller_set
    ),
    `Disjoint unbounded` = n_pct_text(
      sum(dat$fieller_type == "disjoint_unbounded", na.rm = TRUE),
      n_fieller_set
    ),
    `All-real Fieller` = n_pct_text(
      sum(dat$fieller_type == "all_real", na.rm = TRUE),
      n_fieller_set
    ),
    `Fisher p < 0.05` = n_pct_text(
      sum(dat$fisher_sig_05, na.rm = TRUE),
      n_fisher_set
    ),
    `Delta Wald p < 0.05` = n_pct_text(
      sum(dat$delta_wald_sig_05, na.rm = TRUE),
      n_delta_set
    ),
    `Original BACO Wald p < 0.05` = n_pct_text(
      sum(dat$original_baco_wald_sig_05, na.rm = TRUE),
      n_legacy_set
    )
  )
}

table_s8 <- bind_rows(
  summarize_replication_set(
    rep_primary,
    "Primary nested comparisons"
  ),
  summarize_replication_set(
    rep_direct,
    "Direct-count-only sensitivity"
  ),
  summarize_replication_set(
    rep_all_reoriented,
    "All comparisons after adverse reorientation"
  )
)

fts8 <- flextable(table_s8)
fts8 <- style_ft(fts8, font_size = 7.2)
fts8 <- align(
  fts8,
  j = 2:ncol(table_s8),
  align = "center",
  part = "all"
)

fts8 <- set_caption(
  fts8,
  caption = paste0(
    "Table S8. Sensitivity analyses for the external empirical replication."
  )
)

fts8 <- add_footer_lines(
  fts8,
  values = paste0(
    "The three analysis sets contain 21, 19, and 23 unique studies, respectively. ",
    "Primary nested comparisons exclude the two trials whose reported composite ",
    "endpoint was favorable. Direct-count-only sensitivity additionally excludes ",
    "the reconstructed Dequin and Salama comparisons. The all-comparison analysis ",
    "includes Cremer and Temesgen only after explicit reorientation to the adverse ",
    "complement. Percentages use the method-specific available denominator."
  )
)

save_ft(
  fts8,
  file.path(SUPP_DIR, "Table_S8_external_replication_sensitivity.docx")
)

write_csv(
  table_s8,
  file.path(DATA_DIR, "Table_S8_external_replication_sensitivity.csv")
)


# =============================================================================
# SUPPLEMENTARY FIGURES (FINAL MANUSCRIPT)
# =============================================================================

# -----------------------------------------------------------------------------
# FIGURE S1
# Sparse-event type I error profile
# -----------------------------------------------------------------------------

figs1_data <- zero_dat %>%
  select(
    denominator_class,
    expected_total_deaths_target,
    fisher_type1_05,
    midp_type1_05,
    lr_type1_05,
    Delta_type1_05_unconditional,
    Wald_type1_05_unconditional
  ) %>%
  pivot_longer(
    cols = c(
      fisher_type1_05,
      midp_type1_05,
      lr_type1_05,
      Delta_type1_05_unconditional,
      Wald_type1_05_unconditional
    ),
    names_to = "method_raw",
    values_to = "type1"
  ) %>%
  mutate(
    Method = recode(
      method_raw,
      fisher_type1_05 = "Fisher exact",
      midp_type1_05 = "Mid-P",
      lr_type1_05 = "Conditional LR",
      Delta_type1_05_unconditional = "Delta Wald",
      Wald_type1_05_unconditional = "Original BACO Wald"
    ),
    Method = factor(
      Method,
      levels = c(
        "Fisher exact",
        "Mid-P",
        "Delta Wald",
        "Conditional LR",
        "Original BACO Wald"
      )
    ),
    `Critical-endpoint effect` = factor(
      denominator_class,
      levels = c("strong", "weak"),
      labels = c(
        "Strong critical-endpoint effect",
        "Weak critical-endpoint effect"
      )
    )
  ) %>%
  group_by(
    `Critical-endpoint effect`,
    Method,
    expected_total_deaths_target
  ) %>%
  summarise(
    mean_type1 = mean(type1, na.rm = TRUE),
    .groups = "drop"
  )


figs1 <- ggplot(
  figs1_data,
  aes(
    x = expected_total_deaths_target,
    y = mean_type1,
    shape = Method,
    color = Method,
    group = Method
  )
) +
  geom_hline(
    yintercept = 0.05,
    linetype = 3,
    linewidth = 0.5,
    color = "grey40"
  ) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.1) +
  facet_wrap(
    ~`Critical-endpoint effect`,
    nrow = 1
  ) +
  scale_x_continuous(
    breaks = sort(unique(
      figs1_data$expected_total_deaths_target
    ))
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.07))
  ) +
  scale_color_manual(
    values = c(
      "Fisher exact" = "#0072B2",
      "Mid-P" = "#D55E00",
      "Delta Wald" = "#009E73",
      "Conditional LR" = "#CC79A7",
      "Original BACO Wald" = "#E69F00"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Fisher exact" = 16,
      "Mid-P" = 17,
      "Delta Wald" = 15,
      "Conditional LR" = 18,
      "Original BACO Wald" = 8
    )
  ) +
  labs(
    x = "Expected total critical events",
    y = "Empirical type I error",
    color = "Method",
    shape = "Method"
  ) +
  publication_theme()

figs1

save_plot_both(
  figs1,
  "Figure_S1_sparse_event_type1",
  SUPP_DIR,
  width = 8.5,
  height = 5.0
)

write_csv(
  figs1_data,
  file.path(DATA_DIR, "Figure_S1_data.csv")
)


# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# FIGURE S2
# Theoretical relationship between BACO and Delta
# -----------------------------------------------------------------------------

delta_grid <- seq(-0.50, 0.50, by = 0.005)

rrd_values <- c(0.50, 0.95, 1.05, 2.00)

figs2_data <- tidyr::expand_grid(
  Delta = delta_grid,
  rr_d = rrd_values
) %>%
  mutate(
    BACO = 1 + Delta / log(rr_d),
    rr_label = factor(
      rr_d,
      levels = rrd_values,
      labels = c(
        "RRd = 0.50",
        "RRd = 0.95",
        "RRd = 1.05",
        "RRd = 2.00"
      )
    )
  )

figs2 <- ggplot(
  figs2_data,
  aes(
    x = Delta,
    y = BACO,
    color = rr_label,
    linetype = rr_label,
    group = rr_label
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = 3,
    linewidth = 0.45
  ) +
  geom_vline(
    xintercept = 0,
    linetype = 3,
    linewidth = 0.45
  ) +
  geom_line(
    linewidth = 0.75
  ) +
  coord_cartesian(
    ylim = c(-12, 12)
  ) +
  scale_color_manual(
    name = "Critical-endpoint effect",
    values = c(
      "RRd = 0.50" = "#D55E00",
      "RRd = 0.95" = "#984ea3",
      "RRd = 1.05" = "#0072B2",
      "RRd = 2.00" = "#009E73"
    )
  ) +
  scale_linetype_manual(
    name = "Critical-endpoint effect",
    values = c(
      "RRd = 0.50" = "solid",
      "RRd = 0.95" = "dashed",
      "RRd = 1.05" = "dotdash",
      "RRd = 2.00" = "longdash"
    )
  ) +
  labs(
    x = expression(
      Delta == log(RR[c]) - log(RR[d])
    ),
    y = "BACO"
  ) +
  publication_theme()

save_plot_both(
  figs2,
  "Figure_S2_BACO_vs_Delta_relationship",
  SUPP_DIR,
  width = 7.4,
  height = 5.4
)

write_csv(
  figs2_data,
  file.path(DATA_DIR, "Figure_S2_data.csv")
)


# -----------------------------------------------------------------------------

# =============================================================================
# FINAL SUPPLEMENTARY TITLES AND CAPTIONS
# =============================================================================

titles_and_captions_supp <- c(
  "SUPPLEMENT",
  "==========",
  "",
  "Table S1. Final simulation design and inferential procedures.",
  "Table S2. Type I error and estimability across the final Delta-parameterized null simulation grid.",
  "Table S3. Statistical power across the final Delta-parameterized alternative simulation grid.",
  "Table S4. Complete empirical reanalysis of the 130 BACO comparisons.",
  "Table S5. Empirical comparisons with discordant conclusions between the original BACO Wald procedure and Fisher exact testing.",
  "Table S6. Comparisons with evidence of composite-critical effect discrepancy according to Fisher exact testing at alpha = 0.05.",
  "Table S7. Complete external empirical replication of BACO comparisons from COVID-19 randomized trials.",
  "Table S8. Sensitivity analyses for the external empirical replication.",
  "",
  "Figure S1. Empirical type I error under sparse critical events.",
  "Note: Fisher's exact test remains conservative when critical events are very sparse, whereas differences among the regular procedures diminish as the amount of information increases. Results show unconditional rejection probabilities at alpha = 0.05 across expected numbers of critical events.",
  "",
  "Figure S2. Relationship between the BACO ratio and the direct log-effect discrepancy.",
  "Note: For a given Delta, the BACO ratio becomes increasingly extreme as the critical-endpoint relative effect (RRd) approaches 1."
)

writeLines(
  titles_and_captions_supp,
  file.path(SUPP_DIR, "titles_and_captions_supp.txt")
)

cat("\n============================================================\n")
cat("BACO manuscript tables and figures created successfully.\n")
cat("============================================================\n\n")
cat(
  "Output directory:\n  ",
  normalizePath(OUTPUT_DIR, winslash = "/", mustWork = FALSE),
  "\n\n",
  sep = ""
)
cat("MAIN MANUSCRIPT\n")
cat("  Table 1: Table_1_operating_characteristics.docx\n")
cat("  Table 2: Table_2_empirical_reanalysis_external_replication.docx\n")
cat("  Figure 1: Figure_1_type1_error.pdf / .png\n")
cat("  Figure 2: Figure_2_size_adjusted_power_Delta.pdf / .png\n")
cat("  Figure 3: Figure_3_Wald_vs_Fisher.pdf / .png\n\n")
cat("SUPPLEMENT\n")
cat("  Tables S1-S8: individual .docx files\n")
cat("  Figures S1-S2: PDF + 600-dpi PNG\n")
