# =============================================================================
# BACO reproducibility analysis
# Empirical reanalysis and external COVID-19 replication
# =============================================================================
#
# Input
# -----
# data/Supplementary_Data_Studies_standardized.xlsx
#   sheet: Data
#
# This is the same standardized comparison-level dataset intended for public
# availability with the manuscript. The analysis deliberately reconstructs all
# count-based quantities from the six arm-level counts rather than relying on
# the derived columns already stored in the workbook.
#
# Primary methods
# ---------------
#   - BACO point estimate
#   - Fieller confidence set and weak-denominator diagnostic
#   - Delta = log(RR_c) - log(RR_d)
#   - QRR = RR_d / RR_c = exp(-Delta)
#   - probability-ordered two-sided Fisher exact test
#   - Delta Wald test (secondary comparator)
#   - original BACO Wald/Taylor procedure (historical comparator)
#
# No continuity correction is used.
#
# Requires bacoR (the manuscript analyses used version 0.0.9).
# =============================================================================

ROOT_DIR <- getwd()
INPUT_FILE <- file.path(ROOT_DIR, "data", "Supplementary_Data_Studies_standardized.xlsx")
INPUT_SHEET <- "Data"
RESULTS_DIR <- file.path(ROOT_DIR, "results", "empirical")

ALPHA <- 0.05
ALPHA_ORIGINAL <- 0.005

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

required_packages <- c("readxl", "readr", "dplyr", "purrr", "tibble", "bacoR")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(bacoR)
})

if (!file.exists(INPUT_FILE)) {
  stop("Input file not found: ", INPUT_FILE)
}

cat("bacoR version: ", as.character(utils::packageVersion("bacoR")), "\n", sep = "")

# -----------------------------------------------------------------------------
# 1. Read and validate the standardized public dataset
# -----------------------------------------------------------------------------

raw <- readxl::read_excel(INPUT_FILE, sheet = INPUT_SHEET)

required_columns <- c(
  "dataset", "study_id", "comparison_id", "first_author", "doi",
  "source_group", "analysis_set", "data_origin", "primary_analysis",
  "direct_count_sensitivity", "reoriented_for_sensitivity",
  "aggregate_count_valid", "n_treatment", "n_control",
  "composite_treatment", "composite_control",
  "critical_treatment", "critical_control", "total_critical_events"
)
stopifnot(all(required_columns %in% names(raw)))

primary_input <- raw %>%
  filter(dataset == "primary_empirical_reanalysis")

external_input <- raw %>%
  filter(dataset == "external_covid19_replication")

stopifnot(
  nrow(primary_input) == 130,
  n_distinct(primary_input$study_id) == 123,
  sum(primary_input$aggregate_count_valid, na.rm = TRUE) == 129,
  sum(primary_input$source_group == "published_benchmark", na.rm = TRUE) == 97,
  n_distinct(primary_input$study_id[primary_input$source_group == "published_benchmark"]) == 91,
  sum(primary_input$source_group == "additional_restored", na.rm = TRUE) == 33,
  n_distinct(primary_input$study_id[primary_input$source_group == "additional_restored"]) == 32,

  nrow(external_input) == 24,
  n_distinct(external_input$study_id) == 23,
  sum(external_input$primary_analysis, na.rm = TRUE) == 22,
  sum(external_input$direct_count_sensitivity, na.rm = TRUE) == 20,
  sum(external_input$reoriented_for_sensitivity, na.rm = TRUE) == 2,
  sum(external_input$data_origin == "reconstructed", na.rm = TRUE) == 2
)

# The two empirical datasets are independent in the manuscript.
primary_doi <- tolower(trimws(primary_input$doi[!is.na(primary_input$doi)]))
external_doi <- tolower(trimws(external_input$doi[!is.na(external_input$doi)]))
stopifnot(length(intersect(primary_doi, external_doi)) == 0)

# -----------------------------------------------------------------------------
# 2. Original BACO Wald/Taylor comparator
# -----------------------------------------------------------------------------

original_baco_wald <- function(
    n_treatment,
    n_control,
    composite_treatment,
    composite_control,
    critical_treatment,
    critical_control,
    conf_level = 0.95) {

  out_na <- list(
    estimate = NA_real_,
    se = NA_real_,
    ci_lower = NA_real_,
    ci_upper = NA_real_,
    z = NA_real_,
    p = NA_real_
  )

  vals <- c(
    n_treatment, n_control,
    composite_treatment, composite_control,
    critical_treatment, critical_control
  )

  if (any(!is.finite(vals))) return(out_na)
  if (n_treatment <= 0 || n_control <= 0) return(out_na)
  if (composite_treatment <= 0 || composite_control <= 0 ||
      critical_treatment <= 0 || critical_control <= 0) return(out_na)

  rr_c <- (composite_treatment / n_treatment) /
    (composite_control / n_control)
  rr_d <- (critical_treatment / n_treatment) /
    (critical_control / n_control)

  if (!is.finite(rr_c) || !is.finite(rr_d) || rr_c <= 0 || rr_d <= 0) {
    return(out_na)
  }

  beta_c <- log(rr_c)
  beta_d <- log(rr_d)
  tol <- sqrt(.Machine$double.eps)

  if (abs(beta_c) <= tol || abs(beta_d) <= tol) return(out_na)

  var_c <- 1 / composite_treatment + 1 / composite_control -
    1 / n_treatment - 1 / n_control
  var_d <- 1 / critical_treatment + 1 / critical_control -
    1 / n_treatment - 1 / n_control
  cov_cd <- var_c

  theta <- beta_c / beta_d
  var_theta <- theta^2 * (
    var_c / beta_c^2 +
      var_d / beta_d^2 -
      2 * cov_cd / (beta_c * beta_d)
  )

  if (!is.finite(var_theta) || var_theta <= 0) return(out_na)

  se_theta <- sqrt(var_theta)
  zcrit <- qnorm(1 - (1 - conf_level) / 2)
  z_value <- (theta - 1) / se_theta
  p_value <- 2 * pnorm(-abs(z_value))

  list(
    estimate = theta,
    se = se_theta,
    ci_lower = theta - zcrit * se_theta,
    ci_upper = theta + zcrit * se_theta,
    z = z_value,
    p = p_value
  )
}

# -----------------------------------------------------------------------------
# 3. Fieller helpers
# -----------------------------------------------------------------------------

format_endpoint <- function(x, digits = 6) {
  if (is.infinite(x) && x < 0) return("-Inf")
  if (is.infinite(x) && x > 0) return("Inf")
  if (!is.finite(x)) return("NA")
  formatC(x, digits = digits, format = "fg")
}

collapse_fieller <- function(fieller) {
  if (is.null(fieller)) return(NA_character_)
  ints <- fieller$intervals
  if (is.null(ints) || nrow(ints) == 0) return(fieller$set_type)

  pieces <- vapply(seq_len(nrow(ints)), function(i) {
    lo <- ints$lower[i]
    hi <- ints$upper[i]
    left <- if (is.infinite(lo)) "(" else "["
    right <- if (is.infinite(hi)) ")" else "]"
    paste0(left, format_endpoint(lo), ", ", format_endpoint(hi), right)
  }, character(1))

  paste(pieces, collapse = " U ")
}

bounded_limits <- function(fieller) {
  if (is.null(fieller) || !identical(fieller$set_type, "bounded") ||
      nrow(fieller$intervals) != 1) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  c(lower = fieller$intervals$lower[1], upper = fieller$intervals$upper[1])
}

# -----------------------------------------------------------------------------
# 4. Recalculate one comparison from raw counts
# -----------------------------------------------------------------------------

analyse_one <- function(row) {
  base <- tibble(
    dataset = row$dataset,
    study_id = row$study_id,
    comparison_id = row$comparison_id,
    first_author = row$first_author,
    doi = row$doi,
    source_group = row$source_group,
    analysis_set = row$analysis_set,
    data_origin = row$data_origin,
    primary_analysis = row$primary_analysis,
    direct_count_sensitivity = row$direct_count_sensitivity,
    reoriented_for_sensitivity = row$reoriented_for_sensitivity,
    aggregate_count_valid = row$aggregate_count_valid,
    n_treatment = row$n_treatment,
    n_control = row$n_control,
    composite_treatment = row$composite_treatment,
    composite_control = row$composite_control,
    critical_treatment = row$critical_treatment,
    critical_control = row$critical_control,
    total_critical_events = row$total_critical_events
  )

  if (!isTRUE(row$aggregate_count_valid)) {
    return(base %>% mutate(
      rr_composite = NA_real_, rr_critical = NA_real_,
      beta_composite = NA_real_, beta_critical = NA_real_, baco = NA_real_,
      fieller_type = NA_character_, fieller_set = NA_character_,
      fieller_lower = NA_real_, fieller_upper = NA_real_,
      weak_denominator = NA,
      delta = NA_real_, delta_ci_lower = NA_real_, delta_ci_upper = NA_real_,
      delta_wald_p = NA_real_, effect_ratio_c_over_d = NA_real_, qrr = NA_real_,
      fisher_p = NA_real_, fisher_note = "Aggregate-count analysis unavailable.",
      zero_critical_any_arm = NA, zero_critical_both_arms = NA,
      original_baco_wald = NA_real_, original_baco_wald_ci_lower = NA_real_,
      original_baco_wald_ci_upper = NA_real_, original_baco_wald_p = NA_real_
    ))
  }

  res <- bacoR::baco_counts(
    n_treatment = row$n_treatment,
    n_control = row$n_control,
    composite_treatment = row$composite_treatment,
    composite_control = row$composite_control,
    critical_treatment = row$critical_treatment,
    critical_control = row$critical_control,
    conf_level = 0.95,
    midp = FALSE,
    interpret = FALSE
  )

  old <- original_baco_wald(
    row$n_treatment, row$n_control,
    row$composite_treatment, row$composite_control,
    row$critical_treatment, row$critical_control
  )

  flims <- bounded_limits(res$baco$fieller)

  base %>% mutate(
    rr_composite = res$effects$rr_composite,
    rr_critical = res$effects$rr_critical,
    beta_composite = res$effects$beta_composite,
    beta_critical = res$effects$beta_critical,
    baco = res$baco$estimate,
    fieller_type = if (is.null(res$baco$fieller)) NA_character_ else res$baco$fieller$set_type,
    fieller_set = collapse_fieller(res$baco$fieller),
    fieller_lower = unname(flims["lower"]),
    fieller_upper = unname(flims["upper"]),
    weak_denominator = res$diagnostics$weak_denominator,
    delta = res$discrepancy$delta,
    delta_ci_lower = unname(res$discrepancy$delta_ci["lower"]),
    delta_ci_upper = unname(res$discrepancy$delta_ci["upper"]),
    delta_wald_p = res$discrepancy$delta_p_wald,
    effect_ratio_c_over_d = ifelse(is.finite(res$discrepancy$delta),
                                   exp(res$discrepancy$delta), NA_real_),
    qrr = res$discrepancy$qrr,
    fisher_p = res$exact_test$fisher_p,
    fisher_note = if (is.null(res$exact_test$note)) NA_character_ else res$exact_test$note,
    zero_critical_any_arm = res$diagnostics$zero_critical_any_arm,
    zero_critical_both_arms = res$diagnostics$zero_critical_both_arms,
    original_baco_wald = old$estimate,
    original_baco_wald_ci_lower = old$ci_lower,
    original_baco_wald_ci_upper = old$ci_upper,
    original_baco_wald_p = old$p
  )
}

results <- purrr::map_dfr(seq_len(nrow(raw)), function(i) analyse_one(raw[i, ])) %>%
  mutate(
    baco_estimable = is.finite(baco),
    delta_estimable = is.finite(delta),
    fieller_available = !is.na(fieller_type),
    fisher_available = is.finite(fisher_p),
    original_baco_wald_available = is.finite(original_baco_wald_p),

    fisher_sig_05 = fisher_available & fisher_p < ALPHA,
    fisher_sig_005 = fisher_available & fisher_p < ALPHA_ORIGINAL,
    delta_wald_sig_05 = is.finite(delta_wald_p) & delta_wald_p < ALPHA,
    delta_wald_sig_005 = is.finite(delta_wald_p) & delta_wald_p < ALPHA_ORIGINAL,
    original_baco_wald_sig_05 = original_baco_wald_available & original_baco_wald_p < ALPHA,
    original_baco_wald_sig_005 = original_baco_wald_available & original_baco_wald_p < ALPHA_ORIGINAL,

    method_agreement_05 = case_when(
      !fisher_available | !original_baco_wald_available ~ "not_comparable",
      fisher_sig_05 & original_baco_wald_sig_05 ~ "both_significant",
      !fisher_sig_05 & !original_baco_wald_sig_05 ~ "both_nonsignificant",
      fisher_sig_05 & !original_baco_wald_sig_05 ~ "fisher_only",
      TRUE ~ "original_wald_only"
    ),
    method_agreement_005 = case_when(
      !fisher_available | !original_baco_wald_available ~ "not_comparable",
      fisher_sig_005 & original_baco_wald_sig_005 ~ "both_significant",
      !fisher_sig_005 & !original_baco_wald_sig_005 ~ "both_nonsignificant",
      fisher_sig_005 & !original_baco_wald_sig_005 ~ "fisher_only",
      TRUE ~ "original_wald_only"
    )
  )

primary <- results %>% filter(dataset == "primary_empirical_reanalysis")
external <- results %>% filter(dataset == "external_covid19_replication")
external_primary <- external %>% filter(primary_analysis %in% TRUE)
external_direct <- external %>% filter(direct_count_sensitivity %in% TRUE)

# -----------------------------------------------------------------------------
# 5. Freeze the manuscript empirical results
# -----------------------------------------------------------------------------

stopifnot(
  nrow(primary) == 130,
  sum(primary$aggregate_count_valid, na.rm = TRUE) == 129,
  sum(primary$fisher_available, na.rm = TRUE) == 129,
  sum(primary$baco_estimable, na.rm = TRUE) == 111,
  sum(primary$delta_estimable, na.rm = TRUE) == 114,
  sum(primary$fieller_available, na.rm = TRUE) == 114,
  sum(primary$zero_critical_any_arm, na.rm = TRUE) == 15,
  sum(primary$zero_critical_both_arms, na.rm = TRUE) == 6,
  sum(primary$weak_denominator, na.rm = TRUE) == 106,
  sum(primary$fieller_type == "bounded", na.rm = TRUE) == 8,
  sum(primary$fieller_type == "disjoint_unbounded", na.rm = TRUE) == 27,
  sum(primary$fieller_type == "all_real", na.rm = TRUE) == 79,
  sum(primary$fisher_sig_05, na.rm = TRUE) == 5,
  sum(primary$delta_wald_sig_05, na.rm = TRUE) == 8,
  sum(primary$original_baco_wald_sig_05, na.rm = TRUE) == 27,
  sum(primary$fisher_sig_005, na.rm = TRUE) == 1,
  sum(primary$original_baco_wald_sig_005, na.rm = TRUE) == 11,

  nrow(external_primary) == 22,
  n_distinct(external_primary$study_id) == 21,
  sum(external_primary$zero_critical_any_arm, na.rm = TRUE) == 0,
  sum(external_primary$weak_denominator, na.rm = TRUE) == 21,
  sum(external_primary$fieller_type == "bounded", na.rm = TRUE) == 1,
  sum(external_primary$fieller_type == "disjoint_unbounded", na.rm = TRUE) == 5,
  sum(external_primary$fieller_type == "all_real", na.rm = TRUE) == 16,
  sum(external_primary$fisher_sig_05, na.rm = TRUE) == 2,
  sum(external_primary$delta_wald_sig_05, na.rm = TRUE) == 2,
  sum(external_primary$original_baco_wald_sig_05, na.rm = TRUE) == 4,

  nrow(external_direct) == 20,
  n_distinct(external_direct$study_id) == 19,
  sum(external_direct$weak_denominator, na.rm = TRUE) == 19,
  sum(external_direct$fisher_sig_05, na.rm = TRUE) == 1,
  sum(external_direct$delta_wald_sig_05, na.rm = TRUE) == 1,
  sum(external_direct$original_baco_wald_sig_05, na.rm = TRUE) == 3,

  nrow(external) == 24,
  sum(external$weak_denominator, na.rm = TRUE) == 23,
  sum(external$fieller_type == "bounded", na.rm = TRUE) == 1,
  sum(external$fieller_type == "disjoint_unbounded", na.rm = TRUE) == 5,
  sum(external$fieller_type == "all_real", na.rm = TRUE) == 18,
  sum(external$fisher_sig_05, na.rm = TRUE) == 2,
  sum(external$delta_wald_sig_05, na.rm = TRUE) == 2,
  sum(external$original_baco_wald_sig_05, na.rm = TRUE) == 4
)

# Validate selected recalculated quantities against the public workbook itself.
# This is a reproducibility check only; the calculations above use raw counts.
check <- results %>%
  select(dataset, study_id, comparison_id, baco, delta, fisher_p,
         original_baco_wald_p, fieller_type) %>%
  left_join(
    raw %>% select(dataset, study_id, comparison_id,
                   baco_stored = baco,
                   delta_stored = delta,
                   fisher_p_stored = fisher_p,
                   original_baco_wald_p_stored = original_baco_wald_p,
                   fieller_type_stored = fieller_type),
    by = c("dataset", "study_id", "comparison_id")
  )

same_num <- function(x, y, tol = 1e-8) {
  both_na <- is.na(x) & is.na(y)
  both_inf <- is.infinite(x) & is.infinite(y) & sign(x) == sign(y)
  close <- is.finite(x) & is.finite(y) & abs(x - y) <= tol
  both_na | both_inf | close
}

stopifnot(
  all(same_num(check$baco, check$baco_stored)),
  all(same_num(check$delta, check$delta_stored)),
  all(same_num(check$fisher_p, check$fisher_p_stored)),
  all(same_num(check$original_baco_wald_p,
               check$original_baco_wald_p_stored)),
  all((is.na(check$fieller_type) & is.na(check$fieller_type_stored)) |
        check$fieller_type == check$fieller_type_stored)
)

# -----------------------------------------------------------------------------
# 6. Write analysis outputs used by the manuscript table/figure script
# -----------------------------------------------------------------------------

readr::write_csv(primary, file.path(RESULTS_DIR, "empirical_baco_full.csv"), na = "")
readr::write_csv(external, file.path(RESULTS_DIR, "external_replication_baco_full.csv"), na = "")
readr::write_csv(results, file.path(RESULTS_DIR, "combined_empirical_results.csv"), na = "")

summary_lines <- c(
  "BACO empirical reanalysis — reproducibility summary",
  "=================================================",
  paste0("bacoR version: ", as.character(utils::packageVersion("bacoR"))),
  "",
  paste0("Primary dataset: ", nrow(primary), " comparisons; ",
         n_distinct(primary$study_id), " studies"),
  paste0("Valid aggregate counts: ", sum(primary$aggregate_count_valid, na.rm = TRUE)),
  paste0("Weak denominator: ", sum(primary$weak_denominator, na.rm = TRUE), "/",
         sum(primary$fieller_available, na.rm = TRUE)),
  paste0("Fisher / Delta Wald / original BACO Wald rejections at 0.05: ",
         sum(primary$fisher_sig_05, na.rm = TRUE), " / ",
         sum(primary$delta_wald_sig_05, na.rm = TRUE), " / ",
         sum(primary$original_baco_wald_sig_05, na.rm = TRUE)),
  "",
  paste0("External dataset: ", nrow(external), " comparisons; ",
         n_distinct(external$study_id), " studies"),
  paste0("Primary external set: ", nrow(external_primary), " comparisons; ",
         n_distinct(external_primary$study_id), " studies"),
  paste0("Fisher / Delta Wald / original BACO Wald rejections at 0.05: ",
         sum(external_primary$fisher_sig_05, na.rm = TRUE), " / ",
         sum(external_primary$delta_wald_sig_05, na.rm = TRUE), " / ",
         sum(external_primary$original_baco_wald_sig_05, na.rm = TRUE))
)
writeLines(summary_lines, file.path(RESULTS_DIR, "empirical_reanalysis_summary.txt"))

cat("\nEmpirical analyses completed successfully.\n")
cat("Outputs: ", normalizePath(RESULTS_DIR, winslash = "/", mustWork = FALSE), "\n", sep = "")
