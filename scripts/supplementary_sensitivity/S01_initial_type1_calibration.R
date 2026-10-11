# ==============================================================================
# BACO reproducibility analysis — supplementary sensitivity
# Initial Monte Carlo calibration of inference under BACO = 1
#
# Purpose
# -------
# Evaluate whether the original Wald/delta inference for
#
#       BACO = beta_c / beta_d
#
# is correctly calibrated when:
#   1) the critical outcome is sparse; and/or
#   2) the reference effect beta_d = log(RR_d) approaches zero.
#
# The main null hypothesis is:
#
#       H0: BACO = 1
#
# which, whenever beta_d != 0, is equivalent to:
#
#       H0: beta_c - beta_d = 0.
#
# This script compares:
#   A) Original BACO Wald/delta test;
#   B) Direct Wald test of Delta = beta_c - beta_d;
#   C) Fieller confidence-set decision for theta = 1.
#
# IMPORTANT DESIGN FEATURE
# ------------------------
# Under H0 and away from RR = 1:
#
#       RR_c = RR_d = R
#
# so the population BACO is exactly 1.
#
# When R = 1 exactly:
#
#       beta_c = beta_d = 0
#
# and BACO = 0/0 is NOT a defined population parameter.
#
# Therefore, RR = 1 is treated separately as a structural
# "reference-null" scenario and is NOT used to calculate type-I error.
#
# Data-generating structure
# -------------------------
# Each participant belongs to one of three mutually exclusive states:
#
#       0 = no component event
#       1 = nonfatal composite component
#       2 = critical event (death)
#
# with death embedded within the composite outcome.
#
# For the control group:
#       P(death)     = p_d0
#       P(composite) = p_c0
#
# For the experimental group under BACO = 1:
#       P(death)     = R * p_d0
#       P(composite) = R * p_c0
#
# This guarantees:
#       RR_d = RR_c = R
#
# provided all probabilities are valid.
#
# Recommended workflow
# --------------------
# Pilot run:
#       NSIM <- 2000
#
# More stable manuscript run:
#       NSIM <- 10000
#
# Final simulation after methods are frozen:
#       NSIM <- 50000 (or enough to achieve desired Monte Carlo error)
#
# At p = 0.05 and NSIM = 10000, the Monte Carlo SE is approximately:
#       sqrt(.05 * .95 / 10000) = 0.00218
#
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------------------------

required_packages <- c(
  "dplyr",
  "purrr",
  "tibble",
  "tidyr",
  "readr",
  "ggplot2"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", "),
    "\nExample:\ninstall.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "),
    "))"
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
  library(tidyr)
  library(readr)
  library(ggplot2)
})


# ------------------------------------------------------------------------------
# 1. Simulation settings
# ------------------------------------------------------------------------------

ROOT_DIR <- getwd()
OUTPUT_DIR <- file.path(ROOT_DIR, "results", "supplementary_initial_calibration")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

SEED <- 20260902

# Start with 2,000. Increase after checking that the scenario design and
# implementation behave as expected.
NSIM <- 2000

ALPHA_MAIN <- 0.05
ALPHA_STRICT <- 0.005

# Equal allocation is used in this first calibration study.
# Unequal allocation can be added later as a separate factor.
N_PER_ARM <- c(
  100,
  250,
  500,
  1000,
  5000
)

# Baseline probability of the critical endpoint in the control arm.
P_DEATH_CONTROL <- c(
  0.0025,
  0.005,
  0.01,
  0.02,
  0.05
)

# Baseline probability of the composite endpoint in the control arm.
# Always greater than P(death).
P_COMPOSITE_CONTROL <- c(
  0.10,
  0.25
)

# Reference effect. Values close to 1 deliberately stress the ratio denominator.
# R = 1 is included only as a structural null-reference diagnostic scenario.
RR_REFERENCE <- c(
  0.50,
  0.80,
  0.95,
  0.99,
  1.00,
  1.01,
  1.05,
  1.25,
  2.00
)

# Numerical tolerance ONLY for computational equality to zero.
# This is not a scientific threshold for "near-null".
NUMERIC_TOL <- sqrt(.Machine$double.eps)


# ------------------------------------------------------------------------------
# 2. Scenario grid
# ------------------------------------------------------------------------------

scenario_grid <- tidyr::crossing(
  n_per_arm = N_PER_ARM,
  p_death_control = P_DEATH_CONTROL,
  p_composite_control = P_COMPOSITE_CONTROL,
  rr_reference = RR_REFERENCE
) %>%
  mutate(
    p_death_experimental =
      rr_reference * p_death_control,

    p_composite_experimental =
      rr_reference * p_composite_control,

    valid_probability_structure =
      p_death_control >= 0 &
      p_death_control <= p_composite_control &
      p_composite_control <= 1 &
      p_death_experimental >= 0 &
      p_death_experimental <= p_composite_experimental &
      p_composite_experimental <= 1,

    population_beta_reference = log(rr_reference),

    population_baco = ifelse(
      rr_reference == 1,
      NA_real_,
      1
    ),

    scenario_class = ifelse(
      rr_reference == 1,
      "reference_null_BACO_undefined",
      "BACO_equals_1"
    ),

    expected_deaths_control =
      n_per_arm * p_death_control,

    expected_deaths_experimental =
      n_per_arm * p_death_experimental,

    expected_total_deaths =
      expected_deaths_control +
      expected_deaths_experimental,

    expected_composites_control =
      n_per_arm * p_composite_control,

    expected_composites_experimental =
      n_per_arm * p_composite_experimental,

    abs_population_beta_reference =
      abs(population_beta_reference)
  ) %>%
  filter(valid_probability_structure) %>%
  mutate(
    scenario_id = row_number()
  ) %>%
  select(
    scenario_id,
    everything()
  )

write_csv(
  scenario_grid,
  file.path(OUTPUT_DIR, "00_simulation_scenario_grid.csv")
)

cat(
  "\nNumber of valid simulation scenarios: ",
  nrow(scenario_grid),
  "\n",
  sep = ""
)


# ------------------------------------------------------------------------------
# 3. Vectorized Fieller classification
# ------------------------------------------------------------------------------

fieller_vectorized <- function(
    beta_c,
    beta_d,
    var_c,
    var_d,
    cov_cd,
    alpha = 0.05,
    tol = 1e-12
) {

  z <- qnorm(1 - alpha / 2)

  A <- beta_d^2 - z^2 * var_d
  B <- -2 * beta_c * beta_d + 2 * z^2 * cov_cd
  C <- beta_c^2 - z^2 * var_c

  discr <- B^2 - 4 * A * C

  fieller_type <- rep(NA_character_, length(beta_c))

  finite_input <-
    is.finite(beta_c) &
    is.finite(beta_d) &
    is.finite(var_c) &
    is.finite(var_d) &
    is.finite(cov_cd)

  # Conventional quadratic cases.
  fieller_type[
    finite_input &
    abs(A) > tol &
    discr >= 0 &
    A > 0
  ] <- "bounded"

  fieller_type[
    finite_input &
    abs(A) > tol &
    discr >= 0 &
    A < 0
  ] <- "disjoint_unbounded"

  fieller_type[
    finite_input &
    abs(A) > tol &
    discr < 0 &
    A < 0
  ] <- "all_real"

  fieller_type[
    finite_input &
    abs(A) > tol &
    discr < 0 &
    A > 0
  ] <- "empty"

  # Near-linear cases.
  linear_case <-
    finite_input &
    abs(A) <= tol

  fieller_type[
    linear_case &
    abs(B) > tol
  ] <- "half_line"

  fieller_type[
    linear_case &
    abs(B) <= tol &
    C <= 0
  ] <- "all_real"

  fieller_type[
    linear_case &
    abs(B) <= tol &
    C > 0
  ] <- "empty"

  fieller_type[!finite_input] <- "not_available"

  # Does the Fieller confidence set contain theta = 1?
  # Evaluate the quadratic inequality directly at theta = 1.
  q_at_1 <- A + B + C

  contains_1 <- ifelse(
    finite_input,
    q_at_1 <= tol,
    NA
  )

  tibble(
    fieller_type = fieller_type,
    fieller_contains_1 = contains_1
  )
}


# ------------------------------------------------------------------------------
# 4. Simulate one scenario
# ------------------------------------------------------------------------------

simulate_one_scenario <- function(
    scenario_row,
    nsim = 2000,
    alpha_main = 0.05,
    alpha_strict = 0.005,
    base_seed = 20260902
) {

  s <- scenario_row

  n_t <- s$n_per_arm
  n_c <- s$n_per_arm

  p_d0 <- s$p_death_control
  p_c0 <- s$p_composite_control

  p_d1 <- s$p_death_experimental
  p_c1 <- s$p_composite_experimental

  # Mutually exclusive outcome probabilities:
  # none / nonfatal composite / death
  prob_control <- c(
    none = 1 - p_c0,
    nonfatal_composite = p_c0 - p_d0,
    death = p_d0
  )

  prob_experimental <- c(
    none = 1 - p_c1,
    nonfatal_composite = p_c1 - p_d1,
    death = p_d1
  )

  if (
    any(prob_control < 0) ||
    any(prob_experimental < 0) ||
    abs(sum(prob_control) - 1) > 1e-10 ||
    abs(sum(prob_experimental) - 1) > 1e-10
  ) {
    stop("Invalid multinomial probabilities in scenario ", s$scenario_id)
  }

  set.seed(base_seed + s$scenario_id * 1009L)

  control_draws <- rmultinom(
    n = nsim,
    size = n_c,
    prob = prob_control
  )

  experimental_draws <- rmultinom(
    n = nsim,
    size = n_t,
    prob = prob_experimental
  )

  comp_c <-
    control_draws["nonfatal_composite", ] +
    control_draws["death", ]

  comp_t <-
    experimental_draws["nonfatal_composite", ] +
    experimental_draws["death", ]

  death_c <- control_draws["death", ]
  death_t <- experimental_draws["death", ]

  total_deaths <- death_c + death_t

  zero_death_any_arm <-
    death_c == 0 |
    death_t == 0

  zero_death_both_arms <-
    death_c == 0 &
    death_t == 0

  zero_composite_any_arm <-
    comp_c == 0 |
    comp_t == 0

  # Original analytic BACO requires nonzero critical and composite cells.
  estimable_original <-
    !zero_death_any_arm &
    !zero_composite_any_arm

  rr_c <- rep(NA_real_, nsim)
  rr_d <- rep(NA_real_, nsim)

  beta_c <- rep(NA_real_, nsim)
  beta_d <- rep(NA_real_, nsim)

  var_c <- rep(NA_real_, nsim)
  var_d <- rep(NA_real_, nsim)
  cov_cd <- rep(NA_real_, nsim)

  idx <- which(estimable_original)

  if (length(idx) > 0) {

    rr_c[idx] <-
      (comp_t[idx] / n_t) /
      (comp_c[idx] / n_c)

    rr_d[idx] <-
      (death_t[idx] / n_t) /
      (death_c[idx] / n_c)

    beta_c[idx] <- log(rr_c[idx])
    beta_d[idx] <- log(rr_d[idx])

    var_c[idx] <-
      1 / comp_t[idx] +
      1 / comp_c[idx] -
      1 / n_t -
      1 / n_c

    var_d[idx] <-
      1 / death_t[idx] +
      1 / death_c[idx] -
      1 / n_t -
      1 / n_c

    # Original BACO covariance result:
    cov_cd[idx] <- var_c[idx]
  }

  # --------------------------------------------------------------------------
  # A) Original BACO Wald/delta inference
  # --------------------------------------------------------------------------

  baco_hat <- rep(NA_real_, nsim)
  baco_se <- rep(NA_real_, nsim)

  baco_ci_low <- rep(NA_real_, nsim)
  baco_ci_high <- rep(NA_real_, nsim)

  p_baco_wald <- rep(NA_real_, nsim)

  estimable_ratio <-
    estimable_original &
    is.finite(beta_d) &
    abs(beta_d) > NUMERIC_TOL

  ridx <- which(estimable_ratio)

  if (length(ridx) > 0) {

    baco_hat[ridx] <-
      beta_c[ridx] /
      beta_d[ridx]

    # Delta-method variance using the gradient of beta_c / beta_d.
    g1 <- 1 / beta_d[ridx]
    g2 <- -beta_c[ridx] / (beta_d[ridx]^2)

    var_baco <-
      g1^2 * var_c[ridx] +
      g2^2 * var_d[ridx] +
      2 * g1 * g2 * cov_cd[ridx]

    # Numerical protection against tiny negative rounding artifacts.
    var_baco[var_baco < 0 & var_baco > -1e-12] <- 0

    valid_var <- is.finite(var_baco) & var_baco >= 0

    baco_se[ridx[valid_var]] <-
      sqrt(var_baco[valid_var])

    z95 <- qnorm(1 - alpha_main / 2)

    baco_ci_low[ridx[valid_var]] <-
      baco_hat[ridx[valid_var]] -
      z95 * baco_se[ridx[valid_var]]

    baco_ci_high[ridx[valid_var]] <-
      baco_hat[ridx[valid_var]] +
      z95 * baco_se[ridx[valid_var]]

    z_baco <-
      (baco_hat[ridx[valid_var]] - 1) /
      baco_se[ridx[valid_var]]

    p_baco_wald[ridx[valid_var]] <-
      2 * pnorm(-abs(z_baco))
  }

  # --------------------------------------------------------------------------
  # B) Direct contrast:
  #       Delta = beta_c - beta_d
  # --------------------------------------------------------------------------

  delta_hat <- rep(NA_real_, nsim)
  delta_se <- rep(NA_real_, nsim)
  p_delta <- rep(NA_real_, nsim)

  didx <- which(estimable_original)

  if (length(didx) > 0) {

    delta_hat[didx] <-
      beta_c[didx] -
      beta_d[didx]

    var_delta <-
      var_c[didx] +
      var_d[didx] -
      2 * cov_cd[didx]

    var_delta[var_delta < 0 & var_delta > -1e-12] <- 0

    valid_delta_var <-
      is.finite(var_delta) &
      var_delta >= 0

    delta_se[didx[valid_delta_var]] <-
      sqrt(var_delta[valid_delta_var])

    z_delta <-
      delta_hat[didx[valid_delta_var]] /
      delta_se[didx[valid_delta_var]]

    p_delta[didx[valid_delta_var]] <-
      2 * pnorm(-abs(z_delta))
  }

  # --------------------------------------------------------------------------
  # C) Fieller confidence set
  # --------------------------------------------------------------------------

  fieller <- fieller_vectorized(
    beta_c = beta_c,
    beta_d = beta_d,
    var_c = var_c,
    var_d = var_d,
    cov_cd = cov_cd,
    alpha = alpha_main
  )

  # This should agree with the direct contrast test at alpha = 0.05.
  fieller_reject_1 <-
    ifelse(
      is.na(fieller$fieller_contains_1),
      NA,
      !fieller$fieller_contains_1
    )

  # --------------------------------------------------------------------------
  # Denominator diagnostics
  # --------------------------------------------------------------------------

  z95 <- qnorm(1 - alpha_main / 2)

  beta_d_ci_low <-
    beta_d -
    z95 * sqrt(var_d)

  beta_d_ci_high <-
    beta_d +
    z95 * sqrt(var_d)

  denominator_ci_excludes_zero <-
    ifelse(
      is.finite(beta_d_ci_low) &
      is.finite(beta_d_ci_high),
      beta_d_ci_low > 0 |
      beta_d_ci_high < 0,
      NA
    )

  # --------------------------------------------------------------------------
  # Replicate-level data
  # --------------------------------------------------------------------------

  sim <- tibble(
    total_deaths = total_deaths,

    zero_death_any_arm = zero_death_any_arm,
    zero_death_both_arms = zero_death_both_arms,

    estimable_original = estimable_original,
    estimable_ratio = estimable_ratio,

    beta_c = beta_c,
    beta_d = beta_d,

    baco_hat = baco_hat,
    baco_se = baco_se,
    baco_ci_low = baco_ci_low,
    baco_ci_high = baco_ci_high,
    p_baco_wald = p_baco_wald,

    delta_hat = delta_hat,
    delta_se = delta_se,
    p_delta = p_delta,

    fieller_type = fieller$fieller_type,
    fieller_contains_1 = fieller$fieller_contains_1,
    fieller_reject_1 = fieller_reject_1,

    denominator_ci_excludes_zero =
      denominator_ci_excludes_zero
  )

  # --------------------------------------------------------------------------
  # Scenario-level summaries
  # --------------------------------------------------------------------------

  is_defined_baco_population <-
    s$scenario_class == "BACO_equals_1"

  denominator_estimable_n <-
    sum(!is.na(sim$denominator_ci_excludes_zero))

  baco_valid_n <-
    sum(!is.na(sim$p_baco_wald))

  delta_valid_n <-
    sum(!is.na(sim$p_delta))

  fieller_valid_n <-
    sum(!is.na(sim$fieller_contains_1))

  # Conditional rejection rates are evaluated among replicates for which the
  # corresponding method returned a valid result.
  baco_reject_05_cond <-
    if (baco_valid_n > 0) {
      mean(sim$p_baco_wald < alpha_main, na.rm = TRUE)
    } else {
      NA_real_
    }

  baco_reject_005_cond <-
    if (baco_valid_n > 0) {
      mean(sim$p_baco_wald < alpha_strict, na.rm = TRUE)
    } else {
      NA_real_
    }

  delta_reject_05_cond <-
    if (delta_valid_n > 0) {
      mean(sim$p_delta < alpha_main, na.rm = TRUE)
    } else {
      NA_real_
    }

  delta_reject_005_cond <-
    if (delta_valid_n > 0) {
      mean(sim$p_delta < alpha_strict, na.rm = TRUE)
    } else {
      NA_real_
    }

  fieller_reject_05_cond <-
    if (fieller_valid_n > 0) {
      mean(sim$fieller_reject_1, na.rm = TRUE)
    } else {
      NA_real_
    }

  # Unconditional rates count non-estimable replicates as non-rejections.
  # They are useful, but should NOT replace conditional calibration summaries.
  baco_reject_05_uncond <-
    sum(sim$p_baco_wald < alpha_main, na.rm = TRUE) /
    nsim

  delta_reject_05_uncond <-
    sum(sim$p_delta < alpha_main, na.rm = TRUE) /
    nsim

  fieller_reject_05_uncond <-
    sum(sim$fieller_reject_1, na.rm = TRUE) /
    nsim

  # Coverage is meaningful only when population BACO is defined (=1).
  baco_ci_coverage_cond <-
    if (is_defined_baco_population && baco_valid_n > 0) {
      mean(
        sim$baco_ci_low <= 1 &
        sim$baco_ci_high >= 1,
        na.rm = TRUE
      )
    } else {
      NA_real_
    }

  fieller_coverage_cond <-
    if (is_defined_baco_population && fieller_valid_n > 0) {
      mean(
        sim$fieller_contains_1,
        na.rm = TRUE
      )
    } else {
      NA_real_
    }

  # Agreement between Delta test and Fieller decision.
  agreement_idx <-
    !is.na(sim$p_delta) &
    !is.na(sim$fieller_reject_1)

  delta_fieller_agreement <-
    if (any(agreement_idx)) {
      mean(
        (sim$p_delta[agreement_idx] < alpha_main) ==
        sim$fieller_reject_1[agreement_idx]
      )
    } else {
      NA_real_
    }

  # Distributional instability of the BACO point estimate.
  finite_baco <- sim$baco_hat[is.finite(sim$baco_hat)]

  tibble(
    scenario_id = s$scenario_id,

    n_per_arm = s$n_per_arm,

    p_death_control = s$p_death_control,
    p_composite_control = s$p_composite_control,

    rr_reference = s$rr_reference,
    population_beta_reference = s$population_beta_reference,
    abs_population_beta_reference = s$abs_population_beta_reference,

    population_baco = s$population_baco,
    scenario_class = s$scenario_class,

    expected_total_deaths = s$expected_total_deaths,

    nsim = nsim,

    zero_death_any_arm_fraction =
      mean(sim$zero_death_any_arm),

    zero_death_both_arms_fraction =
      mean(sim$zero_death_both_arms),

    original_estimable_fraction =
      mean(sim$estimable_original),

    ratio_estimable_fraction =
      mean(sim$estimable_ratio),

    denominator_ci_excludes_zero_fraction =
      if (denominator_estimable_n > 0) {
        mean(
          sim$denominator_ci_excludes_zero,
          na.rm = TRUE
        )
      } else {
        NA_real_
      },

    fieller_bounded_fraction =
      mean(
        sim$fieller_type == "bounded",
        na.rm = TRUE
      ),

    fieller_disjoint_unbounded_fraction =
      mean(
        sim$fieller_type == "disjoint_unbounded",
        na.rm = TRUE
      ),

    fieller_all_real_fraction =
      mean(
        sim$fieller_type == "all_real",
        na.rm = TRUE
      ),

    baco_wald_type1_05_conditional =
      ifelse(
        is_defined_baco_population,
        baco_reject_05_cond,
        NA_real_
      ),

    baco_wald_type1_005_conditional =
      ifelse(
        is_defined_baco_population,
        baco_reject_005_cond,
        NA_real_
      ),

    baco_wald_rejection_05_unconditional =
      ifelse(
        is_defined_baco_population,
        baco_reject_05_uncond,
        NA_real_
      ),

    delta_type1_05_conditional =
      ifelse(
        is_defined_baco_population,
        delta_reject_05_cond,
        NA_real_
      ),

    delta_type1_005_conditional =
      ifelse(
        is_defined_baco_population,
        delta_reject_005_cond,
        NA_real_
      ),

    delta_rejection_05_unconditional =
      ifelse(
        is_defined_baco_population,
        delta_reject_05_uncond,
        NA_real_
      ),

    fieller_type1_05_conditional =
      ifelse(
        is_defined_baco_population,
        fieller_reject_05_cond,
        NA_real_
      ),

    fieller_rejection_05_unconditional =
      ifelse(
        is_defined_baco_population,
        fieller_reject_05_uncond,
        NA_real_
      ),

    baco_delta_ci_coverage_95_conditional =
      baco_ci_coverage_cond,

    fieller_coverage_95_conditional =
      fieller_coverage_cond,

    delta_fieller_decision_agreement =
      delta_fieller_agreement,

    median_baco_hat =
      if (length(finite_baco) > 0) {
        median(finite_baco)
      } else {
        NA_real_
      },

    q025_baco_hat =
      if (length(finite_baco) > 0) {
        quantile(
          finite_baco,
          0.025,
          names = FALSE
        )
      } else {
        NA_real_
      },

    q975_baco_hat =
      if (length(finite_baco) > 0) {
        quantile(
          finite_baco,
          0.975,
          names = FALSE
        )
      } else {
        NA_real_
      },

    extreme_baco_abs_gt_10_fraction =
      if (length(finite_baco) > 0) {
        mean(abs(finite_baco) > 10)
      } else {
        NA_real_
      },

    extreme_baco_abs_gt_100_fraction =
      if (length(finite_baco) > 0) {
        mean(abs(finite_baco) > 100)
      } else {
        NA_real_
      }
  )
}


# ------------------------------------------------------------------------------
# 5. Run all scenarios
# ------------------------------------------------------------------------------

cat("\nRunning Monte Carlo calibration...\n")
cat("Scenarios: ", nrow(scenario_grid), "\n", sep = "")
cat("Replicates per scenario: ", NSIM, "\n\n", sep = "")

simulation_results <- map_dfr(
  seq_len(nrow(scenario_grid)),
  function(i) {

    if (
      i == 1 ||
      i %% 25 == 0 ||
      i == nrow(scenario_grid)
    ) {
      cat(
        "Scenario ",
        i,
        " / ",
        nrow(scenario_grid),
        "\n",
        sep = ""
      )
    }

    simulate_one_scenario(
      scenario_row = scenario_grid[i, ],
      nsim = NSIM,
      alpha_main = ALPHA_MAIN,
      alpha_strict = ALPHA_STRICT,
      base_seed = SEED
    )
  }
)

write_csv(
  simulation_results,
  file.path(
    OUTPUT_DIR,
    "01_type1_calibration_all_scenarios.csv"
  )
)


# ------------------------------------------------------------------------------
# 6. Primary calibration summaries
# ------------------------------------------------------------------------------

defined_results <- simulation_results %>%
  filter(
    scenario_class == "BACO_equals_1"
  )

overall_summary <- defined_results %>%
  summarise(
    scenarios = n(),

    median_original_estimable_fraction =
      median(
        original_estimable_fraction,
        na.rm = TRUE
      ),

    median_ratio_estimable_fraction =
      median(
        ratio_estimable_fraction,
        na.rm = TRUE
      ),

    median_type1_baco_wald_05 =
      median(
        baco_wald_type1_05_conditional,
        na.rm = TRUE
      ),

    median_type1_delta_05 =
      median(
        delta_type1_05_conditional,
        na.rm = TRUE
      ),

    median_type1_fieller_05 =
      median(
        fieller_type1_05_conditional,
        na.rm = TRUE
      ),

    max_type1_baco_wald_05 =
      max(
        baco_wald_type1_05_conditional,
        na.rm = TRUE
      ),

    max_type1_delta_05 =
      max(
        delta_type1_05_conditional,
        na.rm = TRUE
      ),

    min_delta_fieller_agreement =
      min(
        delta_fieller_decision_agreement,
        na.rm = TRUE
      )
  )

write_csv(
  overall_summary,
  file.path(
    OUTPUT_DIR,
    "02_overall_calibration_summary.csv"
  )
)


# ------------------------------------------------------------------------------
# 7. Identify the most problematic scenarios for the original BACO Wald test
# ------------------------------------------------------------------------------

worst_wald_scenarios <- defined_results %>%
  arrange(
    desc(baco_wald_type1_05_conditional)
  ) %>%
  select(
    scenario_id,
    n_per_arm,
    p_death_control,
    p_composite_control,
    rr_reference,
    abs_population_beta_reference,
    expected_total_deaths,

    original_estimable_fraction,
    ratio_estimable_fraction,

    baco_wald_type1_05_conditional,
    delta_type1_05_conditional,
    fieller_type1_05_conditional,

    baco_delta_ci_coverage_95_conditional,
    fieller_coverage_95_conditional,

    fieller_bounded_fraction,
    denominator_ci_excludes_zero_fraction,

    extreme_baco_abs_gt_10_fraction,
    extreme_baco_abs_gt_100_fraction
  ) %>%
  slice_head(n = 30)

write_csv(
  worst_wald_scenarios,
  file.path(
    OUTPUT_DIR,
    "03_worst_Wald_calibration_scenarios.csv"
  )
)


# ------------------------------------------------------------------------------
# 8. Structural RR = 1 scenarios: BACO is undefined
# ------------------------------------------------------------------------------

reference_null_results <- simulation_results %>%
  filter(
    scenario_class == "reference_null_BACO_undefined"
  ) %>%
  select(
    scenario_id,
    n_per_arm,
    p_death_control,
    p_composite_control,
    rr_reference,
    expected_total_deaths,

    original_estimable_fraction,
    ratio_estimable_fraction,

    fieller_bounded_fraction,
    fieller_disjoint_unbounded_fraction,
    fieller_all_real_fraction,

    denominator_ci_excludes_zero_fraction,

    median_baco_hat,
    q025_baco_hat,
    q975_baco_hat,

    extreme_baco_abs_gt_10_fraction,
    extreme_baco_abs_gt_100_fraction
  )

write_csv(
  reference_null_results,
  file.path(
    OUTPUT_DIR,
    "04_reference_null_RR_equals_1_BACO_undefined.csv"
  )
)


# ------------------------------------------------------------------------------
# 9. Monte Carlo uncertainty
# ------------------------------------------------------------------------------

mcse_05 <- sqrt(
  ALPHA_MAIN *
  (1 - ALPHA_MAIN) /
  NSIM
)

mcse_005 <- sqrt(
  ALPHA_STRICT *
  (1 - ALPHA_STRICT) /
  NSIM
)

mcse_info <- tibble(
  nsim = NSIM,
  nominal_alpha = c(
    ALPHA_MAIN,
    ALPHA_STRICT
  ),
  approximate_Monte_Carlo_SE = c(
    mcse_05,
    mcse_005
  )
)

write_csv(
  mcse_info,
  file.path(
    OUTPUT_DIR,
    "05_Monte_Carlo_standard_error.csv"
  )
)


# ------------------------------------------------------------------------------
# 10. Plot-ready long format
# ------------------------------------------------------------------------------

plot_long <- defined_results %>%
  select(
    scenario_id,
    n_per_arm,
    p_death_control,
    p_composite_control,
    rr_reference,
    abs_population_beta_reference,
    expected_total_deaths,

    baco_wald_type1_05_conditional,
    delta_type1_05_conditional,
    fieller_type1_05_conditional
  ) %>%
  pivot_longer(
    cols = c(
      baco_wald_type1_05_conditional,
      delta_type1_05_conditional,
      fieller_type1_05_conditional
    ),
    names_to = "method",
    values_to = "empirical_type1"
  ) %>%
  mutate(
    method = recode(
      method,
      baco_wald_type1_05_conditional =
        "Original BACO Wald/delta",
      delta_type1_05_conditional =
        "Direct log-effect discrepancy",
      fieller_type1_05_conditional =
        "Fieller"
    )
  )

write_csv(
  plot_long,
  file.path(
    OUTPUT_DIR,
    "06_plot_ready_type1_long.csv"
  )
)


# ------------------------------------------------------------------------------
# 11. Figure A:
#     Type-I error vs distance of reference effect from null
#
# To avoid an unreadable figure, select a representative design:
# p_composite_control = 0.25 and n = 500.
# The full scenario table remains available in the CSV.
# ------------------------------------------------------------------------------

figA_data <- plot_long %>%
  filter(
    n_per_arm == 500,
    p_composite_control == 0.25
  )

pA <- ggplot(
  figA_data,
  aes(
    x = abs_population_beta_reference,
    y = empirical_type1,
    shape = method
  )
) +
  geom_hline(
    yintercept = ALPHA_MAIN,
    linetype = 2
  ) +
  geom_point(size = 2.5) +
  geom_line(
    aes(group = interaction(method, p_death_control)),
    linewidth = 0.5
  ) +
  facet_wrap(
    ~ p_death_control,
    labeller = label_bquote(
      p[D0] == .(p_death_control)
    )
  ) +
  labs(
    x = expression(abs(log(RR[d]))),
    y = "Empirical rejection probability",
    shape = NULL,
    title = "Calibration of BACO inference as the reference effect approaches the null",
    subtitle = paste0(
      "BACO = 1 in the population; ",
      NSIM,
      " Monte Carlo replicates per scenario"
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom"
  )

ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "Figure_A_type1_vs_reference_effect.png"
  ),
  plot = pA,
  width = 11,
  height = 7,
  dpi = 300
)


# ------------------------------------------------------------------------------
# 12. Figure B:
#     Original BACO Wald type-I error by expected death count and weak denominator
# ------------------------------------------------------------------------------

figB_data <- defined_results %>%
  filter(
    p_composite_control == 0.25
  )

pB <- ggplot(
  figB_data,
  aes(
    x = expected_total_deaths,
    y = abs_population_beta_reference,
    size = baco_wald_type1_05_conditional,
    shape = n_per_arm
  )
) +
  geom_point(alpha = 0.75) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "Expected total critical events (log scale)",
    y = expression(abs(log(RR[d])) ~ "(log scale)"),
    size = "Wald rejection\nprobability",
    shape = "N per arm",
    title = "Original BACO Wald calibration depends on both event information and denominator strength",
    subtitle = "Population BACO = 1"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom"
  )

ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "Figure_B_Wald_calibration_event_information.png"
  ),
  plot = pB,
  width = 10,
  height = 7,
  dpi = 300
)


# ------------------------------------------------------------------------------
# 13. Figure C:
#     Non-estimability from sparse critical events
# ------------------------------------------------------------------------------

figC_data <- defined_results %>%
  filter(
    p_composite_control == 0.25,
    rr_reference %in% c(
      0.50,
      0.95,
      1.05,
      2.00
    )
  )

pC <- ggplot(
  figC_data,
  aes(
    x = expected_total_deaths,
    y = ratio_estimable_fraction,
    shape = factor(rr_reference)
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = 2
  ) +
  geom_point(size = 2.5) +
  geom_line(
    aes(
      group = interaction(
        rr_reference,
        n_per_arm
      )
    ),
    linewidth = 0.5
  ) +
  scale_x_log10() +
  labs(
    x = "Expected total critical events (log scale)",
    y = "Fraction of replicates with calculable BACO",
    shape = expression(RR[d]),
    title = "Sparse critical events make the original BACO non-estimable"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom"
  )

ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "Figure_C_BACO_estimability_sparse_events.png"
  ),
  plot = pC,
  width = 9,
  height = 6,
  dpi = 300
)


# ------------------------------------------------------------------------------
# 14. Sanity checks
# ------------------------------------------------------------------------------

# Fieller decision for theta = 1 and the direct Delta test should be nearly
# identical under the same normal-theory covariance structure.

agreement_check <- defined_results %>%
  summarise(
    min_agreement =
      min(
        delta_fieller_decision_agreement,
        na.rm = TRUE
      ),

    median_agreement =
      median(
        delta_fieller_decision_agreement,
        na.rm = TRUE
      ),

    max_absolute_difference_type1 =
      max(
        abs(
          delta_type1_05_conditional -
          fieller_type1_05_conditional
        ),
        na.rm = TRUE
      )
  )

write_csv(
  agreement_check,
  file.path(
    OUTPUT_DIR,
    "07_Delta_Fieller_agreement_check.csv"
  )
)


# ------------------------------------------------------------------------------
# 15. Console output
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("BACO Monte Carlo calibration completed\n")
cat("============================================================\n\n")

cat("Monte Carlo SE at alpha = 0.05: ")
print(mcse_05)

cat("\nMonte Carlo SE at alpha = 0.005: ")
print(mcse_005)

cat("\nOverall calibration summary:\n")
print(overall_summary)

cat("\nDelta vs Fieller agreement:\n")
print(agreement_check)

cat("\nWorst original-Wald scenarios:\n")
print(worst_wald_scenarios %>% slice_head(n = 10))

cat(
  "\nOutputs saved in: ",
  normalizePath(
    OUTPUT_DIR,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)


# ==============================================================================
# INTERPRETATION RULES FOR THIS PILOT
# ==============================================================================
#
# 1) For scenario_class == "BACO_equals_1":
#    - empirical type-I error should be close to nominal alpha;
#    - 95% CI coverage should be close to 0.95.
#
# 2) RR_reference close to 1 is the critical stress region because:
#
#       beta_d = log(RR_d) -> 0.
#
# 3) RR_reference == 1 is NOT a type-I error scenario for BACO because:
#
#       BACO = 0 / 0
#
#    is undefined at the population level.
#    These scenarios are retained only to show the behavior of the sample
#    BACO estimator under a structurally undefined target.
#
# 4) Sparse scenarios must be evaluated using BOTH:
#    - conditional calibration among estimable replicates; and
#    - the fraction of replicates in which the original BACO is estimable.
#
# 5) Do not interpret a well-calibrated conditional test as a complete solution
#    if most simulated datasets cannot produce the estimator.
#
# 6) Bootstrap is deliberately NOT nested inside this first Monte Carlo grid.
#    A nested bootstrap would be computationally expensive and its sparse-cell
#    failure mechanism is already visible empirically. After identifying the
#    critical scenarios, a second targeted simulation can compare percentile,
#    BCa, parametric bootstrap, profile likelihood, and Bayesian inference.
#
# ==============================================================================
