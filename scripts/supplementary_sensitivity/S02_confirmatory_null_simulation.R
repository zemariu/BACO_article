# ==============================================================================
# BACO reproducibility analysis — supplementary sensitivity
# Confirmatory Monte Carlo simulation
#
# Goal
# ----
# Separate four inferential situations under the population null hypothesis
#
#       BACO = 1
#
# into a deliberately controlled 2 x 2 design:
#
#                         Reference effect
#                    strong              weak
#                --------------------------------
# many deaths     | regular         | weak denominator
# few deaths      | sparse          | combined
#
# This confirmatory simulation is designed to distinguish:
#
#   1) sparse-event failure;
#   2) weak-denominator failure;
#   3) their combination;
#   4) regular asymptotic behavior.
#
# Main methods compared
# ---------------------
# A) Original BACO Wald/delta inference:
#
#       theta_hat = beta_c / beta_d
#
# B) Direct discrepancy test:
#
#       Delta = beta_c - beta_d
#
# C) Fieller confidence-set inference for theta = beta_c / beta_d.
#
# Under H0, whenever beta_d != 0:
#
#       RR_c = RR_d = R
#       BACO = 1
#       Delta = 0
#
# A separate RR = 1 block is retained only to illustrate the structural
# 0/0 problem. It is NOT included in type-I error calculations.
#
# Recommended main run:
#       NSIM <- 10000
#
# Optional final run after all methods are frozen:
#       NSIM <- 20000 or 50000
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
    "Install missing packages first: ",
    paste(missing_packages, collapse = ", ")
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
# 1. Global settings
# ------------------------------------------------------------------------------

ROOT_DIR <- getwd()
OUTPUT_DIR <- file.path(ROOT_DIR, "results", "supplementary_confirmatory_null")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

SEED <- 20260902

# Main confirmatory setting.
NSIM <- 10000

ALPHA <- 0.05
ALPHA_STRICT <- 0.005

# A large, common denominator is used deliberately so that event information
# is controlled primarily through the EXPECTED NUMBER OF CRITICAL EVENTS,
# rather than through sample size itself.
#
# Because simulation is performed with multinomial counts, N = 10,000 per arm
# is computationally inexpensive.
N_PER_ARM <- 10000

P_COMPOSITE_CONTROL <- c(
  0.10,
  0.25
)

# Reference effects clearly separated from the null.
RR_STRONG <- c(
  0.50,
  0.70,
  1.50,
  2.00
)

# Reference effects deliberately close to 1.
#
# These values are design points, not proposed diagnostic cutoffs.
RR_WEAK <- c(
  0.90,
  0.95,
  0.975,
  0.99,
  1.01,
  1.025,
  1.05,
  1.10
)

# Expected total critical events across both trial arms.
DEATHS_SPARSE <- c(
  1,
  2,
  3,
  5,
  10
)

DEATHS_INFORMATIVE <- c(
  50,
  100,
  250,
  500
)

# Additional structural block where BACO is undefined in the population.
DEATHS_REFERENCE_NULL <- c(
  5,
  50,
  500
)

NUMERIC_TOL <- sqrt(.Machine$double.eps)


# ------------------------------------------------------------------------------
# 2. Scenario construction
# ------------------------------------------------------------------------------

# Given equal allocation and a target expected total number of deaths:
#
#   E[D_total] = N*p_d0 + N*p_d1
#              = N*p_d0*(1 + RR_d)
#
# therefore:
#
#   p_d0 = E[D_total] / {N*(1 + RR_d)}
#
# Under BACO = 1:
#
#   RR_c = RR_d = R
#
# and:
#
#   p_c1 = R*p_c0
#   p_d1 = R*p_d0

make_block <- function(
    class_name,
    rr_values,
    expected_deaths_values,
    p_composite_values = P_COMPOSITE_CONTROL,
    n_per_arm = N_PER_ARM
) {

  tidyr::crossing(
    scenario_class = class_name,
    rr_reference = rr_values,
    expected_total_deaths_target = expected_deaths_values,
    p_composite_control = p_composite_values,
    n_per_arm = n_per_arm
  ) %>%
    mutate(
      p_death_control =
        expected_total_deaths_target /
        (n_per_arm * (1 + rr_reference)),

      p_death_experimental =
        rr_reference * p_death_control,

      p_composite_experimental =
        rr_reference * p_composite_control,

      population_beta_reference =
        log(rr_reference),

      abs_population_beta_reference =
        abs(population_beta_reference),

      population_baco =
        ifelse(
          rr_reference == 1,
          NA_real_,
          1
        ),

      valid_probability_structure =
        p_death_control >= 0 &
        p_death_control <= p_composite_control &
        p_death_experimental >= 0 &
        p_death_experimental <= p_composite_experimental &
        p_composite_experimental <= 1
    ) %>%
    filter(valid_probability_structure)
}


regular_block <- make_block(
  class_name = "regular",
  rr_values = RR_STRONG,
  expected_deaths_values = DEATHS_INFORMATIVE
)

weak_block <- make_block(
  class_name = "weak_denominator",
  rr_values = RR_WEAK,
  expected_deaths_values = DEATHS_INFORMATIVE
)

sparse_block <- make_block(
  class_name = "sparse",
  rr_values = RR_STRONG,
  expected_deaths_values = DEATHS_SPARSE
)

combined_block <- make_block(
  class_name = "combined_sparse_weak",
  rr_values = RR_WEAK,
  expected_deaths_values = DEATHS_SPARSE
)

reference_null_block <- make_block(
  class_name = "reference_null_BACO_undefined",
  rr_values = 1,
  expected_deaths_values = DEATHS_REFERENCE_NULL
)

scenario_grid <- bind_rows(
  regular_block,
  weak_block,
  sparse_block,
  combined_block,
  reference_null_block
) %>%
  arrange(
    scenario_class,
    p_composite_control,
    expected_total_deaths_target,
    rr_reference
  ) %>%
  mutate(
    scenario_id = row_number()
  ) %>%
  select(
    scenario_id,
    everything()
  )

write_csv(
  scenario_grid,
  file.path(
    OUTPUT_DIR,
    "00_confirmatory_scenario_grid.csv"
  )
)

cat(
  "\nConfirmatory scenarios: ",
  nrow(scenario_grid),
  "\n",
  sep = ""
)

print(
  scenario_grid %>%
    count(
      scenario_class,
      name = "n_scenarios"
    )
)


# ------------------------------------------------------------------------------
# 3. Vectorized Fieller inference
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
  B <- -2 * beta_c * beta_d +
    2 * z^2 * cov_cd
  C <- beta_c^2 - z^2 * var_c

  discriminant <- B^2 - 4 * A * C

  finite_input <-
    is.finite(beta_c) &
    is.finite(beta_d) &
    is.finite(var_c) &
    is.finite(var_d) &
    is.finite(cov_cd)

  fieller_type <- rep(
    "not_available",
    length(beta_c)
  )

  # Standard quadratic cases.
  fieller_type[
    finite_input &
    abs(A) > tol &
    discriminant >= 0 &
    A > 0
  ] <- "bounded"

  fieller_type[
    finite_input &
    abs(A) > tol &
    discriminant >= 0 &
    A < 0
  ] <- "disjoint_unbounded"

  fieller_type[
    finite_input &
    abs(A) > tol &
    discriminant < 0 &
    A < 0
  ] <- "all_real"

  fieller_type[
    finite_input &
    abs(A) > tol &
    discriminant < 0 &
    A > 0
  ] <- "empty"

  # Degenerate/linear cases.
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

  # Direct test of whether theta = 1 belongs to the Fieller set.
  #
  # The Fieller quadratic inequality evaluated at theta=1 is:
  #
  #   A + B + C <= 0
  #
  # which is algebraically equivalent to testing:
  #
  #   beta_c - beta_d = 0
  #
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
# 4. Monte Carlo helper functions
# ------------------------------------------------------------------------------

safe_rate <- function(x) {

  n_valid <- sum(!is.na(x))

  if (n_valid == 0) {
    return(
      c(
        rate = NA_real_,
        n_valid = 0,
        mcse = NA_real_,
        mc_low = NA_real_,
        mc_high = NA_real_
      )
    )
  }

  p <- mean(x, na.rm = TRUE)

  mcse <- sqrt(
    p * (1 - p) /
    n_valid
  )

  c(
    rate = p,
    n_valid = n_valid,
    mcse = mcse,
    mc_low = max(0, p - 1.96 * mcse),
    mc_high = min(1, p + 1.96 * mcse)
  )
}


# ------------------------------------------------------------------------------
# 5. Simulate one confirmatory scenario
# ------------------------------------------------------------------------------

simulate_confirmatory_scenario <- function(
    s,
    nsim = NSIM,
    alpha = ALPHA,
    alpha_strict = ALPHA_STRICT,
    base_seed = SEED
) {

  n_t <- s$n_per_arm
  n_c <- s$n_per_arm

  # Three mutually exclusive states:
  #   none / nonfatal composite / death
  prob_c <- c(
    none =
      1 - s$p_composite_control,

    nonfatal_composite =
      s$p_composite_control -
      s$p_death_control,

    death =
      s$p_death_control
  )

  prob_t <- c(
    none =
      1 - s$p_composite_experimental,

    nonfatal_composite =
      s$p_composite_experimental -
      s$p_death_experimental,

    death =
      s$p_death_experimental
  )

  if (
    any(prob_c < 0) ||
    any(prob_t < 0) ||
    abs(sum(prob_c) - 1) > 1e-10 ||
    abs(sum(prob_t) - 1) > 1e-10
  ) {
    stop(
      "Invalid probability structure in scenario ",
      s$scenario_id
    )
  }

  set.seed(
    base_seed +
    s$scenario_id * 1009L
  )

  draw_c <- rmultinom(
    n = nsim,
    size = n_c,
    prob = prob_c
  )

  draw_t <- rmultinom(
    n = nsim,
    size = n_t,
    prob = prob_t
  )

  composite_c <-
    draw_c["nonfatal_composite", ] +
    draw_c["death", ]

  composite_t <-
    draw_t["nonfatal_composite", ] +
    draw_t["death", ]

  death_c <- draw_c["death", ]
  death_t <- draw_t["death", ]

  total_deaths <- death_c + death_t

  zero_death_any_arm <-
    death_c == 0 |
    death_t == 0

  zero_death_both_arms <-
    death_c == 0 &
    death_t == 0

  zero_composite_any_arm <-
    composite_c == 0 |
    composite_t == 0

  original_estimable <-
    !zero_death_any_arm &
    !zero_composite_any_arm

  rr_c <- rep(NA_real_, nsim)
  rr_d <- rep(NA_real_, nsim)

  beta_c <- rep(NA_real_, nsim)
  beta_d <- rep(NA_real_, nsim)

  var_c <- rep(NA_real_, nsim)
  var_d <- rep(NA_real_, nsim)
  cov_cd <- rep(NA_real_, nsim)

  idx <- which(original_estimable)

  if (length(idx) > 0) {

    rr_c[idx] <-
      (composite_t[idx] / n_t) /
      (composite_c[idx] / n_c)

    rr_d[idx] <-
      (death_t[idx] / n_t) /
      (death_c[idx] / n_c)

    beta_c[idx] <- log(rr_c[idx])
    beta_d[idx] <- log(rr_d[idx])

    var_c[idx] <-
      1 / composite_t[idx] +
      1 / composite_c[idx] -
      1 / n_t -
      1 / n_c

    var_d[idx] <-
      1 / death_t[idx] +
      1 / death_c[idx] -
      1 / n_t -
      1 / n_c

    # Original BACO covariance derivation.
    cov_cd[idx] <- var_c[idx]
  }

  # --------------------------------------------------------------------------
  # Original BACO Wald/delta
  # --------------------------------------------------------------------------

  ratio_estimable <-
    original_estimable &
    is.finite(beta_d) &
    abs(beta_d) > NUMERIC_TOL

  baco_hat <- rep(NA_real_, nsim)
  baco_se <- rep(NA_real_, nsim)
  p_baco <- rep(NA_real_, nsim)

  baco_ci_low <- rep(NA_real_, nsim)
  baco_ci_high <- rep(NA_real_, nsim)

  ridx <- which(ratio_estimable)

  if (length(ridx) > 0) {

    baco_hat[ridx] <-
      beta_c[ridx] /
      beta_d[ridx]

    g1 <- 1 / beta_d[ridx]

    g2 <-
      -beta_c[ridx] /
      beta_d[ridx]^2

    var_baco <-
      g1^2 * var_c[ridx] +
      g2^2 * var_d[ridx] +
      2 * g1 * g2 * cov_cd[ridx]

    # Remove only tiny numerical negatives.
    var_baco[
      var_baco < 0 &
      var_baco > -1e-12
    ] <- 0

    good_var <-
      is.finite(var_baco) &
      var_baco >= 0

    valid_idx <- ridx[good_var]

    baco_se[valid_idx] <-
      sqrt(var_baco[good_var])

    zcrit <- qnorm(
      1 - alpha / 2
    )

    baco_ci_low[valid_idx] <-
      baco_hat[valid_idx] -
      zcrit * baco_se[valid_idx]

    baco_ci_high[valid_idx] <-
      baco_hat[valid_idx] +
      zcrit * baco_se[valid_idx]

    z_baco <-
      (baco_hat[valid_idx] - 1) /
      baco_se[valid_idx]

    p_baco[valid_idx] <-
      2 * pnorm(
        -abs(z_baco)
      )
  }

  # --------------------------------------------------------------------------
  # Direct discrepancy:
  # Delta = beta_c - beta_d
  # --------------------------------------------------------------------------

  delta_hat <- rep(
    NA_real_,
    nsim
  )

  delta_se <- rep(
    NA_real_,
    nsim
  )

  p_delta <- rep(
    NA_real_,
    nsim
  )

  didx <- which(
    original_estimable
  )

  if (length(didx) > 0) {

    delta_hat[didx] <-
      beta_c[didx] -
      beta_d[didx]

    var_delta <-
      var_c[didx] +
      var_d[didx] -
      2 * cov_cd[didx]

    var_delta[
      var_delta < 0 &
      var_delta > -1e-12
    ] <- 0

    good_delta <-
      is.finite(var_delta) &
      var_delta >= 0

    valid_delta_idx <-
      didx[good_delta]

    delta_se[valid_delta_idx] <-
      sqrt(
        var_delta[good_delta]
      )

    z_delta <-
      delta_hat[valid_delta_idx] /
      delta_se[valid_delta_idx]

    p_delta[valid_delta_idx] <-
      2 * pnorm(
        -abs(z_delta)
      )
  }

  # --------------------------------------------------------------------------
  # Fieller
  # --------------------------------------------------------------------------

  fieller <- fieller_vectorized(
    beta_c = beta_c,
    beta_d = beta_d,
    var_c = var_c,
    var_d = var_d,
    cov_cd = cov_cd,
    alpha = alpha
  )

  fieller_reject_1 <-
    ifelse(
      is.na(
        fieller$fieller_contains_1
      ),
      NA,
      !fieller$fieller_contains_1
    )

  # --------------------------------------------------------------------------
  # Denominator strength diagnostics
  # --------------------------------------------------------------------------

  zcrit <- qnorm(
    1 - alpha / 2
  )

  beta_d_se <-
    sqrt(var_d)

  beta_d_ci_low <-
    beta_d -
    zcrit * beta_d_se

  beta_d_ci_high <-
    beta_d +
    zcrit * beta_d_se

  denominator_ci_excludes_zero <-
    ifelse(
      is.finite(beta_d_ci_low) &
      is.finite(beta_d_ci_high),

      beta_d_ci_low > 0 |
      beta_d_ci_high < 0,

      NA
    )

  # --------------------------------------------------------------------------
  # Replicate-level table
  # --------------------------------------------------------------------------

  sim <- tibble(
    total_deaths = total_deaths,

    zero_death_any_arm =
      zero_death_any_arm,

    zero_death_both_arms =
      zero_death_both_arms,

    original_estimable =
      original_estimable,

    ratio_estimable =
      ratio_estimable,

    beta_c = beta_c,
    beta_d = beta_d,

    baco_hat = baco_hat,
    baco_se = baco_se,
    baco_ci_low = baco_ci_low,
    baco_ci_high = baco_ci_high,
    p_baco = p_baco,

    delta_hat = delta_hat,
    delta_se = delta_se,
    p_delta = p_delta,

    fieller_type =
      fieller$fieller_type,

    fieller_contains_1 =
      fieller$fieller_contains_1,

    fieller_reject_1 =
      fieller_reject_1,

    denominator_ci_excludes_zero =
      denominator_ci_excludes_zero
  )

  # --------------------------------------------------------------------------
  # Inferential summaries
  # --------------------------------------------------------------------------

  population_baco_defined <-
    s$rr_reference != 1

  # Conditional type-I error among replicates in which each method is estimable.
  baco_reject_05 <-
    safe_rate(
      ifelse(
        is.na(sim$p_baco),
        NA,
        sim$p_baco < alpha
      )
    )

  baco_reject_005 <-
    safe_rate(
      ifelse(
        is.na(sim$p_baco),
        NA,
        sim$p_baco < alpha_strict
      )
    )

  delta_reject_05 <-
    safe_rate(
      ifelse(
        is.na(sim$p_delta),
        NA,
        sim$p_delta < alpha
      )
    )

  delta_reject_005 <-
    safe_rate(
      ifelse(
        is.na(sim$p_delta),
        NA,
        sim$p_delta < alpha_strict
      )
    )

  fieller_reject_05 <-
    safe_rate(
      sim$fieller_reject_1
    )

  # Coverage.
  baco_coverage <-
    safe_rate(
      ifelse(
        is.na(sim$baco_ci_low) |
        is.na(sim$baco_ci_high),

        NA,

        sim$baco_ci_low <= 1 &
        sim$baco_ci_high >= 1
      )
    )

  fieller_coverage <-
    safe_rate(
      sim$fieller_contains_1
    )

  # Agreement between the direct Delta test and Fieller at theta=1.
  agreement_idx <-
    !is.na(sim$p_delta) &
    !is.na(sim$fieller_reject_1)

  delta_fieller_agreement <-
    if (any(agreement_idx)) {
      mean(
        (sim$p_delta[agreement_idx] < alpha) ==
        sim$fieller_reject_1[agreement_idx]
      )
    } else {
      NA_real_
    }

  finite_baco <-
    sim$baco_hat[
      is.finite(sim$baco_hat)
    ]

  # --------------------------------------------------------------------------
  # Return one scenario-level row
  # --------------------------------------------------------------------------

  tibble(
    scenario_id =
      s$scenario_id,

    scenario_class =
      s$scenario_class,

    n_per_arm =
      s$n_per_arm,

    p_composite_control =
      s$p_composite_control,

    p_death_control =
      s$p_death_control,

    rr_reference =
      s$rr_reference,

    population_beta_reference =
      s$population_beta_reference,

    abs_population_beta_reference =
      s$abs_population_beta_reference,

    expected_total_deaths_target =
      s$expected_total_deaths_target,

    nsim = nsim,

    zero_death_any_arm_fraction =
      mean(sim$zero_death_any_arm),

    zero_death_both_arms_fraction =
      mean(sim$zero_death_both_arms),

    original_estimable_fraction =
      mean(sim$original_estimable),

    ratio_estimable_fraction =
      mean(sim$ratio_estimable),

    denominator_ci_excludes_zero_fraction =
      mean(
        sim$denominator_ci_excludes_zero,
        na.rm = TRUE
      ),

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

    # Original BACO Wald/delta.
    baco_wald_type1_05_conditional =
      ifelse(
        population_baco_defined,
        baco_reject_05["rate"],
        NA_real_
      ),

    baco_wald_type1_05_mcse =
      ifelse(
        population_baco_defined,
        baco_reject_05["mcse"],
        NA_real_
      ),

    baco_wald_type1_05_mc_low =
      ifelse(
        population_baco_defined,
        baco_reject_05["mc_low"],
        NA_real_
      ),

    baco_wald_type1_05_mc_high =
      ifelse(
        population_baco_defined,
        baco_reject_05["mc_high"],
        NA_real_
      ),

    baco_wald_type1_005_conditional =
      ifelse(
        population_baco_defined,
        baco_reject_005["rate"],
        NA_real_
      ),

    # Delta.
    delta_type1_05_conditional =
      ifelse(
        population_baco_defined,
        delta_reject_05["rate"],
        NA_real_
      ),

    delta_type1_05_mcse =
      ifelse(
        population_baco_defined,
        delta_reject_05["mcse"],
        NA_real_
      ),

    delta_type1_005_conditional =
      ifelse(
        population_baco_defined,
        delta_reject_005["rate"],
        NA_real_
      ),

    # Fieller.
    fieller_type1_05_conditional =
      ifelse(
        population_baco_defined,
        fieller_reject_05["rate"],
        NA_real_
      ),

    fieller_type1_05_mcse =
      ifelse(
        population_baco_defined,
        fieller_reject_05["mcse"],
        NA_real_
      ),

    # Coverage.
    baco_wald_coverage_95 =
      ifelse(
        population_baco_defined,
        baco_coverage["rate"],
        NA_real_
      ),

    fieller_coverage_95 =
      ifelse(
        population_baco_defined,
        fieller_coverage["rate"],
        NA_real_
      ),

    # Internal validation.
    delta_fieller_decision_agreement =
      delta_fieller_agreement,

    # BACO point-estimate instability.
    median_baco_hat =
      if (
        length(finite_baco) > 0
      ) {
        median(finite_baco)
      } else {
        NA_real_
      },

    q025_baco_hat =
      if (
        length(finite_baco) > 0
      ) {
        quantile(
          finite_baco,
          0.025,
          names = FALSE
        )
      } else {
        NA_real_
      },

    q975_baco_hat =
      if (
        length(finite_baco) > 0
      ) {
        quantile(
          finite_baco,
          0.975,
          names = FALSE
        )
      } else {
        NA_real_
      },

    abs_BACO_gt_10_fraction =
      if (
        length(finite_baco) > 0
      ) {
        mean(
          abs(finite_baco) > 10
        )
      } else {
        NA_real_
      },

    abs_BACO_gt_100_fraction =
      if (
        length(finite_baco) > 0
      ) {
        mean(
          abs(finite_baco) > 100
        )
      } else {
        NA_real_
      }
  )
}


# ------------------------------------------------------------------------------
# 6. Run confirmatory simulation
# ------------------------------------------------------------------------------

cat(
  "\nRunning confirmatory simulation with ",
  NSIM,
  " replicates/scenario...\n\n",
  sep = ""
)

results <- map_dfr(
  seq_len(
    nrow(scenario_grid)
  ),
  function(i) {

    if (
      i == 1 ||
      i %% 20 == 0 ||
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

    simulate_confirmatory_scenario(
      s = scenario_grid[i, ],
      nsim = NSIM,
      alpha = ALPHA,
      alpha_strict = ALPHA_STRICT,
      base_seed = SEED
    )
  }
)

write_csv(
  results,
  file.path(
    OUTPUT_DIR,
    "01_confirmatory_all_scenarios.csv"
  )
)


# ------------------------------------------------------------------------------
# 7. Class-level summary
# ------------------------------------------------------------------------------

defined_results <- results %>%
  filter(
    scenario_class !=
      "reference_null_BACO_undefined"
  )

class_summary <- defined_results %>%
  group_by(
    scenario_class
  ) %>%
  summarise(
    n_scenarios = n(),

    median_expected_total_deaths =
      median(
        expected_total_deaths_target
      ),

    median_abs_log_RRd =
      median(
        abs_population_beta_reference
      ),

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

    median_denominator_identified_fraction =
      median(
        denominator_ci_excludes_zero_fraction,
        na.rm = TRUE
      ),

    median_fieller_bounded_fraction =
      median(
        fieller_bounded_fraction,
        na.rm = TRUE
      ),

    mean_type1_Wald_05 =
      mean(
        baco_wald_type1_05_conditional,
        na.rm = TRUE
      ),

    median_type1_Wald_05 =
      median(
        baco_wald_type1_05_conditional,
        na.rm = TRUE
      ),

    max_type1_Wald_05 =
      max(
        baco_wald_type1_05_conditional,
        na.rm = TRUE
      ),

    mean_type1_Delta_05 =
      mean(
        delta_type1_05_conditional,
        na.rm = TRUE
      ),

    median_type1_Delta_05 =
      median(
        delta_type1_05_conditional,
        na.rm = TRUE
      ),

    max_type1_Delta_05 =
      max(
        delta_type1_05_conditional,
        na.rm = TRUE
      ),

    mean_type1_Fieller_05 =
      mean(
        fieller_type1_05_conditional,
        na.rm = TRUE
      ),

    median_type1_Fieller_05 =
      median(
        fieller_type1_05_conditional,
        na.rm = TRUE
      ),

    mean_Wald_coverage_95 =
      mean(
        baco_wald_coverage_95,
        na.rm = TRUE
      ),

    mean_Fieller_coverage_95 =
      mean(
        fieller_coverage_95,
        na.rm = TRUE
      ),

    .groups = "drop"
  )

write_csv(
  class_summary,
  file.path(
    OUTPUT_DIR,
    "02_confirmatory_summary_by_class.csv"
  )
)


# ------------------------------------------------------------------------------
# 8. Weak-denominator block: many critical events
# ------------------------------------------------------------------------------

weak_summary <- defined_results %>%
  filter(
    scenario_class ==
      "weak_denominator"
  ) %>%
  group_by(
    rr_reference,
    abs_population_beta_reference,
    expected_total_deaths_target
  ) %>%
  summarise(
    mean_Wald_type1 =
      mean(
        baco_wald_type1_05_conditional,
        na.rm = TRUE
      ),

    mean_Delta_type1 =
      mean(
        delta_type1_05_conditional,
        na.rm = TRUE
      ),

    mean_Fieller_type1 =
      mean(
        fieller_type1_05_conditional,
        na.rm = TRUE
      ),

    mean_Wald_coverage =
      mean(
        baco_wald_coverage_95,
        na.rm = TRUE
      ),

    mean_Fieller_coverage =
      mean(
        fieller_coverage_95,
        na.rm = TRUE
      ),

    mean_Fieller_bounded_fraction =
      mean(
        fieller_bounded_fraction,
        na.rm = TRUE
      ),

    .groups = "drop"
  )

write_csv(
  weak_summary,
  file.path(
    OUTPUT_DIR,
    "03_weak_denominator_many_events.csv"
  )
)


# ------------------------------------------------------------------------------
# 9. Sparse block: strong denominator
# ------------------------------------------------------------------------------

sparse_summary <- defined_results %>%
  filter(
    scenario_class ==
      "sparse"
  ) %>%
  group_by(
    expected_total_deaths_target,
    rr_reference
  ) %>%
  summarise(
    mean_ratio_estimable_fraction =
      mean(
        ratio_estimable_fraction,
        na.rm = TRUE
      ),

    mean_zero_death_any_arm_fraction =
      mean(
        zero_death_any_arm_fraction,
        na.rm = TRUE
      ),

    mean_Wald_type1 =
      mean(
        baco_wald_type1_05_conditional,
        na.rm = TRUE
      ),

    mean_Delta_type1 =
      mean(
        delta_type1_05_conditional,
        na.rm = TRUE
      ),

    mean_Fieller_type1 =
      mean(
        fieller_type1_05_conditional,
        na.rm = TRUE
      ),

    .groups = "drop"
  )

write_csv(
  sparse_summary,
  file.path(
    OUTPUT_DIR,
    "04_sparse_strong_denominator.csv"
  )
)


# ------------------------------------------------------------------------------
# 10. Combined sparse + weak block
# ------------------------------------------------------------------------------

combined_summary <- defined_results %>%
  filter(
    scenario_class ==
      "combined_sparse_weak"
  ) %>%
  group_by(
    expected_total_deaths_target,
    rr_reference
  ) %>%
  summarise(
    mean_ratio_estimable_fraction =
      mean(
        ratio_estimable_fraction,
        na.rm = TRUE
      ),

    mean_Wald_type1 =
      mean(
        baco_wald_type1_05_conditional,
        na.rm = TRUE
      ),

    mean_Delta_type1 =
      mean(
        delta_type1_05_conditional,
        na.rm = TRUE
      ),

    mean_Fieller_type1 =
      mean(
        fieller_type1_05_conditional,
        na.rm = TRUE
      ),

    mean_Fieller_bounded_fraction =
      mean(
        fieller_bounded_fraction,
        na.rm = TRUE
      ),

    .groups = "drop"
  )

write_csv(
  combined_summary,
  file.path(
    OUTPUT_DIR,
    "05_combined_sparse_weak.csv"
  )
)


# ------------------------------------------------------------------------------
# 11. Reference-null structural block
# ------------------------------------------------------------------------------

reference_null_summary <- results %>%
  filter(
    scenario_class ==
      "reference_null_BACO_undefined"
  ) %>%
  select(
    scenario_id,
    p_composite_control,
    expected_total_deaths_target,

    ratio_estimable_fraction,

    denominator_ci_excludes_zero_fraction,

    fieller_bounded_fraction,
    fieller_disjoint_unbounded_fraction,
    fieller_all_real_fraction,

    median_baco_hat,
    q025_baco_hat,
    q975_baco_hat,

    abs_BACO_gt_10_fraction,
    abs_BACO_gt_100_fraction
  )

write_csv(
  reference_null_summary,
  file.path(
    OUTPUT_DIR,
    "06_reference_null_structural_0_over_0.csv"
  )
)


# ------------------------------------------------------------------------------
# 12. Worst original-Wald scenarios
# ------------------------------------------------------------------------------

worst_Wald <- defined_results %>%
  arrange(
    desc(
      baco_wald_type1_05_conditional
    )
  ) %>%
  select(
    scenario_id,
    scenario_class,

    expected_total_deaths_target,
    rr_reference,
    abs_population_beta_reference,

    p_composite_control,

    ratio_estimable_fraction,
    denominator_ci_excludes_zero_fraction,

    fieller_bounded_fraction,

    baco_wald_type1_05_conditional,
    baco_wald_type1_05_mcse,

    delta_type1_05_conditional,
    fieller_type1_05_conditional,

    baco_wald_coverage_95,
    fieller_coverage_95,

    abs_BACO_gt_10_fraction,
    abs_BACO_gt_100_fraction
  ) %>%
  slice_head(
    n = 40
  )

write_csv(
  worst_Wald,
  file.path(
    OUTPUT_DIR,
    "07_worst_Wald_scenarios.csv"
  )
)


# ------------------------------------------------------------------------------
# 13. Internal Delta-Fieller validation
# ------------------------------------------------------------------------------

agreement_check <- defined_results %>%
  summarise(
    minimum_decision_agreement =
      min(
        delta_fieller_decision_agreement,
        na.rm = TRUE
      ),

    median_decision_agreement =
      median(
        delta_fieller_decision_agreement,
        na.rm = TRUE
      ),

    maximum_absolute_difference_in_rejection_rate =
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
    "08_Delta_Fieller_validation.csv"
  )
)


# ------------------------------------------------------------------------------
# 14. Plot data
# ------------------------------------------------------------------------------

type1_long <- defined_results %>%
  select(
    scenario_id,
    scenario_class,

    p_composite_control,
    expected_total_deaths_target,

    rr_reference,
    abs_population_beta_reference,

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
  type1_long,
  file.path(
    OUTPUT_DIR,
    "09_plot_ready_type1_long.csv"
  )
)


# ------------------------------------------------------------------------------
# 15. Figure 1
#     Confirmatory 2 x 2 conceptual design
# ------------------------------------------------------------------------------

figure1_data <- class_summary %>%
  select(
    scenario_class,

    median_ratio_estimable_fraction,
    median_denominator_identified_fraction,
    median_fieller_bounded_fraction,

    mean_type1_Wald_05,
    mean_type1_Delta_05,
    mean_type1_Fieller_05
  )

write_csv(
  figure1_data,
  file.path(
    OUTPUT_DIR,
    "10_figure1_class_summary_data.csv"
  )
)


# ------------------------------------------------------------------------------
# 16. Figure 2
#     Weak denominator isolated with many critical events
# ------------------------------------------------------------------------------

fig2_data <- type1_long %>%
  filter(
    scenario_class ==
      "weak_denominator",
    p_composite_control == 0.25
  )

p2 <- ggplot(
  fig2_data,
  aes(
    x =
      abs_population_beta_reference,

    y =
      empirical_type1,

    shape =
      method
  )
) +
  geom_hline(
    yintercept = ALPHA,
    linetype = 2
  ) +
  geom_point(
    size = 2.4
  ) +
  geom_line(
    aes(
      group =
        interaction(
          method,
          expected_total_deaths_target
        )
    ),
    linewidth = 0.45
  ) +
  facet_wrap(
    ~ expected_total_deaths_target,
    labeller = label_bquote(
      E[D] == .(
        expected_total_deaths_target
      )
    )
  ) +
  labs(
    x =
      expression(
        abs(log(RR[d]))
      ),

    y =
      "Empirical type-I error",

    shape = NULL,

    title =
      "Weak-denominator behavior persists despite many critical events",

    subtitle =
      paste0(
        "Population BACO = 1; ",
        NSIM,
        " replicates per scenario"
      )
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    legend.position =
      "bottom"
  )

ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "Figure_2_weak_denominator_many_events.png"
  ),
  plot = p2,
  width = 11,
  height = 7,
  dpi = 300
)


# ------------------------------------------------------------------------------
# 17. Figure 3
#     Sparse event information with strong denominator
# ------------------------------------------------------------------------------

fig3_data <- sparse_summary

p3 <- ggplot(
  fig3_data,
  aes(
    x =
      expected_total_deaths_target,

    y =
      mean_ratio_estimable_fraction,

    shape =
      factor(rr_reference)
  )
) +
  geom_point(
    size = 2.8
  ) +
  geom_line(
    aes(
      group =
        factor(rr_reference)
    ),
    linewidth = 0.5
  ) +
  scale_x_continuous(
    breaks =
      DEATHS_SPARSE
  ) +
  labs(
    x =
      "Expected total critical events",

    y =
      "Probability that original BACO is calculable",

    shape =
      expression(RR[d]),

    title =
      "Sparse critical events primarily threaten estimability"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    legend.position =
      "bottom"
  )

ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "Figure_3_sparse_estimability.png"
  ),
  plot = p3,
  width = 9,
  height = 6,
  dpi = 300
)


# ------------------------------------------------------------------------------
# 18. Figure 4
#     Wald coverage vs Fieller coverage in the weak-denominator block
# ------------------------------------------------------------------------------

coverage_long <- defined_results %>%
  filter(
    scenario_class ==
      "weak_denominator",
    p_composite_control == 0.25
  ) %>%
  select(
    expected_total_deaths_target,
    rr_reference,
    abs_population_beta_reference,

    baco_wald_coverage_95,
    fieller_coverage_95
  ) %>%
  pivot_longer(
    cols = c(
      baco_wald_coverage_95,
      fieller_coverage_95
    ),
    names_to =
      "method",
    values_to =
      "coverage"
  ) %>%
  mutate(
    method = recode(
      method,

      baco_wald_coverage_95 =
        "Original BACO Wald/delta",

      fieller_coverage_95 =
        "Fieller"
    )
  )

p4 <- ggplot(
  coverage_long,
  aes(
    x =
      abs_population_beta_reference,

    y =
      coverage,

    shape =
      method
  )
) +
  geom_hline(
    yintercept = 0.95,
    linetype = 2
  ) +
  geom_point(
    size = 2.5
  ) +
  geom_line(
    aes(
      group =
        interaction(
          method,
          expected_total_deaths_target
        )
    ),
    linewidth = 0.45
  ) +
  facet_wrap(
    ~ expected_total_deaths_target,
    labeller = label_bquote(
      E[D] == .(
        expected_total_deaths_target
      )
    )
  ) +
  labs(
    x =
      expression(
        abs(log(RR[d]))
      ),

    y =
      "Empirical 95% coverage",

    shape = NULL,

    title =
      "Wald coverage deteriorates as the BACO denominator approaches zero"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    legend.position =
      "bottom"
  )

ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "Figure_4_coverage_weak_denominator.png"
  ),
  plot = p4,
  width = 11,
  height = 7,
  dpi = 300
)


# ------------------------------------------------------------------------------
# 19. Monte Carlo precision table
# ------------------------------------------------------------------------------

nominal_precision <- tibble(
  nsim = NSIM,

  nominal_probability =
    c(
      0.05,
      0.005,
      0.95
    ),

  Monte_Carlo_SE =
    sqrt(
      nominal_probability *
      (1 - nominal_probability) /
      NSIM
    ),

  approximate_95_MC_half_width =
    1.96 * Monte_Carlo_SE
)

write_csv(
  nominal_precision,
  file.path(
    OUTPUT_DIR,
    "11_Monte_Carlo_precision.csv"
  )
)


# ------------------------------------------------------------------------------
# 20. Console summary
# ------------------------------------------------------------------------------

cat(
  "\n============================================================\n"
)

cat(
  "BACO confirmatory simulation complete\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Scenarios by class:\n"
)

print(
  scenario_grid %>%
    count(
      scenario_class
    )
)

cat(
  "\nClass-level results:\n"
)

print(
  class_summary
)

cat(
  "\nDelta-Fieller internal validation:\n"
)

print(
  agreement_check
)

cat(
  "\nMonte Carlo precision:\n"
)

print(
  nominal_precision
)

cat(
  "\nWorst Wald scenarios:\n"
)

print(
  worst_Wald %>%
    slice_head(
      n = 10
    )
)

cat(
  "\nOutputs saved in:\n",
  normalizePath(
    OUTPUT_DIR,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)


# ==============================================================================
# INTERPRETATION FRAMEWORK
# ==============================================================================
#
# REGULAR
# -------
# Many critical events + reference effect clearly away from RR = 1.
#
# Expected:
#   - high estimability;
#   - bounded Fieller sets common;
#   - all methods approach nominal calibration.
#
#
# WEAK DENOMINATOR
# ----------------
# Many critical events + RR_d close to 1.
#
# Key question:
#   Does the original BACO Wald test remain anti-conservative even when
#   sparse-event problems have been removed?
#
# If yes, this isolates a ratio/parameterization problem rather than a
# sparse-data problem.
#
#
# SPARSE
# ------
# Few critical events + RR_d clearly away from 1.
#
# Key question:
#   Does sparsity primarily reduce the probability that the original BACO
#   can be estimated at all?
#
# Conditional calibration must always be reported together with estimability.
#
#
# COMBINED
# --------
# Few critical events + RR_d close to 1.
#
# Expected:
#   - frequent zero-event cells;
#   - low BACO estimability;
#   - severe denominator weakness;
#   - worst overall behavior.
#
#
# REFERENCE NULL
# --------------
# RR_c = RR_d = 1:
#
#       BACO = log(1) / log(1) = 0/0.
#
# There is no defined population BACO parameter here.
# Sample BACO values in this block are diagnostic only.
#
#
# NEXT STEP AFTER CONFIRMATION
# ----------------------------
# If the confirmatory results reproduce the pilot:
#
#   1) freeze the evidence that Wald ratio inference is unreliable;
#   2) retain Fieller for ratio confidence sets;
#   3) retain Delta = beta_c - beta_d for testing H0: BACO = 1;
#   4) move to zero-event estimation:
#        - profile likelihood,
#        - bias reduction / penalization,
#        - Bayesian multinomial inference,
#        - continuity correction only as sensitivity analysis.
#
# ==============================================================================
