# ==============================================================================
# BACO reproducibility analysis — supplementary sensitivity
# Power and size-adjusted power under BACO-parameterized alternatives
#
# PURPOSE
# -------
# The preceding simulations established:
#   1) severe type-I error inflation for the original BACO Wald ratio test;
#   2) correct asymptotic behavior of Delta/Fieller with adequate information;
#   3) valid but conservative Fisher exact inference under sparse/zero events.
#
# The next question is POWER:
#
#       How much power is lost by exact inference,
#       especially when the critical endpoint is sparse?
#
# IMPORTANT: raw power alone is not a fair comparison when methods have
# different type-I error. Therefore this script reports BOTH:
#
#   A) nominal power at p < 0.05;
#   B) SIZE-ADJUSTED POWER, using method- and scenario-specific empirical
#      critical p-value cutoffs estimated under H0: BACO = 1.
#
# This prevents an anti-conservative method from appearing "more powerful"
# merely because it rejects too often under the null.
#
#
# POPULATION PARAMETERIZATION
# ---------------------------
# Let:
#
#       theta = BACO = log(RR_c) / log(RR_d)
#
# For a chosen RR_d and theta:
#
#       RR_c = RR_d^theta
#
# and:
#
#       Delta = log(RR_c) - log(RR_d)
#             = (theta - 1) * log(RR_d).
#
# The conditional critical-event risk ratio among composite-event participants:
#
#       QRR = (D_t/C_t) / (D_c/C_c)
#           = RR_d / RR_c
#           = exp(-Delta).
#
# This is important because the exact test is testing QRR = 1.
#
# A large numerical departure of BACO from 1 near RR_d = 1 can correspond to
# a very small Delta and therefore a genuinely difficult-to-detect difference.
# The script therefore reports power against BOTH theta and Delta.
#
#
# METHODS
# -------
# 1) Fisher exact conditional test
# 2) Fisher probability-ordered mid-P
# 3) Conditional likelihood-ratio test with chi-square(1) reference
# 4) Direct Delta Wald test
# 5) Original BACO Wald/delta ratio test
#
#
# RECOMMENDED PILOT SETTINGS
# --------------------------
# Null calibration:
#       NSIM_NULL <- 10000
#
# Alternative power:
#       NSIM_ALT  <- 5000
#
# After checking the results, increase NSIM_ALT to 10000 or 20000 for the
# manuscript-final run.
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
# 1. Settings
# ------------------------------------------------------------------------------

ROOT_DIR <- getwd()
OUTPUT_DIR <- file.path(ROOT_DIR, "results", "supplementary_baco_power")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

SEED <- 20260902

NSIM_NULL <- 10000
NSIM_ALT <- 5000

ALPHA <- 0.05
ALPHA_STRICT <- 0.005

N_PER_ARM <- 10000

# Control-arm composite risk.
P_COMPOSITE_CONTROL <- c(
  0.10,
  0.25
)

# Expected total critical events across both groups.
EXPECTED_TOTAL_DEATHS <- c(
  1,
  2,
  3,
  5,
  10,
  20,
  50
)

# Strong and weak reference effects.
RR_D_REFERENCE <- c(
  0.50,
  0.95,
  1.05,
  2.00
)

# BACO alternatives.
#
# theta < 0  : reversal of association direction between composite and critical
# theta = 0  : no composite effect while critical effect is present
# 0 < theta < 1 : attenuation
# theta > 1  : amplification
THETA_ALT <- c(
  -1.00,
  0.00,
  0.25,
  0.50,
  0.75,
  1.25,
  1.50,
  2.00,
  3.00
)

NUMERIC_TOL <- sqrt(.Machine$double.eps)


# ------------------------------------------------------------------------------
# 2. Base nuisance-parameter scenarios
# ------------------------------------------------------------------------------

# With equal group sizes:
#
#   E[D_total] = N*p_d0 + N*p_d1
#              = N*p_d0*(1 + RR_d)
#
# Therefore:
#
#   p_d0 = E[D_total] / {N*(1 + RR_d)}

base_grid <- tidyr::crossing(
  n_per_arm = N_PER_ARM,
  p_composite_control = P_COMPOSITE_CONTROL,
  expected_total_deaths_target = EXPECTED_TOTAL_DEATHS,
  rr_d = RR_D_REFERENCE
) %>%
  mutate(
    p_death_control = expected_total_deaths_target /
      (n_per_arm *
        (1 + rr_d)),

    p_death_treatment = rr_d *
      p_death_control,

    denominator_class = ifelse(
      rr_d %in%
        c(
          0.95,
          1.05
        ),
      "weak",
      "strong"
    ),

    base_scenario_id = row_number(),

    null_rr_c = rr_d,

    null_p_composite_treatment = rr_d *
      p_composite_control,

    valid_null = p_death_control >= 0 &
      p_death_control <= p_composite_control &
      p_death_treatment >= 0 &
      p_death_treatment <= null_p_composite_treatment &
      null_p_composite_treatment > 0 &
      null_p_composite_treatment < 1
  ) %>%
  filter(
    valid_null
  )

write_csv(
  base_grid,
  file.path(
    OUTPUT_DIR,
    "00_base_nuisance_scenarios.csv"
  )
)


# ------------------------------------------------------------------------------
# 3. Alternative scenarios
# ------------------------------------------------------------------------------

alt_grid <- tidyr::crossing(
  base_scenario_id = base_grid$base_scenario_id,

  theta = THETA_ALT
) %>%
  left_join(
    base_grid,
    by = "base_scenario_id"
  ) %>%
  mutate(
    rr_c = rr_d^theta,

    p_composite_treatment = rr_c *
      p_composite_control,

    population_Delta = log(
      rr_c
    ) -
      log(
        rr_d
      ),

    abs_population_Delta = abs(
      population_Delta
    ),

    population_conditional_critical_RR = rr_d /
      rr_c,

    alternative_class = case_when(
      theta < 0 ~
        "association_reversal",

      theta == 0 ~
        "null_composite_effect",

      theta > 0 &
        theta < 1 ~
        "attenuation",

      theta > 1 ~
        "amplification",

      TRUE ~
        "other"
    ),

    valid_alternative = p_composite_treatment > 0 &
      p_composite_treatment < 1 &
      p_death_treatment >= 0 &
      p_death_treatment <= p_composite_treatment
  ) %>%
  filter(
    valid_alternative
  ) %>%
  mutate(
    alternative_scenario_id = row_number()
  ) %>%
  select(
    alternative_scenario_id,
    everything()
  )

write_csv(
  alt_grid,
  file.path(
    OUTPUT_DIR,
    "01_power_alternative_scenario_grid.csv"
  )
)

cat(
  "\nBase scenarios: ",
  nrow(
    base_grid
  ),
  "\n",
  sep = ""
)

cat(
  "Valid alternative scenarios: ",
  nrow(
    alt_grid
  ),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------------------------
# 4. Fisher probability-ordered and mid-P tests
# ------------------------------------------------------------------------------

fisher_probability_ordered <- function(
  d_t,
  c_t,
  c_c,
  d_total
) {
  if (
    c_t <= 0 ||
      c_c <= 0
  ) {
    return(
      c(
        exact = NA_real_,
        midp = NA_real_
      )
    )
  }

  if (d_total == 0) {
    return(
      c(
        exact = 1,
        midp = 1
      )
    )
  }

  lo <-
    max(
      0,
      d_total -
        c_c
    )

  hi <-
    min(
      c_t,
      d_total
    )

  support <-
    lo:hi

  probs <- dhyper(
    support,
    m = c_t,
    n = c_c,
    k = d_total
  )

  p_obs <- dhyper(
    d_t,
    m = c_t,
    n = c_c,
    k = d_total
  )

  exact_p <-
    sum(
      probs[
        probs <=
          p_obs *
            (1 +
              1e-12)
      ]
    )

  exact_p <-
    min(
      1,
      exact_p
    )

  mid_p <-
    max(
      0,
      exact_p -
        0.5 *
          p_obs
    )

  c(
    exact = exact_p,

    midp = mid_p
  )
}


exact_p_vector <- function(
  d_t,
  c_t,
  d_c,
  c_c
) {
  n <-
    length(
      d_t
    )

  exact <-
    rep(
      NA_real_,
      n
    )

  midp <-
    rep(
      NA_real_,
      n
    )

  for (i in seq_len(
    n
  )) {
    ans <-
      fisher_probability_ordered(
        d_t = d_t[i],

        c_t = c_t[i],

        c_c = c_c[i],

        d_total = d_t[i] +
          d_c[i]
      )

    exact[i] <-
      ans[
        "exact"
      ]

    midp[i] <-
      ans[
        "midp"
      ]
  }

  tibble(
    Fisher = exact,

    MidP = midp
  )
}


# ------------------------------------------------------------------------------
# 5. Conditional LR test
# ------------------------------------------------------------------------------

binom_loglik <- function(
  d,
  c,
  q
) {
  if (c == 0) {
    return(0)
  }

  if (q <= 0) {
    return(
      ifelse(
        d == 0,
        0,
        -Inf
      )
    )
  }

  if (q >= 1) {
    return(
      ifelse(
        d == c,
        0,
        -Inf
      )
    )
  }

  d *
    log(
      q
    ) +
    (c -
      d) *
      log1p(
        -q
      )
}


conditional_lr_vector <- function(
  d_t,
  c_t,
  d_c,
  c_c
) {
  n <-
    length(
      d_t
    )

  p <-
    rep(
      NA_real_,
      n
    )

  for (i in seq_len(
    n
  )) {
    if (
      c_t[i] <= 0 ||
        c_c[i] <= 0
    ) {
      next
    }

    d_total <-
      d_t[i] +
      d_c[i]

    c_total <-
      c_t[i] +
      c_c[i]

    if (d_total == 0) {
      p[i] <- 1
      next
    }

    q_t <-
      d_t[i] /
      c_t[i]

    q_c <-
      d_c[i] /
      c_c[i]

    q_0 <-
      d_total /
      c_total

    ll_u <-
      binom_loglik(
        d_t[i],
        c_t[i],
        q_t
      ) +
      binom_loglik(
        d_c[i],
        c_c[i],
        q_c
      )

    ll_0 <-
      binom_loglik(
        d_t[i],
        c_t[i],
        q_0
      ) +
      binom_loglik(
        d_c[i],
        c_c[i],
        q_0
      )

    LR <-
      max(
        0,
        2 *
          (ll_u -
            ll_0)
      )

    p[i] <-
      pchisq(
        LR,
        df = 1,
        lower.tail = FALSE
      )
  }

  p
}


# ------------------------------------------------------------------------------
# 6. Delta and original BACO Wald tests
# ------------------------------------------------------------------------------

asymptotic_pvalues <- function(
  n_t,
  n_c,
  comp_t,
  comp_c,
  death_t,
  death_c
) {
  nsim <-
    length(
      comp_t
    )

  p_delta <-
    rep(
      NA_real_,
      nsim
    )

  p_wald <-
    rep(
      NA_real_,
      nsim
    )

  valid <-
    comp_t > 0 &
    comp_c > 0 &
    death_t > 0 &
    death_c > 0

  idx <-
    which(
      valid
    )

  if (
    length(
      idx
    ) ==
      0
  ) {
    return(
      tibble(
        Delta = p_delta,

        BACO_Wald = p_wald
      )
    )
  }

  rr_c <-
    (comp_t[idx] /
      n_t) /
    (comp_c[idx] /
      n_c)

  rr_d <-
    (death_t[idx] /
      n_t) /
    (death_c[idx] /
      n_c)

  beta_c <-
    log(
      rr_c
    )

  beta_d <-
    log(
      rr_d
    )

  var_c <-
    1 /
    comp_t[idx] +
    1 /
      comp_c[idx] -
    1 /
      n_t -
    1 /
      n_c

  var_d <-
    1 /
    death_t[idx] +
    1 /
      death_c[idx] -
    1 /
      n_t -
    1 /
      n_c

  cov_cd <-
    var_c

  # Delta inference.
  delta <-
    beta_c -
    beta_d

  var_delta <-
    var_c +
    var_d -
    2 *
      cov_cd

  good_delta <-
    is.finite(
      var_delta
    ) &
    var_delta > 0

  if (
    any(
      good_delta
    )
  ) {
    z_delta <-
      delta[
        good_delta
      ] /
      sqrt(
        var_delta[
          good_delta
        ]
      )

    p_delta[
      idx[
        good_delta
      ]
    ] <-
      2 *
      pnorm(
        -abs(
          z_delta
        )
      )
  }

  # Original BACO Wald ratio.
  good_ratio <-
    is.finite(
      beta_d
    ) &
    abs(
      beta_d
    ) >
      NUMERIC_TOL

  if (
    any(
      good_ratio
    )
  ) {
    ridx <-
      idx[
        good_ratio
      ]

    bc <-
      beta_c[
        good_ratio
      ]

    bd <-
      beta_d[
        good_ratio
      ]

    vc <-
      var_c[
        good_ratio
      ]

    vd <-
      var_d[
        good_ratio
      ]

    cv <-
      cov_cd[
        good_ratio
      ]

    theta_hat <-
      bc /
      bd

    g1 <-
      1 /
      bd

    g2 <-
      -bc /
      bd^2

    var_theta <-
      g1^2 *
      vc +
      g2^2 *
        vd +
      2 *
        g1 *
        g2 *
        cv

    good_var <-
      is.finite(
        var_theta
      ) &
      var_theta > 0

    if (
      any(
        good_var
      )
    ) {
      final_idx <-
        ridx[
          good_var
        ]

      z_theta <-
        (theta_hat[
          good_var
        ] -
          1) /
        sqrt(
          var_theta[
            good_var
          ]
        )

      p_wald[
        final_idx
      ] <-
        2 *
        pnorm(
          -abs(
            z_theta
          )
        )
    }
  }

  tibble(
    Delta = p_delta,

    BACO_Wald = p_wald
  )
}


# ------------------------------------------------------------------------------
# 7. Simulate one p-value dataset
# ------------------------------------------------------------------------------

simulate_pvalues <- function(
  n_per_arm,
  p_composite_control,
  p_composite_treatment,
  p_death_control,
  p_death_treatment,
  nsim,
  seed
) {
  prob_c <- c(
    none = 1 -
      p_composite_control,

    nonfatal = p_composite_control -
      p_death_control,

    critical = p_death_control
  )

  prob_t <- c(
    none = 1 -
      p_composite_treatment,

    nonfatal = p_composite_treatment -
      p_death_treatment,

    critical = p_death_treatment
  )

  if (
    any(
      prob_c < 0
    ) ||
      any(
        prob_t < 0
      ) ||
      abs(
        sum(
          prob_c
        ) -
          1
      ) >
        1e-10 ||
      abs(
        sum(
          prob_t
        ) -
          1
      ) >
        1e-10
  ) {
    stop(
      "Invalid multinomial probabilities."
    )
  }

  set.seed(
    seed
  )

  draw_c <-
    rmultinom(
      n = nsim,

      size = n_per_arm,

      prob = prob_c
    )

  draw_t <-
    rmultinom(
      n = nsim,

      size = n_per_arm,

      prob = prob_t
    )

  comp_c <-
    draw_c[
      "nonfatal",
    ] +
    draw_c[
      "critical",
    ]

  comp_t <-
    draw_t[
      "nonfatal",
    ] +
    draw_t[
      "critical",
    ]

  death_c <-
    draw_c[
      "critical",
    ]

  death_t <-
    draw_t[
      "critical",
    ]

  exact <-
    exact_p_vector(
      d_t = death_t,

      c_t = comp_t,

      d_c = death_c,

      c_c = comp_c
    )

  p_lr <-
    conditional_lr_vector(
      d_t = death_t,

      c_t = comp_t,

      d_c = death_c,

      c_c = comp_c
    )

  asym <-
    asymptotic_pvalues(
      n_t = n_per_arm,

      n_c = n_per_arm,

      comp_t = comp_t,

      comp_c = comp_c,

      death_t = death_t,

      death_c = death_c
    )

  tibble(
    Fisher = exact$Fisher,

    MidP = exact$MidP,

    Conditional_LR = p_lr,

    Delta = asym$Delta,

    BACO_Wald = asym$BACO_Wald,

    zero_death_any_arm = death_t == 0 |
      death_c == 0,

    zero_death_both_arms = death_t == 0 &
      death_c == 0
  )
}


# ------------------------------------------------------------------------------
# 8. Empirical size-adjusted critical cutoff
# ------------------------------------------------------------------------------

# Non-estimable p-values are treated as p = 1:
# the complete procedure does not reject when it cannot compute the test.
#
# The cutoff is the largest observed p-value threshold whose empirical null
# rejection rate does not exceed alpha.

calibrated_cutoff <- function(
  p_values,
  alpha = ALPHA
) {
  p <-
    ifelse(
      is.na(
        p_values
      ),
      1,
      p_values
    )

  candidates <-
    sort(
      unique(
        p
      )
    )

  cdf <-
    vapply(
      candidates,
      function(x) {
        mean(
          p <= x
        )
      },
      numeric(1)
    )

  ok <-
    which(
      cdf <=
        alpha +
          1e-12
    )

  if (
    length(
      ok
    ) ==
      0
  ) {
    cutoff <-
      0

    achieved_size <-
      0
  } else {
    j <-
      max(
        ok
      )

    cutoff <-
      candidates[j]

    achieved_size <-
      cdf[j]
  }

  c(
    cutoff = cutoff,

    achieved_size = achieved_size
  )
}


# ------------------------------------------------------------------------------
# 9. Null calibration for every nuisance scenario
# ------------------------------------------------------------------------------

METHODS <- c(
  "Fisher",
  "MidP",
  "Conditional_LR",
  "Delta",
  "BACO_Wald"
)

null_calibration_list <-
  vector(
    "list",
    nrow(
      base_grid
    )
  )

cat(
  "\nCalibrating null distributions...\n"
)

for (i in seq_len(
  nrow(base_grid)
)) {
  s <- base_grid[i, ]

  if (
    i == 1 ||
      i %% 10 == 0 ||
      i == nrow(base_grid)
  ) {
    cat(
      "Null scenario ",
      i,
      " / ",
      nrow(base_grid),
      "\n",
      sep = ""
    )
  }

  null_p <-
    simulate_pvalues(
      n_per_arm = s$n_per_arm,

      p_composite_control = s$p_composite_control,

      p_composite_treatment = s$null_p_composite_treatment,

      p_death_control = s$p_death_control,

      p_death_treatment = s$p_death_treatment,

      nsim = NSIM_NULL,

      seed = SEED +
        1009L *
          s$base_scenario_id
    )

  rows <-
    map_dfr(
      METHODS,
      function(m) {
        p_raw <- null_p[[m]]

        p_proc <-
          ifelse(
            is.na(p_raw),
            1,
            p_raw
          )

        cal <-
          calibrated_cutoff(
            p_raw,
            alpha = ALPHA
          )

        tibble(
          base_scenario_id = s$base_scenario_id,

          method = m,

          null_estimable_fraction = mean(
            !is.na(p_raw)
          ),

          nominal_type1_05 = mean(
            p_proc < ALPHA
          ),

          nominal_type1_005 = mean(
            p_proc < ALPHA_STRICT
          ),

          calibrated_cutoff_05 = cal["cutoff"],

          calibrated_null_size_05 = cal["achieved_size"]
        )
      }
    )

  null_calibration_list[[i]] <-
    rows
}

null_calibration <-
  bind_rows(
    null_calibration_list
  ) %>%
  left_join(
    base_grid %>%
      select(
        base_scenario_id,
        denominator_class,
        rr_d,
        p_composite_control,
        expected_total_deaths_target
      ),
    by = "base_scenario_id"
  )

write_csv(
  null_calibration,
  file.path(
    OUTPUT_DIR,
    "02_null_method_calibration.csv"
  )
)


# ------------------------------------------------------------------------------
# 10. Alternative power simulation
# ------------------------------------------------------------------------------

power_list <-
  vector(
    "list",
    nrow(
      alt_grid
    )
  )

cat(
  "\nRunning alternative power simulations...\n"
)

for (i in seq_len(
  nrow(alt_grid)
)) {
  s <- alt_grid[i, ]

  if (
    i == 1 ||
      i %% 25 == 0 ||
      i == nrow(alt_grid)
  ) {
    cat(
      "Alternative scenario ",
      i,
      " / ",
      nrow(alt_grid),
      "\n",
      sep = ""
    )
  }

  alt_p <-
    simulate_pvalues(
      n_per_arm = s$n_per_arm,

      p_composite_control = s$p_composite_control,

      p_composite_treatment = s$p_composite_treatment,

      p_death_control = s$p_death_control,

      p_death_treatment = s$p_death_treatment,

      nsim = NSIM_ALT,

      seed = SEED +
        1000003L +
        2017L *
          s$alternative_scenario_id
    )

  cal_s <-
    null_calibration %>%
    filter(
      base_scenario_id == s$base_scenario_id
    )

  rows <-
    map_dfr(
      METHODS,
      function(m) {
        p_raw <- alt_p[[m]]

        p_proc <-
          ifelse(
            is.na(p_raw),
            1,
            p_raw
          )

        cutoff <-
          cal_s %>%
          filter(
            method == m
          ) %>%
          pull(
            calibrated_cutoff_05
          )

        if (length(cutoff) != 1) {
          stop(
            "Calibration cutoff not uniquely identified."
          )
        }

        tibble(
          alternative_scenario_id = s$alternative_scenario_id,

          base_scenario_id = s$base_scenario_id,

          method = m,

          estimable_fraction = mean(
            !is.na(p_raw)
          ),

          nominal_power_05 = mean(
            p_proc < ALPHA
          ),

          nominal_power_005 = mean(
            p_proc < ALPHA_STRICT
          ),

          size_adjusted_power_05 = mean(
            p_proc <= cutoff
          ),

          calibrated_cutoff_05 = cutoff
        )
      }
    )

  power_list[[i]] <-
    rows
}

power_results <-
  bind_rows(
    power_list
  ) %>%
  left_join(
    alt_grid %>%
      select(
        alternative_scenario_id,
        base_scenario_id,

        denominator_class,
        rr_d,
        rr_c,
        theta,

        population_Delta,
        abs_population_Delta,
        population_conditional_critical_RR,

        alternative_class,

        p_composite_control,
        expected_total_deaths_target
      ),
    by = c(
      "alternative_scenario_id",
      "base_scenario_id"
    )
  )

write_csv(
  power_results,
  file.path(
    OUTPUT_DIR,
    "03_power_all_methods_all_scenarios.csv"
  )
)


# ------------------------------------------------------------------------------
# 11. Summary by method
# ------------------------------------------------------------------------------

summary_by_method <-
  power_results %>%
  group_by(
    method,
    denominator_class
  ) %>%
  summarise(
    n_scenarios = n(),

    median_estimable_fraction = median(
      estimable_fraction
    ),

    median_nominal_power = median(
      nominal_power_05
    ),

    median_size_adjusted_power = median(
      size_adjusted_power_05
    ),

    mean_nominal_power = mean(
      nominal_power_05
    ),

    mean_size_adjusted_power = mean(
      size_adjusted_power_05
    ),

    .groups = "drop"
  )

write_csv(
  summary_by_method,
  file.path(
    OUTPUT_DIR,
    "04_power_summary_by_method.csv"
  )
)


# ------------------------------------------------------------------------------
# 12. Summary by expected critical events and method
# ------------------------------------------------------------------------------

summary_by_events <-
  power_results %>%
  group_by(
    method,
    denominator_class,
    expected_total_deaths_target
  ) %>%
  summarise(
    n_scenarios = n(),

    mean_estimable_fraction = mean(
      estimable_fraction
    ),

    mean_nominal_power = mean(
      nominal_power_05
    ),

    mean_size_adjusted_power = mean(
      size_adjusted_power_05
    ),

    .groups = "drop"
  )

write_csv(
  summary_by_events,
  file.path(
    OUTPUT_DIR,
    "05_power_summary_by_expected_events.csv"
  )
)


# ------------------------------------------------------------------------------
# 13. Summary by BACO alternative
# ------------------------------------------------------------------------------

summary_by_theta <-
  power_results %>%
  group_by(
    method,
    denominator_class,
    theta,
    alternative_class
  ) %>%
  summarise(
    n_scenarios = n(),

    median_abs_Delta = median(
      abs_population_Delta
    ),

    mean_nominal_power = mean(
      nominal_power_05
    ),

    mean_size_adjusted_power = mean(
      size_adjusted_power_05
    ),

    .groups = "drop"
  )

write_csv(
  summary_by_theta,
  file.path(
    OUTPUT_DIR,
    "06_power_summary_by_BACO_alternative.csv"
  )
)


# ------------------------------------------------------------------------------
# 14. Figure 1:
#     Size-adjusted power vs actual log-effect discrepancy
# ------------------------------------------------------------------------------

fig1_data <-
  power_results %>%
  filter(
    p_composite_control == 0.25,
    expected_total_deaths_target %in%
      c(
        3,
        5,
        10,
        20
      )
  ) %>%
  mutate(
    deaths_facet = factor(
      expected_total_deaths_target,
      levels = c(3, 5, 10, 20),
      labels = c(
        "E[D] = 3",
        "E[D] = 5",
        "E[D] = 10",
        "E[D] = 20"
      )
    )
  )

p1 <- ggplot(
  fig1_data,
  aes(
    x = abs_population_Delta,
    y = size_adjusted_power_05,
    shape = method
  )
) +
  geom_line(
    aes(
      group = interaction(
        method,
        rr_d
      )
    ),
    linewidth = 0.35
  ) +
  geom_point(
    size = 2.3
  ) +
  facet_grid(
    denominator_class ~ deaths_facet
  ) +
  labs(
    x = expression(
      abs(Delta) ==
        abs(
          log(RR[c]) -
            log(RR[d])
        )
    ),
    y = "Size-adjusted power",
    shape = NULL,
    title = "Power is driven by the actual composite-critical effect discrepancy",
    subtitle = "Method-specific critical values calibrated to empirical 5% null size"
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    legend.position = "bottom"
  )

ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "Figure_1_size_adjusted_power_vs_Delta.png"
  ),
  plot = p1,
  width = 13,
  height = 7,
  dpi = 300
)

# ------------------------------------------------------------------------------
# 15. Figure 2:
#     Power vs expected event count for selected BACO alternatives
# ------------------------------------------------------------------------------

fig2_data <-
  power_results %>%
  filter(
    p_composite_control == 0.25,

    theta %in%
      c(
        -1,
        0.5,
        1.5,
        2
      )
  )

p2 <- ggplot(
  fig2_data,
  aes(
    x = expected_total_deaths_target,

    y = size_adjusted_power_05,

    shape = method
  )
) +
  geom_point(
    size = 2.4
  ) +
  geom_line(
    aes(
      group = interaction(
        method,
        rr_d
      )
    ),
    linewidth = 0.4
  ) +
  facet_grid(
    theta ~
      denominator_class,

    labeller = label_bquote(
      BACO ==
        .(
          theta
        )
    )
  ) +
  scale_x_continuous(
    breaks = EXPECTED_TOTAL_DEATHS
  ) +
  labs(
    x = "Expected total critical events",

    y = "Size-adjusted power",

    shape = NULL,

    title = "Power increases with critical-event information but depends strongly on reference-effect strength"
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    legend.position = "bottom"
  )

ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "Figure_2_power_vs_expected_events.png"
  ),

  plot = p2,

  width = 11,
  height = 9,
  dpi = 300
)


# ------------------------------------------------------------------------------
# 16. Figure 3:
#     Nominal vs size-adjusted power
# ------------------------------------------------------------------------------

fig3_data <-
  power_results %>%
  filter(
    p_composite_control == 0.25
  )

p3 <- ggplot(
  fig3_data,
  aes(
    x = nominal_power_05,

    y = size_adjusted_power_05,

    shape = method
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2
  ) +
  geom_point(
    alpha = 0.65,
    size = 2.1
  ) +
  facet_wrap(
    ~denominator_class
  ) +
  labs(
    x = "Nominal power at p < 0.05",

    y = "Size-adjusted power",

    shape = NULL,

    title = "Apparent power can be inflated when the test itself is anti-conservative"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    legend.position = "bottom"
  )

ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "Figure_3_nominal_vs_size_adjusted_power.png"
  ),

  plot = p3,

  width = 10,
  height = 6,
  dpi = 300
)


# ------------------------------------------------------------------------------
# 17. Figure 4:
#     Why BACO magnitude alone is misleading near RR_d = 1
# ------------------------------------------------------------------------------

concept_data <-
  alt_grid %>%
  distinct(
    denominator_class,
    rr_d,
    theta,
    population_Delta,
    abs_population_Delta,
    population_conditional_critical_RR
  )

p4 <- ggplot(
  concept_data,
  aes(
    x = theta,

    y = abs_population_Delta,

    shape = factor(
      rr_d
    )
  )
) +
  geom_point(
    size = 2.7
  ) +
  geom_line(
    aes(
      group = rr_d
    ),
    linewidth = 0.45
  ) +
  labs(
    x = "Population BACO",

    y = expression(
      abs(
        Delta
      )
    ),

    shape = expression(
      RR[d]
    ),

    title = "The same BACO departure can imply very different effect discrepancies",

    subtitle = "Near RRd = 1, large BACO values may correspond to small absolute log-effect differences"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    legend.position = "bottom"
  )

ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "Figure_4_BACO_vs_actual_effect_discrepancy.png"
  ),

  plot = p4,

  width = 9,
  height = 6,
  dpi = 300
)


# ------------------------------------------------------------------------------
# 18. Monte Carlo precision
# ------------------------------------------------------------------------------

mc_precision <-
  tibble(
    component = c(
      "null calibration",
      "alternative power"
    ),

    n_sim = c(
      NSIM_NULL,
      NSIM_ALT
    ),

    probability = 0.05,

    approximate_MCSE_at_5pct = sqrt(
      0.05 *
        0.95 /
        n_sim
    )
  )

write_csv(
  mc_precision,
  file.path(
    OUTPUT_DIR,
    "07_Monte_Carlo_precision.csv"
  )
)


# ------------------------------------------------------------------------------
# 19. Console report
# ------------------------------------------------------------------------------

cat(
  "\n============================================================\n"
)

cat(
  "BACO power simulation complete\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Null calibration by method:\n"
)

print(
  null_calibration %>%
    group_by(
      method,
      denominator_class
    ) %>%
    summarise(
      mean_nominal_type1 = mean(
        nominal_type1_05
      ),

      max_nominal_type1 = max(
        nominal_type1_05
      ),

      mean_calibrated_size = mean(
        calibrated_null_size_05
      ),

      .groups = "drop"
    )
)

cat(
  "\nPower summary by method:\n"
)

print(
  summary_by_method
)

cat(
  "\nMonte Carlo precision:\n"
)

print(
  mc_precision
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
# 1) PRIMARY power comparison:
#       use SIZE-ADJUSTED POWER.
#
#    Raw nominal power is misleading for methods whose null rejection rate
#    substantially exceeds 5%.
#
#
# 2) Fisher exact:
#       expected to have the strongest size control but lower power with very
#       sparse data due to discreteness.
#
#
# 3) Mid-P:
#       expected to recover part of the power loss. Its acceptability depends on
#       whether its null calibration remains sufficiently close to 5%.
#
#
# 4) Conditional LR:
#       raw power may appear high in sparse settings because its chi-square
#       approximation can be anti-conservative. Size-adjusted power allows a
#       fairer comparison.
#
#
# 5) Delta:
#       should become efficient as event information increases, but its complete
#       procedure cannot reject when a zero critical cell makes it non-estimable.
#
#
# 6) Original BACO Wald:
#       any apparent power advantage near a weak denominator must be interpreted
#       against its severe type-I error inflation. Size-adjusted power is the
#       relevant comparison.
#
#
# 7) PARAMETER INTERPRETATION:
#
#       Delta = (BACO - 1) * log(RR_d)
#
#    Therefore a fixed BACO departure has very different substantive magnitude
#    depending on RR_d.
#
#    Example:
#       BACO = 2, RR_d = 0.50 -> |Delta| = |log(0.50)| = 0.693
#       BACO = 2, RR_d = 0.95 -> |Delta| = |log(0.95)| = 0.051
#
#    Power should therefore be interpreted against Delta / conditional QRR,
#    not against BACO magnitude alone.
#
#
# NEXT STEP AFTER THIS SIMULATION
# -------------------------------
# If size-adjusted power confirms Fisher as the best calibrated sparse method:
#
#   - freeze the testing framework;
#   - perform a smaller sensitivity analysis for unequal allocation;
#   - optionally study mid-P/Boschloo;
#   - then begin manuscript Methods/Results and package API design.
#
# ==============================================================================
