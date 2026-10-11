# ==============================================================================
# BACO reproducibility analysis
# Final Delta/QRR-parameterized simulation used in the JCE manuscript
#
# WHY THIS SIMULATION?
# --------------------
# Earlier simulations showed that a fixed numerical BACO value does NOT
# represent a fixed substantive discrepancy when RR_d changes:
#
#       BACO = log(RR_c) / log(RR_d)
#
# and
#
#       Delta = log(RR_c) - log(RR_d)
#             = (BACO - 1) * log(RR_d).
#
# Therefore, when RR_d -> 1, even a very large BACO can correspond to a very
# small absolute discrepancy between the composite and critical-outcome effects.
#
# This final sensitivity simulation fixes Delta directly:
#
#       Delta = log(RR_c) - log(RR_d)
#
# so that:
#
#       RR_c = RR_d * exp(Delta)
#
# and the corresponding BACO is:
#
#       BACO = 1 + Delta / log(RR_d),     RR_d != 1.
#
# An equivalent and often more intuitive parameter is:
#
#       QRR = [P(D | C, treatment)] / [P(D | C, control)]
#           = RR_d / RR_c
#           = exp(-Delta).
#
# Thus:
#
#       H0: no composite-critical discrepancy
#       <=> Delta = 0
#       <=> QRR = 1
#       <=> RR_c = RR_d.
#
#
# OBJECTIVES
# ----------
# 1) Confirm that power is primarily driven by Delta/QRR, not by BACO magnitude.
# 2) Evaluate robustness to:
#       - unequal allocation (1:1 and 2:1);
#       - smaller total sample sizes;
#       - different composite risks;
#       - sparse versus more informative critical-event counts;
#       - strong versus weak RR_d.
# 3) Compare:
#       - Fisher exact;
#       - Fisher mid-P;
#       - conditional LR;
#       - Delta Wald;
#       - original BACO Wald.
# 4) Report BOTH nominal power and size-adjusted power.
#
#
# IMPORTANT
# ---------
# Size-adjusted power is used only to compare procedures fairly in simulation.
# The empirically calibrated p-value thresholds are NOT proposed for real-world
# use.
#
#
# DEFAULT COMPUTATION
# -------------------
# This is a moderately large final-sensitivity grid.
#
#   NSIM_NULL <- 5000
#   NSIM_ALT  <- 3000
#
# These replicate counts are the frozen manuscript settings.
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
OUTPUT_DIR <- file.path(ROOT_DIR, "results", "simulation_final_delta")
dir.create(
  OUTPUT_DIR,
  showWarnings = FALSE,
  recursive = TRUE
)

SEED <- 20260902

NSIM_NULL <- 5000
NSIM_ALT <- 3000

ALPHA <- 0.05

# Total randomized sample size.
N_TOTAL <- c(
  500,
  2000,
  10000
)

# Treatment : control allocation ratio.
ALLOCATION_RATIO <- c(
  1, # 1:1
  2 # 2:1
)

# Control-arm composite risk.
P_COMPOSITE_CONTROL <- c(
  0.05,
  0.25,
  0.50
)

# Expected total critical events across both groups.
EXPECTED_TOTAL_DEATHS <- c(
  2,
  5,
  10,
  20
)

# Strong and weak reference effects.
RR_D_REFERENCE <- c(
  0.50,
  0.95,
  1.05,
  2.00
)

# Fixed effect discrepancies.
#
# Delta < 0 => RR_c < RR_d
# Delta > 0 => RR_c > RR_d
DELTA_ALT <- c(
  -0.50,
  -0.25,
  -0.10,
  0.10,
  0.25,
  0.50
)

NUMERIC_TOL <- sqrt(.Machine$double.eps)

METHODS <- c(
  "Fisher",
  "MidP",
  "Conditional_LR",
  "Delta",
  "BACO_Wald"
)


# ------------------------------------------------------------------------------
# 2. Build base nuisance-parameter grid
# ------------------------------------------------------------------------------

base_grid <- tidyr::crossing(
  n_total = N_TOTAL,
  allocation_ratio = ALLOCATION_RATIO,
  p_composite_control = P_COMPOSITE_CONTROL,
  expected_total_deaths_target = EXPECTED_TOTAL_DEATHS,
  rr_d = RR_D_REFERENCE
) %>%
  mutate(
    # Treatment : control = r : 1
    n_control = round(
      n_total /
        (1 +
          allocation_ratio)
    ),

    n_treatment = n_total -
      n_control,

    allocation_label = ifelse(
      allocation_ratio == 1,
      "1:1",
      "2:1"
    ),

    denominator_class = ifelse(
      rr_d %in%
        c(
          0.95,
          1.05
        ),
      "weak",
      "strong"
    ),

    # Expected total deaths:
    #
    # E[D] = N_c*p_d0 + N_t*p_d1
    #      = p_d0 * (N_c + N_t*RR_d)
    #
    # => p_d0 = E[D] / (N_c + N_t*RR_d)
    p_death_control = expected_total_deaths_target /
      (n_control +
        n_treatment *
          rr_d),

    p_death_treatment = rr_d *
      p_death_control,

    # Under H0: Delta = 0, RR_c = RR_d.
    null_rr_c = rr_d,

    null_p_composite_treatment = rr_d *
      p_composite_control,

    expected_composite_control = n_control *
      p_composite_control,

    expected_composite_treatment_null = n_treatment *
      null_p_composite_treatment,

    valid_null = n_control > 0 &
      n_treatment > 0 &
      p_death_control >= 0 &
      p_death_control <= p_composite_control &
      p_death_treatment >= 0 &
      p_death_treatment <= null_p_composite_treatment &
      null_p_composite_treatment > 0 &
      null_p_composite_treatment < 1
  ) %>%
  filter(
    valid_null
  ) %>%
  mutate(
    base_scenario_id = row_number()
  ) %>%
  select(
    base_scenario_id,
    everything()
  )

write_csv(
  base_grid,
  file.path(
    OUTPUT_DIR,
    "00_base_scenario_grid.csv"
  )
)


# ------------------------------------------------------------------------------
# 3. Build Delta-parameterized alternative grid
# ------------------------------------------------------------------------------

alt_grid <- tidyr::crossing(
  base_scenario_id = base_grid$base_scenario_id,

  population_Delta = DELTA_ALT
) %>%
  left_join(
    base_grid,
    by = "base_scenario_id"
  ) %>%
  mutate(
    # Fixed discrepancy parameterization.
    rr_c = rr_d *
      exp(
        population_Delta
      ),

    p_composite_treatment = rr_c *
      p_composite_control,

    # Equivalent conditional critical-risk ratio.
    population_QRR = exp(
      -population_Delta
    ),

    # BACO implied by the fixed Delta.
    population_BACO = 1 +
      population_Delta /
        log(
          rr_d
        ),

    abs_population_BACO_minus_1 = abs(
      population_BACO -
        1
    ),

    abs_population_Delta = abs(
      population_Delta
    ),

    expected_composite_treatment = n_treatment *
      p_composite_treatment,

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
    "01_delta_parameterized_alternative_grid.csv"
  )
)

cat(
  "\nBase null scenarios: ",
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
# 4. Fisher probability-ordered p-value + mid-P
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

  # No critical events => no information against the null.
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

  probs <-
    dhyper(
      support,
      m = c_t,
      n = c_c,
      k = d_total
    )

  p_obs <-
    dhyper(
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
# 5. Conditional LR
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

  # Delta Wald.
  delta_hat <-
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
      delta_hat[
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
# 7. Simulate p-values for one scenario
# ------------------------------------------------------------------------------

simulate_pvalues <- function(
  n_treatment,
  n_control,
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
      "Invalid multinomial probability structure."
    )
  }

  set.seed(
    seed
  )

  draw_c <-
    rmultinom(
      n = nsim,

      size = n_control,

      prob = prob_c
    )

  draw_t <-
    rmultinom(
      n = nsim,

      size = n_treatment,

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

  lr_p <-
    conditional_lr_vector(
      d_t = death_t,

      c_t = comp_t,

      d_c = death_c,

      c_c = comp_c
    )

  asym <-
    asymptotic_pvalues(
      n_t = n_treatment,

      n_c = n_control,

      comp_t = comp_t,

      comp_c = comp_c,

      death_t = death_t,

      death_c = death_c
    )

  tibble(
    Fisher = exact$Fisher,

    MidP = exact$MidP,

    Conditional_LR = lr_p,

    Delta = asym$Delta,

    BACO_Wald = asym$BACO_Wald,

    zero_death_any_arm = death_t == 0 |
      death_c == 0,

    zero_death_both_arms = death_t == 0 &
      death_c == 0,

    composite_zero_any_arm = comp_t == 0 |
      comp_c == 0
  )
}


# ------------------------------------------------------------------------------
# 8. Empirical size-adjustment helper
# ------------------------------------------------------------------------------

calibrated_cutoff <- function(
  p_values,
  alpha = ALPHA
) {
  # Non-estimable procedure = no rejection.
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

  empirical_size <-
    vapply(
      candidates,
      function(x) {
        mean(
          p <= x
        )
      },
      numeric(1)
    )

  valid <-
    which(
      empirical_size <=
        alpha +
          1e-12
    )

  if (
    length(
      valid
    ) ==
      0
  ) {
    return(
      c(
        cutoff = 0,
        achieved_size = 0
      )
    )
  }

  j <-
    max(
      valid
    )

  c(
    cutoff = candidates[j],

    achieved_size = empirical_size[j]
  )
}


# ------------------------------------------------------------------------------
# 9. Null calibration
# ------------------------------------------------------------------------------

cat(
  "\nRunning null calibration...\n"
)

null_calibration_list <-
  vector(
    "list",
    nrow(
      base_grid
    )
  )

for (i in seq_len(
  nrow(base_grid)
)) {
  s <-
    base_grid[i, ]

  if (
    i == 1 ||
      i %% 25 == 0 ||
      i == nrow(base_grid)
  ) {
    cat(
      "Null ",
      i,
      " / ",
      nrow(base_grid),
      "\n",
      sep = ""
    )
  }

  null_p <-
    simulate_pvalues(
      n_treatment = s$n_treatment,

      n_control = s$n_control,

      p_composite_control = s$p_composite_control,

      p_composite_treatment = s$null_p_composite_treatment,

      p_death_control = s$p_death_control,

      p_death_treatment = s$p_death_treatment,

      nsim = NSIM_NULL,

      seed = SEED +
        1009L *
          s$base_scenario_id
    )

  null_calibration_list[[i]] <-
    map_dfr(
      METHODS,
      function(method_name) {
        p_raw <-
          null_p[[method_name]]

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

          method = method_name,

          null_estimable_fraction = mean(
            !is.na(p_raw)
          ),

          nominal_type1 = mean(
            p_proc < ALPHA
          ),

          calibrated_cutoff = cal["cutoff"],

          calibrated_size = cal["achieved_size"]
        )
      }
    )
}

null_calibration <-
  bind_rows(
    null_calibration_list
  ) %>%
  left_join(
    base_grid %>%
      select(
        base_scenario_id,
        n_total,
        n_treatment,
        n_control,
        allocation_label,
        denominator_class,
        rr_d,
        p_composite_control,
        expected_total_deaths_target,
        expected_composite_control,
        expected_composite_treatment_null
      ),
    by = "base_scenario_id"
  )

write_csv(
  null_calibration,
  file.path(
    OUTPUT_DIR,
    "02_null_calibration_final_sensitivity.csv"
  )
)


# ------------------------------------------------------------------------------
# 10. Alternative power simulation
# ------------------------------------------------------------------------------

cat(
  "\nRunning Delta-parameterized alternatives...\n"
)

power_list <-
  vector(
    "list",
    nrow(
      alt_grid
    )
  )

for (i in seq_len(
  nrow(alt_grid)
)) {
  s <-
    alt_grid[i, ]

  if (
    i == 1 ||
      i %% 50 == 0 ||
      i == nrow(alt_grid)
  ) {
    cat(
      "Alternative ",
      i,
      " / ",
      nrow(alt_grid),
      "\n",
      sep = ""
    )
  }

  alt_p <-
    simulate_pvalues(
      n_treatment = s$n_treatment,

      n_control = s$n_control,

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

  power_list[[i]] <-
    map_dfr(
      METHODS,
      function(method_name) {
        p_raw <-
          alt_p[[method_name]]

        p_proc <-
          ifelse(
            is.na(p_raw),
            1,
            p_raw
          )

        cutoff <-
          cal_s %>%
          filter(
            method == method_name
          ) %>%
          pull(
            calibrated_cutoff
          )

        if (length(cutoff) != 1) {
          stop(
            "Could not identify unique size-adjustment cutoff."
          )
        }

        tibble(
          alternative_scenario_id = s$alternative_scenario_id,

          base_scenario_id = s$base_scenario_id,

          method = method_name,

          estimable_fraction = mean(
            !is.na(p_raw)
          ),

          nominal_power = mean(
            p_proc < ALPHA
          ),

          size_adjusted_power = mean(
            p_proc <= cutoff
          ),

          calibrated_cutoff = cutoff
        )
      }
    )
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

        n_total,
        n_treatment,
        n_control,
        allocation_ratio,
        allocation_label,

        denominator_class,
        rr_d,
        rr_c,

        population_Delta,
        abs_population_Delta,
        population_QRR,
        population_BACO,
        abs_population_BACO_minus_1,

        p_composite_control,
        p_composite_treatment,

        expected_total_deaths_target,
        expected_composite_control,
        expected_composite_treatment
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
    "03_power_final_delta_sensitivity.csv"
  )
)


# ------------------------------------------------------------------------------
# 11. Key summaries
# ------------------------------------------------------------------------------

summary_type1 <-
  null_calibration %>%
  group_by(
    method,
    allocation_label,
    denominator_class
  ) %>%
  summarise(
    n_scenarios = n(),

    mean_type1 = mean(
      nominal_type1
    ),

    max_type1 = max(
      nominal_type1
    ),

    mean_calibrated_size = mean(
      calibrated_size
    ),

    median_estimable_fraction = median(
      null_estimable_fraction
    ),

    .groups = "drop"
  )

write_csv(
  summary_type1,
  file.path(
    OUTPUT_DIR,
    "04_type1_summary_allocation_denominator.csv"
  )
)


summary_power_delta <-
  power_results %>%
  group_by(
    method,
    population_Delta,
    population_QRR,
    denominator_class
  ) %>%
  summarise(
    n_scenarios = n(),

    median_population_BACO = median(
      population_BACO
    ),

    mean_estimable_fraction = mean(
      estimable_fraction
    ),

    mean_nominal_power = mean(
      nominal_power
    ),

    mean_size_adjusted_power = mean(
      size_adjusted_power
    ),

    .groups = "drop"
  )

write_csv(
  summary_power_delta,
  file.path(
    OUTPUT_DIR,
    "05_power_summary_by_Delta_QRR.csv"
  )
)


summary_power_design <-
  power_results %>%
  group_by(
    method,
    n_total,
    allocation_label,
    p_composite_control,
    expected_total_deaths_target
  ) %>%
  summarise(
    n_scenarios = n(),

    mean_estimable_fraction = mean(
      estimable_fraction
    ),

    mean_size_adjusted_power = mean(
      size_adjusted_power
    ),

    .groups = "drop"
  )

write_csv(
  summary_power_design,
  file.path(
    OUTPUT_DIR,
    "06_power_summary_by_design.csv"
  )
)


# ------------------------------------------------------------------------------
# 12. Fisher-focused robustness summary
# ------------------------------------------------------------------------------

fisher_robustness <-
  power_results %>%
  filter(
    method == "Fisher"
  ) %>%
  group_by(
    population_Delta,
    population_QRR,
    expected_total_deaths_target,
    n_total,
    allocation_label
  ) %>%
  summarise(
    n_scenarios = n(),

    mean_size_adjusted_power = mean(
      size_adjusted_power
    ),

    min_size_adjusted_power = min(
      size_adjusted_power
    ),

    max_size_adjusted_power = max(
      size_adjusted_power
    ),

    .groups = "drop"
  )

write_csv(
  fisher_robustness,
  file.path(
    OUTPUT_DIR,
    "07_Fisher_robustness_summary.csv"
  )
)


# ------------------------------------------------------------------------------
# 13. Figure 1:
#     Same Delta, very different BACO values
# ------------------------------------------------------------------------------

concept_data <-
  alt_grid %>%
  distinct(
    rr_d,
    denominator_class,
    population_Delta,
    population_QRR,
    population_BACO
  )

p1 <- ggplot(
  concept_data,
  aes(
    x = population_Delta,

    y = population_BACO,

    shape = factor(
      rr_d
    )
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = 2
  ) +
  geom_vline(
    xintercept = 0,
    linetype = 2
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
    x = expression(
      Delta ==
        log(RR[c]) -
          log(RR[d])
    ),

    y = "Implied population BACO",

    shape = expression(
      RR[d]
    ),

    title = "A fixed effect discrepancy can imply radically different BACO values",

    subtitle = "The BACO scale becomes increasingly stretched as RRd approaches 1"
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
    "Figure_1_fixed_Delta_implied_BACO.png"
  ),

  plot = p1,

  width = 9,
  height = 6,
  dpi = 300
)


# ------------------------------------------------------------------------------
# 14. Figure 2:
#     Fisher power versus Delta under varying design conditions
# ------------------------------------------------------------------------------

fig2_data <-
  power_results %>%
  filter(
    method == "Fisher"
  ) %>%
  group_by(
    population_Delta,
    expected_total_deaths_target,
    n_total,
    allocation_label,
    denominator_class
  ) %>%
  summarise(
    mean_size_adjusted_power = mean(
      size_adjusted_power
    ),

    .groups = "drop"
  ) %>%
  mutate(
    deaths_facet = factor(
      expected_total_deaths_target,
      levels = EXPECTED_TOTAL_DEATHS,
      labels = paste0(
        "E[D] = ",
        EXPECTED_TOTAL_DEATHS
      )
    )
  )

p2 <- ggplot(
  fig2_data,
  aes(
    x = population_Delta,

    y = mean_size_adjusted_power,

    shape = factor(
      n_total
    )
  )
) +
  geom_line(
    aes(
      group = interaction(
        n_total,
        allocation_label,
        denominator_class
      )
    ),
    linewidth = 0.40
  ) +
  geom_point(
    size = 2.3
  ) +
  facet_grid(
    allocation_label ~
      deaths_facet
  ) +
  labs(
    x = expression(
      Delta
    ),

    y = "Fisher size-adjusted power",

    shape = "Total N",

    title = "Exact-test power is driven by discrepancy magnitude and critical-event information",

    subtitle = "Rows compare 1:1 versus 2:1 allocation; points average across composite risks and RRd classes"
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
    "Figure_2_Fisher_power_design_robustness.png"
  ),

  plot = p2,

  width = 13,
  height = 7,
  dpi = 300
)


# ------------------------------------------------------------------------------
# 15. Figure 3:
#     Fisher power by composite risk
# ------------------------------------------------------------------------------

fig3_data <-
  power_results %>%
  filter(
    method == "Fisher",
    n_total ==
      max(
        N_TOTAL
      )
  ) %>%
  group_by(
    population_Delta,
    p_composite_control,
    expected_total_deaths_target,
    allocation_label
  ) %>%
  summarise(
    mean_size_adjusted_power = mean(
      size_adjusted_power
    ),

    .groups = "drop"
  ) %>%
  mutate(
    composite_label = factor(
      p_composite_control,
      levels = P_COMPOSITE_CONTROL,
      labels = paste0(
        "P(C0) = ",
        P_COMPOSITE_CONTROL
      )
    )
  )

p3 <- ggplot(
  fig3_data,
  aes(
    x = expected_total_deaths_target,

    y = mean_size_adjusted_power,

    shape = factor(
      population_Delta
    )
  )
) +
  geom_line(
    aes(
      group = interaction(
        population_Delta,
        allocation_label
      )
    ),
    linewidth = 0.40
  ) +
  geom_point(
    size = 2.3
  ) +
  facet_grid(
    allocation_label ~
      composite_label
  ) +
  scale_x_continuous(
    breaks = EXPECTED_TOTAL_DEATHS
  ) +
  labs(
    x = "Expected total critical events",

    y = "Fisher size-adjusted power",

    shape = expression(
      Delta
    ),

    title = "Composite-event prevalence has limited impact once the critical-event process is fixed"
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
    "Figure_3_Fisher_power_composite_risk.png"
  ),

  plot = p3,

  width = 12,
  height = 7,
  dpi = 300
)


# ------------------------------------------------------------------------------
# 16. Figure 4:
#     Method comparison by Delta
# ------------------------------------------------------------------------------

fig4_data <-
  power_results %>%
  group_by(
    method,
    population_Delta,
    denominator_class,
    expected_total_deaths_target
  ) %>%
  summarise(
    mean_size_adjusted_power = mean(
      size_adjusted_power
    ),

    .groups = "drop"
  ) %>%
  filter(
    expected_total_deaths_target %in%
      c(
        5,
        20
      )
  ) %>%
  mutate(
    deaths_facet = factor(
      expected_total_deaths_target,
      levels = c(
        5,
        20
      ),
      labels = c(
        "E[D] = 5",
        "E[D] = 20"
      )
    )
  )

p4 <- ggplot(
  fig4_data,
  aes(
    x = population_Delta,

    y = mean_size_adjusted_power,

    shape = method
  )
) +
  geom_line(
    aes(
      group = method
    ),
    linewidth = 0.45
  ) +
  geom_point(
    size = 2.5
  ) +
  facet_grid(
    denominator_class ~
      deaths_facet
  ) +
  labs(
    x = expression(
      Delta
    ),

    y = "Mean size-adjusted power",

    shape = NULL,

    title = "After size calibration, power differences between valid procedures are modest"
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
    "Figure_4_method_comparison_fixed_Delta.png"
  ),

  plot = p4,

  width = 11,
  height = 7,
  dpi = 300
)


# ------------------------------------------------------------------------------
# 17. Monte Carlo precision
# ------------------------------------------------------------------------------

mc_precision <-
  tibble(
    component = c(
      "null calibration",
      "alternative power"
    ),

    nsim = c(
      NSIM_NULL,
      NSIM_ALT
    ),

    MCSE_at_probability_0_05 = sqrt(
      0.05 *
        0.95 /
        nsim
    ),

    MCSE_at_probability_0_50 = sqrt(
      0.50 *
        0.50 /
        nsim
    )
  )

write_csv(
  mc_precision,
  file.path(
    OUTPUT_DIR,
    "08_Monte_Carlo_precision.csv"
  )
)


# ------------------------------------------------------------------------------
# 18. Console summary
# ------------------------------------------------------------------------------

cat(
  "\n============================================================\n"
)

cat(
  "BACO final Delta-parameterized sensitivity simulation complete\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Type-I error summary:\n"
)

print(
  summary_type1
)

cat(
  "\nPower by Delta / QRR:\n"
)

print(
  summary_power_delta
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
# INTERPRETATION TARGETS
# ==============================================================================
#
# The final methodological narrative will be strongest if the simulation shows:
#
# 1) A FIXED Delta gives broadly comparable power across strong/weak RR_d once
#    critical-event information and design are held constant.
#
# 2) The implied BACO can differ dramatically for the SAME Delta:
#
#       BACO = 1 + Delta/log(RR_d).
#
#    This demonstrates why BACO magnitude alone should not be interpreted as
#    absolute discrepancy magnitude.
#
# 3) Fisher exact retains valid size under:
#       - sparse events;
#       - zero cells;
#       - smaller N;
#       - unequal allocation;
#       - different composite risks.
#
# 4) Mid-P may recover small amounts of power but can exceed nominal size.
#
# 5) Delta is useful with adequate nonzero event information but remains
#    unavailable when a critical cell is zero.
#
# 6) Original BACO Wald remains unsuitable because:
#       - it is parameterization-sensitive;
#       - it is weak-denominator unstable;
#       - it can be severely anti-conservative;
#       - apparent power can be driven by inflated type-I error.
#
#
# AFTER THIS SCRIPT
# -----------------
# If these targets are confirmed, stop expanding the simulation grid.
#
# Next steps:
#
#   A) freeze the methodological recommendations;
#   B) draft the Methods and Results of the methodological manuscript;
#   C) define the package API:
#
#          baco()
#          baco_test()
#          baco_fieller()
#          baco_diagnostics()
#
#   D) implement unit tests using the empirical 123-study dataset and selected
#      simulation edge cases.
#
# ==============================================================================
