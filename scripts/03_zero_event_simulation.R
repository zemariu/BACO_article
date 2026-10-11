# ==============================================================================
# BACO reproducibility analysis
# Sparse/zero-event Monte Carlo validation used in the JCE manuscript
#
# This simulation validates the key algebraic reformulation:
#
#   H0: BACO = 1
#   <=> RR_c = RR_d
#   <=> P(critical | composite, treatment)
#       =
#       P(critical | composite, control)
#
# Hence, among participants with the composite endpoint, the null can be tested
# using the 2x2 table:
#
#                         critical     noncritical composite
#   treatment                D_t              C_t - D_t
#   control                  D_c              C_c - D_c
#
# Methods compared
# ----------------
# 1) Fisher exact conditional test (probability-ordered two-sided p-value)
# 2) Two-sided mid-P version of Fisher's test (sensitivity comparator)
# 3) Conditional binomial likelihood-ratio test (chi-square approximation)
#    -- algebraically equivalent to the full joint-multinomial profile LR test
#       of theta = 1, but vastly simpler computationally.
# 4) Direct asymptotic Delta test:
#       Delta = log(RR_c) - log(RR_d)
# 5) Original BACO Wald/delta ratio test
#
# The focus is TYPE-I ERROR and ESTIMABILITY under sparse/zero critical events.
#
# IMPORTANT
# ---------
# - Fisher exact is available with one zero critical cell.
# - If both arms have zero critical events, Fisher p = 1 by construction:
#   the data contain no information about differential critical composition.
# - Mid-P is not exact in the strict coverage sense; it is included because
#   Fisher can be conservative with very sparse data.
# - The LR chi-square approximation is expected to be anti-conservative in
#   some extremely sparse scenarios and is included to quantify that behavior.
#
# Recommended first run:
#       NSIM <- 10000
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
OUTPUT_DIR <- file.path(ROOT_DIR, "results", "simulation_zero_event")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

SEED <- 20260902

NSIM <- 10000

ALPHA <- 0.05
ALPHA_STRICT <- 0.005

N_PER_ARM <- 10000

# Composite risks span uncommon to common composite endpoints.
P_COMPOSITE_CONTROL <- c(
  0.05,
  0.10,
  0.25,
  0.50
)

# Expected total critical events across both arms.
# This deliberately concentrates on the region where zero cells are common.
EXPECTED_TOTAL_DEATHS <- c(
  0.5,
  1,
  2,
  3,
  5,
  10,
  20
)

# Strong and weak reference effects.
# Under H0, RR_c = RR_d = R, so population BACO = 1.
RR_REFERENCE <- c(
  0.50,
  0.95,
  1.05,
  2.00
)

NUMERIC_TOL <- sqrt(.Machine$double.eps)


# ------------------------------------------------------------------------------
# 2. Scenario grid
# ------------------------------------------------------------------------------

scenario_grid <- tidyr::crossing(
  n_per_arm = N_PER_ARM,
  p_composite_control = P_COMPOSITE_CONTROL,
  expected_total_deaths_target = EXPECTED_TOTAL_DEATHS,
  rr_reference = RR_REFERENCE
) %>%
  mutate(
    # E[D_total] = N*p_d0*(1 + RR_d)
    p_death_control =
      expected_total_deaths_target /
      (
        n_per_arm *
        (1 + rr_reference)
      ),

    p_death_treatment =
      rr_reference *
      p_death_control,

    p_composite_treatment =
      rr_reference *
      p_composite_control,

    population_beta_reference =
      log(
        rr_reference
      ),

    population_BACO = 1,

    valid_probability_structure =
      p_death_control >= 0 &
      p_death_control <= p_composite_control &
      p_death_treatment >= 0 &
      p_death_treatment <= p_composite_treatment &
      p_composite_treatment <= 1
  ) %>%
  filter(
    valid_probability_structure
  ) %>%
  mutate(
    scenario_id =
      row_number(),

    denominator_class =
      ifelse(
        rr_reference %in%
          c(
            0.95,
            1.05
          ),
        "weak",
        "strong"
      )
  ) %>%
  select(
    scenario_id,
    everything()
  )

write_csv(
  scenario_grid,
  file.path(
    OUTPUT_DIR,
    "00_zero_event_test_scenario_grid.csv"
  )
)


# ------------------------------------------------------------------------------
# 3. Probability-ordered two-sided Fisher p-value
# ------------------------------------------------------------------------------

# For fixed:
#   C_t = treatment composite events
#   C_c = control composite events
#   D   = total critical events
#
# under H0 q_t = q_c, the treatment critical count has:
#
#   D_t | C_t, C_c, D ~ Hypergeometric(C_t, C_c, D).
#
# The two-sided Fisher p-value used here matches the probability-ordering
# principle used by fisher.test(): sum probabilities no greater than the
# probability of the observed table.

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

  # No critical events: there is no information against H0.
  if (
    d_total == 0
  ) {
    return(
      c(
        exact = 1,
        midp = 1
      )
    )
  }

  lo <- max(
    0,
    d_total - c_c
  )

  hi <- min(
    c_t,
    d_total
  )

  support <- lo:hi

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

  # Floating-point tolerance analogous to standard exact-test implementations.
  exact_p <- sum(
    probs[
      probs <=
        p_obs *
        (1 + 1e-12)
    ]
  )

  exact_p <- min(
    1,
    exact_p
  )

  # Probability-ordered two-sided mid-P:
  # exact p minus half the observed-table probability.
  mid_p <- max(
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


# ------------------------------------------------------------------------------
# 4. Vectorized exact p-values
# ------------------------------------------------------------------------------

exact_p_vector <- function(
    d_t,
    c_t,
    d_c,
    c_c
) {

  n <- length(
    d_t
  )

  exact <- rep(
    NA_real_,
    n
  )

  midp <- rep(
    NA_real_,
    n
  )

  for (
    i in seq_len(
      n
    )
  ) {

    ans <- fisher_probability_ordered(
      d_t =
        d_t[i],

      c_t =
        c_t[i],

      c_c =
        c_c[i],

      d_total =
        d_t[i] +
        d_c[i]
    )

    exact[i] <-
      ans["exact"]

    midp[i] <-
      ans["midp"]
  }

  tibble(
    p_fisher =
      exact,

    p_midp =
      midp
  )
}


# ------------------------------------------------------------------------------
# 5. Conditional binomial LR test
# ------------------------------------------------------------------------------

binom_loglik <- function(
    d,
    c,
    q
) {

  if (
    c == 0
  ) {
    return(0)
  }

  if (
    q <= 0
  ) {
    return(
      ifelse(
        d == 0,
        0,
        -Inf
      )
    )
  }

  if (
    q >= 1
  ) {
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
    (c - d) *
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

  n <- length(
    d_t
  )

  LR <- rep(
    NA_real_,
    n
  )

  p <- rep(
    NA_real_,
    n
  )

  for (
    i in seq_len(
      n
    )
  ) {

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

    # If no critical events occur, unrestricted and null models coincide.
    if (
      d_total == 0
    ) {
      LR[i] <- 0
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

    LR[i] <-
      max(
        0,
        2 *
          (
            ll_u -
            ll_0
          )
      )

    p[i] <-
      pchisq(
        LR[i],
        df = 1,
        lower.tail = FALSE
      )
  }

  tibble(
    LR_conditional =
      LR,

    p_LR_conditional =
      p
  )
}


# ------------------------------------------------------------------------------
# 6. Delta and original BACO Wald tests
# ------------------------------------------------------------------------------

asymptotic_tests_vector <- function(
    n_t,
    n_c,
    comp_t,
    comp_c,
    death_t,
    death_c
) {

  nsim <- length(
    comp_t
  )

  p_delta <- rep(
    NA_real_,
    nsim
  )

  p_wald <- rep(
    NA_real_,
    nsim
  )

  baco <- rep(
    NA_real_,
    nsim
  )

  valid <-
    comp_t > 0 &
    comp_c > 0 &
    death_t > 0 &
    death_c > 0

  idx <- which(
    valid
  )

  if (
    length(
      idx
    ) == 0
  ) {
    return(
      tibble(
        BACO =
          baco,

        p_Delta =
          p_delta,

        p_BACO_Wald =
          p_wald
      )
    )
  }

  rr_c <-
    (comp_t[idx] / n_t) /
    (comp_c[idx] / n_c)

  rr_d <-
    (death_t[idx] / n_t) /
    (death_c[idx] / n_c)

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

  # Delta.
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
      delta[good_delta] /
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

  # Original BACO ratio Wald.
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

    covv <-
      cov_cd[
        good_ratio
      ]

    theta <-
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
      covv

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

      theta_good <-
        theta[
          good_var
        ]

      se_good <-
        sqrt(
          var_theta[
            good_var
          ]
        )

      baco[
        final_idx
      ] <-
        theta_good

      p_wald[
        final_idx
      ] <-
        2 *
        pnorm(
          -abs(
            (
              theta_good -
              1
            ) /
            se_good
          )
        )
    }
  }

  tibble(
    BACO =
      baco,

    p_Delta =
      p_delta,

    p_BACO_Wald =
      p_wald
  )
}


# ------------------------------------------------------------------------------
# 7. Monte Carlo summary helper
# ------------------------------------------------------------------------------

rate_summary <- function(
    p_values,
    alpha = ALPHA,
    total_n
) {

  valid <-
    !is.na(
      p_values
    )

  n_valid <-
    sum(
      valid
    )

  n_reject <-
    sum(
      p_values[
        valid
      ] <
        alpha
    )

  conditional <-
    if (
      n_valid > 0
    ) {
      n_reject /
        n_valid
    } else {
      NA_real_
    }

  unconditional <-
    n_reject /
    total_n

  mcse_unconditional <-
    sqrt(
      unconditional *
      (
        1 -
        unconditional
      ) /
      total_n
    )

  c(
    n_valid =
      n_valid,

    estimable_fraction =
      n_valid /
      total_n,

    rejection_conditional =
      conditional,

    rejection_unconditional =
      unconditional,

    mcse_unconditional =
      mcse_unconditional
  )
}


# ------------------------------------------------------------------------------
# 8. Simulate one scenario
# ------------------------------------------------------------------------------

simulate_one_scenario <- function(
    s,
    nsim = NSIM,
    seed = SEED
) {

  n_t <-
    s$n_per_arm

  n_c <-
    s$n_per_arm

  probs_c <- c(
    none =
      1 -
      s$p_composite_control,

    nonfatal =
      s$p_composite_control -
      s$p_death_control,

    critical =
      s$p_death_control
  )

  probs_t <- c(
    none =
      1 -
      s$p_composite_treatment,

    nonfatal =
      s$p_composite_treatment -
      s$p_death_treatment,

    critical =
      s$p_death_treatment
  )

  set.seed(
    seed +
    1009L *
    s$scenario_id
  )

  draw_c <- rmultinom(
    n = nsim,
    size = n_c,
    prob = probs_c
  )

  draw_t <- rmultinom(
    n = nsim,
    size = n_t,
    prob = probs_t
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

  total_deaths <-
    death_t +
    death_c

  zero_any <-
    death_t == 0 |
    death_c == 0

  zero_both <-
    death_t == 0 &
    death_c == 0

  exact <- exact_p_vector(
    d_t = death_t,
    c_t = comp_t,
    d_c = death_c,
    c_c = comp_c
  )

  lr <- conditional_lr_vector(
    d_t = death_t,
    c_t = comp_t,
    d_c = death_c,
    c_c = comp_c
  )

  asym <- asymptotic_tests_vector(
    n_t = n_t,
    n_c = n_c,
    comp_t = comp_t,
    comp_c = comp_c,
    death_t = death_t,
    death_c = death_c
  )

  r_fisher_05 <-
    rate_summary(
      exact$p_fisher,
      alpha =
        ALPHA,
      total_n =
        nsim
    )

  r_midp_05 <-
    rate_summary(
      exact$p_midp,
      alpha =
        ALPHA,
      total_n =
        nsim
    )

  r_lr_05 <-
    rate_summary(
      lr$p_LR_conditional,
      alpha =
        ALPHA,
      total_n =
        nsim
    )

  r_delta_05 <-
    rate_summary(
      asym$p_Delta,
      alpha =
        ALPHA,
      total_n =
        nsim
    )

  r_wald_05 <-
    rate_summary(
      asym$p_BACO_Wald,
      alpha =
        ALPHA,
      total_n =
        nsim
    )

  r_fisher_005 <-
    rate_summary(
      exact$p_fisher,
      alpha =
        ALPHA_STRICT,
      total_n =
        nsim
    )

  r_midp_005 <-
    rate_summary(
      exact$p_midp,
      alpha =
        ALPHA_STRICT,
      total_n =
        nsim
    )

  r_lr_005 <-
    rate_summary(
      lr$p_LR_conditional,
      alpha =
        ALPHA_STRICT,
      total_n =
        nsim
    )

  r_delta_005 <-
    rate_summary(
      asym$p_Delta,
      alpha =
        ALPHA_STRICT,
      total_n =
        nsim
    )

  r_wald_005 <-
    rate_summary(
      asym$p_BACO_Wald,
      alpha =
        ALPHA_STRICT,
      total_n =
        nsim
    )

  tibble(
    scenario_id =
      s$scenario_id,

    denominator_class =
      s$denominator_class,

    rr_reference =
      s$rr_reference,

    p_composite_control =
      s$p_composite_control,

    expected_total_deaths_target =
      s$expected_total_deaths_target,

    nsim =
      nsim,

    observed_mean_total_deaths =
      mean(
        total_deaths
      ),

    zero_death_any_arm_fraction =
      mean(
        zero_any
      ),

    zero_death_both_arms_fraction =
      mean(
        zero_both
      ),

    # alpha = .05
    fisher_type1_05 =
      r_fisher_05[
        "rejection_unconditional"
      ],

    fisher_mcse_05 =
      r_fisher_05[
        "mcse_unconditional"
      ],

    midp_type1_05 =
      r_midp_05[
        "rejection_unconditional"
      ],

    lr_type1_05 =
      r_lr_05[
        "rejection_unconditional"
      ],

    Delta_estimable_fraction =
      r_delta_05[
        "estimable_fraction"
      ],

    Delta_type1_05_conditional =
      r_delta_05[
        "rejection_conditional"
      ],

    Delta_type1_05_unconditional =
      r_delta_05[
        "rejection_unconditional"
      ],

    Wald_estimable_fraction =
      r_wald_05[
        "estimable_fraction"
      ],

    Wald_type1_05_conditional =
      r_wald_05[
        "rejection_conditional"
      ],

    Wald_type1_05_unconditional =
      r_wald_05[
        "rejection_unconditional"
      ],

    # alpha = .005
    fisher_type1_005 =
      r_fisher_005[
        "rejection_unconditional"
      ],

    midp_type1_005 =
      r_midp_005[
        "rejection_unconditional"
      ],

    lr_type1_005 =
      r_lr_005[
        "rejection_unconditional"
      ],

    Delta_type1_005_unconditional =
      r_delta_005[
        "rejection_unconditional"
      ],

    Wald_type1_005_unconditional =
      r_wald_005[
        "rejection_unconditional"
      ]
  )
}


# ------------------------------------------------------------------------------
# 9. Run all scenarios
# ------------------------------------------------------------------------------

cat(
  "\nRunning zero-event exact-test simulation...\n"
)

cat(
  "Scenarios: ",
  nrow(
    scenario_grid
  ),
  "\n",
  sep = ""
)

cat(
  "Replicates/scenario: ",
  NSIM,
  "\n\n",
  sep = ""
)

results <- map_dfr(
  seq_len(
    nrow(
      scenario_grid
    )
  ),
  function(i) {

    if (
      i == 1 ||
      i %% 10 == 0 ||
      i ==
        nrow(
          scenario_grid
        )
    ) {
      cat(
        "Scenario ",
        i,
        " / ",
        nrow(
          scenario_grid
        ),
        "\n",
        sep = ""
      )
    }

    simulate_one_scenario(
      scenario_grid[
        i,
      ],
      nsim =
        NSIM,
      seed =
        SEED
    )
  }
)

write_csv(
  results,
  file.path(
    OUTPUT_DIR,
    "01_zero_event_test_simulation_results.csv"
  )
)


# ------------------------------------------------------------------------------
# 10. Summary by expected critical-event count
# ------------------------------------------------------------------------------

summary_by_deaths <- results %>%
  group_by(
    expected_total_deaths_target,
    denominator_class
  ) %>%
  summarise(
    n_scenarios =
      n(),

    mean_zero_any =
      mean(
        zero_death_any_arm_fraction
      ),

    mean_zero_both =
      mean(
        zero_death_both_arms_fraction
      ),

    mean_Fisher_type1_05 =
      mean(
        fisher_type1_05
      ),

    max_Fisher_type1_05 =
      max(
        fisher_type1_05
      ),

    mean_midP_type1_05 =
      mean(
        midp_type1_05
      ),

    max_midP_type1_05 =
      max(
        midp_type1_05
      ),

    mean_LR_type1_05 =
      mean(
        lr_type1_05
      ),

    max_LR_type1_05 =
      max(
        lr_type1_05
      ),

    mean_Delta_estimable =
      mean(
        Delta_estimable_fraction
      ),

    mean_Delta_type1_cond =
      mean(
        Delta_type1_05_conditional,
        na.rm = TRUE
      ),

    mean_Wald_estimable =
      mean(
        Wald_estimable_fraction
      ),

    mean_Wald_type1_cond =
      mean(
        Wald_type1_05_conditional,
        na.rm = TRUE
      ),

    .groups =
      "drop"
  )

write_csv(
  summary_by_deaths,
  file.path(
    OUTPUT_DIR,
    "02_summary_by_expected_deaths.csv"
  )
)


# ------------------------------------------------------------------------------
# 11. Summary by method
# ------------------------------------------------------------------------------

long_results <- results %>%
  select(
    scenario_id,
    denominator_class,
    rr_reference,
    p_composite_control,
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
    names_to =
      "method",
    values_to =
      "type1"
  ) %>%
  mutate(
    method = recode(
      method,

      fisher_type1_05 =
        "Fisher exact",

      midp_type1_05 =
        "Fisher mid-P",

      lr_type1_05 =
        "Conditional LR",

      Delta_type1_05_unconditional =
        "Delta Wald",

      Wald_type1_05_unconditional =
        "Original BACO Wald"
    )
  )

write_csv(
  long_results,
  file.path(
    OUTPUT_DIR,
    "03_plot_ready_type1_long.csv"
  )
)


# ------------------------------------------------------------------------------
# 12. Figure 1: Type-I error across sparse event information
# ------------------------------------------------------------------------------

fig1_data <- long_results %>%
  filter(
    p_composite_control ==
      0.25
  )

p1 <- ggplot(
  fig1_data,
  aes(
    x =
      expected_total_deaths_target,

    y =
      type1,

    shape =
      method
  )
) +
  geom_hline(
    yintercept =
      ALPHA,
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
          rr_reference
        )
    ),
    linewidth = 0.45
  ) +
  facet_wrap(
    ~ denominator_class
  ) +
  scale_x_continuous(
    breaks =
      EXPECTED_TOTAL_DEATHS
  ) +
  labs(
    x =
      "Expected total critical events",

    y =
      "Empirical type-I error",

    shape =
      NULL,

    title =
      "Exact conditional testing controls false-positive error under sparse BACO settings",

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
  filename =
    file.path(
      OUTPUT_DIR,
      "Figure_1_type1_sparse_exact_methods.png"
    ),

  plot =
    p1,

  width = 11,
  height = 7,
  dpi = 300
)


# ------------------------------------------------------------------------------
# 13. Figure 2: Estimability of asymptotic methods
# ------------------------------------------------------------------------------

fig2_data <- results %>%
  filter(
    p_composite_control ==
      0.25
  )

p2 <- ggplot(
  fig2_data,
  aes(
    x =
      expected_total_deaths_target,

    y =
      Delta_estimable_fraction,

    shape =
      factor(
        rr_reference
      )
  )
) +
  geom_point(
    size = 2.5
  ) +
  geom_line(
    aes(
      group =
        rr_reference
    ),
    linewidth = 0.45
  ) +
  facet_wrap(
    ~ denominator_class
  ) +
  scale_x_continuous(
    breaks =
      EXPECTED_TOTAL_DEATHS
  ) +
  labs(
    x =
      "Expected total critical events",

    y =
      "Probability asymptotic Delta/BACO is calculable",

    shape =
      expression(
        RR[d]
      ),

    title =
      "Exact conditional inference avoids the zero-cell estimability barrier"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    legend.position =
      "bottom"
  )

ggsave(
  filename =
    file.path(
      OUTPUT_DIR,
      "Figure_2_estimability_vs_exact.png"
    ),

  plot =
    p2,

  width = 10,
  height = 6,
  dpi = 300
)


# ------------------------------------------------------------------------------
# 14. Figure 3: Fisher exact vs mid-P vs LR only
# ------------------------------------------------------------------------------

fig3_data <- long_results %>%
  filter(
    method %in%
      c(
        "Fisher exact",
        "Fisher mid-P",
        "Conditional LR"
      ),
    p_composite_control ==
      0.25
  )

p3 <- ggplot(
  fig3_data,
  aes(
    x =
      expected_total_deaths_target,

    y =
      type1,

    shape =
      method
  )
) +
  geom_hline(
    yintercept =
      ALPHA,
    linetype = 2
  ) +
  geom_point(
    size = 2.7
  ) +
  geom_line(
    aes(
      group =
        interaction(
          method,
          rr_reference
        )
    ),
    linewidth = 0.5
  ) +
  facet_wrap(
    ~ denominator_class
  ) +
  scale_x_continuous(
    breaks =
      EXPECTED_TOTAL_DEATHS
  ) +
  labs(
    x =
      "Expected total critical events",

    y =
      "Empirical type-I error",

    shape =
      NULL,

    title =
      "Trade-off between exactness and conservatism in sparse BACO testing"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    legend.position =
      "bottom"
  )

ggsave(
  filename =
    file.path(
      OUTPUT_DIR,
      "Figure_3_exact_midP_LR_comparison.png"
    ),

  plot =
    p3,

  width = 10,
  height = 6,
  dpi = 300
)


# ------------------------------------------------------------------------------
# 15. Worst-case calibration table
# ------------------------------------------------------------------------------

worst_cases <- results %>%
  arrange(
    desc(
      lr_type1_05
    )
  ) %>%
  select(
    scenario_id,
    denominator_class,
    rr_reference,
    p_composite_control,
    expected_total_deaths_target,

    zero_death_any_arm_fraction,

    fisher_type1_05,
    midp_type1_05,
    lr_type1_05,

    Delta_estimable_fraction,
    Delta_type1_05_conditional,

    Wald_estimable_fraction,
    Wald_type1_05_conditional
  ) %>%
  slice_head(
    n = 40
  )

write_csv(
  worst_cases,
  file.path(
    OUTPUT_DIR,
    "04_worst_sparse_calibration_cases.csv"
  )
)


# ------------------------------------------------------------------------------
# 16. Monte Carlo precision
# ------------------------------------------------------------------------------

mc_precision <- tibble(
  nsim =
    NSIM,

  nominal_probability =
    c(
      0.05,
      0.005
    ),

  MCSE =
    sqrt(
      nominal_probability *
      (
        1 -
        nominal_probability
      ) /
      NSIM
    )
)

write_csv(
  mc_precision,
  file.path(
    OUTPUT_DIR,
    "05_Monte_Carlo_precision.csv"
  )
)


# ------------------------------------------------------------------------------
# 17. Console
# ------------------------------------------------------------------------------

cat(
  "\n============================================================\n"
)

cat(
  "Sparse/zero-event BACO test simulation complete\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Summary by expected critical events:\n"
)

print(
  summary_by_deaths
)

cat(
  "\nWorst conditional-LR calibration scenarios:\n"
)

print(
  worst_cases %>%
    slice_head(
      n = 10
    )
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
# EXPECTED INTERPRETATION
# ==============================================================================
#
# If the empirical findings generalize:
#
# 1) Fisher exact:
#       - unconditional type-I error <= nominal level;
#       - conservative with very small critical-event counts;
#       - available even when one death cell is zero.
#
# 2) Fisher mid-P:
#       - closer to nominal 5%;
#       - improved power/calibration trade-off;
#       - may be mildly anti-conservative in selected scenarios and therefore
#         should remain a sensitivity method unless simulation strongly supports it.
#
# 3) Conditional LR:
#       - same null hypothesis as the full profile-LR test at theta = 1;
#       - likely anti-conservative when event counts are extremely low;
#       - should become well calibrated as information increases.
#
# 4) Delta:
#       - asymptotically well calibrated;
#       - unavailable with a zero death cell;
#       - conditional calibration alone is insufficient because estimability
#         can be very low.
#
# 5) Original BACO Wald:
#       - inherits both the zero-cell estimability problem and the previously
#         demonstrated weak-denominator anti-conservatism.
#
# Potential publication recommendation after validation:
#
#   Test H0: BACO = 1
#       sparse/zero critical events -> Fisher exact
#       adequate event information -> Delta/Fieller
#
#   Estimate BACO ratio / CI
#       use BACO point estimate when defined;
#       use Fieller confidence sets;
#       report weak-denominator/unbounded-set diagnostics;
#       do NOT use continuity correction as a default solution.
#
# ==============================================================================
