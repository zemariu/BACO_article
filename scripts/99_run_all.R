# =============================================================================
# BACO JCE reproducibility workflow
# =============================================================================
# Run this script from the repository root.
#
# By default, the empirical analyses are recomputed from the public workbook,
# while the manuscript tables/figures use the frozen final simulation CSVs
# supplied in frozen_results/simulation/. Set RUN_LONG_SIMULATIONS <- TRUE to
# regenerate the final Monte Carlo outputs before creating tables and figures.
# =============================================================================

ROOT_DIR <- getwd()

RUN_LONG_SIMULATIONS <- FALSE
RUN_SUPPLEMENTARY_SENSITIVITY <- FALSE

source(file.path(ROOT_DIR, "scripts", "01_empirical_reanalysis.R"))

if (RUN_LONG_SIMULATIONS) {
  source(file.path(ROOT_DIR, "scripts", "02_final_delta_simulation.R"))
  source(file.path(ROOT_DIR, "scripts", "03_zero_event_simulation.R"))
}

if (RUN_SUPPLEMENTARY_SENSITIVITY) {
  source(file.path(
    ROOT_DIR, "scripts", "supplementary_sensitivity",
    "S01_initial_type1_calibration.R"
  ))
  source(file.path(
    ROOT_DIR, "scripts", "supplementary_sensitivity",
    "S02_confirmatory_null_simulation.R"
  ))
  source(file.path(
    ROOT_DIR, "scripts", "supplementary_sensitivity",
    "S03_baco_parameterized_power.R"
  ))
}

source(file.path(ROOT_DIR, "scripts", "04_make_tables_figures.R"))

cat("\nReproducibility workflow completed.\n")
