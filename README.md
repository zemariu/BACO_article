# Reproducibility files for the BACO methodological study

**Manuscript:** da Silva, JMN. *Wald inference for the BACO ratio is poorly calibrated when the critical effect is near the null: a simulation and empirical study*.

This repository contains the curated analysis code and the standardized comparison-level dataset used for the final manuscript. The public workflow is intentionally limited to analyses reported in the manuscript or Supplementary Material.

## Software

The analyses were developed in **R 4.5.1**. The empirical analysis requires the open-source **`bacoR`** package, available at [https://github.com/zemariu/bacoR](https://github.com/zemariu/bacoR). Other required packages are checked by each script and include `readxl`, `readr`, `dplyr`, `purrr`, `tibble`, `tidyr`, `ggplot2`, `scales`, `stringr`, `flextable`, and `officer`.

## Repository structure

The repository is organised as follows:

```text
BACO_article/
├── README.md
├── .gitignore
├── scripts/
│   ├── 01_empirical_reanalysis.R
│   ├── 02_final_delta_simulation.R
│   ├── 03_zero_event_simulation.R
│   ├── 04_make_tables_figures.R
│   ├── 99_run_all.R
│   └── supplementary_sensitivity/
│       ├── S01_initial_type1_calibration.R
│       ├── S02_confirmatory_null_simulation.R
│       └── S03_baco_parameterized_power.R
├── frozen_results/
│   └── simulation/
├── data/
│   └── Supplementary_Data_Studies_standardized.xlsx
└── results/
```

## Public input data

`data/Supplementary_Data_Studies_standardized.xlsx`

The workbook contains one standardized `Data` sheet for both empirical datasets, plus `README` and `Codebook` sheets. Each row is a randomized treatment/intervention-control comparison. The same variable names are used for both datasets. Article-level descriptive information already available in the published source reviews/trial reports is not duplicated; first author and DOI are retained for linkage.

The empirical script recalculates all inferential quantities from the six arm-level counts (`n_treatment`, `n_control`, `composite_treatment`, `composite_control`, `critical_treatment`, `critical_control`) rather than using stored derived values.

## Final analysis scripts

1. `scripts/01_empirical_reanalysis.R`  
   Recalculates the primary empirical reanalysis and external COVID-19 replication from the standardized public workbook using `bacoR`. It reproduces BACO, Fieller sets, weak-denominator diagnostics, Delta, QRR, Fisher exact testing, Delta Wald inference, and the original BACO Wald/Taylor comparator. It also reproduces the prespecified external sensitivity sets.

2. `scripts/02_final_delta_simulation.R`  
   Final Delta/QRR-parameterized Monte Carlo experiment used for the main null calibration and power results. Frozen manuscript settings are 5,000 null and 3,000 alternative replicates per scenario.

3. `scripts/03_zero_event_simulation.R`  
   Separate 10,000-replicate sparse/zero-critical-event experiment comparing Fisher exact, mid-P, conditional LR, Delta Wald, and original BACO Wald procedures.

4. `scripts/04_make_tables_figures.R`  
   Produces the final manuscript Tables 1-2, Figures 1-3, Supplementary Tables S1-S8, and the two supplementary figures used in the current manuscript: Figure S1 (sparse-event type I error) and Figure S2 (BACO-Delta relationship).

5. `scripts/99_run_all.R`  
   Convenience workflow. By default it recomputes the empirical analyses and uses the frozen final simulation CSVs, avoiding a lengthy Monte Carlo rerun. Set `RUN_LONG_SIMULATIONS <- TRUE` to regenerate the final simulation outputs first.

## Supplementary sensitivity scripts

The following scripts correspond to simulation stages described in the Supplementary Methods but are not required to recreate the main manuscript tables and figures:

- `scripts/supplementary_sensitivity/S01_initial_type1_calibration.R` — initial pilot calibration of sparse-event and near-null-denominator behavior.
- `scripts/supplementary_sensitivity/S02_confirmatory_null_simulation.R` — confirmatory separation of regular, sparse, weak-denominator, and combined settings.
- `scripts/supplementary_sensitivity/S03_baco_parameterized_power.R` — initial power experiment parameterized by BACO, used to demonstrate why fixed BACO alternatives confound discrepancy magnitude with denominator strength.

## Frozen simulation outputs

`frozen_results/simulation/` contains the exact CSV inputs used by the final table/figure script:

- `02_null_calibration_final_sensitivity.csv`
- `03_power_final_delta_sensitivity.csv`
- `01_zero_event_test_simulation_results.csv`

These files allow exact regeneration of reported tables and figures without repeating computationally intensive simulations. They are outputs, not additional empirical data sources.

## Recommended workflow

From the repository root:

```r
source("scripts/99_run_all.R")
```

For a full Monte Carlo rerun, edit `scripts/99_run_all.R` and set:

```r
RUN_LONG_SIMULATIONS <- TRUE
```

The main analysis outputs are written under `results/`, and publication-ready tables/figures under `manuscript_outputs/`.


## Reproducibility checks

The empirical script contains frozen checks for the manuscript datasets and principal results, including:

- primary dataset: 130 comparisons from 123 studies; 129 with valid aggregate counts;
- external dataset: 24 comparisons from 23 studies; 22 comparisons from 21 studies in the primary external analysis;
- availability/estimability of BACO, Delta, Fieller, and Fisher procedures;
- weak-denominator and Fieller-set counts;
- rejection counts at alpha = 0.05 and 0.005;
- external direct-count-only and adverse-reorientation sensitivity analyses;
- no DOI overlap between the two empirical datasets.

The script also verifies selected recalculated quantities against the derived columns distributed in the public workbook.

## Contact
For questions, reproducibility issues, or technical problems, please contact:

José Mário Nunes da Silva, PhD

Email: [zemariu@hotmail.com](zemariu@hotmail.com)

## Licence
The source code in this repository may be distributed under the MIT License.

The MIT License applies to the original repository code and documentation. It does not override the terms of use, attribution requirements, or redistribution restrictions of the original mortality, population, contextual, or spatial datasets.
