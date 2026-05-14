# Simulation Study

This directory contains the simulation study scripts and results supporting the following
publication:

Bobak, C. A., Mohan, D., Murphy, M. A., Barnato, A. E., & O'Malley, A. J. (2026). Gaming
the system: evaluating spillover in a video game intervention for advanced care planning
using physician social networks. *Health Services and Outcomes Research Methodology*, 1–24.
https://doi.org/10.1007/s10742-025-00370-9

---

## Overview

The simulation framework mirrors the extended SW-CRT model from the paper. Each replicate
begins with a panel of step-indexed physician networks (one row-stochastic adjacency matrix
per trial step), a physician-hospital mapping, and each physician's first treated step.
Patient outcomes are generated from a logistic mixed-effects model with hospital random
intercepts, with parameters set to match the empirical data structure. Both a
network-adjusted model and a naïve stepped-wedge model (ignoring spillover) are fit to
each replicate, enabling evaluation of bias, RMSE, coverage, Type I error, and power
across a range of simulation conditions.

The scripts are intended to be run in sequence. `create_network_bundle.R` prepares the
network input; `simulation_block_code.R` defines all simulation functions;
`run_sim_code.Rmd` executes the simulation scenarios and writes results to CSV.
`simulation_code.R` is an earlier draft of the simulation functions retained for
reference.

---

## Script Descriptions

### 1. `create_network_bundle.R`

Converts the step-wise physician network objects (produced by `Create_Spillover_Networks.R`
in the parent HSORM Scripts directory) into the `network_panel` format required by the
simulation functions. Key operations include:

- Aligning vertex sets across all trial steps to a common NPI universe
- Converting `igraph` adjacency objects to row-stochastic sparse `dgCMatrix` format via
  the `Matrix` package, with isolated-node handling
- Assembling a named list containing `W_list` (step-indexed adjacency matrices),
  `npi`, `trial_step`, and `hospital` vectors
- Saving the bundle as `network_panel.rds` for use in simulation runs

**Dependencies:** `igraph`, `dplyr`, `Matrix`

---

### 2. `simulation_block_code.R`

Defines all functions used by the simulation runner. Key components include:

- **`phys_exposure()`** -- Computes physician-level spillover exposure under four
  data-generating mechanisms: proportional (default), threshold (any treated neighbor),
  count of treated peers, and log-sum of treated peers. Supports both strict-prior and
  non-strict-prior definitions of intervention status.
- **`sample_team()`** -- Samples a care team for each simulated patient stay, with
  configurable team size and a probability parameter controlling cross-hospital physician
  mixing.
- **`safe_glmer()`** -- Fits a GLMM via `lme4::glmer` with automatic fallback to a
  second optimizer on convergence failure, returning `NULL` on persistent failure rather
  than erroring.
- **`sim_one()`** -- Simulates a single dataset and fits both the network-adjusted model
  (`y ~ Post * Expose + stepF + (1 | hosp)`) and the naïve model (`y ~ Post + stepF +
  (1 | hosp)`). Accepts all key design and truth parameters as arguments.
- **`run_setting()`** -- Runs `R` replicates under a fixed parameter setting and
  computes summary statistics (mean estimate, bias, empirical SE, mean model SE, RMSE,
  coverage, Type I error, power, sign error rate, and convergence failure rate) for each
  target coefficient. Optionally writes results to CSV.

**Dependencies:** `Matrix`, `lme4`, `lmtest`, `sandwich`

---

### 3. `run_sim_code.Rmd`

Executes the simulation scenarios reported in the paper using the functions defined in
`simulation_block_code.R`. Three scenario sets are run with R = 1,000 replicates each:

- **Set A** -- Spillover sweep: `beta2` varied across {-1, 0, 1} with fixed `beta3`
- **Set B** -- Interaction sweep: `beta3` varied across {-1, 0, 1} with fixed positive
  spillover (`beta2 = log(2.8)`)
- **Set C** -- Misspecification scenario: data generated under a binary threshold
  exposure mechanism while the analysis model assumes the proportional form

Results from all three sets are combined and written to `simulation_results.csv`.

**Dependencies:** `Matrix`, `lme4`, `lmtest`, `sandwich`

---

### 4. `simulation_code.R`

An earlier prototype of the simulation functions, retained for reference. Implements
the core `sim_one()` and `sim_many()` functions with a simpler exposure engine
(proportional only) and no convergence fallback. Includes `make_blocky_W_list()`, a
self-contained toy network generator useful for local testing without access to the
secure compute environment.

**Dependencies:** `Matrix`, `lme4`, `lmtest`, `sandwich`

---

### 5. `simulation_results.csv`

Pre-computed simulation results for all scenarios described above, as reported in the
paper. Columns include `term`, `model`, `truth`, `mean_est`, `bias`, `emp_se`,
`mean_se`, `rmse`, `coverage`, `type1`, `power`, `sign_error`, `n_effective`,
`conv_fail_rate_net`, `conv_fail_rate_naive`, and `setting`.

---

## Pipeline Summary

```
Create_Spillover_Networks.R  (HSORM Scripts)
        ↓
  phys_network_list.RData + npi_trial_data.csv
        ↓
create_network_bundle.R
        ↓
  network_panel.rds
        ↓
simulation_block_code.R  (function definitions)
        ↓
run_sim_code.Rmd
        ↓
  simulation_results.csv
```

---

## Notes

- All scripts that read data reference paths within the secure computing environment
  and will require path updates to run in other environments.
- `simulation_code.R` includes `make_blocky_W_list()`, which generates a synthetic
  block-structured network and can be used to test the simulation pipeline locally
  without access to the ACP trial data.
- The default truth parameters match those used in the paper:
  `beta1 = log(1.15)`, `beta2 = log(2.8)`, `beta3 = log(0.88)`,
  baseline ACP probability ≈ 0.22, hospital RE SD = 0.35.
- R = 1,000 replicates per scenario were used for the published results; Monte Carlo
  standard errors for proportions are reported as √p(1-p)/R.

---

## Citation

Bobak, C. A., Mohan, D., Murphy, M. A., Barnato, A. E., & O'Malley, A. J. (2026). Gaming
the system: evaluating spillover in a video game intervention for advanced care planning
using physician social networks. *Health Services and Outcomes Research Methodology*, 1–24.
https://doi.org/10.1007/s10742-025-00370-9
