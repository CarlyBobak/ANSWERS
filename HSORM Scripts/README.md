# HSORM Scripts

This directory contains the statistical analysis scripts for the following publication:

Bobak, C. A., Mohan, D., Murphy, M. A., Barnato, A. E., & O'Malley, A. J. (2026). Gaming
the system: evaluating spillover in a video game intervention for advanced care planning
using physician social networks. *Health Services and Outcomes Research Methodology*, 1–24.
https://doi.org/10.1007/s10742-025-00370-9

---

## Overview

Three scripts implement the end-to-end analysis pipeline: network construction from
physician billing data, aggregation of physician-level spillover exposure to the patient
level, and regression modeling to estimate intervention and spillover effects on advance
care planning (ACP) documentation. The scripts are intended to be run in sequence.

---

## Script Descriptions

### 1. `Create_Spillover_Networks.R`

Constructs physician social networks from Medicare billing encounter data and computes a
time-varying spillover exposure metric for each trial physician across each stepped-wedge
trial step. Key operations include:

- Filtering encounter data to the trial period and relevant physician NPIs
- Building a bipartite physician-by-patient-stay incidence matrix, processed in chunks
  to manage memory for large datasets
- Projecting the bipartite matrix to a weighted physician-by-physician adjacency matrix
  via shared patient counts
- Constructing `igraph` network objects and computing step-wise row-stochastic adjacency
  matrices
- Computing physician-level spillover exposure vectors for each trial step via matrix
  multiplication of the normalized adjacency matrix and a binary intervention status vector
- Writing physician-level spillover effects to CSV for downstream use

**Dependencies:** `dplyr`, `Matrix`, `igraph`

---

### 2. `Make_Spillover_Variables.R`

Aggregates the physician-level spillover exposure computed in script (1) to the patient
level. For each patient stay, the script identifies all physicians encountered during the
hospitalization and sums their spillover exposure values at the corresponding trial step.
The resulting patient-level spillover variable (`physWPost`) is appended to the analysis
dataset and written to CSV for use in downstream modeling.

**Dependencies:** base R

---

### 3. `Spillover_Models.Rmd`

Fits a series of mixed-effect logistic regression models using the patient-level spillover
exposure variable constructed in script (2), with ACP billing as the binary outcome.
Model specifications progress from a traditional stepped-wedge model (no spillover) through
network-adjusted models with and without intervention-by-spillover interaction terms.
Results are formatted and written to CSV files for reporting.

**Dependencies:** `lme4`, `nlme`, `geepack`, `Gmisc`, `dplyr`, `ggplot2`, `ggpubr`,
`broom.mixed`

---

## Pipeline Summary

```
Create_Spillover_Networks.R
        ↓
  physician_level_spillover_effects_rw_strict_prior.csv
        ↓
Make_Spillover_Variables.R
        ↓
  acp_data_with_spill_strict_prior_stochastic.csv
        ↓
Spillover_Models.Rmd
        ↓
  model output CSVs
```

---

## Notes

- Input data (encounter records and trial metadata) reside in a secure computing
  environment and are not included in this repository.
- The chunked matrix multiplication in `Create_Spillover_Networks.R` is designed for
  datasets with hundreds of thousands of patient stays; chunk boundaries may need
  adjustment depending on available memory.
- The spillover variable is scaled so that a unit change corresponds to a 10 percentage
  point increase in the proportion of peer physicians who have been intervened on,
  facilitating comparability with the direct intervention effect estimate.

---

## Citation

Bobak, C. A., Mohan, D., Murphy, M. A., Barnato, A. E., & O'Malley, A. J. (2026). Gaming
the system: evaluating spillover in a video game intervention for advanced care planning
using physician social networks. *Health Services and Outcomes Research Methodology*, 1–24.
https://doi.org/10.1007/s10742-025-00370-9
