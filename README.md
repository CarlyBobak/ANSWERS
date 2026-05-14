# ANSWERS: Adjustment for Network Spillover in Stepped-Wedge-based Experimental Rollouts and Studies

**ANSWERS** is an R-based toolkit for quantifying and adjusting for spillover effects in stepped-wedge clinical trials using physician social networks. It supports network-informed exposure calculation and provides outputs ready for use in regression models that account for contamination due to social proximity.

This package helps trialists *get ANSWERS* when treatment effects aren't contained by design.

---

## What It Does

- Computes weighted adjacency matrices for each trial step
- Calculates physician-level exposure to spillover using pre/post intervention data
- Supports summary of exposure at the **patient level** for those treated by multiple physicians
- Returns step-wise exposure vectors for integration into downstream statistical models

---

## Core Method

Spillover is estimated via a matrix multiplication of a normalized adjacency matrix and an indicator vector reflecting intervention status. Optional steps allow aggregation across multiple treating physicians per patient using customizable summary functions (e.g., sum, mean, max).

---

## Installation

Coming soon as an R package. For now, clone this repo:

```r
git clone https://github.com/CarlyBobak/ANSWERS.git
```

---

## Repository Contents

This repository contains both the in-development ANSWERS R package source and archived analysis scripts from prior publications.

### HSORM Scripts

Analysis scripts associated with Bobak et al. (2026), published in *Health Services and Outcomes Research Methodology*:

- **`Create_Spillover_Networks.R`** -- Generates physician social networks from billing data and computes a time-varying spillover metric for each physician across trial steps.
- **`Make_Spillover_Variables.R`** -- Aggregates physician-level spillover exposure to the patient level based on encounters during a hospital stay. Output is formatted for use in downstream statistical modeling.
- **`Spillover_Models.Rmd`** -- Contains analysis code for fitting models using the spillover variables and formatting model results for interpretation and reporting.
- **`Simulation Study`** -- Contains code and results from a simulation study evaluating the ANSWERS methodology.


### Causal Networks Scripts

Analysis scripts associated with O'Malley, Bobak & Barnato (2026), published in *Social Networks*:

- **`plots4paper.R`** -- Generates figures for the publication.
- **`retention_models_fe.Rmd`** -- Retention models with fixed effects specification.
- **`retention_models_potout.Rmd`** -- Retention models under potential outcomes framework.
- **`retention_models_tinv_fe.Rmd`** -- Time-invariant retention models with fixed effects.
- **`retention_models_tinv_potout.Rmd`** -- Time-invariant retention models under potential outcomes framework.

---

## Citations

If you use this repository or the methods herein, please cite the relevant publications:

Bobak, C. A., Mohan, D., Murphy, M. A., Barnato, A. E., & O'Malley, A. J. (2026). Gaming the system: evaluating spillover in a video game intervention for advanced care planning using physician social networks. *Health Services and Outcomes Research Methodology*, 1-24. https://doi.org/10.1007/s10742-025-00370-9

O'Malley, A. J., Bobak, C. A., & Barnato, A. E. (2026). Causal inference for intervention spillover in a stepped wedge cluster-randomized trial: Lessons from a physician network. *Social Networks*, 86, 431-445. https://doi.org/10.1016/j.socnet.2026.04.017

---

## Collaborate with Us

We are actively looking for collaborators running **stepped-wedge** or other **cluster randomized trials** to pilot the **ANSWERS** framework. If you're interested in exploring how network-informed spillover adjustment could enhance your trial analysis, we'd love to hear from you!

Reach out to **Carly Bobak** at <Carly.A.Bobak@dartmouth.edu> to learn more or discuss a potential collaboration.
