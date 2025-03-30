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

```
git clone https://github.com/CarlyBobak/ANSWER.git
```

---

## Repository Contents

- **`Create_Spillover_Networks.R`**  
  Generates physician social networks from billing data and computes a time-varying spillover metric for each physician across trial steps.

- **`Make_Spillover_Variables.R`**  
  Aggregates physician-level spillover exposure to the patient level based on encounters during a hospital stay. Output is formatted for use in downstream statistical modeling.

- **`Spillover_Models.Rmd`**  
  Contains analysis code for fitting models using the spillover variables and formatting model results for interpretation and reporting.

---

## Citation

A pre-print can be found here:

---

## Collaborate with Us

## 🤝 Collaborate with Us

We are actively looking for collaborators running **stepped-wedge** or other **cluster randomized trials** to pilot the **ANSWER** framework. If you're interested in exploring how network-informed spillover adjustment could enhance your trial analysis, we'd love to hear from you!

Reach out to **Carly Bobak** at [Carly.A.Bobak@dartmouth.edu](mailto:Carly.A.Bobak@dartmouth.edu) to learn more or discuss a potential collaboration.


