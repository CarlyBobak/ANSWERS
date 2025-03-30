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

---

## Citation

A pre-print can be found here:

---

## Contact

We are actively seeking additional stepped-wedge and other cluster randomized trials to pilot ANSWERS. Please don't hesitate to contact Carly Bobak at Carly.A.Bobak@dartmouth.edu if you're interested in trying out ANSWERS.

