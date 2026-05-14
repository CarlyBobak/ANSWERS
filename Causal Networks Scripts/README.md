# Causal Networks Scripts

This directory contains the statistical analysis scripts for the following publication:

O'Malley, A. J., Bobak, C. A., & Barnato, A. E. (2026). Causal inference for intervention
spillover in a stepped wedge cluster-randomized trial: Lessons from a physician network.
*Social Networks*, 86, 431–445. https://doi.org/10.1016/j.socnet.2026.04.017

---

## Overview

Five scripts -- four R Markdown files and one standard R script -- together generate all
statistical results reported in the main text and supplementary material of the article.
Scripts (1)--(4) correspond exactly to the four model specifications reported in the paper.
Script (5) generates selected tables and figures. The first sections of scripts (1)--(4)
follow a common structure; causal estimands based on potential outcomes are implemented
only for the mixed-effect model specifications.

---

## Script Descriptions

### 1. `retention_models_tinv_potout.Rmd`

Estimates mixed-effect logistic regression models and potential outcome-defined causal
estimands (average direct effect, average indirect effect, and related quantities) reported
in the main text. The physician network is fixed at baseline (time-invariant specification).
Regression parameter estimates, standard errors, 95% confidence intervals, and p-values for
both model coefficients and causal estimands are written to separate output files.

### 2. `retention_models_potout.Rmd`

Analogous to (1), but uses the time-varying network specification in which the physician
network is evaluated in the prior time period rather than fixed at baseline.

### 3. `retention_models_tinv_fe.Rmd`

Analogous to (1) in using the time-invariant baseline network, but differs in two respects:
(i) hospital effects are treated as fixed effects rather than random effects, and (ii)
Generalized Estimating Equations (GEE) are used to adjust for clustering by hospital,
yielding population-average (marginal) estimates. Results under both the fixed-effect and
GEE specifications are written to separate output files. Potential outcomes-based causal
estimands are not computed under this specification.

### 4. `retention_models_fe.Rmd`

Analogous to (3) but uses the time-varying network specification in which the physician
network is evaluated in the prior time period.

### 5. `plots4paper.R`

Generates tables and figures involving the key study variables, including visualizations of
the diffusion of direct and indirect intervention effects across study steps and calendar
months. The user must manually toggle between the time-invariant baseline network and the
time-varying network specifications.

---

## Notes

- Scripts (1)--(4) share a common data loading and preprocessing structure in their opening
  sections.
- Potential outcomes-based causal estimands (average direct effect, average indirect effect,
  etc.) are implemented only for the mixed-effect model specifications in scripts (1) and (2).
- All analyses use mixed-effect logistic regression via `lme4`/`glmer` (scripts 1--2) or
  fixed-effect logistic regression and GEE via `geepack` (scripts 3--4).

---

## Dependencies

```r
library(lme4)
library(nlme)
library(geepack)   # scripts 3-4 only
library(Gmisc)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(broom.mixed)
```

---

## Citation

O'Malley, A. J., Bobak, C. A., & Barnato, A. E. (2026). Causal inference for intervention
spillover in a stepped wedge cluster-randomized trial: Lessons from a physician network.
*Social Networks*, 86, 431–445. https://doi.org/10.1016/j.socnet.2026.04.017