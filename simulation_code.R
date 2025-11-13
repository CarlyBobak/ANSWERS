## ---- deps (secure compute friendly)
suppressPackageStartupMessages({
  library(Matrix)
  library(lme4)        # for GLMM fits (preferred)
  library(lmtest)      # for coeftest with CRSE (fallback)
  library(sandwich)    # for cluster-robust vcov (fallback)
})

## ---- exposure (strict prior by default)
compute_exposure <- function(W, trial_step, step, strict_prior = TRUE) {
  npi <- rownames(W)
  ts  <- trial_step[npi]
  Post_phys <- as.numeric(if (strict_prior) ts < step else ts <= step)
  as.numeric(W %*% Post_phys)  # names follow rownames(W)
}

## ---- tiny helper: sample a care team for a patient
sample_team <- function(npi_by_hosp, hosp, K = 2, p_cross = 0.05) {
  # choose 1 from focal hospital, and with p_cross add 1 from other hospitals
  team <- sample(npi_by_hosp[[hosp]], size = 1)
  while (length(team) < K) {
    if (runif(1) < p_cross) {
      other_h <- sample(setdiff(names(npi_by_hosp), hosp), 1)
      team <- c(team, sample(npi_by_hosp[[other_h]], 1))
    } else {
      team <- c(team, sample(npi_by_hosp[[hosp]], 1))
    }
  }
  unique(team)
}

## ---- single replicate
sim_one <- function(network_panel,
                    N_per_hosp_step = 400,   # patients per hospital per step
                    K_phys_per_stay = 2,     # team size (avg)
                    p_cross_hosp = 0.05,     # small cross-hospital care mixing
                    # truth parameters (log-odds scale)
                    intercept = qlogis(0.22),
                    step_FE = rep(0, 6),     # optional fixed effects for steps 0..5
                    beta1 = log(1.15),       # direct effect (Post)
                    beta2 = log(2.8),        # spillover (Expose)
                    beta3 = log(0.88),       # interaction Post*Expose
                    tau_hosp = 0.35,         # SD of hospital RE
                    strict_prior = TRUE,
                    use_glmm = TRUE,
                    seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  W_list     <- network_panel$W_list              # list of dgCMatrix, steps 0..5
  npi        <- network_panel$npi
  trial_step <- network_panel$trial_step          # named by NPI (first treated step)
  hospital   <- network_panel$hospital            # named by NPI
  stopifnot(length(W_list) == 6)                  # steps 0..5
  
  # map hospital -> NPIs
  hosp_by_npi <- hospital[npi]
  npi_by_hosp <- split(npi, hosp_by_npi)
  H <- length(npi_by_hosp)
  hosp_levels <- names(npi_by_hosp)
  
  # hospital random intercepts
  b_h <- rnorm(H, 0, tau_hosp); names(b_h) <- hosp_levels
  
  # precompute physician exposures per step (strict prior)
  phys_expo <- lapply(0:5, function(s) {
    W <- W_list[[s + 1]]
    ex <- compute_exposure(W, trial_step, step = s, strict_prior = strict_prior)
    names(ex) <- rownames(W)
    ex
  })
  
  # simulate patients across steps and hospitals
  rows <- vector("list", H * 6 * N_per_hosp_step)
  r_id <- 0L
  
  for (s in 0:5) {
    for (h in hosp_levels) {
      # patient-level Post: cluster-level rollout (hospital turns on at min trial_step among NPIs in that hospital)
      hosp_ts <- suppressWarnings(min(trial_step[npi_by_hosp[[h]]], na.rm = TRUE))
      hosp_ts <- if (is.infinite(hosp_ts)) 99L else as.integer(hosp_ts)
      Post_patient <- as.integer(hosp_ts <= s)  # treated from its step onward
      
      for (i in seq_len(N_per_hosp_step)) {
        # sample team
        K <- max(1, round(rnorm(1, mean = K_phys_per_stay, sd = 0.3)))
        team <- sample_team(npi_by_hosp, h, K, p_cross_hosp)
        # patient exposure = sum of team exposures (you can switch to mean if desired)
        Ex_team <- sum(phys_expo[[s + 1]][team], na.rm = TRUE)
        
        eta <- intercept +
          step_FE[s + 1] +
          b_h[h] +
          beta1 * Post_patient +
          beta2 * Ex_team +
          beta3 * (Post_patient * Ex_team)
        
        p <- plogis(eta)
        y <- rbinom(1, 1, p)
        
        r_id <- r_id + 1L
        rows[[r_id]] <- list(
          y = y,
          step = s,
          hosp = h,
          Post = Post_patient,
          Expose = Ex_team
        )
      }
    }
  }
  
  dat <- do.call(rbind.data.frame, rows)
  dat$hosp <- factor(dat$hosp, levels = hosp_levels)
  # Step FE as factors (0..5) to mimic your analysis
  dat$stepF <- factor(dat$step)
  
  # Fit models
  fit_naive <- NULL
  fit_net   <- NULL
  if (use_glmm) {
    # GLMMs with hospital RE and step FEs
    fit_naive <- suppressWarnings(
      glmer(y ~ Post + stepF + (1 | hosp),
            data = dat, family = binomial(), control = glmerControl(optimizer = "bobyqa"))
    )
    fit_net <- suppressWarnings(
      glmer(y ~ Post * Expose + stepF + (1 | hosp),
            data = dat, family = binomial(), control = glmerControl(optimizer = "bobyqa"))
    )
  } else {
    # GLM with cluster-robust SE by hospital (fallback)
    fit_naive_glm <- glm(y ~ Post + stepF, data = dat, family = binomial())
    fit_net_glm   <- glm(y ~ Post * Expose + stepF, data = dat, family = binomial())
    V_naive <- vcovCL(fit_naive_glm, cluster = ~ hosp)
    V_net   <- vcovCL(fit_net_glm,   cluster = ~ hosp)
    fit_naive <- list(glm = fit_naive_glm, V = V_naive)
    fit_net   <- list(glm = fit_net_glm,   V = V_net)
  }
  
  list(data = dat, fit_naive = fit_naive, fit_net = fit_net,
       truth = c(beta1 = beta1, beta2 = beta2, beta3 = beta3))
}

## ---- extract estimates (+ CRSE support)
get_coef <- function(fit, term, use_glmm = TRUE) {
  if (use_glmm) {
    cf <- suppressWarnings(summary(fit)$coefficients)
    if (!(term %in% rownames(cf))) return(c(NA, NA))
    c(est = unname(cf[term, "Estimate"]), se = unname(cf[term, "Std. Error"]))
  } else {
    # fit is list(glm, V); use coeftest to pull robust SE
    z <- suppressWarnings(coeftest(fit$glm, vcov. = fit$V))
    if (!(term %in% rownames(z))) return(c(NA, NA))
    c(est = unname(z[term, 1]), se = unname(z[term, 2]))
  }
}

## ---- replicate wrapper with bias/coverage/power
sim_many <- function(R = 500, network_panel, ...,
                     use_glmm = TRUE,
                     target = c("Post", "Expose", "Post:Expose"),
                     null = c(0, 0, 0),
                     alpha = 0.05, seed = 123) {
  set.seed(seed)
  target <- match.arg(target, several.ok = TRUE)
  out <- vector("list", R)
  
  for (r in seq_len(R)) {
    s1 <- sim_one(network_panel = network_panel, use_glmm = use_glmm, seed = sample.int(1e8, 1), ...)
    fits <- list(naive = s1$fit_naive, net = s1$fit_net)
    
    # pull from network-adjusted model
    ests <- list()
    for (trm in target) {
      term_name <- switch(trm,
                          "Post" = "Post",
                          "Expose" = "Expose",
                          "Post:Expose" = "Post:Expose")
      ce <- get_coef(fits$net, term_name, use_glmm = use_glmm)
      ests[[trm]] <- ce
    }
    out[[r]] <- ests
  }
  
  # aggregate
  agg <- lapply(target, function(trm) {
    est <- sapply(out, function(e) e[[trm]]["est"])
    se  <- sapply(out, function(e) e[[trm]]["se"])
    tru <- switch(trm,
                  "Post" = 0 + 0,      # we’ll compute bias relative to truth injected via sim_one(...)
                  "Expose" = 0 + 0,
                  "Post:Expose" = 0 + 0)
    # we can’t know the per-run truth here—so return summary of estimates and SEs;
    # in practice, run per "setting" (fixed truths), compute stats below:
    c( mean_est = mean(est, na.rm = TRUE),
       mean_se  = mean(se,  na.rm = TRUE),
       sd_est   = sd(est,   na.rm = TRUE),
       prop_na  = mean(!is.finite(est)) )
  })
  names(agg) <- target
  do.call(rbind, agg)
}
