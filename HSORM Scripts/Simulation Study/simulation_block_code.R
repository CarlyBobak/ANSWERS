## =========================
## Minimal deps
## =========================
suppressPackageStartupMessages({
  library(Matrix)
  library(lme4)      # preferred
  library(lmtest)    # fallback CRSE
  library(sandwich)  # fallback CRSE
})

## =========================
## Exposure engines (strict prior)
## =========================
# W: row-stochastic dgCMatrix (weights allowed)
# trial_step: named integer vector; value = first treated step (1..5), Inf if never treated
# step: 0..5; strict prior uses neighbors with trial_step < step
phys_exposure <- function(W, trial_step, step,
                          truth = c("proportion","threshold","count","logsum"),
                          strict_prior = TRUE) {
  truth <- match.arg(truth)
  npi <- rownames(W)
  ts  <- trial_step[npi]
  Post <- as.numeric(if (strict_prior) ts < step else ts <= step)
  
  if (truth == "proportion") {
    ex <- as.numeric(W %*% Post)                    # fraction treated (row-stochastic)
  } else {
    B <- as(Matrix::drop0(W) > 0, "dgCMatrix")      # unweighted adjacency
    cnt <- as.numeric(B %*% Post)                   # # treated neighbors
    if (truth == "threshold") ex <- as.numeric(cnt > 0)
    else if (truth == "count") ex <- cnt
    else ex <- log1p(cnt)                           # "logsum"
  }
  names(ex) <- npi
  ex
}

## =========================
## Team sampling utility
## =========================
sample_team <- function(npi_by_hosp, hosp, K = 2, p_cross = 0.05) {
  team <- sample(npi_by_hosp[[hosp]], 1)
  while (length(team) < K) {
    if (runif(1) < p_cross) {
      other <- sample(setdiff(names(npi_by_hosp), hosp), 1)
      team <- c(team, sample(npi_by_hosp[[other]], 1))
    } else {
      team <- c(team, sample(npi_by_hosp[[hosp]], 1))
    }
  }
  unique(team)
}

## =========================
## Safe GLMM fit (retry once)
## =========================
safe_glmer <- function(formula, data) {
  fit <- try(suppressWarnings(
    glmer(formula, data=data, family = binomial(),
          control = glmerControl(optimizer = "bobyqa"))
  ), silent=TRUE)
  ok <- !(inherits(fit, "try-error") || !is.null(fit@optinfo$conv$lme4$messages))
  if (!ok) {
    fit <- try(suppressWarnings(
      glmer(formula, data=data, family = binomial(),
            control = glmerControl(optimizer = "nlminbwrap"))
    ), silent=TRUE)
    ok <- !(inherits(fit, "try-error") || !is.null(fit@optinfo$conv$lme4$messages))
    if (!ok) return(NULL)
  }
  fit
}

## =========================
## Single replicate simulator
## =========================
sim_one <- function(network_panel,
                    N_per_hosp_step = 400,
                    K_phys_per_stay = 2,
                    p_cross_hosp = 0.05,
                    # truth (log-odds)
                    intercept = qlogis(0.22),
                    step_FE = rep(0, 6),        # length 6 for steps 0..5
                    beta1 = log(1.15),          # Post
                    beta2 = log(2.8),           # Expose
                    beta3 = log(0.88),          # Post:Expose
                    tau_hosp = 0.35,
                    # exposure truth / aggregation
                    exposure_truth = c("proportion","threshold","count","logsum"),
                    team_agg = c("sum","mean"),
                    strict_prior = TRUE,
                    use_glmm = TRUE,
                    seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  exposure_truth <- match.arg(exposure_truth)
  team_agg <- match.arg(team_agg)
  
  W_list     <- network_panel$W_list     # list of length 6 (steps 0..5)
  npi        <- network_panel$npi
  trial_step <- network_panel$trial_step # named by NPI
  hospital   <- network_panel$hospital   # named by NPI
  stopifnot(length(W_list) == 6)
  
  # map hospital -> NPIs
  hosp_by_npi <- hospital[npi]
  npi_by_hosp <- split(npi, hosp_by_npi)
  H <- length(npi_by_hosp)
  hosp_levels <- names(npi_by_hosp)
  
  # hospital RE
  b_h <- rnorm(H, 0, tau_hosp); names(b_h) <- hosp_levels
  
  # precompute physician exposures per step according to the TRUTH
  phys_expo <- lapply(0:5, function(s) {
    ex <- phys_exposure(W_list[[s+1]], trial_step, s,
                        truth = exposure_truth, strict_prior = strict_prior)
    names(ex) <- rownames(W_list[[s+1]])
    ex
  })
  
  # simulate patients
  rows <- vector("list", H * 6 * N_per_hosp_step)
  r_id <- 0L
  
  for (s in 0:5) {
    for (h in hosp_levels) {
      # patient-level Post from hospital rollout (min NPI trial step within hospital)
      hosp_ts <- suppressWarnings(min(trial_step[npi_by_hosp[[h]]], na.rm = TRUE))
      hosp_ts <- if (is.infinite(hosp_ts)) 99L else as.integer(hosp_ts)
      Post_patient <- as.integer(hosp_ts <= s)
      
      for (i in seq_len(N_per_hosp_step)) {
        K <- max(1, round(rnorm(1, mean = K_phys_per_stay, sd = 0.3)))
        team <- sample_team(npi_by_hosp, h, K, p_cross_hosp)
        ex_vec <- phys_expo[[s+1]][team]
        Ex_team <- if (team_agg == "sum") sum(ex_vec, na.rm = TRUE) else mean(ex_vec, na.rm = TRUE)
        
        eta <- intercept + step_FE[s+1] + b_h[h] +
          beta1 * Post_patient + beta2 * Ex_team + beta3 * (Post_patient * Ex_team)
        y <- rbinom(1, 1, plogis(eta))
        
        r_id <- r_id + 1L
        rows[[r_id]] <- list(y = y, step = s, hosp = h, Post = Post_patient, Expose = Ex_team)
      }
    }
  }
  
  dat <- do.call(rbind.data.frame, rows)
  dat$hosp  <- factor(dat$hosp, levels = hosp_levels)
  dat$stepF <- factor(dat$step)
  
  if (use_glmm) {
    fit_net   <- safe_glmer(y ~ Post * Expose + stepF + (1 | hosp), data = dat)
    fit_naive <- safe_glmer(y ~ Post + stepF + (1 | hosp), data = dat)
  } else {
    fit_net   <- glm(y ~ Post * Expose + stepF, data = dat, family = binomial())
    fit_naive <- glm(y ~ Post + stepF, data = dat, family = binomial())
  }
  
  list(data = dat, fit_net = fit_net, fit_naive = fit_naive,
       truth = c(beta1 = beta1, beta2 = beta2, beta3 = beta3))
}

## =========================
## Coef/SE extractors
## =========================
coef_se <- function(fit, term, use_glmm = TRUE, cluster = NULL) {
  if (use_glmm) {
    if (is.null(fit)) return(c(NA_real_, NA_real_))
    cf <- suppressWarnings(summary(fit)$coefficients)
    if (!(term %in% rownames(cf))) return(c(NA_real_, NA_real_))
    return(c(est = unname(cf[term, "Estimate"]), se = unname(cf[term, "Std. Error"])))
  } else {
    V <- if (!is.null(cluster)) sandwich::vcovCL(fit, cluster = cluster) else vcov(fit)
    z <- suppressWarnings(coeftest(fit, vcov. = V))
    if (!(term %in% rownames(z))) return(c(NA_real_, NA_real_))
    return(c(est = unname(z[term, 1]), se = unname(z[term, 2])))
  }
}

## =========================
## One setting -> many replicates
## =========================
run_setting <- function(R = 300,
                        network_panel,
                        # pass-through to sim_one
                        exposure_truth = "proportion",
                        team_agg = "sum",
                        beta1 = log(1.15), beta2 = log(2.8), beta3 = log(0.88),
                        intercept = qlogis(0.22),
                        step_FE = rep(0, 6),
                        tau_hosp = 0.35,
                        N_per_hosp_step = 400,
                        K_phys_per_stay = 2,
                        p_cross_hosp = 0.05,
                        strict_prior = TRUE,
                        use_glmm = TRUE,
                        alpha = 0.05,
                        seed = 12345,
                        csv_path = NULL) {
  
  set.seed(seed)
  
  # storage
  terms_net   <- c("Post","Expose","Post:Expose")
  terms_naive <- c("Post (naive)")
  est_net <- se_net <- matrix(NA_real_, nrow = R, ncol = length(terms_net),
                              dimnames = list(NULL, terms_net))
  est_nv  <- se_nv  <- matrix(NA_real_, nrow = R, ncol = length(terms_naive),
                              dimnames = list(NULL, terms_naive))
  fail_net <- fail_nv <- logical(R)
  seeds   <- sample.int(.Machine$integer.max, R)
  
  # helpers
  get_cs <- function(fit, term, use_glmm, cluster_hosp) {
    if (use_glmm) {
      if (is.null(fit)) return(c(NA_real_, NA_real_))
      cf <- suppressWarnings(summary(fit)$coefficients)
      if (!(term %in% rownames(cf))) return(c(NA_real_, NA_real_))
      c(est = unname(cf[term, "Estimate"]), se = unname(cf[term, "Std. Error"]))
    } else {
      V <- sandwich::vcovCL(fit, cluster = cluster_hosp)
      z <- suppressWarnings(lmtest::coeftest(fit, vcov. = V))
      if (!(term %in% rownames(z))) return(c(NA_real_, NA_real_))
      c(est = unname(z[term, 1]), se = unname(z[term, 2]))
    }
  }
  mcse_p <- function(p, n) ifelse(is.na(p), NA_real_, sqrt(p*(1-p)/n))
  
  t0 <- proc.time()["elapsed"]
  
  for (r in seq_len(R)) {
    s1 <- sim_one(network_panel = network_panel,
                  N_per_hosp_step = N_per_hosp_step,
                  K_phys_per_stay = K_phys_per_stay,
                  p_cross_hosp = p_cross_hosp,
                  intercept = intercept,
                  step_FE = step_FE,
                  beta1 = beta1, beta2 = beta2, beta3 = beta3,
                  tau_hosp = tau_hosp,
                  exposure_truth = exposure_truth,
                  team_agg = team_agg,
                  strict_prior = strict_prior,
                  use_glmm = use_glmm,
                  seed = seeds[r])
    
    # NETWORK-ADJUSTED
    if (use_glmm && is.null(s1$fit_net)) {
      fail_net[r] <- TRUE
    } else {
      csP  <- get_cs(s1$fit_net, "Post",        use_glmm, ~ s1$data$hosp)
      csE  <- get_cs(s1$fit_net, "Expose",      use_glmm, ~ s1$data$hosp)
      csPE <- get_cs(s1$fit_net, "Post:Expose", use_glmm, ~ s1$data$hosp)
      est_net[r,] <- c(csP["est"], csE["est"], csPE["est"])
      se_net[r, ] <- c(csP["se"],  csE["se"],  csPE["se"])
    }
    
    # NAÏVE (Post only)
    if (use_glmm && is.null(s1$fit_naive)) {
      fail_nv[r] <- TRUE
    } else {
      csPn <- get_cs(s1$fit_naive, "Post", use_glmm, ~ s1$data$hosp)
      est_nv[r,] <- csPn["est"]; se_nv[r,] <- csPn["se"]
    }
  }
  
  elapsed <- as.numeric(proc.time()["elapsed"] - t0)
  truth <- c(Post = beta1, Expose = beta2, "Post:Expose" = beta3, "Post (naive)" = beta1)
  
  # summaries (network)
  zcrit <- qnorm(1 - alpha/2)
  pack <- function(est_mat, se_mat, term_names) {
    R_eff <- apply(est_mat, 2, function(x) sum(is.finite(x)))
    cover <- (est_mat - zcrit * se_mat <= matrix(truth[term_names], nrow = R, ncol = length(term_names), byrow = TRUE)) &
      (est_mat + zcrit * se_mat >= matrix(truth[term_names], nrow = R, ncol = length(term_names), byrow = TRUE))
    reject <- (abs(est_mat / se_mat) >= zcrit)
    
    out <- lapply(seq_along(term_names), function(j) {
      trm <- term_names[j]
      est <- est_mat[, j]; se <- se_mat[, j]; n <- R_eff[j]
      tval <- truth[trm]
      # stats
      mean_est <- mean(est, na.rm = TRUE)
      bias     <- mean(est - tval, na.rm = TRUE)
      emp_se   <- sd(est, na.rm = TRUE)
      mean_se  <- mean(se, na.rm = TRUE)
      rmse     <- sqrt(mean((est - tval)^2, na.rm = TRUE))
      covg     <- mean(cover[, j], na.rm = TRUE)
      rej      <- mean(reject[, j], na.rm = TRUE)
      type1    <- if (abs(tval) < 1e-12) rej else NA_real_
      power    <- if (abs(tval) >= 1e-12) rej else NA_real_
      # MCSEs for proportions
      mcse_cov <- mcse_p(covg, n)
      mcse_t1  <- if (!is.na(type1)) mcse_p(type1, n) else NA_real_
      mcse_pow <- if (!is.na(power)) mcse_p(power, n) else NA_real_
      # sign error (only when truth ≠ 0)
      sign_err <- if (abs(tval) >= 1e-12) mean(sign(est) != sign(tval), na.rm = TRUE) else NA_real_
      
      data.frame(
        term = trm, truth = tval, mean_est, bias, emp_se, mean_se, rmse,
        coverage = covg, coverage_mcse = mcse_cov,
        type1 = type1, type1_mcse = mcse_t1,
        power = power, power_mcse = mcse_pow,
        sign_error = sign_err,
        n_effective = n,
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, out)
  }
  
  tab_net   <- pack(est_net, se_net, terms_net)
  tab_naive <- pack(est_nv,  se_nv,  terms_naive)
  
  out <- rbind(
    transform(tab_net,  model = "network-adjusted"),
    transform(tab_naive, model = "naive")
  )
  out$conv_fail_rate_net   <- mean(fail_net)
  out$conv_fail_rate_naive <- mean(fail_nv)
  
  attr(out, "elapsed_sec") <- elapsed
  attr(out, "R") <- R
  attr(out, "setting") <- list(exposure_truth = exposure_truth, team_agg = team_agg,
                               beta1 = beta1, beta2 = beta2, beta3 = beta3,
                               N_per_hosp_step = N_per_hosp_step)
  
  if (!is.null(csv_path)) {
    utils::write.csv(out, csv_path, row.names = FALSE)
  }
  out
}

## =========================
## Tiny blocky toy network (Matrix-only)  ---- use locally
## =========================
make_blocky_W_list <- function(H = 6, G_per_h = 25, Tsteps = 6,
                               p_within = 0.12, p_between = 0.02,
                               rollout = c(1,1,2,2,3,3)) {
  N <- H * G_per_h
  npi <- sprintf("N%03d", 1:N)
  hospital <- rep(paste0("H", 1:H), each = G_per_h)
  # base 0/1 adjacency
  A <- Matrix(0, N, N, sparse = TRUE)
  # within
  for (h in 1:H) {
    idx <- ((h-1)*G_per_h + 1):(h*G_per_h)
    for (i in idx) for (j in idx) if (i < j && runif(1) < p_within) { A[i,j] <- 1; A[j,i] <- 1 }
  }
  # between
  for (i in 1:(N-1)) for (j in (i+1):N)
    if (hospital[i] != hospital[j] && runif(1) < p_between) { A[i,j] <- 1; A[j,i] <- 1 }
  diag(A) <- 0
  # row-stochastic
  W <- A; rs <- rowSums(W); iso <- which(rs == 0); if (length(iso)) W[cbind(iso, iso)] <- 1
  W <- W / rowSums(W); rownames(W) <- colnames(W) <- npi
  W_list <- replicate(Tsteps, W, simplify = FALSE)
  # hospital rollout steps; NPIs inherit hospital's step (1..5)
  hr <- setNames(rep_len(rollout, H), paste0("H", 1:H))
  ts <- setNames(as.integer(hr[hospital]), npi)
  list(W_list = W_list, npi = npi, trial_step = ts, hospital = setNames(hospital, npi))
}










