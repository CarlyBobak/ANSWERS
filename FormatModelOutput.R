variable_mapping <- list(
  "(Intercept)" = "Baseline",
  "gameTime" = "Game Time",
  "chosp_step2" = "Step of the trial (Two)",
  "chosp_step3" = "Step of the trial (Three)",
  "chosp_step4" = "Step of the trial (Four)",
  "chosp_step5" = "Step of the trial (Five)",
  "physWPost" = "Physician Spillover",
  "ctrialmonth2" = "Month of the Trial (Two)",
  "ctrialmonth3" = "Month of the Trial (Three)",
  "ctrialmonth4" = "Month of the Trial (Four)",
  "ctrialmonth5" = "Month of the Trial (Five)",
  "ctrialmonth6" = "Month of the Trial (Six)",
  "ctrialmonth7" = "Month of the Trial (Seven)",
  "ctrialmonth8" = "Month of the Trial (Eight)",
  "ctrialmonth9" = "Month of the Trial (Nine)",
  "ctrialmonth10" = "Month of the Trial (Ten)",
  "ctrialmonth11" = "Month of the Trial (Eleven)",
  "ctrialday" = "Trial Day",
  "cbregionsouth" = "Region (South)",
  "cbregionwest" = "Region (West)",
  "cageg75" = "Patient age (75–84 years)",
  "cageg85" = "Patient age (≥ 85 years)",
  "ccovid_statusYes" = "Covid status (Covid positive)",
  "cad_cd" = "Coronary artery disease",
  "dem_cd" = "Dementia",
  "heart_cd" = "Heart failure",
  "diab_cd" = "Diabetes",
  "lung_cd" = "COPD",
  "ccharl_sc1" = "Charlson score (One)",
  "ccharl_sc2" = "Charlson score (Two)",
  "ccharl_sc3" = "Charlson score (Three)",
  "ccharl_sc4" = "Charlson score (Four)",
  "ccharl_sc5" = "Charlson score (Five)",
  "n_w2mds" = "Practice size",
  "covid_stays" = "Proportion of hospitalizations with COVID",
  "Q2_Y20_Adj_ACP" = "ACP rate in Quarter 2 of 2020 (effect of a 10 percentage-point increase)",
  "Delta_ACP" = "Change in ACP rate between Q2 of 2019 and Q1 of 2020",
  "csqY" = "Would you be surprised if the patient died in the next year (Yes)",
  "phase" = "Intervention"
)


format_model_output_logistic <- function(model, variable_mapping) {
  # Exponentiate estimates for odds ratios
  tidy_model <- broom::tidy(model) %>%
    mutate(estimate = exp(estimate))
  
  # Dynamically replace term components based on mapping
  tidy_model$term <- sapply(tidy_model$term, function(term) {
    components <- unlist(strsplit(term, ":"))
    mapped_components <- sapply(components, function(component) {
      if (component %in% names(variable_mapping)) {
        return(variable_mapping[[component]])
      } else {
        return(component)
      }
    })
    return(paste(mapped_components, collapse = " x "))
  })
  
  # Calculate and exponentiate 95% CI for odds ratios
  ci <- exp(confint(model, method="Wald", level=0.95))
  ci_df <- as.data.frame(ci, row.names = NULL)
  names(ci_df) <- c("lower_ci", "upper_ci")
  
  # Merge tidy model output with CIs
  tidy_model <- merge(tidy_model, ci_df, by = "row.names", all.x = TRUE)
  names(tidy_model)[1] <- "term"
  
  # Include formatted output with Odds Ratios, CIs, and p-values
  tidy_model <- tidy_model %>%
    mutate(formatted_output = paste0("Odds Ratio: ", round(estimate, 3),
                                     ", 95% CI: [", round(lower_ci, 3), ", ", round(upper_ci, 3), "]",
                                     ", p-value: ", format.pval(p.value, digits = 3, eps = 0.001)))
  
  # Define priority variables and their interactions
  priority_vars <- c("Intervention", "Game Time", "Physician Spillover", "Step of the trial", "Proportion of hospitalizations with COVID")
  
  # Sort dataframe by priority, placing priority variables and their interactions at the top
  tidy_model <- tidy_model %>%
    arrange(!Term %in% priority_vars, Term)
  
  # Select and rename columns for publication
  tidy_model <- tidy_model %>%
    select(Term = term, Odds_Ratio = estimate, Lower_CI = lower_ci, Upper_CI = upper_ci, P_Value = p.value, Formatted = formatted_output)
  
  return(tidy_model)
}
