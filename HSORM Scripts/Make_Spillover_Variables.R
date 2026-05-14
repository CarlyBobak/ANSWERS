# Load physician-level spillover effects from CSV
phys_trial_effects_rw <- read.csv(
  file.path(outdir_phys, "physician_level_spillover_effects_rw_strict_prior.csv"),
  row.names = 1
)

# Initialize column to hold patient-level aggregate exposure
data$physWPost <- NA

# Loop over each row (i.e., patient stay) in the dataset
for (i in seq_len(nrow(data))) {
  
  # Identify hospital and trial step for the current patient stay
  hospid <- data$HospitalID[i]
  step <- as.numeric(data$trial_step[i])
  stay <- data$patient_stay_id[i]
  
  # Find all unique NPIs (physicians) encountered by this patient during their stay
  encountered_npi <- unique(
    as.character(trial_encounters[trial_encounters$patient_stay_id == stay, "npi"])
  )
  
  # Extract and sum the spillover values from physicians for this step (step + 1 = column index)
  # NA values are ignored in the summation
  data$physWPost[i] <- sum(
    phys_trial_effects_rw[phys_trial_effects_rw$NPI %in% encountered_npi, step + 1],
    na.rm = TRUE
  )
}

# Save updated dataset with patient-level spillover exposure
write.csv(
  data,
  file.path(outdir_spillover, "acp_data_with_spill_strict_prior_stochastic.csv"),
  row.names = FALSE
)
