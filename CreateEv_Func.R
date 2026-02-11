# ------------------------------------------------------------------------------
# Function: create_ev
#
# Purpose:
#   Generate event timing (EV) files for fMRI simulations and store them locally
#   Two sets of EV files are created:
#     1) Data EVs: used to generate BOLD signal (with misspecifications)
#     2) Analysis EVs: used for model fitting
#   
# Inputs:
#   t        : Vector of trial types (e.g., 1 or 2)
#   n_trial  : Total number of trials
#   isi      : Vector of inter-stimulus intervals (seconds)
#   amp      : Signal amplitude (scalar or vector)
#   mis_dt   : Duration misspecification (seconds; scalar or vector)
#   mis_lat  : Latency misspecification (seconds; scalar or vector)
#   
# Output:
#   A list containing:
#     1) ev_data    : data frame with true event timings for signal generation
#     2) ev_analysis: data frame with assumed event timings for model fitting
#
# Intermediate output:
#   Writes EV files to "./datasets/analysis/custom_timing_files/" and 
#   "./datasets/data/custom_timing_files/"
#
# Additional notes:
#   1) For model fitting (Analysis EVs), the assumed amplitude is 1,
#      the assumed duration is 2, and the assumed latency is 0.
#      When no misspecifications in duration and latency are introduced in the 
#      data, mis_dt is 2 and mis_lat is 0 (default setting).
#   2) EV files follow the standard format in FSL. Each event timing file
#      contains 3 columns: onset timing, duration, and weight.
# 
# ------------------------------------------------------------------------------

create_ev <- function(t, n_trial, isi,  amp = 1, mis_dt = 2, mis_lat = 0){
  #Create directories to store EVs
  if (!dir.exists("./datasets")) {
    dir.create("./datasets", recursive = TRUE)
  }
  if (!dir.exists("./datasets/analysis")) {
    dir.create("./datasets/analysis", recursive = TRUE)
  }
  if (!dir.exists("./datasets/data")) {
    dir.create("./datasets/data", recursive = TRUE)
  }
  
  #Create Analysis EVs
  ev_analysis <- data.frame(timing = rep(0, n_trial), duration = 2, weight = 1, trial_type = t)
  ev_analysis$timing[1] <- round(isi[1],0)
  for (j in 2:n_trial){
    #onset timing is defined as the cumulative ISI
    ev_analysis$timing[j] <- round(ev_analysis$timing[j-1] + isi[j],0)
  }
  #Store Analysis EVs locally
  if (!dir.exists("./datasets/analysis/custom_timing_files")) {
    dir.create("./datasets/analysis/custom_timing_files", recursive = TRUE)
  }
  write.table(subset(ev_analysis, trial_type == 1)[, c("timing", "duration", "weight")], file = "./datasets/analysis/custom_timing_files/ev1.txt", sep = "\t", col.names = FALSE, row.names = FALSE)
  write.table(subset(ev_analysis, trial_type == 2)[, c("timing", "duration", "weight")], file = "./datasets/analysis/custom_timing_files/ev2.txt", sep = "\t", col.names = FALSE, row.names = FALSE)
  
  #Create Data EVs
  ev_data <- ev_analysis
  ev_data$timing <- ev_data$timing + round(mis_lat,0)
  ev_data$duration <- round(mis_dt,0)
  ev_data$weight <- amp
  #Store Data EVs locally
  if (!dir.exists("./datasets/data/custom_timing_files")) {
    dir.create("./datasets/data/custom_timing_files", recursive = TRUE)
  }
  write.table(subset(ev_data, trial_type == 1)[, c("timing", "duration", "weight")], file = "./datasets/data/custom_timing_files/ev1.txt", sep = "\t", col.names = FALSE, row.names = FALSE)
  write.table(subset(ev_data, trial_type == 2)[, c("timing", "duration", "weight")], file = "./datasets/data/custom_timing_files/ev2.txt", sep = "\t", col.names = FALSE, row.names = FALSE)
  write.table(ev_data[, c("timing", "duration", "weight")], file = "./datasets/data/custom_timing_files/ev.txt", sep = "\t", col.names = FALSE, row.names = FALSE)
  
  return(list(ev_data = ev_data, ev_analysis = ev_analysis))
}






