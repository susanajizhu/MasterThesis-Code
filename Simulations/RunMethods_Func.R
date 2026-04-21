# ------------------------------------------------------------------------------
# Function: run_methods
#
# Purpose:
#   Run single-trial estimation methods for fMRI data:
#   1) Least Squares-All (LS-A)
#   2) Least Squares-Separate (LS-S) 
#   3) Extended Least Squares-Separate (ELS-S)
#   4) Finite Bold Response-Separate (FS)
#   The function fits general linear models (GLM) for each method and
#   extracts single-trial beta estimates
#
# Inputs:
#   n_trial         : Number of trials
#   functional_data : Filename of the 4D functional NIfTI image 
#                     (default: 'simulated_filtered_func.nii.gz')
#
# Output:
#   A list containing:
#     1) stmat      : Matrix of single-trial beta estimates
#                     (rows = trials, columns = methods)
#     2) out_LA     : output from lm object from fitting LS-A model
#     3) out_LS     : List of outputs from lm objects from fitting LS-S model
#     4) out_ELS    : List of outputs from lm objects from fitting ELS-S model
#     5) out_FS     : List of outputs from lm objects from fitting FS model
#
# Additional notes:
#     1) For the ELS-S method, the beta coefficient associated with the
#        canonical HRF is selected for each trial as the activation estimate
#     2) For the FS method, the maximum beta value across stimulus regressors
#        is selected for each trial as the activation estimate
#
# ------------------------------------------------------------------------------

#Packages
library(fmristat)
library(niftiR6)

run_methods <- function(n_trial, functional_data = 'simulated_filtered_func.nii.gz'){
  #Read data and create fMRI designs
  tims_LA <- fmristat::read_FSL_featdir('./datasets/analysis/', functional_data = functional_data)
  tims_LS <- fmristat::read_FSL_featdir('./datasets/analysis/', functional_data = functional_data)
  tims_ELS <- fmristat::read_FSL_featdir('./datasets/analysis/', functional_data = functional_data)
  tims_FS <- fmristat::read_FSL_featdir('./datasets/analysis/', functional_data = functional_data)
  
  data <- niftiR6::readNifti(tims_LA$get_link()) #Same for all methods
  y <- data[1,1,1,] #Extract BOLD signal
  
  #Create design matrices
  tims_LA$trial_model <- 'LA+S'
  tims_LA$make_design_matrix()
  
  tims_LS$trial_model <- 'LS+S'
  tims_LS$make_design_matrix()
  
  tims_ELS$trial_model <- 'LS+T+D+S'
  tims_ELS$make_design_matrix()
  
  tims_FS$trial_model <- 'FS+S'
  tims_FS$make_design_matrix()
  
  #Store results
  stmat <- matrix(NA, nrow=n_trial, ncol = 4) #ncol = number of methods
  colnames(stmat) <- c('LS-A', 'LS-S', 'ELS-S', 'FS')
  
  out_LA <- vector("list", 1) 
  out_LS <- vector("list", length(tims_LS$get_design_matrix()))
  out_ELS <- vector("list", length(tims_ELS$get_design_matrix()))
  out_FS <- vector("list", length(tims_FS$get_design_matrix()))
  
  #Run LS-A Method
  X_LA <- tims_LA$get_design_matrix()[[1]]
  out_LA <- lm(y ~ 1 + X_LA)
  stmat[,1] <- coef(out_LA)[2:(n_trial+1)]
  #Run LS-S Method
  for(i in 1:length(tims_LS$get_design_matrix())) {
    X_LS <- tims_LS$get_design_matrix()[[i]]
    out_LS[[i]] <- lm(y ~ 1 + X_LS)
    stmat[i,2] <- coef(out_LS[[i]])[2]
  }
  #Run ELS-S Method
  for(j in 1:length(tims_ELS$get_design_matrix())) {
    X_ELS <- tims_ELS$get_design_matrix()[[j]]
    out_ELS[[j]] <- lm(y ~ 1 + X_ELS)
    stmat[j,3] <- coef(out_ELS[[j]])[2]
    #out_coef <- coef(out_ELS[[j]])
    #stmat[j, 3] <- sign(out_coef[2])*sqrt(out_coef[2]^2 + out_coef[3]^2 + out_coef[4]^2)
  }
  #Run FS Method
  for(k in 1:length(tims_FS$get_design_matrix())) {
    X_FS <- tims_FS$get_design_matrix()[[k]]
    out_FS[[k]] <- lm(y ~ 1 + X_FS)
    stmat[k,4] <- max(coef(out_FS[[k]])[2:9], na.rm = TRUE)
  }
  
  return(list(stmat = stmat,
              out_LA = out_LA,
              out_LS = out_LS, 
              out_ELS = out_ELS,
              out_FS = out_FS))
}