# ------------------------------------------------------------------------------
# Script: AnalysisFlaker.R
#
# Purpose:
#   Implement four single trial estimation methods across all voxels of a single
#   run from a Flaker experiment to obtain trial-wise activation estimates.
#   Methods:
#   1) Least squares-All (LS-A)
#   2) Least squares-Separate (LS-s)
#   3) Extended Least squares-Separate (ELS-S)
#   4) Finite BOLD Response-Separate (FS)
#   Additionally, the voxel showing the maximum z-statistic in the COPE 5 contrast
#   (incongruent > congruent) is identified, and the full model outputs are saved
#   separately for further visualization.
#
# Outputs:
#   - results-run1.RData: List containing voxel coordinates and trial-wise 
#     beta estimates for all non-empty voxels
#   - output_max-run1.RData: Outputs from lm objects for the voxel with maximum 
#     activation
#
# Additional notes:
#   - Confound regressors are included for LS-A, LS-S, and ELS-S methods but 
#     excluded for FS method to avoid multicollinearity
#   - The zstat5_anatomical.nii.gz file contains the z-statistics of COPE 5 but
#     in the subject/individual space (not standard space)
#
# ------------------------------------------------------------------------------

#Load packages
library(fmristat)
library(niftiR6)

sub <- 'sub-01'
run <- 1
path <- './data/'
cv <- read.table(paste0(path,sub,'/run-',run,'.feat/confoundevs.txt')) #Confound regressors
n_trial <- 24

#Find voxel of maximum activation based on z-stat
c5 <- niftiR6::readNifti(paste0(path,sub,'/run-',run,'.feat/stats/zstat5_anatomical.nii.gz'))
max_zstat <- max(c5[]) #Get maximum z-stat
maxloc <- which(c5[]== max_zstat,arr.ind=T) #Find location of max z-stat

#Create fMRI design for each method
#LS-A
tims_LA <- fmristat::read_FSL_featdir(paste0(path,sub,'/run-',run,'.feat/'))
tims_LA$trial_model <- 'LA+S'
tims_LA$make_design_matrix()
tims_LA$set_confound_matrix(cv[,1:6])
#LS-S
tims_LS <- fmristat::read_FSL_featdir(paste0(path,sub,'/run-',run,'.feat/'))
tims_LS$trial_model <- 'LS+S'
tims_LS$make_design_matrix()
tims_LS$set_confound_matrix(cv[,1:6])
#ELS-S
tims_ELS <- fmristat::read_FSL_featdir(paste0(path,sub,'/run-',run,'.feat/'))
tims_ELS$trial_model <- 'LS+T+D+S'
tims_ELS$make_design_matrix()
tims_ELS$set_confound_matrix(cv[,1:6])
#FS
tims_FS <- fmristat::read_FSL_featdir(paste0(path,sub,'/run-',run,'.feat/'))
tims_FS$trial_model <- 'FS+S'
tims_FS$make_design_matrix()

#get data
data <- niftiR6::readNifti(tims_LA$get_link()) #Same for all methods

#Store all estimates
results <- list()
output_max <- list()
idx <- 1

#Loop over all voxels
for (x in 1:data$dims[2]){
  for (y in 1:data$dims[3]){
    for (z in 1:data$dims[4]){
      if (var(data[x,y,z,]) == 0){ #Skip empty voxels
        next
        } else {
          cat('Non-empty voxel:', idx, '\n')
        
          signal <- data[x,y,z,] #Extract BOLD signal
          
          #Apply methods
          stmat <- matrix(NA, nrow=n_trial, ncol = 4) #ncol = number of methods
          colnames(stmat) <- c('LS-A', 'LS-S', 'ELS-S', 'FS')
          
          out_LA <- vector("list", 1) 
          out_LS <- vector("list", length(tims_LS$get_design_matrix()))
          out_ELS <- vector("list", length(tims_ELS$get_design_matrix()))
          out_FS <- vector("list", length(tims_FS$get_design_matrix()))
          
          #LS-A
          X_LA <- tims_LA$get_design_matrix()[[1]]
          out_LA <- lm(signal ~ 1 + X_LA)
          stmat[,1] <- coef(out_LA)[2:(n_trial+1)]
          #LS-S 
          for(i in 1:length(tims_LS$get_design_matrix())) {
            X_LS <- tims_LS$get_design_matrix()[[i]]
            out_LS[[i]] <- lm(signal ~ 1 + X_LS)
            stmat[i,2] <- coef(out_LS[[i]])[2]
          }
          #ELS-S
          for(j in 1:length(tims_ELS$get_design_matrix())) {
            X_ELS <- tims_ELS$get_design_matrix()[[j]]
            out_ELS[[j]] <- lm(signal ~ 1 + X_ELS)
            stmat[j,3] <- coef(out_ELS[[j]])[2]
          }
          #FS
          for(k in 1:length(tims_FS$get_design_matrix())) {
            X_FS <- tims_FS$get_design_matrix()[[k]]
            out_FS[[k]] <- lm(signal ~ 1 + X_FS)
            stmat[k,4] <- max(coef(out_FS[[k]])[2:9], na.rm = TRUE)
          }
          
          #Store results
          results[[idx]] <- list(x = x, y = y, z = z, betas = stmat)
          idx <- idx + 1
          
          #Store full model outputs for the voxel with maximum z-statistic
          if ((maxloc[1] == x) & (maxloc[2] == y) & (maxloc[3] == z)){
            output_max <- list(signal = signal,
                               x = x, 
                               y = y, 
                               z = z,
                               output_LA = out_LA, 
                               output_LS = out_LS,
                               output_ELS = out_ELS,
                               output_FS = out_FS)
          }
        }
    }
  }
}

#Save all estimation results
save(results, file = './Results/results-run1.RData')

#Save outputs from the selected voxel
save(output_max, file = './Results/output_max-run1.RData')







