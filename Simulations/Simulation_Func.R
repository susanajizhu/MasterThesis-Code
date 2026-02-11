# ------------------------------------------------------------------------------
# Function: simulation
#
# Purpose:
#   Run a complete fMRI simulation study comparing single-trial estimation
#   methods under varying conditions of ISI, amplitude, latency and duration
#   For each iteration, the function:
#     1) Generates event timings and simulated BOLD signal
#     2) Runs LS-A, LS-S, ELS-S, and FS methods
#     3) Computes RMSE and Pearson's correlation coefficients of beta estimates
#     5) Saves all relevant outputs to RData files
#
# Inputs:
#   N      : Number of iterations per simulation
#   isi_a  : Lower bound of ISI sampled from U[isi_a, isi_b] (seconds)
#   isi_b  : Upper bound of ISI sampled from U[isi_a, isi_b] (seconds)
#   dur_a  : Lower bound of duration sampled from U[dur_a, dur_b] (seconds)
#   dur_b  : Upper bound of duration sampled from U[dur_a, dur_b] (seconds)
#   lat_a  : Lower bound of latency sampled from U[lat_a, lat_b] (seconds)
#   lat_b  : Upper bound of latency sampled from U[lat_a, lat_b] (seconds)
#   amp_a  : Lower bound of amplitude sampled from U[amp_a, amp_b]
#   amp_b  : Upper bound of amplitude sampled from U[amp_a, amp_b]
#
# Output:
#   A list containing:
#     1) y.all        : List of simulated BOLD time series
#     2) beta.all     : List of true beta values
#     3) stmat.all    : Array of estimated beta values
#                      (trials × methods × iterations)
#     4) rmse         : Matrix of RMSE values (iterations × methods)
#     5) corr         : Matrix of correlation values (iterations × methods)
#     6) out_LA.all   : List of outputs from lm objects from LS-A models
#     7) out_LS.all   : List of outputs from lm objects from LS-S models
#     8) out_ExtLS.all: List of outputs from lm objects from ELS-S models
#     9) out_FS.all   : List of outputs from lm objects from FS models   
#
# Additional notes:
#   - EV files are generated using create_ev()
#   - Methods are implemented using run_methods()
#   - A template FEAT directory exists at './real-data/sub-01/run-1.feat/'
#   - A simulated functional data is written to './datasets/analysis/'
#
# ------------------------------------------------------------------------------

#Load packages
library(fmristat)
library(niftiR6)

simulation <- function(N, isi_a, isi_b, dur_a, dur_b, lat_a, lat_b, amp_a, amp_b){
  #Set fixed parameters
  TR = 2
  n_trial = 24 #Number of trials
  signal_noise_ratio = 5 #Signal to noise ratio (SNR)
  sampling_rate = 20
  
  #Store RMSE
  rmse <- matrix(NA, nrow = N, ncol = 4) #ncol = number of methods
  colnames(rmse) <- c('LS-A', 'LS-S', 'ELS-S', 'FS')
  rownames(rmse) <- paste0('Iter ', 1:N)
  
  #Store Pearson's correlation coefficients
  corr <- matrix(NA, nrow = N, ncol = 4) #ncol = number of methods
  colnames(corr) <- c('LS-A', 'LS-S', 'ELS-S', 'FS')
  rownames(corr) <- paste0('Iter ', 1:N)
  
  #Store beta estimates
  stmat.all <- array(NA, dim = c(n_trial, 4, N)) #Second dimension is the number of methods
  dimnames(stmat.all) <- list(1:n_trial, c('LS-A', 'LS-S', 'ELS-S', 'FS'), paste0('Iteration: ', 1:N))
  
  #Store output results from lm
  out_LA.all <- vector('list', N)
  out_LS.all <- vector('list', N)
  out_ExtLS.all <- vector('list', N)
  out_FS.all <- vector('list', N)
  
  #Store BOLD signal
  y.all <- vector('list', N)
  
  #Store true betas
  beta.all <- vector('list', N)
  
  #Set seed for each iteration for reproducibility
  seeds <- 1:N
  
  # ========================= Simulation ======================================
  start.time <- Sys.time()
  
  for (s in 1:N){
    set.seed(seeds[s])
    cat('Iteration: ', s,'\n')
    
    # ====================================================
    #             Part 1: Generate data         
    # ====================================================
    
    #Set conditions
    isi <- runif(n_trial, min = isi_a, max = isi_b)
    duration <- runif(n_trial, min = dur_a, max = dur_b)
    latency <- runif(n_trial, min = lat_a, max = lat_b)
    amp <- runif(n_trial, min = amp_a, max = amp_b)
    
    #Set parameters
    n = round(sum(isi),0) + 20 #total stimulus time + offset of 20 seconds
    num_volumes = n/TR
    
    #Sample trial type
    t = rep(c(1,2), each = n_trial/2) 
    t = sample(t)
    
    #Create EV files with create_ev()
    ev_files <- create_ev(t = t,
                          n_trial = n_trial,
                          isi = isi,
                          amp = amp,
                          mis_dt = duration,
                          mis_lat = latency)
    
    ev_data <- ev_files$ev_data
    
    #Create new fmri design
    fmri_des <- fmristat::fmri_design$new()
    fmri_des$set_time_vitals(Hz = sampling_rate, TR = TR, nvol = num_volumes)
    fmri_des$trial_model <- 'SD+S' #Standard regressor
    
    #Event timings
    EV <- "./datasets/data/custom_timing_files/ev.txt"
    timings <- fmristat::ev_timings$new()
    timings$set_link(EV)
    fmri_des$ev_timings[[1]] <- timings$read_FSL_timings(EV)
    
    #Create design matrix
    fmri_des$make_design_matrix()
    
    X <- fmri_des$get_design_matrix()[[1]]
    e <- rnorm(num_volumes, mean = 0, sd = mean(X[X>0])/signal_noise_ratio)
    y <- X + e
    y.all[[s]] <- y #Store BOLD signal from each iteration
    
    #Store time series data in a new clone file
    path = './real-data/sub-01/run-1.feat/'
    tims <- fmristat::read_FSL_featdir(path)
    dat <- niftiR6::readNifti(tims$get_link())
    
    dat2 = dat$clone()
    dat2$setData(array(y, dim = c(1,1,1,num_volumes)))
    dat2$changeName('simulated_filtered_func')
    dat2$writeData()
    file.rename(file.path(path,'simulated_filtered_func.nii.gz'),
                file.path('./datasets/analysis/simulated_filtered_func.nii.gz'))
    
    # ====================================================
    #             Part 2: Run methods         
    # ====================================================
    
    cat('Running methods...','\n')
    
    #Methods: LS-A, LS-S, ELS-S, FS
    results <- run_methods(n_trial = n_trial)
    stmat <- results$stmat
    out_LA <- results$out_LA
    out_LS <- results$out_LS
    out_ExtLS <- results$out_ExtLS
    out_FS <- results$out_FS
    
    cat('Beta estimates from Iteration', s, ':','\n')
    print(stmat)
    v.beta <- c(subset(ev_data, trial_type == 1)$weight, subset(ev_data, trial_type == 2)$weight)
    beta.true <- matrix(rep(v.beta, 4), nrow = n_trial, ncol = 4)
    
    beta.all[[s]] <- v.beta #Store true beta values from each iteration
    
    stmat.all[,,s] <- stmat #Store estimates from each iteration
    out_LA.all[[s]] <- out_LA #Store output results from each iteration
    out_LS.all[[s]] <- out_LS
    out_ExtLS.all[[s]] <- out_ExtLS
    out_FS.all[[s]] <- out_FS
    
    #Calculate and store RMSE
    for (col in 1:ncol(stmat)){
      rmse[s, col] <- round(sqrt(mean((stmat[,col] - beta.true[,col])^2)),3)
    }
    
    #Calculate and store Pearson's correlation coefficients
    if (sd(amp) != 0){ #Calculate correlation only if amplitude is not constant
      for (col in 1:ncol(stmat)){
        corr[s,col] = round(cor(stmat[,col], beta.true[,col]),3) 
      }
    }
    
    #End of iteration
    cat('Iteration ', s, ' finished','\n')
  }
  
  #Save information in RData file
  cat('Saving data in RData file...', '\n')
  
  if (!dir.exists("./Results")) {
    dir.create("./Results", recursive = TRUE)
  }
  
  save(y.all,
       beta.all,
       rmse,
       corr,
       stmat.all, 
       out_LA.all, 
       out_LS.all, 
       out_ExtLS.all,
       out_FS.all,
       file = (paste0('./Results/','isi', isi_a, '_', isi_b, '-', 
                      'dur', dur_a, '_', dur_b, '-',
                      'lat', lat_a, '_', lat_b, '-',
                      'amp', amp_a, '_', amp_b,'.RData')))
  
  end.time <- Sys.time()
  print(end.time - start.time)

  return(list(y.all = y.all,
              beta.all = beta.all,
              stmat.all = stmat.all, 
              rmse = rmse,
              corr = corr,
              out_LA.all = out_LA.all, 
              out_LS.all = out_LS.all,
              out_ExtLS.all = out_ExtLS.all, 
              out_FS.all = out_FS.all))
}


  









