# ------------------------------------------------------------------------------
# Script: PlotResults.R
#
# Purpose:
#   Visualize results from the fMRI simulations with boxplots.
#   For each amplitude condition, the script:
#     1) Loads simulation results from the './Results/' directory
#     2) Extracts simulation conditions from filenames
#     3) Extracts RMSE and correlation values
#     4) Formats results into dataframes with 4 columns: ISI, latency, duration,
#        and rmse/correlation value
#     5) Generates and saves boxplots comparing estimation methods in the
#        './Results/' directory
#
# Amplitude conditions:
#   amp1_1: Fixed amplitude at 1
#   amp2_4: Amplitude sampled from U[2,4]
#   amp4_6: Amplitude sampled from U[4,6]
#
# Outputs:
#   For each amplitude condition, the script saves:
#     - RMSE boxplots: 'RMSE-<amplitude>.png'
#     - Correlation boxplots (if applicable): 'Corr-<amplitude>.png'
#
# Notes:
#   - Correlation plots are omitted when amplitude is fixed at 1 because
#     Pearson's correlation can not be calculated when no variability is present
#   - The figures use a color-blind friendly palette (Okabe-Ito palette)
#
# ------------------------------------------------------------------------------

#Load packages
library(ggplot2)
library(ggokabeito)
library(stringr)

path = './Results/'
N = 100 #Number of iterations
methods = c('LS-A', 'LS-S', 'ELS-S', 'FS')
amplitudes = c('amp1_1', 'amp2_4', 'amp4_6')

for (amp in amplitudes){ #Iterate over amplitude conditions
  #List all files
  results <- list.files(path, pattern = amp)
  #Store correlations
  corr.l <- vector('list', length = length(results))
  #Store RMSE
  rmse.l <- vector('list', length = length(results))
  
  for (f in 1:length(results)){
    #Load Rdata file
    load(paste0(path, results[f]))
    #Extract file conditions
    if (str_length(results[f]) == 33){ #ISI ~ U(2,4), U(4,6), U(6,8)
      isi <- substr(results[f], 4, 6)
      dur <- substr(results[f], 11, 13)
      lat <- substr(results[f], 18, 20)
    } else { #ISI ~ U[8,10]
      isi <- substr(results[f], 4, 7)
      dur <- substr(results[f], 12, 14)
      lat <- substr(results[f], 19, 21)
    }
    #RMSE in dataframe format
    rmse.l[[f]] <- data.frame(ISI = isi, latency = lat, duration = dur,
                              method = c(rep(methods[1], N), rep(methods[2], N), 
                                         rep(methods[3], N), rep(methods[4], N)), 
                              rmse = c(rmse[,1], rmse[,2], rmse[,3], rmse[,4]))
    #Correlations in dataframe format
    if (amp != 'amp1_1'){
      corr.l[[f]] <- data.frame(ISI = isi, latency = lat, duration = dur,
                                method = c(rep(methods[1], N), rep(methods[2], N), 
                                           rep(methods[3], N), rep(methods[4], N)), 
                                corr = c(corr[,1], corr[,2], corr[,3], corr[,4]))
    }
  }
  #Bind all data frames into one
  rmse.all <- do.call(rbind, rmse.l)
  if (amp != "amp1_1") {
    corr.all <- do.call(rbind, corr.l)
  }
  #Format data frames
  rmse.all$ISI <- sub("(\\d+)_(\\d+)", "\\1s-\\2s", rmse.all$ISI)
  rmse.all$latency <- sub("(\\d+)_(\\d+)", "lat: U(\\1,\\2)", rmse.all$latency)
  rmse.all$duration <- sub("(\\d+)_(\\d+)", "dur: U(\\1,\\2)", rmse.all$duration)
  rmse.all$method <- factor(rmse.all$method, levels = methods, ordered = TRUE)
  
  if (amp != 'amp1_1'){
    corr.all$ISI <- sub("(\\d+)_(\\d+)", "\\1s-\\2s", corr.all$ISI)
    corr.all$latency <- sub("(\\d+)_(\\d+)", "lat: U(\\1,\\2)", corr.all$latency)
    corr.all$duration <- sub("(\\d+)_(\\d+)", "dur: U(\\1,\\2)", corr.all$duration)
    corr.all$method <- factor(corr.all$method, levels = methods, ordered = TRUE)
  }
  #Plot RMSE
  ggplot(rmse.all, aes(x = ISI, y = rmse, fill = method)) +
    geom_boxplot(outlier.size = 1) +
    facet_grid(duration ~ latency) +
    scale_fill_okabe_ito() +
    labs(x = "Interstimulus interval (ISI) in seconds",
         y = "Root Mean Square Error (RMSE)",
         fill = "method",
         title = paste0("Simulation results: Amplitude ~ U(", 
                        substr(amp, 4, 4), ",", substr(amp, 6, 6), ")")) +
    theme_bw()
  ggsave(path = path, filename = paste0('RMSE-',amp,'.png'),
         width = 12, height = 6, units = "in", dpi = 300)
  #Plot correlations
  if (amp != 'amp1_1'){
    ggplot(corr.all, aes(x = ISI, y = corr, fill = method)) +
      geom_boxplot(outlier.size = 1) +
      facet_grid(duration ~ latency) +
      scale_fill_okabe_ito() +
      labs(x = "Interstimulus interval (ISI) in seconds",
           y = "Correlation",
           fill = "method",
           title = paste0("Simulation results: Amplitude ~ U(", 
                          substr(amp, 4, 4), ",", substr(amp, 6, 6), ")")) +
      theme_bw()
    ggsave(path = path, filename = paste0('Corr-',amp,'.png'),
           width = 12, height = 6, units = "in", dpi = 300)
  }
}













