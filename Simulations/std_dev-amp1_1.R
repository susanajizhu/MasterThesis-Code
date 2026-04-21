# ------------------------------------------------------------------------------
# Script: RawStmat-amp1_1.R
#
# Purpose:
#   Visualize raw estimates for fixed amplitude at 1 with boxplots and calculate
#   the standard deviation.
#   The script:
#     1) Loads results of fixed amplitude at 1 from the './Results/' directory
#     2) Extracts simulation conditions from filenames
#     3) Extracts raw estimate values
#     4) Formats results into dataframes with 4 columns: ISI, latency, duration,
#        and estimate value
#     5) Generates and saves boxplots comparing estimation methods in the
#        './Results/' directory
#     6) Calculates standard deviations and saves in a latex table
#
# Outputs:
#     - Raw estimates boxplot: 'RawStmat-amp1_1.png'
#     - Latex table with standard deviations of estimates: 'std_dev-amp1_1.tex'
#
# Notes:
#   - The figures use a color-blind friendly palette (Okabe-Ito palette)
#
# ------------------------------------------------------------------------------

#Load packages
library(ggplot2)
library(ggokabeito)
library(stringr)
library(dplyr)

path = './Results/'
N = 100 #Number of iterations
n_trial = 24 #Number of trials
methods = c('LS-A', 'LS-S', 'ELS-S', 'FS')
amp = 'amp1_1'

#List all files
results <- list.files(path, pattern = amp)

#Store estimates
stmat.l <- vector('list', length = length(results))

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
  
  LSA <- as.vector(stmat.all[,1,]) #Second dimension of the array is method
  LSS <- as.vector(stmat.all[,2,])
  ELSS <- as.vector(stmat.all[,3,])
  FS <- as.vector(stmat.all[,4,])

  #Store all estimates in a data frame
  stmat.l[[f]] <- data.frame(ISI = isi, latency = lat, duration = dur,
                             method = rep(methods, each = N * n_trial),
                             stmat = c(LSA, LSS, ELSS, FS))
}

#Bind all data frames into one
df_stmat <- do.call(rbind, stmat.l)

#Format data frames
df_stmat$ISI <- sub("(\\d+)_(\\d+)", "\\1s-\\2s", df_stmat$ISI)
df_stmat$latency <- sub("(\\d+)_(\\d+)", "\\1s-\\2s", df_stmat$latency)
df_stmat$duration <- sub("(\\d+)_(\\d+)", "\\1s-\\2s", df_stmat$duration)
df_stmat$method <- factor(df_stmat$method, levels = methods, ordered = TRUE)

#Plot estimates
ggplot(df_stmat, aes(x = ISI, y = stmat, fill = method)) +
  geom_boxplot(outlier.size = 1) +
  facet_grid(duration ~ latency) +
  geom_hline(yintercept = 1,  linetype = "dashed", color = 'red') + 
  scale_fill_okabe_ito() +
  labs(x = "Interstimulus interval (ISI) in seconds",
       y = "Raw estimated amplitudes",
       fill = "method",
       title = paste0("Variability in estimated amplitudes")) +
  theme_bw()
ggsave(path = path, filename = paste0('RawStmat-',amp,'.png'),
       width = 12, height = 6, units = "in", dpi = 300)

#Calculate standard deviations
sd_stmat <- df_stmat %>%
  group_by(ISI, latency, duration, method) %>%
  summarise(sd_stmat = round(sd(stmat),3), .groups = 'drop')

print(xtable(sd_stmat, type = "latex"), file = "std_dev-amp1_1.tex")


