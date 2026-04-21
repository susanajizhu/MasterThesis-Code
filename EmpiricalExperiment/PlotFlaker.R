# ------------------------------------------------------------------------------
# Script: PlotFlaker.R
#
# Purpose:
#   Visualize the observed BOLD signal and the predicted BOLD signal for 
#   the voxel with maximum activation (highest z-statistic) in the contrast of
#   incongruent > congruent (COPE 5). 
#   Event/trial onsets are displayed as vertical reference lines.
#
# Outputs:
#   - Two figures showing observed vs predicted BOLD signals for each
#     estimation method (LS-A, LS-S, ELS-S, FS):
#     a) One plot displays full time series information
#     b) Second plot displays time series around a window of 14 seconds after 
#        stimulus onset
#
# Additional notes:
#   - The script reads the 'output_max-run1.RData' file generated from
#     AnalysisFlaker.R
# ------------------------------------------------------------------------------

#Read file
load("./Results/output_max-run1.RData")

#Extract BOLD signal
bold_max <- output_max$signal

#Read custom timing files
ev1 <- data.frame(read.delim("./data/sub-01/run-1.feat/custom_timing_files/ev1.txt", header = FALSE))
colnames(ev1) <- c('onsets', 'duration', 'weight')
ev2 <- data.frame(read.delim("./data/sub-01/run-1.feat/custom_timing_files/ev2.txt", header = FALSE))
colnames(ev2) <- c('onsets', 'duration', 'weight')
ev <- rbind(ev1, ev2)
ev$onsets_TR <- ev$onsets/2 #Select onset times, transform seconds to TR

#Plot results: Full times series
png("./Results/FlakerPlot_all.png", width = 12, height = 12, units = 'in', res = 300)
par(
  mfrow = c(4, 1),
  mar = c(4, 4, 3, 1),
  oma = c(10, 4, 2, 4),
  cex.axis = 1.3,
  cex.lab = 1.3,
  cex.main = 1.6
)

#LS-A
plot(bold_max,type='l',lwd=1,col='darkblue',bty='n',ylab='BOLD signal',xlab='Time (TR)', main = 'LS-A')
lines(predict(output_max$output_LA), col = 'darkorange', lwd = 1)
abline(v = ev$onsets_TR, lty = 1, col = "grey")

#LS-S
plot(bold_max,type='l',lwd=1,col='darkblue',bty='n', ylab='BOLD signal',xlab='Time (TR)', main='LS-S')
for (i in 1:length(output_max$output_LS)){
  lines(predict(output_max$output_LS[[i]]), col = 'darkorange', lwd = 1)
}
abline(v = ev$onsets_TR, lty = 1, col = "grey")

#ELS-S
plot(bold_max,type='l',lwd=1,col='darkblue',bty='n',ylab='BOLD signal',xlab='Time (TR)', main = 'ELS-S')
for (i in 1:length(output_max$output_ELS)){
  lines(predict(output_max$output_ELS[[i]]), col = 'darkorange', lwd = 1)
}
abline(v = ev$onsets_TR, lty = 1, col = "grey")

#FS
plot(bold_max,type='l',lwd=1,col='darkblue',bty='n',ylab='BOLD signal',xlab='Time (TR)', main = 'FS')
for (i in 1:length(output_max$output_FS)){
  lines(predict(output_max$output_FS[[i]]), col = 'darkorange', lwd = 1)
}
abline(v = ev$onsets_TR, lty = 1, col = "grey")

legend("bottom",
       inset = -0.80,
       legend = c("Observed BOLD signal", "Predicted BOLD signal"),
       col = c("darkblue", "darkorange"),
       lwd = 2,
       horiz = TRUE,
       bty = "n",
       xpd = NA,
       cex = 1.5)

dev.off()

#Plot results: Window of 14 seconds around time series
png("./Results/FlakerPlot.png", width = 12, height = 12, units = 'in', res = 300)
par(
  mfrow = c(4, 1),
  mar = c(4, 4, 3, 1),
  oma = c(10, 4, 2, 4),
  cex.axis = 1.3,
  cex.lab = 1.3,
  cex.main = 1.6
)

#LS-A
plot(bold_max,type='l',lwd=1,col='darkblue',bty='n',ylab='BOLD signal',xlab='Time (TR)', main = 'LS-A')
lines(predict(output_max$output_LA), col = 'darkorange', lwd = 1)
abline(v = ev$onsets_TR, lty = 1, col = "grey")

#LS-S
plot(bold_max,type='l',lwd=1,col='darkblue',bty='n', ylab='BOLD signal',xlab='Time (TR)', main='LS-S')
for (i in 1:length(output_max$output_LS)){
  pred <- predict(output_max$output_LS[[i]])
  start <- ev$onsets_TR[i] #Find onset time of current trial
  end <- ev$onsets_TR[i] + 7 #Window of 14 seconds after trial onset
  inter <- start:end
  pred_int <- pred
  pred_int[-inter] <- NA
  lines(pred_int, col='darkorange', lwd=1)
}
abline(v = ev$onsets_TR, lty = 1, col = "grey")

#ELS-S
plot(bold_max,type='l',lwd=1,col='darkblue',bty='n',ylab='BOLD signal',xlab='Time (TR)', main = 'ELS-S')
for (i in 1:length(output_max$output_ELS)){
  pred <- predict(output_max$output_ELS[[i]])
  start <- ev$onsets_TR[i] #Find onset time of current trial
  end <- ev$onsets_TR[i] + 7 #Window of 14 seconds after trial onset
  inter <- start:end
  pred_int <- pred
  pred_int[-inter] <- NA
  lines(pred_int, col='darkorange', lwd=1)
}
abline(v = ev$onsets_TR, lty = 1, col = "grey")

#FS
plot(bold_max,type='l',lwd=1,col='darkblue',bty='n',ylab='BOLD signal',xlab='Time (TR)', main = 'FS')
for (i in 1:length(output_max$output_FS)){
  pred <- predict(output_max$output_FS[[i]])
  start <- ev$onsets_TR[i] #Find onset time of current trial
  end <- ev$onsets_TR[i] + 7 #Window of 14 seconds after trial onset
  inter <- start:end
  pred_int <- pred
  pred_int[-inter] <- NA
  lines(pred_int, col='darkorange', lwd=1)
}
abline(v = ev$onsets_TR, lty = 1, col = "grey")

legend("bottom",
       inset = -0.80,
       legend = c("Observed BOLD signal", "Predicted BOLD signal"),
       col = c("darkblue", "darkorange"),
       lwd = 2,
       horiz = TRUE,
       bty = "n",
       xpd = NA,
       cex = 1.5)

dev.off()

