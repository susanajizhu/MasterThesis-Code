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
#   - A figure showing observed vs predicted BOLD signals for each
#     estimation method (LS-A, LS-S, ELS-S, FS)
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
ev_TR <- ev$onsets/2 #Select onset times and transform seconds to TR

#Plot results
png("./Results/FlakerPlot.png", width = 12, height = 12, units = 'in', res = 300)
par(mfrow = c(4, 1), mar = c(4, 4, 3, 1), oma = c(10, 4, 2, 4))

#LS-A
plot(bold_max,type='l',lwd=2,col='darkblue',bty='n',ylab='BOLD signal',xlab='Time (TR)', main = 'LS-A')
lines(predict(output_max$output_LA), col = 'darkorange', lwd = 1)
abline(v = ev_TR, lty = 1, col = "grey")
#LS-S
plot(bold_max,type='l',lwd=2,col='darkblue',bty='n',ylab='BOLD signal',xlab='Time (TR)', main = 'LS-S')
for (i in 1:length(output_max$output_LS)){
  lines(predict(output_max$output_LS[[i]]), col = 'darkorange', lwd = 1)
}
abline(v = ev_TR, lty = 1, col = "grey")
#ELS-S
plot(bold_max,type='l',lwd=2,col='darkblue',bty='n',ylab='BOLD signal',xlab='Time (TR)', main = 'ELS-S')
for (i in 1:length(output_max$output_ELS)){
  lines(predict(output_max$output_ELS[[i]]), col = 'darkorange', lwd = 1)
}
abline(v = ev_TR, lty = 1, col = "grey")
#FS
plot(bold_max,type='l',lwd=2,col='darkblue',bty='n',ylab='BOLD signal',xlab='Time (TR)', main = 'FS')
for (i in 1:length(output_max$output_FS)){
  lines(predict(output_max$output_FS[[i]]), col = 'darkorange', lwd = 1)
}
abline(v = ev_TR, lty = 1, col = "grey")

legend("bottom",
       inset = -0.80,
       legend = c("Observed BOLD signal", "Predicted BOLD signal"),
       col = c("darkblue", "darkorange"),
       lwd = 2,
       horiz = TRUE,
       bty = "n",
       xpd = NA)

dev.off()

