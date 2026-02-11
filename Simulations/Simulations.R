# ------------------------------------------------------------------------------
# Script: Simulations.R
#
# Purpose:
#   Run fMRI simulation experiments under varying conditions of ISI, amplitude, 
#   latency, and duration using single-trial estimation methods.
#   This script:
#     1) Loads required packages (fmristat, niftiR6)
#     2) Loads required functions (create_ev(), run_methods(), simulation())
#     3) Runs the function simulation() according to the conditions defined by
#        the user
#
# Outputs:
#   Each simulation generates an .RData file in the directory "./Results/" with:
#     - Simulated BOLD signals
#     - True beta values (true amplitudes)
#     - Estimated beta values (estimated amplitudes)
#     - Root mean square errors
#     - Pearson's correlation coefficients
#     - Outputs from the lm objects from all single trial estimation methods
#
# ------------------------------------------------------------------------------

#Load packages
library(fmristat)
library(niftiR6)

#Load functions
source("CreateEv_Func.R") #create_ev()
source("RunMethods_Func.R") #run_methods()
source("Simulation_Func.R") #simulation()

# ================= Define conditions and run simulation ======================

#Examples of conditions
c1 <- simulation(N = 100,
                 isi_a = 8, isi_b = 10,
                 dur_a = 2, dur_b = 2,
                 lat_a = 0, lat_b = 0,
                 amp_a = 1, amp_b = 1)

c2 <- simulation(N = 100,
                 isi_a = 8, isi_b = 10,
                 dur_a = 2, dur_b = 2,
                 lat_a = 0, lat_b = 0,
                 amp_a = 2, amp_b = 4)

c3 <- simulation(N = 100,
                 isi_a = 8, isi_b = 10,
                 dur_a = 2, dur_b = 2,
                 lat_a = 0, lat_b = 0,
                 amp_a = 4, amp_b = 6)










