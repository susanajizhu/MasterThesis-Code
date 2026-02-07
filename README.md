# Statistics and Data Science: Single-trial estimation methods in fMRI data analysis
This repository contains the code developed for my Master's thesis about single trial estimation methods for fMRI data. The methods implemented and compared are:
1) Least Squares-All (LS-A)
2) Least Squares-Separate (LS-S)
3) Extended Least Squares-Separate (ELS-S)
4) Finite Bold Response-Separate (FS)

## Repository structure
A) Simulations
    1) CreateEv_Func.R: Generate event timing (EV) files
    2) RunMethods_Func.R: Run LS-A, LS-S, ELS-S, and FS methods 
    3) Simulation_Func.R: Full simulation process (makes use of create_ev() and run_methods() functions)
    4) Simulations.R: Run the simulations according to the conditions specified
    5) PlotResults.R: Plot RMSE and Pearson's correlation values
    6) "real-data" directory: contains an example of a real data set (same data was used for the empirical studies)
B) Empirical experiment
    1) AnalysisFlaker.R: Run LS-A, LS-S, ELS-S, and FS methods
    2) PlotFlaker.R: Plot observed BOLD signal vs predicted signal for one voxel
    3) "data" directory: contains the data set used for the analysis

## How to Run the Code
### 1. Requirements
The following R packages are required:
- fmristat
- niftiR6
- ggplot2
- ggokabeito
- stringr
### 2. Steps
A) Simulation studies
    1) Set conditions (ISI, amplitude, latency, and duration) and run the Simulations.R script 
    2) Run PlotResults.R script to see the RMSE and correlation values with boxplots
B) Empirical experiment
    1) Run AnalysisFlaker.R script to obtain estimates and output results
    2) Run PlotFlaker.R to compare observed BOLD signal and predicted BOLD signal by each method
