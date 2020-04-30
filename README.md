# thav
Thresholded Adaptive Validation: Tuning the Graphical Lasso for Graph Recovery

This repository provides the implementations and simulations described in the paper Thresholded Adaptive Validation: Tuning the Graphical Lasso for Graph Recovery by Mike Laszkiewicz, Asja Fischer, and Johannes Lederer.

# Usage
The file (R/generate_precmatrix.R) contains the implementation of the precision matrix generation method used in the paper. 
The file R/av_estimation.R contains the implementations of the adaptive validation (AV) estimator, av_glasso, and the corresponding thresholded adaptive validation (thAV) estimator, thAV.estimator.
The implementation of the simulations can be found in R/def_simulations.R

# Simulations
The file simulations/simulation_study.Rmd provides the code for the execution of the simulations. 

The files simulations/adaptation_threshold.ipynb , simulations/visualization_weights.ipynb , and example_graphs.R generate plots and graphs, which were used in the mentioned Paper and can be found in the directory plots/ .
