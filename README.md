# thav
Thresholded Adaptive Validation: Tuning the Graphical Lasso for Graph Recovery

This repository provides the implementations and simulations described in the paper Thresholded Adaptive Validation: Tuning the Graphical Lasso for Graph Recovery by Mike Laszkiewicz, Asja Fischer, and Johannes Lederer.

# Getting started

## Installation 
```R
# install.packages("devtools")
# library("devtools")
devtools::install_github('UUPharmacometrics/xpose')
## Further requirements:
# install.packages("MASS") # to sample from multivariate normal
# install.packages("igraph") # to generate graphs
# install.packages("stargazer") # to produce latex outputs
# intsall.packages("huge") # contains adjacency matrix generation and StARS + RIC estimation
# install.packages("matrixcalc") # to compute eigenvalues
# install.packages("glasso") # to compute the graphical lasso

## or simply
# install.packages(c("MASS", "igraph", "stargazer", "huge", "matrixcalc", "glasso"))
```
## Load thav.glasso and requirements
```R
library("thav.glasso", "MASS", "igraph", "stargazer", "huge", "matrixcalc", "glasso")
```

# Usage
## Generate Precision matrix (and layout for graphical presentation using igraph)
```R
graph <- generateGraph(d=200, graph="random")
theta <- graph[[1]]
layout <- graph[[2]]
```

## Visualize the graph
```R
# settings for visualization
igraph_options(vertex.size       =2,
               vertex.label.cex  = 1,
               edge.width        = 1 ,
               vertex.label      ="",
               vertex.color      ="black",
               vertex.frame.color=NA,
               edge.color        ="coral1")

# plot the graph
adjacency_matrix <- apply(theta, c(1,2), function(x) x!=0)
network <- graph.adjacency(adjacency_matrix, mode="undirected", weighted=TRUE, diag=FALSE)
plot.igraph(network, layout=layout)
```

## Apply the AV and thAV
```R
# generate Data
data <- mvrnorm(300, mu=rep(0, 200), Sigma=solve(theta))

# calculate the AV
av <- av_glasso(data)
av_estimator <- av$Thetahat
av_tuningparameter <- av$TuningParameter

# calculate the thAV
thav <- thAV.estimator(data)
```

## Plot thAV
```R
adjacency_matrix_thav <- apply(thav, c(1,2), function(x) x!=0)
network_thav <- graph.adjacency(adjacency_matrix_thav, mode="undirected", weighted=TRUE, diag=FALSE)
plot.igraph(network_thav, layout=layout)
```

# Simulations
The file `R/generate_precmatrix.R` contains the implementation of the precision matrix generation method used in the paper. 
The file `R/av_estimation.R` contains the implementations of the adaptive validation (AV) estimator, av_glasso, and the corresponding thresholded adaptive validation (thAV) estimator, thAV.estimator.
The implementation of the simulations can be found in `R/def_simulations.R`

The file `simulations/simulation_study.Rmd` provides the code for the execution of the simulations. 

The files `simulations/adaptation_threshold.ipynb` , `simulations/visualization_weights.ipynb` , and `example_graphs.R` generate plots and graphs, which were used in the mentioned Paper and can be found in the directory `plots/`.
