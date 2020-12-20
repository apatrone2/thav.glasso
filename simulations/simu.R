require("devtools") # to install thav.glasso
library("devtools")
install("../../thav.glasso") # contains the estimator and simulation methods
require("thav.glasso")
require("MASS") # to sample from multivariate normal
require("igraph") # to generate graphs
require("stargazer") # to produce latex outputs
require("huge") # contains adjacency matrix generation and StARS + RIC estimation
require("matrixcalc") # to compute eigenvalues
require("glasso") # to compute the graphical lasso
require("scalreg")


set.seed(18)
comparisonf1(300, 200, latex=TRUE, graph="random")
