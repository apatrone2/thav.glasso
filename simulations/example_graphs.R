require("devtools")
install("../thav.glasso")
require("thav.glasso")
require("MASS") # to sample from multivariate normal
require("igraph") # to generate graphs
require("huge") # to generate adjacency matrices
require("matrixcalc") # to calculate eigenvalues
require("glasso") # to compute the graphical lasso

igraph_options(vertex.size       =2,
               vertex.label.cex  = 1,
               edge.width        = 1 ,
               vertex.label      ="",
               vertex.color      ="black",
               vertex.frame.color=NA,
               edge.color        ="coral1")

d <- 200
set.seed(12)
theta_2d <- generateGraph(d, prob=2/d)[[1]]
network_2d <- graph.adjacency(apply(theta_2d, c(1,2), function(x) x!=0), mode="undirected", weighted=TRUE, diag=FALSE)
layout_2d <- layout_nicely(network_2d)
png("../plots/density2d.png")
par(mar = c(0, 0, 0, 0) + 0.1)
plot.igraph(network_2d, layout=layout_2d )
dev.off()

theta_3d <- generateGraph(d, prob=3/d)[[1]]
network_3d <- graph.adjacency(apply(theta_3d, c(1,2), function(x) x!=0), mode="undirected", weighted=TRUE, diag=FALSE)
layout_3d <- layout_nicely(network_3d)
png("../plots/density3d.png")
par(mar = c(0, 0, 0, 0) + 0.1)
plot.igraph(network_3d, layout=layout_3d )
dev.off()

theta_4d <- generateGraph(d, prob=4/d)[[1]]
network_4d <- graph.adjacency(apply(theta_4d, c(1,2), function(x) x!=0), mode="undirected", weighted=TRUE, diag=FALSE)
layout_4d <- layout_nicely(network_4d)
png("../plots/density4d.png")
par(mar = c(0, 0, 0, 0) + 0.1)
plot.igraph(network_4d, layout=layout_4d )
dev.off()

