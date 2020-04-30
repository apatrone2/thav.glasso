install("..\\..\\thav.glasso")
require("thav.glasso")
require("MASS") # to sample from multivariate normal
require("igraph") # to generate graphs
require("stargazer") # to produce latex outputs
require("huge")
require("matrixcalc")
require("glasso")
require("hdi")

set.seed(52)
data("riboflavin")
data.ribo <- matrix(riboflavin$x, ncol=4088)
colnames( data.ribo) <- colnames( riboflavin$x)
variances <- diag( var(data.ribo))
variances <- sort( variances, decreasing=TRUE)
largest100 <- names(variances[1:100])

data.v100 <- data.ribo[, largest100]

dim( data.v100)

genes.scaled <- huge.npn( data.v100)

'av <- av_glasso( genes.scaled, seq(0.02, 0.2, length.out = 100), C=0.3)
thAV <- threshold(av)'
thav <- thAV.estimator(genes.scaled)

network.thAV <- graph.adjacency( thav, mode="undirected", weighted=TRUE, diag=FALSE)
network.layout <- layout.graphopt( network.thAV)
png("..\\plots\\ribo_thAV.png")
par(mar=c(0, 0, 0, 0))
plot.igraph( network.thAV, layout=network.layout)
dev.off()

stars <- huge.select( huge( genes.scaled, method = "mb", nlambda=30), criterion="stars", stars.thres=0.05)#$opt.icov
network.mb <- graph.adjacency( stars$refit, mode="undirected", weighted=TRUE, diag=FALSE)
network.layout <- layout.graphopt( network.mb)
png("..\\plots\\ribo_mb.png")
par(mar=c(0, 0, 0, 0))
plot.igraph( network.mb, layout=network.layout)
dev.off()
