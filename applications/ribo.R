require("devtools")
install("..\\..\\thav.glasso")
load_libraries()
require("hdi")

igraph_options(vertex.size       =2,
               vertex.label.cex  = 1,
               edge.width        = 1 ,
               vertex.label      ="",
               vertex.color      ="black",
               vertex.frame.color=NA,
               edge.color        ="coral1")

# Select 100 most-varying features
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

thav <- thAV.estimator(genes.scaled, mute=FALSE)
  
network.thAV <- graph.adjacency( thav, mode="undirected", weighted=TRUE, diag=FALSE)
network.layout <- layout.graphopt( network.thAV)
png("../plots/ribo_thAV.png", width=1900, height=1600, res=200) # saves the plot
par(mar=c(0, 0, 0, 0))
plot.igraph( network.thAV, layout=network.layout)
dev.off() # saves the plot

#lambda = 0.5
thav2 <- thAV.estimator(genes.scaled, lambda=0.5)

network.thAV2 <- graph.adjacency(thav2, mode="undirected", weighted=TRUE, diag=FALSE)
network.layout2 <- layout.graphopt( network.thAV2)
png("../plots/ribo_thAV2.png", width=1900, height=1600, res=200) # saves the plot
par(mar=c(0, 0, 0, 0))
plot.igraph( network.thAV2, layout=network.layout2)
dev.off()

stars <- huge.select( huge( genes.scaled, method = "mb", nlambda=30), criterion="stars", stars.thres=0.05)
network.mb <- graph.adjacency( stars$refit, mode="undirected", weighted=TRUE, diag=FALSE)
network.layout <- layout.graphopt( network.mb)
#png("..\\plots\\ribo_mb.png") # saves the plot
par(mar=c(0, 0, 0, 0))
plot.igraph( network.mb, layout=network.layout)
#dev.off() # saves the plot
