library(devtools)
install_github("zdk123/SpiecEasi")
install("..\\..\\thav.glasso")
library("SpiecEasi")
library("thav.glasso")
require("phyloseq")
require("igraph")
require("huge")
require("glasso")

# https://github.com/zdk123/SpiecEasi
set.seed(132)
data('amgut2.filt.phy')
se.mb.amgut2 <- spiec.easi(amgut2.filt.phy, method='mb', lambda.min.ratio=1e-2,
                           nlambda=20, icov.select.params=list(rep.num=50))
ig2.mb <- adj2igraph(se.mb.amgut2$refit$stars,  vertex.attr=list(name=taxa_names(amgut2.filt.phy)))
plot_network(ig2.mb, amgut2.filt.phy, type='taxa', color="Rank4", layout.method=layout.graphopt)

# our method:
thav <- thAV.estimator(huge.npn(se.mb.amgut2$est$data), threshold=3)
adj_thav <- apply( thav, c(1,2), function(x) if(x!=0){ return(1)} else{return(0)})
ig.thav <- adj2igraph( adj_thav, vertex.attr=list(name=taxa_names(amgut2.filt.phy)))

plot_network(ig.thav, amgut2.filt.phy, type="taxa", color="Rank4")

# delete isolated vertices
plot_network( delete.vertices(ig.thav, degree(ig.thav)==0), amgut2.filt.phy, type="taxa", color="Rank4")#, layout.method=layout.graphopt)
plot_network( delete.vertices(ig2.mb, degree(ig2.mb)==0), amgut2.filt.phy, type="taxa", color="Rank4")#, layout.method=layout.graphopt)


# calcuate the similar graph
similar_adj <- function( theta1, theta2)
{
  d <- dim(theta1)[1]
  similar <- matrix( rep(0, d^2), ncol=d)
  for(j in 1:(d - 1))
  {
    for(i in (j + 1):d)
    {
      if( theta1[i,j] != 0 && theta2[i,j] != 0)
      {
        similar[i,j] <- similar[j, i] <- 1 
      }
    }
  }
  return( similar)
}

sim_comp <- similar_adj(thav, se.mb.amgut2$refit$stars)
sim_graph <- adj2igraph( sim_comp, vertex.attr=list(name=taxa_names(amgut2.filt.phy)))

# plots
set.seed(12)
plot_network(delete.vertices(ig.thav, degree(ig.thav)==0), amgut2.filt.phy, type="taxa",
              color="Rank3", shape="Rank3")
set.seed(44)
plot_network(delete.vertices(ig2.mb, degree(ig2.mb)==0), amgut2.filt.phy, type="taxa",
              color="Rank3", shape="Rank3")
set.seed(4)
plot_network(delete.vertices(sim_graph, degree(sim_graph)==0), amgut2.filt.phy, type="taxa",
              color="Rank3", shape="Rank3")


# compare edges of the thav and SPIEC-EASI:
paste0("thav has ", length(which(thav[upper.tri(thav)] != 0)), " edges")
paste0("similar graph has ", length(which(sim_comp[upper.tri(sim_comp)] != 0)), " edges")
paste0("SPIEC-EASI graph has ", length(which(se.mb.amgut2$refit$stars[upper.tri(se.mb.amgut2$refit$stars)] != 0)), " edges")


# threshold = 2
set.seed(492)
thav2 <- thAV.estimator(huge.npn(se.mb.amgut2$est$data), lambda=2)
adj_thav2 <- apply(thav2, c(1,2), function(x) if(x!=0){ return(1)} else{return(0)})
ig.thav2 <- adj2igraph(adj_thav2, vertex.attr=list(name=taxa_names(amgut2.filt.phy)))

plot_network(delete.vertices(ig.thav2, degree(ig.thav2)==0), amgut2.filt.phy, type="taxa", color="Rank3", shape="Rank3")#, layout.method = layout.graphopt)

sim_comp <- similar_adj(thav2, se.mb.amgut2$refit$stars)
sim_graph <- adj2igraph( sim_comp, vertex.attr=list(name=taxa_names(amgut2.filt.phy)))

# compare edges of the thav and SPIEC-EASI:
paste0("thav2 has ", length(which(thav2[upper.tri(thav2)] != 0)), " edges")
paste0("similar graph has ", length(which(sim_comp[upper.tri(sim_comp)] != 0)), " edges")
paste0("SPIEC-EASI graph has ", length(which(se.mb.amgut2$refit$stars[upper.tri(se.mb.amgut2$refit$stars)] != 0)), " edges")
