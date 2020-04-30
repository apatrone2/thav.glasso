#' Generates a precision matrix list consisting of a precision matrix and a igraph layout
#'
#' @param d Dimension of the precision matrix 
#' @param graph Graph type
#' @param rangemin, rangemax The generation of the entries is based on sampling from a uniform distribution from [-rangemax, -rangemin] and [-rangemin, rangemax]
#' @param prob Probability that is used to generate the Gilbert graph, which defines the adjacency matrix of the graph
#' @return list containing the precision matrix and an igraph layout
#' @examples 
#' rd <- generateGraph(100)
#' network_rd <- graph.adjacency(rd[[1]], mode="undirected", weighted=TRUE, diag=FALSE)
#' plot.igraph(network_rd, layout=rd[[2]])
#' 
#' sf <- generateGraph(100, graph="scale-free")
#' network_sf <- graph.adjacency(sf[[1]], mode="undirected", weighted=TRUE, diag=FALSE)
#' plot.igraph(network_sf, layout=sf[[2]])
#' @export
generateGraph <- function(d, graph="random", rangemin=0.5, rangemax=0.9, prob=NULL)
{
  data <- huge.generator(d=d, prob=prob, graph=graph, verbose=FALSE)
  adj <- data$theta
  network <- graph.adjacency(round(adj, 4), mode="undirected", weighted=TRUE, 
                             diag=FALSE)
  network.layout <- layout_nicely(network)
  
  for(j in 1:(d - 1))
  {
    for(i in (j+1):d)
    {
      if (adj[i,j] != 0)
      {
        random <- sample(c(0,1), 1)
        if( random == 0)
        {
          adj[i,j] <- runif(1, min=rangemin, max=rangemax)
          adj[j,i] <- adj[i,j]
        }
        else
        {
          adj[i,j] <- -runif(1, min=rangemin, max=rangemax)
          adj[j,i] <- adj[i,j]
        }
      }
    }
  }
  diag(adj) <- rep(1, d)
  if( is.positive.definite(matrix(adj, ncol=d))==FALSE)
  {
    adj <- adj - diag( floor( 10 * eigen(adj)$values[d]) / 10, d)
  }
  
  
  return( list(round(adj, 4) / adj[1,1], network.layout))
}

