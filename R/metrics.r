#' Compute the infloss on the off-diagonals of two matrices
#'
#' @param A Matrix
#' @param B Matrix with the same size as A
#' @return Maximal off-diagonal difference between A and B
#' @export
infloss <- function(A, B)
{
  dif <- abs(A - B)
  diag(dif) <- rep(0, dim(A)[1])
  return(max(dif))
}

#' Compute precision of a precision matrix estimator
#'
#' @param theta True precision matrix
#' @param theta.estimate Estimator of theta
#' @return Precision of theta.estimate with theta being the target
#' @export
precision <- function(theta, theta.estimate)
{
  adj.theta <- apply(theta, c(1,2), function(x) if (x!=0){return(1)} else {return(0)})
  diag(adj.theta) <- rep(0, dim(theta)[1])
  adj.estimate <- apply(theta.estimate, c(1,2), function(x) if (x!=0){return(1)} else {return(0)})
  diag(adj.estimate) <- rep(0, dim(theta)[1])
  num.edges.estimate <- sum(adj.estimate)/2
  num.edges.same <- sum (apply(adj.theta - adj.estimate, c(1,2), function(x) if(x!=0) {return(1)} else{return(0)}))/2
  return(num.edges.same/ num.edges.estimate)
}

#' Compute recall of a precision matrix estimator
#'
#' @param theta True precision matrix
#' @param theta.estimate Estimator of theta
#' @return Recall of theta.estimate with theta being the target
#' @export
recall <- function(theta, theta.estimate)
{
  adj.theta <- apply(theta, c(1,2), function(x) if (x!=0){return(1)} else {return(0)})
  diag(adj.theta) <- rep(0, dim(theta)[1])
  adj.estimate <- apply(theta.estimate, c(1,2), function(x) if (x!=0){return(1)} else {return(0)})
  diag(adj.estimate) <- rep(0, dim(theta)[1])
  num.edges.theta <- sum(adj.theta)/2
  num.edges.same <- sum (apply(adj.theta - adj.estimate, c(1,2), function(x) if(x!=0) {return(1)} else{return(0)}))/2
  return(num.edges.same/ num.edges.theta)
}

#' Compute performance of a precision matrix estimator
#'
#' @param theta True precision matrix
#' @param theta.estimate Estimator of theta
#' @param decimal Number of decimals such that entries of absolute size 1e-decimal are set to zero
#' @return f1 F_1-score of theta.estimate with theta being the target
#' @return recall Recall of theta.estimate with theta being the target
#' @return precision Precision of theta.estimate with theta being the target
#' @export
f1score <- function(theta, theta.estimate, decimal=4)
{
  #we don't want to include edges with values like e-17 which occure due to the 
  #precision matrix generation algorithm. Thats why we round...
  theta <- round(theta, decimal)
  theta.estimate <- round(theta.estimate, decimal)
  
  #generate the adjacency matrices. diagonal entries are set to 0
  adj.theta <- apply(theta, c(1,2), function(x) if (x!=0){return(1)} else {return(0)})
  diag(adj.theta) <- rep(0, dim(theta)[1])
  #we use a little trick here: set values of the adjacency matrix to 2. 
  #we do this so that we can count the number of same edges in an easy way:
  #if both edges are nonexistent -> adj.Theta - adj.estimate = 0 
  #if one of them is existent -> adj.Theta - adj.estimate \in \{ 1, -2\}
  #if both are existent -> adj.Theta - adj.estimate = -1
  adj.estimate <- apply(theta.estimate, c(1,2), function(x) if (x!=0){return(2)} else {return(0)})
  diag(adj.estimate) <- rep(0, dim(theta)[1])
  
  # every edge occurs exactly 2 times
  num.edges.theta <- sum(adj.theta)/2
  num.edges.estimate <- sum(adj.estimate/2)/2
  
  num.edges.same <- sum(apply(adj.theta - adj.estimate, c(1,2), function(x) if(x==-1) {return(1)} else{return(0)})) / 2
  
  recall <- num.edges.same / num.edges.theta
  if( num.edges.estimate > 0)
  {
    precision <- num.edges.same / num.edges.estimate
    f1 <- 2 * (precision * recall)/ (precision + recall)
  }
  else
  {
    precision <- 0
    f1 <- 0
  }
  
  return(list("f1"=f1, "recall"=recall, "precision"=precision))
}
