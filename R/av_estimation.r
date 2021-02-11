#' Loads the necessary libraries, which are: huge, MASS, thav.glasso, scalreg, igraph, matrixcalc, genscore, scio, stargazer, glasso
#'
#' @export
load_libraries <- function()
{
  libraries <- c("huge", "MASS", "thav.glasso", "scalreg", "igraph", "matrixcalc", "genscore", "scio", "stargazer", "glasso")
  lapply(libraries, library, character.only=TRUE)
}

#' Turns a covariance matrix into the correlation matrix
#'
#' @param theta Covariance matrix
#' @return Correlation matrix
#' @export
normalize_theta <- function(theta)
{
  for(j in 1:dim(theta)[1])
  {
    for(i in 1:dim(theta)[1])
    {
      theta[i,j] <- theta[i,j] / (sqrt(theta[i,i]) * sqrt(theta[j, j]))
    }
  }
  return(theta)
}

#' Computes the adaptive validation graphical lasso
#' 
#' @param data Data matrix of size n x d 
#' @param C constant C that is used to calibrate the AV regularization parameter
#' @return The selected regularization parameter
#' @return The resulting graphical lasso estimate
#' @examples 
#' theta <- generateGraph(d=100)[[1]]
#' Sigma <- solve(theta)
#' data <- scale(mvrnorm(n=200, mu=rep(0, 100), Sigma=Sigma))
#' av <- av_glasso(data)
#' av_tp <- av$TuningParameter
#' av_est <- av$Thetahat
#' @export
av_glasso <- function(data, C=0.7) 
{
  empvar <- var(data)
  
  # Automatically set max(seq_r) <- tp that sets all offdiagonal entries to zero
  max_r <- max( abs(empvar[upper.tri(empvar)]))
  seq_r <- seq(0.05, max_r, length.out = 40)
  
  r <- max(seq_r)
  estimates <- list()
  
  #calculate glasso[max(r)] and setup the warm start
  glasso <- glasso(empvar, rho=r, penalize.diagonal=FALSE)
  warm <- list(glasso$wi, glasso$w)
  estimates[[length(seq_r)]] <- normalize_theta(warm[[1]]) #warm[[1]] #
  
  counter_r <- length(seq_r) #used to fill the list estimates with corresponding estimates
  while(r > min(seq_r))
  {
    #compute glasso[r] using a warm start
    counter_r <- counter_r - 1
    r <- seq_r[counter_r]
    newwarm <- glasso(empvar, rho=r, penalize.diagonal=FALSE, start="warm", wi.init=warm[[1]], w.init=warm[[2]])
    estimates[[counter_r]] <- normalize_theta(newwarm$wi) #newwarm$wi #
    #warm <- list(normalize_theta(newwarm$wi), normalize_theta(newwarm$w))
    warm <- list(newwarm$wi, newwarm$w)
    
    counter_r_dash <- length(seq_r)
    r_dash <- seq_r[counter_r_dash]
    while(r_dash > r)
    {
      #if l(...) > Cr_dash + Cr
      if(avtest(estimates, seq_r, counter_r, counter_r_dash, C)) 
      {
        return(list(TuningParameter=seq_r[counter_r + 1], Thetahat=estimates[[counter_r + 1]]))
      }
      counter_r_dash <- counter_r_dash - 1
      r_dash <- seq_r[counter_r_dash]
    }
  }
  estimates[[1]] <- normalize_theta(glasso(empvar, rho=min(seq_r), penalize.diagonal=FALSE, start="warm", wi.init=warm[[1]], w.init=warm[[2]])$wi)
  return(list(TuningParameter=seq_r[1], Thetahat=estimates[[1]])) 
}

#' Computes the test that is used in the AV calibration scheme
#' 
#' @param estimates List of the graphical lasso estimators so far
#' @param seq_r List of possible regularization paraeters
#' @param counter_r position of the current r in r_seq (see algorithm)
#' @param counter_r_dash position of the current r_dash in seq_r (see algorithm)
#' @param C constant C that is used for the AV calibration
#' @return TRUE if l(...) > Cr_dash + Cr, else FALSE
#' @export
avtest <- function(estimates, seq_r, counter_r, counter_r_dash, C)
{
  infloss(estimates[[counter_r]], estimates[[counter_r_dash]]) > C * (seq_r[counter_r] + seq_r[counter_r_dash])
}

#' Thresholds the AV estimator
#' 
#' @param av List AV containing the regularization parameter and the estimator, which is obtained by applying the av_glasso function
#' @param C Constant C that was used to calibrate the AV
#' @param lambda Threshold factor lambda: We threshold the AV solution by lambda*C*av$TuningParameter
#' @return The thresholded AV estimator
#' @export
threshold <- function( av, C=0.7, lambda=1)
{
  return( apply( av$Thetahat, c(1,2), function(x) if( abs(x) <= lambda*C*av$TuningParameter){ return(0)} else{ return(x)}))
}

#' Computes the thAV estimator
#' 
#' @param data Data matrix of size n x d 
#' @param C constant C that is used to calibrate the AV regularization parameter
#' @param lambda Threshold factor lambda: We threshold the AV solution by lambda*C*av$TuningParameter
#' @param mute If FALSE, the method prints the selected AV regularization parameter. The default value is TRUE.
#' @return The thresholded AV estimator
#' @examples 
#' theta <- generateGraph(d=100)[[1]]
#' Sigma <- solve(theta)
#' data <- scale(mvrnorm(n=200, mu=rep(0, 100), Sigma=Sigma))
#' thav <- thAV.estimator(data)
#' @export
thAV.estimator <- function( data, C=0.7, lambda=1, mute=TRUE)
{
  av <- av_glasso(data, C) 
  
  if(!mute)
  {
    print(paste0("Selected Regularization Parameter: ", round(av$TuningParameter,2)))
  }
  return( threshold(av, C, lambda))
}


#' Computes the rSME via eBIC
#'
#' @param data Data matrix of size n x d
#' @param gamma Parameter that is used in the eBIC criterion. Usually gamma=0.5 or gamma=1. The default value is set to 1.
#' @return rSME estimator that is tuned via eBIC
#' @export
score_ebic <- function(data, gamma=1)
{
  d <- dim(data)[2]
  domain <- make_domain("R", p=d)
  elts <- get_elts(NULL, data, "gaussian", domain)
  max_r <- lambda_max(elts, "symmetric")
  seq_r <- seq(0.05, max_r, length.out = 40)
  best_eBIC <- 10^5
  for(j in 1:length(seq_r))
  {
    result <- get_results(elts, "symmetric", seq_r[j])
    if(eBIC(result, elts, gammas=gamma)[1] < best_eBIC)
    {
      ebic_est <- result$K
      best_eBIC <- eBIC(result, elts, gammas=gamma)[1]
    }
  }
  return(ebic_est)
}

#' Computes the adaptive validation regularized score matching estimator rSME (AV)
#' 
#' @param data Data matrix of size n x d
#' @param C constant C that is used to calibrate the AV regularization parameter
#' @return The selected regularization parameter
#' @return The resulting AV regularized score matching estimator rSME (AV)
#' @export
av_rsme <- function(data, C=2)
{
  empvar <- var(data)
  seq_r <- seq(5e-4, 0.3, length.out = 40)
  #initialization
  r <- max(seq_r)
  
  domain <- make_domain("R", p=dim(empvar)[1])
  estimates <- estimate(data, "gaussian", domain, scale="sd", return_raw=TRUE, lambda1s=seq_r,
                        verbose=FALSE, verbosetext=FALSE)$raw_estimate
  estimates <- rev(estimates) # genscore implementation returns a list of estimates 
                              # where the first entry corresponds to the largest TP estimate
  # normalize each theta
  for(j in 1:40)
  {
    estimates[[j]] <- normalize_theta(estimates[[j]])
  }
  counter_r <- length(seq_r) #used to fill the list estimates with corresponding estimates
  while(r != min(seq_r))
  {
    counter_r <- counter_r - 1
    r <- seq_r[counter_r]

    counter_r_dash <- length(seq_r)
    r_dash <- seq_r[counter_r_dash]
    while(r_dash > r)
    {
      #if l(...) > Cr_dash + Cr
      if(avtest(estimates, seq_r, counter_r, counter_r_dash, C))
      {
        return(list(TuningParameter=seq_r[counter_r + 1], Thetahat=estimates[[counter_r + 1]]))
      }
      counter_r_dash <- counter_r_dash - 1
      r_dash <- seq_r[counter_r_dash]
    }
  }
  return(list(TuningParameter=seq_r[1], Thetahat=estimates[[1]]))
}

#' Computes the trace of a matrix
#'
#' @param matrix Matrix
#' @return Trace of matrix
#' @export
tr <- function(matrix)
{
  return(sum(diag(matrix)))
}

#' Computation of the Bregman loss as defined in "Weidong Liu, Xi Luo, Fast and adaptive sparse precision matrix estimation in high dimensions, Journal of Multivariate Analysis, Volume 135, 2015, Pages 153-162"
#'
#' @param sigma Covariance matrix Sigma
#' @param theta Precision matrix Theta
#' @return Bregman loss
bregman_loss <- function(sigma, theta)
{
  return(tr(sigma %*% theta) - log(det(theta)))
}

#' Computation of the Sparse Columnwise Inverse Operator (SCIO) estimator tuned via the Bregman loss: We pick the regularization parameter that minimizes the Bregman loss.
#'
#' @param data Data matrix of size n x d 
#' @return SCIO estimator
#' @export 
scio_bregman <- function(data)
{
  empvar <- var(data)
  seq_r <- seq(0.05, max(abs(empvar[upper.tri(empvar)])), length.out = 40)
  
  bregman <- rep(0, length(seq_r))
  for(j in 1:length(seq_r))
  {
    bregman[j] <- bregman_loss(empvar, scio(empvar, lambda=seq_r[j])$w)
  }
  best_r <- seq_r[which.min(bregman)]
  return(scio(empvar, lambda=best_r)$w)
}

#' Computation of the adaptive validation SCIO
#' 
#' @param data Data matrix of size n x d
#' @param C constant C that is used to calibrate the AV regularization parameter
#' @return The selected regularization parameter
#' @return The resulting AV-SCIO estimate
#' @export
av_scio <- function(data, C=0.7)
{
  empvar <- var(data)
  max_r <- max(abs(empvar[upper.tri(empvar)]))
  seq_r <- seq(0.05, max_r, length.out = 40)
  #initialization
  r <- max(seq_r)
  estimates <- list()
  
  
  estimates[[length(seq_r)]] <- normalize_theta(scio(empvar, lambda=r)$w)
  
  counter_r <- length(seq_r) #used to fill the list estimates with corresponding estimates
  while(r != min(seq_r))
  {
    counter_r <- counter_r - 1
    r <- seq_r[counter_r]
    estimates[[counter_r]] <- normalize_theta(scio(empvar, lambda=r)$w)
    
    
    counter_r_dash <- length(seq_r)
    r_dash <- seq_r[counter_r_dash]
    while(r_dash > r)
    {
      #if l(...) > Cr_dash + Cr
      if(avtest(estimates, seq_r, counter_r, counter_r_dash, C))
      {
        return(list(TuningParameter=seq_r[counter_r + 1], Thetahat=estimates[[counter_r + 1]]))
      }
      counter_r_dash <- counter_r_dash - 1
      r_dash <- seq_r[counter_r_dash]
    }
  }
  return(list(TuningParameter=seq_r[1], Thetahat=estimates[[1]]))
}


#' Computes the reweighted l1 regularization to compute the graphical lasso with power law regularization (Liu & Ihler(2011)).
#' 
#' @param est Estimation of Theta from the previous iteration
#' @param alpha Regularization parameter
#' @return The reweighted penalties as a dxd matrix. 
#' @export
penalty_matrix <- function(est, alpha)
{
  d <- dim(est)[1]
  eps <- diag(est) 
  lambda <- matrix(rep(0, d*d), ncol=d)
  for(j in 1:(d-1))
  {
    for(i in j:d)
    {
      if(i==j)
      {
        lambda[i, i] <- 0
      }
      else
      {
        lambda[i, j] <- lambda[j, i] <- alpha * ( 1/(sum(abs(est[i, ])) + eps[i]) + 1/(sum(abs(est[j, ])) + eps[j]))
      }
    }
  }
  return(lambda)
}

#' Computes the graphical lasso with power law regularization sf_glasso (Liu & Ihler(2011)).
#' 
#' @param data Data matrix of size n x d
#' @param alpha Regularization parameter
#' @return Graphical lasso with power law regularization
#' @export
sf_glasso <- function(data, alpha)
{
  d <- dim(data)[2]
  est <- diag(d)
  empvar <- var(data)
  for(j in 1:3)
  {
    rho <- penalty_matrix(est, alpha)
    est <- glasso(empvar, rho=rho)$wi
  }
  return(est)
}

#' Computes the graphical lasso with power law regularization and regularization parameter tuned via adaptive validation sf_glasso (AV).
#' 
#' @param data Data matrix of size n x d
#' @param C Constant C used to tune the thAV. Default value is 1.5
#' @return Graphical lasso with power law regularization and regularization parameter tuned via adaptive validation
#' @export 
av_sf_glasso <- function(data, C=1.5)
{
  empvar <- var(data)
  max_r <- 0.8#max(abs(empvar[upper.tri(empvar)]))
  seq_r <- seq(0.01, max_r, length.out = 40)
  #initialization
  r <- max(seq_r)
  estimates <- list()
  
  estimates[[length(seq_r)]] <- normalize_theta(sf_glasso(data, r))
  
  counter_r <- length(seq_r) #used to fill the list estimates with corresponding estimates
  while(r != min(seq_r))
  {
    counter_r <- counter_r - 1
    r <- seq_r[counter_r]
    estimates[[counter_r]] <- normalize_theta(sf_glasso(data, r))
    
    counter_r_dash <- length(seq_r)
    r_dash <- seq_r[counter_r_dash]
    while(r_dash > r)
    {
      #if l(...) > Cr_dash + Cr
      if(avtest(estimates, seq_r, counter_r, counter_r_dash, C))
      {
        return(list(TuningParameter=seq_r[counter_r + 1], Thetahat=estimates[[counter_r + 1]]))
      }
      counter_r_dash <- counter_r_dash - 1
      r_dash <- seq_r[counter_r_dash]
    }
  }
  return(list(TuningParameter=seq_r[1], Thetahat=estimates[[1]]))
}
