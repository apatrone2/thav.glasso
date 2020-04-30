#' Compute the Adaptive Validation Graphical Lasso
#' 
#' @param data Data matrix of size $n\times d$ 
#' @param seq_r Sequence of considered regularization parameters
#' @param C constant C that is used to calibrate the AV regularization parameter
#' @return The selected regularization parameter
#' @return The resulting graphical lasso estimate
#' @export
av_glasso <- function(data, seq_r= 0.2 * 0.95^(40:0), C=0.3)
{
  #initialization
  r <- max(seq_r)
  estimates <- list()
  empvar <- var(data)
  
  #calculate glasso[max(r)] and setup the warm start
  glasso <- glasso(empvar, rho=r, penalize.diagonal=FALSE)
  warm <- list(glasso$wi, glasso$w)
  estimates[[length(seq_r)]] <- warm[[1]]
  
  counter_r <- length(seq_r) #used to fill the list estimates with corresponding estimates
  while(r != min(seq_r))
  {
    #compute glasso[r] using a warm start
    counter_r <- counter_r - 1
    newwarm <- glasso(empvar, rho=seq_r[ counter_r], penalize.diagonal=FALSE, start="warm", wi.init=warm[[1]], w.init=warm[[2]])
    estimates[[counter_r]] <- newwarm$wi
    warm <- list(newwarm$wi, newwarm$w)
    
    counter_r_dash <- length(seq_r)
    r_dash <- seq_r[counter_r_dash]
    while(r_dash > r)
    {
      #if l(...) > Cr_dash + Cr
      if(avtest(estimates, seq_r, counter_r, counter_r_dash, C))
      {
        return(list(TuningParameter=seq_r[counter_r + 1], Thetahat=estimates[[counter_r]]))
      }
      counter_r_dash <- counter_r_dash - 1
      r_dash <- seq_r[counter_r_dash]
    }
    #counter_r <- counter_r - 1
    r <- seq_r[counter_r - 1]
  }
  return(list(TuningParameter=r, Thetahat=estimates[[counter_r]]))
}

#' Computes the test that is used in the AV calibration scheme
#' 
#' @param estimates List of the graphical lasso estimators so far
#' @param seq_r List of possible regularization paraeters
#' @param counter_r position of the current r in r_seq (see algorithm)
#' @param counter_r_dash position of the current r_dash in seq_r (see algorithm)
#' @param C constant C that is used for the AV calibration
#' @return TRUE if $l(...) > Cr_dash + Cr$, else FALSE
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
threshold <- function( av, C=0.3, lambda=3)
{
  return( apply( av$Thetahat, c(1,2), function(x) if( abs(x) <= lambda*C*av$TuningParameter){ return(0)} else{ return(x)}))
}

#' Compute the thAV estimator
#' 
#' @param data Data matrix of size $n\times d$ 
#' @param seq_r Sequence of considered regularization parameters
#' @param C constant C that is used to calibrate the AV regularization parameter
#' @param lambda Threshold factor lambda: We threshold the AV solution by lambda*C*av$TuningParameter
#' @return The thresholded AV estimator
#' @export
thAV.estimator <- function( data, seq_r=0.2 * 0.95^(40:0), C=0.3, lambda=3)
{
  av <- av_glasso(data, seq_r, C)
  return( threshold(av, C, lambda))
}
