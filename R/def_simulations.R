#' Returns the standard deviation of a vector
#'
#' @param x Vector
#' @result Standard deviation of vector x
#' @export
std <- function(x)
{
  return(sqrt(var(x, na.rm=TRUE)))
}

#' Concatenate mean and std to one frame
#'
#' @param mean_matrix Matrix consisting of means
#' @param std_matrix Matrix consisting of means
#' @result Matrix that consists of the mean values with corresponding standard deviations in brackets
#' @export
conc_matrix <- function(mean_matrix, std_matrix)
{
  size <- dim(mean_matrix)
  concatenated_matrix <- matrix(rep(0, size[1] * size[2]), ncol=size[2])
  for(i in 1:size[1])
  {
    for(j in 1:size[2])
    {
      concatenated_matrix[i, j] <- paste0(round(mean_matrix[i, j], 2), " (", round(std_matrix[i, j], 2), ")")
    }
  }
  colnames(concatenated_matrix) <- colnames(mean_matrix)
  rownames(concatenated_matrix) <- rownames(mean_matrix)
  return(concatenated_matrix)
}

#' Computes a simulation study of the resulting $F_1$-score, Precision and Recall using the oracle, thAV, StARS, and RIC
#' 
#' @param n Number of samples
#' @param d Number of dimensions
#' @param seq_r Set of possible regularization parameters for the thAV
#' @param C Constant C used to tune the thAV
#' @param graph Type of graph. If graph=NULL, this function does the simulation study for random and scale-free graphs, subsequentially
#' @param lambda Threshold factor used to threshold the AV. We thresholdung by lambda*C*r.
#' @param num_reps Number of repetitions of the recovery task for each recovery method
#' @param latex If TRUE, this function returns a latex output via the stargazer package
#' @return A Matrix that consists of the mean of each performance measure ($F_1$-score, Precision, and Recall) and the standard deviation in brackets
#' @export
comparisonf1 <- function(n, d, seq_r=0.2 * 0.95^(40:0), C=0.3, graph=NULL, lambda=3, num_reps=50, latex=FALSE)
{
  if( is.null(graph))
  {
    f1_results_rand <- matrix(rep(0, num_reps*4), ncol=num_reps)
    prec_results_rand <- matrix(rep(0, num_reps*4), ncol=num_reps)
    recall_results_rand <- matrix(rep(0, num_reps*4), ncol=num_reps)
    
    f1_results_sf <- matrix(rep(0, num_reps*4), ncol=num_reps)
    prec_results_sf <- matrix(rep(0, num_reps*4), ncol=num_reps)
    recall_results_sf <- matrix(rep(0, num_reps*4), ncol=num_reps)
    
    rownames(f1_results_rand) <- rownames(f1_results_sf) <- rownames(prec_results_rand) <- rownames(prec_results_sf) <- rownames(recall_results_rand) <- rownames(recall_results_sf) <- c("oracle", "thAV", "StARS", "RIC")
    colnames(f1_results_rand) <- colnames(f1_results_sf) <- colnames(prec_results_rand) <- colnames(prec_results_sf) <- colnames(recall_results_rand) <- colnames(recall_results_sf) <- 1:num_reps
    # count the times when ric estimation is NULL
    ric_calc_rand <- rep(TRUE, num_reps)
    ric_calc_sf <- rep(TRUE, num_reps)
    # random graph #
    ################
    ################
    
    pb <- txtProgressBar(min = 0, max = 2 * num_reps, style = 3)
    for(j in 1:num_reps)
    {
      theta <- generateGraph(d, graph="random")[[1]]
      
      data <- mvrnorm(n, mu=rep(0, d), Sigma=solve(theta))
      
      #oracle generation
      seq_r_oracle <- seq(0.1, 0.8, length.out = 40)
      glasso_list <- glassopath(var(data), rholist=seq_r_oracle, penalize.diagonal=FALSE, trace=0)$wi
      f1_list <- rep(0, length( seq_r_oracle))
      for(k in 1:length(seq_r_oracle))
      {
        f1_list[k] <- f1score(theta, glasso_list[, , k])$f1
      }
      oracle_estimate <- glasso_list[, , which.max(f1_list)]
      oracle_results <- f1score(theta, oracle_estimate)
      f1_results_rand[1, j] <- oracle_results$f1
      prec_results_rand[1, j] <- oracle_results$precision
      recall_results_rand[1, j] <- oracle_results$recall
      
      remove(f1_list)
      remove(glasso_list)
      remove(oracle_estimate)
      remove(oracle_results)
      
      #thAV
      thav <- thAV.estimator(data, seq_r=seq_r, C=C, lambda=lambda)
      thav_results <- f1score(theta, thav)
      f1_results_rand[2, j] <- thav_results$f1
      prec_results_rand[2, j] <- thav_results$precision
      recall_results_rand[2, j] <- thav_results$recall
      
      remove(thav)
      remove(thav_results)
      
      huge_path <- huge(data, method="glasso", verbose=FALSE)
      #stars 
      stars <- huge.select(huge_path, criterion="stars", verbose = FALSE)$opt.icov
      stars_results <- f1score(theta, stars)
      f1_results_rand[3, j] <- stars_results$f1
      prec_results_rand[3, j] <- stars_results$precision
      recall_results_rand[3, j] <- stars_results$recall
      
      remove(stars)
      remove(stars_results)
      
      #ric
      ric <- huge.select(huge_path, criterion="ric", verbose=FALSE)$opt.icov
      if(is.null(ric))
      {
        ric_calc_rand[j] <- FALSE
      }
      
      if(ric_calc_rand[j])
      {
        ric_results <- f1score(theta, ric)
        f1_results_rand[4, j] <- ric_results$f1
        prec_results_rand[4, j] <- ric_results$precision
        recall_results_rand[4, j] <- ric_results$recall
        remove(ric_results)
      }
      else
      {  
        f1_results_rand[4, j] <- NA
        prec_results_rand[4, j] <- NA
        recall_results_rand[4, j] <- NA
      }
      remove(ric)
      remove(huge_path)
      setTxtProgressBar(pb, j)
    }
    
    
    average_results_rand <- matrix(rep(0, 4*3), ncol=3)
    std_results_rand <- matrix(rep(0, 4*3), ncol=3)
    colnames(average_results_rand) <- c("$F_1$", "Precision", "Recall")
    rownames(average_results_rand) <- c("oracle", "thAV", "StARS", "RIC")
    
    average_results_rand[, 1] <- apply(f1_results_rand, 1, function(x) mean(x, na.rm=TRUE))
    average_results_rand[, 2] <- apply(prec_results_rand, 1, function(x) mean(x, na.rm=TRUE))
    average_results_rand[, 3] <- apply(recall_results_rand, 1, function(x) mean(x, na.rm=TRUE))
    
    std_results_rand[, 1] <- apply(f1_results_rand, 1, function(x) std(x))
    std_results_rand[, 2] <- apply(prec_results_rand, 1, function(x) std(x))
    std_results_rand[, 3] <- apply(recall_results_rand, 1, function(x) std(x))
    
    # scale-free graph #
    ####################
    ####################
    for(j in 1:num_reps)
    {
      network <- generateGraph(d, graph="scale-free")
      theta <- network[[1]]
      layout <- network[[2]]
      
      data <- mvrnorm(n, mu=rep(0, d), Sigma=solve(theta))
      
      #oracle generation
      seq_r_oracle <- seq(0.1, 0.8, length.out = 40)
      glasso_list <- glassopath(var(data), rholist=seq_r_oracle, penalize.diagonal=FALSE, trace=0)$wi
      f1_list <- rep(0, length( seq_r_oracle))
      for(k in 1:length(seq_r_oracle))
      {
        f1_list[k] <- f1score(theta, glasso_list[, , k])$f1
      }
      oracle_estimate <- glasso_list[, , which.max(f1_list)]
      oracle_results <- f1score(theta, oracle_estimate)
      f1_results_sf[1, j] <- oracle_results$f1
      prec_results_sf[1, j] <- oracle_results$precision
      recall_results_sf[1, j] <- oracle_results$recall
      
      remove(glasso_list)
      remove(oracle_estimate)
      remove(oracle_results)
      
      #thAV
      thav <- thAV.estimator(data, seq_r=seq_r, C=C, lambda=lambda)
      thav_results <- f1score(theta, thav)
      f1_results_sf[2, j] <- thav_results$f1
      prec_results_sf[2, j] <- thav_results$precision
      recall_results_sf[2, j] <- thav_results$recall
      
      remove(thav)
      remove(thav_results)
      
      huge_path <- huge(data, method="glasso", verbose=FALSE)
      #stars 
      stars <- huge.select(huge_path, criterion="stars", verbose = FALSE)$opt.icov
      stars_results <- f1score(theta, stars)
      f1_results_sf[3, j] <- stars_results$f1
      prec_results_sf[3, j] <- stars_results$precision
      recall_results_sf[3, j] <- stars_results$recall
      
      remove(stars)
      remove(stars_results)
      
      #ric
      ric <- huge.select(huge_path, criterion="ric", verbose=FALSE)$opt.icov
      if( is.null(ric))
      {
        ric_calc_sf[j] <- FALSE
      }
      
      if( ric_calc_sf[j])
      {
        ric_results <- f1score(theta, ric)
        f1_results_sf[4, j] <- ric_results$f1
        prec_results_sf[4, j] <- ric_results$precision
        recall_results_sf[4, j] <- ric_results$recall
        remove(ric_results)
      }
      else
      {  
        f1_results_sf[4, j] <- NA
        prec_results_sf[4, j] <- NA
        recall_results_sf[4, j] <- NA
      }
      remove(ric)
      remove(huge_path)
      setTxtProgressBar(pb, num_reps + j)
    }
    
    
    average_results_sf <- matrix(rep(0, 4*3), ncol=3)
    std_results_sf <- matrix( rep(0, 4*3), ncol=3)
    colnames(average_results_sf) <- c("$F_1$", "Precision", "Recall")
    rownames(average_results_sf) <- c("oracle", "thAV", "StARS", "RIC")
    
    average_results_sf[, 1] <- apply(f1_results_sf, 1, function(x) mean(x, na.rm=TRUE))
    average_results_sf[, 2] <- apply(prec_results_sf, 1, function(x) mean(x, na.rm=TRUE))
    average_results_sf[, 3] <- apply(recall_results_sf, 1, function(x) mean(x, na.rm=TRUE))
    
    std_results_sf[, 1] <- apply(f1_results_sf, 1, function(x) std(x))
    std_results_sf[, 2] <- apply(prec_results_sf, 1, function(x) std(x))
    std_results_sf[, 3] <- apply(recall_results_sf, 1, function(x) std(x))  
    
    close(pb)
    summary_rand <- conc_matrix(average_results_rand, std_results_rand)
    summary_sf <- conc_matrix(average_results_sf, std_results_sf)
    summary <- cbind(summary_rand, summary_sf)
    ric_stability <- rbind(ric_calc_rand, ric_calc_sf)
    rownames(ric_stability) <- c("random", "scale-free")
    print("number of times ric was computed:")
    print(apply(ric_stability, 1, function(x) sum(x)))
    if( latex)
    {
      stargazer(summary)
    }
    else
    {
      return(summary)
    }
  }
  else
  {
    f1_results <- matrix(rep(0, num_reps*4), ncol=num_reps)
    prec_results <- matrix(rep(0, num_reps*4), ncol=num_reps)
    recall_results <- matrix(rep(0, num_reps*4), ncol=num_reps)
    
    rownames(f1_results) <- rownames(prec_results) <- rownames(recall_results) <- c("oracle", "thAV", "StARS", "RIC")
    colnames(f1_results) <- colnames(prec_results) <- colnames(recall_results) <- 1:num_reps
    # count the times when ric estimation is NULL
    ric_calc <- rep(TRUE, num_reps)
    
    pb <- txtProgressBar(min = 0, max = num_reps, style = 3)
    for(j in 1:num_reps)
    {
      network <- generateGraph(d, graph=graph)
      theta <- network[[1]]
      layout <- network[[2]]
      
      data <- mvrnorm(n, mu=rep(0, d), Sigma=solve(theta))
      
      #oracle generation
      seq_r_oracle <- seq(0.1, 0.8, length.out = 40)
      glasso_list <- glassopath(var(data), rholist=seq_r_oracle, penalize.diagonal=FALSE, trace=0)$wi
      f1_list <- rep(0, length( seq_r_oracle))
      for(k in 1:length(seq_r_oracle))
      {
        f1_list[k] <- f1score(theta, glasso_list[, , k])$f1
      }
      oracle_estimate <- glasso_list[, , which.max(f1_list)]
      oracle_results <- f1score(theta, oracle_estimate)
      f1_results[1, j] <- oracle_results$f1
      prec_results[1, j] <- oracle_results$precision
      recall_results[1, j] <- oracle_results$recall
      
      remove(glasso_list)
      remove(oracle_estimate)
      remove(oracle_results)
      
      #thAV
      thav <- thAV.estimator(data, seq_r=seq_r, C=C, lambda=lambda)
      thav_results <- f1score(theta, thav)
      f1_results[2, j] <- thav_results$f1
      prec_results[2, j] <- thav_results$precision
      recall_results[2, j] <- thav_results$recall
      
      remove(thav)
      remove(thav_results)
      
      huge_path <- huge(data, method="glasso", verbose=FALSE)
      #stars 
      stars <- huge.select(huge_path, criterion="stars", verbose = FALSE)$opt.icov
      stars_results <- f1score( theta, stars)
      f1_results[3, j] <- stars_results$f1
      prec_results[3, j] <- stars_results$precision
      recall_results[3, j] <- stars_results$recall
      
      remove(stars)
      remove(stars_results)
      
      #ric
      ric <- huge.select(huge_path, criterion="ric", verbose=FALSE)$opt.icov
      if( is.null(ric))
      {
        ric_calc[j] <- FALSE
      }
      
      if( ric_calc[j])
      {
        ric_results <- f1score(theta, ric)
        f1_results[4, j] <- ric_results$f1
        prec_results[4, j] <- ric_results$precision
        recall_results[4, j] <- ric_results$recall
        remove(ric_results)
      }
      else
      {  
        f1_results[4, j] <- NA
        prec_results[4, j] <- NA
        recall_results[4, j] <- NA
      }
      remove(ric)
      remove(huge_path)
      setTxtProgressBar(pb, j)
    }
    
    average_results <- matrix( rep(0, 4*3), ncol=3)
    std_results <- matrix( rep(0, 4*3), ncol=3)
    colnames(average_results) <- c("$F_1$", "Precision", "Recall")
    rownames(average_results) <- c("oracle", "thAV", "StARS", "RIC")
    
    average_results[, 1] <- apply(f1_results, 1, function(x) mean(x, na.rm=TRUE))
    average_results[, 2] <- apply(prec_results, 1, function(x) mean(x, na.rm=TRUE))
    average_results[, 3] <- apply(recall_results, 1, function(x) mean(x, na.rm=TRUE))
    
    std_results[, 1] <- apply(f1_results, 1, function(x) std(x))
    std_results[, 2] <- apply(prec_results, 1, function(x) std(x))
    std_results[, 3] <- apply(recall_results, 1, function(x) std(x))
    
    close(pb)
    summary <- conc_matrix(average_results, std_results)
    print("number of times ric was computed:")
    print(sum(ric_calc))
    if( latex)
    {
      stargazer(summary)
    }
    else
    {
      return( summary)
    }
  }
}

#' Computes a simulation study of the resulting $F_1$-score, Precision and Recall using only thAV
#' 
#' @param n_seq Sequence of number of samples ((n_seq[j], d_seq[j]) corresponds to one setting)
#' @param d_seq Sequence of number of dimensions ((n_seq[j], d_seq[j]) corresponds to one setting)
#' @param seq_r Set of possible regularization parameters for the thAV
#' @param C Constant C used to tune the thAV
#' @param graph Type of graph. If graph=NULL, this function does the simulation study for random and scale-free graphs, subsequentially
#' @param lambda Threshold factor used to threshold the AV. We thresholdung by lambda*C*r.
#' @param num_reps Number of repetitions of the recovery task for each recovery method
#' @param latex If TRUE, this function returns a latex output via the stargazer package
#' @return A Matrix that consists of the mean of each performance measure ($F_1$-score, Precision, and Recall) and the standard deviation in brackets
#' @export
simulation_thavf1 <- function(n_seq, d_seq, seq_r=0.2 * 0.95^(40:0), C=0.3, graph="random", lambda=3, num_reps=50, latex=FALSE)
{
  num_simus <- length(n_seq)
  f1_results <- precision_results <- recall_results <- matrix( rep(0, num_reps * num_simus), ncol=num_reps)
  pb <- txtProgressBar(min = 0, max = num_reps * num_simus, style = 3)
  
  mean_results <- std_results <- matrix( rep(0, num_simus * 3), ncol=3)
  colnames(mean_results) <- colnames(std_results) <- c("F_1", "Precision", "Recall")
  
  for(simu in 1:num_simus)
  {
    n <- n_seq[simu]
    d <- d_seq[simu]
    
    for(j in 1:num_reps)
    {
      theta <- generateGraph(d, graph=graph)[[1]]
      data <- mvrnorm(n, mu=rep(0, d), Sigma=solve(theta))
      
      thav <- thAV.estimator(data, seq_r=seq_r, C=C, lambda=lambda)
      score <- f1score(theta, thav)
      f1_results[simu, j] <- score$f1
      precision_results[simu, j] <- score$precision
      recall_results[simu, j] <- score$recall
      
      setTxtProgressBar(pb, j + (simu - 1) * num_reps)
    }
    mean_results[simu, 1] <- mean(f1_results[simu, ])
    mean_results[simu, 2] <- mean(precision_results[simu, ])
    mean_results[simu, 3] <- mean(recall_results[simu, ])
    
    std_results[simu, 1] <- std(f1_results[simu, ])
    std_results[simu, 2] <- std(precision_results[simu, ])
    std_results[simu, 3] <- std(recall_results[simu, ])
  }
  
  summary <- conc_matrix(mean_results, std_results)
  if( latex)
  {
    stargazer(summary)
  }
  else
  {
    return(summary)
  }
}
#' Compares the resulting graphs for several estimation methods
#' 
#' @param n Number of samples
#' @param d Number of dimensions
#' @param type Graph type
#' @return Saves plots and $F_1$-scores of the true graph, the oracle, the thAV, the thAV2 (lambda=2), the StARS, and the RIC
#' @export
graph_comparison <- function(n, d, type="random")
{
  graph <- generateGraph(d, graph=type)
  theta <- graph[[1]]
  layout <- graph[[2]]
  
  data <- mvrnorm(n, mu=rep(0, d), Sigma=solve(theta))
  
  seq_r_oracle <- seq(0.1, 0.8, length.out = 40)
  path <- glassopath(var(data), rholist=seq_r_oracle, penalize.diagonal=FALSE, trace=0)$wi
  # find best solution
  f1_path <- rep(0, length(seq_r_oracle))
  for(k in 1:length(seq_r_oracle))
  {
    f1_path[k] <- f1score(theta, path[, , k])$f1
  }
  f1_oracle <- f1_path[which.max(f1_path)]
  
  oracle <- path[, ,which.max(f1_path)]
  thav <- thAV.estimator(data)
  huge_path <- huge(data, method="glasso", verbose=FALSE)
  stars <- huge.select(huge_path, criterion="stars", verbose=FALSE)$opt.icov
  #ric <- huge.select(huge_path, criterion="ric", verbose=FALSE)$opt.icov
  thav2 <- thAV.estimator(data, lambda=2)
  
  path_plots <- paste0("../plots/", type, "/")

  print(paste0("f1-score oracle: ", round(f1score(theta, oracle)$f1, 2)))
  print(paste0("f1-score thav: ", round(f1score(theta, thav)$f1, 2)))
  print(paste0("f1-score thav2: ", round(f1score(theta, thav2)$f1, 2)))
  print(paste0("f1-score stars: ", round(f1score(theta, stars)$f1, 2)))
'  if( is.null( ric))
  {
    print("f1-score ric: 0")
  }
  else
  {
    print(paste0("f1-score ric: ", round(f1score(theta, ric)$f1, 2)))
    ric_network <- graph.adjacency(ric, mode="undirected", weighted=TRUE, diag=FALSE)
    png(paste0(path_plots, "ric_graphn", n, "d", d, ".png"))
    par(mar=c(0, 0, 0, 0) + 0.1)
    plot.igraph( ric_network, layout=layout)
    dev.off()
  }'
  
  
  true_network <- graph.adjacency(theta, mode="undirected", weighted=TRUE, diag=FALSE)
  oracle_network <- graph.adjacency(oracle, mode="undirected", weighted=TRUE, diag=FALSE)
  thav_network <- graph.adjacency(thav, mode="undirected", weighted=TRUE, diag=FALSE)
  thav2_network <- graph.adjacency(thav2, mode="undirected", weighted=TRUE, diag=FALSE)
  stars_network <- graph.adjacency(stars, mode="undirected", weighted=TRUE, diag=FALSE)
  
  
  png(paste0(path_plots, "true_graphn", n, "d", d, ".png"))
  par(mar=c(0, 0, 0, 0) + 0.1)
  plot.igraph( true_network, layout=layout)
  dev.off()
  png(paste0(path_plots, "oracle_graphn", n, "d", d, ".png"))
  par(mar=c(0, 0, 0, 0) + 0.1)
  plot.igraph( oracle_network, layout=layout)
  dev.off()
  png(paste0(path_plots, "thav_graphn", n, "d", d, ".png"))
  par(mar=c(0, 0, 0, 0) + 0.1)
  plot.igraph( thav_network, layout=layout)
  dev.off()
  png(paste0(path_plots, "thav2_graphn", n, "d", d, ".png"))
  par(mar=c(0, 0, 0, 0) + 0.1)
  plot.igraph( thav2_network, layout=layout)
  dev.off()
  png(paste0(path_plots, "stars_graphn", n, "d", d, ".png"))
  par(mar=c(0, 0, 0, 0) + 0.1)
  plot.igraph( stars_network, layout=layout)
  dev.off()
}

#' Compare the similarity of thAV estimates based on different C.
#' 
#' @param n Number of samples
#' @param d Number of dimensions
#' @param seq_C Sequence of different C used to be compared
#' @param graph Graph type
#' @param num_reps Number of repetitions of the recovery task for each recovery method
#' @param latex If TRUE, this function returns a latex output via the stargazer package
#' @return Matrix with average (and standard deviation) of F_1-score between thAV estimates corresponding to the values in seq_C. In addition, this method also prints the average F_1-score between a thAV (based on a specific C) and the true graph.
#' @export
comparison_similarity <- function(n, d, seq_C=c(0.2, 0.3, 0.4, 0.5), graph="random", num_reps=50, latex=FALSE)
{
  pb <- txtProgressBar(min = 0, max = length(seq_C) * (length(seq_C) - 1) / 2, style = 3)
  dif_f1 <- matrix(rep(0, length(seq_C)^2 + length(seq_C)), ncol=length(seq_C))
  colnames(dif_f1) <- seq_C 
  rownames(dif_f1) <- c(seq_C, "f1")
  seq_r <- seq(0.03, 1.5, length.out = 80)
  
  for(i in 1:(length(seq_C) - 1))
  {
    for(j in (i + 1):length(seq_C))
    {
      dif_f1_over_reps <- rep(0, num_reps)
      for(rep in 1:num_reps)
      {
        theta <- generateGraph(d, graph=graph)[[1]]
        #theta <- graph[[1]]
        data <- mvrnorm(n, mu=rep(0, d), Sigma=solve(theta))
        
        thav1 <- thAV.estimator(data, C=seq_C[i], seq_r=seq_r)
        thav2 <- thAV.estimator(data, C=seq_C[j], seq_r=seq_r)
        
        score <- f1score(thav1, thav2)$f1
        if(is.nan(score))
        {
          dif_f1_over_reps[rep] <- 0
        }
        else
        {
          dif_f1_over_reps[rep] <- f1score(thav1, thav2)$f1
        }
      }
      setTxtProgressBar(pb, i + j)
      dif_f1[i, j] <- dif_f1[j, i] <- paste0(round(mean(dif_f1_over_reps), 2), " (", round(std(dif_f1_over_reps), 2), ")")
    }
  }
  remove(thav1)
  remove(thav2)
  close(pb)
  
  #calculation of f1 scores between estimation and true graph
  pb <- txtProgressBar(min = 0, max = num_reps, style = 3)
  diag( dif_f1) <- rep(1, length(seq_C))
  for(i in 1:length(seq_C))
  {
    general_f1 <- rep(0, num_reps)
    for(rep in 1:num_reps)
    {
      theta <- generateGraph(d, graph=graph)[[1]]
      data <- mvrnorm( n, mu=rep(0, d), Sigma=solve(theta))
      
      thav <- thAV.estimator(data, C=seq_C[i])
      general_f1[rep] <- f1score(theta, thav)$f1
      setTxtProgressBar(pb, rep)
    }
    dif_f1[length(seq_C) + 1, i] <- paste0(round(mean(general_f1), 2), "(", round(std(general_f1), 2), ")")
  }
  remove(thav)
  close(pb)
  
  if( latex)
  {
    stargazer( dif_f1)
  }
  else
  {
    return( dif_f1) 
  }
}


#' Compare the resulting graphs of thAV with different C
#' 
#' @param n Number of samples
#' @param d Number of dimensions
#' @param seq_C Sequence of different C used to be compared
#' @param graph Graph type
#' @param add_titel If TRUE, all plots have a title
#' @return Plots of the different thAV graphs that are computet using C from seq_C. In addition, this method prints a table with the f_1 score between the truth and each thAV estimate.
#' @export
visualization_differentC <- function(n, d, seq_C=c(0.2, 0.3, 0.4, 0.5), graph="random", add_title=TRUE)
{
  network <- generateGraph(d, graph=graph)
  theta <- network[[1]]
  layout <- network[[2]]
  
  #true graph
  network.true <- graph.adjacency(theta, mode="undirected", weighted=TRUE, diag=FALSE)
  if(add_title)
  {
    plot.igraph(network.true, layout=layout, main="True graph")
  }
  else
  {
    plot.igraph(network.true, layout=layout)
  }
  
  data <- mvrnorm(n, mu=rep(0, d), Sigma=solve(theta))
  f1_table <- matrix(rep(0, length(seq_C)), ncol=length(seq_C))
  colnames(f1_table) <- seq_C
  rownames(f1_table) <- "f1"
  for(j in 1:length(seq_C))
  {
    thav <- thAV.estimator(data, C=seq_C[j])
    f1_table[1, j] <- round(f1score( theta, thav)$f1, 2)
    
    
    network.thav <- graph.adjacency(thav, mode="undirected", weighted=TRUE, diag=FALSE)
    if(add_title)
    {
      plot.igraph(network.thav, layout=layout, main= paste0("thAV with C=", seq_C[j], ", F1=", f1_table[1, j]))
    }
    else
    {
      plot.igraph( network.thav, layout=layout)
    }
  }
  return(f1_table)
}


#' Computes the performance (F_1, Precision, and Recall) of thAV using several thresholds, and of StARS and RIC.
#' 
#' @param n Number of samples
#' @param d Number of dimensions
#' @param C Constant C used to tune the thAV
#' @param graph Type of graph. If graph=NULL, this function does the simulation study for random and scale-free graphs, subsequentially
#' @return Saves .csv-files that will be used by adaptation_threshold.ipynb to generate a plot
#' @export
adaptation_threshold <- function( n, d, C=0.3, graph="random")
{
  auxiliaries <- matrix(rep(0, 3), ncol=3)
  colnames(auxiliaries) <- c("avtp", "starsf1", "ricf1")
  thresholds <- seq(0.01, 0.2, length.out = 40)
  f1scores <- matrix(rep(0, 4 * length(thresholds)), ncol=4)
  
  theta <- generateGraph(d, graph=graph)[[1]]
  sigma <- solve(theta)
  data <- mvrnorm( n, mu=rep(0, d), Sigma=sigma)
  
  colnames(f1scores) <- c("threshold", "f1score", "precision", "recall")
  f1scores[, 1] <- thresholds
  av <- av_glasso(data, C=C)
  for(j in 1:length(thresholds))
  {
    performance <- f1score(theta, apply(av$Thetahat, c(1,2), function(x) if(abs(x)<thresholds[j]) {return(0)} else{ return(x)}))
    f1scores[j, 2] <- performance$f1
    f1scores[j, 3] <- performance$precision
    f1scores[j, 4] <- performance$recall
  }
  
  auxiliaries[1, 1] <- av$TuningParameter
  stars <- huge.select( huge(data, method="glasso", verbose=FALSE), criterion="stars", verbose=FALSE)$opt.icov
  auxiliaries[1, 2] <- f1score(theta, stars)$f1
  
  #ric
  ric <- huge.select(huge(data, method="glasso", verbose=FALSE), criterion="ric", verbose=FALSE)$opt.icov
  ric.calculated <- TRUE
  counter <- 1
  while(is.null(ric))
  {
    print(paste0("RIC failed to compute: try ", counter))
    counter <- counter + 1
    ric <- huge.select(huge(data, method="glasso", verbose=FALSE), criterion="ric", verbose=FALSE)$opt.icov
    if(counter > 5)
    {
      print("RIC failed for the 5th time: we set all the f1score to 0")
      ric.calculated <- FALSE
      break()
    }
  }
  if(ric.calculated)
  {
    auxiliaries[1, 3] <- f1score( theta, ric)$f1
  }
  else
  {
    auxiliaries[1, 3] <- 0
  }
  
  if(graph == "random")
  {
    write.csv(f1scores, paste0("data//f1thresholded_C", C, "n",n,"d", d, ".csv"),row.names = FALSE)
    write.csv(auxiliaries, paste0("data//aux_f1thresholded_C", C, "n", n, "d", d, ".csv"), row.names=FALSE)
  }
  else
  {
    write.csv(f1scores, paste0("data//sff1thresholded_C", C, "n",n,"d", d, ".csv"),row.names = FALSE)
    write.csv(auxiliaries, paste0("data//aux_sff1thresholded_C", C, "n", n, "d", d, ".csv"), row.names=FALSE)
  }
}

#' Simulation study of the recovery performance of thAV, StARS, and RIC if we change the density of a random graph
#' 
#' @param n Number of samples
#' @param d Number of dimensions
#' @param seq_density Sequence with values for p, which is the probability to connect 2 edges of a Gilbert graph (this is the way how the random graph is generated).
#' @param C Constant C used to tune the thAV
#' @param seq_r Set of possible regularization parameters for the thAV
#' @param graph Type of graph. If graph=NULL, this function does the simulation study for random and scale-free graphs, subsequentially
#' @param num_reps Number of repetitions of the recovery task for each recovery method
#' @param latex If TRUE, this function returns a latex output via the stargazer package
#' @return Matrix that measures average (and standard deviation) of the F_1-score between the true graph and an estimation (thAV, StARS, and RIC)
#' @export
comparison_density <- function( n, d, seq_density, C=0.3, seq_r=0.2 * 0.95^(40:0), graph="random", num_reps=50, latex=FALSE)
{
  pb <- txtProgressBar(min = 0, max = length(seq_density) * num_reps, style = 3)
  average_f1 <- matrix(rep(0, length(seq_density) * 3), ncol=length(seq_density))
  std_f1 <- matrix(rep(0, length(seq_density) * 3), ncol=length(seq_density))
  ric_calc <- matrix(rep(TRUE, length(seq_density) * num_reps), ncol=num_reps)
  rownames(ric_calc) <- seq_density
  colnames(ric_calc) <- 1:num_reps
  for( dens in 1:length(seq_density))
  {
    f1.thav <- rep(0, num_reps)
    f1.stars <- rep(0, num_reps)
    f1.ric <- rep(0, num_reps)
    
    for( j in 1:num_reps)
    {
      theta <- generateGraph(d, graph=graph, prob=seq_density[dens])[[1]]
      data <- mvrnorm(n, mu=rep(0, d), Sigma=solve( theta))
      thAV <- thAV.estimator(data, seq_r, C=C)
      f1.thav[j] <- f1score(theta, thAV)$f1
      remove(thAV)
      
      huge_path <- huge(data, method="glasso", verbose=FALSE)
      stars <- huge.select(huge_path, criterion="stars", verbose=FALSE)$opt.icov
      f1.stars[j] <- f1score(theta, stars)$f1
      remove(stars)
      
      ric <- huge.select(huge_path, criterion="ric", verbose=FALSE)$opt.icov
      if( is.null(ric))
      {
        ric_calc[dens, j] <- FALSE
        f1.ric[j] <- 0
      }
      else
      {
        f1.ric[j] <- f1score(theta, ric)$f1
      }
      remove(ric)
      setTxtProgressBar(pb, (dens - 1) * num_reps + j)
    }
    average_f1[1, dens] <- mean(f1.thav)
    average_f1[2, dens] <- mean(f1.stars)
    average_f1[3, dens] <- mean(f1.ric)
    
    std_f1[1, dens] <- std(f1.thav)
    std_f1[2, dens] <- std(f1.stars)
    std_f1[3, dens] <- std(f1.ric)
    
    results <- conc_matrix(round(average_f1, 2), round(std_f1, 2))
  }
  colnames(results) <- round(seq_density, 3)
  rownames(results) <- c("thAV", "stars", "ric")
  
  close(pb)
  print("Number of times RIC was computed in each setting:")
  print(apply(ric_calc, 1, function(x) sum(x)))
  if(latex)
  {
    stargazer(results)
  }
  else
  {
    return(results)
  }
}


#' Simulation study about the recovery of the num_top most significant precision matrix entries
#'
#' @param n_seq Sequence of sample sizes used for the simulation
#' @param d_seq Sequence of dimensionaly of the precision matrix used for the simulation
#' @param num_top Amount of entries that are going to be analyzed
#' @param num_reps Number of repetitions which are run during the simulation
#' @param graph Type of the graph
#' @param seq_r Sequence of possible regularization parameters for the thAV
#' @param lambda Factor used for the threshold lambda*C*r
#' @param latex If TRUE, returns a latex output
#' @return Mean and standard deviation of the proportional amount of top_edges that got recovered
#' @export
recovery_topvalues <- function(n_seq, d_seq, num_top=50, num_reps=50, graph="random", seq_r=0.2 * 0.95^(40:0), threshold=3, C=0.3, latex=FALSE)
{
  simu_avg <- matrix( rep(0, length(n_seq) * 3), ncol=length(n_seq))
  simu_std <- matrix( rep(0, length(n_seq) * 3), ncol=length(n_seq))
  
  rownames(simu_avg) <- rownames(simu_std) <- c("thAV", "StARS", "RIC")
  names_col <- rep("", length(n_seq))
  pb <- txtProgressBar(min = 0, max = length(n_seq) * num_reps, style = 3)
  for(j in 1:length(n_seq))
  {
    names_col[j] <- paste0("n=", n_seq[j], ", d=", d_seq[j])
  }
  colnames(simu_avg) <- colnames(simu_std) <- names_col
  
  for( num_setting in 1:length(n_seq))
  {
    d <- d_seq[num_setting]
    n <- n_seq[num_setting]
    score_thav <- rep(0, num_reps)
    score_stars <- rep(0, num_reps)
    score_ric <- rep(0, num_reps)
    counter_ricnull <- 0
    for(j in 1:num_reps)
    {
      theta <- generateGraph(d, graph=graph)[[1]]
      data <- mvrnorm(n, mu=rep(0, d), Sigma=solve(theta))
      
      av <- av_glasso(data, seq_r, C)
      av_tp <- av$TuningParameter
      av_estimate <- av$Thetahat
      
      thav <- threshold(av, C, threshold)
      
      huge_path <- huge(data, method="glasso", verbose=FALSE)
      stars <- huge.select(huge_path, criterion="stars", verbose=FALSE)$opt.icov
      ric <- huge.select(huge_path, criterion="ric", verbose=FALSE)$opt.icov
      
      top_values <- sort( abs(theta[ upper.tri(theta)]), decreasing=TRUE )
      top_coords <- matrix( rep(0, 2*num_top), ncol=2)
      top_thav <- rep(0, num_top)
      top_stars <- rep(0, num_top)
      top_ric <- rep(0, num_top)
      counter_thav <- 0
      counter_stars <- 0
      counter_ric <- 0
      for( k in 1:num_top)
      {
        top_coords[k, ] <- which( matrix(abs(theta), ncol=d)==top_values[k], arr.ind=TRUE)[1,]
        top_thav[k] <- thav[top_coords[k, 1], top_coords[k, 2]]
        top_stars[k] <- stars[top_coords[k, 1], top_coords[k, 2]]
        if( top_stars[k] != 0)
        {
          counter_stars <- counter_stars + 1
        }
        if( !is.null(ric))
        {
          top_ric[k] <- ric[top_coords[k, 1], top_coords[k, 2]]
          if( top_ric[k] != 0)
          {
            counter_ric <- counter_ric + 1
          }
        }
        else
        {
          score_ric[j] <- NA
          counter_ricnull <- counter_ricnull + 1/num_top # we add 1/num_top because we do this                                                             # num_top times: once for each top_value
        }
        if( top_thav[k] != 0)
        {
          counter_thav <- counter_thav + 1
        }
        
      }
      score_thav[j] <- counter_thav / num_top
      score_stars[j] <- counter_stars / num_top
      if( !is.null(score_ric[j]))
      {
        score_ric[j] <- counter_ric / num_top
      }
      
      setTxtProgressBar(pb, (num_setting - 1) * num_reps + j)
    }
    print(paste0("ric failed in setting n=", n, ", d=", d, " ", counter_ricnull ," times."))
    print(score_thav)
    simu_avg[1, num_setting] <- round(mean(score_thav), 2)
    simu_std[1, num_setting] <- round(std(score_thav), 2)
    
    simu_avg[2, num_setting] <- round(mean(score_stars), 2)
    simu_std[2, num_setting] <- round(std(score_stars), 2)
    
    simu_avg[3, num_setting] <- round(mean(score_ric, na.rm=TRUE), 2)
    simu_std[3, num_setting] <- round(std(score_ric), 2)
    
    remove(huge_path)
    remove(theta)
    remove(data)
    remove(ric)
    remove(stars)
    remove(thav)
    remove(av)
    remove(av_estimate)
  }
  summary <- conc_matrix(simu_avg, simu_std)
  close(pb)
  if( latex)
  {
    stargazer(summary)
  }
  else
  {
    return(summary)
  }
}

#' Exemplary recovery in F_1-score of the top edges for several settings
#'
#' @param n_seq Sequence of sample sizes used for the simulation
#' @param d_seq Sequence of dimensionaly of the precision matrix used for the simulation
#' @param num_top Amount of entries that are going to be analyzed
#' @param graph Type of the graph
#' @param latex If TRUE, returns a latex output
#' @return Table including the largest, the num_top'th largest absolute value of Theta, the reference value 6Cr, and the proportional amount of recovered top edges for several settings
#' @export 
thav_topvalues <- function(n_seq, d_seq, num_top=50, graph="random", latex=TRUE)
{
  summary <- matrix( rep(0, 4 * length(n_seq)), ncol=length(n_seq))
  rownames(summary) <- c("top1", paste0("top",num_top), "6Cr", "occurence%")
  for(num_setting in 1:length(n_seq))
  {
    d <- d_seq[num_setting]
    n <- n_seq[num_setting]
    
    theta <- generateGraph(d, graph=graph)[[1]]
    data <- mvrnorm(n, mu=rep(0, d), Sigma=solve(theta))
    
    av <- av_glasso(data)
    av_tp <- av$TuningParameter
    thav <- threshold(av)
    
    top_values <- sort( abs(theta[ upper.tri(theta)]), decreasing=TRUE)
    top_1 <- top_values[1]
    top_num_top <- top_values[num_top]
    top_coords <- matrix(rep(0, 2*50), ncol=2)
    counter_thav <- 0
    for(j in 1:50)
    {
      top_coords[j, ] <- which( matrix(abs(theta), ncol=d)==top_values[j], arr.ind=TRUE)[1, ]
      if( thav[ top_coords[j, 1], top_coords[j, 2]] != 0)
      {
        counter_thav <- counter_thav + 1
      }
    }
    summary[, num_setting] <- c(top_values[1], top_values[num_top], 6 * 0.3 * av_tp, counter_thav / num_top)
  }
  stargazer( round(summary, 2))
}

#' Compute the performance of the thAV for large-scale examples
#' 
#' @param d Number of dimensions
#' @param n Number of samples
#' @param seq_r Set of possible regularization parameters for the thAV
#' @param C Constant C used to tune the thAV
#' @param graph Type of graph. If graph=NULL, this function does the simulation study for random and scale-free graphs, subsequentially
#' @param lambda Threshold factor used to threshold the AV. We thresholdung by lambda*C*r.
#' @param num_reps Number of repetitions of the recovery task for each recovery method
#' @param latex If TRUE, this function returns a latex output via the stargazer package
#' @return A Matrix that consists of the mean of each performance measure ($F_1$-score, Precision, and Recall) and the standard deviation in brackets.
#' @export
performance_thav <- function(d, n=500, seq_r=0.2 * 0.95^(40:0), C=0.3, graph="random", lambda=3, num_reps=20, latex=TRUE)
{
  performance <- matrix(rep(0, 3 * num_reps), ncol=3)
  colnames(performance) <- c("F_1", "Precision", "Recall")
  pb <- txtProgressBar(min=0, max=num_reps ,style=3)
  for(j in 1:num_reps)
  {
    theta <- generateGraph(d, graph=graph)[[1]]
    data <- mvrnorm(n, mu=rep(0, d), Sigma=solve(theta))
    
    thav <- thAV.estimator(data, seq_r=seq_r, C=C, lambda=lambda)
    this_performance <- f1score(theta, thav)
    performance[j, ] <- c(this_performance$f1, this_performance$precision, this_performance$recall)
    setTxtProgressBar(pb, j)
  }
  avg <- round(apply(performance, 2, function(x) mean(x)), 2)
  std <- round(apply(performance, 2, function(x) std(x)), 2)
  summary <- matrix(paste0( avg, " (", std, ")"),ncol=3)
  colnames(summary) <- colnames(performance)
  if(latex)
  {
    stargazer(summary)
  }
  else
  {
    print(summary)
  }
}

#' Investigate the empirical distribution of the absolute values of the graphs edges
#'
#' @param d Number of dimensions
#' @param graph Graph type
#' @param num_reps Number of repetitions
#' @return Saves a .csv-file that can be read by visualizations_weights.ipynb to plot a histogram of the weight distribution
#' @export
weight_dist <- function(d, graph="random", num_reps=50)
{
  weights <- c()
  for( rep in 1:num_reps)
  {
    theta <- generateGraph(d, graph=graph)[[1]]
    
    upper_weights <- as.vector(theta[upper.tri(theta)])
    weights <- c(weights, abs(upper_weights[which(upper_weights != 0)]))
  }
  file_name <- paste0("..\\data\\weights_", graph, "d", d, ".csv")
  write.csv(weights, file_name, row.names=FALSE)
}
