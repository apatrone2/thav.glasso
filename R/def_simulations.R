#' Returns the standard deviation of a vector
#'
#' @param x Vector
#' @return Standard deviation of vector x
#' @export
std <- function(x)
{
  return(sqrt(var(x, na.rm=TRUE)))
}

#' Concatenate mean and std to one frame
#'
#' @param mean_matrix Matrix consisting of means
#' @param std_matrix Matrix consisting of standard deviations
#' @return Matrix (of strings) that consists of the mean values with corresponding standard deviations in brackets
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

#' Computes parts of the simulation study that compares F1-score, Precision and Recall using the oracle, StARS, scaled Lasso, TIGER, rSME using eBIC with gamma=1, rSME using eBIC with gamma=0.5, SCIO using Cross-Validation, SCIO using the bregman loss, thAV-glasso, AV-rSME, and AV-SCIO.
#' We split the simulation into several parts because of memory issues.
#' 
#' @param n Number of samples
#' @param d Number of dimensions
#' @param part Number that describes the part of the simulation study
#' @param num_reps Number of repetitions (in this part) of the recovery task for each recovery method
#' @param graph Type of graph. 
#' @param C_thav Constant C used to tune the thAV
#' @param C_score Constant C used to tune the AV-rSME
#' @param C_scio Constant C used to tune the AV-SCIO
#' @return Writes a .txt-file that can be read by the simu_f1 function to summarize the simulation study. 
#' @export
simu_f1_split <- function(n, d, part, num_reps=10, graph="random", C_thav=0.7)
{
  f1 <- precision <- recall <- matrix(rep(0, num_reps*9), ncol=num_reps)
  
  pb <- txtProgressBar(min = 0, max = num_reps * 9, style = 3)
  counter_progress <- 0
  for(j in 1:num_reps)
  {
    theta <- generateGraph(d, graph=graph)[[1]]
    data <- scale(mvrnorm(n, mu=rep(0, d), Sigma=solve(theta)))
    
    
    #oracle generation
    seq_r_oracle <- seq(0.1, 0.8, length.out = 40)
    glasso_list <- glassopath(var(data), rholist=seq_r_oracle, penalize.diagonal=FALSE, trace=0)$wi
    f1_list <- rep(0, length( seq_r_oracle))
    for(k in 1:length(seq_r_oracle))
    {
      f1_list[k] <- f1score(theta, glasso_list[, , k])$f1
    }
    oracle_estimate <- glasso_list[, , which.max(f1_list)]
    counter_progress <- counter_progress + 1
    setTxtProgressBar(pb, counter_progress)
    
    thav <- thAV.estimator(data, C=C_thav)
    counter_progress <- counter_progress + 1
    setTxtProgressBar(pb, counter_progress)
    
    score_ebic1 <- score_ebic(data, gamma=1)
    counter_progress <- counter_progress + 1
    setTxtProgressBar(pb, counter_progress)
    
    score_ebic2 <- score_ebic(data, gamma=0.5)
    counter_progress <- counter_progress + 1
    setTxtProgressBar(pb, counter_progress)
    
    scio_cv <- scio.cv(data)$w
    counter_progress <- counter_progress + 1
    setTxtProgressBar(pb, counter_progress)
    
    scio_bregman_est <- scio_bregman(data)
    counter_progress <- counter_progress + 1
    setTxtProgressBar(pb, counter_progress)
    
    stars <- huge.select(huge(data, method="glasso", verbose=FALSE), criterion="stars", verbose=FALSE)$opt.icov
    counter_progress <- counter_progress + 1
    setTxtProgressBar(pb, counter_progress)
    
    scaled_lasso <- scalreg(data)$precision
    counter_progress <- counter_progress + 1
    setTxtProgressBar(pb, counter_progress)
    
    tiger <- huge.tiger(data, verbose=FALSE)
    tiger_tp <- sqrt(log(d) / n)
    tiger_est <- tiger$icov[[which.min(abs(tiger$lambda - tiger_tp))]]
    counter_progress <- counter_progress + 1
    setTxtProgressBar(pb, counter_progress)
    
    performance_oracle <- f1score(theta, oracle_estimate)
    performance_thav <- f1score(theta, thav)
    performance_score_ebic1 <- f1score(theta, score_ebic1)
    performance_score_ebic2 <- f1score(theta, score_ebic2)
    performance_scio_cv <- f1score(theta, scio_cv)
    performance_scio_bregman <- f1score(theta, scio_bregman_est)
    performance_stars <- f1score(theta, stars)
    performance_scaled_lasso <- f1score(theta, scaled_lasso)
    performance_tiger <- f1score(theta, tiger_est)
    
    f1[1, j] <- performance_oracle$f1 
    f1[2, j] <- performance_stars$f1 
    f1[3, j] <- performance_scaled_lasso$f1
    f1[4, j] <- performance_tiger$f1
    f1[5, j] <- performance_score_ebic1$f1 
    f1[6, j] <- performance_score_ebic2$f1 
    f1[7, j] <- performance_scio_cv$f1 
    f1[8, j] <- performance_scio_bregman$f1
    f1[9, j] <- performance_thav$f1
    
    
    precision[1, j] <- performance_oracle$precision 
    precision[2, j] <- performance_stars$precision
    precision[3, j] <- performance_scaled_lasso$precision
    precision[4, j] <- performance_tiger$precision
    precision[5, j] <- performance_score_ebic1$precision
    precision[6, j] <- performance_score_ebic2$precision
    precision[7, j] <- performance_scio_cv$precision
    precision[8, j] <- performance_scio_bregman$precision
    precision[9, j] <- performance_thav$precision
    
    recall[1, j] <- performance_oracle$recall 
    recall[2, j] <- performance_stars$recall
    recall[3, j] <- performance_scaled_lasso$recall
    recall[4, j] <- performance_tiger$recall
    recall[5, j] <- performance_score_ebic1$recall 
    recall[6, j] <- performance_score_ebic2$recall 
    recall[7, j] <- performance_scio_cv$recall
    recall[8, j] <- performance_scio_bregman$recall
    recall[9, j] <- performance_thav$recall
  }
  
  PATH_f1 <- paste0("../data/simu_results/comparison_f1/f1n", n, "d", d , graph, part,".txt")
  PATH_prec <- paste0("../data/simu_results/comparison_f1/precn", n, "d", d , graph, part,".txt")
  PATH_recall <- paste0("../data/simu_results/comparison_f1/recalln", n, "d", d , graph, part,".txt")
  write.csv(f1, PATH_f1, row.names=FALSE)
  write.csv(precision, PATH_prec, row.names=FALSE)
  write.csv(recall, PATH_recall, row.names=FALSE)
  
  close(pb)
}

#' Reads the .txt files created by the simu_f1_split function to summarize the resulting F1-score, Precision and Recall of the oracle, StARS, scaled Lasso, TIGER, rSME using eBIC with gamma=1, rSME using eBIC with gamma=0.5, SCIO using Cross-Validation, SCIO using the bregman loss, thAV-glasso, AV-rSME, and AV-SCIO. It is necessary to run simu_f1_spit beforehand.
#' 
#' @param n Number of samples
#' @param d Number of dimensions
#' @param graph Type of graph, "random" or "scale-free"
#' @param num_parts Number of parts in which we split the simulation study
#' @return A tex-ready code for a table that consists of the mean of each performance measure (F1-score, Precision, and Recall) and the standard deviation in brackets. The result is saved as a .txt file.
#' @export
simu_f1 <- function(n, d, graph="random", num_parts=3)
{
  for(part in 1:num_parts)
  {
    PATH_f1 <- paste0("../data/simu_results/comparison_f1/f1n", n, "d", d, graph, part,".txt")
    PATH_prec <- paste0("../data/simu_results/comparison_f1/precn", n, "d", d, graph, part, ".txt")
    PATH_recall <- paste0("../data/simu_results/comparison_f1/recalln", n, "d", d, graph, part, ".txt")
    
    readf1_results <- as.matrix(read.csv(PATH_f1), header=TRUE)[1:9,]
    readprec_results <- as.matrix(read.csv(PATH_prec), header=TRUE)[1:9, ]
    readrecall_results <- as.matrix(read.csv(PATH_recall), header=TRUE)[1:9, ]
    
    if(part==1)
    {
      f1_results <- readf1_results
      prec_results <- readprec_results
      recall_results <- readrecall_results
    }
    else
    {
      f1_results <- cbind(f1_results, readf1_results)
      prec_results <- cbind(prec_results, readprec_results)
      recall_results <- cbind(recall_results, readrecall_results)
    }
  }
  
  f1_results <- f1_results[, 1:9]
  prec_results <- prec_results[, 1:9]
  recall_results <- recall_results[, 1:9]
  
  rownames(f1_results) <- rownames(prec_results) <- rownames(recall_results) <- c("oracle", "StARS", "scaled lasso", "TIGER", "score_ebic1", "score_ebic0.5", "scio_cv", "scio_bregman", "thAV")
  
  average_results <- matrix( rep(0, 9*3), ncol=3)
  std_results <- matrix( rep(0, 9*3), ncol=3)
  colnames(average_results) <- c("F1", "Precision", "Recall")
  rownames(average_results) <- c("oracle", "StARS", "scaled lasso", "TIGER", "score_ebic1", "score_ebic0.5", "scio_cv", "scio_bregman", "thAV")
  
  average_results[, 1] <- apply(f1_results, 1, function(x) mean(x, na.rm=TRUE))
  average_results[, 2] <- apply(prec_results, 1, function(x) mean(x, na.rm=TRUE))
  average_results[, 3] <- apply(recall_results, 1, function(x) mean(x, na.rm=TRUE))
  
  std_results[, 1] <- apply(f1_results, 1, function(x) std(x))
  std_results[, 2] <- apply(prec_results, 1, function(x) std(x))
  std_results[, 3] <- apply(recall_results, 1, function(x) std(x))

  summary <- conc_matrix(average_results, std_results)

  PATH = paste0("../data/simu_results/comparison_f1/n", n, "d", d , graph, ".txt")
  file_results <- file(PATH)
  writeLines(stargazer(summary), file_results)
  close(file_results)
}


#' Runs a simulation study of the resulting F1-score, Precision and Recall using only thAV. 
#' 
#' @param n_seq Sequence of number of samples ((n_seq[j], d_seq[j]) corresponds to one setting)
#' @param d_seq Sequence of number of dimensions ((n_seq[j], d_seq[j]) corresponds to one setting)
#' @param estimator Chosen estimator. Must be either "thav_glasso", "thav_rsme", or "av_scio"
#' @param graph Type of graph
#' @param num_reps Number of repetitions of the recovery task for each recovery method
#' @param latex If TRUE, this function returns a latex output using the stargazer package
#' @return A tex-ready code for a table that consists of the mean of each performance measure (F1-score, Precision, and Recall) and the standard deviation in brackets. The result is saved as a .txt file.
#' @export
simu_av_f1 <- function(n_seq, d_seq, estimator="thav_glasso", graph="random", num_reps=10, latex=TRUE)
{
  num_simus <- length(n_seq)
  f1_results <- precision_results <- recall_results <- matrix( rep(0, num_reps * num_simus), ncol=num_reps)
  pb <- txtProgressBar(min = 0, max = num_reps * num_simus, style = 3)
  
  mean_results <- std_results <- matrix( rep(0, num_simus * 3), ncol=3)
  colnames(mean_results) <- colnames(std_results) <- c("F1", "Precision", "Recall")
  
  for(simu in 1:num_simus)
  {
    n <- n_seq[simu]
    d <- d_seq[simu]
    
    for(j in 1:num_reps)
    {
      theta <- generateGraph(d, graph=graph)[[1]]
      data <- scale(mvrnorm(n, mu=rep(0, d), Sigma=solve(theta)))
      
      if(estimator == "thav_glasso")
      {
          est <- thAV.estimator(data)
      }
      else if(estimator=="thav_rsme")
      {
          est <- av_rsme(data, C=0.5)
          est_th <- threshold(est, C=0.5, lambda=3)
      }
      else if(estimator=="av_scio")
      {
          est <- av_scio(data, C=0.7)$Thetahat
      }
      
      score <- f1score(theta, est)
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
    print(paste0("Employed Estimator: ", estimator))
    stargazer(summary)
  }
  else
  {
    print(paste0("Employed Estimator: ", estimator))
    return(summary)
  }
}

#' Compares the resulting graphs for several estimation methods (true graph, oracle, StARS, scaled lasso, TIGER, rSME (eBIC), SCIO (CV), SCIO (bregman), and thAV)
#' 
#' @param n Number of samples
#' @param d Number of dimensions
#' @param type Graph type
#' @return Saves plots of the graphs and returns F1-scores
#' @export
graph_comparison <- function(n, d, type="random")
{
  graph <- generateGraph(d, graph=type)
  theta <- graph[[1]]
  layout <- graph[[2]]
  
  data <- scale(mvrnorm(n, mu=rep(0, d), Sigma=solve(theta)))
  path_plots <- paste0("../plots/", type, "/n", n, "d", d, "/")
  
  ### True graph:
  true_network <- graph.adjacency(theta, mode="undirected", weighted=TRUE, diag=FALSE)
  pdf(paste0(path_plots, "true_graphn", n, "d", d, ".pdf"))#, width=1000, height=800, res=200)
  par(mar=c(0, 0, 0, 0) + 0.1)
  plot.igraph( true_network, layout=layout)
  dev.off()
  remove(true_network)
  
  ### oracle estimator:
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
  print(paste0("f1-score oracle: ", round(f1score(theta, oracle)$f1, 2)))
  oracle_network <- graph.adjacency(oracle, mode="undirected", weighted=TRUE, diag=FALSE)
  pdf(paste0(path_plots, "oracle_graphn", n, "d", d, ".pdf"))#, width=1000, height=800, res=200)
  par(mar=c(0, 0, 0, 0) + 0.1)
  plot.igraph( oracle_network, layout=layout)
  dev.off()
  remove(path)
  remove(f1_oracle)
  remove(oracle_network)
  
  
  ### thav:
  thav <- thAV.estimator(data)
  print(paste0("f1-score thav: ", round(f1score(theta, thav)$f1, 2)))
  thav_network <- graph.adjacency(thav, mode="undirected", weighted=TRUE, diag=FALSE)
  pdf(paste0(path_plots, "thav_graphn", n, "d", d, ".pdf"))#, width=1000, height=800, res=200)
  par(mar=c(0, 0, 0, 0) + 0.1)
  plot.igraph( thav_network, layout=layout)
  dev.off()
  remove(thav)
  remove(thav_network)
  
  ### stars:
  huge_path <- huge(data, method="glasso", verbose=FALSE)
  stars <- huge.select(huge_path, criterion="stars", verbose=FALSE)$opt.icov
  print(paste0("f1-score stars: ", round(f1score(theta, stars)$f1, 2)))
  stars_network <- graph.adjacency(stars, mode="undirected", weighted=TRUE, diag=FALSE)
  pdf(paste0(path_plots, "stars_graphn", n, "d", d, ".pdf"))#, width=1000, height=800, res=200)
  par(mar=c(0, 0, 0, 0) + 0.1)
  plot.igraph( stars_network, layout=layout)
  dev.off()
  remove(huge_path)
  remove(stars)
  remove(stars_network)
  
  ### scaled lasso:
  scaled_lasso <- scalreg(data)$precision
  print(paste0("f1-score scaled lasso:", round(f1score(theta, scaled_lasso)$f1, 2)))
  scaled_lasso_network <- graph.adjacency(scaled_lasso, mode="undirected", weighted=TRUE, diag=FALSE)
  pdf(paste0(path_plots, "scaledlasso_graphn", n, "d", d, ".pdf"))#, width=1000, height=800, res=200)
  par(mar=c(0, 0, 0, 0) + 0.1)
  plot.igraph( scaled_lasso_network, layout=layout)
  dev.off()
  remove(scaled_lasso_network)
  remove(scaled_lasso)
  
  ### tiger:
  tiger <- huge.tiger(data, verbose=FALSE)
  tiger_tp <- sqrt(2) * sqrt(log(d)/(2*n))
  tiger_est <- tiger$icov[[which.min(abs(tiger$lambda - tiger_tp))]]
  print(paste0("f1-score tiger:", round(f1score(theta, tiger_est)$f1, 2)))
  tiger_network <- graph.adjacency(tiger_est, mode="undirected", weight=TRUE, diag=FALSE)
  pdf(paste0(path_plots, "tiger_graphn", n, "d", d, ".pdf"))#, width=1000, height=800, res=200)
  par(mar=c(0, 0, 0, 0) + 0.1)
  plot.igraph( tiger_network, layout=layout)
  dev.off()
  remove(tiger)
  remove(tiger_network)
  remove(tiger_est)
  
  
  ### scio(cv):
  scio_cv <- scio.cv(data)$w
  print(paste0("f1-score scio (CV):", round(f1score(theta, scio_cv)$f1, 2)))
  scio_cv_network <- graph.adjacency(scio_cv, mode="undirected", weight=TRUE, diag=FALSE)
  pdf(paste0(path_plots, "sciocv_graphn", n, "d", d, ".pdf"))#, width=1000, height=800, res=200)
  par(mar=c(0, 0, 0, 0) + 0.1)
  plot.igraph(scio_cv_network, layout=layout)
  dev.off()
  remove(scio_cv)
  remove(scio_cv_network)
  
  ### scio (bregman):
  scio_bregman_est <- scio_bregman(data)
  print(paste0("f1-score scio (Bregman):", round(f1score(theta, scio_bregman_est)$f1, 2)))
  scio_bregman_network <- graph.adjacency(scio_bregman_est, mode="undirected", weight=TRUE, diag=FALSE)
  pdf(paste0(path_plots, "sciobregman_graphn", n, "d", d, ".pdf"))#, width=1000, height=800, res=200)
  par(mar=c(0, 0, 0, 0) + 0.1)
  plot.igraph( scio_bregman_network, layout=layout)
  dev.off()
  remove(scio_bregman_est)
  remove(scio_bregman_network)
  
  ### rSME:
  score_ebic1 <- score_ebic(data, gamma=1)
  print(paste0("f1-score score (eBIC):", round(f1score(theta, score_ebic1)$f1, 2)))
  score_network <- graph.adjacency(score_ebic1, mode="undirected", weight=TRUE, diag=FALSE)
  pdf(paste0(path_plots, "score_graphn", n, "d", d, ".pdf"))#, width=1000, height=800, res=200)
  par(mar=c(0, 0, 0, 0) + 0.1)
  plot.igraph( score_network, layout=layout)
  dev.off()
  remove(score_ebic1)
  remove(score_network)
  
  ### thav (rSME):
  thav_rsme <- threshold(av_rsme(data, C=0.5), C=0.5, lambda=2)
  print(paste0("f1-score thav (rSME): ", round(f1score(theta, thav_rsme)$f1, 2)))
  rsme_network <- graph.adjacency(thav_rsme, mode="undirected", weighted=TRUE, diag=FALSE)
  pdf(paste0(path_plots, "thav_rsme_graphn", n, "d", d, ".pdf"))#, width=1000, height=800, res=200)
  par(mar=c(0, 0, 0, 0) + 0.1)
  plot.igraph( rsme_network, layout=layout)
  dev.off()
  remove(thav_rsme)
  remove(rsme_network)
  
  ### thav (sf):
  thav_sf <- threshold(av_sf_glasso(data, C=1.5), C=1.5, lambda=1.5)
  print(paste0("f1-score thav (sf-glasso): ", round(f1score(theta, thav_sf)$f1, 2)))
  sf_network <- graph.adjacency(thav_sf, mode="undirected", weighted=TRUE, diag=FALSE)
  pdf(paste0(path_plots, "thav_sf_graphn", n, "d", d, ".pdf"))#, width=1000, height=800, res=200)
  par(mar=c(0, 0, 0, 0) + 0.1)
  plot.igraph( sf_network, layout=layout)
  dev.off()
}

#' Compares the similarity of thAV estimates based on different C.
#' 
#' @param n Number of samples
#' @param d Number of dimensions
#' @param seq_C Sequence of different C used to be compared
#' @param graph Graph type
#' @param num_reps Number of repetitions of the recovery task for each recovery method
#' @param latex If TRUE, this function returns a latex output via the stargazer package
#' @return Matrix with average (and standard deviation) of F1-score between thAV estimates corresponding to the values in seq_C. In addition, this method also prints the average F1-score between a thAV (based on a specific C) and the true graph. The results are saved in a .txt file.
#' @export
comparison_similarity <- function(n, d, seq_C=c(0.5, 0.6, 0.7, 0.8), graph="random", num_reps=25, latex=FALSE)
{
  pb <- txtProgressBar(min = 0, max = length(seq_C) * (length(seq_C) - 1) / 2 + num_reps, style = 3)
  dif_f1 <- matrix(rep(0, length(seq_C)^2 + length(seq_C)), ncol=length(seq_C))
  colnames(dif_f1) <- seq_C 
  rownames(dif_f1) <- c(seq_C, "f1")
  
  for(i in 1:(length(seq_C) - 1))
  {
    for(j in (i + 1):length(seq_C))
    {
      dif_f1_over_reps <- rep(0, num_reps)
      for(rep in 1:num_reps)
      {
        theta <- generateGraph(d, graph=graph)[[1]]
        data <- scale(mvrnorm(n, mu=rep(0, d), Sigma=solve(theta)))
        
        thav1 <- thAV.estimator(data, C=seq_C[i], lambda=1)
        thav2 <- thAV.estimator(data, C=seq_C[j], lambda=1)
        
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
  
  #calculation of f1 scores between estimation and true graph
  diag( dif_f1) <- rep(1, length(seq_C))
  for(i in 1:length(seq_C))
  {
    general_f1 <- rep(0, num_reps)
    for(rep in 1:num_reps)
    {
      theta <- generateGraph(d, graph=graph)[[1]]
      data <- scale(mvrnorm( n, mu=rep(0, d), Sigma=solve(theta)))
      
      thav <- thAV.estimator(data, C=seq_C[i], lambda=1)
      general_f1[rep] <- f1score(theta, thav)$f1
      setTxtProgressBar(pb, length(seq_C) * (length(seq_C) - 1) / 2 + rep)
    }
    dif_f1[length(seq_C) + 1, i] <- paste0(round(mean(general_f1), 2), "(", round(std(general_f1), 2), ")")
  }
  remove(thav)
  close(pb)
  
  if( latex)
  {
    PATH = paste0("../data/simu_results/comparison_similarity/n", n, "d", d , graph, ".txt")
    file_results <- file(PATH)
    writeLines(stargazer(dif_f1), file_results)
    print(paste0("The results were saved in ", PATH))
    close(file_results)
  }
  else
  {
    return(dif_f1) 
  }
}


#' Compares the resulting graphs of thAV with different C
#' 
#' @param n Number of samples
#' @param d Number of dimensions
#' @param seq_C Sequence of different C used to be compared
#' @param graph Graph type
#' @param add_titlw If TRUE, all plots have a title
#' @return Plots of the different thAV graphs that are computet using C from seq_C. In addition, this method prints a table with the F1 score between the truth and each thAV estimate.
#' @export
visualization_differentC <- function(n, d, seq_C=c(0.5, 0.6, 0.7, 0.8), graph="random", add_title=FALSE)
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
  
  data <- scale(mvrnorm(n, mu=rep(0, d), Sigma=solve(theta)))
  
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
      pdf(paste0("../plots/", graph, "/n", n, "d",d, "/different_C/C", seq_C[j], ".pdf"))#, width=1000, height=1000, res=200)
      par(mar=c(0, 0, 0, 0) + 0.1)
      plot.igraph( network.thav, layout=layout)
      dev.off()
    }
  }
  return(f1_table)
}


#' Computes the performance (F1, Precision, and Recall) of thAV using different thresholds. The generated .txt-files are used for the adaptation_threshold.ipynb to create a plot.
#' 
#' @param n Number of samples
#' @param d Number of dimensions
#' @param C Constant C used to tune the thAV
#' @param graph Type of graph
#' @return Saves .txt-files that will be used by adaptation_threshold.ipynb to generate a plot
#' @export
adaptation_threshold <- function(n, d, seq_C=c(0.5, 0.6, 0.7, 0.8), graph="random")
{
  theta <- generateGraph(d, graph=graph)[[1]]
  sigma <- solve(theta)
  data <- scale(mvrnorm( n, mu=rep(0, d), Sigma=sigma))
  for(C in seq_C)
  {
    auxiliaries <- matrix(rep(0, 5), ncol=5)
    colnames(auxiliaries) <- c("avtp", "starsf1", "ricf1", "scaledlassof1", "tigerf1")
    thresholds <- seq(0.01, 0.2, length.out = 40)
    f1scores <- matrix(rep(0, 4 * length(thresholds)), ncol=4)
    
    
    
    colnames(f1scores) <- c("threshold", "f1score", "precision", "recall")
    f1scores[, 1] <- thresholds
    av <- av_glasso(data, C=C)
    print(paste0("Using C=", C, ", we select the regularization parameter r=", round(av$TuningParameter, 2),"."))
    thav <- threshold(av, C=C)
    print(paste0("The thAV (C=", C, ") F1-score is: ", round(f1score(theta, thav)$f1, 2)))
    for(j in 1:length(thresholds))
    {
      performance <- f1score(theta, apply(av$Thetahat, c(1,2), function(x) if(abs(x)<thresholds[j]) {return(0)} else{ return(x)}))
      f1scores[j, 2] <- performance$f1
      f1scores[j, 3] <- performance$precision
      f1scores[j, 4] <- performance$recall
    }
    
    auxiliaries[1, 1] <- av$TuningParameter
    
    if(graph == "random")
    {
      write.csv(f1scores, paste0("../data/f1thresholded_C", C, "n",n,"d", d, ".csv"),row.names = FALSE)
      write.csv(auxiliaries, paste0("../data/aux_f1thresholded_C", C, "n", n, "d", d, ".csv"), row.names=FALSE)
    }
    else
    {
      write.csv(f1scores, paste0("../data/sff1thresholded_C", C, "n",n,"d", d, ".csv"),row.names = FALSE)
      write.csv(auxiliaries, paste0("../data/aux_sff1thresholded_C", C, "n", n, "d", d, ".csv"), row.names=FALSE)
    }
  }
}

#' Computes the performance (F1, Precision, and Recall) of the thresholded Maximum-Likelihood estimator using different thresholds. The generated .txt-files are used for the adaptation_threshold.ipynb to create a plot.
#' 
#' @param n Number of samples
#' @param d Number of dimensions
#' @param graph Type of graph
#' @return Saves .txt-files that will be used by adaptation_threshold.ipynb to generate a plot
#' @export
unregularized_threshold <- function(n, d, graph="random")
{
  theta <- generateGraph(d, graph=graph)[[1]]
  sigma <- solve(theta)
  data <- scale(mvrnorm( n, mu=rep(0, d), Sigma=sigma))
  
  thresholds <- seq(0.01, 1, length.out = 40)
  
  f1scores <- matrix(rep(0, 4 * length(thresholds)), ncol=4)
  
  
  
  colnames(f1scores) <- c("threshold", "f1score", "precision", "recall")
  f1scores[, 1] <- thresholds
  
  empvar <- var(data)
  mle <- solve(empvar)
  mle <- normalize_theta(mle)
  
  for(j in 1:length(thresholds))
  {
    performance <- f1score(theta, apply(mle, c(1,2), function(x) if(abs(x)<thresholds[j]) {return(0)} else{ return(x)}))
    f1scores[j, 2] <- performance$f1
    f1scores[j, 3] <- performance$precision
    f1scores[j, 4] <- performance$recall
  }
  
  if(graph == "random")
  {
    write.csv(f1scores, paste0("../data/f1thresholded_mlen",n,"d", d, ".csv"),row.names = FALSE)
  }
  else
  {
    write.csv(f1scores, paste0("../data/sff1thresholded_mlen" ,n,"d", d, ".csv"),row.names = FALSE)
  }
}

#' Simulation study of the recovery performance of thAV, StARS, scaled lasso, and TIGER for graphs with varying density.
#' We split the simulation study into several parts due to memory issues.
#' 
#' @param n Number of samples
#' @param d Number of dimensions
#' @param density Value for p, which is the probability to connect 2 edges of a Gilbert graph (this is the way how the random graph is generated).
#' @param part Number that describes the part of the simulation study
#' @param graph Type of graph. If graph=NULL, this function does the simulation study for random and scale-free graphs, subsequentially
#' @param num_reps Number of repetitions of the recovery task for each recovery method
#' @return Writes a .txt-file that can be read by the comparisonf1 function to summarize the simulation study. 
#' @export
comparison_density_split <- function( n, d, seq_density, part, graph="random", num_reps=10)
{
  pb <- txtProgressBar(min = 0, max = length(seq_density) * num_reps, style = 3)
  counter_pb <- 0 
  for(density in seq_density)
  {
    average_f1 <- matrix(rep(0, 7), ncol=1)
    std_f1 <- matrix(rep(0, 7), ncol=1)
    f1_results <- matrix(rep(0, num_reps * 7), ncol=num_reps)
    
    for( j in 1:num_reps)
    {
      theta <- generateGraph(d, graph=graph, prob=density)[[1]]
      data <- scale(mvrnorm(n, mu=rep(0, d), Sigma=solve( theta)))
      thAV <- thAV.estimator(data)
      f1_results[7, j] <- f1score(theta, thAV)$f1
      remove(thAV)
      
      huge_path <- huge(data, method="glasso", verbose=FALSE)
      stars <- huge.select(huge_path, criterion="stars", verbose=FALSE)$opt.icov
      f1_results[1, j] <- f1score(theta, stars)$f1
      remove(stars)
      
      scaled_lasso <- scalreg(data)$precision
      f1_results[2, j] <- f1score(theta, scaled_lasso)$f1
      
      remove(scaled_lasso)
      
      #tiger
      tiger <- huge.tiger(data, verbose=FALSE)
      tiger_tp <- sqrt(log(d) / n)
      tiger_est <- tiger$icov[[which.min(abs(tiger$lambda - tiger_tp))]]
      f1_results[3, j] <- f1score(theta, tiger_est)$f1
      
      remove(tiger)
      remove(tiger_est)
      
      
      score <- score_ebic(data)
      f1_results[4, j] <- f1score(theta, score)$f1
      remove(score)
      
      scio_cv <- scio.cv(data)$w
      f1_results[5, j] <- f1score(theta, scio_cv)$f1
      remove(scio_cv)
      
      scio_breg <- scio_bregman(data)
      f1_results[6, j] <- f1score(theta, scio_breg)$f1
      
      
      counter_pb <- counter_pb + 1
      setTxtProgressBar(pb, counter_pb)
    }
    
    close(pb)
    PATH <- paste0("../data/simu_results/comparison_density/f1n", n, "d", d , graph, part, "density", round(density, 3), ".txt")
    write.csv(f1_results, PATH, row.names=FALSE)
    print(paste0("The results were saved in ", PATH)) 
  }
}
 
#' Reads the .txt files created by the comparison_density_split function to summarize the resulting F1-score, Precision and Recall using the oracle, thAV, StARS, scaled lasso, and TIGER.
#' 
#' @param n Number of samples
#' @param d Number of dimensions
#' @param seq_density Sequence with values for p, which is the probability to connect 2 edges of a Gilbert graph (this is the way how the random graph is generated).
#' @param graph Type of graph. If graph=NULL, this function does the simulation study for random and scale-free graphs, subsequentially
#' @param num_parts Number of parts in which we split the simulation study
#' @return Matrix that measures average (and standard deviation) of the F1-score between the true graph and an estimation (thAV, StARS, and RIC)
#' @export
comparison_density <- function(n, d, seq_density, graph, num_parts)
{
  f1_avg <- matrix(rep(0, length(seq_density) * 7), ncol=length(seq_density))
  rownames(f1_avg) <- c("StARS", "scaled lasso", "TIGER", "score (eBIC)", "SCIO (CV)", "SCIO (Bregman)", "thAV")
  colnames(f1_avg) <- seq_density
  f1_std <- f1_avg
  for(j in 1:length(seq_density))
  {
    density <- seq_density[j]
    print(density)
    for(part in 1:num_parts)
    {
      PATH <- paste0("../data/simu_results/comparison_density/f1n", n, "d", d , graph, part, "density", round(density, 3), ".txt")
      readf1_results <- as.matrix(read.csv(PATH), header=TRUE)
      
      if(part==1)
      {
        f1_results <- readf1_results
      }
      else
      {
        f1_results <- cbind(f1_results, readf1_results)
      }
    }
    f1_avg[, j] <- apply(f1_results, 1, function(x) mean(x))
    f1_std[, j] <- apply(f1_results, 1, function(x) std(x))
  }
  
  summary <- conc_matrix(f1_avg, f1_std)
  
  PATH = paste0("../data/simu_results/comparison_density/n", n, "d", d , graph, ".txt")
  file_results <- file(PATH)
  writeLines(stargazer(summary), file_results)
  print(paste0("The results were saved in ", PATH))
  close(file_results)
}

#' Simulation study about the recovery of the most significant precision matrix entries
#' We split the simulation into several parts due to memory issues.
#'
#' @param n Sample size
#' @param d Size of the graph
#' @param part Number that describes the part of the simulation study
#' @param perc_top Percentage of most significant precision matrix entries that are taken into consideration
#' @param num_reps Number of repetitions which are run during the simulation
#' @param graph Type of the graph
#' @return Writes a file that is read by the recovery_topvalues function to summarize the experiment
#' @export
recovery_topvalues_split <- function(n, d, part, perc_top=0.2, num_reps=10, graph="random")
{
  pb <- txtProgressBar(min = 0, max = num_reps, style = 3)
  
  recovery <- matrix(rep(0, 7 * num_reps), ncol=num_reps)
  score_thav <- score_stars <- score_scaledlasso <- score_tiger <- score_rsme <- score_sciocv <- score_sciobreg <- rep(0, num_reps)
  for(j in 1:num_reps)
  {
    theta <- generateGraph(d, graph=graph)[[1]]
    data <- scale(mvrnorm(n, mu=rep(0, d), Sigma=solve(theta)))
    
    thav <- thAV.estimator(data)[upper.tri(theta)]
    
    huge_path <- huge(data, method="glasso", verbose=FALSE)
    stars <- huge.select(huge_path, criterion="stars", verbose=FALSE)$opt.icov[upper.tri(theta)]
    
    scaled_lasso <- scalreg(data)$precision[upper.tri(theta)]
    
    tiger <- huge.tiger(data, verbose=FALSE)
    tiger_tp <- sqrt(2) * sqrt(log(d)/(2*n))
    tiger_est <- tiger$icov[[which.min(abs(tiger$lambda - tiger_tp))]][upper.tri(theta)]
    
    rsme_ebic1 <- score_ebic(data, gamma=1)[upper.tri(theta)]
    
    scio_cv <- scio.cv(data)$w[upper.tri(theta)]
    
    scio_bregman_est <- scio_bregman(data)[upper.tri(theta)]
    
    top_values <- sort( abs(theta[ upper.tri(theta)]), decreasing=TRUE )
    num_edges <- sum(theta[upper.tri(theta)]!=0)
    
    num_top <- as.integer(perc_top * num_edges)
    lower_bound <- top_values[num_top]
    
    top_coords <- which(abs(theta[upper.tri(theta)]) >= lower_bound)
    top_thav <- top_stars <- top_scaledlasso <- top_tiger <- top_rsme <- top_sciocv <- top_sciobreg <- rep(0, num_top)
    counter_thav <- counter_stars <- counter_scaledlasso <- counter_tiger <- counter_rsme <- counter_sciocv <- counter_sciobreg <- 0
    
    for( k in 1:num_top)
    {
      top_thav[k] <- thav[top_coords[k]]
      top_stars[k] <- stars[top_coords[k]]
      top_scaledlasso[k] <- scaled_lasso[top_coords[k]]
      top_tiger[k] <- tiger_est[top_coords[k]]
      top_rsme[k] <- rsme_ebic1[top_coords[k]]
      top_sciocv[k] <- scio_cv[top_coords[k]]
      top_sciobreg[k] <- scio_bregman_est[top_coords[k]]
      
      if( top_stars[k] != 0)
      {
        counter_stars <- counter_stars + 1
      }
      if( top_thav[k] != 0)
      {
        counter_thav <- counter_thav + 1
      }
      if( top_scaledlasso[k] != 0)
      {
        counter_scaledlasso <- counter_scaledlasso + 1
      }
      if( top_tiger[k] != 0)
      {
        counter_tiger <- counter_tiger + 1
      }
      if( top_rsme[k] != 0)
      {
        counter_rsme <- counter_rsme + 1
      }
      if( top_sciocv[k] != 0)
      {
        counter_sciocv <- counter_sciocv + 1
      }
      if( top_sciobreg[k] != 0)
      {
        counter_sciobreg <- counter_sciobreg + 1
      }
    }
    
    score_thav[j] <- counter_thav / num_top
    score_stars[j] <- counter_stars / num_top
    score_scaledlasso[j] <- counter_scaledlasso / num_top
    score_tiger[j] <- counter_tiger / num_top  
    score_rsme[j] <- counter_rsme / num_top
    score_sciocv[j] <- counter_sciocv / num_top
    score_sciobreg[j] <- counter_sciobreg / num_top
    
    remove(huge_path)
    remove(theta)
    remove(data)
    remove(stars)
    remove(thav)
    remove(scaled_lasso)
    remove(tiger)
    remove(scio_cv)
    remove(scio_bregman_est)
    remove(rsme_ebic1)
    
    setTxtProgressBar(pb, j)
  }
  recovery[1, ] <- score_stars
  recovery[2, ] <- score_scaledlasso
  recovery[3, ] <- score_tiger
  recovery[4, ] <- score_rsme
  recovery[5, ] <- score_sciocv
  recovery[6, ] <- score_sciobreg
  recovery[7, ] <- score_thav
  
  PATH <- paste0("../data/simu_results/recovery_topedges/n", n, "d", d , graph, part, ".txt")
  write.csv(recovery, PATH, row.names=FALSE)
  print(paste0("The results were saved in ", PATH)) 
  close(pb)
}

#' Simulation study about the recovery of the most significant precision matrix entries.
#' This function concatinates the results of recovery_topvalues_split and returns a summary.
#'
#' @param n_seq Sequence of sample size
#' @param d_seq Sequence of dimensions 
#' @param num_parts Sequence of number of parts, in each setting n_seq[j], d_seq[j]
#' @param graph Type of the graph
#' @return Mean and standard deviation of the proportional amount of most significant edges that got recovered
#' @export
recovery_topvalues <- function(n_seq, d_seq, num_parts=c(2,3,2), graph="random")
{
  recovery_avg <- recovery_std <- matrix(rep(0, length(n_seq) * 7), ncol=length(n_seq))
  rownames(recovery_avg) <- rownames(recovery_std) <- c("StARS", "scaled lasso", "TIGER", "rSME", "SCIO (CV)", "SCIO (bregman)", "thAV")
  
  for(j in 1:length(n_seq))
  {
    n <- n_seq[j]
    d <- d_seq[j]
    for(part in 1:num_parts[j])
    {
      PATH <- paste0("../data/simu_results/recovery_topedges/n", n, "d", d , graph, part, ".txt")
      read_results <- as.matrix(read.csv(PATH), header=TRUE)
      
      if(part==1)
      {
        results <- read_results
      }
      else
      {
        results <- cbind(results, read_results)
      }
    }
    recovery_avg[, j] <- apply(results, 1, function(x) mean(x))
    recovery_std[, j] <- apply(results, 1, function(x) std(x))
  }
  
  summary <- conc_matrix(recovery_avg, recovery_std)
  
  PATH = paste0("../data/simu_results/recovery_topedges/results", graph, ".txt")
  file_results <- file(PATH)
  writeLines(stargazer(summary), file_results)
  print(paste0("The results were saved in ", PATH))
  close(file_results)
}
  
#' Exemplary recovery in F1-score of the top edges for several settings
#'
#' @param n_seq Sequence of sample sizes used for the simulation
#' @param d_seq Sequence of dimensionaly of the precision matrix used for the simulation
#' @param perc_top Percentage of most significant precision matrix entries that are taken into consideration
#' @param graph Type of the graph
#' @param latex If TRUE, returns a latex output
#' @return Table including the largest, the num_top'th largest absolute value of Theta, the reference value 6Cr, and the proportional amount of recovered top edges for several settings
#' @export 
thav_topvalues <- function(n_seq, d_seq, perc_top=0.2, graph="random", latex=TRUE)
{
  summary <- matrix( rep(0, 4 * length(n_seq)), ncol=length(n_seq))
  rownames(summary) <- c("top1", paste0("top ",perc_top, "%"), "4Cr", "occurence%")
  for(num_setting in 1:length(n_seq))
  {
    d <- d_seq[num_setting]
    n <- n_seq[num_setting]
    
    theta <- generateGraph(d, graph=graph)[[1]]
    data <- scale(mvrnorm(n, mu=rep(0, d), Sigma=solve(theta)))
    
    av <- av_glasso(data)
    av_tp <- av$TuningParameter
    thav <- threshold(av, C=0.7, lambda=1)[upper.tri(theta)]
    
    top_values <- sort( abs(theta[ upper.tri(theta)]), decreasing=TRUE )
    num_edges <- sum(theta[upper.tri(theta)]!=0)
    
    num_top <- as.integer(perc_top * num_edges)
    lower_bound <- top_values[num_top]
    top_coords <- which(abs(theta[upper.tri(theta)]) >= lower_bound)
    
    counter_thav <- 0
    for(j in 1:num_top)
    {
      if( thav[top_coords[j]] != 0)
      {
        counter_thav <- counter_thav + 1
      }
    }
    summary[, num_setting] <- c(top_values[1], top_values[num_top], 4 * 0.7 * av_tp, counter_thav / num_top)
  }
  stargazer( round(summary, 2))
}

#' Computes the performance of the thAV for large-scale examples
#' 
#' @param d Number of dimensions
#' @param n Number of samples
#' @param seq_r Set of possible regularization parameters for the thAV
#' @param C Constant C used to tune the thAV
#' @param graph Type of graph. If graph=NULL, this function does the simulation study for random and scale-free graphs, subsequentially
#' @param lambda Threshold factor used to threshold the AV. We thresholdung by lambda*C*r.
#' @param num_reps Number of repetitions of the recovery task for each recovery method
#' @param latex If TRUE, this function returns a latex output via the stargazer package
#' @return A Matrix that consists of the mean of each performance measure (F1-score, Precision, and Recall) and the standard deviation in brackets.
#' @export
performance_thav <- function(d, n=500, seq_r=seq(0.05, 0.4, length.out = 40), C=0.7, graph="random", lambda=1, num_reps=20, latex=TRUE)
{
  performance <- matrix(rep(0, 3 * num_reps), ncol=3)
  colnames(performance) <- c("F1", "Precision", "Recall")
  pb <- txtProgressBar(min=0, max=num_reps ,style=3)
  for(j in 1:num_reps)
  {
    theta <- generateGraph(d, graph=graph)[[1]]
    data <- scale(mvrnorm(n, mu=rep(0, d), Sigma=solve(theta)))
    
    thav <- thAV.estimator(data, C=C, lambda=lambda)
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
    PATH = paste0("../data/simu_results/largescale_experiments/n", n, "d", d , graph, ".txt")
    file_results <- file(PATH)
    writeLines(stargazer(summary), file_results)
    print(paste0("The results were saved in ", PATH))
    close(file_results)
  }
  else
  {
    print(summary)
  }
}

#' Investigates the empirical distribution of the absolute values of the graphs edges
#'
#' @param d Number of dimensions
#' @param graph Graph type
#' @param num_reps Number of repetitions
#' @return Saves a .txt file that can be read by visualizations_weights.ipynb to plot a histogram of the weight distribution
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
  print(paste0("The results were saved in ", file_name))
}

#' Computes F1-score, Precision and Recall using the thAV-version of the rSME estimator or the glasso with power law regularization (proposed by Liu & Ihler (2011). 
#' The results are averaged over several repetitions.
#' 
#' @param n Number of samples
#' @param d Number of dimensions
#' @param estimator Type of estiamtor. Either thav_rsme or thav_sf
#' @param C Constant C used to tune the thAV. Default value is 1.5
#' @param lambda Thresholding parameter lambda. Default value is 1.5
#' @param num_reps Number of repetitions of the recovery task. Default value is 10.
#' @param graph Type of graph. One can choose between "scale-free" and "random". Default value is "scale-free".
#' @return Writes a .txt-file that can be read by the comparisonf1 function to summarize the simulation study. 
#' @export
simu_thav_others <- function(n, d, estimator, C=1.5, lambda=1.5, num_reps=10, graph="scale-free")
{
  f1 <- precision <- recall <- rep(0, num_reps)
  pb <- txtProgressBar(min = 0, max = num_reps, style = 3)
  for(j in 1:num_reps)
  {
    theta <- generateGraph(d, graph=graph)[[1]]
    data <- scale(mvrnorm(n, mu=rep(0, d), Sigma=solve(theta)))
    
    if(estimator == "thav_sf")
    {
      av_sf <- av_sf_glasso(data, C=C)
      thav <- threshold(av_sf, C=C, lambda=lambda)
    }
    if(estimator == "thav_rsme")
    {
      rsme <- av_rsme(data, C=C)
      thav <- threshold(rsme, C=C, lambda=lambda)
    }
    
    
    performance <- f1score(theta, thav)
    f1[j] <- performance$f1
    precision[j] <- performance$precision
    recall[j] <- performance$recall
    setTxtProgressBar(pb, j)
  }
  results <- matrix(rep(0, 6), ncol=3)
  colnames(results) <- c("F1", "Precision", "Recall")
  rownames(results) <- c("Mean", "Std")
  results[1, 1] <- round(mean(f1), 2)
  results[2, 1] <- round(std(f1), 2)
  results[1, 2] <- round(mean(precision), 2)
  results[2, 2] <- round(std(precision), 2)
  results[1, 3] <- round(mean(recall), 2)
  results[2, 3] <- round(std(recall), 2)
  close(pb)
  if(estimator=="thav_sf")
  {
    PATH <- paste0("../data/simu_results/comparison_f1/sf_glasso/n", n, "d", d, graph, ".txt")
  }
  if(estimator=="thav_rsme")
  {
    PATH <- paste0("../data/simu_results/comparison_f1/rsme/n", n, "d", d, graph, ".txt")
  }
  file_results <- file(PATH)
  writeLines(stargazer(results), file_results)
  close(file_results)
  print(paste0("The results were saved in ", PATH))
  return(results)
}
