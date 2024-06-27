#' Generate one stepped wedge dataset
#'
#' @param data_type One of c("normal", "binomial"). Type of outcome.
#' @param sigma Numeric; SD of the outcome (ignored if data_type=="binomial")
#' @param tau Numeric. SD of the cluster random effect
#' @param beta_j Vector of length J; calendar time effects
#' @param delta_s Vector of length J-1. Exposure time treatment effects
#' @param n_sequences Integer. Number of sequences
#' @param n_clust_per_seq Integer. Number of clusters per sequence
#' @param n_ind_per_cell Integer. Number of individuals per cluster - time
#'     period cell (K)
#' @param re One of c("cluster", "cluster+time"); whether a cluster intercept
#'     should be included ("cluster") or a cluster intercept plus a
#'     cluster-period intercept ("cluster+time")
#' @return A dataframe representing a stepped wedge dataset
generate_dataset <- function(
  data_type, sigma, tau, beta_j, delta_s, n_sequences, n_clust_per_seq,
  n_ind_per_cell, re
) {
  
  # Misc
  passed_args <- as.list(environment())
  I <- round(n_sequences * n_clust_per_seq)
  J <- round(n_sequences + 1)
  K <- n_ind_per_cell
  i_vec <- j_vec <- k_vec <- y_ij_vec <- x_ij_vec <- s_ij_vec <- rep(NA, I*J*K)
  
  # Dataset checks
  if (length(delta_s)!=(J-1)) { stop("Length of delta_s must equal J-1.") }
  if (length(beta_j)!=J) { stop("Length of beta_j must equal J.") }
  
  clusters <- c(1:I)
  crossover <- rep(c(2:J), each=round(I/(J-1)))
  crossover <- sample(crossover)
  pointer <- 1
  for (i in clusters) {
    
    if (re=="cluster") {
      alpha_i <- rnorm(n=1, mean=0, sd=tau)
    } else if (re=="cluster+time") {
      alpha_i <- rnorm(n=1, mean=0, sd=tau/sqrt(2))
    }
    
    for (j in c(1:J)) {
      
      x_ij <- In(j>=crossover[i])
      s_ij <- round(x_ij*(j-crossover[i]+1))
      
      if (re=="cluster") {
        gamma_ij <- 0
      } else if (re=="cluster+time") {
        gamma_ij <- rnorm(n=1, mean=0, sd=tau/sqrt(2))
      }
      
      if (data_type=="normal") {
        
        if (x_ij==0) {
          mu_ij <- beta_j[j] + alpha_i + gamma_ij
        } else if (x_ij==1) {
          mu_ij <- beta_j[j] + alpha_i + gamma_ij + delta_s[s_ij]
        }
        y_ij <- rnorm(n=K, mean=mu_ij, sd=sigma)
        
      } else if (data_type=="binomial") {
        
        if (x_ij==0) {
          p_ij <- expit(beta_j[j] + alpha_i + gamma_ij)
        } else if (x_ij==1) {
          p_ij <- expit(beta_j[j] + alpha_i + gamma_ij + delta_s[s_ij])
        }
        y_ij <- rbinom(n=K, size=1, prob=p_ij)
        
      }
      
      inds <- c(pointer:(pointer+(K-1)))
      i_vec[inds] <- rep(i, K)
      j_vec[inds] <- rep(j, K)
      k_vec[inds] <- c(1:K)
      y_ij_vec[inds] <- y_ij
      x_ij_vec[inds] <- rep(x_ij, K)
      s_ij_vec[inds] <- rep(s_ij, K)
      pointer <- pointer + K
      
    }
    
  }
  
  dat <- data.frame(
    "i" = i_vec,
    "j" = j_vec,
    "k" = k_vec,
    # "c_i" = c_i_vec, # crossover time (i.e., start time of intervention)
    "y_ij" = y_ij_vec,
    "x_ij" = x_ij_vec,
    "s_ij" = s_ij_vec
  )
  
  attr(dat, "params") <- passed_args
  
  return (dat)
  
}
