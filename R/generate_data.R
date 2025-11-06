#' Generate one stepped wedge dataset
#'
#' @param data_type One of c("normal", "binomial"). Type of outcome.
#' @param sigma Numeric; SD of the outcome (ignored if data_type=="binomial")
#' @param tau Numeric. SD of the cluster random effect
#' @param gamma Numeric. SD of the random time effect
#' @param beta_j Vector of length J; calendar time effects
#' @param delta_s Vector of length J-1. Exposure time treatment effects
#' @param n_sequences Integer. Number of sequences
#' @param n_clust_per_seq Integer. Number of clusters per sequence
#' @param n_ind_per_cell Integer. Number of individuals per cluster - time
#'     period cell (K)
#' @param n_extra_c Integer. Number of extra time points to add to the beginning
#'     of the design (i.e., extra control periods)
#' @param n_extra_t Integer. Number of extra time points to add to the end of
#'     of the design (i.e., extra treatment periods)
#' @return A dataframe representing a stepped wedge dataset
#' @details If `n_extra_c` and/or `n_extra_t` are used, the lengths of `beta_j`
#'     and `delta_s` should be kept the same. The first/last values of `beta_j`
#'     and the last value of `delta_s` will be "extended" to the extra time
#'     points.
generate_dataset <- function(
  data_type, sigma, tau, gamma, beta_j, delta_s, n_sequences, n_clust_per_seq,
  n_ind_per_cell, n_extra_c=0, n_extra_t=0
) {
  
  # Misc
  passed_args <- as.list(environment())
  I <- round(n_sequences * n_clust_per_seq)
  J <- round(n_sequences + 1)
  js <- c((1-n_extra_c):(J+n_extra_t))
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
    
    alpha_i <- rnorm(n=1, mean=0, sd=tau)
    
    for (j in js) {
      
      x_ij <- In(j>=crossover[i])
      s_ij <- round(x_ij*(j-crossover[i]+1))
      
      # Modify indices to account for extra time points
      s_ij_mod <- min(s_ij, length(delta_s))
      j_mod <- max(min(j, J), 1)
      
      gamma_ij <- rnorm(n=1, mean=0, sd=gamma)
      
      if (data_type=="normal") {
        
        if (x_ij==0) {
          mu_ij <- beta_j[j_mod] + alpha_i + gamma_ij
        } else if (x_ij==1) {
          mu_ij <- beta_j[j_mod] + alpha_i + gamma_ij + delta_s[s_ij_mod]
        }
        y_ij <- rnorm(n=K, mean=mu_ij, sd=sigma)
        
      } else if (data_type=="binomial") {
        
        if (x_ij==0) {
          p_ij <- expit(beta_j[j_mod] + alpha_i + gamma_ij)
        } else if (x_ij==1) {
          p_ij <- expit(beta_j[j_mod] + alpha_i + gamma_ij + delta_s[s_ij_mod])
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
    "y_ij" = y_ij_vec,
    "x_ij" = x_ij_vec,
    "s_ij" = s_ij_vec
  )
  
  attr(dat, "params") <- passed_args
  
  return (dat)
  
}
