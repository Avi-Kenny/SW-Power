######################.
##### Estimation #####
######################.

if (cfg$sim_which=="Power") {
  
  #' Run a single simulation
  #'
  #' @return A list with ...
  one_simulation <- function() {
    
    # # Testing
    # L <- list(data_type="normal", sigma=1, tau=1, n_sequences=6, n_clust_per_seq=4, n_ind_per_cell=20)
    
    # Generate data
    batch({
      dat <- generate_dataset(
        data_type = L$data_type,
        sigma = L$sigma,
        tau = L$tau,
        beta_j = rep(0,L$n_sequences+1),
        delta_s = rep(0.2,L$n_sequences),
        n_sequences = L$n_sequences,
        n_clust_per_seq = L$n_clust_per_seq,
        n_ind_per_cell = L$n_ind_per_cell
      )
    })
    
    # Analyze data
    results <- analyze_data(
      dat = dat,
      cal_time = L$model$cal_time,
      exp_time = L$model$exp_time,
      re = L$re,
      n_omit = L$model$n_omit
      # target = "ETATE"
    )
    
    # Conduct Wald-type hypothesis test
    ci <- results$est + c(-1,1)*qnorm(1-0.05/2)*results$se
    reject <- In(ci[1]>0 || ci[2]<0)
    
    # Return results
    return(list(
      reject = reject,
      est = results$est,
      se = results$se,
      true_tate = 0.2
    ))
    
  }
  
}
