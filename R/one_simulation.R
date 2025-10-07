if (cfg$sim_which=="Power") {
  
  ###############################################.
  ##### Power for different analysis models #####
  ###############################################.
  
  one_simulation <- function() {
    
    # # Testing
    # L <- list(data_type="normal", sigma=2, icc=0.1, n_sequences=5, n_clust_per_seq=2, n_ind_per_cell=8, re="cluster+time", estimand="TATE", model=list(cal_time="cat", exp_time="IT", n_omit=0))
    
    # Generate data
    batch({
      if (L$tvte) {
        # delta_s <- seq(0.15, 0.3, length.out=L$n_sequences)
        delta_s <- seq(0.2, 0.4, length.out=L$n_sequences)
      } else {
        delta_s <- rep(0.2,L$n_sequences)
      }
      tau <- L$sigma * sqrt(L$icc/(1-L$icc))
      dat <- generate_dataset(
        data_type = L$data_type,
        sigma = L$sigma,
        tau = tau,
        beta_j = seq(from=0, to=1, length.out=L$n_sequences+1),
        delta_s = delta_s,
        n_sequences = L$n_sequences,
        n_clust_per_seq = L$n_clust_per_seq,
        n_ind_per_cell = L$n_ind_per_cell,
        re = L$re
      )
    })
    
    # Parse estimand
    if (L$estimand=="TATE") {
      estimand_type <- "TATE"
      estimand_time <- c(1,L$n_sequences)
    } else if (L$estimand=="PTE-1") {
      estimand_type <- "PTE"
      estimand_time <- 1
    } else if (L$estimand=="PTE-S") {
      estimand_type <- "PTE"
      estimand_time <- L$n_sequences
    }
    
    # Parse random effects
    if (L$re=="cluster+time") {
      rand_eff <- c("clust", "time")
    } else if (L$re=="cluster") {
      rand_eff <- c("clust")
    }
    
    # # Analyze data
    # results_old <- analyze_data(
    #   dat = dat,
    #   cal_time = L$model$cal_time,
    #   exp_time = L$model$exp_time,
    #   re = L$re,
    #   estimand_type = estimand_type,
    #   estimand_time = estimand_time,
    #   return_curve = F
    # )
    
    # Read in data
    suppressMessages({
      sw_dat <- steppedwedge::load_data(
        time = "j",
        cluster_id = "i",
        treatment = "x_ij",
        outcome = "y_ij",
        data = dat
      )
    })
    
    # Analyze data
    results <- steppedwedge::analyze(
      dat = sw_dat,
      estimand_type = estimand_type,
      estimand_time = estimand_time,
      exp_time = L$model$exp_time,
      cal_time = L$model$cal_time,
      re = rand_eff
    )
    
    # !!!!! TESTING
    if (F) {
      
      res_ETI <- analyze_data(
        dat = dat,
        cal_time = "linear",
        exp_time = "cat",
        re = "cluster",
        n_omit = 0,
        return_curve = T
      )
      res_NCS <- analyze_data(
        dat = dat,
        cal_time = "linear",
        exp_time = "NCS",
        re = "cluster",
        n_omit = 0,
        return_curve = T
      )
      
      df_plot <- data.frame(
        x = rep(c(1:attr(dat, "params")$n_sequences), 2),
        y = c(res_ETI$curve, res_NCS$curve),
        which = rep(c("ETI", "NCS"), each=length(res_ETI$curve))
      )
      ggplot(df_plot, aes(x=x, y=y, color=factor(which))) +
        geom_line()
      
    }
    
    # Conduct Wald-type hypothesis test
    ci <- results$te_est + c(-1,1)*qnorm(1-0.05/2)*results$te_se
    reject <- In(ci[1]>0 || ci[2]<0)
    
    # Return results
    return(list(
      reject = reject,
      est = results$te_est,
      se = results$te_se,
      true_tate = 0.2
    ))
    
  }
  
} else if (cfg$sim_which=="Washout") {
  
  ##############################################.
  ##### Power for different washout models #####
  ##############################################.
  
  one_simulation <- function() {
    
    # Generate data
    batch({
      tau <- L$sigma * sqrt(L$icc/(1-L$icc))
      dat <- generate_dataset(
        data_type = L$data_type,
        sigma = L$sigma,
        tau = tau,
        beta_j = seq(from=0, to=1, length.out=L$n_sequences+1),
        delta_s = rep(0.2,L$n_sequences),
        n_sequences = L$n_sequences,
        n_clust_per_seq = L$n_clust_per_seq,
        n_ind_per_cell = L$n_ind_per_cell,
        re = L$re
      )
    })
    
    # Parse estimand
    if (L$estimand=="TATE(0,S)") {
      estimand_time <- c(1,L$n_sequences)
      n_wash <- 0
    } else if (L$estimand=="TATE(1,S)") {
      estimand_time <- c(2,L$n_sequences)
      n_wash <- 1
    } else if (L$estimand=="TATE(2,S)") {
      estimand_time <- c(3,L$n_sequences)
      n_wash <- 2
    } else if (L$estimand=="TATE(3,S)") {
      estimand_time <- c(4,L$n_sequences)
      n_wash <- 3
    }
    
    # Parse random effects
    if (L$re=="cluster+time") {
      rand_eff <- c("clust", "time")
    } else if (L$re=="cluster") {
      rand_eff <- c("clust")
    }
    
    # Prep work
    if (L$model$exp_time=="DCT") {
      
      dat %<>% dplyr::mutate(
        s_ij = pmin(s_ij, n_wash+1),
        ij = as.integer(factor(paste0(i,"-",j)))
      )
      
    } else if (L$model$exp_time=="IT") {
      
      if (n_wash>0) {
        dat %<>% dplyr::filter(
          s_ij==0 | s_ij>n_wash
        )
      }
      
    }
    
    # Read in data
    suppressMessages({
      sw_dat <- steppedwedge::load_data(
        time = "j",
        cluster_id = "i",
        treatment = "x_ij",
        outcome = "y_ij",
        data = dat
      )
    })
    
    # Analyze data
    if (L$model$exp_time=="DCT") {
      
      if (L$model$cal_time!="categorical") {
        stop("DCT model only implemented for categorical calendar time")
      }
      if (L$re!="cluster+time") {
        stop("DCT model only implemented for re=='cluster+time'")
      }
      
      model_pit <- lme4::lmer(
        formula = "y_ij~factor(j)+factor(s_ij)+(1|i)+(1|ij)",
        data = dat
      )
      summ_pit <- summary(model_pit)
      c_name <- paste0("factor(s_ij)", round(n_wash+1))
      results <- list(
        te_est = summ_pit$coefficients[c_name,"Estimate"],
        te_se = summ_pit$coefficients[c_name,"Std. Error"]
      )
      
    } else if (L$model$exp_time=="IT") {
      
      results <- steppedwedge::analyze(
        dat = sw_dat,
        exp_time = "IT",
        cal_time = L$model$cal_time,
        re = rand_eff
      )
      
    } else if (L$model$exp_time=="ETI") {
      
      results <- steppedwedge::analyze(
        dat = sw_dat,
        estimand_type = "TATE",
        estimand_time = estimand_time,
        exp_time = "ETI",
        cal_time = L$model$cal_time,
        re = rand_eff
      )
      
    } else if (L$model$exp_time=="IT-full-temp") {
      
      results <- steppedwedge::analyze(
        dat = sw_dat,
        exp_time = "IT",
        cal_time = L$model$cal_time,
        re = rand_eff
      )
      
    }
    
    # Conduct Wald-type hypothesis test
    ci <- results$te_est + c(-1,1)*qnorm(1-0.05/2)*results$te_se
    reject <- In(ci[1]>0 || ci[2]<0)
    
    # Return results
    return(list(
      reject = reject,
      est = results$te_est,
      se = results$te_se,
      true_tate = 0.2
    ))
    
  }
  
}
