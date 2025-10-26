# Robust standard error
if (T) {
  
  # > head(dat)
  # i j k      y_ij x_ij s_ij
  # 1 1 1 1 -4.159351    0    0
  # 2 1 1 2 -2.047745    0    0
  # 3 1 1 3  1.589721    0    0
  # 4 1 1 4  3.848683    0    0
  # 5 1 1 5 -2.607160    0    0
  # 6 1 1 6  3.042304    0    0
  
  library(lme4)
  library(clubSandwich)
  
  # IT model
  model_lmer <- lmer(
    y_ij ~ factor(j) + x_ij + (1|i),
    data = dat
  )
  sqrt(diag(vcov(model_lmer)))
  robust_var <- clubSandwich::vcovCR(model_lmer, cluster=dat$i, type="CR1")
  sqrt(diag(robust_var))
  
  # ETI model
  model_lmer <- lmer(
    y_ij ~ factor(j) + factor(s_ij) + (1|i),
    data = dat
  )
  robust_var <- clubSandwich::vcovCR(model_lmer, cluster=dat$i, type="CR1")
  sqrt(diag(vcov(model_lmer)))
  sqrt(diag(robust_var))
  
}

# Debugging staircase design
if (T) {
  
  n_clust_per_seq <- 4
  n_sequences <- 4
  n_ind_per_cell <- 10^2
  n_before_and_after <- 1
  design <- swDsn(
    clusters = rep(n_clust_per_seq, n_sequences),
    extra.ctl.time = n_before_and_after-1,
    extra.trt.time = n_before_and_after-1
  )
  n_row <- nrow(design$swDsn)
  n_col <- ncol(design$swDsn)
  n_matrix <- matrix(NA, nrow=n_row, ncol=n_col)
  
  for (row in c(1:n_row)) {
    vec <- design$swDsn[row,]
    first_tx <- which.max(vec==1)
    n_vals <- rep(0, n_col)
    ind_1 <- first_tx-n_before_and_after
    ind_2 <- first_tx+n_before_and_after-1
    n_vals[ind_1:ind_2] <- n_ind_per_cell
    n_matrix[row,] <- n_vals
  }
  
  # print(design$swDsn)
  # print(n_matrix)
  
  power <- swPwr(
    design = design,
    distn = "gaussian",
    n = n_matrix,
    # n = n_ind_per_cell,
    mu0 = 0,
    mu1 = 0.1,
    # H = H,
    sigma = 1,
    icc = 0.05,
    cac = 1,
    alpha = 0.05
  )
  print(paste("Power:", round(power,2)))
  
}

# New simulation to test bizarre issue around extra time points
if (F) {
  
  sim <- new_sim()
  sim %<>% set_levels(
    case = c("base", "add_3_c", "add_3_t")
    # case = c("base", "add_1_c", "add_1_t", "add_3_c", "add_3_t")
  )
  sim %<>% set_config(num_sim=1000)
  sim %<>% set_script(function() {
    
    if (L$case=="base") {
      n_extra_c=0; n_extra_t=0;
    } else if (L$case=="add_1_c") {
      n_extra_c=1; n_extra_t=0;
    } else if (L$case=="add_1_t") {
      n_extra_c=0; n_extra_t=1;
    } else if (L$case=="add_3_c") {
      n_extra_c=3; n_extra_t=0;
    } else if (L$case=="add_3_t") {
      n_extra_c=0; n_extra_t=3;
    }
    
    dat <- generate_dataset(
      data_type = "normal",
      sigma = 1,
      tau = sqrt(0.1/0.9),
      # tau = 0,
      beta_j = rep(0, 5),
      delta_s = rep(0.2, 4),
      n_sequences = 4,
      n_clust_per_seq = 4,
      n_ind_per_cell = 30,
      re = "cluster",
      n_extra_c = n_extra_c,
      n_extra_t = n_extra_t
    )
    
    res <- analyze_data(
      dat = dat,
      cal_time = "cat",
      exp_time = "cat",
      re = "cluster",
      estimand_type = "TATE",
      estimand_time = c(1,4),
      return_curve = F,
      return_ses = T
    )
    
    ci <- res$est + c(-1,1)*qnorm(1-0.05/2)*res$se
    reject <- In(ci[1]>0 || ci[2]<0)
    
    return (list(
      "est" = res$est,
      "se" = res$se,
      "reject" = reject,
      "se_curve_1" = res$curve_se[1],
      "se_curve_4" = res$curve_se[4]
    ))
    
  })
  sim %<>% run()
  
  print(SimEngine::summarize(
    sim,
    list(stat="mean", x="reject", name="power"),
    list(stat="mean", x="se"),
    list(stat="mean", x="se_curve_1"),
    list(stat="mean", x="se_curve_4")
  ))
  #   level_id    case n_reps power
  # 1        1    base   1000 0.430
  # 2        2 add_1_c   1000 0.550
  # 3        3 add_1_t   1000 0.482
  
  # swCRTdesign power estimates
  # Base: 42.6%
  # Add 1 C: 52.5%
  # Add 1 T: 44.2%
  # Add 3 C: 63.3%
  # Add 3 T: 44.5%
  
}

# New simulation to test bizarre issue around extra time points (smaller design)
if (F) {
  
  sim <- new_sim()
  sim %<>% set_levels(
    case = c("base", "add_1_c", "add_1_t")
  )
  sim %<>% set_config(num_sim=1000)
  sim %<>% set_script(function() {
    
    if (L$case=="base") {
      n_extra_c=0; n_extra_t=0;
    } else if (L$case=="add_1_c") {
      n_extra_c=1; n_extra_t=0;
    } else if (L$case=="add_1_t") {
      n_extra_c=0; n_extra_t=1;
    }
    
    dat <- generate_dataset(
      data_type = "normal",
      sigma = 1,
      tau = sqrt(0.1/0.9),
      # tau = 0,
      beta_j = rep(0, 3),
      delta_s = rep(0.2, 2),
      n_sequences = 2,
      n_clust_per_seq = 4,
      n_ind_per_cell = 30,
      re = "cluster",
      n_extra_c = n_extra_c,
      n_extra_t = n_extra_t
    )
    
    res <- analyze_data(
      dat = dat,
      cal_time = "cat",
      exp_time = "cat",
      re = "cluster",
      estimand_type = "TATE",
      estimand_time = c(1,2),
      return_curve = F,
      return_ses = T
    )
    
    ci <- res$est + c(-1,1)*qnorm(1-0.05/2)*res$se
    reject <- In(ci[1]>0 || ci[2]<0)
    
    return (list(
      "est" = res$est,
      "se" = res$se,
      "reject" = reject,
      "se_curve_1" = res$curve_se[1],
      "se_curve_2" = res$curve_se[2]
    ))
    
  })
  sim %<>% run()

  print(SimEngine::summarize(
    sim,
    list(stat="mean", x="reject", name="power"),
    list(stat="mean", x="se"),
    list(stat="mean", x="se_curve_1"),
    list(stat="mean", x="se_curve_2")
  ))
  
  #   level_id    case n_reps  power   mean_se mean_se_curve_1 mean_se_curve_2
  # 1        1    base  10000 0.1672 0.2175115       0.1686696       0.2835353
  # 2        2 add_1_c  10000 0.1935 0.1904643       0.1534062       0.2465464
  # 3        3 add_1_t  10000 0.1606 0.2180981       0.1690114       0.2843349
  
}


# New simulation to test bizarre issue around extra time points (large design)
if (T) {
  
  sim <- new_sim()
  sim %<>% set_levels(
    case = c("base", "add_3_c", "add_3_t")
  )
  sim %<>% set_config(num_sim=1000)
  sim %<>% set_script(function() {
    
    if (L$case=="base") {
      n_extra_c=0; n_extra_t=0;
    } else if (L$case=="add_3_c") {
      n_extra_c=3; n_extra_t=0;
    } else if (L$case=="add_3_t") {
      n_extra_c=0; n_extra_t=3;
    }
    
    dat <- generate_dataset(
      data_type = "normal",
      sigma = 1,
      tau = sqrt(0.1/0.9),
      # tau = 0,
      beta_j = rep(0, 8),
      delta_s = rep(0.2, 7),
      n_sequences = 7,
      n_clust_per_seq = 2,
      n_ind_per_cell = 30,
      re = "cluster",
      n_extra_c = n_extra_c,
      n_extra_t = n_extra_t
    )
    
    res <- analyze_data(
      dat = dat,
      cal_time = "cat",
      exp_time = "cat",
      re = "cluster",
      estimand_type = "TATE",
      estimand_time = c(1,7),
      return_curve = F,
      return_ses = T
    )
    
    ci <- res$est + c(-1,1)*qnorm(1-0.05/2)*res$se
    reject <- In(ci[1]>0 || ci[2]<0)
    
    return (list(
      "est" = res$est,
      "se" = res$se,
      "reject" = reject,
      "se_curve_1" = res$curve_se[1],
      "se_curve_2" = res$curve_se[2]
    ))
    
  })
  sim %<>% run()
  
  print(SimEngine::summarize(
    sim,
    list(stat="mean", x="reject", name="power"),
    list(stat="mean", x="se"),
    list(stat="mean", x="se_curve_1"),
    list(stat="mean", x="se_curve_2")
  ))
  
  #   level_id    case n_reps power    mean_se mean_se_curve_1 mean_se_curve_2
  # 1        1    base   1000 0.544 0.09605261      0.06900267      0.08025247
  # 2        2 add_3_c   1000 0.730 0.07668250      0.06374699      0.07220469
  # 3        3 add_3_t   1000 0.576 0.09230057      0.06871098      0.07849939
  
}

# Testing out time trend function
if (F) {
  
  dat <- generate_dataset(
    mu, tau, theta, n_clusters, n_time_points, n_ind_per_cluster, data_type,
    sigma=NA, delay_model, n_extra_time_points, rte=NA, time_trend="incr"
  )
  
}

# Effect curves
if (F) {
  
  #' Return the effect curve value (% of Tx effect reached as a fn of time steps)
  #'
  #' @param x Value at which to evaluate the effect curve
  #' @param type One of the following character strings:
  #'     - "exp": Hussey & Hughes; ignores delay
  #'     - "spline": "exposure treatment indicators" model
  #' @param params List of parameters, which differs by "type"
  #'     - For "exp", `d` is the rate parameter
  #'     - For "spline", `knots` specifies the knot locations (the first of
  #'       which must be zero) and `slopes` represents the slopes between
  #'       knots. Spline starts at (0,0) and has slope zero after the last knot.
  #'     - For "parabola", paramaters `a`, `b`, and `c` correspond to the parabola
  #'       given by y = a(x^2) + bx + c
  #' @return A number between zero and one, representing the % of the Tx effect
  #'     reached after x time steps
  
  effect_curve <- function(x, type, params) {
    
    p <- params
    
    if (type=="exp") { y <- ifelse(x>6, 1, (1-exp(-x/p$d)) / (1-exp(-6/p$d))) }
    
    if (type=="spline") {
      
      if ((length(p$knots)-1)!=length(p$slopes)) {
        stop("Length of `knots` must equal length of `slopes` plus one")
      }
      
      y <- p$slopes[1] * pmax(0,x)
      
      for (i in 2:length(p$knots)) {
        if (i<length(p$knots)) {
          y <- y + (p$slopes[i] - p$slopes[i-1]) * pmax(0,x-p$knots[i])
        } else {
          y <- y - p$slopes[i-1] * pmax(0,x-p$knots[i])
        }
      }
      
    }
    
    return (y)
    
  }
  
  effect_curve(x, type, params)
  
  # Create theta_s vector (intervention effects) based on "delay_model" function
  theta_ls <- theta * effect_curve(
    x = 1:(n_time_points-1),
    type = delay_model$type,
    params = delay_model$params
  )
  
}

# Levels (original)
if (F) {
  
  # Set global constants
  C <- list(
    mu = log(0.1),
    delay_models = list(
      "Instantaneous" = list(
        type = "spline",
        params = list(knots=c(0,0.1), slopes=10)
      ),
      "Lagged" = list(
        type = "spline",
        params = list(knots=c(0,2,2.1), slopes=c(0,10))
      ),
      "Curved" = list(
        type = "exp",
        params = list(d=1.5)
      ),
      "Partially convex" = list(
        type = "spline",
        params = list(knots=c(0,2,4), slopes=c(0.1,0.4))
      )
    )
  )
  
  if (cfg$run_sims && Sys.getenv("sim_run") %in% c("first", "")) {
    
    level_sets <- list()
    
    # Pitfalls of the immediate treatment (IT) model, estimation of the TATE and
    #   LTE, estimation of the entire effect curve
    # Figures: sim2, sim3
    # !!!!! Formerly `level_set_123`
    level_sets[["primary_1"]] <- list(
      n_clusters = 24,
      n_time_points = 7,
      n_ind_per_cluster = 20,
      theta = 0.5,
      tau = 0.25,
      sigma = 1,
      data_type = "normal",
      method = list(
        "IT" = list(method="IT"),
        "ETI" = list(method="ETI"),
        "NCS (4df)" = list(method="NCS-4df"),
        "MEC" = list(method="MCMC-MON-Stan",
                     enforce="simplex",
                     mcmc=cfg$mcmc)
      ),
      delay_model = C$delay_models,
      n_extra_time_points = 0,
      rte = NA,
      return_extra = list("whole_curve"=list(whole_curve=TRUE))
    )
    
    # Power of Wald-type hypothesis tests
    # Figures: sim4
    # !!!!! Formerly `level_set_4`
    level_sets[["primary_2"]] <- list(
      n_clusters = 24,
      n_time_points = 7,
      n_ind_per_cluster = 20,
      theta = seq(0,0.5,0.1),
      tau = 0.25,
      sigma = 1,
      data_type = "normal",
      method = list(
        "IT" = list(method="IT"),
        "ETI" = list(method="ETI"),
        "NCS (4df)" = list(method="NCS-4df"),
        "MEC" = list(method="MCMC-MON-Stan",
                     enforce="simplex",
                     mcmc=cfg$mcmc)),
      delay_model = C$delay_models,
      n_extra_time_points = 0,
      rte = NA,
      return_extra = list("none"=list())
    )
    
    # Performance of RTE models (data generated with/without random treatment effect)
    # Figures: sim6, sim7
    # !!!!! Formerly `level_set_67`
    level_sets[["primary_3"]] <- list(
      n_clusters = 24,
      n_time_points = 7,
      n_ind_per_cluster = 20,
      theta = 0.5,
      tau = 0.25,
      sigma = 1,
      data_type = "normal",
      method = list(
        "ETI" = list(method="ETI"),
        "ETI (RTE; height)" = list(method="ETI", re="height")
      ),
      delay_model = C$delay_models,
      n_extra_time_points = 0,
      rte = list(
        "none" = list(),
        "height" = list(type="height", nu=1, rho1=-0.2)
      ),
      return_extra = list("rte"=list(rte=TRUE))
    )
    
    # Effect of adding extra time points
    # Figures: sim8
    # !!!!! Formerly `level_set_8`
    level_sets[["primary_4"]] <- list(
      n_clusters = 24,
      n_time_points = 7,
      n_ind_per_cluster = 20,
      theta = 0.5,
      tau = 0.25,
      sigma = 1,
      data_type = "normal",
      method = list("ETI" = list(method="ETI")),
      delay_model = C$delay_models,
      n_extra_time_points = c(0,1,2),
      rte = NA,
      return_extra = list("none"=list())
    )
    
    level_set <- level_sets[[cfg$sim_level_set]]
    
    # if (cfg$sim_level_set=="asdf") { cfg$keep = c(1:3,7:9,16:18,22:24) }
    
  }
  
}