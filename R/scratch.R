
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