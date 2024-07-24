# Set simulation levels
if (cfg$run_sims && Sys.getenv("sim_run") %in% c("first", "")) {
  
  level_sets <- list()
  
  # Level set: ETATE estimation
  # Figures: ...
  level_sets[["Power set 1"]] <- list(
    data_type = "normal",
    sigma = 2,
    # tau = 0.5,
    icc = c(0.01,0.05,0.1),
    n_sequences = c(5, 10, 20, 40),
    # n_clust_per_seq = c(2,4),
    n_clust_per_seq = 2,
    n_ind_per_cell = 8,
    re = "cluster+time",
    estimand = c("TATE", "PTE-1", "PTE-S"),
    model = list(
      "IT" = list(cal_time="cat", exp_time="IT", n_omit=0),
      "ETI, cat cal time" = list(cal_time="cat", exp_time="cat", n_omit=0),
      "ETI, linear cal time" = list(cal_time="linear", exp_time="cat", n_omit=0),
      # "ETI, spline time" = list(cal_time="NCS", exp_time="cat", n_omit=0),
      # "NCS, cat cal time" = list(cal_time="cat", exp_time="NCS", n_omit=0),
      "NCS, linear cal time" = list(cal_time="linear", exp_time="NCS", n_omit=0)
      # "NCS, spline cal time" = list(cal_time="NCS", exp_time="NCS", n_omit=0)
      
      # "ETI, omit last 1" = list(cal_time="cat", exp_time="cat", n_omit=1),
      # "ETI, omit last 2" = list(cal_time="cat", exp_time="cat", n_omit=2)
    )
  )
  
  level_set <- level_sets[[cfg$sim_level_set]]
  
}
