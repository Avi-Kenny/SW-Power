# Set simulation levels
if (cfg$run_sims && Sys.getenv("sim_run") %in% c("first", "")) {
  
  level_sets <- list()
  
  # Level set: ETATE estimation
  # Figures: ...
  level_sets[["Power set 1"]] <- list(
    data_type = "normal",
    sigma = 1,
    tau = 0.25,
    # n_sequences = 6,
    n_sequences = c(6, 10, 14),
    n_clust_per_seq = 2,
    n_ind_per_cell = 10,
    re = "cluster",
    # re = c("cluster", "cluster+time"),
    model = list(
      "IT" = list(cal_time="cat", exp_time="IT", n_omit=0),
      "ETI" = list(cal_time="cat", exp_time="cat", n_omit=0),
      "ETI, linear cal time" = list(cal_time="linear", exp_time="cat", n_omit=0),
      # # "ETI, spline time" = list(cal_time="NCS", exp_time="cat", n_omit=0),
      # # "NCS, spline time" = list(cal_time="NCS", exp_time="NCS", n_omit=0),
      # "ETI, no time" = list(cal_time="none", exp_time="cat", n_omit=0),
      
      
      
      "ETI, omit last 1" = list(cal_time="cat", exp_time="cat", n_omit=1),
      "ETI, omit last 2" = list(cal_time="cat", exp_time="cat", n_omit=2)
    )
  )
  
  level_set <- level_sets[[cfg$sim_level_set]]
  
}
