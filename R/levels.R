# Set simulation levels
if (cfg$run_sims && Sys.getenv("sim_run") %in% c("first", "")) {
  
  level_sets <- list()
  
  # Level set: Comparing analysis models
  # Figures: ...
  level_sets[["Power set 1"]] <- list(
    data_type = "normal",
    sigma = 1, # 1.5
    icc = c(0.01,0.05,0.1),
    cac = 0.75,
    n_sequences = c(6, 12, 18, 24),
    n_clust_per_seq = 4,
    n_ind_per_cell = 5, # 8
    tvte = TRUE,
    estimand = c("TATE", "PTE-1", "PTE-S"),
    model = list(
      # "IT" = list(cal_time="categorical", exp_time="IT", n_omit=0),
      "ETI, cat cal time" = list(cal_time="categorical", exp_time="ETI", n_omit=0),
      "ETI, linear cal time" = list(cal_time="linear", exp_time="ETI", n_omit=0),
      # "ETI, spline time" = list(cal_time="NCS", exp_time="ETI", n_omit=0),
      # "NCS, cat cal time" = list(cal_time="categorical", exp_time="NCS", n_omit=0),
      "NCS, cat cal time" = list(cal_time="categorical", exp_time="NCS", n_omit=0)
      # "NCS, linear cal time" = list(cal_time="linear", exp_time="NCS", n_omit=0)
      # "NCS, spline cal time" = list(cal_time="NCS", exp_time="NCS", n_omit=0)
    )
  )
  
  # Level set: Washout models
  # Figures: ...
  level_sets[["Washout set 1"]] <- list(
    data_type = "normal",
    sigma = 1, # 0.8
    icc = c(0.01,0.05,0.1),
    cac = 0.75,
    n_sequences = c(6, 9, 12),
    n_clust_per_seq = 4,
    n_ind_per_cell = 5, # 8
    tvte = FALSE,
    estimand = c("TATE(1,S)", "TATE(3,S)"),
    model = list(
      # "IT-full-temp" = list(cal_time="categorical", exp_time="IT-full-temp", n_omit=0), # Testing to make sure this gives highest power
      "IT" = list(cal_time="categorical", exp_time="IT", n_omit=0),
      "DCT" = list(cal_time="categorical", exp_time="DCT", n_omit=0),
      "ETI" = list(cal_time="categorical", exp_time="ETI", n_omit=0)
    )
  )
  
  level_set <- level_sets[[cfg$sim_level_set]]
  
}
