if (cfg$sim_run_or_update=="run") {
  
  run_on_cluster(
    
    first = {
      
      # Simulation setup
      sim <- new_sim()
      sim %<>% set_config(
        num_sim = cfg$sim_num,
        parallel = cfg$sim_parallel,
        n_cores = cfg$sim_n_cores,
        stop_at_error = cfg$sim_stop_at_error,
        batch_levels = c("data_type", "sigma", "icc", "n_sequences",
                         "n_clust_per_seq", "n_ind_per_cell", "re"),
        # return_batch_id = T,
        seed = 123,
        packages = cfg$pkgs
      )
      sim <- do.call(set_levels, c(list(sim), level_set))
      # if (!is.null(cfg$keep)) { sim %<>% set_levels(.keep=cfg$keep) }
      
      # Simulation script
      sim %<>% set_script(one_simulation)
      
    },
    
    main = { sim %<>% run() },
    
    last = {},
    
    cluster_config = cluster_config
    
  )
  
} else if (cfg$sim_run_or_update=="update") {
  
  update_sim_on_cluster(
    
    first = {
      sim <- readRDS(paste0(cluster_config$dir,"/sim.rds"))
      sim <- do.call(set_levels, c(list(sim), level_set))
    },
    
    main = { sim %<>% update_sim() },
    
    last = {},
    
    cluster_config = cluster_config
    
  )
  
}
