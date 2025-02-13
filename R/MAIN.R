# Main config
cfg <- list(
  run_sims = T,
  run_it_vs_eti = F,
  run_eti_estimands = F,
  run_double_clusters = F,
  run_extra_time = F,
  run_process = F,
  sim_which = "Washout", # "Power", "Washout"
  sim_level_set = "Washout set 1", # "Power set 1", "Washout set 1"
  sim_run_or_update = "run",
  sim_num = 1000,
  # sim_num = 1, # !!!!!
  sim_parallel = F,
  sim_n_cores = 1000,
  sim_stop_at_error = F
)

# Secondary config
source("R/config.R", local=T)

# Load SimEngine + functions
{
  library(SimEngine)
  source("R/misc_functions.R", local=T)
  source("R/one_simulation.R", local=T)
  source("R/generate_data.R", local=T)
  source("R/analyze_data.R", local=T)
  source("../Misc/swPwr.r") # !!!!! Temporary until swCRTdesign is updated
}

# Set level sets
source("R/levels.R", local=T)

# Run simulation
if (cfg$run_sims) { source("R/run.R", local=T) }

# Run IT vs. ETI comparison (based on power calculator)
if (cfg$run_it_vs_eti) { source("R/it_vs_eti.R", local=T) }

# Run ETI vs. ETI comparison (based on power calculator)
if (cfg$run_eti_estimands) { source("R/eti_estimands.R", local=T) }

# Run comparison of doubling individuals vs. doubling clusters (based on power calculator)
if (cfg$run_double_clusters) { source("R/double_clusters.R", local=T) }

# Run IT vs. ETI comparison (based on power calculator)
if (cfg$run_extra_time) { source("R/extra_time.R", local=T) }

# Run IT vs. ETI comparison (based on power calculator)
# !!!!! Move this code to other folder
# if (cfg$run_unequal_alloc) { source("R/unequal_alloc.R", local=T) }

# Tables and figures
if (cfg$run_process) { source("R/process.R", local=T) }
