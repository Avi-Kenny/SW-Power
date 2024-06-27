# Main config
cfg <- list(
  run_sims = T,
  run_it_vs_eti = F,
  run_process = T,
  sim_which = "Power",
  sim_level_set = "Power set 1",
  sim_run_or_update = "run",
  sim_num = 100,
  sim_parallel = F,
  sim_n_cores = 3,
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
}

# Set level sets
source("R/levels.R", local=T)

# Run simulation
if (cfg$run_sims) { source("R/run.R", local=T) }

# Run IT vs. ETI comparison (based on power calculator)
if (cfg$run_it_vs_eti) { source("R/it_vs_eti.R", local=T) }

# Tables and figures
if (cfg$run_process) { source("R/process.R", local=T) }
