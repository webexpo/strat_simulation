# Script for the interpretation of the simulations

# testing existing simulation and comparison across runs


##### LIBRARIES ####


##### DATA ####

load(file = "F:/Dropbox/GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/all_scenarios_run3a_sim.RDS")
run3a <- simulation_results

load(file = "F:/Dropbox/GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/all_scenarios_run3b_sim.RDS")
run3b <- simulation_results

run1 <- readRDS( "F:/Dropbox/GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/all_scenarios_run1_d.RDS")

run2 <- simulation_results


##### SCENARIOS ####

scenarios <- expand.grid(
  true_gsd        = c(1.5, 2.5, 3.5),
  true_exceedance_perc = c(0.1, 1, 3, 7, 10, 25),
  sample_size      = c(3L, 6L , 9L, 12L),
  stringsAsFactors = FALSE)
