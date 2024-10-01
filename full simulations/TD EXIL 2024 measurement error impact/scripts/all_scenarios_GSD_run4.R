# Script for performing a simulation study with real GSD : 5000 iterations per scenario using STAN

# Script for the second official run of 3  

##### LIBRARIES ####

library(tolerance)
library(parallel)
library(rstan)

##### DATA ####

real_gsds <- readRDS( "created data/real_gsd_values.RDS")

##### SCRIPTS ####

source("full simulations/TD EXIL 2024 measurement error impact/scripts/sample_analysis_functions.R")

source("full simulations/TD EXIL 2024 measurement error impact/scripts/simulation_interpretation_functions.R")

source("other/support_functions.R")

source("other/performance_metrics.R")

source("parameter estimation/bayesian/load.webexpo.SEG.functions.R")

source("parameter estimation/frequentist/percentile.R")

source("parameter estimation/frequentist/exceedance_fraction.R")

source("data generation/SEG data simulation.R")

source("full simulations/TD EXIL 2024 measurement error impact/scripts/one_scenario_functions.R")


##### SIMULATUON SPACE ####

scenarios <- expand.grid(
  true_exceedance_perc = c(0.1, 1, 3, 7, 10, 25),
  sample_size      = c(3L, 6L , 9L, 12L),
  proportion_censored        = c(0, 0.3, 0.6),
  stringsAsFactors = FALSE)

## fixed parameters

n_sim <- 5000

me_cv <- 0.25

n_iterations_gum = 5000  

true_p95 <- 100




##### ANALYSES ####

simulation_results <- vector("list", length = dim(scenarios)[1])

simulated_data_objects <- vector("list", length = dim(scenarios)[1])


## simulated exposure values

for (i in 1:dim(scenarios)[1]) {
  
  ## parameter vectors due to GSDs
  
  true_gsd <- sample( real_gsds[ real_gsds>=quantile(real_gsds,0.025) & real_gsds<=quantile(real_gsds,0.975)] , n_sim , replace = TRUE )
  
  true_gm <- exp( log(true_p95) - qnorm(0.95)*log(true_gsd) )
  
  oel <- exp( qnorm(1 - scenarios$true_exceedance_perc[i]/100, mean = log(true_gm) , sd = log(true_gsd) ) )
  
  loq <- signif(exp(qnorm(scenarios$proportion_censored[i], mean = log(true_gm) , sd = log(true_gsd) )),3)
  
  
  ## simulated exposure values
  
  simulated_data_objects[[i]] <- data.simulation.seg.gen(  true_gm = true_gm, 
                                                           true_gsd = true_gsd,
                                                           sample_size = scenarios$sample_size[i]  ,
                                                           me_cv_inf = me_cv,
                                                           me_cv_sup = me_cv,
                                                           censor_level = loq )
  
}


#saveRDS(simulated_data_objects, file = "C:/jerome/Dropbox/GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/all_scenarios_GSD_run4_data.RDS")

simulated_data_objects <- readRDS(file = "C:/jerome/Dropbox/GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/all_scenarios_GSD_run4_data.RDS")


## stan models

stan.models.list <- list()

stan.models.list <- augment.stan.models.list(stan.models.list, stan.file=code_seg_informedvar)
stan.models.list <- augment.stan.models.list(stan.models.list, stan.file=code_seg_informedvar_lognormal_mecv)
stan.models.list <- augment.stan.models.list(stan.models.list, stan.file=code_seg_informedvar_lognormal_mecvknown)

compiled.models.list(stan.models.list)


## running the parallel function 

simulation_results <-readRDS("C:/jerome/Dropbox/GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/all_scenarios_GSD_run4f_sim.RDS")


i <-44 
  
  true_gsd <- simulated_data_objects[[i]]$true_gsd
  
  true_gm <- simulated_data_objects[[i]]$true_gm
  
  oel <- exp( qnorm(1 - scenarios$true_exceedance_perc[i]/100, mean = log(true_gm) , sd = log(true_gsd) ) )
  
  simulation_results[[i]] <- parallel.function.s( simulated_data_object = simulated_data_objects[[i]] , me_cv = me_cv , 
                                                  n_iterations_gum = n_iterations_gum , n_sim = n_sim , 
                                                  n_clusters = 24, oel = oel , models.list = stan.models.list)
  
  
  print(i)
  
  print(simulation_results[[i]]$time)
  
  print(Sys.time())
  



## saving the simulation results and the simulated data

#saveRDS(simulation_results, file = "C:/jerome/Dropbox/GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/all_scenarios_GSD_run4a_sim.RDS")
#saveRDS(simulation_results, file = "C:/jerome/Dropbox/GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/all_scenarios_GSD_run4b_sim.RDS")
#saveRDS(simulation_results, file = "C:/jerome/Dropbox/GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/all_scenarios_GSD_run4d_sim.RDS")
#saveRDS(simulation_results, file = "C:/jerome/Dropbox/GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/all_scenarios_GSD_run4e_sim.RDS")
#saveRDS(simulation_results, file = "C:/jerome/Dropbox/GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/all_scenarios_GSD_run4f_sim.RDS")
saveRDS(simulation_results, file = "C:/jerome/Dropbox/GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/all_scenarios_GSD_run4g_sim.RDS")




