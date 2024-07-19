# Script for performing a simulation study with all scenarios : only ten iterations per scenario

# test seccessfully run on 2024-07-17

##### LIBRARIES ####

library(tolerance)
library(parallel)

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
  true_gsd        = c(1.5, 2.5, 3.5),
  true_exceedance_perc = c(0.1, 1, 3, 7, 10, 25),
  sample_size      = c(3L, 6L , 9L, 12L),
  stringsAsFactors = FALSE)

## parameters depending on scenario

scenarios$true_p95 <- 100

scenarios$true_gm <- exp( log(scenarios$true_p95) - qnorm(0.95)*log(scenarios$true_gsd) )

scenarios$oel <- exp( qnorm(1 - scenarios$true_exceedance_perc/100, mean = log(scenarios$true_gm) , sd = log(scenarios$true_gsd) ) )

scenarios$proportion_censored <- 0 # need to adress the issue of "all points censored", or ros log(GSD) == 0

scenarios$loq <- exp(qnorm(scenarios$proportion_censored, mean = log(scenarios$true_gm) , sd = log(scenarios$true_gsd) ))


## fixed parameters

n_sim <- 10

me_cv <- 0.25

n_iterations_gum = 50  

sim_quantile = 0.95



##### ANALYSES ####

simulation_results <- vector("list", length = dim(scenarios)[1])

simulated_data_objects <- vector("list", length = dim(scenarios)[1])


## simulated exposure values

for (i in 1:dim(scenarios)[1]) {
  
  simulated_data_objects[[i]] <- data.simulation.seg( true_gsd = scenarios$true_gsd[i], 
                                                      true_gm = scenarios$true_gm[i], 
                                                      n_sim = n_sim, 
                                                      sample_size = scenarios$sample_size[i] ,
                                                      me_cv_inf = me_cv,
                                                      me_cv_sup = me_cv,
                                                      censor_level = scenarios$loq[i] )
  
}

## running the parallel function


for (i in 1:dim(scenarios)[1]) {
  
  simulation_results[[i]] <- parallel.function( simulated_data_object = simulated_data_objects[[i]] , me_cv = me_cv , 
                                                n_iterations_gum = n_iterations_gum , sim_quantile = sim_quantile , n_sim = n_sim , 
                                                n_clusters = 8, oel = scenarios$oel[i])
  
}


#### INTERPREATION ####

simulation_interpretation_results <- vector("list", length = dim(scenarios)[1])

for (i in 1:dim(scenarios)[1]) {
  
  simulation_interpretation_results[[i]] <- list( rmse = rmse.result( results_one_scenario = simulation_results[[i]] , 
                                                                      true_gm = scenarios$true_gm[i] , true_gsd = scenarios$true_gsd[i] , true_p95 = scenarios$true_p95[i] , 
                                                                      true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                                                  
                                                  precision = precision.result( results_one_scenario = simulation_results[[i]] , 
                                                                              true_gm = scenarios$true_gm[i] , true_gsd = scenarios$true_gsd[i] , true_p95 = scenarios$true_p95[i] , 
                                                                              true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                                                  
                                                  bias = bias.result( results_one_scenario = simulation_results[[i]] , 
                                                                      true_gm = scenarios$true_gm[i] , true_gsd = scenarios$true_gsd[i] , true_p95 = scenarios$true_p95[i] , 
                                                                      true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                                                  
                                                  median_error = median.error.result( results_one_scenario = simulation_results[[i]] , 
                                                                                       true_gm = scenarios$true_gm[i] , true_gsd = scenarios$true_gsd[i] , true_p95 = scenarios$true_p95[i] , 
                                                                                       true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                                                  
                                                  rmsle = rmsle.result( results_one_scenario = simulation_results[[i]] , 
                                                                        true_gm = scenarios$true_gm[i] , true_gsd = scenarios$true_gsd[i] , true_p95 = scenarios$true_p95[i] , 
                                                                        true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                                                  
                                                  mad = mad.result( results_one_scenario = simulation_results[[i]] , 
                                                                    true_gm = scenarios$true_gm[i] , true_gsd = scenarios$true_gsd[i] , true_p95 = scenarios$true_p95[i] , 
                                                                    true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                                                  
                                                  coverage = coverage.result( results_one_scenario = simulation_results[[i]] , 
                                                                             true_p95 = scenarios$true_p95[i] , 
                                                                             true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                                                  
                                                  perc_mistake = perc.mistake.result( results_one_scenario = simulation_results[[i]] , 
                                                                                      true_p95 = scenarios$true_p95[i] , 
                                                                                      true_exceedance_perc = scenarios$true_exceedance_perc[i] ,
                                                                                      oel = scenarios$oel[i] ) )
  
}                     
  
  
                                                  

