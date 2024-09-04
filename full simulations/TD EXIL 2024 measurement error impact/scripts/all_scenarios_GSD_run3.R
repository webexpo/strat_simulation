# Script for performing a simulation study with real GSD : 5000 iterations per scenario

#  Because of errors, see run 1 script, will be fully run in STAN

#  Here : all scenarios in STAN 

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


#saveRDS(simulated_data_objects, file = "C:/jerome/Dropbox/GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/all_scenarios_GSD_run3_data.RDS")

simulated_data_objects <- readRDS( file = "F:/Dropbox/GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/all_scenarios_GSD_run3_sim.RDS")


## stan models

stan.models.list <- list()

stan.models.list <- augment.stan.models.list(stan.models.list, stan.file=code_seg_informedvar)
stan.models.list <- augment.stan.models.list(stan.models.list, stan.file=code_seg_informedvar_lognormal_mecv)
stan.models.list <- augment.stan.models.list(stan.models.list, stan.file=code_seg_informedvar_lognormal_mecvknown)

compiled.models.list(stan.models.list)


## running the parallel function 

start_time <- Sys.time()

for (i in 1:dim(scenarios)[1]) {
  
  #for (i in 5:14) {  
  
  true_gsd <- simulated_data_objects[[i]]$true_gsd
  
  true_gm <- simulated_data_objects[[i]]$true_gm
  
  oel <- exp( qnorm(1 - scenarios$true_exceedance_perc[i]/100, mean = log(true_gm) , sd = log(true_gsd) ) )
  
  simulation_results[[i]] <- parallel.function.s( simulated_data_object = simulated_data_objects[[i]] , me_cv = me_cv , 
                                                n_iterations_gum = n_iterations_gum , n_sim = n_sim , 
                                                n_clusters = 6, oel = oel , models.list = stan.models.list)
  
  
  print(i)
  
  print(simulation_results[[i]]$time)
  
}

end_time <- Sys.time()
mytime <- end_time - start_time 
mytime

## saving the simulation results and the simulated data

saveRDS(simulation_results, file = "C:/jerome/Dropbox/GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/all_scenarios_GSD_run3_sim.RDS")




## debug


i <- 5

true_gsd <- simulated_data_objects[[i]]$true_gsd

true_gm <- simulated_data_objects[[i]]$true_gm

oel <- exp( qnorm(1 - scenarios$true_exceedance_perc[i]/100, mean = log(true_gm) , sd = log(true_gsd) ) )


my_X <- vector("list", length = n_sim)

for ( j in 1:n_sim ) { my_X[[j]] <- list( index = j,
                                          oel = oel[j]) }


simulation_result_nonparallel <- vector("list", length = n_sim)

for ( j in 1:n_sim ) { simulation_result_nonparallel[[j]] <- ithpair.function.me.b.s(my_X[[j]]$index, 
                                                                               simulated_data_object = simulated_data_objects[[i]] , 
                                                                               me_cv = me_cv , 
                                                                               #n_iterations_gum = n_iterations_gum , 
                                                                               oel = my_X[[j]]$oel,
                                                                               models.list = stan.models.list)

print(j)

}

4195

true <-simulated_data_objects[[i]]$true[,1406] 

observed <- simulated_data_objects[[i]]$observed[,1406]

my_X[[1406]]$oel

ithpair.function.s(my_X[[1407]]$index, 
                   simulated_data_object = simulated_data_objects[[i]] , 
                   me_cv = me_cv , 
                   n_iterations_gum = n_iterations_gum , 
                   oel = my_X[[j]]$oel,
                   models.list = stan.models.list)


simulation_result_nonparallel <- lapply(X = my_X , function(x){ ithpair.function.s(x$index, 
                                                                                       simulated_data_object = simulated_data_objects[[i]] , 
                                                                                       me_cv = me_cv , 
                                                                                       n_iterations_gum = n_iterations_gum , 
                                                                                       oel = x$oel,
                                                                                       models.list = stan.models.list) } ) 



SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
Chain 1: Rejecting initial value:
  Chain 1:   Error evaluating the log probability at the initial value.
Chain 1: Exception: lognormal_lpdf: Random variable is -0.019719, but must be nonnegative! (in 'string', line 43, column 4 to column 37)
Chain 1: 
Chain 1: Initialization between (-2, 2) failed after 1 attempts. 
Chain 1:  Try specifying initial values, reducing ranges of constrained values, or reparameterizing the model.
[1] "Error : Initialization failed."
[1] "error occurred during calling the sampler; sampling not done"
Stan model 'anon_model' does not contain samples.
Stan model 'anon_model' does not contain samples.


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
  
  
                                                  

