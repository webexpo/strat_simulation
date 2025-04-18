# Script for performing a simulation study with one scenario - test script


##### LIBRARIES ####

library(tolerance)
library(parallel)


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



##### PARAMETERS ####

true_gsd <- 2.5

true_exceedance_perc <- 5.1  ## in percentage

true_p95 <- 100

true_gm <- exp( log(true_p95) - qnorm(0.95)*log(true_gsd) )

oel <- exp( qnorm(1 - true_exceedance_perc/100, mean = log(true_gm) , sd = log(true_gsd) ) )

proportion_censored <- 0.3

loq <- signif(exp(qnorm(proportion_censored, mean = log(true_gm) , sd = log(true_gsd) )),3)

sample_size <- 6

n_sim <- 1000 

me_cv <- 0.25

n_iterations_gum = 5000  

##### ANALYSES ####

###### testing sample generation ####    

test_data_sim <- data.simulation.seg( true_gsd = true_gsd, 
                                      true_gm = true_gm, 
                                      n_sim = n_sim, 
                                      sample_size = sample_size ,
                                      me_cv_inf = me_cv,
                                      me_cv_sup = me_cv,
                                      censor_level = loq )


###### testing ith.pair function ####  

test_ith.pair <- ithpair.function( index = 1 , simulated_data_object = test_data_sim , me_cv = me_cv , 
                                   n_iterations_gum = n_iterations_gum , oel = oel)

test_ith.pair

###### testing parallel function ####  

test_parallel <- parallel.function( simulated_data_object = test_data_sim , me_cv = me_cv , 
                                    n_iterations_gum = n_iterations_gum , n_sim = n_sim , 
                                    n_clusters = 10, oel = rep(oel,n_sim))



test_parallel$time

#data interpretation : mean over the simulation

apply(test_parallel$array, c(1,2), mean)


###### testing performance functions ####      

test_rmse <- rmse.result( results_one_scenario = test_parallel , 
                          true_gm = true_gm , true_gsd = true_gsd , true_p95 = true_p95 , 
                          true_exceedance_perc = true_exceedance_perc )


test_precision <- precision.result( results_one_scenario = test_parallel , 
                                    true_gm = true_gm , true_gsd = true_gsd , true_p95 = true_p95 , 
                                    true_exceedance_perc = true_exceedance_perc )

test_bias <- bias.result( results_one_scenario = test_parallel , 
                          true_gm = true_gm , true_gsd = true_gsd , true_p95 = true_p95 , 
                          true_exceedance_perc = true_exceedance_perc )

test_median_error <- median.error.result( results_one_scenario = test_parallel , 
                                          true_gm = true_gm , true_gsd = true_gsd , true_p95 = true_p95 , 
                                          true_exceedance_perc = true_exceedance_perc )

test_rmsle <- rmsle.result( results_one_scenario = test_parallel , 
                            true_gm = true_gm , true_gsd = true_gsd , true_p95 = true_p95 , 
                            true_exceedance_perc = true_exceedance_perc )

test_mad <- mad.result( results_one_scenario = test_parallel , 
                        true_gm = true_gm , true_gsd = true_gsd , true_p95 = true_p95 , 
                        true_exceedance_perc = true_exceedance_perc )

test_coverage <- coverage.result( results_one_scenario = test_parallel ,
                                  true_p95 = true_p95 , 
                                  true_exceedance_perc = true_exceedance_perc )


test_perc_mistake <- perc.mistake.result( results_one_scenario = test_parallel ,
                                          true_p95 = true_p95 , 
                                          true_exceedance_perc = true_exceedance_perc,
                                          oel=oel)


