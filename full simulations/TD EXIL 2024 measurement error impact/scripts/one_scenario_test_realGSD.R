# Script for performing a simulation study with one scenario - test script


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

##### ANALYSIS 1 no ND####

###### parameters ####

true_exceedance_perc <- 3  ## in percentage

true_p95 <- 100

proportion_censored <- 0

sample_size <- 12

n_sim <- 1000 

me_cv <- 0.25

n_iterations_gum = 5000  

true_gsd <- sample( real_gsds , n_sim , replace = TRUE )

true_gm <- exp( log(true_p95) - qnorm(0.95)*log(true_gsd) )

oel <- exp( qnorm(1 - true_exceedance_perc/100, mean = log(true_gm) , sd = log(true_gsd) ) )

loq <- signif(exp(qnorm(proportion_censored, mean = log(true_gm) , sd = log(true_gsd) )),3)




###### testing sample generation ####    

test_data_sim <- data.simulation.seg.gen(  true_gm = true_gm, 
                                           true_gsd = true_gsd,
                                           sample_size = sample_size  ,
                                           me_cv_inf = me_cv,
                                           me_cv_sup = me_cv,
                                           censor_level = loq )


###### testing parallel function ####  

test_parallel <- parallel.function( simulated_data_object = test_data_sim , me_cv = me_cv , 
                                         n_iterations_gum = n_iterations_gum , n_sim = n_sim , 
                                         n_clusters = 10, oel = oel)

        
        
test_parallel$time

###### testing performance functions ####      
        
test_rmse <- rmse.result( results_one_scenario = test_parallel , 
                          true_gm = NA , true_gsd = NA , true_p95 = true_p95 , 
                          true_exceedance_perc = true_exceedance_perc )
  

test_precision <- precision.result( results_one_scenario = test_parallel , 
                          true_gm = NA , true_gsd = NA , true_p95 = true_p95 , 
                          true_exceedance_perc = true_exceedance_perc )

test_bias <- bias.result( results_one_scenario = test_parallel , 
                                    true_gm = NA , true_gsd = NA , true_p95 = true_p95 , 
                                    true_exceedance_perc = true_exceedance_perc )

test_median_error <- median.error.result( results_one_scenario = test_parallel , 
                          true_gm = NA , true_gsd = NA , true_p95 = true_p95 , 
                          true_exceedance_perc = true_exceedance_perc )

test_rmsle <- rmsle.result( results_one_scenario = test_parallel , 
                          true_gm = NA , true_gsd = NA , true_p95 = true_p95 , 
                          true_exceedance_perc = true_exceedance_perc )

test_mad <- mad.result( results_one_scenario = test_parallel , 
                            true_gm = NA , true_gsd = NA , true_p95 = true_p95 , 
                            true_exceedance_perc = true_exceedance_perc )

test_coverage <- coverage.result( results_one_scenario = test_parallel ,
                                  true_p95 = true_p95 , 
                                  true_exceedance_perc = true_exceedance_perc )


test_perc_mistake <- perc.mistake.result( results_one_scenario = test_parallel ,
                                  true_p95 = true_p95 , 
                                  true_exceedance_perc = true_exceedance_perc,
                                  oel=oel)


##### ANALYSIS 2 ND ####

###### parameters ####

true_exceedance_perc <- 5.1  ## in percentage

true_p95 <- 100

proportion_censored <- 0.30

sample_size <- 6

n_sim <- 1000 

me_cv <- 0.25

n_iterations_gum = 5000  

sim_quantile = 0.975

true_gsd <- sample( real_gsds , n_sim , replace = TRUE )

true_gm <- exp( log(true_p95) - qnorm(0.95)*log(true_gsd) )

oel <- exp( qnorm(1 - true_exceedance_perc/100, mean = log(true_gm) , sd = log(true_gsd) ) )

loq <- signif(exp(qnorm(proportion_censored, mean = log(true_gm) , sd = log(true_gsd) )),3)




###### testing sample generation ####    

test_data_sim2 <- data.simulation.seg.gen(  true_gm = true_gm, 
                                           true_gsd = true_gsd,
                                           sample_size = sample_size  ,
                                           me_cv_inf = me_cv,
                                           me_cv_sup = me_cv,
                                           censor_level = loq )

###### testing parallel function ####  

test_parallel2 <- parallel.function( simulated_data_object = test_data_sim2 , me_cv = me_cv , 
                                    n_iterations_gum = n_iterations_gum , n_sim = n_sim , 
                                    n_clusters = 10, oel = oel)



test_parallel2$time

###### testing performance functions ####      

test_rmse2 <- rmse.result( results_one_scenario = test_parallel2 , 
                          true_gm = NA , true_gsd = NA , true_p95 = true_p95 , 
                          true_exceedance_perc = true_exceedance_perc )


test_precision2 <- precision.result( results_one_scenario = test_parallel2 , 
                                    true_gm = NA , true_gsd = NA , true_p95 = true_p95 , 
                                    true_exceedance_perc = true_exceedance_perc )

test_bias2 <- bias.result( results_one_scenario = test_parallel2 , 
                          true_gm = NA , true_gsd = NA , true_p95 = true_p95 , 
                          true_exceedance_perc = true_exceedance_perc )

test_median_error2 <- median.error.result( results_one_scenario = test_parallel2 , 
                                          true_gm = NA , true_gsd = NA , true_p95 = true_p95 , 
                                          true_exceedance_perc = true_exceedance_perc )

test_rmsle2 <- rmsle.result( results_one_scenario = test_parallel2 , 
                            true_gm = NA , true_gsd = NA , true_p95 = true_p95 , 
                            true_exceedance_perc = true_exceedance_perc )

test_mad2 <- mad.result( results_one_scenario = test_parallel2 , 
                        true_gm = NA , true_gsd = NA , true_p95 = true_p95 , 
                        true_exceedance_perc = true_exceedance_perc )

test_coverage2 <- coverage.result( results_one_scenario = test_parallel2 ,
                                  true_p95 = true_p95 , 
                                  true_exceedance_perc = true_exceedance_perc )


test_perc_mistake2 <- perc.mistake.result( results_one_scenario = test_parallel2 ,
                                          true_p95 = true_p95 , 
                                          true_exceedance_perc = true_exceedance_perc,
                                          oel=oel)

