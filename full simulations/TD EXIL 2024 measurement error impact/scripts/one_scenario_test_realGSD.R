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

true_exceedance_perc <- 1  ## in percentage

true_p95 <- 100

proportion_censored <- 0

sample_size <- 12

n_sim <- 5000 

me_cv <- 0.25

n_iterations_gum = 5000  

sim_quantile = 0.975

#true_gsd <- rep( 2.5 , n_sim )

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
                                         n_iterations_gum = n_iterations_gum , sim_quantile = sim_quantile , n_sim = n_sim , 
                                         n_clusters = 2, oel = oel)

        
        
test_parallel$time

# saving the results on dropbox

saveRDS(list(data = test_data_sim,
             sim = test_parallel),"F:/Dropbox/GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/test_one_scneario_realGSD1.RDS")
    
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

true_exceedance_perc <- 1  ## in percentage

true_p95 <- 100

proportion_censored <- 0.43

sample_size <- 12

n_sim <- 5000 

me_cv <- 0.25

n_iterations_gum = 5000  

sim_quantile = 0.975

#true_gsd <- rep( 2.5 , n_sim )

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
                                    n_iterations_gum = n_iterations_gum , sim_quantile = sim_quantile , n_sim = n_sim , 
                                    n_clusters = 2, oel = oel)



test_parallel2$time

saveRDS(list(data = test_data_sim2,
             sim = test_parallel2),"F:/Dropbox/GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/test_one_scneario_realGSD2.RDS")


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

##### ANALYSIS 3 no ND GSD 2.5####

###### parameters ####

true_exceedance_perc <- 1  ## in percentage

true_p95 <- 100

proportion_censored <- 0

sample_size <- 12

n_sim <- 5000 

me_cv <- 0.25

n_iterations_gum = 5000  

sim_quantile = 0.975

true_gsd <- rep( 2.5 , n_sim )

#true_gsd <- sample( real_gsds , n_sim , replace = TRUE )

true_gm <- exp( log(true_p95) - qnorm(0.95)*log(true_gsd) )

oel <- exp( qnorm(1 - true_exceedance_perc/100, mean = log(true_gm) , sd = log(true_gsd) ) )

loq <- signif(exp(qnorm(proportion_censored, mean = log(true_gm) , sd = log(true_gsd) )),3)




###### testing sample generation ####    

test_data_sim3 <- data.simulation.seg.gen(  true_gm = true_gm, 
                                           true_gsd = true_gsd,
                                           sample_size = sample_size  ,
                                           me_cv_inf = me_cv,
                                           me_cv_sup = me_cv,
                                           censor_level = loq )


###### testing parallel function ####  

test_parallel3 <- parallel.function( simulated_data_object = test_data_sim3 , me_cv = me_cv , 
                                    n_iterations_gum = n_iterations_gum , sim_quantile = sim_quantile , n_sim = n_sim , 
                                    n_clusters = 8, oel = oel)



test_parallel3$time

saveRDS(list(data = test_data_sim3,
             sim = test_parallel3),"F:/Dropbox/GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/test_one_scneario_realGSD3.RDS")

###### testing performance functions ####      

test_rmse3 <- rmse.result( results_one_scenario = test_parallel3 , 
                          true_gm = NA , true_gsd = NA , true_p95 = true_p95 , 
                          true_exceedance_perc = true_exceedance_perc )


test_precision3 <- precision.result( results_one_scenario = test_parallel3 , 
                                    true_gm = NA , true_gsd = NA , true_p95 = true_p95 , 
                                    true_exceedance_perc = true_exceedance_perc )

test_bias3 <- bias.result( results_one_scenario = test_parallel3 , 
                          true_gm = NA , true_gsd = NA , true_p95 = true_p95 , 
                          true_exceedance_perc = true_exceedance_perc )

test_median_error3 <- median.error.result( results_one_scenario = test_parallel3 , 
                                          true_gm = NA , true_gsd = NA , true_p95 = true_p95 , 
                                          true_exceedance_perc = true_exceedance_perc )

test_rmsle3 <- rmsle.result( results_one_scenario = test_parallel3 , 
                            true_gm = NA , true_gsd = NA , true_p95 = true_p95 , 
                            true_exceedance_perc = true_exceedance_perc )

test_mad3 <- mad.result( results_one_scenario = test_parallel3 , 
                        true_gm = NA , true_gsd = NA , true_p95 = true_p95 , 
                        true_exceedance_perc = true_exceedance_perc )

test_coverage3 <- coverage.result( results_one_scenario = test_parallel3 ,
                                  true_p95 = true_p95 , 
                                  true_exceedance_perc = true_exceedance_perc )


test_perc_mistake3 <- perc.mistake.result( results_one_scenario = test_parallel3 ,
                                          true_p95 = true_p95 , 
                                          true_exceedance_perc = true_exceedance_perc,
                                          oel=oel)


##### ANALYSIS 4 ND GSD 2.5####

###### parameters ####

true_exceedance_perc <- 1  ## in percentage

true_p95 <- 100

proportion_censored <- 0.43

sample_size <- 12

n_sim <- 5000 

me_cv <- 0.25

n_iterations_gum = 5000  

sim_quantile = 0.975

true_gsd <- rep( 2.5 , n_sim )

#true_gsd <- sample( real_gsds , n_sim , replace = TRUE )

true_gm <- exp( log(true_p95) - qnorm(0.95)*log(true_gsd) )

oel <- exp( qnorm(1 - true_exceedance_perc/100, mean = log(true_gm) , sd = log(true_gsd) ) )

loq <- signif(exp(qnorm(proportion_censored, mean = log(true_gm) , sd = log(true_gsd) )),3)




###### testing sample generation ####    

test_data_sim4 <- data.simulation.seg.gen(  true_gm = true_gm, 
                                            true_gsd = true_gsd,
                                            sample_size = sample_size  ,
                                            me_cv_inf = me_cv,
                                            me_cv_sup = me_cv,
                                            censor_level = loq )

###### testing parallel function ####  

test_parallel4 <- parallel.function( simulated_data_object = test_data_sim4 , me_cv = me_cv , 
                                     n_iterations_gum = n_iterations_gum , sim_quantile = sim_quantile , n_sim = n_sim , 
                                     n_clusters = 8, oel = oel)



test_parallel4$time

saveRDS(list(data = test_data_sim4,
             sim = test_parallel4),"F:/Dropbox/GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/test_one_scneario_realGSD4.RDS")

###### testing performance functions ####      

test_rmse4 <- rmse.result( results_one_scenario = test_parallel4 , 
                           true_gm = NA , true_gsd = NA , true_p95 = true_p95 , 
                           true_exceedance_perc = true_exceedance_perc )


test_precision4 <- precision.result( results_one_scenario = test_parallel4 , 
                                     true_gm = NA , true_gsd = NA , true_p95 = true_p95 , 
                                     true_exceedance_perc = true_exceedance_perc )

test_bias4 <- bias.result( results_one_scenario = test_parallel4 , 
                           true_gm = NA , true_gsd = NA , true_p95 = true_p95 , 
                           true_exceedance_perc = true_exceedance_perc )

test_median_error4 <- median.error.result( results_one_scenario = test_parallel4 , 
                                           true_gm = NA , true_gsd = NA , true_p95 = true_p95 , 
                                           true_exceedance_perc = true_exceedance_perc )

test_rmsle4 <- rmsle.result( results_one_scenario = test_parallel4 , 
                             true_gm = NA , true_gsd = NA , true_p95 = true_p95 , 
                             true_exceedance_perc = true_exceedance_perc )

test_mad4 <- mad.result( results_one_scenario = test_parallel4 , 
                         true_gm = NA , true_gsd = NA , true_p95 = true_p95 , 
                         true_exceedance_perc = true_exceedance_perc )

test_coverage4 <- coverage.result( results_one_scenario = test_parallel4 ,
                                   true_p95 = true_p95 , 
                                   true_exceedance_perc = true_exceedance_perc )


test_perc_mistake4 <- perc.mistake.result( results_one_scenario = test_parallel4 ,
                                           true_p95 = true_p95 , 
                                           true_exceedance_perc = true_exceedance_perc,
                                           oel=oel)

