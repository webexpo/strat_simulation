# Script for the interpretation of the simulations

# testing existing simulation and comparison across runs

# validation against the initial pilot GUM measurement error_interpretation.R

##### LIBRARIES ####

##### Scripts ####

source("full simulations/TD EXIL 2024 measurement error impact/scripts/simulation_interpretation_functions.R")

source("other/performance_metrics.R")

##### DATA ####

init_path <- "F:/Dropbox/"
init_path <- "C:/jerome/Dropbox/"


load(file = paste( init_path ,"GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/all_scenarios_run3a_sim.RDS", sep=""))
run3a <- simulation_results

load(file = paste( init_path ,"GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/all_scenarios_run3b_sim.RDS", sep=""))
run3b <- simulation_results

run1 <- readRDS(file = paste( init_path ,"GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/all_scenarios_run1_d.RDS", sep=""))

run2 <- load(file = paste( init_path ,"GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/all_scenarios_run2_sim.RDS", sep=""))
run2 <- simulation_results


##### SCENARIOS ####

scenarios <- expand.grid(
  true_gsd        = c(1.5, 2.5, 3.5),
  true_exceedance_perc = c(0.1, 1, 3, 7, 10, 25),
  sample_size      = c(3L, 6L , 9L, 12L),
  stringsAsFactors = FALSE)

scenarios$true_p95 <- 100

scenarios$true_gm <- exp( log(scenarios$true_p95) - qnorm(0.95)*log(scenarios$true_gsd) )

scenarios$oel <- exp( qnorm(1 - scenarios$true_exceedance_perc/100, mean = log(scenarios$true_gm) , sd = log(scenarios$true_gsd) ) )


##### interpreting scenario 23 across all runs : 1% exceedanace, gsd=2.5, sample size = 6

i <- 23

## most recent functions : run 3

run3a_inter <- list( rmse = rmse.result( results_one_scenario = run3a[[i]] , 
                                          true_gm = scenarios$true_gm[i] ,
                                          true_gsd = scenarios$true_gsd[i] ,
                                          true_p95 = scenarios$true_p95[i] ,
                                          true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                                                  
                    precision = precision.result( results_one_scenario = run3a[[i]] , 
                                                  true_gm = scenarios$true_gm[i] , 
                                                  true_gsd = scenarios$true_gsd[i] ,
                                                  true_p95 = scenarios$true_p95[i] ,
                                                  true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                                                  
                    bias = bias.result( results_one_scenario = run3a[[i]] , 
                                        true_gm = scenarios$true_gm[i] , 
                                        true_gsd = scenarios$true_gsd[i] ,
                                        true_p95 = scenarios$true_p95[i] ,
                                        true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                    
                    median_error = median.error.result( results_one_scenario = run3a[[i]] , 
                                                        true_gm = scenarios$true_gm[i] , 
                                                        true_gsd = scenarios$true_gsd[i] ,
                                                        true_p95 = scenarios$true_p95[i] , 
                                                        true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                                                  
                    rmsle = rmsle.result( results_one_scenario = run3a[[i]] , 
                                          true_gm = scenarios$true_gm[i] , 
                                          true_gsd = scenarios$true_gsd[i] , 
                                          true_p95 = scenarios$true_p95[i] , 
                                          true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                                                  
                    mad = mad.result( results_one_scenario = run3a[[i]] , 
                                      true_gm = scenarios$true_gm[i] , 
                                      true_gsd = scenarios$true_gsd[i] , 
                                      true_p95 = scenarios$true_p95[i] , 
                                      true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                                                  
                    coverage = coverage.result( results_one_scenario = run3a[[i]] , 
                                                true_p95 = scenarios$true_p95[i] , 
                                                true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                                                  
                    perc_mistake = perc.mistake.result( results_one_scenario = run3a[[i]] , 
                                                        true_p95 = scenarios$true_p95[i] , 
                                                        true_exceedance_perc = scenarios$true_exceedance_perc[i] ,
                                                        oel = scenarios$oel[i] ) )
  
                     
run3b_inter <- list( rmse = rmse.result( results_one_scenario = run3b[[i]] , 
                                         true_gm = scenarios$true_gm[i] ,
                                         true_gsd = scenarios$true_gsd[i] ,
                                         true_p95 = scenarios$true_p95[i] ,
                                         true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                     
                     precision = precision.result( results_one_scenario = run3b[[i]] , 
                                                   true_gm = scenarios$true_gm[i] , 
                                                   true_gsd = scenarios$true_gsd[i] ,
                                                   true_p95 = scenarios$true_p95[i] ,
                                                   true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                     
                     bias = bias.result( results_one_scenario = run3b[[i]] , 
                                         true_gm = scenarios$true_gm[i] , 
                                         true_gsd = scenarios$true_gsd[i] ,
                                         true_p95 = scenarios$true_p95[i] ,
                                         true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                     
                     median_error = median.error.result( results_one_scenario = run3b[[i]] , 
                                                         true_gm = scenarios$true_gm[i] , 
                                                         true_gsd = scenarios$true_gsd[i] ,
                                                         true_p95 = scenarios$true_p95[i] , 
                                                         true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                     
                     rmsle = rmsle.result( results_one_scenario = run3b[[i]] , 
                                           true_gm = scenarios$true_gm[i] , 
                                           true_gsd = scenarios$true_gsd[i] , 
                                           true_p95 = scenarios$true_p95[i] , 
                                           true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                     
                     mad = mad.result( results_one_scenario = run3b[[i]] , 
                                       true_gm = scenarios$true_gm[i] , 
                                       true_gsd = scenarios$true_gsd[i] , 
                                       true_p95 = scenarios$true_p95[i] , 
                                       true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                     
                     coverage = coverage.result( results_one_scenario = run3b[[i]] , 
                                                 true_p95 = scenarios$true_p95[i] , 
                                                 true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                     
                     perc_mistake = perc.mistake.result( results_one_scenario = run3b[[i]] , 
                                                         true_p95 = scenarios$true_p95[i] , 
                                                         true_exceedance_perc = scenarios$true_exceedance_perc[i] ,
                                                         oel = scenarios$oel[i] ) )


## run 2 : only quantile for the GUM approach is 97.5%

# tweak to have the same format at run3

run2_tweaked <- run2

run2_tweaked[[i]] <- run2[[i]]

run2_tweaked[[i]]$array <- array(NA, dim = c( 8 , 11 , 50000 ) )

for (j in 1:dim(run2[[i]]$array)[3]) {
  
  run2_tweaked[[i]]$array[ , 1:8 , j] <- run2[[i]]$array[ , , j]
  
  run2_tweaked[[i]]$array[ , 9:11 , j] <- run2[[i]]$array[ , c(8,8,8) , j]
  
}



run2_inter <- list( rmse = rmse.result( results_one_scenario = run2_tweaked[[i]] , 
                                         true_gm = scenarios$true_gm[i] ,
                                         true_gsd = scenarios$true_gsd[i] ,
                                         true_p95 = scenarios$true_p95[i] ,
                                         true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                     
                     precision = precision.result( results_one_scenario = run2_tweaked[[i]] , 
                                                   true_gm = scenarios$true_gm[i] , 
                                                   true_gsd = scenarios$true_gsd[i] ,
                                                   true_p95 = scenarios$true_p95[i] ,
                                                   true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                     
                     bias = bias.result( results_one_scenario = run2_tweaked[[i]] , 
                                         true_gm = scenarios$true_gm[i] , 
                                         true_gsd = scenarios$true_gsd[i] ,
                                         true_p95 = scenarios$true_p95[i] ,
                                         true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                     
                     median_error = median.error.result( results_one_scenario = run2_tweaked[[i]] , 
                                                         true_gm = scenarios$true_gm[i] , 
                                                         true_gsd = scenarios$true_gsd[i] ,
                                                         true_p95 = scenarios$true_p95[i] , 
                                                         true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                     
                     rmsle = rmsle.result( results_one_scenario = run2_tweaked[[i]] , 
                                           true_gm = scenarios$true_gm[i] , 
                                           true_gsd = scenarios$true_gsd[i] , 
                                           true_p95 = scenarios$true_p95[i] , 
                                           true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                     
                     mad = mad.result( results_one_scenario = run2_tweaked[[i]] , 
                                       true_gm = scenarios$true_gm[i] , 
                                       true_gsd = scenarios$true_gsd[i] , 
                                       true_p95 = scenarios$true_p95[i] , 
                                       true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                     
                     coverage = coverage.result( results_one_scenario = run2_tweaked[[i]] , 
                                                 true_p95 = scenarios$true_p95[i] , 
                                                 true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                     
                     perc_mistake = perc.mistake.result( results_one_scenario = run2_tweaked[[i]] , 
                                                         true_p95 = scenarios$true_p95[i] , 
                                                         true_exceedance_perc = scenarios$true_exceedance_perc[i] ,
                                                         oel = scenarios$oel[i] ) )



## run 1 : only quantile for the GUM approach is 95%


run1_tweaked <- run1

run1_tweaked[[i]] <- run1[[i]]

run1_tweaked[[i]]$array <- array(NA, dim = c( 8 , 11 , 50000 ) )

for (j in 1:dim(run1[[i]]$array)[3]) {
  
  run1_tweaked[[i]]$array[ , 1:8 , j] <- run1[[i]]$array[ , , j]
  
  run1_tweaked[[i]]$array[ , 9:11 , j] <- run1[[i]]$array[ , c(8,8,8) , j]
  
}



run1_inter <- list( rmse = rmse.result( results_one_scenario = run1_tweaked[[i]] , 
                                        true_gm = scenarios$true_gm[i] ,
                                        true_gsd = scenarios$true_gsd[i] ,
                                        true_p95 = scenarios$true_p95[i] ,
                                        true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                    
                    precision = precision.result( results_one_scenario = run1_tweaked[[i]] , 
                                                  true_gm = scenarios$true_gm[i] , 
                                                  true_gsd = scenarios$true_gsd[i] ,
                                                  true_p95 = scenarios$true_p95[i] ,
                                                  true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                    
                    bias = bias.result( results_one_scenario = run1_tweaked[[i]] , 
                                        true_gm = scenarios$true_gm[i] , 
                                        true_gsd = scenarios$true_gsd[i] ,
                                        true_p95 = scenarios$true_p95[i] ,
                                        true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                    
                    median_error = median.error.result( results_one_scenario = run1_tweaked[[i]] , 
                                                        true_gm = scenarios$true_gm[i] , 
                                                        true_gsd = scenarios$true_gsd[i] ,
                                                        true_p95 = scenarios$true_p95[i] , 
                                                        true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                    
                    rmsle = rmsle.result( results_one_scenario = run1_tweaked[[i]] , 
                                          true_gm = scenarios$true_gm[i] , 
                                          true_gsd = scenarios$true_gsd[i] , 
                                          true_p95 = scenarios$true_p95[i] , 
                                          true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                    
                    mad = mad.result( results_one_scenario = run1_tweaked[[i]] , 
                                      true_gm = scenarios$true_gm[i] , 
                                      true_gsd = scenarios$true_gsd[i] , 
                                      true_p95 = scenarios$true_p95[i] , 
                                      true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                    
                    coverage = coverage.result( results_one_scenario = run1_tweaked[[i]] , 
                                                true_p95 = scenarios$true_p95[i] , 
                                                true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                    
                    perc_mistake = perc.mistake.result( results_one_scenario = run1_tweaked[[i]] , 
                                                        true_p95 = scenarios$true_p95[i] , 
                                                        true_exceedance_perc = scenarios$true_exceedance_perc[i] ,
                                                        oel = scenarios$oel[i] ) )

# RMSE and bias table for all runs - P95 (comparable for 7 methods)

rmse_comp_p95 <- data.frame( run1 = run1_inter$rmse$p95[1:7] , 
                             run2 = run2_inter$rmse$p95[1:7] , 
                             run3a = run3a_inter$rmse$p95[1:7] , 
                             run3b = run3b_inter$rmse$p95[1:7] )


bias_comp_p95 <- data.frame( run1 = run1_inter$bias$p95[1:7] , 
                             run2 = run2_inter$bias$p95[1:7] , 
                             run3a = run3a_inter$bias$p95[1:7] , 
                             run3b = run3b_inter$bias$p95[1:7] )

     # Results :  fair reproductibility and validated against the pilot markdown


# coverage

coverage_comp_ucl70 <- data.frame( run1 = run1_inter$coverage$p95_ucl70[1:7] , 
                                 run2 = run2_inter$coverage$p95_ucl70[1:7] , 
                                 run3a = run3a_inter$coverage$p95_ucl70[1:7] , 
                                 run3b = run3b_inter$coverage$p95_ucl70[1:7] )

        # Results :  fair reproductibility and validated against the pilot markdown

# perc mistake

perc_mistake_comp <- data.frame( run1 = run1_inter$perc_mistake$exceedance_ucl70[1:7] , 
                                     run2 = run2_inter$perc_mistake$exceedance_ucl70[1:7] , 
                                     run3a = run3a_inter$perc_mistake$exceedance_ucl70[1:7] , 
                                     run3b = run3b_inter$perc_mistake$exceedance_ucl70[1:7] )

# Results :  fair reproductibility and validated against the pilot markdown. Warning : perc mistake for  F methods not good for runs 1 and 2 (OEL was stuck)

