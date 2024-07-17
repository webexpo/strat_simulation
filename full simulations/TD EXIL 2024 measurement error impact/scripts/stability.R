# Script for assessing the stability of the simulation 


##### LIBRARIES ####


##### DATA ####

test_repeatability <- readRDS("full simulations/TD EXIL 2024 measurement error impact/data/test_repeatability1.RDS")

##### SCRIPTS ####

source("full simulations/TD EXIL 2024 measurement error impact/scripts/simulation_interpretation.R")
source("other/performance_metrics.R")

##### TEST 1 ####

### parameters

true_gsd <- 2.5

true_exceedance_perc <- 2.5  ## in percentage

true_p95 <- 100

true_gm <- exp( log(true_p95) - qnorm(0.95)*log(true_gsd) )

oel <- exp( qnorm(1 - true_exceedance_perc/100, mean = log(true_gm) , sd = log(true_gsd) ) )

sample_size <- 6

n_sim <- 5000

me_cv <- 0.25

n_iterations_gum = 5000

sim_quantile = 0.95

### analysis

test_rmse <- rmse.result( results_one_scenario = test_repeatability[[1]] , 
                          true_gm = true_gm , true_gsd = true_gsd , true_p95 = true_p95 , 
                          true_exceedance_perc = true_exceedance_perc )


rmse_array <- array( dim=c(8,4,5) )

for (i in 1:5) { rmse_array[,,i] <- as.matrix(rmse.result( results_one_scenario = test_repeatability[[i]] , 
                                                           true_gm = true_gm , true_gsd = true_gsd , true_p95 = true_p95 , 
                                                           true_exceedance_perc = true_exceedance_perc )[,2:5]) }

result_meansd <-apply( rmse_array , c(1,2) , function(x) { paste( signif(mean(x),2) , 
                                                  "(", signif(100*sd(x)/mean(x),2),")", sep="")})
                                                  
                                                  
result_minmax <-apply( rmse_array , c(1,2) , function(x) { paste("[", signif(min(x),2) , "-" , signif(max(x),2) , "]",
                                                  sep="") } )

