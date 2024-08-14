# Script for performing a simulation study with real GSD : 5000 iterations per scenario

# test successfully run on 2024-07-17, then after ND and real GSD update, on 2024-07-31 

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

for (i in 68:dim(scenarios)[1]) { #issue with 37,43,67

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
  
  
  
 ## running the parallel function


    simulation_results[[i]] <- parallel.function( simulated_data_object = simulated_data_objects[[i]] , me_cv = me_cv , 
                                                  n_iterations_gum = n_iterations_gum , n_sim = n_sim , 
                                                  n_clusters = 18, oel = oel)
    
    
    print(i)
    
    print(simulation_results[[i]]$time)
  
}


saveRDS( simulation_results , "C:/jerome/dropbox/temp/simulation_results_real_gsd_1.RDS")
saveRDS( simulated_data_objects , "C:/jerome/dropbox/temp/simulated_data_real_gsd_1.RDS")



#### INTERPRETATION ####

simulation_interpretation_results <- vector("list", length = dim(scenarios)[1])

for (i in 1:dim(scenarios)[1]) {
  
  simulation_interpretation_results[[i]] <- list( rmse = rmse.result( results_one_scenario = simulation_results[[i]] , 
                                                                      true_gm = NA , true_gsd = NA , true_p95 = true_p95 , 
                                                                      true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                                                  
                                                  precision = precision.result( results_one_scenario = simulation_results[[i]] , 
                                                                                true_gm = NA , true_gsd = NA , true_p95 = true_p95 , 
                                                                                true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                                                  
                                                  bias = bias.result( results_one_scenario = simulation_results[[i]] , 
                                                                      true_gm = NA , true_gsd = NA , true_p95 = true_p95 , 
                                                                      true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                                                  
                                                  median_error = median.error.result( results_one_scenario = simulation_results[[i]] , 
                                                                                      true_gm = NA , true_gsd = NA , true_p95 = true_p95 , 
                                                                                      true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                                                  
                                                  rmsle = rmsle.result( results_one_scenario = simulation_results[[i]] , 
                                                                        true_gm = NA , true_gsd = NA , true_p95 = true_p95 , 
                                                                        true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                                                  
                                                  mad = mad.result( results_one_scenario = simulation_results[[i]] , 
                                                                    true_gm = NA , true_gsd = NA , true_p95 = true_p95 , 
                                                                    true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                                                  
                                                  coverage = coverage.result( results_one_scenario = simulation_results[[i]] ,
                                                                              true_p95 = true_p95 , 
                                                                              true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                                                  
                                                  perc_mistake = perc.mistake.result( results_one_scenario = simulation_results[[i]] ,
                                                                                      true_p95 = true_p95 , 
                                                                                      true_exceedance_perc = scenarios$true_exceedance_perc[i],
                                                                                      oel=oel ) )
  
}                     
  
  
                                                  

### debug attempts for the problematic scenarios : low exceedance and censoring, only patterns seems higher OELs
    
    #param and data
    
    n_sim <- 5000
    
    me_cv <- 0.25
    
    n_iterations_gum = 5000  
    
    true_p95 <- 100
    
    scenarios <- expand.grid(
      true_exceedance_perc = c(0.1, 1, 3, 7, 10, 25),
      sample_size      = c(3L, 6L , 9L, 12L),
      proportion_censored        = c(0, 0.3, 0.6),
      stringsAsFactors = FALSE)
    
    simulated_data_objects <- readRDS("F:/Dropbox/temp/mytest.RDS")
    
    
    ## parameter vectors due to GSDs i=31 prob
    
    i <- 31
    
    true_gsd <- simulated_data_objects$true_gsd
    
    true_gm <- simulated_data_objects$true_gm
    
    oel <- exp( qnorm(1 - scenarios$true_exceedance_perc[i]/100, mean = log(true_gm) , sd = log(true_gsd) ) )
    
    loq <- signif(exp(qnorm(scenarios$proportion_censored[i], mean = log(true_gm) , sd = log(true_gsd) )),3)
    
    
    
    
    #calculus
    
    my_X <- vector("list", length = n_sim)
    
    for ( j in 1:n_sim ) { my_X[[j]] <- list( index = j,
                                              oel = oel[j]) }
    
    
    # calculations
    
    myresults <- vector("list", length = n_sim)
    
    for (j in 1:n_sim) { 
      
      
      
      myresults[[j]] <- ithpair.function( my_X[[j]]$index, 
                                            simulated_data_object = simulated_data_objects , 
                                            me_cv = me_cv , 
                                            n_iterations_gum = n_iterations_gum , 
                                            oel = my_X[[j]]$oel) 
      }
    
    ## error for j=129 bad sampling of GSDs (no quantiles used)
    
    true_gm[129]
    true_gsd[129]
    oel[129]
    loq[129]



