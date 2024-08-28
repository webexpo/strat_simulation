# Script for performing a simulation study with real GSD : 5000 iterations per scenario

#  errors encountered with JAGS on certain scenarios (see debug section)

#  solved with STAN

##### LIBRARIES ####

library(tolerance)
library(parallel)
library(rstan)
library(rjags)

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

for (i in 1:dim(scenarios)[1]) { #issue with 37,43,67

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
  
  
####  DEBUG ####                                                 

### debug attempts for the problematic scenarios : low exceedance and censoring, only patterns seems higher OELs

### no error at lower itration (<2000)

#### loop in loop approach to separate iteration did not work

### no error in the measurement model

### no influence if fixing the OEL instead of the P95

### no error in the frequentist approach



####  SEPARATE ANALYSIS #### 
    
###

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


###### generating the data ####

simulated_data_objects <- vector("list", length = 72)

for (i in 1:dim(scenarios)[1]) { 
  
  true_gsd <- sample( real_gsds[ real_gsds>=quantile(real_gsds,0.025) & real_gsds<=quantile(real_gsds,0.975)] , n_sim , replace = TRUE )
  
  true_gm <- exp( log(true_p95) - qnorm(0.95)*log(true_gsd) )

  loq <- signif(exp(qnorm(scenarios$proportion_censored[i], mean = log(true_gm) , sd = log(true_gsd) )),3)
  
  simulated_data_objects[[i]] <- data.simulation.seg.gen(  true_gm = true_gm, 
                                                     true_gsd = true_gsd,
                                                     sample_size = scenarios$sample_size[i]  ,
                                                     me_cv_inf = me_cv,
                                                     me_cv_sup = me_cv,
                                                     censor_level = loq )
  
  
}

###### ideal bayesian ####  
    
    ## error at 37
    
    simulation_results_ideal_b <- vector("list", length = 3)
    
    for (i in c(37,43,67)) { 
    
      
      oel <- exp( qnorm(1 - scenarios$true_exceedance_perc[i]/100, mean = log( simulated_data_objects[[i]]$true_gm) , sd = log(simulated_data_objects[[i]]$true_gsd) ) )
      
      
      simulation_results_ideal_b[[i]] <- parallel.function.ideal.b( simulated_data_object = simulated_data_objects[[i]], 
                                                                    n_sim = n_sim , 
                                                                    n_clusters = 20, oel = oel)
      
      
      print(i)
      
      print(simulation_results_ideal_b[[i]]$time)
      
    }
    

    
###### naive bayesian #### 
    
    # error
    
    simulation_results_naive_b <- vector("list", length = dim(scenarios)[1])
    
    for (i in c(37,43,67)) { 
      
    
      oel <- exp( qnorm(1 - scenarios$true_exceedance_perc[i]/100, mean = log( simulated_data_objects[[i]]$true_gm) , sd = log(simulated_data_objects[[i]]$true_gsd) ) )
      
      
      simulation_results_naive_b[[i]] <- parallel.function.naive.b( simulated_data_object = simulated_data_objects[[i]], 
                                                                    n_sim = n_sim , 
                                                                    n_clusters = 20, oel = oel)
      
      
      print(i)
      
      print(simulation_results_naive_b[[i]]$time)
      
    }    

    
###### me bayesian ####  
    
    #no error
    
    simulation_results_me_b <- vector("list", length = dim(scenarios)[1])
    
    for (i in c(37,43,67)) { 
      
    oel <- exp( qnorm(1 - scenarios$true_exceedance_perc[i]/100, mean = log( simulated_data_objects[[i]]$true_gm) , sd = log(simulated_data_objects[[i]]$true_gsd) ) )
      
      
     simulation_results_me_b[[i]] <- parallel.function.me.b( simulated_data_object = simulated_data_object, 
                                                                    n_sim = n_sim ,
                                                                    me_cv = me_cv,
                                                                    n_clusters = 4, oel = oel)
      
      
      print(i)
      
      print(simulation_results_me_b[[i]]$time)
      
    }    
###### ideal bayesian webexpo####  14h for 1 scenario, not realistic
    

    simulation_results_ideal_b_w <- vector("list", length = dim(scenarios)[1])
    
    for (i in c(37,43,67)) { 
      
    oel <- exp( qnorm(1 - scenarios$true_exceedance_perc[i]/100, mean = log( simulated_data_objects[[i]]$true_gm) , sd = log(simulated_data_objects[[i]]$true_gsd) ) )
      
      
    simulation_results_ideal_b_w[[i]] <- parallel.function.ideal.b.w( simulated_data_object = simulated_data_object, 
                                                              n_sim = n_sim ,
                                                              n_clusters = 4, oel = oel)
      
      
      print(i)
      
      print(simulation_results_ideal_b_w[[i]]$time)
      
    }     
    

    ###### bayesian stan ####
    
        ## prep
    
                # compiling models
                
                stan.models.list <- list()
    
                stan.models.list <- augment.stan.models.list(stan.models.list, stan.file=code_seg_informedvar)
                stan.models.list <- augment.stan.models.list(stan.models.list, stan.file=code_seg_informedvar_lognormal_mecv)
                stan.models.list <- augment.stan.models.list(stan.models.list, stan.file=code_seg_informedvar_lognormal_mecvknown)
                
                compiled.models.list(stan.models.list)
                
                # testing results on similar datasets
                
                test_data <- simulated_data_objects[[17]]$true[,1]
    
                test_naive_jags <- ithpair.function.ideal.b(index = 1 , 
                                                            simulated_data_object = simulated_data_objects[[17]] , 
                                                            oel = exp( qnorm(1 - scenarios$true_exceedance_perc[17]/100, mean = log( simulated_data_objects[[17]]$true_gm) , sd = log(simulated_data_objects[[17]]$true_gsd) ) )[1]
                )
                
                test_naive_stan <- ithpair.function.ideal.b.s(index = 1 , 
                                                            simulated_data_object = simulated_data_objects[[17]] , 
                                                            oel = exp( qnorm(1 - scenarios$true_exceedance_perc[17]/100, mean = log( simulated_data_objects[[17]]$true_gm) , sd = log(simulated_data_objects[[17]]$true_gsd) ) )[1],
                                                            models.list = stan.models.list
                )
                
                
                # testing parallel functions for JAGS and STAN for one scenario - NAIVE
                
                i <- 17
                
                oel <- exp( qnorm(1 - scenarios$true_exceedance_perc[i]/100, mean = log( simulated_data_objects[[i]]$true_gm) , sd = log(simulated_data_objects[[i]]$true_gsd) ) )
                
                simulation_results_ideal_b <- parallel.function.ideal.b( simulated_data_object = simulated_data_objects[[i]], 
                                                                              n_sim = 100 , 
                                                                              n_clusters = 8, oel = oel)
                
                
                print(simulation_results_ideal_b$time) # 7 secs
                

                simulation_results_ideal_b_s <- parallel.function.ideal.b.s( simulated_data_object = simulated_data_objects[[i]], 
                                                                         n_sim = 100 , 
                                                                         n_clusters = 8, oel = oel,
                                                                         models.list = stan.models.list)
                
                
                print(simulation_results_ideal_b_s$time) # 17 secs
                
                      
                # testing parallel functions for JAGS and STAN for one scenario - ME
                
                simulation_results_me_b <- parallel.function.me.b( simulated_data_object = simulated_data_objects[[i]], 
                                                                          me_cv = me_cv,
                                                                         n_sim = 100 , 
                                                                         n_clusters = 8, oel = oel)
                
                
                print(simulation_results_me_b$time) # 11 secondes
                
                
                simulation_results_me_b_s <- parallel.function.me.b.s( simulated_data_object = simulated_data_objects[[i]], 
                                                                         me_cv = me_cv,
                                                                             n_sim = 100 , 
                                                                             n_clusters = 8, oel = oel,
                                                                             models.list = stan.models.list)
                
                
                print(simulation_results_me_b_s$time) # 80 secondes
                
                # STAN FUNCTION ON THE PROBLEMATIC SCENARIOS : no error
        
          
        simulation_results_ideal_b_s <- vector("list", length = 3)
        
        
        for (i in c(37,43,67)) { 
          
          oel <- exp( qnorm(1 - scenarios$true_exceedance_perc[i]/100, mean = log( simulated_data_objects[[i]]$true_gm) , sd = log(simulated_data_objects[[i]]$true_gsd) ) )
          
          
          simulation_results_ideal_b_s[[i]] <- parallel.function.ideal.b.s( simulated_data_object = simulated_data_objects[[i]], 
                                                                        n_sim = n_sim , 
                                                                        n_clusters = 20, oel = oel)
          
          
          print(i)
          
          print(simulation_results_ideal_b_s[[i]]$time)
          
        }
        
    
   
 
    