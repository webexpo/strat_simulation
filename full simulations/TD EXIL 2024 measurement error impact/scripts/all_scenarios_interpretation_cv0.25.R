# Script for the interpretation of the simulations - all scenarios - error CV = 0.25


##### LIBRARIES ####



##### Scripts ####

source("full simulations/TD EXIL 2024 measurement error impact/scripts/simulation_interpretation_functions.R")

source("other/performance_metrics.R")


##### parameters of the simulation ####
    
    scenarios <- expand.grid(
      true_gsd        = c(1.5, 2.5, 3.5),
      true_exceedance_perc = c(0.1, 1, 3, 7, 10, 25),
      sample_size      = c(3L, 6L , 9L, 12L),
      stringsAsFactors = FALSE)
    
    ##parameters depending on scenario##
    
    scenarios$true_p95 <- 100
    
    scenarios$true_gm <- exp( log(scenarios$true_p95) - qnorm(0.95)*log(scenarios$true_gsd) )
    
    scenarios$oel <- exp( qnorm(1 - scenarios$true_exceedance_perc/100, mean = log(scenarios$true_gm) , sd = log(scenarios$true_gsd) ) )
    
    scenarios$proportion_censored <- 0 # need to adress the issue of "all points censored", or ros log(GSD) == 0
    
    scenarios$loq <- exp(qnorm(scenarios$proportion_censored, mean = log(scenarios$true_gm) , sd = log(scenarios$true_gsd) ))

    
    ## fixed parameters
    
    n_sim <- 5000
    
    me_cv <- 0.25
    
    n_iterations_gum = 5000  

##### DATA ####

init_path <- "F:/Dropbox/"
init_path <- "C:/jerome/Dropbox/"


load(file = paste( init_path ,"GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/all_scenarios_run3a_sim.RDS", sep=""))
run3a <- simulation_results

load(file = paste( init_path ,"GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/all_scenarios_run3b_sim.RDS", sep=""))
run3b <- simulation_results

run3c <- readRDS(paste( init_path ,"GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/all_scenarios_run3c_sim.RDS", sep=""))

##### SCENARIOS ####

scenarios <- expand.grid(
  true_gsd        = c(1.5, 2.5, 3.5),
  true_exceedance_perc = c(0.1, 1, 3, 7, 10, 25),
  sample_size      = c(3L, 6L , 9L, 12L),
  stringsAsFactors = FALSE)

scenarios$true_p95 <- 100

scenarios$true_gm <- exp( log(scenarios$true_p95) - qnorm(0.95)*log(scenarios$true_gsd) )

scenarios$oel <- exp( qnorm(1 - scenarios$true_exceedance_perc/100, mean = log(scenarios$true_gm) , sd = log(scenarios$true_gsd) ) )


##### FUNCTION ####


#' function which creates the interpretation of the simulations for a single scenario and one run 
#'
#' @param index index of the scenario of interest
#' @param simulation_result object of the results of a full run of the simulation 
#'
#' @return lits of results with all performance metrics calculated across the iterationz
#'


one_run_summary <- function( index , simulation_result ) {
  
  
  run_inter <- list( rmse = rmse.result( results_one_scenario = simulation_result[[index]] , 
                                         true_gm = scenarios$true_gm[index] ,
                                         true_gsd = scenarios$true_gsd[index] ,
                                         true_p95 = scenarios$true_p95[index] ,
                                         true_exceedance_perc = scenarios$true_exceedance_perc[index] ),
                     
                     precision = precision.result( results_one_scenario = simulation_result[[index]] , 
                                                   true_gm = scenarios$true_gm[index] , 
                                                   true_gsd = scenarios$true_gsd[index] ,
                                                   true_p95 = scenarios$true_p95[index] ,
                                                   true_exceedance_perc = scenarios$true_exceedance_perc[index] ),
                     
                     bias = bias.result( results_one_scenario = simulation_result[[index]] , 
                                         true_gm = scenarios$true_gm[index] , 
                                         true_gsd = scenarios$true_gsd[index] ,
                                         true_p95 = scenarios$true_p95[index] ,
                                         true_exceedance_perc = scenarios$true_exceedance_perc[index] ),
                     
                     median_error = median.error.result( results_one_scenario = simulation_result[[index]] , 
                                                         true_gm = scenarios$true_gm[index] , 
                                                         true_gsd = scenarios$true_gsd[index] ,
                                                         true_p95 = scenarios$true_p95[index] , 
                                                         true_exceedance_perc = scenarios$true_exceedance_perc[index] ),
                     
                     rmsle = rmsle.result( results_one_scenario = simulation_result[[index]] , 
                                           true_gm = scenarios$true_gm[index] , 
                                           true_gsd = scenarios$true_gsd[index] , 
                                           true_p95 = scenarios$true_p95[index] , 
                                           true_exceedance_perc = scenarios$true_exceedance_perc[index] ),
                     
                     mad = mad.result( results_one_scenario = simulation_result[[index]] , 
                                       true_gm = scenarios$true_gm[index] , 
                                       true_gsd = scenarios$true_gsd[index] , 
                                       true_p95 = scenarios$true_p95[index] , 
                                       true_exceedance_perc = scenarios$true_exceedance_perc[index] ),
                     
                     coverage = coverage.result( results_one_scenario = simulation_result[[index]] , 
                                                 true_p95 = scenarios$true_p95[i] , 
                                                 true_exceedance_perc = scenarios$true_exceedance_perc[i] ),
                     
                     perc_mistake = perc.mistake.result( results_one_scenario = simulation_result[[index]] , 
                                                         true_p95 = scenarios$true_p95[i] , 
                                                         true_exceedance_perc = scenarios$true_exceedance_perc[i] ,
                                                         oel = scenarios$oel[i] )
  )
  
return(run_inter)
  
  
  
}
  

#### Time calculation ####


#### calculating the total computing time for the 3 runs

run3a_time <- sum( sapply( run3a , function(x) as.numeric( x$time , units = "hours") ) )

run3b_time <- sum( sapply( run3b , function(x) as.numeric( x$time , units = "hours") ) )

run3c_time <- sum( sapply( run3c , function(x) as.numeric( x$time , units = "hours") ) )



#### for one scenario : mean and sd(%) of all performance metrics

### organisation : 

## scenario i

    # performance metric a

              #exposure metric 1 : dataframe with lines for the 7 methods and columns for reapeat 1, 2, 3, mean and CV%
    
              #exposure metric 2

    # performance metric b


## calculations

results_me_0.25 <- vector("list", length = dim(scenarios)[1])


i <- 23


# summary of the runs

summary_repa <-one_run_summary( index=i , simulation_result=run3a )

summary_repb <-one_run_summary( index=i , simulation_result=run3b )

summary_repc <-one_run_summary( index=i , simulation_result=run3c )  

# rmse

results_me_0.25[[i]]$rmse <- vector("list", length = 5)


    #gm
    
    results_me_0.25[[i]]$rmse$gm <- data.frame( method = summary_repa$rmse$method ,
                                                rep1 = summary_repa$rmse$gm , 
                                                rep2 = summary_repb$rmse$gm , 
                                                rep3 = summary_repc$rmse$gm  )
    
    
            # ading the mean
            
            results_me_0.25[[i]]$rmse$gm$mean <- rowMeans( results_me_0.25[[i]]$rmse$gm[ , 2:4 ] )
        
            # adding the sd
            
            results_me_0.25[[i]]$rmse$gm$sd <- apply( results_me_0.25[[i]]$rmse$gm[ , 2:4 ] , 1 , sd )*100 / results_me_0.25[[i]]$rmse$gm$mean
    
    #gsd
            
    results_me_0.25[[i]]$rmse$gsd <- data.frame( method = summary_repa$rmse$method ,
                                                rep1 = summary_repa$rmse$gsd , 
                                                rep2 = summary_repb$rmse$gsd , 
                                                rep3 = summary_repc$rmse$gsd  )
    
    
            # ading the mean
            
            results_me_0.25[[i]]$rmse$gsd$mean <- rowMeans( results_me_0.25[[i]]$rmse$gsd[ , 2:4 ] )
        
            # adding the sd
            
            results_me_0.25[[i]]$rmse$gsd$sd <- apply( results_me_0.25[[i]]$rmse$gsd[ , 2:4 ] , 1 , sd )*100 / results_me_0.25[[i]]$rmse$gsd$mean        
