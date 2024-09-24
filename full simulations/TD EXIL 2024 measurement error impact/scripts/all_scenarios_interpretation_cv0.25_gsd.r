# Script for the interpretation of the simulations - all scenarios - error CV = 0.25


##### LIBRARIES ####



##### SCRIPTS ####

source("full simulations/TD EXIL 2024 measurement error impact/scripts/simulation_interpretation_functions.R")

source("other/performance_metrics.R")


##### DATA ####

#init_path <- "F:/Dropbox/"
init_path <- "C:/jerome/Dropbox/"


run3 <- readRDS(paste( init_path ,"GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/all_scenarios_GSD_run3f_sim.RDS", sep=""))

run4 <- readRDS(paste( init_path ,"GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/all_scenarios_GSD_run4e_sim.RDS", sep=""))

run5 <- readRDS(paste( init_path ,"GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/all_scenarios_GSD_run5c_sim.RDS", sep=""))

##### SCENARIOS ####

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




#### RESULTS ####

###### computing time across all scenarios ######
    
    # computing time for each scenario in run 3
    
    computing_time3_h <- vector("numeric", length = dim(scenarios)[1])
    
    for ( i in 1:dim(scenarios)[1] ) {
      
      computing_time3_h[i] <- as.numeric(run3[[i]]$time, units = "hours") 
      
    }
    
    summary(computing_time3_h )
    sum( computing_time3_h )
    # computing time for each scenario in run 4
    
    computing_time4_h <- vector("numeric", length = dim(scenarios)[1])
    
    for ( i in 1:dim(scenarios)[1] ) {
      
      if (!is.null(run4[[i]]$time)) computing_time4_h[i] <- as.numeric(run4[[i]]$time, units = "hours", na.rm = FALSE) 
      
      else computing_time4_h[i] <- NA
    }
    
    summary(computing_time4_h )
    sum( computing_time4_h )
    
    
    # computing time for each scenario in run 5
    
    computing_time5_h <- vector("numeric", length = dim(scenarios)[1])
    
    for ( i in 1:dim(scenarios)[1] ) {
      
      computing_time5_h[i] <- as.numeric(run5[[i]]$time, units = "hours") 
      
    }
    
    summary(computing_time5_h )
    sum( computing_time5_h )
    
    
    
    
###### summary across all scenarios ######
    

#### for one scenario : mean and sd(%) of all performance metrics

### organisation : 

## scenario i

    # performance metric a

              #exposure metric 1 : dataframe with lines for the 7 methods and columns for reapeat 1, 2, 3, mean and CV%
    
              #exposure metric 2

    # performance metric b


## calculations

results_me_0.125 <- vector("list", length = dim(scenarios)[1])

for ( i in 1:dim(scenarios)[1] ) {
  
  # summary of each run
  
  summary_repa <-one_run_summary( index=i , simulation_result=run4a )
  
  summary_repb <-one_run_summary( index=i , simulation_result=run4b )
  
  summary_repc <-one_run_summary( index=i , simulation_result=run4c )  
  
  # summary across the 3 runs
  
  results_me_0.125[[i]] <- three_run_summary( summary_repa , summary_repb , summary_repc )
  
}


