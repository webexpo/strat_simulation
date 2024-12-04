# Script for the interpretation of the simulations - all scenarios - error CV = 0.125


##### LIBRARIES ####



##### SCRIPTS ####

source("full simulations/TD EXIL 2024 measurement error impact/scripts/simulation_interpretation_functions.R")

source("other/performance_metrics.R")


##### DATA ####

#init_path <- "F:/Dropbox/"
init_path <- "C:/jerome/Dropbox/"


run1_data <- readRDS(file = paste( init_path ,"GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/all_scenarios_GSD_run6a_data.RDS", sep=""))

run2_data <- readRDS(file = paste( init_path ,"GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/all_scenarios_GSD_run6b_data.RDS", sep=""))

run3_data <- readRDS(file = paste( init_path ,"GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/all_scenarios_GSD_run6c_data.RDS", sep=""))

run1 <- readRDS(paste( init_path ,"GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/all_scenarios_GSD_run6a2_sim.RDS", sep=""))

run2 <- readRDS(paste( init_path ,"GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/all_scenarios_GSD_run6b_sim.RDS", sep=""))

run3 <- readRDS(paste( init_path ,"GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/all_scenarios_GSD_run6c_sim.RDS", sep=""))

##### SCENARIOS ####

    scenarios <- expand.grid(
      true_exceedance_perc = c(0.1, 1, 3, 7, 10, 25),
      sample_size      = c(3L, 6L , 9L, 12L),
      proportion_censored        = c(0, 0.3, 0.6),
      stringsAsFactors = FALSE)
    
    scenarios$true_p95 <- 100


    ## fixed parameters
    
    n_sim <- 5000
    
    me_cv <- 0.25
    
    n_iterations_gum = 5000  
    
    true_p95 <- 100

    ## variable parameters
    
    for (i in 1:dim(scenarios)[1]) {
      
      run1_data[[i]]$oel <- exp( qnorm(1 - rep(scenarios$true_exceedance_perc[i],n_sim)/100, mean = log(run1_data[[i]]$true_gm) , sd = log(run1_data[[i]]$true_gsd) ) )
      
      run2_data[[i]]$oel <- exp( qnorm(1 - rep(scenarios$true_exceedance_perc[i],n_sim)/100, mean = log(run2_data[[i]]$true_gm) , sd = log(run2_data[[i]]$true_gsd) ) )
      
      run3_data[[i]]$oel <- exp( qnorm(1 - rep(scenarios$true_exceedance_perc[i],n_sim)/100, mean = log(run3_data[[i]]$true_gm) , sd = log(run3_data[[i]]$true_gsd) ) )
      
    }
    

    


#### RESULTS ####

###### computing time across all scenarios ######
    
    # computing time for each scenario in run 3
    
    computing_time1_h <- vector("numeric", length = dim(scenarios)[1])
    
    for ( i in 1:dim(scenarios)[1] ) {
      
      computing_time1_h[i] <- as.numeric(run1[[i]]$time, units = "hours") 
      
    }
    
    summary(computing_time1_h )
    sum( computing_time1_h )
    # computing time for each scenario in run 4
    
    computing_time2_h <- vector("numeric", length = dim(scenarios)[1])
    
    for ( i in 1:dim(scenarios)[1] ) {
      
      if (!is.null(run2[[i]]$time)) computing_time2_h[i] <- as.numeric(run2[[i]]$time, units = "hours", na.rm = FALSE) 
      
      else computing_time2_h[i] <- NA
    }
    
    summary(computing_time2_h )
    sum( computing_time2_h )
    
    
    # computing time for each scenario in run 5
    
    computing_time3_h <- vector("numeric", length = dim(scenarios)[1])
    
    for ( i in 1:dim(scenarios)[1] ) {
      
      computing_time3_h[i] <- as.numeric(run3[[i]]$time, units = "hours") 
      
    }
    
    summary(computing_time3_h )
    sum( computing_time3_h )
    
    
    
    
###### summary across all scenarios ######
    

#### for one scenario : mean and sd(%) of all performance metrics

### organisation : 

## scenario i

    # performance metric a

              #exposure metric 1 : dataframe with lines for the 7 methods and columns for reapeat 1, 2, 3, mean and CV%
    
              #exposure metric 2

    # performance metric b


## calculations

results_me_0.125_gsd <- vector("list", length = dim(scenarios)[1])

for ( i in 1:dim(scenarios)[1] ) {
  
  # summary of each run
  
  summary_repa <-one_run_summary_gsd( index=i , simulation_result=run1 , simulation_data = run1_data )
  
  summary_repb <-one_run_summary_gsd( index=i , simulation_result=run2 , simulation_data = run2_data)
  
  summary_repc <-one_run_summary_gsd( index=i , simulation_result=run3 , simulation_data = run3_data)  
  

  
  # summary across the 3 runs
  
  results_me_0.125_gsd[[i]] <- three_run_summary( summary_repa , summary_repb , summary_repc )
  
}

#### EXPORT ####

saveRDS( list( scenarios = scenarios,
               parameters = list( n_sim = n_sim, me_cv = me_cv, n_iterations_gum = n_iterations_gum, true_p95 = true_p95 ),
               results = results_me_0.125_gsd ), file = paste( init_path ,"GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/aggregated results/cv0.125_realgsd.RDS", sep=""))    
    
    