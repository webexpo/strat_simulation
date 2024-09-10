# Script for the interpretation of the simulations - all scenarios - error CV = 0.25


##### LIBRARIES ####



##### Scripts ####

source("full simulations/TD EXIL 2024 measurement error impact/scripts/simulation_interpretation_functions.R")

source("other/performance_metrics.R")


##### DATA ####

#init_path <- "F:/Dropbox/"
init_path <- "C:/jerome/Dropbox/"


run4a <- readRDS(paste( init_path ,"GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/all_scenarios_run4a_sim.RDS", sep=""))

run4b <- readRDS(paste( init_path ,"GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/all_scenarios_run4b3_sim.RDS", sep=""))

run4c <- readRDS(paste( init_path ,"GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/all_scenarios_run4c_sim.RDS", sep=""))

##### SCENARIOS ####

scenarios <- expand.grid(
  true_gsd        = c(1.5, 2.5, 3.5),
  true_exceedance_perc = c(0.1, 1, 3, 7, 10, 25),
  sample_size      = c(3L, 6L , 9L, 12L),
  stringsAsFactors = FALSE)

scenarios$true_p95 <- 100

scenarios$true_gm <- exp( log(scenarios$true_p95) - qnorm(0.95)*log(scenarios$true_gsd) )

scenarios$oel <- exp( qnorm(1 - scenarios$true_exceedance_perc/100, mean = log(scenarios$true_gm) , sd = log(scenarios$true_gsd) ) )



#### RESULTS ####


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


