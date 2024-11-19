# Script for the interpretation of the simulations - all scenarios - error CV = 0.25


##### LIBRARIES ####



##### Scripts ####

source("full simulations/TD EXIL 2024 measurement error impact/scripts/simulation_interpretation_functions.R")

source("other/performance_metrics.R")




##### DATA ####

#init_path <- "F:/Dropbox/"
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

## fixed parameters

n_sim <- 5000

me_cv <- 0.25

n_iterations_gum = 5000  

true_p95 <- 100

#### RESULTS ####

###### computing time across all scenarios ######

# computing time for each scenario in run a

computing_timea_h <- vector("numeric", length = dim(scenarios)[1])

for ( i in 1:dim(scenarios)[1] ) {
  
  computing_timea_h[i] <- as.numeric(run3a[[i]]$time, units = "hours") 
  
}

summary(computing_timea_h )
sum( computing_timea_h )
# computing time for each scenario in run b

computing_timeb_h <- vector("numeric", length = dim(scenarios)[1])

for ( i in 1:dim(scenarios)[1] ) {
  
  if (!is.null(run3b[[i]]$time)) computing_timeb_h[i] <- as.numeric(run3b[[i]]$time, units = "hours", na.rm = FALSE) 
  
  else computing_timeb_h[i] <- NA
}

summary(computing_timeb_h )
sum( computing_timeb_h )


# computing time for each scenario in run 5

computing_timec_h <- vector("numeric", length = dim(scenarios)[1])

for ( i in 1:dim(scenarios)[1] ) {
  
  computing_timec_h[i] <- as.numeric(run3c[[i]]$time, units = "hours") 
  
}

summary(computing_timec_h )
sum( computing_timec_h )


#total : 

sum( computing_timea_h )+sum( computing_timeb_h )+sum( computing_timec_h ) #290 hours
#### for one scenario : mean and sd(%) of all performance metrics

### organisation : 

## scenario i

    # performance metric a

              #exposure metric 1 : dataframe with lines for the 7 methods and columns for reapeat 1, 2, 3, mean and CV%
    
              #exposure metric 2

    # performance metric b


## calculations

results_me_0.25 <- vector("list", length = dim(scenarios)[1])

for ( i in 1:dim(scenarios)[1] ) {
  
  # summary of each run
  
  summary_repa <-one_run_summary( index=i , simulation_result=run3a )
  
  summary_repb <-one_run_summary( index=i , simulation_result=run3b )
  
  summary_repc <-one_run_summary( index=i , simulation_result=run3c )  
  
  # summary across the 3 runs
  
  results_me_0.25[[i]] <- three_run_summary( summary_repa , summary_repb , summary_repc )
  
}


#### EXPORT ####

saveRDS( list( scenarios = scenarios,
               parameters = list( n_sim = n_sim, me_cv = me_cv, n_iterations_gum = n_iterations_gum, true_p95 = true_p95 ),
               results = results_me_0.25 ), file = paste( init_path ,"GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/aggregated results/cv0.25.RDS", sep=""))    

