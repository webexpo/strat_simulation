#' function that estimates the compliance probability based on the expostats bayesian model  / a frequentist model
#' 
#' compliance is defined as : the X% upper confidence limit on the 95th percentile is below the OEL
#' 
#' inputs include sample size, true GSD, true 95th percentile relative to OEL, level of confidence 
#' 
#' These functions assume no non detects

##### Libraries #####  

library(rjags)

##### External scripts #####

source("data generation/load.webexpo.data_gen.functions.R")
source("parameter estimation/bayesian/load.webexpo.SEG.functions.R")

##### Fonctions #####

decision_brent_bayesian <- function( raw_data = c("0.39", "0.59", "0.93", "0.166", "0.86", "0.10"),
                                     confidence_level_perc = 70 ) {
  
  
  
  
  ### bayesian analysis
  
  hmc <- Webexpo.seg.globalbayesian.jags( data.sample  =  raw_data , n.iter = 25000 , oel = 1)  
  
  estimated_p95_ucl <-  perc( mu.chain    = hmc$mu.chain , 
                              sigma.chain = hmc$sigma.chain , 
                              target_perc = 95,
                              probacred   = confidence_level_perc*2 - 100  )$ucl
  
  
  #### output metrics - risk decision
  
  if ( estimated_p95_ucl < 1) { decision <- "compliant" } else { decision <- "not compliant" }
  
  ##### final return
  
  return( decision )
  
}


simulation_brent_bayesian <- function( sample_size           = 10,
                                       ratio_p95overoel      = 1,
                                       gsd_value             = 2.5,
                                       n_sim                 = 500,
                                       confidence_level_perc = 70 ) {
  
  ## intermediate calculations
  
  gm_value <- exp( ( log( ratio_p95overoel) - qnorm(0.95)*log( gsd_value ) ) )
  
  ## data generation
  
  simulated_values <-  webexpo.seg.gener.LN( n            = sample_size * n_sim,
                                             no.censoring = TRUE,
                                             gm           = gm_value,
                                             gsd          = gsd_value)
  
  simulated_data_matrix <- matrix( data = simulated_values , nrow = sample_size )
  
  ## Looping across the n.sim iterations
  
  simulation_result <- apply( X = simulated_data_matrix ,
                              MARGIN = 2,
                              FUN = decision_brent_bayesian ,
                              confidence_level_perc = confidence_level_perc ,
                              simplify = TRUE)
  
  # result
   
  decision_perc_compliant <- 100*length(simulation_result[simulation_result=="compliant"])/length(simulation_result)
  
  
  return(decision_perc_compliant)
  
}

##### inputs #####

sample_size <- 10

ratio_p95overoel <- 1

gsd_value <- 2.5

n_sim <- 1000

confidence_level_perc <- 95

##### Example ##### ( ~ 5 mins for n_sim = 10000 )

simulation_brent_bayesian( sample_size           = sample_size,
                           ratio_p95overoel      = ratio_p95overoel,
                           gsd_value             = gsd_value,
                           n_sim                 = n_sim,
                           confidence_level_perc = confidence_level_perc )