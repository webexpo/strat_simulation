#' functions that estimates the compliance probability based on the expostats bayesian model  or a frequentist model
#' 
#' compliance is defined as : the X% upper confidence limit on the 95th percentile is below the OEL
#' 
#' inputs include sample size, true GSD, true 95th percentile relative to OEL, and level of confidence 
#' 
#' These functions assume there is no non detect

##### Libraries #####  

library(rjags)

##### External scripts #####

source("data generation/load.webexpo.data_gen.functions.R")
source("parameter estimation/bayesian/load.webexpo.SEG.functions.R")
source("parameter estimation/frequentist/percentile.R")

##### Fonctions #####

#' Compliance function using the Expostats Bayesian model for a single sample
#' 
#' 
#' @param raw_data vector of values to be analysed ( character format, see default )
#' @param confidence_level_perc Level of confidence in % for the upper confidence limit 
#' 
#' @return decision : either "compliant" or "not compliant"

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

#' Compliance function using a frequentist equation for a single sample
#' 
#' 
#' @param raw_data vector of values to be analysed ( character format, see default )
#' @param confidence_level_perc Level of confidence in % for the upper confidence limit 
#' 
#' @return decision : either "compliant" or "not compliant"


decision_brent_frequentist <- function( raw_data = c("0.39", "0.59", "0.93", "0.166", "0.86", "0.10"),
                                     confidence_level_perc = 70 ) {
  
  
  ### frequentist equation
  
  estimated_p95_ucl <-  fun.perc( as.numeric( raw_data ) , alpha= 1-confidence_level_perc/100,perc=0.95)$uc 
    
  #### output metrics - risk decision
  
  if ( estimated_p95_ucl < 1) { decision <- "compliant" } else { decision <- "not compliant" }
  
  ##### final return
  
  return( decision )
  
}


#' Compliance probability function ( proportion of "compliant" across many repetitions )using the Expostats Bayesian model
#' 
#' 
#' @param sample_size sample size for the  strategy
#' @param ratio_p95overoel True 95th percentile relative to the OEL 
#' @param gsd_value True geometric standard deviation
#' @param n_sim Number of repetitions of the simulation 
#' @param confidence_level_perc Level of confidence in % for the upper confidence limit  
#' 
#' @return proportion of the repetitions for which the decision was "compliant" in %


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



#' Compliance probability function ( proportion of "compliant" across many repetitions )using a frequientist model
#' 
#' 
#' @param sample_size sample size for the  strategy
#' @param ratio_p95overoel True 95th percentile relative to the OEL 
#' @param gsd_value True geometric standard deviation
#' @param n_sim Number of repetitions of the simulation 
#' @param confidence_level_perc Level of confidence in % for the upper confidence limit  
#' 
#' @return proportion of the repetitions for which the decision was "compliant" in %


simulation_brent_frequentist <- function( sample_size           = 10,
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
                              FUN = decision_brent_frequentist ,
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

confidence_level_perc <- 80

##### Example ##### ( ~ 5 mins for n_sim = 1000 for the Bayesian model, ideally, n_sim should be closer to 5000 even more for good reproducibility ). the frequentist function is instantaneous up to n_sim=100 000

simulation_brent_bayesian( sample_size           = sample_size,
                           ratio_p95overoel      = ratio_p95overoel,
                           gsd_value             = gsd_value,
                           n_sim                 = n_sim,
                           confidence_level_perc = confidence_level_perc )

simulation_brent_frequentist( sample_size           = sample_size,
                           ratio_p95overoel      = ratio_p95overoel,
                           gsd_value             = gsd_value,
                           n_sim                 = n_sim,
                           confidence_level_perc = confidence_level_perc )

