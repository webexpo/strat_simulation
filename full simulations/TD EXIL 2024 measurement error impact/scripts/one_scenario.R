# Script for peforming a simulation study with one scenario


##### LIBRARIES ####

library(tolerance)
library(parallel)

##### DATA ####

real_gsds <- readRDS( "created data/real_gsd_values.RDS")


##### SCRIPTS ####

source("full simulations/TD EXIL 2024 measurement error impact/scripts/sample_analysis.R")

source("full simulations/TD EXIL 2024 measurement error impact/scripts/simulation_interpretation.R")

source("other/support_functions.R")

source("other/performance_metrics.R")

source("parameter estimation/bayesian/load.webexpo.SEG.functions.R")

source("parameter estimation/frequentist/percentile.R")

source("parameter estimation/frequentist/exceedance_fraction.R")

source("data generation/SEG data simulation.R")

##### FUNCTIONS ####

#' function which analyses the ith pair of samples ( clean and dirtied ) in the simulation 
#'
#' @param index index of the simulation 
#' @param simulated_data_object list containing the simulated data for a single scenario 
#' @param me.cv value for the measurement error as standard error CV, not in % and not expanded uncertainty CV
#' @param oel occupational exposure limit
#' @param n_iterations number of iteration for the GUM approach
#' @param sim_quantile quantile selection for the summaries of metrics across iterations in addition to mean and median
#'
#' @return matrix of results with one column per approach
#'


ithpair.function <- function( index , simulated_data_object , me_cv , n_iterations_gum = 10000 , sim_quantile = 0.95, oel ) {
  
  # data preparation
  
  true_data <- simulated_data_object$true[,index]
  
  observed_data <- simulated_data_object$observed[,index]
  
  ## analysis
  
  ideal_b = expostats.naive( true_data , oel)
  ideal_f = frequentist.naive( true_data , oel)
  naive_b = expostats.naive( observed_data , oel)
  naive_f = frequentist.naive( observed_data , oel)
  me_b = expostats.me( observed_data , oel , me_cv)
  me_f = frequentist.me( observed_data , oel , me_cv , n_iterations_gum , sim_quantile ) # is a list of 3 vectors (mean, median,quantile)
  
  ## results as a single matrix
  
  one_column <- c("gm_est","gsd_est","p95_est","p95_70ucl","p95_95ucl","F_est","F_70ucl","F_95ucl")
  
  results <- matrix( c( ideal_b , ideal_f , naive_b , naive_f , me_b , me_f$mean , me_f$median , me_f$quantile ) , nrow = length(one_column) )

  colnames(results) <- c("ideal_b","ideal_f","naive_b","naive_f","me_b","me_f_mean","me_f_median","me_f_quantile")
  rownames(results) <- one_column
  
  return(results)
  
}


#' function which uses parallel computing to perform the simulation for one data generation obkect (one scenario)
#'
#' @param index index of the simulation 
#' @param simulated_data_object list containing the simulated data for a single scenario 
#' @param me.cv value for the measurement error as standard error CV, not in % and not expanded uncertainty CV
#' @param oel occupational exposure limit
#' @param n_iterations number of iteration for the GUM approach
#' @param sim_quantile quantile selection for the summaries of metrics across iterations in addition to mean and median
#'
#' @return matrix of results with one column per approach
#'


parallel.function <- function( simulated_data_object , me_cv , n_iterations_gum = 10000 , sim_quantile = 0.95 , n_sim , n_clusters = 10, oel ) {
  
  # compteur de temps initialisé
  start_time <- Sys.time()
  
  
  #procédure parallele
  
  # creating the clusters ( use detectCores() to count the cores available, choose this number minus 2 or 4 below) 
  
  cl <- makeCluster(n_clusters)
  
  # libraries and scripts to be used in each cluster
  
  clusterEvalQ(cl, library(rjags))
  clusterEvalQ(cl, library(tolerance))
  
  # sending objects to clusters
  clusterExport( cl , "jags.model.informedvar" , envir=environment())
  clusterExport( cl , "webexpo.seg.datapreparation" , envir=environment())
  clusterExport( cl , "Webexpo.seg.globalbayesian.jags" , envir=environment())
  clusterExport( cl , "fun.jags.informedvar" , envir=environment())
  clusterExport( cl , "expostats.naive" , envir=environment())
  clusterExport( cl , "expostats.me" , envir=environment())
  clusterExport( cl , "frequentist.naive" , envir=environment())
  clusterExport( cl , "frequentist.me" , envir=environment())
  clusterExport( cl , "ithpair.function" , envir=environment())
  clusterExport( cl , "fun.NdExpo.lognorm" , envir=environment())
  clusterExport( cl , "fun.perc.en689" , envir=environment())
  clusterExport( cl , "fun.frac.dep" , envir=environment())
  
  
  clusterExport( cl , "sample_size" , envir=environment())
  clusterExport( cl , "oel" , envir=environment())
  clusterExport( cl , "me_cv" , envir=environment())
  clusterExport( cl , "n_iterations_gum" , envir=environment())
  clusterExport( cl , "sim_quantile" , envir=environment())
  clusterExport( cl , "test_data_sim" , envir=environment())
  # calculations
  
  simulation_result_parallel <- parLapply(cl, X = as.list(1:n_sim) , function(x){ ithpair.function(x, 
                                                                                                   simulated_data_object = simulated_data_object , 
                                                                                                   me_cv = me_cv , 
                                                                                                   n_iterations_gum = n_iterations_gum , 
                                                                                                   sim_quantile = sim_quantile,
                                                                                                   oel = oel) } ) 
  
  # recommendation from the net: close the clusters
  stopCluster(cl)
  
  # making an array of the results
  
  simulation_result_parallel_array <- array( data = unlist(simulation_result_parallel) , dim = c( 8 , 8 , n_sim ) )
  
  # estimation of computing time ( 9 min on my computer for 5000 iterations)
  end_time <- Sys.time()
  mytime <- end_time - start_time 
  
  # results
  
  results <- list( array = simulation_result_parallel_array , time = mytime )
  
  return(results)
  
}
  


##### PARAMETERS ####

true_gsd <- 2.5

true_exceedance_perc <- 2.5  ## in percentage

true_p95 <- 100

true_gm <- exp( log(true_p95) - qnorm(0.95)*log(true_gsd) )

oel <- exp( qnorm(1 - true_exceedance_perc/100, mean = log(true_gm) , sd = log(true_gsd) ) )

proportion_censored <- 0 # need to adress the issue of "all points censored", or ros log(GSD) == 0

loq <- exp(qnorm(proportion_censored, mean = log(true_gm) , sd = log(true_gsd) ))

sample_size <- 6

n_sim <- 5000

me_cv <- 0.25

expostats_sample <- c("28.9","19.4","<5.5","149.9","26.42","56.1")

n_iterations_gum = 5000

sim_quantile = 0.95

##### ANALYSES ####

###### testing sample generation ####    

test_data_sim <- data.simulation.seg( true_gsd = true_gsd, 
                             true_gm = true_gm, 
                             n_sim = n_sim, 
                             sample_size = sample_size ,
                             me_cv_inf = me_cv,
                             me_cv_sup = me_cv,
                             censor_level = loq )


###### testing sample analysis functions ####  

test_expostats_naive <- expostats.naive( mysample = expostats_sample , oel)

test_expostats_me <- expostats.me( mysample = expostats_sample , oel , me_cv)

test_frequentist_naive <- frequentist.naive( mysample = expostats_sample , oel  )

test_frequentist_me <- frequentist.me( mysample = expostats_sample , oel , me_cv , n_iterations = n_iterations_gum , sim_quantile )

View(data.frame( naive_b = test_expostats_naive , naive_f=test_frequentist_naive , me_b=test_expostats_me , me_f=test_frequentist_me$mean))


###### testing ith.pair function ####  

test_ith.pair <- ithpair.function( index = 1 , simulated_data_object = test_data_sim , me_cv = me_cv , 
                                    n_iterations_gum = n_iterations_gum , sim_quantile = sim_quantile , oel = oel)

test_ith.pair

###### testing parallel function ####  

test_parallel <- parallel.function( simulated_data_object = test_data_sim , me_cv = me_cv , 
                                         n_iterations_gum = n_iterations_gum , sim_quantile = sim_quantile , n_sim = n_sim , 
                                         n_clusters = 10, oel = oel)

        
        
test_parallel$time

#data interpretation : mean over the simulation
        
apply(test_parallel$array, c(1,2), mean)
        
        
###### testing performance functions ####      
        
test_rmse <- rmse.result( results_one_scenario = test_parallel , 
                          true_gm = true_gm , true_gsd = true_gsd , true_p95 = true_p95 , 
                          true_exceedance_perc = true_exceedance_perc )
  

test_precision <- precision.result( results_one_scenario = test_parallel , 
                          true_gm = true_gm , true_gsd = true_gsd , true_p95 = true_p95 , 
                          true_exceedance_perc = true_exceedance_perc )

test_bias <- bias.result( results_one_scenario = test_parallel , 
                                    true_gm = true_gm , true_gsd = true_gsd , true_p95 = true_p95 , 
                                    true_exceedance_perc = true_exceedance_perc )

test_median_error <- median.error.result( results_one_scenario = test_parallel , 
                          true_gm = true_gm , true_gsd = true_gsd , true_p95 = true_p95 , 
                          true_exceedance_perc = true_exceedance_perc )

test_rmsle <- rmsle.result( results_one_scenario = test_parallel , 
                          true_gm = true_gm , true_gsd = true_gsd , true_p95 = true_p95 , 
                          true_exceedance_perc = true_exceedance_perc )

test_mad <- mad.result( results_one_scenario = test_parallel , 
                            true_gm = true_gm , true_gsd = true_gsd , true_p95 = true_p95 , 
                            true_exceedance_perc = true_exceedance_perc )

test_coverage <- coverage.result( results_one_scenario = test_parallel ,
                                  true_p95 = true_p95 , 
                                  true_exceedance_perc = true_exceedance_perc )


test_perc_mistake <- perc.mistake.result( results_one_scenario = test_parallel ,
                                  true_p95 = true_p95 , 
                                  true_exceedance_perc = true_exceedance_perc,
                                  oel=oel)

###### testing repeatability ####

test_repeatability <- vector("list", length = 5)

for (i in 1:5) {
  
  test_repeatability[[i]] <- parallel.function( simulated_data_object = test_data_sim , me_cv = me_cv , 
                                         n_iterations_gum = n_iterations_gum , sim_quantile = sim_quantile , n_sim = n_sim , 
                                         n_clusters = 10, oel = oel)
  
}

test_repeatability[[1]]$time
test_repeatability[[5]]$time

saveRDS(test_repeatability, "full simulations/TD EXIL 2024 measurement error impact/data/test_repeatability1.RDS")

