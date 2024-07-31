# Functions for performing a simulation study with one scenario



##### FUNCTIONS ####

#' function which analyses the ith pair of samples ( clean and dirtied ) in the simulation 
#'
#' @param index index of the simulation 
#' @param simulated_data_object list containing the simulated data for a single scenario 
#' @param me.cv value for the measurement error as standard error CV, not in % and not expanded uncertainty CV
#' @param oel occupational exposure limit
#' @param n_iterations number of iteration for the GUM approach
#'
#' @return matrix of results with one column per approach
#'


ithpair.function <- function( index , simulated_data_object , me_cv , n_iterations_gum = 10000 , oel ) {
  
  # data preparation
  
  true_data <- simulated_data_object$true[,index]
  
  observed_data <- simulated_data_object$observed[,index]
  
  ## analysis
  
  ideal_b = expostats.naive( true_data , oel)
  ideal_f = frequentist.naive( true_data , oel)
  naive_b = expostats.naive( observed_data , oel)
  naive_f = frequentist.naive( observed_data , oel)
  me_b = expostats.me( observed_data , oel , me_cv)
  me_f = frequentist.me( observed_data , oel , me_cv , n_iterations_gum  ) # is a list of 6 vectors (mean, median,quantile)
  
  ## results as a single matrix
  
  one_column <- c("gm_est","gsd_est","p95_est","p95_70ucl","p95_95ucl","F_est","F_70ucl","F_95ucl")
  
  results <- matrix( c( ideal_b , ideal_f , naive_b , naive_f , me_b , me_f$mean , me_f$median , me_f$q2.5 , me_f$q5, me_f$q95 , me_f$q97.5) , nrow = length(one_column) )

  colnames(results) <- c("ideal_b","ideal_f","naive_b","naive_f","me_b","me_f_mean","me_f_median","me_f_q2.5","me_f_q5","me_f_q95","me_f_q97.5")
  rownames(results) <- one_column
  
  return(results)
  
}


#' function which uses parallel computing to perform the simulation for one data generation object (one scenario)
#'
#' @param simulated_data_object list containing the simulated data for a single scenario 
#' @param me.cv value for the measurement error as standard error CV, not in % and not expanded uncertainty CV
#' @param oel occupational exposure limit : vector across iterations
#' @param n_iterations number of iteration for the GUM approach
#' @param sim_quantile quantile selection for the summaries of metrics across iterations in addition to mean and median
#'
#' @return matrix of results with one column per approach
#'


parallel.function <- function( simulated_data_object , me_cv , n_iterations_gum = 10000 , n_sim , n_clusters = 10, oel ) {
  
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
  
  
  clusterExport( cl , "oel" , envir=environment())
  clusterExport( cl , "me_cv" , envir=environment())
  clusterExport( cl , "n_iterations_gum" , envir=environment())
  clusterExport( cl , "simulated_data_object" , envir=environment())
  
  # list for the LApply function
  
  my_X <- vector("list", length = n_sim)
    
  for ( i in 1:n_sim ) { my_X[[i]] <- list( index = i,
                                            oel = oel[i]) }
  

  # calculations
  
  simulation_result_parallel <- parLapply(cl, X = my_X , function(x){ ithpair.function(x$index, 
                                                                                       simulated_data_object = simulated_data_object , 
                                                                                       me_cv = me_cv , 
                                                                                       n_iterations_gum = n_iterations_gum , 
                                                                                       oel = x$oel) } ) 

  # recommendation from the net: close the clusters
  stopCluster(cl)
  
  # making an array of the results
  
  simulation_result_parallel_array <- array( data = unlist(simulation_result_parallel) , dim = c( 8 , 11 , n_sim ) )
  
  # estimation of computing time ( 9 min on my computer for 5000 iterations)
  end_time <- Sys.time()
  mytime <- end_time - start_time 
  
  # results
  
  results <- list( array = simulation_result_parallel_array , time = mytime )
  
  return(results)
  
}
  

