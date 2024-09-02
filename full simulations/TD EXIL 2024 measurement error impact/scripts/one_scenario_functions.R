# Functions for performing a simulation study with one scenario



##### FUNCTIONS ####

###### grouped analysis  #####

#' function which analyses the ith pair of samples ( clean and dirtied ) in the simulation 
#'
#' @param index index of the simulation 
#' @param simulated_data_object list containing the simulated data for a single scenario 
#' @param me_cv value for the measurement error as standard error CV, not in % and not expanded uncertainty CV
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


#' function which analyses the ith pair of samples ( clean and dirtied ) in the simulation : STAN version
#'
#' @param index index of the simulation 
#' @param simulated_data_object list containing the simulated data for a single scenario 
#' @param me_cv value for the measurement error as standard error CV, not in % and not expanded uncertainty CV
#' @param oel occupational exposure limit
#' @param n_iterations number of iteration for the GUM approach
#' @models.list llist object containing the models for the STAN approach
#'
#' @return matrix of results with one column per approach
#'


ithpair.function.s <- function( index , simulated_data_object , me_cv , n_iterations_gum = 10000 , oel ,
                                models.list) {
  
  # data preparation
  
  true_data <- simulated_data_object$true[,index]
  
  observed_data <- simulated_data_object$observed[,index]
  
  ## analysis
  
  ideal_b = expostats.naive.s( mysample = true_data , oel = oel , models.list = models.list)
  ideal_f = frequentist.naive( mysample = true_data , oel = oel)
  naive_b = expostats.naive.s( mysample = observed_data , oel = oel, models.list = models.list)
  naive_f = frequentist.naive( mysample = observed_data , oel = oel)
  me_b = expostats.me.s( mysample = observed_data , oel = oel , models.list = models.list , me_cv = me_cv)
  me_f = frequentist.me( mysample = observed_data , oel = oel , me_cv = me_cv , n_iterations_gum = n_iterations_gum ) # is a list of 6 vectors (mean, median,quantile)
  
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
#' @param me_cv value for the measurement error as standard error CV, not in % and not expanded uncertainty CV
#' @param oel occupational exposure limit : vector across iterations
#' @param n_iterations number of iteration for the GUM approach


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
  

#' function which uses parallel computing to perform the simulation for one data generation object (one scenario) : STAN VERSION
#'
#' @param simulated_data_object list containing the simulated data for a single scenario 
#' @param me_cv value for the measurement error as standard error CV, not in % and not expanded uncertainty CV
#' @param oel occupational exposure limit : vector across iterations
#' @param n_iterations number of iteration for the GUM approach
#' @models.list llist object containing the models for the STAN approach

#'
#' @return matrix of results with one column per approach
#'


parallel.function.s <- function( simulated_data_object , me_cv , n_iterations_gum = 10000 , n_sim , n_clusters = 10, oel , models.list) {
  
  # compteur de temps initialisé
  start_time <- Sys.time()
  
  
  #procédure parallele
  
  # creating the clusters ( use detectCores() to count the cores available, choose this number minus 2 or 4 below) 
  
  cl <- makeCluster(n_clusters)
  
  # libraries and scripts to be used in each cluster
  
  clusterEvalQ(cl, library(rstan))
  clusterEvalQ(cl, library(tolerance))
  
  # sending objects to clusters
  clusterExport( cl , "Webexpo.seg.globalbayesian.stan" , envir=environment())
  clusterExport( cl , "webexpo.seg.datapreparation" , envir=environment())
  
  clusterExport( cl , "any.me" , envir=environment())
  clusterExport( cl , "extracted.nodes" , envir=environment())
  clusterExport( cl , "webexpo.stan.inits" , envir=environment())
  clusterExport( cl , "webexpo.stan.model" , envir=environment())
  
  clusterExport( cl , "SEG.informedvar.stan" , envir=environment())
  clusterExport( cl , "models.list" , envir=environment())
  clusterExport( cl , "expostats.naive.s" , envir=environment())
  clusterExport( cl , "expostats.me.s" , envir=environment())
  
  clusterExport( cl , "frequentist.naive" , envir=environment())
  clusterExport( cl , "frequentist.me" , envir=environment())
  clusterExport( cl , "ithpair.function.s" , envir=environment())
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
  
  simulation_result_parallel <- parLapply(cl, X = my_X , function(x){ ithpair.function.s(x$index, 
                                                                                       simulated_data_object = simulated_data_object , 
                                                                                       me_cv = me_cv , 
                                                                                       n_iterations_gum = n_iterations_gum , 
                                                                                       oel = x$oel,
                                                                                       models.list = models.list) } ) 
  
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

###### separate analysis  #####

#' function which analyses the ith pair of samples ( clean and dirtied ) in the simulation for case "IDEAL - BAYESIAN" 
#'
#' @param index index of the simulation 
#' @param simulated_data_object list containing the simulated data for a single scenario 
#' @param oel occupational exposure limit
#'
#' @return vector of results 
#'



ithpair.function.ideal.b <- function( index , simulated_data_object , oel ) {
  
  # data preparation
  
  true_data <- simulated_data_object$true[,index]
  
   ## analysis
  
  ideal_b = expostats.naive( true_data , oel)
  
  ## results 
  
  results <- ideal_b 
  
  return(results)
  
}

#' function which analyses the ith pair of samples ( clean and dirtied ) in the simulation for case "IDEAL - BAYESIAN"   : Webexpo
#'
#' @param index index of the simulation 
#' @param simulated_data_object list containing the simulated data for a single scenario 
#' @param oel occupational exposure limit
#'
#' @return vector of results 
#'



ithpair.function.ideal.b.w <- function( index , simulated_data_object , oel ) {
  
  # data preparation
  
  true_data <- simulated_data_object$true[,index]
  
  ## analysis
  
  ideal_b = expostats.naive.w( true_data , oel)
  
  ## results 
  
  results <- ideal_b 
  
  return(results)
  
}

#' function which analyses the ith pair of samples ( clean and dirtied ) in the simulation for case "IDEAL - BAYESIAN"   : STAN Webexpo
#'
#' @param index index of the simulation 
#' @param simulated_data_object list containing the simulated data for a single scenario 
#' @param oel occupational exposure limit
#'
#' @return vector of results 
#'



ithpair.function.ideal.b.s <- function( index , simulated_data_object , oel ,  models.list) {
  
  # data preparation
  
  true_data <- simulated_data_object$true[,index]
  
  ## analysis
  
  ideal_b = expostats.naive.s( true_data , oel , models.list)
  
  ## results 
  
  results <- ideal_b 
  
  return(results)
  
}

#' function which analyses the ith pair of samples ( clean and dirtied ) in the simulation for case "IDEAL - FREQUENTIST" 
#'
#' @param index index of the simulation 
#' @param simulated_data_object list containing the simulated data for a single scenario 
#' @param oel occupational exposure limit
#'
#' @return vector of results 
#'



ithpair.function.ideal.f <- function( index , simulated_data_object , oel ) {
  
  # data preparation
  
  true_data <- simulated_data_object$true[,index]
  
  ## analysis
  
  ideal_f = frequentist.naive( true_data , oel)
  
  ## results 
  
  results <- ideal_f 
  
  return(results)
  
}


#' function which analyses the ith pair of samples ( clean and dirtied ) in the simulation for case "NAIVE - BAYESIAN" 
#'
#' @param index index of the simulation 
#' @param simulated_data_object list containing the simulated data for a single scenario 
#' @param oel occupational exposure limit
#'
#' @return vector of results 
#'

ithpair.function.naive.b <- function( index , simulated_data_object , oel ) {
  
  # data preparation
  
  observed_data <- simulated_data_object$observed[,index]
  
  ## analysis
  
  naive_b = expostats.naive( observed_data , oel)
  
  ## results 
  
  results <- naive_b 
  
  return(results)
  
}

#' function which analyses the ith pair of samples ( clean and dirtied ) in the simulation for case "NAIVE - BAYESIAN"  STAN
#'
#' @param index index of the simulation 
#' @param simulated_data_object list containing the simulated data for a single scenario 
#' @param oel occupational exposure limit
#'
#' @return vector of results 
#'

ithpair.function.naive.b.s <- function( index , simulated_data_object , oel , models.list) {
  
  # data preparation
  
  observed_data <- simulated_data_object$observed[,index]
  
  ## analysis
  
  naive_b = expostats.naive.s( observed_data , oel , models.list)
  
  ## results 
  
  results <- naive_b 
  
  return(results)
  
}


#' function which analyses the ith pair of samples ( clean and dirtied ) in the simulation for case "NAIVE - FREQUENTIST" 
#'
#' @param index index of the simulation 
#' @param simulated_data_object list containing the simulated data for a single scenario 
#' @param oel occupational exposure limit
#'
#' @return vector of results 
#'



ithpair.function.naive.f <- function( index , simulated_data_object , oel ) {
  
  # data preparation
  
  observed_data <- simulated_data_object$observed[,index]
  
  ## analysis
  
  naive_f = frequentist.naive( observed_data , oel)
  
  ## results 
  
  results <- naive_f 
  
  return(results)
  
}


#' function which analyses the ith pair of samples ( clean and dirtied ) in the simulation for case "ME - BAYESIAN" 
#'
#' @param index index of the simulation 
#' @param simulated_data_object list containing the simulated data for a single scenario 
#' @param me_cv value for the measurement error as standard error CV, not in % and not expanded uncertainty CV
#' @param oel occupational exposure limit

#'
#' @return vector of results 
#'



ithpair.function.me.b <- function( index , simulated_data_object , me_cv , oel ) {
  
  # data preparation
  
  observed_data <- simulated_data_object$observed[,index]
  
  ## analysis
  
  me_b = expostats.me( observed_data , oel , me_cv)
 
  ## results 
  
  results <- me_b 
  
  return(results)
  
}

#' function which analyses the ith pair of samples ( clean and dirtied ) in the simulation for case "ME - BAYESIAN"  STAN
#'
#' @param index index of the simulation 
#' @param simulated_data_object list containing the simulated data for a single scenario 
#' @param me_cv value for the measurement error as standard error CV, not in % and not expanded uncertainty CV
#' @param oel occupational exposure limit

#'
#' @return vector of results 
#'



ithpair.function.me.b.s <- function( index , simulated_data_object , me_cv , oel , models.list) {
  
  # data preparation
  
  observed_data <- simulated_data_object$observed[,index]
  
  ## analysis
  
  me_b = expostats.me.s( observed_data , oel , me_cv , models.list)
  
  ## results 
  
  results <- me_b 
  
  return(results)
  
}

#' function which analyses the ith pair of samples ( clean and dirtied ) in the simulation for case "ME - FREQUENTIST" 
#'
#' @param index index of the simulation 
#' @param simulated_data_object list containing the simulated data for a single scenario 
#' @param me_cv value for the measurement error as standard error CV, not in % and not expanded uncertainty CV
#' @param oel occupational exposure limit
#' @param n_iterations_gum number of iteration for the GUM approach
#'
#' @return matrix of results (4 quantiles of the GUM approach * 8 parameters)
#'



ithpair.function.me.f <- function( index , simulated_data_object , me_cv , n_iterations_gum = 10000 , oel ) {
  
  # data preparation
  
  observed_data <- simulated_data_object$observed[,index]
  
  ## analysis
  
  me_f = frequentist.me( observed_data , oel , me_cv , n_iterations_gum  ) # is a list of 6 vectors (mean, median,quantile)
  
  ## results 
  
  one_column <- c("gm_est","gsd_est","p95_est","p95_70ucl","p95_95ucl","F_est","F_70ucl","F_95ucl")
  
  results <- matrix( c( me_f$q2.5 , me_f$q5, me_f$q95 , me_f$q97.5) , nrow = length(one_column) )
  
  colnames(results) <- c("me_f_q2.5","me_f_q5","me_f_q95","me_f_q97.5")
  rownames(results) <- one_column
  
  
  return(results)
  
}


#' function which uses parallel computing to perform the simulation for one data generation object (one scenario): approach IDEAL BAYESIAN
#'
#' @param simulated_data_object list containing the simulated data for a single scenario 
#' @param oel occupational exposure limit : vector across iterations
#' @param n_sim number of simulation iteration
#' @param n_clusters number of clusters for parallel computing

#'
#' @return matrix of results with one column per approach
#'


parallel.function.ideal.b <- function( simulated_data_object , n_sim , n_clusters = 10, oel ) {
  
  # compteur de temps initialisé
  start_time <- Sys.time()
  
  
  #procédure parallele
  
  # creating the clusters ( use detectCores() to count the cores available, choose this number minus 2 or 4 below) 
  
  cl <- makeCluster(n_clusters)
  
  # libraries and scripts to be used in each cluster
  
  clusterEvalQ(cl, library(rjags))

  # sending objects to clusters
  clusterExport( cl , "jags.model.informedvar" , envir=environment())
  clusterExport( cl , "webexpo.seg.datapreparation" , envir=environment())
  clusterExport( cl , "Webexpo.seg.globalbayesian.jags" , envir=environment())
  clusterExport( cl , "fun.jags.informedvar" , envir=environment())
  clusterExport( cl , "expostats.naive" , envir=environment())
  clusterExport( cl , "ithpair.function.ideal.b" , envir=environment())

  clusterExport( cl , "oel" , envir=environment())
  clusterExport( cl , "simulated_data_object" , envir=environment())
  
  # list for the LApply function
  
  my_X <- vector("list", length = n_sim)
  
  for ( i in 1:n_sim ) { my_X[[i]] <- list( index = i,
                                            oel = oel[i]) }
  
  
  # calculations
  
  simulation_result_parallel <- parLapply(cl, X = my_X , function(x){ ithpair.function.ideal.b(x$index, 
                                                                                       simulated_data_object = simulated_data_object , 
                                                                                       oel = x$oel) } ) 
  
  # recommendation from the net: close the clusters
  stopCluster(cl)
  
  # making an matrix of the results
  
  simulation_result_parallel_matrix <- matrix( data = unlist(simulation_result_parallel) , nrow =  8  , ncol = n_sim ) 
  
  # estimation of computing time ( 9 min on my computer for 5000 iterations)
  end_time <- Sys.time()
  mytime <- end_time - start_time 
  
  # results
  
  results <- list( matrix = simulation_result_parallel_matrix , time = mytime )
  
  return(results)
  
}

#' function which uses parallel computing to perform the simulation for one data generation object (one scenario): approach IDEAL BAYESIAN : WEBEXPO
#'
#' @param simulated_data_object list containing the simulated data for a single scenario 
#' @param oel occupational exposure limit : vector across iterations
#' @param n_sim number of simulation iteration
#' @param n_clusters number of clusters for parallel computing

#'
#' @return matrix of results with one column per approach
#'


parallel.function.ideal.b.w <- function( simulated_data_object , n_sim , n_clusters = 10, oel ) {
  
  # compteur de temps initialisé
  start_time <- Sys.time()
  
  
  #procédure parallele
  
  # creating the clusters ( use detectCores() to count the cores available, choose this number minus 2 or 4 below) 
  
  cl <- makeCluster(n_clusters)
  
  # libraries and scripts to be used in each cluster
  
  clusterEvalQ(cl, library(rjags))
  
  # sending objects to clusters
  clusterExport( cl , "Webexpo.seg.globalbayesian.mcgill" , envir=environment())
  clusterExport( cl , "webexpo.seg.datapreparation" , envir=environment())
  clusterExport( cl , "data.summary" , envir=environment())
  clusterExport( cl , "dens.gen.icdf" , envir=environment())
  clusterExport( cl , list("any.me","cond.values","Default.inits","empty.matrices","logp.from.logpcum",
                           "logPhi.quadratic.approx.coeff","lphi","me.gen","me.gen.object","mu.truncatedData.gen",
                           "mu.truncatedData.gen.object","out.logout.moments","out.sample","pnorm.logpcum",
                           "polynomial.product.coeff","quadratic.solution","random.pow","real.cubic.roots",
                           "renamed.me.parm","rgamma.truncated","rnorm.censored","runif.logp",
                           "sigma.gen.object","sigma.truncatedData.gen","sigma.truncatedData.gen.object",
                           "sqrt.invertedGamma.gen","take.log","truevalue.gen","truevalue.gen.object",
                           "truevalues.gen","Xinverse.y","y.gen","y.gen.inits") , envir=environment())
  
  
  
  clusterExport( cl , "SEG.informedvar" , envir=environment())
  
  #clusterExport( cl , "cond.values" , envir=environment())
  #clusterExport( cl , "cond.values" , envir=environment())
  
  
  clusterExport( cl , "expostats.naive.w" , envir=environment())
  clusterExport( cl , "ithpair.function.ideal.b.w" , envir=environment())
  
  clusterExport( cl , "oel" , envir=environment())
  clusterExport( cl , "simulated_data_object" , envir=environment())
  
  # list for the LApply function
  
  my_X <- vector("list", length = n_sim)
  
  for ( i in 1:n_sim ) { my_X[[i]] <- list( index = i,
                                            oel = oel[i]) }
  
  
  # calculations
  
  simulation_result_parallel <- parLapply(cl, X = my_X , function(x){ ithpair.function.ideal.b.w(x$index, 
                                                                                               simulated_data_object = simulated_data_object , 
                                                                                               oel = x$oel) } ) 
  
  # recommendation from the net: close the clusters
  stopCluster(cl)
  
  # making an matrix of the results
  
  simulation_result_parallel_matrix <- matrix( data = unlist(simulation_result_parallel) , nrow =  8  , ncol = n_sim ) 
  
  # estimation of computing time ( 9 min on my computer for 5000 iterations)
  end_time <- Sys.time()
  mytime <- end_time - start_time 
  
  # results
  
  results <- list( matrix = simulation_result_parallel_matrix , time = mytime )
  
  return(results)
  
}


#' function which uses parallel computing to perform the simulation for one data generation object (one scenario): approach IDEAL BAYESIAN : STAN WEBEXPO
#'
#' @param simulated_data_object list containing the simulated data for a single scenario 
#' @param oel occupational exposure limit : vector across iterations
#' @param n_sim number of simulation iteration
#' @param n_clusters number of clusters for parallel computing

#'
#' @return matrix of results with one column per approach
#'


parallel.function.ideal.b.s <- function( simulated_data_object , n_sim , n_clusters = 10, oel, models.list ) {
  
  # compteur de temps initialisé
  start_time <- Sys.time()
  
  
  #procédure parallele
  
  # creating the clusters ( use detectCores() to count the cores available, choose this number minus 2 or 4 below) 
  
  cl <- makeCluster(n_clusters)
  
  # libraries and scripts to be used in each cluster
  
  clusterEvalQ(cl, library(rstan))
  
  # sending objects to clusters
  
  clusterExport( cl , "Webexpo.seg.globalbayesian.stan" , envir=environment())
  clusterExport( cl , "webexpo.seg.datapreparation" , envir=environment())
  
  clusterExport( cl , "any.me" , envir=environment())
  clusterExport( cl , "extracted.nodes" , envir=environment())
  clusterExport( cl , "webexpo.stan.inits" , envir=environment())
  clusterExport( cl , "webexpo.stan.model" , envir=environment())

  clusterExport( cl , "SEG.informedvar.stan" , envir=environment())
  clusterExport( cl , "models.list" , envir=environment())
  
  clusterExport( cl , "expostats.naive.s" , envir=environment())
  clusterExport( cl , "ithpair.function.ideal.b.s" , envir=environment())
  
  clusterExport( cl , "oel" , envir=environment())
  clusterExport( cl , "simulated_data_object" , envir=environment())
  
  # list for the LApply function
  
  my_X <- vector("list", length = n_sim)
  
  for ( i in 1:n_sim ) { my_X[[i]] <- list( index = i,
                                            oel = oel[i]) }
  
  
  # calculations
  
  simulation_result_parallel <- parLapply(cl, X = my_X , function(x){ ithpair.function.ideal.b.s(x$index, 
                                                                                                 simulated_data_object = simulated_data_object , 
                                                                                                 oel = x$oel,
                                                                                                 models.list = models.list) } ) 
  
  # recommendation from the net: close the clusters
  stopCluster(cl)
  
  # making an matrix of the results
  
  simulation_result_parallel_matrix <- matrix( data = unlist(simulation_result_parallel) , nrow =  8  , ncol = n_sim ) 
  
  # estimation of computing time ( 9 min on my computer for 5000 iterations)
  end_time <- Sys.time()
  mytime <- end_time - start_time 
  
  # results
  
  results <- list( matrix = simulation_result_parallel_matrix , time = mytime )
  
  return(results)
  
}

#' function which uses parallel computing to perform the simulation for one data generation object (one scenario): approach NAIVE BAYESIAN
#'
#' @param simulated_data_object list containing the simulated data for a single scenario 
#' @param oel occupational exposure limit : vector across iterations
#' @param n_sim number of simulation iteration
#' @param n_clusters number of clusters for parallel computing



parallel.function.naive.b <- function( simulated_data_object , n_sim , n_clusters = 10, oel ) {
  
  # compteur de temps initialisé
  start_time <- Sys.time()
  
  
  #procédure parallele
  
  # creating the clusters ( use detectCores() to count the cores available, choose this number minus 2 or 4 below) 
  
  cl <- makeCluster(n_clusters)
  
  # libraries and scripts to be used in each cluster
  
  clusterEvalQ(cl, library(rjags))
  
  # sending objects to clusters
  clusterExport( cl , "jags.model.informedvar" , envir=environment())
  clusterExport( cl , "webexpo.seg.datapreparation" , envir=environment())
  clusterExport( cl , "Webexpo.seg.globalbayesian.jags" , envir=environment())
  clusterExport( cl , "fun.jags.informedvar" , envir=environment())
  clusterExport( cl , "expostats.naive" , envir=environment())
  clusterExport( cl , "ithpair.function.naive.b" , envir=environment())
  
  clusterExport( cl , "oel" , envir=environment())
  clusterExport( cl , "simulated_data_object" , envir=environment())
  
  # list for the LApply function
  
  my_X <- vector("list", length = n_sim)
  
  for ( i in 1:n_sim ) { my_X[[i]] <- list( index = i,
                                            oel = oel[i]) }
  
  
  # calculations
  
  simulation_result_parallel <- parLapply(cl, X = my_X , function(x){ ithpair.function.naive.b(x$index, 
                                                                                               simulated_data_object = simulated_data_object , 
                                                                                               oel = x$oel) } ) 
  
  # recommendation from the net: close the clusters
  stopCluster(cl)
  
  # making an matrix of the results
  
  simulation_result_parallel_matrix <- matrix( data = unlist(simulation_result_parallel) , nrow =  8  , ncol = n_sim ) 
  
  # estimation of computing time ( 9 min on my computer for 5000 iterations)
  end_time <- Sys.time()
  mytime <- end_time - start_time 
  
  # results
  
  results <- list( matrix = simulation_result_parallel_matrix , time = mytime )
  
  return(results)
  
}


#' function which uses parallel computing to perform the simulation for one data generation object (one scenario): approach NAIVE BAYESIAN : STAN WEBEXPO
#'
#' @param simulated_data_object list containing the simulated data for a single scenario 
#' @param oel occupational exposure limit : vector across iterations
#' @param n_sim number of simulation iteration
#' @param n_clusters number of clusters for parallel computing

#'
#' @return matrix of results with one column per approach
#'


parallel.function.naive.b.s <- function( simulated_data_object , n_sim , n_clusters = 10, oel, models.list ) {
  
  # compteur de temps initialisé
  start_time <- Sys.time()
  
  
  #procédure parallele
  
  # creating the clusters ( use detectCores() to count the cores available, choose this number minus 2 or 4 below) 
  
  cl <- makeCluster(n_clusters)
  
  # libraries and scripts to be used in each cluster
  
  clusterEvalQ(cl, library(rstan))
  
  # sending objects to clusters
  
  clusterExport( cl , "Webexpo.seg.globalbayesian.stan" , envir=environment())
  clusterExport( cl , "webexpo.seg.datapreparation" , envir=environment())
  
  clusterExport( cl , "any.me" , envir=environment())
  clusterExport( cl , "extracted.nodes" , envir=environment())
  clusterExport( cl , "webexpo.stan.inits" , envir=environment())
  clusterExport( cl , "webexpo.stan.model" , envir=environment())
  
  clusterExport( cl , "SEG.informedvar.stan" , envir=environment())
  clusterExport( cl , "models.list" , envir=environment())
  
  clusterExport( cl , "expostats.naive.s" , envir=environment())
  clusterExport( cl , "ithpair.function.naive.b.s" , envir=environment())
  
  clusterExport( cl , "oel" , envir=environment())
  clusterExport( cl , "simulated_data_object" , envir=environment())
  
  # list for the LApply function
  
  my_X <- vector("list", length = n_sim)
  
  for ( i in 1:n_sim ) { my_X[[i]] <- list( index = i,
                                            oel = oel[i]) }
  
  
  # calculations
  
  simulation_result_parallel <- parLapply(cl, X = my_X , function(x){ ithpair.function.naive.b.s(x$index, 
                                                                                                 simulated_data_object = simulated_data_object , 
                                                                                                 oel = x$oel,
                                                                                                 models.list = models.list) } ) 
  
  # recommendation from the net: close the clusters
  stopCluster(cl)
  
  # making an matrix of the results
  
  simulation_result_parallel_matrix <- matrix( data = unlist(simulation_result_parallel) , nrow =  8  , ncol = n_sim ) 
  
  # estimation of computing time ( 9 min on my computer for 5000 iterations)
  end_time <- Sys.time()
  mytime <- end_time - start_time 
  
  # results
  
  results <- list( matrix = simulation_result_parallel_matrix , time = mytime )
  
  return(results)
  
}



#' function which uses parallel computing to perform the simulation for one data generation object (one scenario): approach ME BAYESIAN
#'
#' @param simulated_data_object list containing the simulated data for a single scenario 
#' @param oel occupational exposure limit : vector across iterations
#' @param me_cv value for the measurement error as standard error CV, not in % and not expanded uncertainty CV
#' @param n_sim number of simulation iteration
#' @param n_clusters number of clusters for parallel computing



parallel.function.me.b <- function( simulated_data_object , me_cv , n_sim , n_clusters = 10, oel ) {
  
  # compteur de temps initialisé
  start_time <- Sys.time()
  
  
  #procédure parallele
  
  # creating the clusters ( use detectCores() to count the cores available, choose this number minus 2 or 4 below) 
  
  cl <- makeCluster(n_clusters)
  
  # libraries and scripts to be used in each cluster
  
  clusterEvalQ(cl, library(rjags))
  
  # sending objects to clusters
  clusterExport( cl , "jags.model.informedvar" , envir=environment())
  clusterExport( cl , "webexpo.seg.datapreparation" , envir=environment())
  clusterExport( cl , "Webexpo.seg.globalbayesian.jags" , envir=environment())
  clusterExport( cl , "fun.jags.informedvar" , envir=environment())
  clusterExport( cl , "expostats.me" , envir=environment())
  clusterExport( cl , "ithpair.function.me.b" , envir=environment())
  
  clusterExport( cl , "oel" , envir=environment())
  clusterExport( cl , "me_cv" , envir=environment())
  clusterExport( cl , "simulated_data_object" , envir=environment())
  
  # list for the LApply function
  
  my_X <- vector("list", length = n_sim)
  
  for ( i in 1:n_sim ) { my_X[[i]] <- list( index = i,
                                            oel = oel[i]) }
  
  
  # calculations
  
  simulation_result_parallel <- parLapply(cl, X = my_X , function(x){ ithpair.function.me.b(x$index, 
                                                                                               simulated_data_object = simulated_data_object , 
                                                                                               me_cv=me_cv,
                                                                                               oel = x$oel) } ) 
  
  # recommendation from the net: close the clusters
  stopCluster(cl)
  
  # making an matrix of the results
  
  simulation_result_parallel_matrix <- matrix( data = unlist(simulation_result_parallel) , nrow =  8  , ncol = n_sim ) 
  
  # estimation of computing time ( 9 min on my computer for 5000 iterations)
  end_time <- Sys.time()
  mytime <- end_time - start_time 
  
  # results
  
  results <- list( matrix = simulation_result_parallel_matrix , time = mytime )
  
  return(results)
  
}


#' function which uses parallel computing to perform the simulation for one data generation object (one scenario): approach ME BAYESIAN : STAN WEBEXPO
#'
#' @param simulated_data_object list containing the simulated data for a single scenario 
#' @param oel occupational exposure limit : vector across iterations
#' @param n_sim number of simulation iteration
#' @param n_clusters number of clusters for parallel computing

#'
#' @return matrix of results with one column per approach
#'


parallel.function.me.b.s <- function( simulated_data_object , me_cv, n_sim , n_clusters = 10, oel, models.list ) {
  
  # compteur de temps initialisé
  start_time <- Sys.time()
  
  
  #procédure parallele
  
  # creating the clusters ( use detectCores() to count the cores available, choose this number minus 2 or 4 below) 
  
  cl <- makeCluster(n_clusters)
  
  # libraries and scripts to be used in each cluster
  
  clusterEvalQ(cl, library(rstan))
  
  # sending objects to clusters
  
  clusterExport( cl , "Webexpo.seg.globalbayesian.stan" , envir=environment())
  clusterExport( cl , "webexpo.seg.datapreparation" , envir=environment())
  
  clusterExport( cl , "any.me" , envir=environment())
  clusterExport( cl , "extracted.nodes" , envir=environment())
  clusterExport( cl , "webexpo.stan.inits" , envir=environment())
  clusterExport( cl , "webexpo.stan.model" , envir=environment())
  
  clusterExport( cl , "SEG.informedvar.stan" , envir=environment())
  clusterExport( cl , "models.list" , envir=environment())
  
  clusterExport( cl , "expostats.me.s" , envir=environment())
  clusterExport( cl , "ithpair.function.me.b.s" , envir=environment())
  
  clusterExport( cl , "oel" , envir=environment())
  clusterExport( cl , "me_cv" , envir=environment())
  clusterExport( cl , "simulated_data_object" , envir=environment())
  
  # list for the LApply function
  
  my_X <- vector("list", length = n_sim)
  
  for ( i in 1:n_sim ) { my_X[[i]] <- list( index = i,
                                            oel = oel[i]) }
  
  
  # calculations
  
  simulation_result_parallel <- parLapply(cl, X = my_X , function(x){ ithpair.function.me.b.s(x$index, 
                                                                                                 simulated_data_object = simulated_data_object , 
                                                                                                 me_cv=me_cv,
                                                                                                 oel = x$oel,
                                                                                                 models.list = models.list) } ) 
  
  # recommendation from the net: close the clusters
  stopCluster(cl)
  
  # making an matrix of the results
  
  simulation_result_parallel_matrix <- matrix( data = unlist(simulation_result_parallel) , nrow =  8  , ncol = n_sim ) 
  
  # estimation of computing time ( 9 min on my computer for 5000 iterations)
  end_time <- Sys.time()
  mytime <- end_time - start_time 
  
  # results
  
  results <- list( matrix = simulation_result_parallel_matrix , time = mytime )
  
  return(results)
  
}


#' function which uses parallel computing to perform the simulation for one data generation object (one scenario) IDEAL FREQUENTIST
#'
#' @param simulated_data_object list containing the simulated data for a single scenario 
#' @param oel occupational exposure limit : vector across iterations
#' @param n_sim number of iteration for the simulation

#'
#' @return matrix of results with one column per approach
#'


parallel.function.ideal.f <- function( simulated_data_object , n_sim , n_clusters = 10, oel ) {
  
  # compteur de temps initialisé
  start_time <- Sys.time()
  
  
  #procédure parallele
  
  # creating the clusters ( use detectCores() to count the cores available, choose this number minus 2 or 4 below) 
  
  cl <- makeCluster(n_clusters)
  
  # libraries and scripts to be used in each cluster
  
  clusterEvalQ(cl, library(tolerance))
  
  # sending objects to clusters
  clusterExport( cl , "frequentist.naive" , envir=environment())
  clusterExport( cl , "ithpair.function.ideal.f" , envir=environment())
  clusterExport( cl , "fun.NdExpo.lognorm" , envir=environment())
  clusterExport( cl , "fun.perc.en689" , envir=environment())
  clusterExport( cl , "fun.frac.dep" , envir=environment())
  
  
  clusterExport( cl , "oel" , envir=environment())
  clusterExport( cl , "simulated_data_object" , envir=environment())
  
  # list for the LApply function
  
  my_X <- vector("list", length = n_sim)
  
  for ( i in 1:n_sim ) { my_X[[i]] <- list( index = i,
                                            oel = oel[i]) }
  
  
  # calculations
  
  simulation_result_parallel <- parLapply(cl, X = my_X , function(x){ ithpair.function.ideal.f(x$index, 
                                                                                       simulated_data_object = simulated_data_object , 
                                                                                       oel = x$oel) } ) 
  
  # recommendation from the net: close the clusters
  stopCluster(cl)
  
  # making an matrix of the results
  
  simulation_result_parallel_matrix <- matrix( data = unlist(simulation_result_parallel) , nrow =  8  , ncol = n_sim ) 
  
  # estimation of computing time ( 9 min on my computer for 5000 iterations)
  end_time <- Sys.time()
  mytime <- end_time - start_time 
  
  # results
  
  results <- list( matrix = simulation_result_parallel_matrix , time = mytime )
  
  
  return(results)
  
  
}


#' function which uses parallel computing to perform the simulation for one data generation object (one scenario) NAIVE FREQUENTIST
#'
#' @param simulated_data_object list containing the simulated data for a single scenario 
#' @param oel occupational exposure limit : vector across iterations
#' @param n_sim number of iteration for the simulation

#'
#' @return matrix of results with one column per approach
#'


parallel.function.naive.f <- function( simulated_data_object , n_sim , n_clusters = 10, oel ) {
  
  # compteur de temps initialisé
  start_time <- Sys.time()
  
  
  #procédure parallele
  
  # creating the clusters ( use detectCores() to count the cores available, choose this number minus 2 or 4 below) 
  
  cl <- makeCluster(n_clusters)
  
  # libraries and scripts to be used in each cluster
  
  clusterEvalQ(cl, library(tolerance))
  
  # sending objects to clusters
  clusterExport( cl , "frequentist.naive" , envir=environment())
  clusterExport( cl , "ithpair.function.naive.f" , envir=environment())
  clusterExport( cl , "fun.NdExpo.lognorm" , envir=environment())
  clusterExport( cl , "fun.perc.en689" , envir=environment())
  clusterExport( cl , "fun.frac.dep" , envir=environment())
  
  
  clusterExport( cl , "oel" , envir=environment())
  clusterExport( cl , "simulated_data_object" , envir=environment())
  
  # list for the LApply function
  
  my_X <- vector("list", length = n_sim)
  
  for ( i in 1:n_sim ) { my_X[[i]] <- list( index = i,
                                            oel = oel[i]) }
  
  
  # calculations
  
  simulation_result_parallel <- parLapply(cl, X = my_X , function(x){ ithpair.function.naive.f(x$index, 
                                                                                               simulated_data_object = simulated_data_object , 
                                                                                               oel = x$oel) } ) 
  
  # recommendation from the net: close the clusters
  stopCluster(cl)
  
  # making an matrix of the results
  
  simulation_result_parallel_matrix <- matrix( data = unlist(simulation_result_parallel) , nrow =  8  , ncol = n_sim ) 
  
  # estimation of computing time ( 9 min on my computer for 5000 iterations)
  end_time <- Sys.time()
  mytime <- end_time - start_time 
  
  # results
  
 results <- list( matrix = simulation_result_parallel_matrix , time = mytime )
  
  
  return(results)
  
  
}


#' function which uses parallel computing to perform the simulation for one data generation object (one scenario) ME FREQUENTIST
#'
#' @param simulated_data_object list containing the simulated data for a single scenario 
#' @param me_cv value for the measurement error as standard error CV, not in % and not expanded uncertainty CV
#' @param oel occupational exposure limit : vector across iterations
#' @param n_iterations_gum number of iteration for the GUM approach
#' @param n_sim number of iteration for the simulation

#'
#' @return matrix of results with one column per approach
#'


parallel.function.me.f <- function( simulated_data_object , me_cv , n_iterations_gum = 10000 , n_sim , n_clusters = 10, oel ) {
  
  # compteur de temps initialisé
  start_time <- Sys.time()
  
  
  #procédure parallele
  
  # creating the clusters ( use detectCores() to count the cores available, choose this number minus 2 or 4 below) 
  
  cl <- makeCluster(n_clusters)
  
  # libraries and scripts to be used in each cluster
  
  clusterEvalQ(cl, library(tolerance))
  
  # sending objects to clusters
  clusterExport( cl , "frequentist.me" , envir=environment())
  clusterExport( cl , "ithpair.function.me.f" , envir=environment())
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
  
  simulation_result_parallel <- parLapply(cl, X = my_X , function(x){ ithpair.function.me.f(x$index, 
                                                                                       simulated_data_object = simulated_data_object , 
                                                                                       me_cv = me_cv , 
                                                                                       n_iterations_gum = n_iterations_gum , 
                                                                                       oel = x$oel) } ) 
  
  # recommendation from the net: close the clusters
  stopCluster(cl)
  
  # making an array of the results
  
  simulation_result_parallel_array <- array( data = unlist(simulation_result_parallel) , dim = c( 8 , 4 , n_sim ) )
  
  # estimation of computing time ( 9 min on my computer for 5000 iterations)
  end_time <- Sys.time()
  mytime <- end_time - start_time 
  
  # results
  
  results <- list( array = simulation_result_parallel_array , time = mytime )

  return(results)
  
}


