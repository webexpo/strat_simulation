#' Script for simulating lognormal data with or without censoring, and with and without a measurement error

# WARNING :  random numbers generated below 0.001 will be replaced by 0.001. simulation is designed for 95th percentile around 100, GM around 10. Avoid values below.

##### DATA ####


real_gsds <- readRDS( "created data/real_gsd_values.RDS")


#####  FUNCTIONS



#' simulating lognormal data for a single exposure group with or without censoring, and with or without a measurement error
#'
#' @param true_gsd 
#' @param true_gm 
#' @param sample_size 
#' @param me_cv_inf lower limit for the CV of the measurement error
#' @param me_cv_sup upper limit for the CV of the measurement error
#' @param censor_level LOQ
#' @param n_sim number of simulations 
#' @return list of 2 (character) matrices of generated data (one column per iteration), one matrix of true values, one matrix of observed values based on the coefficient of variation of measurement error
#'


data.simulation.seg <- function(  true_gsd = 2.5, 
                                  true_gm = 30, 
                                  n_sim = 50, 
                                  sample_size = 9 ,
                                  me_cv_inf = 0.15,
                                  me_cv_sup = 0.15,
                                  censor_level = 0 )  {
  
  ###### data generation
  
  
  data_vector <- exp(rnorm( n = sample_size*n_sim , mean = log(true_gm) , sd = log(true_gsd) ) )           
  
  data_vector_me <- pmax( rep(0.1,sample_size*n_sim)  , rnorm( sample_size*n_sim , data_vector , data_vector*runif( sample_size*n_sim , me_cv_inf ,  me_cv_sup ) ) ) # truncation at 0.1 to avoind negative numbers
  
  
  ###### censoring
  
  data_vector[ data_vector < censor_level ] <- paste("<", censor_level , sep="")
  
  data_vector_me[ data_vector_me < censor_level ] <- paste("<", censor_level , sep="")
  
  #formatting pre treatment
  
  data_matrix <- matrix( data = data_vector , nrow = sample_size , ncol = n_sim)
  
  data_matrix_me <- matrix( data = data_vector_me , nrow = sample_size , ncol = n_sim)
  
  return(list(true = data_matrix, 
              observed = data_matrix_me,
              true_gsd = rep(true_gsd,n_sim)
              ))
  
}



#' simulating lognormal data for a single exposure group with or without censoring, and with or without a measurement error, using the "REAL GSD" approach
#'
#' @param true_gm 
#' @param sample_size 
#' @param me_cv_inf lower limit for the CV of the measurement error
#' @param me_cv_sup upper limit for the CV of the measurement error
#' @param censor_level LOQ
#' @param n_sim number of simulations 
#' @return list of 2 (character) matrices of generated data (one column per iteration), one matrix of true values, one matrix of observed values based on the coefficient of variation of measurement error
#'
#'WARNING : This function requires object "real_gsds" to be loaded in the environment from "created data/real_gsd_values.RDS"

data.simulation.seg.realgsd <- function(  true_gm = 30, 
                                  n_sim = 50, 
                                  sample_size = 9 ,
                                  me_cv_inf = 0.15,
                                  me_cv_sup = 0.15,
                                  censor_level = 0 )  {
  
  
  ###### sampling GSD values within the 95% I of the real GSD values
  
  gsd <- sample( real_gsds[ real_gsds>=quantile(real_gsds,0.025) & real_gsds<=quantile(real_gsds,0.975)] , n_sim , replace = TRUE )
  
  ###### data generation
  
  data_matrix <- sapply( gsd , function(x) { exp( rnorm( n = sample_size , mean = log(true_gm) , sd = log(x) ) ) }  )
  
  data_matrix_me <- apply( data_matrix, 2 , function(x) { pmax( rep(0.001,sample_size)  , rnorm( sample_size , x , x*runif( sample_size , me_cv_inf ,  me_cv_sup ) ) )  }  )
  

  ###### censoring
  
  data_matrix[ data_matrix < censor_level ] <- paste("<", censor_level , sep="")
  
  data_matrix_me[ data_matrix_me < censor_level ] <- paste("<", censor_level , sep="")
  
  #formatting pre treatment
  
  return(list(true = data_matrix, 
              observed = data_matrix_me,
              true_gsd = gsd))
  
}



#' simulating lognormal data for a single exposure group with or without censoring, and with or without a measurement error, using a generic approach with vectors of input values
#'
#' @param true_gm  # vector of true GM values across the simulations
#' @param sample_size # vector of sample sizes across the simulations
#' @param me_cv_inf lower limit for the CV of the measurement error
#' @param me_cv_sup upper limit for the CV of the measurement error
#' @param censor_level LOQ # vector of LOQ values across the simulations
#' @param n_sim number of simulations 
#' @return list of 2 (character) matrices of generated data (one column per iteration), one matrix of true values, one matrix of observed values based on the coefficient of variation of measurement error
#'
#'WARNING : This function requires object "real_gsds" to be loaded in the environment from "created data/real_gsd_values.RDS"

data.simulation.seg.gen <- function(      true_gm = rep(30,10), 
                                          true_gsd = rep(2.5,10),
                                          sample_size = 9 ,
                                          me_cv_inf = 0.15,
                                          me_cv_sup = 0.15,
                                          censor_level = rep(9,10) )  {
  
  

  ###### data generation
  
  data_matrix <- exp( mapply(FUN = rnorm, n = sample_size , mean = log(true_gm) , sd = log(true_gsd) ) )
  
  data_matrix_me <- apply( data_matrix, 2 , function(x) { pmax( rep(0.001,sample_size)  , rnorm( sample_size , x , x*runif( sample_size , me_cv_inf ,  me_cv_sup ) ) )  }  )
  
  
  ###### censoring
  
  for ( i in 1:length(true_gm) ) {
    
    concentrations <- data_matrix[,i]
    
    concentrations[ concentrations < censor_level[i] ] <- paste("<", signif(censor_level[i],3) , sep="")
    
    data_matrix[,i] <- concentrations
    
    concentrations_me <- data_matrix_me[,i]
    
    concentrations_me[ concentrations_me < censor_level[i] ] <- paste("<", signif(censor_level[i],3) , sep="")
    
    data_matrix_me[,i] <- concentrations_me
    
      }
  
  
  #formatting pre treatment
  
  return(list(true = data_matrix, 
              observed = data_matrix_me,
              true_gsd = true_gsd,
              true_gm = true_gm))
  
}
