#' Script for simulating lognormal data with or without censoring, and with or without a measurement error


#' simulating lognormal data for a single exposure group with or without censoring, and with or without a measurement error
#'
#' @param true_gsd 
#' @param true_gm 
#' @param sample_size 
#' @param me_cv_inf lower limit for the CV of the measurement error
#' @param me_cv_sup upper limit for the CV of the measurement error
#' @param prop_censor proportion of censorship
#' @param n_sim number of simulations 
#' @return list of 2 matrices of generated data (one column per iteration), one matrix of true values, one matrix of observed values based on the coefficient of variation of measurement error
#'


data.simulation.seg <- function(  true_gsd = 2.5, 
                                  true_gm = 30, 
                                  n_sim = 50, 
                                  sample_size = 9 ,
                                  me_cv_inf = 0.15,
                                  me_cv_sup = 0.15,
                                  prop_censor = 0.1 )  {
  
  ###### data generation
  
  
  data_vector <- exp(rnorm( n = sample_size*n_sim , mean = log(true_gm) , sd = log(true_gsd) ) )           
  
  data_vector_me <- pmax( rep(0.1,sample_size*n_sim)  , rnorm( sample_size*n_sim , data_vector , data_vector*me.cv ) )
  
  #formatting pre treatment
  
  data_matrix <- matrix( data = data_vector , nrow = sample_size , ncol = n_sim)
  
  data_matrix_me <- matrix( data = data_vector_me , nrow = sample_size , ncol = n_sim)
  
  return(list(true = data_matrix, observed = data_matrix_me))
  
}

