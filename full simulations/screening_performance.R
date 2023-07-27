#' exploring the performance of screening tests in industrial hygiene measurement data interpretation

#" question : given a truth, what are the chances of the screening test yielding a "PASS"

#' principle of a screening test studied here : PASS if the maximum of a small number of measurements is lower than a fraction of the OEL


##### Data #####
  
   ## GSDs from the Kromhout et al's paper

    real_gsds <- readRDS( "created data/real_gsd_values.RDS")

    ## support functions
    
    source("other/support_functions.R")
    
##### Functions #####

    #' Performance of a screening test given a true exposure profile described by true exceedance and true GSD
    #'
    #' @param screening_threshold Fraction of the OEL used for the screening
    #' @param sample_size number of measurements in the screening procedure
    #' @param true_exceedance True exceedance fraction of the OEL
    #' @param true_gsd True geometric standard deviation
    #' @param n_simulation monte carlo simulation size,  100 000 remains fast
    #'
    #' @return proportion in % of PASS results
    #'

    
    screening.performance.sim <- function(screening_threshold = 0.1,
                                          sample_size         = 3,
                                          true_exceedance     = 0.05,
                                          true_gsd            = 3.67,
                                          n_simulation        = 100000 ) { 
      
      
      
      # generating random data
      
      true_gm <- exp(mu.from.f( f = true_exceedance , sig = log(true_gsd) , oel = 100 ))
      
      random_data <- matrix(  exp( rnorm( n_simulation*sample_size , log(true_gm) , log(true_gsd)) ) , nrow = sample_size )
      
      
      ## Looping across the n.sim iterations
      
      simulation_results <- apply( X = random_data ,
                                   MARGIN = 2,
                                   FUN = function(x){  
                                     
                                     if (max(x)<(100*screening_threshold)) return("pass") else return("fail")
                                     
                                   } ,
                                   simplify = TRUE)
      
      ## proportion of pass
      
      return(100*length(simulation_results[simulation_results=="pass"])/length(simulation_results))
      
      
    }
    
    
      
    #' Performance of a screening test given a true exposure profile described by true exceedance 
    #' GSDs are sampled for each iteration of the simulation among the dataset reported by Krouomout et al
    #'
    #' @param screening_threshold Fraction of the OEL used for the screening
    #' @param sample_size number of measurements in the screening procedure
    #' @param true_exceedance True exceedance fraction of the OEL
    #' @param n_simulation monte carlo simulation size,  100 000 remains fast
    #'
    #' @return proportion in % of PASS results
    #'
    
    
    screening.performance.sim.unknowngsd <- function(screening_threshold = 0.1,
                                                     sample_size         = 3,
                                                     true_exceedance     = 0.05,
                                                     n_simulation        = 100000 ) { 
      
      
      
      # generating random data
      
      gsd_values <- sample( real_gsds , n_simulation , replace = TRUE )
      
      random_data <- sapply( gsd_values , function(x) { exp( rnorm( n = sample_size , 
                                                                    mean = mu.from.f( f = true_exceedance , sig = log(x) , oel = 100 ) ,
                                                                    sd = log(x) ) ) }  )
      
      ## Looping across the n.sim iterations
      
      simulation_results <- apply( X = random_data ,
                                   MARGIN = 2,
                                   FUN = function(x){  
                                     
                                     if (max(x)<(100*screening_threshold)) return("pass") else return("fail")
                                     
                                   } ,
                                   simplify = TRUE)
      
      ## proportion of pass
      
      return(100*length(simulation_results[simulation_results=="pass"])/length(simulation_results))
      
      
    }
      
##### Input parameters

    screening_threshold <- 0.1
    
    sample_size <- 3
    
    true_exceedance <- 0.10
    
    true_gsd <- 2.5
    
    n_simulation <- 100000


##### Examples #####


screening.performance.sim(screening_threshold = screening_threshold,
                          sample_size         = sample_size,
                          true_exceedance     = true_exceedance,
                          true_gsd            = true_gsd,
                          n_simulation        = n_simulation )


screening.performance.sim.unknowngsd( screening_threshold = screening_threshold,
                                      sample_size         = sample_size,
                                      true_exceedance     = true_exceedance,
                                      n_simulation        = n_simulation )



