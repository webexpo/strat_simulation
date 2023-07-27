#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(utils)

shinyServer(function(input, output) {
  
  
  
  ##### DATA #####
  
  real_gsds <-readRDS(file="data/real_gsd_values.RDS")


  ##### FUNCTIONS ######
  
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
  
  
  #' estimation of log(geometric mean) from GSD and exceedance fraction
  #' 
  #'
  #' @param f : exceedance fraction in proportion, NOT in percentage
  #' @param sig : log-transformed gsd
  #' @param oel : occupational exposure limit
  
  #'
  #' @return log-geometric mean
  #'
  
  
  mu.from.f <-function(f,sig,oel)
    
    #f : exceedance fraction in proportion, NOT in percentage
    #sig : log-transformed gsd
    #oel : occupational exposure limit
    
  {
    
    if (f>=1)return(log(oel)-sig*qnorm(1-0.995))
    
    else return(log(oel)-sig*qnorm(1-f))
    
  }	
  
  
##### SERVER STUFF #####
  
  
  known_gsd_result <- reactive( screening.performance.sim(screening_threshold = input$screening_threshold,
                                                          sample_size         = input$sample_size,
                                                          true_exceedance     = input$true_exceedance/100,
                                                          true_gsd            = input$true_gsd,
                                                          n_simulation        = 100000 ) )
  
  unknown_gsd_result <- reactive( screening.performance.sim.unknowngsd(screening_threshold = input$screening_threshold,
                                                                        sample_size         = input$sample_size,
                                                                        true_exceedance     = input$true_exceedance/100,
                                                                        n_simulation        = 100000 ) )
  
  

  output$known_gsd_result  <- renderText( format( known_gsd_result() , digits=1))
  output$unknown_gsd_result  <- renderText( format( unknown_gsd_result() , digits=1))
  
})


