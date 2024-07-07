
#'   Support functions for sampling strat simulations



################### estimation of log(geometric mean) from GSD and exceedance fraction


#' Calculation of ln(GM) as a function of the exceedance fraction, ln(GSD) and the OEL
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
        
        
        
 
######## transfer between 95th percentile and exceedance fraction
       

#' Calculation of 95th percentile as a function of the exceedance fraction, the GSD and the OEL
#'
#' @param frac : exceedance fraction in proportion, NOT in percentage
#' @param gsd : log-transformed gsd 
#' @param oel : occupational exposure limit 
#'
#' @return value of the 95th percentile
#'
       
       
             
        p95fromF <- function( oel = 1 , gsd = 1.5 , frac = 0.05) {
          
          
          return( oel*exp( log(gsd) * ( qnorm(0.95) - qnorm(1-frac) ) ) )  
          
          
        }
        
#' Calculation of exceedance fraction as a function of the 95th percentile, the GSD and the OEL
#'
#' @param p95 : 95th percentile
#' @param gsd : log-transformed gsd 
#' @param oel : occupational exposure limit 
#'
#' @return exceedance fraction in proportion, NOT in percentage
#'
        
        
        Ffromp95 <- function( oel = 1 , gsd = 1.5 , p95 = 0.05) {
          
          
          return( 1-pnorm( (qnorm(0.95)*log(gsd) - log(p95/oel))/log(gsd) )   )  
          
          
        }
        
     
  