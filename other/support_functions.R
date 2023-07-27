#########################################################################
#
#
#   Support functions for the revisited pilot study for Santé Canada project
#
#
###########################################################################


################### estimation of log(geometric mean) from GSD and exceedance fraction

        mu.from.f <-function(f,sig,oel)
          
          #f : exceedance fraction in proportion, NOT in percentage
          #sig : log-transformed gsd
          #oel : occupational exposure limit
          
        {
          
          if (f>=1)return(log(oel)-sig*qnorm(1-0.995))
          
          else return(log(oel)-sig*qnorm(1-f))
          
        }	
        
        
        
 
######## transfer betweem 95th perc and exceedance        
             
        p95fromF <- function( oel = 1 , gsd = 1.5 , frac = 0.05) {
          
          
          return( oel*exp( log(gsd) * ( qnorm(0.95) - qnorm(1-frac) ) ) )  
          
          
        }
        
        
        Ffromp95 <- function( oel = 1 , gsd = 1.5 , p95 = 0.05) {
          
          
          return( 1-pnorm( (qnorm(0.95)*log(gsd) - log(p95/oel))/log(gsd) )   )  
          
          
        }
        
        Ffromp95(oel = 1 , gsd = 2.5 , p95 = 0.5)
        