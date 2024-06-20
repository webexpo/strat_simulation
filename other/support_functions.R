
#'   Support functions for sampling strat simulations



################### estimation of log(geometric mean) from GSD and exceedance fraction

        mu.from.f <-function(f,sig,oel)
          
          #f : exceedance fraction in proportion, NOT in percentage
          #sig : log-transformed gsd
          #oel : occupational exposure limit
          
        {
          
          if (f>=1)return(log(oel)-sig*qnorm(1-0.995))
          
          else return(log(oel)-sig*qnorm(1-f))
          
        }	
        
        
        
 
######## transfer between 95th perc and exceedance        
             
        p95fromF <- function( oel = 1 , gsd = 1.5 , frac = 0.05) {
          
          
          return( oel*exp( log(gsd) * ( qnorm(0.95) - qnorm(1-frac) ) ) )  
          
          
        }
        
        
        Ffromp95 <- function( oel = 1 , gsd = 1.5 , p95 = 0.05) {
          
          
          return( 1-pnorm( (qnorm(0.95)*log(gsd) - log(p95/oel))/log(gsd) )   )  
          
          
        }
        
        Ffromp95(oel = 1 , gsd = 2.5 , p95 = 0.5)
        
        
        
        compute_rmse = function(x, theta){
          
          mean.x = mean(x)
          
          N = length(x)
          
          if(is.null(N))
            stop("x must be a vector of minimum length 1")
          
          rmse = sqrt((mean.x - theta) ^ 2 + sum((x - mean.x) ^ 2) / (N - 1))
          return(rmse)
        }
        
        compute_relrmse = function(x, theta){
          
          mean.x = mean(x)
          
          N = length(x)
          
          if(is.null(N))
            stop("x must be a vector of minimum length 1")
          
          rmse = 100*sqrt((mean.x - theta) ^ 2 + sum((x - mean.x) ^ 2) / (N - 1))/theta
          return(rmse)
        }
        
        
        compute_bias = function(x, theta){
          
          bias = mean(x) - theta 
          
          return(bias)
        }
        
        
        compute_relbias = function(x, theta){
          
          bias = 100 * (mean(x) - theta) / theta
          
          return(bias)
        }
        
        
        compute_relprecision = function(x, theta){
          
          sd.x = sd(x)
          
          rp = 100 * sd.x / theta
          
          return(rp)
        }
        