# Script describing functions used to analyse one sample with all approaches


###### DATA ####

##### FUNCTIONS ####

#' Naive bayesian analysis using the expostats prior
#'
#' @param mysample sample to be analysed
#' @param oel occupational exposure limit
#'
#' @return point estimates for GM, GSD, exceedance and P95, and 70 and 95 UCL for P95 and exceedance
#'

expostats.naive <- function( mysample , oel) {
  
  mcmc <- suppressWarnings( Webexpo.seg.globalbayesian.jags( data.sample = mysample ,
                                           is.lognormal = TRUE , 
                                           error.type = "none" ,
                                           me.range = c(0.3,0.3) , 
                                           oel = oel ,
                                           prior.model = "informedvar",
                                           n.iter = 25000))
  
  results <- c( gm_est = exp(median(mcmc$mu.chain)),
                gsd_est = exp(median(mcmc$sigma.chain)),
                p95_est = median(exp(mcmc$mu.chain+qnorm(0.95)*mcmc$sigma.chain)),
                p95_70ucl = quantile(exp(mcmc$mu.chain+qnorm(0.95)*mcmc$sigma.chain),0.7),
                p95_95ucl = quantile(exp(mcmc$mu.chain+qnorm(0.95)*mcmc$sigma.chain),0.95),
                F_est = unname(quantile(100*(1 - pnorm( log(oel) , mean = mcmc$mu.chain , sd = mcmc$sigma.chain )),0.5)),
                F_70ucl = unname(quantile(100*(1 - pnorm( log(oel) , mean = mcmc$mu.chain , sd = mcmc$sigma.chain )),0.7)),
                F_95ucl = unname(quantile(100*(1 - pnorm( log(oel) , mean = mcmc$mu.chain , sd = mcmc$sigma.chain )),0.95)))
  
  
  return(results)
  
}

#' Naive bayesian analysis using the expostats prior and the webexpo algorithm
#'
#' @param mysample sample to be analysed
#' @param oel occupational exposure limit
#'
#' @return point estimates for GM, GSD, exceedance and P95, and 70 and 95 UCL for P95 and exceedance
#'

expostats.naive.w <- function( mysample , oel) {
  
  mcmc <- suppressWarnings( Webexpo.seg.globalbayesian.mcgill( data.sample = mysample ,
                                                             is.lognormal = TRUE , 
                                                             error.type = "none" ,
                                                             me.range = c(0.3,0.3) , 
                                                             oel = oel ,
                                                             prior.model = "informedvar",
                                                             n.iter = 25000))
  
  results <- c( gm_est = exp(median(mcmc$mu.chain)),
                gsd_est = exp(median(mcmc$sigma.chain)),
                p95_est = median(exp(mcmc$mu.chain+qnorm(0.95)*mcmc$sigma.chain)),
                p95_70ucl = quantile(exp(mcmc$mu.chain+qnorm(0.95)*mcmc$sigma.chain),0.7),
                p95_95ucl = quantile(exp(mcmc$mu.chain+qnorm(0.95)*mcmc$sigma.chain),0.95),
                F_est = unname(quantile(100*(1 - pnorm( log(oel) , mean = mcmc$mu.chain , sd = mcmc$sigma.chain )),0.5)),
                F_70ucl = unname(quantile(100*(1 - pnorm( log(oel) , mean = mcmc$mu.chain , sd = mcmc$sigma.chain )),0.7)),
                F_95ucl = unname(quantile(100*(1 - pnorm( log(oel) , mean = mcmc$mu.chain , sd = mcmc$sigma.chain )),0.95)))
  
  
  return(results)
  
}


#' Naive bayesian analysis using the expostats prior and the webexpo STAN algorithm
#'
#' @param mysample sample to be analysed
#' @param oel occupational exposure limit
#'
#' @return point estimates for GM, GSD, exceedance and P95, and 70 and 95 UCL for P95 and exceedance
#'

expostats.naive.s <- function( mysample , oel , models.list) {
  
  mcmc <- suppressWarnings( Webexpo.seg.globalbayesian.stan( data.sample = mysample ,
                                                               is.lognormal = TRUE , 
                                                               error.type = "none" ,
                                                               me.range = c(0.3,0.3) , 
                                                               oel = oel ,
                                                               prior.model = "informedvar",
                                                               n.iter = 25000,
                                                               models.list=models.list))
  
  results <- c( gm_est = exp(median(mcmc$mu.chain)),
                gsd_est = exp(median(mcmc$sigma.chain)),
                p95_est = median(exp(mcmc$mu.chain+qnorm(0.95)*mcmc$sigma.chain)),
                p95_70ucl = quantile(exp(mcmc$mu.chain+qnorm(0.95)*mcmc$sigma.chain),0.7),
                p95_95ucl = quantile(exp(mcmc$mu.chain+qnorm(0.95)*mcmc$sigma.chain),0.95),
                F_est = unname(quantile(100*(1 - pnorm( log(oel) , mean = mcmc$mu.chain , sd = mcmc$sigma.chain )),0.5)),
                F_70ucl = unname(quantile(100*(1 - pnorm( log(oel) , mean = mcmc$mu.chain , sd = mcmc$sigma.chain )),0.7)),
                F_95ucl = unname(quantile(100*(1 - pnorm( log(oel) , mean = mcmc$mu.chain , sd = mcmc$sigma.chain )),0.95)))
  
  
  return(results)
  
}


#' Measurement error bayesian analysis using the expostats prior
#'
#' @param mysample sample to be analysed
#' @param oel occupational exposure limit
#' @param me_cv coefficient of variation in proportion
#'
#' @return point estimates for GM, GSD, exceedance and P95, and 70 and 95 UCL for P95 and exceedance
#'

expostats.me <- function( mysample , oel , me_cv) {
  
  mcmc <- suppressWarnings( Webexpo.seg.globalbayesian.jags( data.sample = mysample ,
                                           is.lognormal = TRUE , 
                                           error.type = "CV" ,
                                           me.range = c(me_cv,me_cv) , 
                                           oel = oel ,
                                           prior.model = "informedvar",
                                           n.iter = 25000) )
  
  
  results <- c( gm_est = exp(median(mcmc$mu.chain)),
               gsd_est = exp(median(mcmc$sigma.chain)),
               p95_est = median(exp(mcmc$mu.chain+qnorm(0.95)*mcmc$sigma.chain)),
               p95_70ucl = quantile(exp(mcmc$mu.chain+qnorm(0.95)*mcmc$sigma.chain),0.7),
               p95_95ucl = quantile(exp(mcmc$mu.chain+qnorm(0.95)*mcmc$sigma.chain),0.95),
               F_est = unname(quantile(100*(1 - pnorm( log(oel) , mean = mcmc$mu.chain , sd = mcmc$sigma.chain )),0.5)),
               F_70ucl = unname(quantile(100*(1 - pnorm( log(oel) , mean = mcmc$mu.chain , sd = mcmc$sigma.chain )),0.7)),
               F_95ucl = unname(quantile(100*(1 - pnorm( log(oel) , mean = mcmc$mu.chain , sd = mcmc$sigma.chain )),0.95)))
  

  return(results)
  
}


#' Measurement error bayesian analysis using the expostats prior
#'
#' @param mysample sample to be analysed
#' @param oel occupational exposure limit
#' @param me_cv coefficient of variation in proportion
#'
#' @return point estimates for GM, GSD, exceedance and P95, and 70 and 95 UCL for P95 and exceedance
#'

expostats.me.s <- function( mysample , oel , me_cv , models.list) {
  
  mcmc <- suppressWarnings( Webexpo.seg.globalbayesian.stan( data.sample = mysample ,
                                                             is.lognormal = TRUE , 
                                                             error.type = "CV" ,
                                                             me.range = c(me_cv-0.001,me_cv+0.001) , 
                                                             oel = oel ,
                                                             prior.model = "informedvar",
                                                             n.iter = 25000,
                                                             models.list=models.list) )
  
  
  results <- c( gm_est = exp(median(mcmc$mu.chain)),
                gsd_est = exp(median(mcmc$sigma.chain)),
                p95_est = median(exp(mcmc$mu.chain+qnorm(0.95)*mcmc$sigma.chain)),
                p95_70ucl = quantile(exp(mcmc$mu.chain+qnorm(0.95)*mcmc$sigma.chain),0.7),
                p95_95ucl = quantile(exp(mcmc$mu.chain+qnorm(0.95)*mcmc$sigma.chain),0.95),
                F_est = unname(quantile(100*(1 - pnorm( log(oel) , mean = mcmc$mu.chain , sd = mcmc$sigma.chain )),0.5)),
                F_70ucl = unname(quantile(100*(1 - pnorm( log(oel) , mean = mcmc$mu.chain , sd = mcmc$sigma.chain )),0.7)),
                F_95ucl = unname(quantile(100*(1 - pnorm( log(oel) , mean = mcmc$mu.chain , sd = mcmc$sigma.chain )),0.95)))
  
  
  return(results)
  
}




#' Naive frequentist analysis 
#'
#' @param mysample sample to be analysed
#' @param oel occupational exposure limit
#'
#' @return point estimates for GM, GSD, exceedance and P95, and 70 and 95 UCL for P95 and exceedance
#'

frequentist.naive <- function( mysample , oel  ) {
  
  # ROS as implemented in ND expo (if less than 3 detects or less than 5 total : LOQ/2)
  
  mysample.ros <- fun.NdExpo.lognorm( mysample )$data$xfin
  
  # parameter estimates
  
  # if variability is not estimable :  NA if there is only one unique value (e.g. all NDs imputed as LD/2)
  
  if ( length(unique(mysample.ros)) == 1 ) return( c( gm_est = NA , gsd_est = NA , p95_est = NA , p95_70ucl = NA , p95_95ucl = NA , F_est = NA , F_70ucl = NA , F_95ucl = NA ))
  
  else {    
    
    gm_est <- exp(mean(log(mysample.ros)))
    
    gsd_est <- exp(sd(log(mysample.ros)))
    
    p95_est <- fun.perc.en689(mysample.ros,alpha=0.05,perc=0.95)$est
    
    p95_70ucl <- fun.perc.en689(mysample.ros,alpha=0.30,perc=0.95)$uc
    
    p95_95ucl <- fun.perc.en689(mysample.ros,alpha=0.05,perc=0.95)$uc
    
    F_est <- fun.frac.dep( mysample.ros , gam = 0.95, L = oel , logx = TRUE, wpnt = FALSE)$fe 
    
    F_70ucl <- fun.frac.dep( mysample.ros , gam = 0.70, L = oel , logx = TRUE, wpnt = FALSE)$fe.UCL
    
    F_95ucl <- fun.frac.dep( mysample.ros , gam = 0.95, L = oel , logx = TRUE, wpnt = FALSE)$fe.UCL
    
    results <- c( gm_est = gm_est,
                  gsd_est = gsd_est,
                  p95_est = p95_est,
                  p95_70ucl = p95_70ucl,
                  p95_95ucl = p95_95ucl,
                  F_est = F_est,
                  F_70ucl = F_70ucl,
                  F_95ucl = F_95ucl)
    
    
    
    return(results) }
  
}



#' Measurement error frequentist analysis using GUM approach as implemented by Scheffers and Emond 
#'
#' @param mysample sample to be analysed
#' @param oel occupational exposure limit
#' @param me_cv coefficient of variation in proportion
#' @param n_iterations number of iteration for the GUM approach
#' 
#' @return mean, median and 4 quantiles summarizing point estimates for GM, GSD, exceedance and P95, and 70 and 95 UCL for P95 and exceedance across iterations, as 

frequentist.me <- function( mysample , oel , me_cv , n_iterations_gum = 10000  ) {
  
  # preliminary calculations
  
  sample_size <- length(mysample)
  
  
  # ROS as implemented in ND expo (if less than 3 detects or less than 5 total : LOQ/2)
  
  mysample.ros <- fun.NdExpo.lognorm( mysample )$data$xfin
  
  # if variability is not estimable :  NA if there is only one unique value (e.g. all NDs imputed as LD/2)
  
  if ( length(unique(mysample.ros)) == 1 ) {  
    
    results <- list( mean = c( gm_est = NA , gsd_est = NA , p95_est = NA , p95_70ucl = NA , p95_95ucl = NA , F_est = NA , F_70ucl = NA , F_95ucl = NA ),
                     median = c( gm_est = NA , gsd_est = NA , p95_est = NA , p95_70ucl = NA , p95_95ucl = NA , F_est = NA , F_70ucl = NA , F_95ucl = NA ),
                     q2.5 = c( gm_est = NA , gsd_est = NA , p95_est = NA , p95_70ucl = NA , p95_95ucl = NA , F_est = NA , F_70ucl = NA , F_95ucl = NA ),
                     q5 = c( gm_est = NA , gsd_est = NA , p95_est = NA , p95_70ucl = NA , p95_95ucl = NA , F_est = NA , F_70ucl = NA , F_95ucl = NA ),
                     q95 = c( gm_est = NA , gsd_est = NA , p95_est = NA , p95_70ucl = NA , p95_95ucl = NA , F_est = NA , F_70ucl = NA , F_95ucl = NA ),
                     q97.5 = c( gm_est = NA , gsd_est = NA , p95_est = NA , p95_70ucl = NA , p95_95ucl = NA , F_est = NA , F_70ucl = NA , F_95ucl = NA ))    
    
    return(results) }
  
  else { 
    
    
    # GUM simulation approach (replicating the sample, then adding random noise to each replicate)
    
    data_matrix <- as.matrix( replicate( n_iterations_gum , mysample.ros ))
    
    data_matrix_me <- apply( data_matrix, 2 , function(x) { pmax( rep(0.001,sample_size)  , rnorm( sample_size , x , x*me_cv ) )  }  )
    
    
    # parameter estimates across iterations
    
    result_matrix <- apply( data_matrix_me , 2 , function(x) { 
      
      
      gm_est <- exp(mean(log(x)))
      
      gsd_est <- exp(sd(log(x)))
      
      p95_est <- fun.perc.en689(x,alpha=0.05,perc=0.95)$est
      
      p95_70ucl <- fun.perc.en689(x,alpha=0.30,perc=0.95)$uc
      
      p95_95ucl <- fun.perc.en689(x,alpha=0.05,perc=0.95)$uc
      
      F_est <- fun.frac.dep( x , gam = 0.95, L = oel , logx = TRUE, wpnt = FALSE)$fe 
      
      F_70ucl <- fun.frac.dep( x , gam = 0.70, L = oel , logx = TRUE, wpnt = FALSE)$fe.UCL
      
      F_95ucl <- fun.frac.dep( x , gam = 0.95, L = oel , logx = TRUE, wpnt = FALSE)$fe.UCL
      
      return( c( gm_est , gsd_est , p95_est , p95_70ucl , p95_95ucl , F_est , F_70ucl , F_95ucl ) )
      
    } )
    
    
    # summaries across iterations
    
    results <- list( mean = apply( result_matrix , 1 , mean ),
                     median = apply( result_matrix , 1 , median ),
                     q2.5 = apply( result_matrix , 1 , function(x) { quantile(x,0.025) } ),
                     q5 = apply( result_matrix , 1 , function(x) { quantile(x,0.05) } ),
                     q95 = apply( result_matrix , 1 , function(x) { quantile(x,0.95) } ),
                     q97.5 = apply( result_matrix , 1 , function(x) { quantile(x,0.975) } ))    
    
    names(results$mean) <- c("gm_est","gsd_est","p95_est","p95_70ucl","p95_95ucl","F_est","F_70ucl","F_95ucl")
    names(results$median) <- c("gm_est","gsd_est","p95_est","p95_70ucl","p95_95ucl","F_est","F_70ucl","F_95ucl")
    names(results$q2.5) <- c("gm_est","gsd_est","p95_est","p95_70ucl","p95_95ucl","F_est","F_70ucl","F_95ucl")
    names(results$q5) <- c("gm_est","gsd_est","p95_est","p95_70ucl","p95_95ucl","F_est","F_70ucl","F_95ucl")
    names(results$q95) <- c("gm_est","gsd_est","p95_est","p95_70ucl","p95_95ucl","F_est","F_70ucl","F_95ucl")
    names(results$q97.5) <- c("gm_est","gsd_est","p95_est","p95_70ucl","p95_95ucl","F_est","F_70ucl","F_95ucl")
    
    return(results) }
  
}






