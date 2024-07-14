#' Function for estimating ta percentile from a vector of lognormal exposure values
#' 
#' estimation methods for the 95th percentile tolerance limit as presented in [Selvin et al](https://www.tandfonline.com/doi/abs/10.1080/15298668791384445)
#' The k values to calculate a (100-100alpha)% upper confidence limit on the 95th percentile (UTL95,100-100alpha) was determined based on Zwillinger, D.; Kokoska, S. (2000) Standard Probability and Statistics Tables and Formulae. Chapman & Hall / CRC, Boca Raton, FL
#' Pp 176
#' 
#' 
#' @param x vector of values to be analysed ( at least 2 values)
#' @param alpha alpha value associated with the probability for the confidence limit ( alpha=0.05 corresponds to a 95% upper confidence limit)
#' @param perc percentile of interest (as a proportion not %).
#' @param logx if TRUE function assume lognormal distribution, normal if FALSE

#' @return est : percentile value
#' @return uc : upper confidence limit
#' @return k : k factor for the calculation of the upper confidence limit (see Selvin et al)


fun.perc <-function(x,alpha=0.05,perc=0.95) {
  
  muy <-mean(log(x))
  
  sigy <-sd(log(x))
  
  n<-length(x)
  
  perc.est <-exp(muy+sigy*qnorm(perc))
  
  u<-1-(qnorm(1-alpha)^2)/(2*(n-1))
  
  v <-qnorm(perc)^2-(qnorm(1-alpha)^2)/n
  
  k <-(qnorm(perc)+sqrt(qnorm(perc)^2-u*v))/(u)
  
  perc.uc <-exp(muy+k*sigy)
  
  perc.lc <-exp(muy-k*sigy)
  
  
  return(list(est=perc.est,uc=perc.uc,lc=perc.lc,k=k))
  
}


#' Function for estimating ta percentile from a vector of lognormal exposure values
#' 
#' estimation methods for the 95th percentile tolerance limit as presented in [Selvin et al](https://www.tandfonline.com/doi/abs/10.1080/15298668791384445)
#' The k values to calculate a (100-100alpha)% upper confidence limit on the 95th percentile (UTL95,100-100alpha) was determined with the R package named ["tolerance"](https://cran.r-project.org/web/packages/tolerance/index.html)
#' 
#' 
#' @param x vector of values to be analysed ( at least 2 values)
#' @param alpha alpha value associated with the probability for the confidence limit ( alpha=0.05 corresponds to a 95% upper confidence limit)
#' @param perc percentile of interest (as a proportion not %).
#' @param logx if TRUE function assume lognormal distribution, normal if FALSE

#' @return est : percentile value
#' @return uc : upper confidence limit
#' @return k : k factor for the calculation of the upper confidence limit (see Selvin et al)


fun.perc.en689 <-function(x,alpha=0.05,perc=0.95) {

  muy <-mean(log(x))
  
  sigy <-sd(log(x))
  
  n<-length(x)
  
  perc.est <-exp(muy+sigy*qnorm(perc))
  
  k <- K.factor(n, alpha = alpha, P = perc, side = 1,
                method = c("EXACT"), m = 100)
  
  perc.uc <-exp(muy+k*sigy)
  
  perc.lc <-exp(muy-k*sigy)
  
  
  return(list(est=perc.est,uc=perc.uc,lc=perc.lc,k=k))
  
}

