#' Function for estimating ta percentile from a vector of lognormal exposure values
#' 
#' Based on Zwillinger, D.; Kokoska, S. (2000) Standard Probability and Statistics Tables and Formulae. Chapman & Hall / CRC, Boca Raton, FL
#' Pp 176
#' 
#' 
#' @param x vector of values to be analysed ( at least 2 values)
#' @param alpha alpha value associated with the probability for the confidence limit ( alpha=0.05 corresponds to a 95% upper confidence limit
#' @param perc percentile of interest (in probability).
#' @param logx if TRUE function assume lognormal distribution, normal if FALSE
#' @param wpnt internal parameter, do not use.
#' 
#' @return est : percentile demande
#' @return uc : limite de confiance superieure
#' @return k : parametre interne de calcul



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
