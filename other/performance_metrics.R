#'   functions for calculating performance metrics from a vector of estimates and the true value of the parameter

#' Root mean squared error (RMSE)
#' https://en.wikipedia.org/wiki/Root_mean_square_deviation
#'
#' @param x : vector of parameter estimates across simulations
#' @param theta : true value of the parameter 

#' @return RMSE
#'


compute.rmse = function(x, theta){
  
  mean.x = mean(x)
  
  N = length(x)
  
  if(is.null(N))
    stop("x must be a vector of minimum length 1")
  
  rmse = sqrt((mean.x - theta) ^ 2 + sum((x - mean.x) ^ 2) / (N - 1))
  return(rmse)
}


#' Normalized, or relative, Root mean squared error (NRMSE)
#' https://en.wikipedia.org/wiki/Root_mean_square_deviation
#'
#' @param x : vector of parameter estimates across simulations
#' @param theta : true value of the parameter 

#' @return normalized RMSE
#'


compute.relrmse = function(x, theta){
  
  mean.x = mean(x)
  
  N = length(x)
  
  if(is.null(N))
    stop("x must be a vector of minimum length 1")
  
  rmse = 100*sqrt((mean.x - theta) ^ 2 + sum((x - mean.x) ^ 2) / (N - 1))/theta
  return(rmse)
}


#' Bias
#' https://en.wikipedia.org/wiki/Bias_of_an_estimator
#'
#' @param x : vector of parameter estimates across simulations
#' @param theta : true value of the parameter 

#' @return bias
#'


compute.bias = function(x, theta){
  
  bias = mean(x) - theta 
  
  return(bias)
}


#' Normalised , or relative Bias
#' https://en.wikipedia.org/wiki/Bias_of_an_estimator
#'
#' @param x : vector of parameter estimates across simulations
#' @param theta : true value of the parameter 

#' @return bias
#'


compute.normalized.bias = function(x, theta){
  
  nbias = 100 * (mean(x) - theta) / theta
  
  return(nbias)
}


#' Precision
#' https://en.wikipedia.org/wiki/Accuracy_and_precision
#'
#' @param x : vector of parameter estimates across simulations
#' @param theta : true value of the parameter 

#' @return bias
#'

compute.precision = function(x, theta){
  
  error = x-theta
  
  p = sd(error)
  
  return(p)
}


#' Normalized or relative Precision
#' https://en.wikipedia.org/wiki/Accuracy_and_precision
#'
#' @param x : vector of parameter estimates across simulations
#' @param theta : true value of the parameter 

#' @return bias
#'


compute.normalized.precision = function(x, theta){
  
  error = x-theta
  
  np = 100*sd(error)/theta
  
  return(np)
}


#' Median error
#' 
#'
#' @param x : vector of parameter estimates across simulations
#' @param theta : true value of the parameter 

#' @return Median error
#'

compute.median.error = function(x, theta){
  
  me = median( x - theta )
  
  return(me)
}


#' Normalized or relative Median error
#' 
#'
#' @param x : vector of parameter estimates across simulations
#' @param theta : true value of the parameter 

#' @return Median error
#'

compute.normalized.median_error = function(x, theta){
  
  me = median( x - theta )

  nme <- 100*me/theta
    
  return(nme)
}



#' function for the median absolute deviation (mad), a robust measure of precision
#' https://en.wikipedia.org/wiki/Median_absolute_deviation
#'
#' @param x : vector of parameter estimates across simulations
#' @param theta : true value of the parameter 

#' @return mad
#'


compute.mad = function(x, theta){
  
  mad = median(abs(x - median(x)))
  
  return(mad)
}



#' function for the rmsle, a robust measure of accuracy
#' https://medium.com/analytics-vidhya/root-mean-square-log-error-rmse-vs-rmlse-935c6cc1802a
#'
#' @param x : vector of parameter estimates across simulations
#' @param theta : true value of the parameter 

#' @return rmsle
#'

compute.rmsle = function(x, theta){
  
  rmsle = sqrt(mean((log(x + 1) - log(theta + 1))^2))
  
  return(rmsle)
}
