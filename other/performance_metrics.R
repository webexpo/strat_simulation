#'   functions for calculating performance metrics from a vector of estimates and the true value of the parameter

#' Root mean squared error (RMSE)
#' https://en.wikipedia.org/wiki/Root_mean_square_deviation
#'
#' @param x : vector of parameter estimates across simulations
#' @param theta : true value of the parameter 

#' @return RMSE, NA values in x are excluded from the calculation, if all values are NA, return 0
#'


compute.rmse = function(x, theta){
  
  ## management of non estimable values
  
  n_na <- sum(is.na(x))
  
  # if all values are NA, return 0
  
  if (n_na == length(x)) return(0) else {
    
    # if not all values are NA, remove NA values
    
    x <- x[!is.na(x)]
    
    # calculations
    
    mean.x = mean(x)
    
    N = length(x)
    
    rmse = sqrt((mean.x - theta) ^ 2 + sum((x - mean.x) ^ 2) / (N - 1))
    
    return(rmse)
    
  }
}


#' Normalized, or relative, Root mean squared error (NRMSE)
#' https://en.wikipedia.org/wiki/Root_mean_square_deviation
#'
#' @param x : vector of parameter estimates across simulations
#' @param theta : true value of the parameter 

#' @return normalized RMSE, NA values in x are excluded from the calculation, if all values are NA, return 0
#'


compute.relrmse = function(x, theta){
  
  ## management of non estimable values
  
  n_na <- sum(is.na(x))
  
  # if all values are NA, return 0
  
  if (n_na == length(x)) return(0) else {
    
    # if not all values are NA, remove NA values
    
    x <- x[!is.na(x)]
    
    
    # calculations
    
    mean.x = mean(x)
    
    N = length(x)
    
    rmse = 100*sqrt((mean.x - theta) ^ 2 + sum((x - mean.x) ^ 2) / (N - 1))/theta
    
    return(rmse)
    
  }
}


#' Bias
#' https://en.wikipedia.org/wiki/Bias_of_an_estimator
#'
#' @param x : vector of parameter estimates across simulations
#' @param theta : true value of the parameter 

#' @return bias, NA values in x are excluded from the calculation, if all values are NA, return 0
#'


compute.bias = function(x, theta){
  
  ## management of non estimable values
  
  n_na <- sum(is.na(x))
  
  # if all values are NA, return 0
  
  if (n_na == length(x)) return(0) else {
    
    # if not all values are NA, remove NA values
    
    x <- x[!is.na(x)]
    
    # calculations
    
    bias = mean(x) - theta 
    
    return(bias)
    
  }
}


#' Normalised , or relative Bias
#' https://en.wikipedia.org/wiki/Bias_of_an_estimator
#'
#' @param x : vector of parameter estimates across simulations
#' @param theta : true value of the parameter 

#' @return normalized bias, NA values in x are excluded from the calculation, if all values are NA, return 0
#'


compute.normalized.bias = function(x, theta){
  
  ## management of non estimable values
  
  n_na <- sum(is.na(x))
  
  # if all values are NA, return 0
  
  if (n_na == length(x)) return(0) else {
    
    # if not all values are NA, remove NA values
    
    x <- x[!is.na(x)]
    
    # calculations
    
    nbias = 100 * (mean(x) - theta) / theta
    
    return(nbias)
    
  }
}


#' Precision
#' https://en.wikipedia.org/wiki/Accuracy_and_precision
#'
#' @param x : vector of parameter estimates across simulations
#' @param theta : true value of the parameter 

#' @return precision, NA values in x are excluded from the calculation, if all values are NA, return 0
#'

compute.precision = function(x, theta){
  
  ## management of non estimable values
  
  n_na <- sum(is.na(x))
  
  # if all values are NA, return 0
  
  if (n_na == length(x)) return(0) else {
    
    # if not all values are NA, remove NA values
    
    x <- x[!is.na(x)]
    
    # calculations
    
    error = x-theta
    
    p = sd(error)
    
    return(p)
    
  }
}


#' Normalized or relative Precision
#' https://en.wikipedia.org/wiki/Accuracy_and_precision
#'
#' @param x : vector of parameter estimates across simulations
#' @param theta : true value of the parameter 

#' @return normalized precision, NA values in x are excluded from the calculation, if all values are NA, return 0
#'


compute.normalized.precision = function(x, theta){
  
  ## management of non estimable values
  
  n_na <- sum(is.na(x))
  
  # if all values are NA, return 0
  
  if (n_na == length(x)) return(0) else {
    
    # if not all values are NA, remove NA values
    
    x <- x[!is.na(x)]
    
    # calculations
    
    error = x-theta
    
    np = 100*sd(error)/theta
    
    return(np)
    
  }
}


#' Median error
#' 
#'
#' @param x : vector of parameter estimates across simulations
#' @param theta : true value of the parameter 

#' @return Median error, NA values in x are excluded from the calculation, if all values are NA, return 0
#'

compute.median.error = function(x, theta){
  
  ## management of non estimable values
  
  n_na <- sum(is.na(x))
  
  # if all values are NA, return 0
  
  if (n_na == length(x)) return(0) else {
    
    # if not all values are NA, remove NA values
    
    x <- x[!is.na(x)]
    
    # calculations
    
    me = median( x - theta )
    
    return(me)
    
  }
}


#' Normalized or relative Median error
#' 
#'
#' @param x : vector of parameter estimates across simulations
#' @param theta : true value of the parameter 

#' @return Normalized median error, NA values in x are excluded from the calculation, if all values are NA, return 0
#'

compute.normalized.median_error = function(x, theta){
  
  ## management of non estimable values
  
  n_na <- sum(is.na(x))
  
  # if all values are NA, return 0
  
  if (n_na == length(x)) return(0) else {
    
    # if not all values are NA, remove NA values
    
    x <- x[!is.na(x)]
    
    # calculations
    
    me = median( x - theta )
    
    nme <- 100*me/theta
    
    return(nme)
    
  }
}



#' function for the median absolute deviation (mad), a robust measure of precision
#' https://en.wikipedia.org/wiki/Median_absolute_deviation
#'
#' @param x : vector of parameter estimates across simulations
#' @param theta : true value of the parameter 

#' @return mad, NA values in x are excluded from the calculation, if all values are NA, return 0
#'


compute.mad = function(x, theta){
  
  ## management of non estimable values
  
  n_na <- sum(is.na(x))
  
  # if all values are NA, return 0
  
  if (n_na == length(x)) return(0) else {
    
    # if not all values are NA, remove NA values
    
    x <- x[!is.na(x)]
    
    # calculations
    
    mad = median(abs(x - median(x)))
    
    return(mad)
    
  }
}



#' function for the rmsle, a robust measure of accuracy
#' https://medium.com/analytics-vidhya/root-mean-square-log-error-rmse-vs-rmlse-935c6cc1802a
#'
#' @param x : vector of parameter estimates across simulations
#' @param theta : true value of the parameter 

#' @return rmsle, NA values in x are excluded from the calculation, if all values are NA, return 0
#'

compute.rmsle = function(x, theta){
  
  ## management of non estimable values
  
  n_na <- sum(is.na(x))
  
  # if all values are NA, return 0
  
  if (n_na == length(x)) return(0) else {
    
    # if not all values are NA, remove NA values
    
    x <- x[!is.na(x)]
    
    # calculations
    
    rmsle = sqrt(mean((log(x + 1) - log(theta + 1))^2))
    
    return(rmsle)
    
  }
}

