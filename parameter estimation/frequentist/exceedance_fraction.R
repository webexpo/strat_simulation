#' Function for estimating the exceedance fraction from a vector of lognormal exposure values
#' 
#' Based on Wild et al. (1996) Environmetrics 7(3) 247-259
#' 
#' 
#' @param x vector of values to be analysed ( at least 2 values)
#' @param gam probability for one sided confidence limits (default 0.95) : LCL and UCL form an IC with probability [100-200*(1-gam)%]
#' @param L  exposure limit
#' @param logx if TRUE function assume lognormal distribution, normal if FALSE
#' @param wpnt internal parameter, do not use.
#' 
#' @return fe : exceedance fraction in %
#' @return fe.LCL : lower confidence limit
#' @return fe.UCL : upper confidence limit
#' @return L : exposure limit
#' @return Logx : initial choice of distribution

fun.frac.dep <-function (x, gam = 0.95, L = 100 , logx = TRUE, wpnt = FALSE) 
{
  if (!wpnt) 
    options(warn = -1)
  if (is.na(L) || L <= 0) 
    stop("Value of L is missing or <= 0")
  n <- length(x)
  if (logx) {
    if (any(x <= 0)) 
      stop("all data values must be positive")
    y <- log(x)
    LL <- log(L)
  }
  else {
    y <- x
    LL <- L
  }
  yb <- mean(y)
  sd <- sd(y)
  del <- function(ncp, tv = t0, df = n - 1, eps = cv) {
    pt(tv, df, ncp) - eps
  }
  u <- (LL - yb)/sd
  t0 <- sqrt(n) * u
  cv <- gam
  dap <- t0 - qnorm(gam) * (1 + t0^2/(2 * (n - 1)))^0.5
  cd <- dap * c(-2, 2)
  while (del(cd[1], t0, n - 1, cv) * del(cd[2], t0, n - 1, 
                                         cv) >= 0) {
    cd <- cd * 2
  }
  u2 <- uniroot(del, cd)$root
  cv <- 1 - gam
  dap <- t0 - qnorm(cv) * (1 + t0^2/(2 * (n - 1)))^0.5
  cd <- dap * c(-2, 2)
  while (del(cd[1], t0, n - 1, cv) * del(cd[2], t0, n - 1, 
                                         cv) >= 0) {
    cd <- cd * 2
  }
  u1 <- uniroot(del, cd)$root
  out <- c(u1, u2)/sqrt(n)
  out <- 100 * c(1 - pnorm(u), 1 - pnorm(out))
  if (!wpnt) 
    options(warn = 0)
  out <- list(fe = out[1], fe.LCL = out[2], fe.UCL = out[3], 
              L = L, gam = gam, Logx = logx)
  out
}
