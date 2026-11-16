#' Function for estimating the arithmetic mean from a vector of lognormal exposure values
#'
#' Point estimate: minimum variance unbiased estimator (MVUE), Finney's series g_n(s^2/2)
#' Confidence limits: Land's "exact" method, with the C factors computed from the
#' Hewett & Ganser (1997) approximation equations, as implemented in IHSTAT macro-free Classic v1.01 (AIHA)
#'
#' References
#' - Finney DJ (1941) On the distribution of a variate whose logarithm is normally distributed.
#'   J R Stat Soc Suppl 7:155-161
#' - Land CE (1975) Tables of confidence limits for linear functions of the normal mean and variance.
#'   In: Selected Tables in Mathematical Statistics, Vol. III. American Mathematical Society
#' - Hewett P, Ganser GH (1997) Simple procedures for calculating confidence intervals around the
#'   sample mean and exceedance fraction derived from lognormally distributed data.
#'   Appl Occup Environ Hyg 12(2):132-142
#'   MVUE: Equations 2-3; Land's limits: Equations 4-5; C factors: Appendix A
#'   (valid for 3 <= n <= 1001 and 0.01 <= s(ln x) <= 4; vs. Land's published C factors,
#'   error -3 to +2% for the LCL C factor and -3 to +5% for the UCL C factor)
#'
#' Notes
#' - The C approximation equations exist ONLY for one-sided 95% limits (LCL95 / UCL95);
#'   unlike fun.frac.dep / fun.perc there is no confidence-level argument.
#' - Valid for n >= 3. As in IHSTAT, s(ln x) is floored at 0.01 and capped at 4 when computing C.
#' - IHSTAT truncates the MVUE series at 5 terms (nterms = 5); default here sums to convergence.
#'
#' @param x vector of values to be analysed (at least 3 values)
#' @param L exposure limit (optional; if provided, limits are also returned as % of L)
#' @param nterms number of series terms for the MVUE (NULL = until convergence; 5 = IHSTAT)
#'
#' @return mvue : MVUE of the arithmetic mean
#' @return mle : MLE-type estimate exp(mean + s^2/2) (s with n-1), used as anchor for Land's limits
#' @return lcl : Land's 95% one-sided lower confidence limit
#' @return ucl : Land's 95% one-sided upper confidence limit
#' @return C.lcl, C.ucl : Land's C factors
#' @return mvue.pct, ucl.pct : MVUE and UCL in % of L (if L provided)

fun.am <- function(x, L = NA, nterms = NULL) {

  if (any(x <= 0)) stop("all data values must be positive")

  n <- length(x)
  if (n < 3) stop("at least 3 values required (C approximation equations)")

  muy  <- mean(log(x))
  sigy <- sd(log(x))

  ## MVUE : exp(muy) * g_n(s^2/2), Finney's series
  t <- sigy^2 / 2
  g <- 1
  k <- 1
  repeat {
    term <- (n - 1)^(2 * k - 1) * t^k /
      (factorial(k) * n^k * prod(n + 2 * seq_len(k - 1) - 1))
    g <- g + term
    if ((!is.null(nterms) && k >= nterms) ||
        (is.null(nterms) && abs(term) < 1e-12 * g) || k >= 100) break
    k <- k + 1
  }
  mvue <- exp(muy) * g

  ## Land's limits : exp( ln(MLE) + C * s / sqrt(n-1) )
  mle <- exp(muy + sigy^2 / 2)

  s1 <- max(sigy, 0.01)   # used in the final formula
  s2 <- min(s1, 4)        # used for C
  m  <- n - 2

  # C for UCL95
  a <- 0.76766658; b <- 3.8716869; c <- 0.80598919; d <- 6.0321019
  e <- 0.89998154; f <- 2.012669;  g2 <- 0.21978875; h <- 0.41575588; i <- 0.29258276
  F1 <- s2 * (i + 1 / m^c)
  F2 <- b + d / m^c
  F3 <- F1 * (1 - e * exp(-f * F1))
  F4 <- 1 + g2 * exp(-h * F1)
  C.ucl <- 1.645 + a / m + F2 * F3 / F4

  # C for LCL95
  a <- -0.85033767; b <- -0.5258052; c <- 0.92416176; d <- 3.3298209
  e <- 0.94348568;  f <- 1.3213281;  g2 <- 0.8155562; h <- -1.018148; i <- 0.25248895
  F1 <- a / m^c + b
  F2 <- d / m^f + e
  F3 <- g2 / m^i + h
  F4 <- -0.6226 / m^0.2426 + 1.147
  C.lcl <- -0.74295 / m - 1.64765 + s2 * (F3 + F2 * exp(F1 * s2^F4))

  lcl <- exp(log(mle) + C.lcl * s1 / sqrt(n - 1))
  ucl <- exp(log(mle) + C.ucl * s1 / sqrt(n - 1))

  out <- list(mvue = mvue, mle = mle, lcl = lcl, ucl = ucl,
              C.lcl = C.lcl, C.ucl = C.ucl)
  if (!is.na(L)) {
    out$mvue.pct <- 100 * mvue / L
    out$ucl.pct  <- 100 * ucl / L
  }
  out
}
