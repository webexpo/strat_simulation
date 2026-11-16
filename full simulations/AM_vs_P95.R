#################################################################################
#
#
#   Comparing two compliance criteria for lognormal exposures :
#
#     95% UCL on the arithmetic mean  vs.  70% UCL on the 95th percentile (EN 689)
#
#   Ratio UCL95(AM) / UCL70(P95) by sample size (n) and GSD
#     - frequentist : Land / Hewett & Ganser 1997 (AM), exact K tolerance factor (P95)
#     - Bayesian    : Expostats informedvar model (JAGS)
#
#   Ratio > 1 : the 95% UCL on the AM is more conservative
#   Ratio < 1 : the 70% UCL on the P95 is more conservative
#
#   Each n x GSD cell uses one "ideal" sample : normal quantiles, rescaled so that
#   the sample GM = 1 and the sample GSD = the target GSD exactly.
#   The frequentist ratio depends only on n and the sample GSD (not on the GM).
#
##################################################################################


# libraries ---------------------------------------------

    library(tolerance)

    ##  frequentist functions
    source("parameter estimation/frequentist/arithmetic_mean.R")
    source("parameter estimation/frequentist/percentile.R")

    ##  Webexpo functions
    source("parameter estimation/bayesian/load.webexpo.SEG.functions.R")


# scenarios -----------------------------------------------

    n.values   <- c(3, 4, 6, 10, 20, 50)

    gsd.values <- c(1.5, 2, 2.5, 3, 3.5, 4)

    oel <- 1   # arbitrary : neither the AM nor the P95 limits depend on the OEL


# functions -----------------------------------------------

    # ideal sample : sample GM = 1 and sample GSD = gsd exactly

    fun.ideal.sample <- function(n, gsd) {

      z <- qnorm(ppoints(n))

      z <- (z - mean(z)) / sd(z)

      exp(log(gsd) * z)

    }


    # frequentist ratio : Land 95% UCL on AM / 70% UCL on P95 (exact K)

    fun.ratio.freq <- function(x) {

      ucl.am  <- fun.am(x)$ucl

      ucl.p95 <- fun.perc.en689(x, alpha = 0.30, perc = 0.95)$uc

      c(ucl.am = ucl.am, ucl.p95 = ucl.p95, ratio = ucl.am / ucl.p95)

    }


    # Bayesian ratio : Expostats informedvar (JAGS), one-sided upper limits
    #   probacred = 90 -> 95% one-sided upper limit ; probacred = 40 -> 70% one-sided upper limit

    fun.ratio.bayes <- function(x, oel) {

      mcmc <- Webexpo.seg.globalbayesian.jags(data.sample = as.character(x),
                                              is.lognormal = TRUE,
                                              error.type = "none",
                                              oel = oel)

      ucl.am <- all.numeric(mu.chain = mcmc$mu.chain,
                            sigma.chain = mcmc$sigma.chain,
                            probacred = 90,
                            oel = oel,
                            frac_threshold = 5,
                            target_perc = 95)$am$ucl

      ucl.p95 <- all.numeric(mu.chain = mcmc$mu.chain,
                             sigma.chain = mcmc$sigma.chain,
                             probacred = 40,
                             oel = oel,
                             frac_threshold = 5,
                             target_perc = 95)$perc$ucl

      c(ucl.am = unname(ucl.am), ucl.p95 = unname(ucl.p95), ratio = unname(ucl.am / ucl.p95))

    }


# frequentist table -----------------------------------------------

    res.freq <- expand.grid(n = n.values, gsd = gsd.values)

    res.freq <- cbind(res.freq,
                      t(mapply(function(n, gsd) fun.ratio.freq(fun.ideal.sample(n, gsd)),
                               res.freq$n, res.freq$gsd)))

    table.freq <- round(xtabs(ratio ~ n + gsd, res.freq), 2)

    table.freq


# Bayesian table -----------------------------------------------

    set.seed(20261116)

    res.bayes <- expand.grid(n = n.values, gsd = gsd.values)

    res.bayes <- cbind(res.bayes,
                       t(mapply(function(n, gsd) fun.ratio.bayes(fun.ideal.sample(n, gsd), oel = oel),
                                res.bayes$n, res.bayes$gsd)))

    table.bayes <- round(xtabs(ratio ~ n + gsd, res.bayes), 2)

    table.bayes
