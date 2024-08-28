#' loads the Webexpo functions necessary to run the Expostats Bayesian models relying on JAGS
#' @requirement needs JAGS to be installed on the computer and the rjags library

devtools::source_url("https://github.com/lhimp/scripts/raw/master/chemin.R")

library(rjags)

# Webexpo models  -----------------------------------------------------------------

## data preparation

chemin(
  fileName = "webexpo.seg.dataprep.R",
  relPath = c("DATA PREPARATION", "SEG ANALYSIS"),
  githubCredentials = list(userName = "webexpo", repoName = "webexpo_r_lib")
)



# JAGS models EXPOSTATS

chemin(
  fileName = "webexpo.seg.mainbayesian.R",
  relPath = c("jags models", "SEG ANALYSIS"),
  githubCredentials = list(userName = "webexpo", repoName = "webexpo_r_lib")
)

chemin(
  fileName = "webexpo.seg.informedvarbayesian.R",
  relPath = c("jags models", "SEG ANALYSIS"),
  githubCredentials = list(userName = "webexpo", repoName = "webexpo_r_lib")
)

chemin(
  fileName = "webexpo.seg.informedvarbayesian.models.R",
  relPath = c("jags models", "SEG ANALYSIS"),
  githubCredentials = list(userName = "webexpo", repoName = "webexpo_r_lib")
)

# WEBEXPO R models EXPOSTATS

chemin(
  fileName = "webexpo.seg.mainbayesian.mcgill.R",
  relPath = c("R MODELS", "WRAPPING FUNCTIONS"),
  githubCredentials = list(userName = "webexpo", repoName = "webexpo_r_lib")
)


chemin(
  fileName = "data-summary.R",
  relPath = c("R MODELS", "McGILL FUNCTIONS"),
  githubCredentials = list(userName = "webexpo", repoName = "webexpo_r_lib")
)

chemin(
  fileName = "fcts.R",
  relPath = c("R MODELS", "McGILL FUNCTIONS"),
  githubCredentials = list(userName = "webexpo", repoName = "webexpo_r_lib")
)

chemin(
  fileName = "model-SEG-informedvar.R",
  relPath = c("R MODELS", "McGILL FUNCTIONS"),
  githubCredentials = list(userName = "webexpo", repoName = "webexpo_r_lib")
)

chemin(
  fileName = "model-SEG-riskband.R",
  relPath = c("R MODELS", "McGILL FUNCTIONS"),
  githubCredentials = list(userName = "webexpo", repoName = "webexpo_r_lib")
)

chemin(
  fileName = "model-SEG-uninformative.R",
  relPath = c("R MODELS", "McGILL FUNCTIONS"),
  githubCredentials = list(userName = "webexpo", repoName = "webexpo_r_lib")
)

# JAGS models UNINFORMATIVE

chemin(
  fileName = "webexpo.seg.uninformativebayesian.R",
  relPath = c("jags models", "SEG ANALYSIS"),
  githubCredentials = list(userName = "webexpo", repoName = "webexpo_r_lib")
)

chemin(
  fileName = "webexpo.seg.uninformativebayesian.models.R",
  relPath = c("jags models", "SEG ANALYSIS"),
  githubCredentials = list(userName = "webexpo", repoName = "webexpo_r_lib")
)


chemin(
  fileName = "webexpo.between.mainbayesian.R",
  relPath = c("jags models", "BETWEEN WORKER ANALYSIS"),
  githubCredentials = list(userName = "webexpo", repoName = "webexpo_r_lib")
)


## STAN WEBEXPO functions and models ( informedvar and uninformative, lognormal, no ME and ME as CV )

source("https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/McGILL%20FUNCTIONS/SEG-informedVar-stan.R")
source("https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/McGILL%20FUNCTIONS/SEG-uninformative-stan.R")
source("https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/McGILL%20FUNCTIONS/stan-fcts.R")
source("https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/WRAPPING%20FUNCTIONS/webexpo.seg.mainbayesian.stan.R")

#uninformative models code
code_seg_uninformative <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/seg_uninformative.stan'
code_seg_uninformative_lognormal_mecv <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/seg_uninformative_lognormal_mecv.stan'
code_seg_uninformative_lognormal_mecvknown <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/seg_uninformative_lognormal_mecvknown.stan'

#informedvar models code
code_seg_informedvar <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/seg_informedvar.stan'
code_seg_informedvar_lognormal_mecv <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/seg_informedvar_lognormal_mecv.stan'
code_seg_informedvar_lognormal_mecvknown <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/seg_informedvar_lognormal_mecvknown.stan'


## Data interpretation

chemin(
  fileName = "webexpo.seg.interpretation.R",
  relPath = c("RESULT INTERPRETATION", "SEG ANALYSIS"),
  githubCredentials = list(userName = "webexpo", repoName = "webexpo_r_lib")
) 

# ROS

chemin(
  fileName = "function NDexpo_EN.R",
  relPath = c("NDEXPO"),
  githubCredentials = list(userName = "webexpo", repoName = "webexpo_r_lib")
) 
