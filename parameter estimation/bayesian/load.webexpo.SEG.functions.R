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



## Data interpretation

chemin(
  fileName = "webexpo.seg.interpretation.R",
  relPath = c("RESULT INTERPRETATION", "SEG ANALYSIS"),
  githubCredentials = list(userName = "webexpo", repoName = "webexpo_r_lib")
) 

