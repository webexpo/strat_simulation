#' loads the Webexpo functions necessary to generate exposure levels

devtools::source_url("https://github.com/lhimp/scripts/raw/master/chemin.R")


# Webexpo functions  -----------------------------------------------------------------

## data generation for one SEG

chemin(
  fileName = "webexpo.seg.randomgeneration.R",
  relPath = c("RANDOM SAMPLE GENERATION"),
  githubCredentials = list(userName = "webexpo", repoName = "webexpo_r_lib")
)



# data generation for variability between workers within one SEG

chemin(
  fileName = "webexpo.between.randomgeneration.R",
  relPath = c("RANDOM SAMPLE GENERATION"),
  githubCredentials = list(userName = "webexpo", repoName = "webexpo_r_lib")
)

