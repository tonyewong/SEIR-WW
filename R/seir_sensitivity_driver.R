##==============================================================================
## seir_sensitivity_driver.R
##
## This is a driver script to run all three of the sensitivity analysis cases
## described in the Wong et al. manuscript. The first, RIT control case, is the
## scenario that the main text focuses on. The other two are for the
## supplemental experiments for a large university and a small college.
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================

rm(list=ls())

if (Sys.info()["user"]=="tony") {setwd("/Users/tony/codes/SEIR-WW/R")}
if (Sys.info()["user"]=="aewsma") {setwd("/Users/aewsma/codes/SEIR-WW/R")}

source("seir_sensitivity.R")

# RIT control case
sens <- "infections"
n_sample <- 10000
n_bootstrap <- 1000
second <- TRUE
pri <- "mid"
population <- 12500
appen <- paste("sens",sens,"_ns",n_sample,"_nb",n_bootstrap,"_second",second,"_priors",pri,"_",sep="")
seir_sensitivity(n_sample=n_sample, n_bootstrap=n_bootstrap, sens=sens, second=second, pri=pri, population=population, appen=appen)

# large campus scenario
sens <- "infections"
n_sample <- 10000
n_bootstrap <- 1000
second <- TRUE
pri <- "high"
population <- 20000
appen <- paste("sens",sens,"_ns",n_sample,"_nb",n_bootstrap,"_second",second,"_priors",pri,"_",sep="")
seir_sensitivity(n_sample=n_sample, n_bootstrap=n_bootstrap, sens=sens, second=second, pri=pri, population=population, appen=appen)

# small campus scenario
sens <- "infections"
n_sample <- 10000
n_bootstrap <- 1000
second <- TRUE
pri <- "small"
population <- 3000
appen <- paste("sens",sens,"_ns",n_sample,"_nb",n_bootstrap,"_second",second,"_priors",pri,"_",sep="")
seir_sensitivity(n_sample=n_sample, n_bootstrap=n_bootstrap, sens=sens, second=second, pri=pri, population=population, appen=appen)


##==============================================================================
## End
##==============================================================================
