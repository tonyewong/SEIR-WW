##==============================================================================
## priors_mid.R
##
## Prior distributions for all the parameters. If you want to add new distributions,
## then you will need to add the appropriate CDF/PPF scaling to the `scale_parameters.R`
## file.
##
## This version is meant to represent the "control" case, for a university with
## about 12,500 individiuals on-campus and relatively low risk (e.g., Upstate
## New York, as of this writing).
##
##  parameter_names = 1d array of parameter names
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================

priors_func <- function(parameter_names) {

  n_parameters <- length(parameter_names)
  priors <- vector("list", length=n_parameters)
  names(priors) <- parameter_names

  for (pp in parameter_names) {

    # initialize list element for each parameter
    priors[[pp]] <- vector("list", 3)
    names(priors[[pp]]) <- c("distribution", "params", "parnames")

    # set the fields depending on which parameters are present - not assuming specific parameters are here

    if (pp == "initial_infected") {
      priors[[pp]]$distribution <- "binomial"
      priors[[pp]]$parnames <- c("size", "p")
      priors[[pp]]$params <- c(12500, 10/12500)

    } else if (pp == "Rt") {
      priors[[pp]]$distribution <- "triangle"
      priors[[pp]]$parnames <- c("a", "b", "c")
      priors[[pp]]$params <- c(0.8, 2.5, 1.1)

    } else if (pp == "days_to_incubation") {
      priors[[pp]]$distribution <- "triangle"
      priors[[pp]]$parnames <- c("a", "b", "c")
      priors[[pp]]$params <- c(3, 12, 5) # Paltiel et al originally had mean of 3, updated from Petersen et al 2020

    } else if (pp == "time_to_recovery") {
      priors[[pp]]$distribution <- "triangle"
      priors[[pp]]$parnames <- c("a", "b", "c")
      priors[[pp]]$params <- c(10, 21, 14) # limit from https://www.cdc.gov/coronavirus/2019-ncov/hcp/duration-isolation.html

    } else if (pp == "pct_advancing_to_symptoms") {
      priors[[pp]]$distribution <- "triangle"
      priors[[pp]]$parnames <- c("a", "b", "c")
      priors[[pp]]$params <- c(5, 50, 30) # lower limit and center from Paltiel dashboard, and defaults

    } else if (pp == "symptom_case_fatality_ratio") {
      priors[[pp]]$distribution <- "triangle"
      priors[[pp]]$parnames <- c("a", "b", "c")
      priors[[pp]]$params <- c(0, 0.01, 0.0005) # limits and center from Paltiel dashboard, and defaults

    } else if (pp == "frequency_of_screening") {
      priors[[pp]]$distribution <- "neg_binomial"
      priors[[pp]]$parnames <- c("size", "p")
      priors[[pp]]$params <- c(42, 1/8)   # to obtain about 300 tests/week

    } else if (pp == "frequency_exogenous_shocks") {
      priors[[pp]]$distribution <- "binomial"
      priors[[pp]]$parnames <- c("size", "p")
      priors[[pp]]$params <- c(14, 0.5)

    } else if (pp == "new_infections_per_shock") {
      priors[[pp]]$distribution <- "neg_binomial"
      priors[[pp]]$parnames <- c("size", "p")
      priors[[pp]]$params <- c(5, 0.25)   # gives 5-95% range of about 5-30 infections

    } else if (pp == "test_sensitivity") {
      priors[[pp]]$distribution <- "triangle"
      priors[[pp]]$parnames <- c("a", "b", "c")
      priors[[pp]]$params <- c(0.7, 0.9, 0.8)

    } else if (pp == "test_specificity") {
      priors[[pp]]$distribution <- "triangle"
      priors[[pp]]$parnames <- c("a", "b", "c")
      priors[[pp]]$params <- c(0.95, 1, 0.98)

    } else if (pp == "frequency_wastewater_samples") {
      priors[[pp]]$distribution <- "disc_uniform"
      priors[[pp]]$parnames <- c("min", "max")
      priors[[pp]]$params <- c(1, 8)

    } else if (pp == "frac_wastewater") {
      priors[[pp]]$distribution <- "triangle"
      priors[[pp]]$parnames <- c("a", "b", "c")
      priors[[pp]]$params <- c(0.2, 0.8, 0.55)

    } else if (pp == "lag_screening_wastewater") {
      priors[[pp]]$distribution <- "disc_uniform"
      priors[[pp]]$parnames <- c("min", "max")
      priors[[pp]]$params <- c(1, 5)

    } else if (pp == "frequency_of_screening_wastewater") {
      priors[[pp]]$distribution <- "disc_uniform"
      priors[[pp]]$parnames <- c("min", "max")
      priors[[pp]]$params <- c(1, 7)

    } else if (pp == "threshold_wastewater") {
      priors[[pp]]$distribution <- "triangle"
      priors[[pp]]$parnames <- c("a", "b", "c")
      priors[[pp]]$params <- c(1, 5, 2)  # not a problem to have non-integer threshold; interpret as an amount of viral RNA needed to be present

    } else if (pp == "building_size") {
      priors[[pp]]$distribution <- "neg_binomial"
      priors[[pp]]$parnames <- c("size", "p")
      priors[[pp]]$params <- c(5, 0.01)   # gives 5-95% range of about 200-900 individuals

    } else if (pp == "frac_noncompliance") {
      priors[[pp]]$distribution <- "triangle"
      priors[[pp]]$parnames <- c("a", "b", "c")
      priors[[pp]]$params <- c(0, 0.1, 0.01)

    } else {
      print(paste("ERROR: parameter",pp,"not recognized"))

    }

  }

  return (priors)

}

##==============================================================================
## End
##==============================================================================
