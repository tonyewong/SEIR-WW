##==============================================================================
## sobol_wrapper.R
##
## This file contains wrapper functions to call the model `seir_model.R` within
## the Sobol' sensitivity analysis. Importantly, they will call the model, and
## compute the sensitivity output of interest (total infections, in the Wong et
## al manuscript), but you can set this to something else if you'd like.
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================

seir_model_wrapper <- function(parameters, forcings, parameter_names, forcing_names, cycles.per.day, n.cycle, sens) {

  # call model
  df <- seir_model(parameters, forcings, parameter_names, forcing_names, cycles.per.day, n.cycle)

  # compute diagnostic things you might want as output
  number_tested <- sum(df[2:nrow(df),"Persons Tested"], na.rm = TRUE) + as.numeric(df[nrow(df),"Wastewater Screening Count"])
  number_confirmatory_tests <- sum(df[2:nrow(df),"Total TPs"], na.rm = TRUE) + sum(df[2:nrow(df),"Total FPs"], na.rm = TRUE)
  average_iu_census <- sum(df[2:nrow(df),"Total"], na.rm=TRUE)/sum(!is.na(df[2:nrow(df),"Total"])) # the denominator is counting up elements that are not NA
  average_pct_isolated <- 1-((sum(df[2:nrow(df),"False Positive"], na.rm=TRUE)/sum(!is.na(df[2:nrow(df),"False Positive"])))/average_iu_census)
  testing_cost <- number_tested*forcings[1,"test_cost"] + number_confirmatory_tests*forcings[1,"confirmatory_test_cost"]
  infections <- df[nrow(df),"Cumulative Infections"]

  # total number of students who have gone into isolation/quarantine
  total_isolation_headcount <- as.numeric(df[nrow(df),"Symptoms"]+df[nrow(df),"TP"]+df[nrow(df),"FP"])

  # some options for the sensitivity metric used:
  #   total number of persons tested
  #   = number_tested
  # or
  #   average isolation pool size
  #   = average_iu_census
  # or
  #   total number of infections
  #   = infections
  # or
  #   total cost
  #   = testing_cost

  if (sens=="infections") {
    output_of_interest <- as.numeric(infections)
  } else if (sens=="number_tested") {
    output_of_interest <- as.numeric(number_tested)
  } else if (sens=="average_isolation") {
    output_of_interest <- as.numeric(average_iu_census)
  } else if (sens=="total_isolation") {
    output_of_interest <- as.numeric(total_isolation_headcount)
  } else if (sens=="total_cost") {
    output_of_interest <- as.numeric(testing_cost)
  } else {
    output_of_interest <- as.numeric(infections) # default to infections
  }

  return(output_of_interest)
}


## Sobol' wrapper - assumes uniform distributions on parameters
seir_sobol <- function(parameters_matrix, forcings, parameter_names, forcing_names, cycles.per.day, n.cycle, sens) {
  output <- sapply(1:nrow(parameters_matrix), function(ii) {
                  seir_model_wrapper(matrix(parameters_matrix[ii,], nr=1, nc=length(parameter_names)), forcings, parameter_names, forcing_names, cycles.per.day, n.cycle, sens)})
  output_centered <- output - mean(output)
  return(output_centered)
}

##==============================================================================
## End
##==============================================================================
