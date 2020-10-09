##==============================================================================
## balance_check.R
##
## Run a simulation to check that the number of indivudals in the overall
## population remains the same throughout the simulation.
##
## Requires the input file `input/balance_check_set_parameters.csv`, because
## this check will be performed on all of the simulations from that file.
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================

rm(list=ls())

library(scales)
library(tidyverse)

##==============================================================================
## Set stuff

if (Sys.info()["user"]=="tony") {setwd("/Users/tony/codes/SEIR-WW/R")}
if (Sys.info()["user"]=="aewsma") {setwd("/Users/aewsma/codes/SEIR-WW/R")}

# set up model configuration
field_to_plot <- "Total"
N_day_running_total <- 14
cycles.per.day <- 3
n.days <- 100
n.cycle <- cycles.per.day*n.days

source("seir_model.R")

# set up inputs
inputs_ensemble <-read.csv(file="../input/balance_check_set_parameters.csv", header=FALSE)
n.parameters <- nrow(inputs_ensemble)
n.ensemble <- ncol(inputs_ensemble) - 1
parameter_names <- inputs_ensemble[,1]
parameters_ensemble <- as.matrix(t(inputs_ensemble[,2:(n.ensemble+1)]))
colnames(parameters_ensemble) <- parameter_names
forcing_names <- forcings <- NULL # all are parameters here

# run the actual simulation
model_out <- vector('list', length=n.ensemble)
running_max_Nday <- running_max_total <- running_max_infections <- running_max_quarantine <- 0
for (ii in 1:n.ensemble) {
  df <- seir_model(parameters=parameters_ensemble[ii,], forcings=forcings,
                   parameter_names=parameter_names, forcing_names=forcing_names,
                   cycles.per.day=cycles.per.day, n.cycle=n.cycle)
  number_tested <- sum(df[2:nrow(df),"Persons Tested"], na.rm = TRUE) + as.numeric(df[nrow(df),"Wastewater Screening Count"])
  number_confirmatory_tests <- sum(df[2:nrow(df),"Total TPs"], na.rm = TRUE) + sum(df[2:nrow(df),"Total FPs"], na.rm = TRUE)
  average_iu_census <- sum(df[2:nrow(df),"Total"], na.rm=TRUE)/sum(!is.na(df[2:nrow(df),"Total"])) # the denominator is counting up elements that are not NA
  average_pct_isolated <- 1-((sum(df[2:nrow(df),"False Positive"], na.rm=TRUE)/sum(!is.na(df[2:nrow(df),"False Positive"])))/average_iu_census)
  testing_cost <- number_tested*parameters_ensemble[ii,"test_cost"] + number_confirmatory_tests*parameters_ensemble[ii,"confirmatory_test_cost"]
  infections <- as.numeric(df[nrow(df),"Cumulative Infections"])
  total_isolation_headcount <- as.numeric(df[nrow(df),"Symptoms"]+df[nrow(df),"TP"]+df[nrow(df),"FP"])
  total_infections_Nday <- rep(0,n.cycle+1)
  model_out[[ii]] <- mat <- as.matrix(df)
  for (tt in 1:(N_day_running_total*cycles.per.day)) {
    total_infections_Nday[tt] <- sum(mat[1:tt, "New Infections"])
  }
  for (tt in (N_day_running_total*cycles.per.day+1):(n.cycle+1)) {
    total_infections_Nday[tt] <- total_infections_Nday[tt-1] + mat[tt, "New Infections"] - mat[(tt-N_day_running_total*cycles.per.day), "New Infections"]
  }
  #total_infections_Nday[1:(N_day_running_total*cycles.per.day)] <- 0
  model_out[[ii]] <- cbind(model_out[[ii]], total_infections_Nday)
  colnames(model_out[[ii]])[ncol(model_out[[ii]])] <- "Total Count Nday"
  print(paste("Total infections:",round(infections),", Total dead:",round(model_out[[ii]][n.cycle+1,"Dead"],3),", Max. ",N_day_running_total,"-day:",round(max(model_out[[ii]][,"Total Count Nday"]),1),", Number tested:",round(number_tested),", Total cost:",round(testing_cost)))
  running_max_Nday <- max(running_max_Nday, max(model_out[[ii]][,"Total Count Nday"]))
  running_max_total <- max(running_max_total, max(model_out[[ii]][,field_to_plot]))
  running_max_infections <- max(running_max_infections, max(model_out[[ii]][,"Cumulative Infections"]))
  running_max_quarantine <- max(running_max_quarantine, max(model_out[[ii]][,"False Positive"]+model_out[[ii]][,"True Positive"]+model_out[[ii]][,"Symptoms"]))
}
days <- (0:n.cycle)/cycles.per.day

##==============================================================================

for (ii in 1:n.ensemble) {
  model_census <- model_out[[ii]][,"Susceptible"] + model_out[[ii]][,"FP"] +
                  model_out[[ii]][,"TP"]          + model_out[[ii]][,"Exposed"] +
                  model_out[[ii]][,"Asympt"]      + model_out[[ii]][,"Symptoms"] +
                  model_out[[ii]][,"Recovered"]   + model_out[[ii]][,"Dead"]
  max_imbalance <- max(abs(model_census - parameters_ensemble[ii,"population_size"]))
  print(paste(ii,"maximum imbalance:",max_imbalance))
}

##==============================================================================
## End
##==============================================================================
