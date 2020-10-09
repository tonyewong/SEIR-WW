##==============================================================================
## seir_ensemble.R
##
## Run a set of SEIR model simulations. This function is not used in the results
## from the Wong et al. manuscript, but is provided to allow testing of the
## model to compare infection rates, expected numbers of deaths, numbers of
## screening tests used, etc. between model configurations.
## * Requires the input file `input/ensemble_set_parameters.csv`
##
## If you want to run your own test cases, you can modify the parameter values
## in the input `ensemble_set_parameters.csv` file. Then, run this script
## interactively in order to noodle around with the results.
##
## Here, no effective difference between "forcings" and "parameters". That only
## makes a difference in the sensitivity analysis, where forcings are taken to
## be the same from simulation to simulation, whereas the parameters are taken
## to be uncertain and vary across simulations according to the distributions
## given in `priors.R`.
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
inputs_ensemble <-read.csv(file="../input/ensemble_set_parameters.csv", header=FALSE)
n.parameters <- nrow(inputs_ensemble)
n.ensemble <- ncol(inputs_ensemble) - 1
parameter_names <- inputs_ensemble[,1]
parameters_ensemble <- as.matrix(t(inputs_ensemble[,2:(n.ensemble+1)]))
colnames(parameters_ensemble) <- parameter_names
forcing_names <- forcings <- NULL # all are parameters here

# run the actual simulation
model_out <- vector('list', length=n.ensemble)
running_max_Nday <- running_max_total <- 0
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
}
days <- (0:n.cycle)/cycles.per.day

# make a plot of the total known (or FP) infections at each time
par(mfrow=c(3,1), mai=c(.5,.5,.2,.2), las=0)
if (n.ensemble > 8) {
  print("Not generating plot. It would be too cluttered to be useful.")
} else {
  plot_colors <- c("black","steelblue","firebrick","seagreen","darkorchid4","goldenrod","hotpink2","coral")
  ii <- 1
  # detailed plot only for the first simulation set
  plot(days, model_out[[ii]][,"Recovered"], type='l', lwd=1.5, col='seagreen', xlab='', ylab='')
  lines(days, model_out[[ii]][,"Cumulative Infections"], type='l', lwd=1.5, col='black')
  lines(days, model_out[[ii]][,"Exposed"], type='l', lwd=1.5, col='orange')
  lines(days, model_out[[ii]][,"Symptoms"], type='l', lwd=1.5, col='firebrick')
  mtext(side=1, text="Day", line=2)
  mtext(side=2, text="Count", line=2.2)
  grid()
  axes <- par("usr")
  legend(60,axes[4],c("Infec","Recov","Expos","Sympt"), lty=1, lwd=1.5, col=c("black","seagreen","orange","firebrick"))
  # plot of 14-day infection count
  plot(days, model_out[[ii]][,"Total Count Nday"], type='l', lwd=1.5, col=plot_colors[ii], ylim=c(0,running_max_Nday*1.25), xlab="", ylab="")
  mtext(side=1, text="Day", line=2)
  mtext(side=2, text=paste("Infections over previous",N_day_running_total,"days"), line=2.2)
  if (n.ensemble > 1) {
    for (ii in 2:n.ensemble) {lines(days, model_out[[ii]][,"Total Count Nday"], lwd=1.5, col=plot_colors[ii])}
    legend(n.days*0.65, 1.3*running_max_Nday, legend=1:n.ensemble, lty=1, lwd=1.5, col=plot_colors[1:n.ensemble], bty="n")
  }
  lines(c(0, max(days)), c(100,100), lty=2, col="gray50", lwd=1.5)
  # plot of something else you wanted (`field_to_plot` set above)
  ii <- 1
  plot(days, model_out[[ii]][,field_to_plot], type='l', lwd=1.5, col=plot_colors[ii], ylim=c(0,running_max_total*1.1), xlab="", ylab="")
  mtext(side=1, text="Day", line=2)
  mtext(side=2, text=field_to_plot, line=2.2)
  lines(days, model_out[[ii]][,"Wastewater Screening"]*50, col=plot_colors[ii], lty=3)
  if (n.ensemble > 1) {
    for (ii in 2:n.ensemble) {lines(days, model_out[[ii]][,field_to_plot], lwd=1.5, col=plot_colors[ii]); lines(days, model_out[[ii]][,"Wastewater Screening"]*50, col=plot_colors[ii], lty=3)}
  }
}

# save relevant parameters and configuration information
#save(model_out, parameters_ensemble, parameter_names, days, file = "../output/ensemble_results.RData")

##==============================================================================
## End
##==============================================================================
