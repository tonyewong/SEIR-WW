##==============================================================================
## plot_sample_size.R
##
## Plotting estimate of reliability based on subsamples of increasing length, to
## verify that the number of samples used provides a stable estimate of these
## probabilities.
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================

rm(list=ls())

##
## Running this exactly like the actual `seir_ensemble_risk.R` calculations proceed
##

##==============================================================================
## Set stuff for model configuration

n.ensemble <- 5000 # how many simulations for each case (this is higher for creating this figure)
cycles.per.day <- 3
n.days <- 100
n.cycle <- cycles.per.day*n.days
N_day_running_total <- 14

failure_cases <- expand.grid(new_infections_per_shock=c(15, 30),
                             Rt=c(1.1, 1.5),
                             frac_noncompliance=c(0.01, 0.1))

screening_cases <- expand.grid(frequency_of_screening=c(999999),
                               frequency_of_screening_wastewater=c(4),
                               lag_screening_wastewater=c(1))

##==============================================================================

library(scales)
library(tidyverse)
library(lhs)
library(triangle)
library(extraDistr)

if (Sys.info()["user"]=="tony") {setwd("/Users/tony/codes/SEIR-WW/R")}
if (Sys.info()["user"]=="aewsma") {setwd("/Users/aewsma/codes/SEIR-WW/R")}

source("seir_model.R")
source("seir_wrapper.R")

# set up input forcings (scenario-based, not sampled uncertainties)
# -- comment out and add to `parameter_names` below if you want any of these in the analysis as uncertainties
forcing_names <- c("exogenous_shocks"
                   ,"frequency_exogenous_shocks"
                   #,"test_sensitivity", "test_specificity"
                   ,"test_cost", "time_to_return_fps"
                   ,"confirmatory_test_cost", "population_size"
                   ,"wastewater_sampling"
                   ,"frequency_of_screening", "frequency_of_screening_wastewater", "lag_screening_wastewater"
                   ,"new_infections_per_shock", "Rt", "frac_noncompliance"
)
forcings <- mat.or.vec(nr=1, nc=length(forcing_names))
colnames(forcings) <- forcing_names
#forcings[1,"test_sensitivity"] <- 0.8
#forcings[1,"test_specificity"] <- 0.98
forcings[1,"exogenous_shocks"] <- 1             # 0 for no, 1 for yes
forcings[1,"frequency_exogenous_shocks"] <- 7   # days between shocks (outside infections)
forcings[1,"test_cost"] <- 25                   # dollars
forcings[1,"time_to_return_fps"] <- 1           # days to return a false positive to the Uninfected pool
forcings[1,"confirmatory_test_cost"] <- 100     # dollars
forcings[1,"population_size"] <- 12500          # size of total population
forcings[1,"wastewater_sampling"] <- 1          # 0 for no, 1 for yes
#forcings[1,"frequency_of_screening_wastewater"] <- 2     # time (days) to screen `building_size` individuals after wastewater positive
#forcings[1,"lag_screening_wastewater"] <- 2     # time between wastewater sample return and screenings

# set parameters
parameter_names <- c("initial_infected", "days_to_incubation", "time_to_recovery"#, "Rt"
                     #,"frequency_exogenous_shocks", "new_infections_per_shock"
                     ,"test_sensitivity", "test_specificity"
                     ,"pct_advancing_to_symptoms", "symptom_case_fatality_ratio"
                     #,"frequency_of_screening"
                     ,"frequency_wastewater_samples", "frac_wastewater", "threshold_wastewater"#, "lag_screening_wastewater"#, "frequency_of_screening_wastewater"
                     ,"building_size"
)
n_parameters <- length(parameter_names)

# initialize output
model_out <- vector("list", length=nrow(failure_cases))
total_experiments <- nrow(failure_cases)*nrow(screening_cases)
cnt <- 0

for (i_failure in c(1,6)) {
  forcings[1,"new_infections_per_shock"] <- failure_cases[i_failure,"new_infections_per_shock"]
  forcings[1,"Rt"] <- failure_cases[i_failure,"Rt"]
  forcings[1,"frac_noncompliance"] <- failure_cases[i_failure,"frac_noncompliance"]
  # initialize matrix for this failure case output
  model_out[[i_failure]] <- vector("list", length=2); names(model_out[[i_failure]]) <- c("forcings", "output")
  model_out[[i_failure]]$forcings <- forcings[1,c("new_infections_per_shock","Rt","frac_noncompliance")] # save, just in case there is any confusion
  #model_out[[i_failure]]$output <- rep(NA, nrow(screening_cases))
  model_out[[i_failure]]$output <- vector("list", length=nrow(screening_cases))

  for (i_screening in 1:nrow(screening_cases)) {
    forcings[1,"frequency_of_screening"] <- screening_cases[i_screening,"frequency_of_screening"]
    forcings[1,"frequency_of_screening_wastewater"] <- screening_cases[i_screening,"frequency_of_screening_wastewater"]
    forcings[1,"lag_screening_wastewater"] <- screening_cases[i_screening,"lag_screening_wastewater"]

    cnt <- cnt+1

    # sample and scale the parameters
    set.seed(2019) # <-- get the same sample each time so result is conditioned on it and reproducible
    sample_lhs <- randomLHS(n.ensemble, n_parameters)
    source("priors_mid.R")
    source("scale_parameters.R")
    priors <- priors_func(parameter_names)
    parameters_lhs <- scale_parameters(sample_lhs, parameter_names, priors)
    colnames(parameters_lhs) <- parameter_names

    # run the simulations
    tbeg <- proc.time()
    #pb <- txtProgressBar(min=0,max=n.ensemble,initial=0,style=3)
    max_perceived <- max_actual <- all_infections <- rep(NA, n.ensemble)
    for (ii in 1:n.ensemble) {
      df <- seir_model(parameters=parameters_lhs[ii,], forcings=forcings,
                       parameter_names=parameter_names, forcing_names=forcing_names,
                       cycles.per.day=cycles.per.day, n.cycle=n.cycle)
      number_tested <- sum(df[2:nrow(df),"Persons Tested"], na.rm = TRUE) + as.numeric(df[nrow(df),"Wastewater Screening Count"])
      number_confirmatory_tests <- sum(df[2:nrow(df),"Total TPs"], na.rm = TRUE) + sum(df[2:nrow(df),"Total FPs"], na.rm = TRUE)
      average_iu_census <- sum(df[2:nrow(df),"Total"], na.rm=TRUE)/sum(!is.na(df[2:nrow(df),"Total"])) # the denominator is counting up elements that are not NA
      average_pct_isolated <- 1-((sum(df[2:nrow(df),"False Positive"], na.rm=TRUE)/sum(!is.na(df[2:nrow(df),"False Positive"])))/average_iu_census)
      testing_cost <- number_tested*forcings[1,"test_cost"] + number_confirmatory_tests*forcings[1,"confirmatory_test_cost"]
      infections <- as.numeric(df[nrow(df),"Cumulative Infections"])
      total_isolation_headcount <- as.numeric(df[nrow(df),"Symptoms"]+df[nrow(df),"TP"]+df[nrow(df),"FP"])

      mat <- as.matrix(df)

      # new FP going into time steps 2:(n.cycle+1)
      new_false_positives <- diff(mat[,"False Positive"])
      new_false_positives[new_false_positives < 0] <- 0

      # new TP going into time steps 2:(n.cycle+1) computed as the change, plus the outgoing
      # "fluxes" in time steps 1:n.cycle (since change = flux_in - flux_out), and we
      # want flux_in
      rho <- 1/(parameters_lhs[ii,"time_to_recovery"]*cycles.per.day)
      sigma <- rho*(parameters_lhs[ii,"pct_advancing_to_symptoms"]/100/(1-parameters_lhs[ii,"pct_advancing_to_symptoms"]/100))
      new_true_positives <- diff(mat[,"True Positive"]) + mat[1:n.cycle,"True Positive"]*rho + mat[1:n.cycle,"True Positive"]*sigma
      new_true_positives[new_true_positives < 0] <- 0

      # new S going into time steps 2:(n.cycle+1) computed as incoming formerly
      # asymptomatic cases only. the former TPs were already counted previously
      new_symptomatic <- mat[1:n.cycle,"Asympt"]*sigma*(1-forcings[1,"frac_noncompliance"])

      # add them all up to get the total number of new (perceived) infections
      new_infections <- c(0,new_false_positives + new_true_positives + new_symptomatic)

      # take N-day running total
      total_infections_Nday <- rep(0,n.cycle+1)
      for (tt in 1:(N_day_running_total*cycles.per.day)) {
        total_infections_Nday[tt] <- sum(new_infections[1:tt])
      }
      for (tt in (N_day_running_total*cycles.per.day+1):(n.cycle+1)) {
        total_infections_Nday[tt] <- total_infections_Nday[tt-1] + new_infections[tt] - new_infections[tt-(N_day_running_total*cycles.per.day)]
      }
      max_perceived[ii] <- max(total_infections_Nday)

      # add up the actual infections too
      actual_infections_Nday <- rep(0,n.cycle+1)
      for (tt in 1:(N_day_running_total*cycles.per.day)) {
        actual_infections_Nday[tt] <- sum(mat[1:tt, "New Infections"])
      }
      for (tt in (N_day_running_total*cycles.per.day+1):(n.cycle+1)) {
        actual_infections_Nday[tt] <- actual_infections_Nday[tt-1] + mat[tt, "New Infections"] - mat[(tt-N_day_running_total*cycles.per.day), "New Infections"]
      }
      max_actual[ii] <- max(actual_infections_Nday)

      # save all infections
      all_infections[ii] <- infections

      #setTxtProgressBar(pb, ii)
    }
    #close(pb)
    tend <- proc.time()
    #print(paste("Took",round((tend-tbeg)[3]/60,1),"minutes"))
    model_out[[i_failure]]$output[[i_screening]] <- c(length(which(max_perceived <= 100))/length(max_perceived),
                                                      length(which(max_actual <= 100))/length(max_actual),
                                                      median(all_infections))
  }
}
days <- (0:n.cycle)/cycles.per.day



est1 <- est2 <- est3 <- rep(NA, n.ensemble)
for (ii in 1:n.ensemble) {
  est1[ii] <- length(which(max_perceived[1:ii] <= 100))/ii
  est2[ii] <- length(which(max_actual[1:ii] <= 100))/ii
  est3[ii] <- median(all_infections[1:ii])
}

panel_labels <- c(expression(bold('a')),
                  expression(bold('b')),
                  expression(bold('c')))

pdf('../figures/sample_size.pdf',width=3.5,height=7,colormodel='cmyk', pointsize=11)
par(mfrow=c(3,1), mai=c(.6,.6,.2,.1))
plot(est1, type='l', main="Perceived reliability infections over 14-days < 100", ylab="Reliability", xlab="Sample size")
mtext(side=3, text=panel_labels[1], line=0, cex=1, adj=0);
plot(est2, type='l', main="Actual reliability infections over 14-days < 100", ylab="Reliability", xlab="Sample size")
mtext(side=3, text=panel_labels[2], line=0, cex=1, adj=0);
plot(est3, type='l', main="Actual total infections", ylab="Count", xlab="Sample size")
mtext(side=3, text=panel_labels[3], line=0, cex=1, adj=0);
dev.off()


##==============================================================================
## End
##==============================================================================
