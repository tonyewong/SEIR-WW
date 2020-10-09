##==============================================================================
## seir_model.R
##
## This function runs the actual SEIR model with the optional wastewater
## surveillance module. Returns `df` a data frame with the subpopulation sizes
## for the different model compartments, numbers of infections, etc.
##
## Some key parameters:
##  wastewater_sampling = 0 if no wastewater sampling, 1 if doing wastewater sampling
##  frequency_wastewater_testing = number of days between wastewater testing results returns
##  frac_wastewater = proportion of asymptomatic individuals identifiable through wastewater sampling - what fraction of the overall population is identifiable through waterwater sampling?
##  frequency_of_screening_wastewater = frequency of screening after positive detection in wastewater - how many days does it take to test everyone in groups identifiable through wastewater sampling?
##  threshold_wastewater = how many asymptomatic people must be present to detect a "hit" in wastewater?
##
## Function inputs:
##  parameters = 1d array of parameter values
##  forcings = 1d array of forcing input values
##  parameter_names = 1d array of parameter names (order matches `parameters`)
##  forcing_names = 1d array of forcing input names (order matches `forcings`)
##  cycles.per.day = number of model time steps per calendar day
##  n.cycle = number of model time steps in the entire simulation
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================

seir_model <- function(parameters, forcings, parameter_names, forcing_names, cycles.per.day, n.cycle) {

  # round things up as `input` list, as in the dashboard code. not done to be
  # efficient, but rather as a convenience so the (working) dashboard code is
  # minimally changed relative to what was done above.
  input_names <- c("initial_infected", "population_size"
                  ,"Rt", "exogenous_shocks", "frequency_exogenous_shocks", "new_infections_per_shock"
                  ,"days_to_incubation", "time_to_recovery", "pct_advancing_to_symptoms", "symptom_case_fatality_ratio"
                  ,"frequency_of_screening", "test_sensitivity", "test_specificity", "test_cost", "time_to_return_fps", "confirmatory_test_cost"
                  ,"wastewater_sampling", "frequency_wastewater_samples", "frac_wastewater", "frequency_of_screening_wastewater", "threshold_wastewater", "lag_screening_wastewater"
                  ,"building_size", "frac_noncompliance"
                  )
  n_input <- length(input_names)
  input <- vector('list', length=n_input)
  names(input) <- input_names
  for (nn in input_names) {
    if (nn %in% parameter_names) {      input[[nn]] <- parameters[match(nn, parameter_names)]
    } else if (nn %in% forcing_names) { input[[nn]] <- forcings[match(nn, forcing_names)]
    } else {                            print(paste("ERROR: Input field",nn,"is not matched."))}
  }
  input$initial_susceptible <- input$population_size - input$initial_infected

  num.exogenous.shocks <- input$exogenous_shocks # 0 for no, 1 for yes
  num.wastewater.sampling <- input$wastewater_sampling

  frequency.exogenous.shocks <- 1+cycles.per.day*input$frequency_exogenous_shocks          # adding 1 to keep away from 0 (continuous shocks would be... problematic)
  frequency.wastewater.samples <- 1+cycles.per.day*input$frequency_wastewater_samples      # this is for how often we get WW sample results
  lag.screening.wastewater <- cycles.per.day*input$lag_screening_wastewater                # how many time steps after wastewater return do we start receiving test results back?
  cycles.per.test <- input$frequency_of_screening*cycles.per.day                           # this is for the regular covid surveillance testing
  cycles.per.test.wastewater <- input$frequency_of_screening_wastewater*cycles.per.day     # this is for the actual covid testing after getting a WW hit
  rho <- 1/(input$time_to_recovery*cycles.per.day)
  sigma <- rho*(input$pct_advancing_to_symptoms/100/(1-input$pct_advancing_to_symptoms/100))
  beta <- input$Rt*(rho+sigma)
  delta <- (input$symptom_case_fatality_ratio/(1-input$symptom_case_fatality_ratio))*rho
  theta <- 1/(input$days_to_incubation*cycles.per.day)
  mu <- 1/(cycles.per.day*input$time_to_return_fps)

  # First time step of the model
  mat <- matrix(c(0,input$initial_susceptible,0,0,input$initial_infected,0,0,0,0), nrow = 1)
  colnames(mat) <- c("i","U","FP","E","A","S","TP","R","D")
  mat <-
    rbind(
      mat,
      c(1,
        max(0,mat[1,2]*(1-beta*(mat[1,5]/(mat[1,2]+mat[1,5]+mat[1,4])))+mat[1,3]*mu),
        #max(0,mat[1,2]*(1-beta*(mat[1,5]/(mat[1,2]+mat[1,5]+mat[1,4]+mat[1,8])))+mat[1,3]*mu),
        max(0,mat[1,3]*(1-mu)),
        max(0,mat[1,4]*(1-theta)+ beta*(mat[1,2]*mat[1,5]/(mat[1,2]+mat[1,5]+mat[1,4]))),
        #max(0,mat[1,4]*(1-theta)+ beta*(mat[1,2]*mat[1,5]/(mat[1,2]+mat[1,5]+mat[1,4]+mat[1,8]))),
        max(0,mat[1,5]*(1-sigma*(1-input$frac_noncompliance)-rho)+mat[1,4]*theta),
        max(0,mat[1,6]*(1-delta-rho)+(mat[1,5]*(1-input$frac_noncompliance)+mat[1,7])*sigma),
        0,
        max(0,mat[1,8]+(mat[1,5]+mat[1,6]+mat[1,7])*rho),
        max(0,delta*mat[1,6]+mat[1,9]))
    )

  if (frequency.exogenous.shocks==0) {
    superspreader.event <- rep(0, n.cycle+1)
  } else {
    superspreader.event <- 0
    superspreader.event <- c(superspreader.event,
                              (1:n.cycle %% frequency.exogenous.shocks == 0)*num.exogenous.shocks)
  }

  if (frequency.wastewater.samples==0) {
    wastewater.sample <- rep(0, n.cycle+1)
  } else {
    wastewater.sample <- 0
    wastewater.sample <- c(wastewater.sample,
                              (1:n.cycle %% frequency.wastewater.samples == 0)*num.wastewater.sampling)
  }

  wastewater.screening <- rep(0, n.cycle+1)      # this will be a 1 if we have triggered screenings based on the wastewater sampling results
  wastewater.screening.count <- c(0,0)           # this is cumulative number of wastewater-based screening tests conducted

  # Subsequent time steps
  for (i in 2:n.cycle) {
    mat <- rbind(mat,
                c(i,
                  max(0,mat[i,2]*(1-beta*(mat[i,5]/(mat[i,2]+mat[i,5]+mat[i,4])))+mat[i,3]*mu-mat[i-1,2]*(1-input$test_specificity)/cycles.per.test-superspreader.event[i+1]*input$new_infections_per_shock),
                  #max(0,mat[i,2]*(1-beta*(mat[i,5]/(mat[i,2]+mat[i,5]+mat[i,4]+mat[i,8])))+mat[i,3]*mu-mat[i-1,2]*(1-input$test_specificity)/cycles.per.test-superspreader.event[i+1]*input$new_infections_per_shock),
                  max(0,mat[i,3]*(1-mu)+mat[i-1,2]*(1-input$test_specificity)/cycles.per.test),
                  max(0,mat[i,4]*(1-theta)+beta*(mat[i,2]*mat[i,5]/(mat[i,2]+mat[i,5]+mat[i,4]))+superspreader.event[i+1]*input$new_infections_per_shock),
                  #max(0,mat[i,4]*(1-theta)+beta*(mat[i,2]*mat[i,5]/(mat[i,2]+mat[i,5]+mat[i,4]+mat[i,8]))+superspreader.event[i+1]*input$new_infections_per_shock),
                  max(0,mat[i,5]*(1-sigma*(1-input$frac_noncompliance)-rho)+mat[i,4]*theta-mat[i-1,5]*(1-input$frac_noncompliance)*input$test_sensitivity/cycles.per.test),
                  max(0,mat[i,6]*(1-delta-rho)+(mat[i,5]*(1-input$frac_noncompliance)+mat[i,7])*sigma),
                  max(0,mat[i,7]*(1-sigma-rho)+mat[i-1,5]*(1-input$frac_noncompliance)*input$test_sensitivity/cycles.per.test),
                  max(0,mat[i,8]+(mat[i,5]+mat[i,6]+mat[i,7])*rho),
                  max(0,delta*mat[i,6]+mat[i,9]))
                  )

    # testing as a result of a wastewater sample return
    # conditions that trigger wastewater-based screenings:
    # (1) get a sample back AND asymptomatic population sufficiently high
    # (2) continued screenings from previous trigger
    # Could implement probability of detection instead of the strict inequality
    # Only adding more wastewater screenings if we aren't already doing them.
    new_tests_wastewater <- 0
    #print(paste(mat[i,],superspreader.event[i+1],frequency.exogenous.shocks,num.exogenous.shocks))
    if (i > frequency.wastewater.samples) {
        if ( (wastewater.sample[i]==1) &&
             ((mat[i-frequency.wastewater.samples,4]+mat[i-frequency.wastewater.samples,5])*input$frac_wastewater >= input$threshold_wastewater)) {
            # triggered testing - make sure all `building_size` individuals are tested by setting wastewater.screening = 1
            # for however many time steps it will take to test everyone.
            if (i <= n.cycle+2-lag.screening.wastewater-cycles.per.test.wastewater) {
                wastewater.screening[(i+lag.screening.wastewater):(i+lag.screening.wastewater+cycles.per.test.wastewater-1)] <- 1
            } else if ( (i >  n.cycle+2-lag.screening.wastewater-cycles.per.test.wastewater) &&
                        (i <= n.cycle+1-lag.screening.wastewater) ) {
                wastewater.screening[(i+lag.screening.wastewater):(n.cycle+1)] <- 1
            } # no other cases because screening results wouldn't come back before end of the simulation
        }
    }
    if ( wastewater.screening[i+1]==1 ) {
        # * only test `building_size` number of people
        # * and `threshold_wastewater` is the minimum number of individuals who must be sick to show covid in WW from that building
        # * this assumes that `frac_wastewater` of all exposed+asymptomatic cases occur in wastewater surveilled buildings,
        #   and if we notice a spike, all are concentrated in a group of size `building_size` (possibly multiple buildings in real life)
        # * as in original model, E assumed to invariably test negative
        # * assuming that at time step i we got the wastewater return and start getting screening results (if we do screenings)
        #   `lag.screening.wastewater` time steps later
        # * building pop = (A[t]+E[t])*frac_wastewater + (U and R at proportions equal to general population)
        frac_U <- mat[i,2]/(mat[i,2]+mat[i,8])  # as with the regular testing, only grab FP from U and TP from A, leaving E alone
        U_ww <- ( max(0, input$building_size - (mat[i,4]+mat[i,5])*input$frac_wastewater)*frac_U * (1 - input$test_specificity) / cycles.per.test.wastewater) * num.wastewater.sampling
        A_ww <- ( mat[i,5] * input$frac_wastewater * (1-input$frac_noncompliance) * input$test_sensitivity / cycles.per.test.wastewater) * num.wastewater.sampling
        new_tests_wastewater <- input$building_size/cycles.per.test.wastewater
    } else {
        U_ww <- A_ww <- 0
        new_tests_wastewater <- 0
    }
    mat[i+1,2] <- mat[i+1,2] - U_ww  # these are false positives, so remove from U...
    mat[i+1,3] <- mat[i+1,3] + U_ww  # ... and add them to FP
    mat[i+1,5] <- mat[i+1,5] - A_ww  # these are true positives, so remove from A...
    mat[i+1,7] <- mat[i+1,7] + A_ww  # ... and add them to TP
    # building_size individuals are all tested; U_ww and A_ww are just the ones who were false or true positives
    wastewater.screening.count <- c(wastewater.screening.count, wastewater.screening.count[i] + new_tests_wastewater)
  }
  mat <- cbind(mat, superspreader.event, wastewater.sample, wastewater.screening, wastewater.screening.count)

  # create data frame for analysis
  names.df <- c("Cycle","Susceptible","FP","Exposed","Asympt","Symptoms","TP","Recovered","Dead","Superspreader Event","Wastewater Sample","Wastewater Screening","Wastewater Screening Count")
  df <-
    mat %>%
    as_tibble() %>%
    rename_all(~names.df) %>%
    mutate(`Persons Tested` = (lag(Susceptible,1,NA)+lag(Exposed,1,NA)+lag(Asympt,1,NA))/cycles.per.test,
            `Total TPs` = lag(Asympt,2,NA)*input$test_sensitivity/cycles.per.test,
            `Total FPs` = lag(Susceptible,2,NA)*(1-input$test_specificity)/cycles.per.test,
            `Total TNs` = lag(Susceptible,2,NA)*input$test_specificity/cycles.per.test,
            `Total FNs` = lag(Exposed,2,NA)+lag(Asympt,2,NA)*(1-input$test_sensitivity)/cycles.per.test) %>%
    mutate(Day = Cycle/cycles.per.day,
            `True Positive` = TP,
            Symptoms = Symptoms,
            `False Positive` = FP,
            Total = TP+Symptoms+FP) %>%
    mutate(`New Infections` = lag(Asympt,1,NA)*beta*lag(Susceptible,1,NA)/(lag(Susceptible,1,NA)+lag(Exposed,1,NA)+lag(Asympt,1,NA)),
            `New Infections` = ifelse(Cycle>1,
                                      `New Infections`+pmin(`Superspreader Event`*input$new_infections_per_shock,lag(Susceptible,1,NA)),
                                      `New Infections`),
            `New Infections` = ifelse(is.na(`New Infections`),0,`New Infections`),
            `Cumulative Infections` = cumsum(`New Infections`),
            `%Cumulative Infections` = `Cumulative Infections`/input$initial_susceptible)

  return(df)

}

# Notes:
# * Cumulative Infections = mat[301,"Exposed"]+mat[301,"Asympt"]+mat[301,"TP"]+mat[301,"Symptoms"]+mat[301,"Recovered"]+mat[301,"Dead"]
# * this counting of Persons Tested is not included the WW screenings
# * Total = FP + TP + Symptomatic, so that's the total current number of cases we know about at that time

##==============================================================================
## End
##==============================================================================
