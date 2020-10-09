##==============================================================================
## seir_sensitivity.R
##
## Sobol' sensitivity analysis of SEIR model of Paltiel et al. 2020, with the
## wastewater sampling module described in Wong et al. manuscript. The essence
## of the Sobol' analysis is to decompose the variance in a model output of
## interest into portions attributable to different model input parameters, as
## well as the interactions among these input parameters.
##
## The total numbers of model evaluations requires are:
## -- If you estimate only first-order and total sensitivity indices: total of n_sample*(n_parameter + 2) model evaluations
## -- With 2nd-order sensitivity indices: total of n_sample*(2*n_parameter + 2) model evaluations
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================

seir_sensitivity <- function(n_sample, n_bootstrap=0, sens="infections", second=FALSE, pri="mid", population=12500, appen="") {

##==============================================================================

  library(scales)
  library(tidyverse)
  library(lhs)
  library(sensitivity)
  library(triangle)
  library(extraDistr)

  # set up input forcings (scenario-based, not sampled uncertainties)
  # -- comment out and add to `parameter_names` below if you want any of these in the sensitivity analysis as uncertainties
  forcing_names <- c("exogenous_shocks"
                     #,"frequency_exogenous_shocks", "new_infections_per_shock"
                     #,"test_sensitivity", "test_specificity"
                     ,"test_cost", "time_to_return_fps"
                     ,"confirmatory_test_cost", "population_size"
                     ,"wastewater_sampling"
                     #,"building_size"
                     )
  forcings <- mat.or.vec(nr=1, nc=length(forcing_names))
  colnames(forcings) <- forcing_names
  forcings[1,"exogenous_shocks"] <- 1             # 0 for no, 1 for yes
  #forcings[1,"frequency_exogenous_shocks"] <- 7   # days between shocks (outside infections)
  #forcings[1,"new_infections_per_shock"] <- 15
  #forcings[1,"test_sensitivity"] <- 0.8
  #forcings[1,"test_specificity"] <- 0.98
  forcings[1,"test_cost"] <- 25                   # dollars
  forcings[1,"time_to_return_fps"] <- 1           # days to return a false positive to the Uninfected pool
  forcings[1,"confirmatory_test_cost"] <- 100     # dollars
  forcings[1,"population_size"] <- population     # size of total population
  #forcings[1,"building_size"] <- 500              # number of occupants of a representative building
  forcings[1,"wastewater_sampling"] <- 1          # 0 for no, 1 for yes

  # set up model configuration
  cycles.per.day <- 3
  n.days <- 100
  n.cycle <- cycles.per.day*n.days

  # sample parameters (need 2 data frames)
  parameter_names <- c("initial_infected", "Rt", "days_to_incubation", "time_to_recovery"
                       ,"frequency_exogenous_shocks", "new_infections_per_shock"
                       ,"test_sensitivity", "test_specificity"
                       ,"pct_advancing_to_symptoms", "symptom_case_fatality_ratio"
                       ,"frequency_of_screening"
                       ,"frequency_wastewater_samples", "frac_wastewater", "frequency_of_screening_wastewater", "threshold_wastewater", "lag_screening_wastewater"
                       ,"building_size"
                       ,"frac_noncompliance"
                       )
  n_parameters <- length(parameter_names)

  source("seir_model.R")
  source("sobol_wrapper.R")

  set.seed(2019)
  sample_lhs1 <- randomLHS(n_sample, n_parameters)
  sample_lhs2 <- randomLHS(n_sample, n_parameters)

  # scale the parameters
  source("scale_parameters.R")
  if (pri=="mid") {source("priors_mid.R")
  } else if (pri=="high") {source('priors_high.R')
  } else if (pri=="small") {source('priors_small.R')
  } else {print("ERROR: unrecognized priors.")}
  priors <- priors_func(parameter_names)

  parameters_lhs1 <- scale_parameters(sample_lhs1, parameter_names, priors)
  parameters_lhs2 <- scale_parameters(sample_lhs2, parameter_names, priors)

  if (second) {
    scheme <- "B"
    n_simulations_total <- n_sample*(2*n_parameters+2)
  } else {
    scheme <- "A"
    n_simulations_total <- n_sample*(n_parameters+2)
  }

##==============================================================================

  # run the actual sobol
  print(paste("Starting",n_simulations_total,"simulations total..."))
  t.out <- system.time(s.out <- sobolSalt(model=seir_sobol,
                               parameters_lhs1,
                               parameters_lhs2,
                               scheme=scheme,
                               nboot=n_bootstrap,
                               forcings=forcings, parameter_names=parameter_names,
                               forcing_names=forcing_names, cycles.per.day=cycles.per.day,
                               n.cycle=n.cycle, sens=sens))
  rownames(s.out$S) <- rownames(s.out$T) <- parameter_names
  print(paste("... done. Took",t.out[3]/60,"minutes."))

  # Check convergence - is maximum confidence interval width < 10% of the highest
  # total sensitivity index?
  max_sens_ind <- 0.1*max(s.out$T[,1])
  max_conf_int <- max(max(s.out$S[,5]-s.out$S[,4]), max(s.out$T[,5]-s.out$T[,4]))
  print(paste('max. conf int=',max_conf_int, ' want less than:  0.1 * max. sensitivity index=',max_sens_ind, sep=''))

##==============================================================================

  # first-order and total sensitivity index results
  headers.1st.tot <- matrix(c('Parameter', 'S1', 'S1_conf_low', 'S1_conf_high',
                              'ST', 'ST_conf_low', 'ST_conf_high'), nrow=1)
  output.1st.tot  <- data.frame(cbind( parameter_names, s.out$S[,c(1,4,5)], s.out$T[,c(1,4,5)]))
  colnames(output.1st.tot) <- headers.1st.tot
  s.out$output.1st.tot <- output.1st.tot

  # second-order sensitivities
  if ("S2" %in% names(s.out)) {
    S2.names <- NULL
    for (i in 1:(n_parameters-1)) {
      for (j in (i+1):n_parameters) {
        S2.names <- rbind(S2.names, c(parameter_names[i],parameter_names[j]))
      }
    }
    s.out$S2.names <- S2.names
    headers.2nd <- matrix(c('Parameter_1', 'Parameter_2', 'S2', 'S2_conf_low',
                            'S2_conf_high'), nrow=1)
    output.2nd <- data.frame(cbind( s.out$S2.names,s.out$S2[,c(1,4,5)] ))
    colnames(output.2nd) <- headers.2nd
    s.out$output.2nd <- output.2nd
  }

  today=Sys.Date(); today=format(today,format="%d%b%Y")
  filename.sobol <- paste('../output/sobol_',sens,'_',appen,today,'.rds', sep='')
  saveRDS(s.out, filename.sobol)

##==============================================================================

## Generate a radial convergence plot

  ## Libraries
  library(RColorBrewer) # good color palettes
  library(graphics)     # used when plotting polygons
  library(plotrix)      # used when plotting circles

  ## Functions in other files
  source('sobol_functions.R')

  # Set up
  sig.cutoff <- 0.01
  plotdir <- '../figures/'
  n_params <- ncol(s.out$X)

  # First-order and total sensitivity indices
  s1st <- s.out$output.1st.tot
  parnames.sobol <- s1st[,1]

  # Import second-order indices
  s2_table <- s.out$output.2nd

  # Convert second-order to upper-triangular matrix
  s2 <- matrix(nrow=n_params, ncol=n_params, byrow=FALSE)
  s2[1:(n_params-1), 2:n_params] = upper.diag(s2_table$S2)
  s2 <- as.data.frame(s2)
  colnames(s2) <- rownames(s2) <- parnames.sobol

  # Convert confidence intervals to upper-triangular matrix
  s2_conf_low <- matrix(nrow=n_params, ncol=n_params, byrow=FALSE)
  s2_conf_high <- matrix(nrow=n_params, ncol=n_params, byrow=FALSE)
  s2_conf_low[1:(n_params-1), 2:n_params] = upper.diag(s2_table$S2_conf_low)
  s2_conf_high[1:(n_params-1), 2:n_params] = upper.diag(s2_table$S2_conf_high)

  s2_conf_low <- as.data.frame(s2_conf_low)
  s2_conf_high <- as.data.frame(s2_conf_high)
  colnames(s2_conf_low) <- rownames(s2_conf_low) <- parnames.sobol
  colnames(s2_conf_high) <- rownames(s2_conf_high) <- parnames.sobol



  ##==============================================================================
  # Determine which indices are statistically significant

  # S1 & ST: confidence lower bound greater than sig.cutoff
  s1st1 <- stat_sig_s1st(s1st
                         ,method="congtr"
                         ,greater=sig.cutoff
                         ,sigCri='either')

  # S2: confidence lower bound greater than sig.cutoff
  # For Supplemental results for high-risk and small-college scenarios, using
  # just method='gtr'; can use more samples to narrow the CIs for those.
  if ((pri=="small") | (pri=="high")) {
    method <- "gtr"
  } ese {
    method <- "congtr"
  }

  s2_sig1 <- stat_sig_s2(s2
                         ,s2_conf_low
                         ,s2_conf_high
                         ,method=method
                         ,greater=sig.cutoff
  )

  ##==============================================================================

  # yields parnames_calib, which are the sensitive parameters
  ind_sensit <- which(s1st1$sig>0)
  ind_insens <- 1:n_params; ind_insens <- ind_insens[-ind_sensit]

  name_list1 <- list('Sensitive' = parnames.sobol[ind_sensit]
                     ,'Insensitive' = parnames.sobol[ind_insens]
  )

  # add Parameter symbols to plot

  if(TRUE) {
    # Use the nicer parameter symbols
#> parnames.sobol
# [1] "initial_infected"                  "Rt"                                "days_to_incubation"
# [4] "time_to_recovery"                  "frequency_exogenous_shocks"        "new_infections_per_shock"
# [7] "test_sensitivity"                  "test_specificity"                  "pct_advancing_to_symptoms"
#[10] "symptom_case_fatality_ratio"       "frequency_of_screening"            "frequency_wastewater_samples"
#[13] "frac_wastewater"                   "frequency_of_screening_wastewater" "threshold_wastewater"
#[16] "lag_screening_wastewater"          "building_size"                     "frac_noncompliance"
    name_symbols <- c('Initial\ninfected', expression('R'['t']), expression('T'['inc'])
                      , expression('T'['rec']), expression('T'['exo']), expression('N'['exo'])
                      , 'Se', 'Sp', expression('p'['symptoms'])
                      , expression('f'['fatal']), expression('T'['s']), expression('T'['ww'])
                      , expression('f'['ww']), expression('T'['s,ww']), 'W'
                      , expression('T'['lag']), expression('N'['building']), expression('f'['nc'])
                      )
    new_name_symbols <- c(name_symbols[ind_sensit], name_symbols[ind_insens])
  } else {
    new_name_symbols <- c(parnames.sobol[ind_sensit], parnames.sobol[ind_insens])
  }

  # defining list of colors for each group
  col_list1 <- list("Sensitive"    = "darkorange3"
                    ,"Insensitive" = "black"
  )

  # using function to assign variables and colors based on group
  s1st1 <- gp_name_col(name_list1
                       ,col_list1
                       ,s1st1)

  s1st1$symbols <- new_name_symbols


  # Things to change order of to group by sensitivity:
  # 1) s1st1
  # 2) s2
  # 3) s2_sig1

  s1st1_rearr <- rbind(s1st1[ind_sensit,], s1st1[ind_insens,])
  s1st1_rearr$symbols <- new_name_symbols # it was right already!
  s2_rearr <- rbind(s2[ind_sensit,],s2[ind_insens,]) # rearrange the rows...
  s2_rearr <- cbind(s2_rearr[,ind_sensit],s2_rearr[,ind_insens]) # ... and the columns
  s2_sig1_rearr <- rbind(s2_sig1[ind_sensit,], s2_sig1[ind_insens,])
  s2_sig1_rearr <- cbind(s2_sig1_rearr[,ind_sensit], s2_sig1_rearr[,ind_insens])

  ##==============================================================================
  ## Make the plot

  plot.filename <- paste(plotdir,'sobol_spider_',appen,today,sep='')

  plotRadCon(df=s1st1_rearr
             ,s2=s2_rearr
             ,plotS2=TRUE
             ,radSc = 2
             ,scaling=.47
             ,widthSc = 0.4
             ,s2_sig=s2_sig1_rearr
             ,filename = plot.filename
             ,plotType = 'EPS'
             ,gpNameMult=100
             ,varNameMult=1.4
             ,RingThick=0.14
             ,legLoc = "bottomcenter"
             ,cex = 1
             ,rt_names = 0
             ,s1_col = 'darkorange'
             ,st_col = 'gray30'
             ,line_col ='gray60'
             ,STthick = 0.4
             ,legFirLabs=c(.4,2.5), legTotLabs=c(.10,.50), legSecLabs=c(.01,.05)
             ,legend_egg=TRUE
  )

}

##==============================================================================
## End
##==============================================================================
