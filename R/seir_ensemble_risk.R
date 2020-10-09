##==============================================================================
## seir_ensemble_risk.R
##
## Risk analysis using the of SEIR model of Paltiel et al 2020.
## We proceed by creating a large latin hypercube ensemble of simulations and
## seeing how those simulations perform in a nominal "control" scenario with
## Rt = 1.1, Nexo = 15 (new exposures per weekly shock), and fnc = 1% of the
## population is noncompliant with quarantine/isolation procedures. The
## `failure_cases` array creates the set of experiments where those assumptions
## are broken by increasing Rt to 1.5, increasing Nexo to 30, and increasing fnc
## to 10%. This is done singly, and with all combinations among these increases.
##
## For each failure case (combination of Rt, Nexo, and fnc parameters), we
## consider all combinations of T_s (frequency of traditional individual
## screenings), T_s,ww (frequency of wastewater surveillance once it is
## initiated), and T_lag (the lag before initiating wastewater-triggered
## surveillance testing). The combinations to consider are given by
## `screening_cases`.
##
## For each combination of failure case and screening case, we run an ensemble
## 3000 model simulations, where the other model parameters are drawn from their
## prior distributions. For each of these ensembles, we calculate the
## reliability of maintaining fewer than 100 infections across any 2-week
## period. This is done using the perceived infections (false positives, true
## positives, and symptomatic cases) and the actual infections (asymptomatic
## cases, true positives, and symptomatic cases).
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================

rm(list=ls())

##==============================================================================
## Set stuff for model configuration
##===================================

if (Sys.info()["user"]=="tony") {setwd("/Users/tony/codes/SEIR-WW/R")}
if (Sys.info()["user"]=="aewsma") {setwd("/Users/aewsma/codes/SEIR-WW/R")}

n.ensemble <- 3000 # how many simulations for each case
cycles.per.day <- 3
n.days <- 100
n.cycle <- cycles.per.day*n.days
N_day_running_total <- 14

failure_cases <- expand.grid(new_infections_per_shock=c(15, 30),
                             Rt=c(1.1, 1.5),
                             frac_noncompliance=c(0.01, 0.1))

screening_cases <- expand.grid(frequency_of_screening=c(999999),
                               frequency_of_screening_wastewater=c(1,2,3,4,5,6,7,8),
                               lag_screening_wastewater=c(1))

# set up input forcings (scenario-based, not sampled uncertainties)
# -- comment out and add to `parameter_names` below if you want any of these in the analysis as uncertainties
forcing_names <- c("exogenous_shocks"
                   ,"frequency_exogenous_shocks"
                   ,"test_cost", "time_to_return_fps"
                   ,"confirmatory_test_cost", "population_size"
                   ,"wastewater_sampling"
                   ,"frequency_of_screening", "frequency_of_screening_wastewater", "lag_screening_wastewater"
                   ,"new_infections_per_shock", "Rt", "frac_noncompliance"
                   )
forcings <- mat.or.vec(nr=1, nc=length(forcing_names))
colnames(forcings) <- forcing_names
forcings[1,"exogenous_shocks"] <- 1             # 0 for no, 1 for yes
forcings[1,"frequency_exogenous_shocks"] <- 7   # days between shocks (outside infections)
forcings[1,"test_cost"] <- 25                   # dollars
forcings[1,"time_to_return_fps"] <- 1           # days to return a false positive to the Uninfected pool
forcings[1,"confirmatory_test_cost"] <- 100     # dollars
forcings[1,"population_size"] <- 12500          # size of total population
forcings[1,"wastewater_sampling"] <- 1          # 0 for no, 1 for yes
#forcings[1,"frequency_of_screening_wastewater"] <- 2     # time (days) to screen `building_size` individuals after wastewater positive
#forcings[1,"lag_screening_wastewater"] <- 2     # time between wastewater sample return and screenings

# set up parameters
parameter_names <- c("initial_infected", "days_to_incubation", "time_to_recovery"#, "Rt"
                     #,"frequency_exogenous_shocks", "new_infections_per_shock"
                     ,"test_sensitivity", "test_specificity"
                     ,"pct_advancing_to_symptoms", "symptom_case_fatality_ratio"
                     #,"frequency_of_screening"
                     ,"frequency_wastewater_samples", "frac_wastewater", "threshold_wastewater"#, "lag_screening_wastewater"#, "frequency_of_screening_wastewater"
                     ,"building_size"
                     )
n_parameters <- length(parameter_names)
##==============================================================================



##==============================================================================
## Run the analysis
##==================

library(scales)
library(tidyverse)
library(lhs)
library(triangle)
library(extraDistr)

source("seir_model.R")
source("seir_wrapper.R")

# initialize output
model_out <- vector("list", length=nrow(failure_cases))
total_experiments <- nrow(failure_cases)*nrow(screening_cases)
cnt <- 0
today <- Sys.Date(); today <- format(today,format="%d%b%Y")
filename_output <- paste("../output/risk_results_",today,".RData", sep="")

for (i_failure in 1:nrow(failure_cases)) {
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
    print(paste("Starting case",cnt,"/",total_experiments))

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
    model_out[[i_failure]]$output[[i_screening]] <- c(length(which(max_perceived < 100))/length(max_perceived),
                                                      length(which(max_actual < 100))/length(max_actual),
                                                      median(all_infections))
  }
  # save relevant parameters and configuration information
  save(list=c("model_out","screening_cases", "failure_cases"), file=filename_output)
}
days <- (0:n.cycle)/cycles.per.day
##==============================================================================



##==============================================================================
## Analyze
##=========

#load("../output/risk_results_25Sep2020.RData")
load("../output/risk_results_05Oct2020.RData")
frequency_of_screening <- unique(screening_cases[,"frequency_of_screening"])
frequency_of_screening_wastewater <- unique(screening_cases[,"frequency_of_screening_wastewater"])
lag_screening_wastewater <- unique(screening_cases[,"lag_screening_wastewater"])
mat_reliability_per <- mat_reliability_act <- mat_infections <- vector("list", length=8)
for (k in 1:8) {
  mat_dims <- c(length(frequency_of_screening_wastewater), length(frequency_of_screening), length(lag_screening_wastewater))
  mat_dimnames <- vector("list", length=length(mat_dims))
  names(mat_dimnames) <- c("frequency_of_screening_wastewater", "frequency_of_screening", "lag_screening_wastewater")
  mat_dimnames$frequency_of_screening_wastewater <- frequency_of_screening_wastewater
  mat_dimnames$frequency_of_screening <- frequency_of_screening
  mat_dimnames$lag_screening_wastewater <- lag_screening_wastewater
  mat_reliability_per[[k]] <- mat_reliability_act[[k]] <- mat_infections[[k]] <- array(dim=mat_dims, dimnames=mat_dimnames)
  for (i in 1:length(frequency_of_screening_wastewater)) {
    for (j in 1:length(frequency_of_screening)) {
      for (l in 1:length(lag_screening_wastewater)) {
        tmp <- model_out[[k]]$output[[which(screening_cases[,"frequency_of_screening_wastewater"]==frequency_of_screening_wastewater[i] &
                                            screening_cases[,"frequency_of_screening"]==frequency_of_screening[j] &
                                            screening_cases[,"lag_screening_wastewater"]==lag_screening_wastewater[l])]]
        mat_reliability_per[[k]][i,j,l] <- tmp[1]
        mat_reliability_act[[k]][i,j,l] <- tmp[2]
        mat_infections[[k]][i,j,l] <- tmp[3]
      }
    }
  }
  # If want relative to the control case, normalize:
  #if (k > 1) {
  #    mat[[k]] <- mat[[k]] - mat[[1]]
  #}
}
save(list=c("model_out","screening_cases", "failure_cases","mat_reliability_per","mat_reliability_act","mat_infections"), file=filename_output)

## write output to CSV files for SOM
for (k in 1:8) {
    filename_mat <- paste("../output/infections_failurecase",k,".csv", sep="")
    write.csv(apply(t(round(mat_infections[[k]][,1,])), 2, rev), file=filename_mat)
}

##==============================================================================



##==============================================================================
## make a plot
##=============

library(RColorBrewer)
library(pracma)

# plotting details
n_splits <- 1000
splits <- linspace(0, 1, n_splits+1)

frequency_of_screening_wastewater <- as.numeric(dimnames(mat_reliability_act[[1]])[[1]]) # (more general than rownames(mat_reliability_act[[1]][,1,]))
midpoints_freq <- 0.5*(frequency_of_screening_wastewater[2:length(frequency_of_screening_wastewater)]+frequency_of_screening_wastewater[1:(length(frequency_of_screening_wastewater)-1)])

# pick a color palette
cols <- colorRampPalette(brewer.pal(8, "PiYG"))(n_splits)
cols <- colorRampPalette(c("firebrick","orangered","orange","steelblue"))(n_splits)
cols <- colorRampPalette(c("firebrick","orangered","orange","white"))(n_splits)


name_symbols <- c('Initial\ninfected', expression('R'['t']), expression('T'['inc'])
                  , expression('T'['rec']), expression('T'['exo']), expression('N'['exo'])
                  , 'Se', 'Sp', expression('p'['symtoms'])
                  , expression('f'['fatal']), expression('T'['s']), expression('T'['ww'])
                  , expression('f'['ww']), expression('T'['s,ww']), 'W'
                  , expression('T'['lag']), expression('N'['building']), expression('f'['nc'])
                  )

failure_labels <- c(expression("R"["t"]*"=1.1, "*"f"["nc"]*"=0.01, "*"N"["exo"]*"=15") #"Control",
                    ,expression("R"["t"]*"=1.1, "*"f"["nc"]*"=0.01, "*bold("N"["exo"]*"=30"))#"Higher shocks",
                    ,expression(bold("R"["t"]*"=1.5")*", "*"f"["nc"]*"=0.01, "*"N"["exo"]*"=15")#"Higher Rt",
                    ,expression(bold("R"["t"]*"=1.5")*", "*"f"["nc"]*"=0.01, "*bold("N"["exo"]*"=30"))#"Higher shocks & Rt",
                    ,expression("R"["t"]*"=1.1, "*bold("f"["nc"]*"=0.1")*", "*"N"["exo"]*"=15")#"Higher noncompliance",
                    ,expression("R"["t"]*"=1.1, "*bold("f"["nc"]*"=0.1")*", "*bold("N"["exo"]*"=30"))#"Higher shocks & noncomp.",
                    ,expression(bold("R"["t"]*"=1.5")*", "*bold("f"["nc"]*"=0.1")*", "*"N"["exo"]*"=15")#"Higher Rt & noncomp.",
                    ,expression(bold("R"["t"]*"=1.5")*", "*bold("f"["nc"]*"=0.1")*", "*bold("N"["exo"]*"=30"))#"Higher shocks, Rt & noncomp.")
                  )

panel_labels <- c(expression(bold('a')),
                  expression(bold('b')),
                  expression(bold('c')),
                  expression(bold('d')),
                  expression(bold('e')),
                  expression(bold('f')),
                  expression(bold('g')),
                  expression(bold('h')),
                  "")

## version of the figure with 8 panels (one per failure case) and bars representing
## the reliabilty under each of the wastewater-trigger screening time-scales.
## put the actual and detected/perceived reliabilities next to each other for
## each screening time case

# grab only the 1-day lag case
reliability_per <- reliability_act <- infections <- mat.or.vec(nr=nrow(failure_cases), nc=dim(mat_reliability_per[[1]])[1])
for (i_failure in 1:nrow(failure_cases)) {
  reliability_per[i_failure,] <- mat_reliability_per[[i_failure]][,1,1]
  reliability_act[i_failure,] <- mat_reliability_act[[i_failure]][,1,1]
  infections[i_failure,] <- mat_infections[[i_failure]][,1,1]
}



pdf('../figures/reliability.pdf',width=5.5,height=7, colormodel='cmyk', pointsize=11)
par(mfrow=c(4,2), mai=c(.55,.6,.25,.1), las=1)
dx <- 0.2 # width of boxes
dx_sep <- 0.02 # separation between boxes
idx_failure_cases_reorder <- c(1, 3,5,2, 7,4,6, 8)
cnt <- 0
for (i_failure in idx_failure_cases_reorder) {
    cnt <- cnt + 1
    plot(c(0,0),c(0,0), xlim=c(0.5,8.5), ylim=c(0,1), xlab="", ylab="", main="", xaxs='i', yaxs='i', xaxt="n", yaxt="n")
    grid()
    for (i_screening in 1:length(frequency_of_screening_wastewater)) {
        polygon(c(i_screening-dx,i_screening-dx_sep,i_screening-dx_sep,i_screening-dx),
                c(0,0,reliability_per[i_failure,i_screening],reliability_per[i_failure,i_screening]),
                col="gray30", border=NA)
        polygon(c(i_screening+dx,i_screening+dx_sep,i_screening+dx_sep,i_screening+dx),
                c(0,0,reliability_act[i_failure,i_screening],reliability_act[i_failure,i_screening]),
                col="darkorange", border=NA)
    }
    mtext(side=1, text="Time to screening results (days)", line=2.5, cex=0.9)
    mtext(side=2, text="Reliability", line=3, cex=0.9, las=0)
    mtext(side=3, text=failure_labels[i_failure], line=0.25, las=0, cex=0.95)
    mtext(side=3, text=panel_labels[cnt], line=0, cex=.9, adj=0)
    axis(1, at=1:length(frequency_of_screening_wastewater), labels=frequency_of_screening_wastewater, cex.axis=1.3)
    axis(2, at=seq(0,1,by=0.2), cex.axis=1.3)
    if (cnt==2) {legend(0.45,1.06, c("Using perceived infections only","Using actual infections"), col=c("gray30","darkorange"), pch=c(15,15), bty="n", cex=1)}
}
dev.off()





## these versions of the figures include the lag parameter, so they plot the
## full, gory matrix. that isn't really necessary, because the lag parameter
## doesn't really do anything.

lag_screening_wastewater <- dimnames(mat_reliability_act[[1]])[[3]] # (more general than rownames(mat_reliability_act[[1]][,1,]))
midpoints_lag <- 0.5*(lag_screening_wastewater[2:length(lag_screening_wastewater)]+lag_screening_wastewater[1:(length(lag_screening_wastewater)-1)])

# perceived reliability (based on FP+TP+S, the quarantine population)
pdf('../figures/reliability_matrix_perceived.pdf',width=7,height=7, colormodel='cmyk', pointsize=11)
par(mfrow=c(3,3), mai=c(.55,.6,.25,.1), las=1)
cnt <- 0
idx_failure_cases_reorder <- c(5,3,2,
                               7,6,4,
                               8,1,9)
for (i_failure in idx_failure_cases_reorder) {
    cnt <- cnt + 1
    if (i_failure==9) {
        #n_colbar <- 100
        #image(z=matrix(rev(linspace(0,1,n_colbar)), nrow=1, ncol=n_colbar))
        z=matrix(1:length(cols),nrow=1)
        x=1
        y=seq(0,1,len=length(cols))
        par(mai=c(.55,.6,.25,1.5), las=1)
        image(x,y,z,col=cols,axes=FALSE,xlab="",ylab="",xlim=c(0,5))
        axis(2)
        mtext(side=4, "reliability of < 100\ \ \ \ninfections over\nany 14-day period", las=1, adj=-0.1, padj=-0.5)
        par(mai=c(.55,.6,.25,.1), las=1)
    } else {
        image(z=mat_reliability_per[[i_failure]][,1,], x=frequency_of_screening_wastewater, y=lag_screening_wastewater, breaks=splits, col=cols, xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(0.5, 9))
        for (ii in midpoints_freq) {lines(c(ii,ii),c(0,100),lty=1,col="gray50",lwd=1)}
        for (ii in midpoints_lag) {lines(c(0,100),c(ii,ii),lty=1,col="gray50",lwd=1)}
        mtext(side=1, text="Time to screen (days)", line=2.5)
        mtext(side=2, text="Lag to screen (days)", line=2.7, las=0)
        mtext(side=3, text=failure_labels[i_failure], line=0.5, las=0)
        mtext(side=3, text=panel_labels[cnt], line=0.5, cex=1, adj=-0.38);
        axis(1, at=frequency_of_screening_wastewater, cex.axis=1.3)
        axis(2, at=lag_screening_wastewater, cex.axis=1.3)
        # highest reliability
        reliability_max <- max(mat_reliability_per[[i_failure]][,1,])
        idx <- which(mat_reliability_per[[i_failure]][,1,] == reliability_max, arr.ind=TRUE)
        lag_max <- lag_screening_wastewater[idx[1,"lag_screening_wastewater"]]
        freq_max <- frequency_of_screening_wastewater[idx[1,"frequency_of_screening_wastewater"]]
        #text(freq_max, lag_max, round(reliability_max,2))
        for (i_freq in 1:length(frequency_of_screening_wastewater)) {
            x_freq <- frequency_of_screening_wastewater[i_freq]
            for (i_lag in lag_screening_wastewater) {
                y_lag <- lag_screening_wastewater[i_lag]
                if (x_freq==2) {
                    text(x_freq+.25, y_lag, round(mat_reliability_per[[i_failure]][i_freq,1,i_lag],2), cex=0.9)
                } else {
                    text(x_freq, y_lag, round(mat_reliability_per[[i_failure]][i_freq,1,i_lag],2), cex=0.9)
                }
            }
        }
    }
}
dev.off()


# actual reliability (based on A+TP+S, the actual people who are infected)
pdf('../figures/reliability_matrix_actual.pdf',width=7,height=7, colormodel='cmyk', pointsize=11)
par(mfrow=c(3,3), mai=c(.55,.6,.25,.1), las=1)
cnt <- 0
for (i_failure in idx_failure_cases_reorder) {
    cnt <- cnt + 1
    if (i_failure==9) {
        #n_colbar <- 100
        #image(z=matrix(rev(linspace(0,1,n_colbar)), nrow=1, ncol=n_colbar))
        z=matrix(1:length(cols),nrow=1)
        x=1
        y=seq(0,1,len=length(cols))
        par(mai=c(.55,.6,.25,1.5), las=1)
        image(x,y,z,col=cols,axes=FALSE,xlab="",ylab="",xlim=c(0,5))
        axis(2)
        mtext(side=4, "reliability of < 100\ \ \ \ninfections over\nany 14-day period", las=1, adj=-0.1, padj=-0.5)
        par(mai=c(.55,.6,.25,.1), las=1)
    } else {
        image(z=mat_reliability_act[[i_failure]][,1,], x=frequency_of_screening_wastewater, y=lag_screening_wastewater, breaks=splits, col=cols, xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(0.5, 9))
        for (ii in midpoints_freq) {lines(c(ii,ii),c(0,100),lty=1,col="gray50",lwd=1)}
        for (ii in midpoints_lag) {lines(c(0,100),c(ii,ii),lty=1,col="gray50",lwd=1)}
        mtext(side=1, text="Time to screen (days)", line=2.5)
        mtext(side=2, text="Lag to screen (days)", line=2.7, las=0)
        mtext(side=3, text=failure_labels[i_failure], line=0.5, las=0)
        mtext(side=3, text=panel_labels[cnt], line=0.5, cex=1, adj=-0.38);
        axis(1, at=frequency_of_screening_wastewater, cex.axis=1.3)
        axis(2, at=lag_screening_wastewater, cex.axis=1.3)
        # highest reliability
        reliability_max <- max(mat_reliability_act[[i_failure]][,1,])
        idx <- which(mat_reliability_act[[i_failure]][,1,] == reliability_max, arr.ind=TRUE)
        lag_max <- lag_screening_wastewater[idx[1,"lag_screening_wastewater"]]
        freq_max <- frequency_of_screening_wastewater[idx[1,"frequency_of_screening_wastewater"]]
        #text(freq_max, lag_max, round(reliability_max,2))
        for (i_freq in 1:length(frequency_of_screening_wastewater)) {
            x_freq <- frequency_of_screening_wastewater[i_freq]
            for (i_lag in lag_screening_wastewater) {
                y_lag <- lag_screening_wastewater[i_lag]
                if (x_freq==2) {
                    text(x_freq+.25, y_lag, round(mat_reliability_act[[i_failure]][i_freq,1,i_lag],2), cex=0.9)
                } else {
                    text(x_freq, y_lag, round(mat_reliability_act[[i_failure]][i_freq,1,i_lag],2), cex=0.9)
                }
            }
        }
    }
}
dev.off()

##==============================================================================



##==============================================================================
## potentially useful code:
##==========================

if(FALSE) {

mat1 <- mat_reliability_per[[1]]; for (i in 2:nrow(failure_cases)) {mat1 <- cbind(mat1, mat_reliability_per[[i]])}
mat2 <- mat_reliability_act[[1]]; for (i in 2:nrow(failure_cases)) {mat2 <- cbind(mat2, mat_reliability_act[[i]])}
mat3 <- mat_infections[[1]]; for (i in 2:nrow(failure_cases)) {mat3 <- cbind(mat3, mat_infections[[i]])}

model_out <- readRDS("risk_results_07Sep2020.rds")
igood <- which(model_out[,"max_14day"] < 100)
ibad <- which(model_out[,"max_14day"] >= 100)

plot(model_out[ibad,"frequency_of_screening_wastewater"], model_out[ibad,"total_infections"], pch=16, col=rgb(.8,.1,.1,.25), type='p', ylim=c(0,1000))
points(model_out[igood,"frequency_of_screening_wastewater"], model_out[igood,"total_infections"], pch=16, col=rgb(.1,.1,.8,.25))
}

if(FALSE) {
  est1 <- est2 <- est3 <- rep(NA, n.ensemble)
  for (ii in 1:n.ensemble) {
    est1[ii] <- length(which(max_perceived[1:ii] <= 100))/ii
    est2[ii] <- length(which(max_actual[1:ii] <= 100))/ii
    est3[ii] <- median(all_infections[1:ii])
  }
  par(mfrow=c(3,1))
  plot(est1, type='l', main="perceived reliability inf < 100")
  plot(est2, type='l', main="actual reliability of inf < 100")
  plot(est3, type='l', main="actual total infections")
}
##==============================================================================



##==============================================================================
## End
##==============================================================================
