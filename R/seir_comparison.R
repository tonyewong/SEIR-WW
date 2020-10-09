##==============================================================================
## seir_comparison.R
##
## Run a set of simulations to compare control, moderate and high-risk scenarios
## using both wastewater-triggered screening tests and traditional individual
## screening tests.
## * Requires the input file `input/comparison_set_parameters.csv`
## * Generates Figure 2 from the Wong et al manuscript.
##
## If you want to run your own test cases, you can modify the parameter values
## in the input `comparison_set_parameters.csv` file. Then, run this script
## interactively in order to noodle around with the results.
##
## Here, no effective difference between "forcings" and "parameters". That only
## makes a difference in the sensitivity and risk analyses, where forcings are
## taken to be the same from simulation to simulation, whereas the parameters
## are taken to be uncertain and vary across simulations according to the
## distributions given in `priors_XXX.R`.
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
inputs_ensemble <-read.csv(file="../input/comparison_set_parameters.csv", header=FALSE)
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

# make a plot of the total known (or FP) infections at each time
plot_colors <- c("black","steelblue","plum4","firebrick","seagreen2","goldenrod","orchid2","coral")
panel_labels <- c(expression(bold('a')),
                  expression(bold('b')),
                  expression(bold('c')),
                  expression(bold('d')),
                  expression(bold('e')),
                  expression(bold('f')),
                  expression(bold('g')),
                  expression(bold('h')),
                  "")
days_extended <- (0:(n.cycle+30))/cycles.per.day
days_to_plot <- seq(from=10,to=n.days,by=10)
idx_to_plot <- match(days_to_plot, days)

pdf('../figures/comparison.pdf',width=6,height=7, colormodel='cmyk', pointsize=11)
par(mfrow=c(3,2), mai=c(.3,.8,.6,.02), las=0)
# plot of control simulations actual (A+S+TP) and perceived (FP+S+TP) infections
# --- PANEL A --- ctrl-WW perceived infections
mm <- 1
plot(days, model_out[[mm]][,"False Positive"]+model_out[[mm]][,"True Positive"]+model_out[[mm]][,"Symptoms"],
     type='l', lwd=1.5, lty=5, col=plot_colors[1], ylim=c(0,120), xlab="", ylab="", xaxs='i', yaxs='i', main="")
mtext(side=1, text="Day", line=2.2, cex=0.8)
mtext(side=2, text="Count", line=2.2, cex=0.8)
mtext(side=2, text="Control", line=4.6, cex=.9)
mtext(side=3, text=panel_labels[1], line=0, cex=1, adj=0);
mtext(side=3, text="Current actual and perceived\ninfection count", line=1.9, cex=.9);
# ctrl-WW actual infections
lines(days, model_out[[mm]][,"Asympt"]+model_out[[mm]][,"True Positive"]+model_out[[mm]][,"Symptoms"],
      type='l', lwd=1.5, lty=1, col=plot_colors[1])
# ctrl-trad perceived infections
lines(days, model_out[[mm+1]][,"False Positive"]+model_out[[mm+1]][,"True Positive"]+model_out[[mm+1]][,"Symptoms"],
      type='l', lwd=1.5, lty=5, col=plot_colors[5])
# ctrl-trad actual infections
lines(days, model_out[[mm+1]][,"Asympt"]+model_out[[mm+1]][,"True Positive"]+model_out[[mm+1]][,"Symptoms"],
      type='l', lwd=1.5, lty=1, col=plot_colors[5])
# add points for exogenous shocks
points(days[which(model_out[[mm]][,"Superspreader Event"]==1)],rep(1,sum(model_out[[mm]][,"Superspreader Event"])),pch=17,col="darkorange",cex=1.5)
text(65,10,"Exogenous infections")
legend(-1, 123, legend=c("Perceived infections, wastewater surv.","Actual infections, wastewater surv.","Perceived infections, traditional surv.","Actual infections, traditional surv."),
       lty=c(5,1,5,1), col=plot_colors[c(1,1,5,5)], lwd=1.5, bty="n")
# --- PANEL B --- plot of breakdown of isolation unit individuals
par(mai=c(.3,.32,.6,.5))
plot(days, model_out[[mm]][,"Asympt"]+model_out[[mm]][,"True Positive"]+model_out[[mm]][,"Symptoms"],
     type='l', lwd=1.5, lty=1, col="black", xlim=c(0,105), ylim=c(0,120), xlab="", ylab="", xaxs='i', yaxs='i', main="")
lines(days, model_out[[mm]][,"False Positive"]+model_out[[mm]][,"True Positive"]+model_out[[mm]][,"Symptoms"],
      type='l', lwd=1.5, lty=5, col="black")
mtext(side=1, text="Day", line=2.2, cex=0.8)
mtext(side=3, text=panel_labels[2], line=0, cex=1, adj=0);
mtext(side=3, text="Breakdown of actual and\nperceived infections", line=1.9, cex=.9);
for (ii in 1:length(days_to_plot)) {
    idx <- idx_to_plot[ii]
    ## ACTUAL infections
    # Symptomatic
    polygon(x=days_extended[c(idx_to_plot[ii]-6,idx_to_plot[ii],idx_to_plot[ii],idx_to_plot[ii]-6)],
            y=c(model_out[[mm]][idx,"Symptoms"],model_out[[mm]][idx,"Symptoms"],0,0),
            col=plot_colors[2], border=NA)
    # True positive
    polygon(x=days_extended[c(idx_to_plot[ii]-6,idx_to_plot[ii],idx_to_plot[ii],idx_to_plot[ii]-6)],
            y=c(model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"],model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"],
                model_out[[mm]][idx,"Symptoms"],model_out[[mm]][idx,"Symptoms"]),
            col=plot_colors[6], border=NA)
    # Asymptomatic
    polygon(x=days_extended[c(idx_to_plot[ii]-6,idx_to_plot[ii],idx_to_plot[ii],idx_to_plot[ii]-6)],
            y=c(model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"]+model_out[[mm]][idx,"Asympt"],
                model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"]+model_out[[mm]][idx,"Asympt"],
                model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"],
                model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"]),
            col=plot_colors[4], border=NA)
    ## PERCEIVED infections
    # Symptomatic
    polygon(x=days_extended[c(idx_to_plot[ii]+7,idx_to_plot[ii]+1,idx_to_plot[ii]+1,idx_to_plot[ii]+7)],
            y=c(model_out[[mm]][idx,"Symptoms"],model_out[[mm]][idx,"Symptoms"],0,0),
            col=plot_colors[2], border=NA)
    # True positive
    polygon(x=days_extended[c(idx_to_plot[ii]+7,idx_to_plot[ii]+1,idx_to_plot[ii]+1,idx_to_plot[ii]+7)],
            y=c(model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"],model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"],
                model_out[[mm]][idx,"Symptoms"],model_out[[mm]][idx,"Symptoms"]),
            col=plot_colors[6], border=NA)
    # Asymptomatic
    polygon(x=days_extended[c(idx_to_plot[ii]+7,idx_to_plot[ii]+1,idx_to_plot[ii]+1,idx_to_plot[ii]+7)],
            y=c(model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"]+model_out[[mm]][idx,"False Positive"],
                model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"]+model_out[[mm]][idx,"False Positive"],
                model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"],
                model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"]),
            col=plot_colors[7], border=NA)
}
legend(-1, 123, legend=c("Symptomatic","True positive","Asymptomatic","False positive"),
       pch=c(15,15,15,15), col=plot_colors[c(2,6,4,7)], bty="n")
# add arrows to denote actual and perceived infections
text(60,90,"Actual"); arrows(x0=65, y0=85.5, x1=68.7, y1=62.5, length=0.05, code=2)
text(81,83,"Perceived"); arrows(x0=75, y0=78, x1=71.5, y1=42.5, length=0.05, code=2)
# same, but for moderate scenario
mm <- 3
# --- PANEL C --- ctrl-WW perceived infections
par(mai=c(.4,.8,.5,.02))
plot(days, model_out[[mm]][,"False Positive"]+model_out[[mm]][,"True Positive"]+model_out[[mm]][,"Symptoms"],
     type='l', lwd=1.5, lty=5, col=plot_colors[1], ylim=c(0,180), xlab="", ylab="", xaxs='i', yaxs='i', main="")
mtext(side=1, text="Day", line=2.2, cex=0.8)
mtext(side=2, text="Count", line=2.2, cex=0.8)
mtext(side=2, text="Moderate", line=4.6, cex=.9)
mtext(side=3, text=panel_labels[3], line=0, cex=1, adj=0);
# ctrl-WW actual infections
lines(days, model_out[[mm]][,"Asympt"]+model_out[[mm]][,"True Positive"]+model_out[[mm]][,"Symptoms"],
      type='l', lwd=1.5, lty=1, col=plot_colors[1])
# ctrl-trad perceived infections
lines(days, model_out[[mm+1]][,"False Positive"]+model_out[[mm+1]][,"True Positive"]+model_out[[mm+1]][,"Symptoms"],
      type='l', lwd=1.5, lty=5, col=plot_colors[5])
# ctrl-trad actual infections
lines(days, model_out[[mm+1]][,"Asympt"]+model_out[[mm+1]][,"True Positive"]+model_out[[mm+1]][,"Symptoms"],
      type='l', lwd=1.5, lty=1, col=plot_colors[5])
# add points for exogenous shocks
points(days[which(model_out[[mm]][,"Superspreader Event"]==1)],rep(1.5,sum(model_out[[mm]][,"Superspreader Event"])),pch=17,col="darkorange",cex=1.5)
text(65,15,"Exogenous infections")
#legend(0, 112, legend=c("perceived infections with wastewater surv.","Actual infections with wastewater surv.","perceived infections with traditional surv.","Actual infections with traditional surv."),
#       lty=c(1,2,1,2), col=plot_colors[c(2,2,5,5)], lwd=1.5, bty="n")
# --- PANEL D --- plot of breakdown of isolation unit individuals
par(mai=c(.4,.32,.5,.5))
plot(days, model_out[[mm]][,"Asympt"]+model_out[[mm]][,"True Positive"]+model_out[[mm]][,"Symptoms"],
     type='l', lwd=1.5, lty=1, col="black", xlim=c(0,105), ylim=c(0,180), xlab="", ylab="", xaxs='i', yaxs='i', main="")
lines(days, model_out[[mm]][,"False Positive"]+model_out[[mm]][,"True Positive"]+model_out[[mm]][,"Symptoms"],
      type='l', lwd=1.5, lty=5, col="black")
mtext(side=1, text="Day", line=2.2, cex=0.8)
mtext(side=3, text=panel_labels[4], line=0, cex=1, adj=0);
for (ii in 1:length(days_to_plot)) {
    idx <- idx_to_plot[ii]
    ## ACTUAL infections
    # Symptomatic
    polygon(x=days_extended[c(idx_to_plot[ii]-6,idx_to_plot[ii],idx_to_plot[ii],idx_to_plot[ii]-6)],
            y=c(model_out[[mm]][idx,"Symptoms"],model_out[[mm]][idx,"Symptoms"],0,0),
            col=plot_colors[2], border=NA)
    # True positive
    polygon(x=days_extended[c(idx_to_plot[ii]-6,idx_to_plot[ii],idx_to_plot[ii],idx_to_plot[ii]-6)],
            y=c(model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"],model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"],
                model_out[[mm]][idx,"Symptoms"],model_out[[mm]][idx,"Symptoms"]),
            col=plot_colors[6], border=NA)
    # Asymptomatic
    polygon(x=days_extended[c(idx_to_plot[ii]-6,idx_to_plot[ii],idx_to_plot[ii],idx_to_plot[ii]-6)],
            y=c(model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"]+model_out[[mm]][idx,"Asympt"],
                model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"]+model_out[[mm]][idx,"Asympt"],
                model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"],
                model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"]),
            col=plot_colors[4], border=NA)
    ## PERCEIVED infections
    # Symptomatic
    polygon(x=days_extended[c(idx_to_plot[ii]+7,idx_to_plot[ii]+1,idx_to_plot[ii]+1,idx_to_plot[ii]+7)],
            y=c(model_out[[mm]][idx,"Symptoms"],model_out[[mm]][idx,"Symptoms"],0,0),
            col=plot_colors[2], border=NA)
    # True positive
    polygon(x=days_extended[c(idx_to_plot[ii]+7,idx_to_plot[ii]+1,idx_to_plot[ii]+1,idx_to_plot[ii]+7)],
            y=c(model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"],model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"],
                model_out[[mm]][idx,"Symptoms"],model_out[[mm]][idx,"Symptoms"]),
            col=plot_colors[6], border=NA)
    # Asymptomatic
    polygon(x=days_extended[c(idx_to_plot[ii]+7,idx_to_plot[ii]+1,idx_to_plot[ii]+1,idx_to_plot[ii]+7)],
            y=c(model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"]+model_out[[mm]][idx,"False Positive"],
                model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"]+model_out[[mm]][idx,"False Positive"],
                model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"],
                model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"]),
            col=plot_colors[7], border=NA)
}
# same, but for severe scenario
mm <- 5
# --- PANEL E --- ctrl-WW perceived infections
par(mai=c(.5,.8,.4,.02))
plot(days, model_out[[mm]][,"False Positive"]+model_out[[mm]][,"True Positive"]+model_out[[mm]][,"Symptoms"],
     type='l', lwd=1.5, lty=5, col=plot_colors[1], ylim=c(0,900), xlab="", ylab="", xaxs='i', yaxs='i', main="")
mtext(side=1, text="Day", line=2.2, cex=0.8)
mtext(side=2, text="Count", line=2.2, cex=0.8)
mtext(side=2, text="Severe", line=4.6, cex=0.9)
mtext(side=3, text=panel_labels[5], line=0, cex=1, adj=0);
# ctrl-WW actual infections
lines(days, model_out[[mm]][,"Asympt"]+model_out[[mm]][,"True Positive"]+model_out[[mm]][,"Symptoms"],
      type='l', lwd=1.5, lty=1, col=plot_colors[1])
# ctrl-trad perceived infections
lines(days, model_out[[mm+1]][,"False Positive"]+model_out[[mm+1]][,"True Positive"]+model_out[[mm+1]][,"Symptoms"],
      type='l', lwd=1.5, lty=5, col=plot_colors[5])
# ctrl-trad actual infections
lines(days, model_out[[mm+1]][,"Asympt"]+model_out[[mm+1]][,"True Positive"]+model_out[[mm+1]][,"Symptoms"],
      type='l', lwd=1.5, lty=1, col=plot_colors[5])
# add points for exogenous shocks
points(days[which(model_out[[mm]][,"Superspreader Event"]==1)],rep(7.5,sum(model_out[[mm]][,"Superspreader Event"])),pch=17,col="darkorange",cex=1.5)
text(65,75,"Exogenous infections")
#legend(0, 112, legend=c("perceived infections with wastewater surv.","Actual infections with wastewater surv.","perceived infections with traditional surv.","Actual infections with traditional surv."),
#       lty=c(1,2,1,2), col=plot_colors[c(2,2,5,5)], lwd=1.5, bty="n")
# --- PANEL F --- plot of breakdown of isolation unit individuals
par(mai=c(.5,.32,.4,.5))
plot(days, model_out[[mm]][,"Asympt"]+model_out[[mm]][,"True Positive"]+model_out[[mm]][,"Symptoms"],
     type='l', lwd=1.5, lty=1, col="black", xlim=c(0,105), ylim=c(0,900), xlab="", ylab="", xaxs='i', yaxs='i', main="")
lines(days, model_out[[mm]][,"False Positive"]+model_out[[mm]][,"True Positive"]+model_out[[mm]][,"Symptoms"],
      type='l', lwd=1.5, lty=5, col="black")
mtext(side=1, text="Day", line=2.2, cex=0.8)
mtext(side=3, text=panel_labels[6], line=0, cex=1, adj=0);
for (ii in 1:length(days_to_plot)) {
    idx <- idx_to_plot[ii]
    ## ACTUAL infections
    # Symptomatic
    polygon(x=days_extended[c(idx_to_plot[ii]-6,idx_to_plot[ii],idx_to_plot[ii],idx_to_plot[ii]-6)],
            y=c(model_out[[mm]][idx,"Symptoms"],model_out[[mm]][idx,"Symptoms"],0,0),
            col=plot_colors[2], border=NA)
    # True positive
    polygon(x=days_extended[c(idx_to_plot[ii]-6,idx_to_plot[ii],idx_to_plot[ii],idx_to_plot[ii]-6)],
            y=c(model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"],model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"],
                model_out[[mm]][idx,"Symptoms"],model_out[[mm]][idx,"Symptoms"]),
            col=plot_colors[6], border=NA)
    # Asymptomatic
    polygon(x=days_extended[c(idx_to_plot[ii]-6,idx_to_plot[ii],idx_to_plot[ii],idx_to_plot[ii]-6)],
            y=c(model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"]+model_out[[mm]][idx,"Asympt"],
                model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"]+model_out[[mm]][idx,"Asympt"],
                model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"],
                model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"]),
            col=plot_colors[4], border=NA)
    ## PERCEIVED infections
    # Symptomatic
    polygon(x=days_extended[c(idx_to_plot[ii]+7,idx_to_plot[ii]+1,idx_to_plot[ii]+1,idx_to_plot[ii]+7)],
            y=c(model_out[[mm]][idx,"Symptoms"],model_out[[mm]][idx,"Symptoms"],0,0),
            col=plot_colors[2], border=NA)
    # True positive
    polygon(x=days_extended[c(idx_to_plot[ii]+7,idx_to_plot[ii]+1,idx_to_plot[ii]+1,idx_to_plot[ii]+7)],
            y=c(model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"],model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"],
                model_out[[mm]][idx,"Symptoms"],model_out[[mm]][idx,"Symptoms"]),
            col=plot_colors[6], border=NA)
    # Asymptomatic
    polygon(x=days_extended[c(idx_to_plot[ii]+7,idx_to_plot[ii]+1,idx_to_plot[ii]+1,idx_to_plot[ii]+7)],
            y=c(model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"]+model_out[[mm]][idx,"False Positive"],
                model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"]+model_out[[mm]][idx,"False Positive"],
                model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"],
                model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"]),
            col=plot_colors[7], border=NA)
}

dev.off()


##==============================================================================


days_to_plot <- seq(from=20,to=n.days,by=20)
idx_to_plot <- match(days_to_plot, days)

pdf('../figures/comparison_both.pdf',width=6,height=7, colormodel='cmyk', pointsize=11)
par(mfrow=c(3,2), mai=c(.3,.8,.6,.02), las=0)
# plot of control simulations actual (A+S+TP) and perceived (FP+S+TP) infections
# --- PANEL A --- ctrl-WW perceived infections
mm <- 1
plot(days, model_out[[mm]][,"False Positive"]+model_out[[mm]][,"True Positive"]+model_out[[mm]][,"Symptoms"],
     type='l', lwd=1.5, lty=5, col=plot_colors[1], ylim=c(0,120), xlab="", ylab="", xaxs='i', yaxs='i', main="")
mtext(side=1, text="Day", line=2.2, cex=0.8)
mtext(side=2, text="Count", line=2.2, cex=0.8)
mtext(side=2, text="Control", line=4.6, cex=.9)
mtext(side=3, text=panel_labels[1], line=0, cex=1, adj=0);
mtext(side=3, text="Current actual and perceived\ninfection count", line=1.9, cex=.9);
# ctrl-WW actual infections
lines(days, model_out[[mm]][,"Asympt"]+model_out[[mm]][,"True Positive"]+model_out[[mm]][,"Symptoms"],
      type='l', lwd=1.5, lty=1, col=plot_colors[1])
# ctrl-trad perceived infections
lines(days, model_out[[mm+1]][,"False Positive"]+model_out[[mm+1]][,"True Positive"]+model_out[[mm+1]][,"Symptoms"],
      type='l', lwd=1.5, lty=5, col=plot_colors[5])
# ctrl-trad actual infections
lines(days, model_out[[mm+1]][,"Asympt"]+model_out[[mm+1]][,"True Positive"]+model_out[[mm+1]][,"Symptoms"],
      type='l', lwd=1.5, lty=1, col=plot_colors[5])
# add points for exogenous shocks
points(days[which(model_out[[mm]][,"Superspreader Event"]==1)],rep(1,sum(model_out[[mm]][,"Superspreader Event"])),pch=17,col="darkorange",cex=1.5)
text(65,10,"Exogenous infections")
legend(-1, 123, legend=c("Perceived infections, wastewater surv.","Actual infections, wastewater surv.","Perceived infections, traditional surv.","Actual infections, traditional surv."),
       lty=c(5,1,5,1), col=plot_colors[c(1,1,5,5)], lwd=1.5, bty="n")
# --- PANEL B --- plot of breakdown of isolation unit individuals
par(mai=c(.3,.32,.6,.5))
plot(days, model_out[[mm]][,"Asympt"]+model_out[[mm]][,"True Positive"]+model_out[[mm]][,"Symptoms"],
     type='l', lwd=1.5, lty=1, col="black", xlim=c(0,107), ylim=c(0,120), xlab="", ylab="", xaxs='i', yaxs='i', main="")
lines(days, model_out[[mm]][,"False Positive"]+model_out[[mm]][,"True Positive"]+model_out[[mm]][,"Symptoms"],
      type='l', lwd=1.5, lty=5, col="black")
lines(days, model_out[[mm+1]][,"Asympt"]+model_out[[mm+1]][,"True Positive"]+model_out[[mm+1]][,"Symptoms"],
      type='l', lwd=1.5, lty=1, col=plot_colors[5])
lines(days, model_out[[mm+1]][,"False Positive"]+model_out[[mm+1]][,"True Positive"]+model_out[[mm+1]][,"Symptoms"],
      type='l', lwd=1.5, lty=5, col=plot_colors[5])
mtext(side=1, text="Day", line=2.2, cex=0.8)
mtext(side=3, text=panel_labels[2], line=0, cex=1, adj=0);
mtext(side=3, text="Breakdown of actual and\nperceived infections", line=1.9, cex=.9);
for (ii in 1:length(days_to_plot)) {
    idx <- idx_to_plot[ii]
    ### WASTEWATER SURVEILLANCE
    ## ACTUAL infections
    # Symptomatic
    polygon(x=days_extended[c(idx_to_plot[ii]+1,idx_to_plot[ii]+7,idx_to_plot[ii]+7,idx_to_plot[ii]+1)],
            y=c(model_out[[mm]][idx,"Symptoms"],model_out[[mm]][idx,"Symptoms"],0,0),
            col=plot_colors[2], border=NA)
    # True positive
    polygon(x=days_extended[c(idx_to_plot[ii]+1,idx_to_plot[ii]+7,idx_to_plot[ii]+7,idx_to_plot[ii]+1)],
            y=c(model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"],model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"],
                model_out[[mm]][idx,"Symptoms"],model_out[[mm]][idx,"Symptoms"]),
            col=plot_colors[6], border=NA)
    # Asymptomatic
    polygon(x=days_extended[c(idx_to_plot[ii]+1,idx_to_plot[ii]+7,idx_to_plot[ii]+7,idx_to_plot[ii]+1)],
            y=c(model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"]+model_out[[mm]][idx,"Asympt"],
                model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"]+model_out[[mm]][idx,"Asympt"],
                model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"],
                model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"]),
            col=plot_colors[4], border=NA)
    ## PERCEIVED infections
    # Symptomatic
    polygon(x=days_extended[c(idx_to_plot[ii]+9,idx_to_plot[ii]+15,idx_to_plot[ii]+15,idx_to_plot[ii]+9)],
            y=c(model_out[[mm]][idx,"Symptoms"],model_out[[mm]][idx,"Symptoms"],0,0),
            col=plot_colors[2], border=NA)
    # True positive
    polygon(x=days_extended[c(idx_to_plot[ii]+9,idx_to_plot[ii]+15,idx_to_plot[ii]+15,idx_to_plot[ii]+9)],
            y=c(model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"],model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"],
                model_out[[mm]][idx,"Symptoms"],model_out[[mm]][idx,"Symptoms"]),
            col=plot_colors[6], border=NA)
    # Asymptomatic
    polygon(x=days_extended[c(idx_to_plot[ii]+9,idx_to_plot[ii]+15,idx_to_plot[ii]+15,idx_to_plot[ii]+9)],
            y=c(model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"]+model_out[[mm]][idx,"False Positive"],
                model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"]+model_out[[mm]][idx,"False Positive"],
                model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"],
                model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"]),
            col=plot_colors[7], border=NA)
    #
    ### WASTEWATER SURVEILLANCE
    ## ACTUAL infections
    # Symptomatic
    polygon(x=days_extended[c(idx_to_plot[ii]-7,idx_to_plot[ii]-1,idx_to_plot[ii]-1,idx_to_plot[ii]-7)],
            y=c(model_out[[mm+1]][idx,"Symptoms"],model_out[[mm+1]][idx,"Symptoms"],0,0),
            col=plot_colors[2], border=NA)
    # True positive
    polygon(x=days_extended[c(idx_to_plot[ii]-7,idx_to_plot[ii]-1,idx_to_plot[ii]-1,idx_to_plot[ii]-7)],
            y=c(model_out[[mm+1]][idx,"Symptoms"]+model_out[[mm+1]][idx,"True Positive"],model_out[[mm+1]][idx,"Symptoms"]+model_out[[mm+1]][idx,"True Positive"],
                model_out[[mm+1]][idx,"Symptoms"],model_out[[mm+1]][idx,"Symptoms"]),
            col=plot_colors[6], border=NA)
    # Asymptomatic
    polygon(x=days_extended[c(idx_to_plot[ii]-7,idx_to_plot[ii]-1,idx_to_plot[ii]-1,idx_to_plot[ii]-7)],
            y=c(model_out[[mm+1]][idx,"Symptoms"]+model_out[[mm+1]][idx,"True Positive"]+model_out[[mm+1]][idx,"Asympt"],
                model_out[[mm+1]][idx,"Symptoms"]+model_out[[mm+1]][idx,"True Positive"]+model_out[[mm+1]][idx,"Asympt"],
                model_out[[mm+1]][idx,"Symptoms"]+model_out[[mm+1]][idx,"True Positive"],
                model_out[[mm+1]][idx,"Symptoms"]+model_out[[mm+1]][idx,"True Positive"]),
            col=plot_colors[4], border=NA)
    ## PERCEIVED infections
    # Symptomatic
    polygon(x=days_extended[c(idx_to_plot[ii]-15,idx_to_plot[ii]-9,idx_to_plot[ii]-9,idx_to_plot[ii]-15)],
            y=c(model_out[[mm+1]][idx,"Symptoms"],model_out[[mm+1]][idx,"Symptoms"],0,0),
            col=plot_colors[2], border=NA)
    # True positive
    polygon(x=days_extended[c(idx_to_plot[ii]-15,idx_to_plot[ii]-9,idx_to_plot[ii]-9,idx_to_plot[ii]-15)],
            y=c(model_out[[mm+1]][idx,"Symptoms"]+model_out[[mm+1]][idx,"True Positive"],model_out[[mm+1]][idx,"Symptoms"]+model_out[[mm+1]][idx,"True Positive"],
                model_out[[mm+1]][idx,"Symptoms"],model_out[[mm+1]][idx,"Symptoms"]),
            col=plot_colors[6], border=NA)
    # Asymptomatic
    polygon(x=days_extended[c(idx_to_plot[ii]-15,idx_to_plot[ii]-9,idx_to_plot[ii]-9,idx_to_plot[ii]-15)],
            y=c(model_out[[mm+1]][idx,"Symptoms"]+model_out[[mm+1]][idx,"True Positive"]+model_out[[mm+1]][idx,"False Positive"],
                model_out[[mm+1]][idx,"Symptoms"]+model_out[[mm+1]][idx,"True Positive"]+model_out[[mm+1]][idx,"False Positive"],
                model_out[[mm+1]][idx,"Symptoms"]+model_out[[mm+1]][idx,"True Positive"],
                model_out[[mm+1]][idx,"Symptoms"]+model_out[[mm+1]][idx,"True Positive"]),
            col=plot_colors[7], border=NA)
}
legend(-1, 123, legend=c("Symptomatic","True positive","Asymptomatic","False positive"),
       pch=c(15,15,15,15), col=plot_colors[c(2,6,4,7)], bty="n")
# add arrows to denote actual and perceived infections
text(84,109,"Actual"); arrows(x0=83.4, y0=105, x1=78.5, y1=61.2, length=0.04, code=2); arrows(x0=83.4, y0=105, x1=81.2, y1=62.5, length=0.04, code=2)
text(95,93,"Perceived"); arrows(x0=88, y0=89, x1=84, y1=43.5, length=0.04, code=2); arrows(x0=88, y0=89, x1=77, y1=74, length=0.04, code=2)
text(55, 99, "Trad."); arrows(x0=57.4, y0=95, x1=57.4, y1=72, length=0.05, angle=90, code=2)
text(65, 87, "WW"); arrows(x0=62.6, y0=83, x1=62.6, y1=63, length=0.05, angle=90, code=2)# same, but for moderate scenario
mm <- 3
# --- PANEL C --- ctrl-WW perceived infections
par(mai=c(.4,.8,.5,.02))
plot(days, model_out[[mm]][,"False Positive"]+model_out[[mm]][,"True Positive"]+model_out[[mm]][,"Symptoms"],
     type='l', lwd=1.5, lty=5, col=plot_colors[1], ylim=c(0,180), xlab="", ylab="", xaxs='i', yaxs='i', main="")
mtext(side=1, text="Day", line=2.2, cex=0.8)
mtext(side=2, text="Count", line=2.2, cex=0.8)
mtext(side=2, text="Moderate", line=4.6, cex=.9)
mtext(side=3, text=panel_labels[3], line=0, cex=1, adj=0);
# ctrl-WW actual infections
lines(days, model_out[[mm]][,"Asympt"]+model_out[[mm]][,"True Positive"]+model_out[[mm]][,"Symptoms"],
      type='l', lwd=1.5, lty=1, col=plot_colors[1])
# ctrl-trad perceived infections
lines(days, model_out[[mm+1]][,"False Positive"]+model_out[[mm+1]][,"True Positive"]+model_out[[mm+1]][,"Symptoms"],
      type='l', lwd=1.5, lty=5, col=plot_colors[5])
# ctrl-trad actual infections
lines(days, model_out[[mm+1]][,"Asympt"]+model_out[[mm+1]][,"True Positive"]+model_out[[mm+1]][,"Symptoms"],
      type='l', lwd=1.5, lty=1, col=plot_colors[5])
# add points for exogenous shocks
points(days[which(model_out[[mm]][,"Superspreader Event"]==1)],rep(1.5,sum(model_out[[mm]][,"Superspreader Event"])),pch=17,col="darkorange",cex=1.5)
text(65,15,"Exogenous infections")
#legend(0, 112, legend=c("perceived infections with wastewater surv.","Actual infections with wastewater surv.","perceived infections with traditional surv.","Actual infections with traditional surv."),
#       lty=c(1,2,1,2), col=plot_colors[c(2,2,5,5)], lwd=1.5, bty="n")
# --- PANEL D --- plot of breakdown of isolation unit individuals
par(mai=c(.4,.32,.5,.5))
plot(days, model_out[[mm]][,"Asympt"]+model_out[[mm]][,"True Positive"]+model_out[[mm]][,"Symptoms"],
     type='l', lwd=1.5, lty=1, col="black", xlim=c(0,107), ylim=c(0,180), xlab="", ylab="", xaxs='i', yaxs='i', main="")
lines(days, model_out[[mm]][,"False Positive"]+model_out[[mm]][,"True Positive"]+model_out[[mm]][,"Symptoms"],
      type='l', lwd=1.5, lty=5, col="black")
lines(days, model_out[[mm+1]][,"Asympt"]+model_out[[mm+1]][,"True Positive"]+model_out[[mm+1]][,"Symptoms"],
      type='l', lwd=1.5, lty=1, col=plot_colors[5])
lines(days, model_out[[mm+1]][,"False Positive"]+model_out[[mm+1]][,"True Positive"]+model_out[[mm+1]][,"Symptoms"],
      type='l', lwd=1.5, lty=5, col=plot_colors[5])
mtext(side=1, text="Day", line=2.2, cex=0.8)
mtext(side=3, text=panel_labels[4], line=0, cex=1, adj=0);
for (ii in 1:length(days_to_plot)) {
    idx <- idx_to_plot[ii]
    ### WASTEWATER SURVEILLANCE
    ## ACTUAL infections
    # Symptomatic
    polygon(x=days_extended[c(idx_to_plot[ii]+1,idx_to_plot[ii]+7,idx_to_plot[ii]+7,idx_to_plot[ii]+1)],
            y=c(model_out[[mm]][idx,"Symptoms"],model_out[[mm]][idx,"Symptoms"],0,0),
            col=plot_colors[2], border=NA)
    # True positive
    polygon(x=days_extended[c(idx_to_plot[ii]+1,idx_to_plot[ii]+7,idx_to_plot[ii]+7,idx_to_plot[ii]+1)],
            y=c(model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"],model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"],
                model_out[[mm]][idx,"Symptoms"],model_out[[mm]][idx,"Symptoms"]),
            col=plot_colors[6], border=NA)
    # Asymptomatic
    polygon(x=days_extended[c(idx_to_plot[ii]+1,idx_to_plot[ii]+7,idx_to_plot[ii]+7,idx_to_plot[ii]+1)],
            y=c(model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"]+model_out[[mm]][idx,"Asympt"],
                model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"]+model_out[[mm]][idx,"Asympt"],
                model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"],
                model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"]),
            col=plot_colors[4], border=NA)
    ## PERCEIVED infections
    # Symptomatic
    polygon(x=days_extended[c(idx_to_plot[ii]+9,idx_to_plot[ii]+15,idx_to_plot[ii]+15,idx_to_plot[ii]+9)],
            y=c(model_out[[mm]][idx,"Symptoms"],model_out[[mm]][idx,"Symptoms"],0,0),
            col=plot_colors[2], border=NA)
    # True positive
    polygon(x=days_extended[c(idx_to_plot[ii]+9,idx_to_plot[ii]+15,idx_to_plot[ii]+15,idx_to_plot[ii]+9)],
            y=c(model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"],model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"],
                model_out[[mm]][idx,"Symptoms"],model_out[[mm]][idx,"Symptoms"]),
            col=plot_colors[6], border=NA)
    # Asymptomatic
    polygon(x=days_extended[c(idx_to_plot[ii]+9,idx_to_plot[ii]+15,idx_to_plot[ii]+15,idx_to_plot[ii]+9)],
            y=c(model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"]+model_out[[mm]][idx,"False Positive"],
                model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"]+model_out[[mm]][idx,"False Positive"],
                model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"],
                model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"]),
            col=plot_colors[7], border=NA)
    #
    ### WASTEWATER SURVEILLANCE
    ## ACTUAL infections
    # Symptomatic
    polygon(x=days_extended[c(idx_to_plot[ii]-7,idx_to_plot[ii]-1,idx_to_plot[ii]-1,idx_to_plot[ii]-7)],
            y=c(model_out[[mm+1]][idx,"Symptoms"],model_out[[mm+1]][idx,"Symptoms"],0,0),
            col=plot_colors[2], border=NA)
    # True positive
    polygon(x=days_extended[c(idx_to_plot[ii]-7,idx_to_plot[ii]-1,idx_to_plot[ii]-1,idx_to_plot[ii]-7)],
            y=c(model_out[[mm+1]][idx,"Symptoms"]+model_out[[mm+1]][idx,"True Positive"],model_out[[mm+1]][idx,"Symptoms"]+model_out[[mm+1]][idx,"True Positive"],
                model_out[[mm+1]][idx,"Symptoms"],model_out[[mm+1]][idx,"Symptoms"]),
            col=plot_colors[6], border=NA)
    # Asymptomatic
    polygon(x=days_extended[c(idx_to_plot[ii]-7,idx_to_plot[ii]-1,idx_to_plot[ii]-1,idx_to_plot[ii]-7)],
            y=c(model_out[[mm+1]][idx,"Symptoms"]+model_out[[mm+1]][idx,"True Positive"]+model_out[[mm+1]][idx,"Asympt"],
                model_out[[mm+1]][idx,"Symptoms"]+model_out[[mm+1]][idx,"True Positive"]+model_out[[mm+1]][idx,"Asympt"],
                model_out[[mm+1]][idx,"Symptoms"]+model_out[[mm+1]][idx,"True Positive"],
                model_out[[mm+1]][idx,"Symptoms"]+model_out[[mm+1]][idx,"True Positive"]),
            col=plot_colors[4], border=NA)
    ## PERCEIVED infections
    # Symptomatic
    polygon(x=days_extended[c(idx_to_plot[ii]-15,idx_to_plot[ii]-9,idx_to_plot[ii]-9,idx_to_plot[ii]-15)],
            y=c(model_out[[mm+1]][idx,"Symptoms"],model_out[[mm+1]][idx,"Symptoms"],0,0),
            col=plot_colors[2], border=NA)
    # True positive
    polygon(x=days_extended[c(idx_to_plot[ii]-15,idx_to_plot[ii]-9,idx_to_plot[ii]-9,idx_to_plot[ii]-15)],
            y=c(model_out[[mm+1]][idx,"Symptoms"]+model_out[[mm+1]][idx,"True Positive"],model_out[[mm+1]][idx,"Symptoms"]+model_out[[mm+1]][idx,"True Positive"],
                model_out[[mm+1]][idx,"Symptoms"],model_out[[mm+1]][idx,"Symptoms"]),
            col=plot_colors[6], border=NA)
    # Asymptomatic
    polygon(x=days_extended[c(idx_to_plot[ii]-15,idx_to_plot[ii]-9,idx_to_plot[ii]-9,idx_to_plot[ii]-15)],
            y=c(model_out[[mm+1]][idx,"Symptoms"]+model_out[[mm+1]][idx,"True Positive"]+model_out[[mm+1]][idx,"False Positive"],
                model_out[[mm+1]][idx,"Symptoms"]+model_out[[mm+1]][idx,"True Positive"]+model_out[[mm+1]][idx,"False Positive"],
                model_out[[mm+1]][idx,"Symptoms"]+model_out[[mm+1]][idx,"True Positive"],
                model_out[[mm+1]][idx,"Symptoms"]+model_out[[mm+1]][idx,"True Positive"]),
            col=plot_colors[7], border=NA)
}
# same, but for severe scenario
mm <- 5
# --- PANEL E --- ctrl-WW perceived infections
par(mai=c(.5,.8,.4,.02))
plot(days, model_out[[mm]][,"False Positive"]+model_out[[mm]][,"True Positive"]+model_out[[mm]][,"Symptoms"],
     type='l', lwd=1.5, lty=5, col=plot_colors[1], ylim=c(0,900), xlab="", ylab="", xaxs='i', yaxs='i', main="")
mtext(side=1, text="Day", line=2.2, cex=0.8)
mtext(side=2, text="Count", line=2.2, cex=0.8)
mtext(side=2, text="Severe", line=4.6, cex=0.9)
mtext(side=3, text=panel_labels[5], line=0, cex=1, adj=0);
# ctrl-WW actual infections
lines(days, model_out[[mm]][,"Asympt"]+model_out[[mm]][,"True Positive"]+model_out[[mm]][,"Symptoms"],
      type='l', lwd=1.5, lty=1, col=plot_colors[1])
# ctrl-trad perceived infections
lines(days, model_out[[mm+1]][,"False Positive"]+model_out[[mm+1]][,"True Positive"]+model_out[[mm+1]][,"Symptoms"],
      type='l', lwd=1.5, lty=5, col=plot_colors[5])
# ctrl-trad actual infections
lines(days, model_out[[mm+1]][,"Asympt"]+model_out[[mm+1]][,"True Positive"]+model_out[[mm+1]][,"Symptoms"],
      type='l', lwd=1.5, lty=1, col=plot_colors[5])
# add points for exogenous shocks
points(days[which(model_out[[mm]][,"Superspreader Event"]==1)],rep(7.5,sum(model_out[[mm]][,"Superspreader Event"])),pch=17,col="darkorange",cex=1.5)
text(65,75,"Exogenous infections")
#legend(0, 112, legend=c("perceived infections with wastewater surv.","Actual infections with wastewater surv.","perceived infections with traditional surv.","Actual infections with traditional surv."),
#       lty=c(1,2,1,2), col=plot_colors[c(2,2,5,5)], lwd=1.5, bty="n")
# --- PANEL F --- plot of breakdown of isolation unit individuals
par(mai=c(.5,.32,.4,.5))
plot(days, model_out[[mm]][,"Asympt"]+model_out[[mm]][,"True Positive"]+model_out[[mm]][,"Symptoms"],
     type='l', lwd=1.5, lty=1, col="black", xlim=c(0,107), ylim=c(0,900), xlab="", ylab="", xaxs='i', yaxs='i', main="")
lines(days, model_out[[mm]][,"False Positive"]+model_out[[mm]][,"True Positive"]+model_out[[mm]][,"Symptoms"],
      type='l', lwd=1.5, lty=5, col="black")
lines(days, model_out[[mm+1]][,"Asympt"]+model_out[[mm+1]][,"True Positive"]+model_out[[mm+1]][,"Symptoms"],
      type='l', lwd=1.5, lty=1, col=plot_colors[5])
lines(days, model_out[[mm+1]][,"False Positive"]+model_out[[mm+1]][,"True Positive"]+model_out[[mm+1]][,"Symptoms"],
      type='l', lwd=1.5, lty=5, col=plot_colors[5])
mtext(side=1, text="Day", line=2.2, cex=0.8)
mtext(side=3, text=panel_labels[6], line=0, cex=1, adj=0);
for (ii in 1:length(days_to_plot)) {
    idx <- idx_to_plot[ii]
    ### WASTEWATER SURVEILLANCE
    ## ACTUAL infections
    # Symptomatic
    polygon(x=days_extended[c(idx_to_plot[ii]+1,idx_to_plot[ii]+7,idx_to_plot[ii]+7,idx_to_plot[ii]+1)],
            y=c(model_out[[mm]][idx,"Symptoms"],model_out[[mm]][idx,"Symptoms"],0,0),
            col=plot_colors[2], border=NA)
    # True positive
    polygon(x=days_extended[c(idx_to_plot[ii]+1,idx_to_plot[ii]+7,idx_to_plot[ii]+7,idx_to_plot[ii]+1)],
            y=c(model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"],model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"],
                model_out[[mm]][idx,"Symptoms"],model_out[[mm]][idx,"Symptoms"]),
            col=plot_colors[6], border=NA)
    # Asymptomatic
    polygon(x=days_extended[c(idx_to_plot[ii]+1,idx_to_plot[ii]+7,idx_to_plot[ii]+7,idx_to_plot[ii]+1)],
            y=c(model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"]+model_out[[mm]][idx,"Asympt"],
                model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"]+model_out[[mm]][idx,"Asympt"],
                model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"],
                model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"]),
            col=plot_colors[4], border=NA)
    ## PERCEIVED infections
    # Symptomatic
    polygon(x=days_extended[c(idx_to_plot[ii]+9,idx_to_plot[ii]+15,idx_to_plot[ii]+15,idx_to_plot[ii]+9)],
            y=c(model_out[[mm]][idx,"Symptoms"],model_out[[mm]][idx,"Symptoms"],0,0),
            col=plot_colors[2], border=NA)
    # True positive
    polygon(x=days_extended[c(idx_to_plot[ii]+9,idx_to_plot[ii]+15,idx_to_plot[ii]+15,idx_to_plot[ii]+9)],
            y=c(model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"],model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"],
                model_out[[mm]][idx,"Symptoms"],model_out[[mm]][idx,"Symptoms"]),
            col=plot_colors[6], border=NA)
    # Asymptomatic
    polygon(x=days_extended[c(idx_to_plot[ii]+9,idx_to_plot[ii]+15,idx_to_plot[ii]+15,idx_to_plot[ii]+9)],
            y=c(model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"]+model_out[[mm]][idx,"False Positive"],
                model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"]+model_out[[mm]][idx,"False Positive"],
                model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"],
                model_out[[mm]][idx,"Symptoms"]+model_out[[mm]][idx,"True Positive"]),
            col=plot_colors[7], border=NA)
    #
    ### WASTEWATER SURVEILLANCE
    ## ACTUAL infections
    # Symptomatic
    polygon(x=days_extended[c(idx_to_plot[ii]-7,idx_to_plot[ii]-1,idx_to_plot[ii]-1,idx_to_plot[ii]-7)],
            y=c(model_out[[mm+1]][idx,"Symptoms"],model_out[[mm+1]][idx,"Symptoms"],0,0),
            col=plot_colors[2], border=NA)
    # True positive
    polygon(x=days_extended[c(idx_to_plot[ii]-7,idx_to_plot[ii]-1,idx_to_plot[ii]-1,idx_to_plot[ii]-7)],
            y=c(model_out[[mm+1]][idx,"Symptoms"]+model_out[[mm+1]][idx,"True Positive"],model_out[[mm+1]][idx,"Symptoms"]+model_out[[mm+1]][idx,"True Positive"],
                model_out[[mm+1]][idx,"Symptoms"],model_out[[mm+1]][idx,"Symptoms"]),
            col=plot_colors[6], border=NA)
    # Asymptomatic
    polygon(x=days_extended[c(idx_to_plot[ii]-7,idx_to_plot[ii]-1,idx_to_plot[ii]-1,idx_to_plot[ii]-7)],
            y=c(model_out[[mm+1]][idx,"Symptoms"]+model_out[[mm+1]][idx,"True Positive"]+model_out[[mm+1]][idx,"Asympt"],
                model_out[[mm+1]][idx,"Symptoms"]+model_out[[mm+1]][idx,"True Positive"]+model_out[[mm+1]][idx,"Asympt"],
                model_out[[mm+1]][idx,"Symptoms"]+model_out[[mm+1]][idx,"True Positive"],
                model_out[[mm+1]][idx,"Symptoms"]+model_out[[mm+1]][idx,"True Positive"]),
            col=plot_colors[4], border=NA)
    ## PERCEIVED infections
    # Symptomatic
    polygon(x=days_extended[c(idx_to_plot[ii]-15,idx_to_plot[ii]-9,idx_to_plot[ii]-9,idx_to_plot[ii]-15)],
            y=c(model_out[[mm+1]][idx,"Symptoms"],model_out[[mm+1]][idx,"Symptoms"],0,0),
            col=plot_colors[2], border=NA)
    # True positive
    polygon(x=days_extended[c(idx_to_plot[ii]-15,idx_to_plot[ii]-9,idx_to_plot[ii]-9,idx_to_plot[ii]-15)],
            y=c(model_out[[mm+1]][idx,"Symptoms"]+model_out[[mm+1]][idx,"True Positive"],model_out[[mm+1]][idx,"Symptoms"]+model_out[[mm+1]][idx,"True Positive"],
                model_out[[mm+1]][idx,"Symptoms"],model_out[[mm+1]][idx,"Symptoms"]),
            col=plot_colors[6], border=NA)
    # Asymptomatic
    polygon(x=days_extended[c(idx_to_plot[ii]-15,idx_to_plot[ii]-9,idx_to_plot[ii]-9,idx_to_plot[ii]-15)],
            y=c(model_out[[mm+1]][idx,"Symptoms"]+model_out[[mm+1]][idx,"True Positive"]+model_out[[mm+1]][idx,"False Positive"],
                model_out[[mm+1]][idx,"Symptoms"]+model_out[[mm+1]][idx,"True Positive"]+model_out[[mm+1]][idx,"False Positive"],
                model_out[[mm+1]][idx,"Symptoms"]+model_out[[mm+1]][idx,"True Positive"],
                model_out[[mm+1]][idx,"Symptoms"]+model_out[[mm+1]][idx,"True Positive"]),
            col=plot_colors[7], border=NA)
}

dev.off()


##==============================================================================


if (FALSE) {
# plot of quarantine population size
ii <- 1
plot(days, model_out[[ii]][,"False Positive"]+model_out[[ii]][,"True Positive"]+model_out[[ii]][,"Symptoms"],
     type='l', lwd=1.5, col=plot_colors[ii], ylim=c(0,running_max_total*1.1), xlab="", ylab="")
mtext(side=1, text="Day", line=2)
mtext(side=2, text="Quarantine (count)", line=2.2)
lines(days, model_out[[ii]][,"Wastewater Screening"]*50, col=plot_colors[ii], lty=3)
if (n.ensemble > 1) {
  for (ii in 2:n.ensemble) {lines(days, model_out[[ii]][,"False Positive"]+model_out[[ii]][,"True Positive"]+model_out[[ii]][,"Symptoms"],
                                  lwd=1.5, col=plot_colors[ii]); lines(days, model_out[[ii]][,"Wastewater Screening"]*50, col=plot_colors[ii], lty=3)}
}

# plot of quarantine population breakdown
ii <- 1
plot(days, model_out[[ii]][,"False Positive"]/(model_out[[ii]][,"False Positive"]+model_out[[ii]][,"True Positive"]+model_out[[ii]][,"Symptoms"]),
     type='l', lwd=1.5, col=plot_colors[ii], ylim=c(0,1), xlab="", ylab="")
mtext(side=1, text="Day", line=2)
mtext(side=2, text="False positives in quarantine (fraction)", line=2.2)
lines(days, model_out[[ii]][,"Wastewater Screening"]*50, col=plot_colors[ii], lty=3)
if (n.ensemble > 1) {
  for (ii in 2:n.ensemble) {lines(days, model_out[[ii]][,"False Positive"]/(model_out[[ii]][,"False Positive"]+model_out[[ii]][,"True Positive"]+model_out[[ii]][,"Symptoms"]),
                                  lwd=1.5, col=plot_colors[ii]); lines(days, model_out[[ii]][,"Wastewater Screening"]*50, col=plot_colors[ii], lty=3)}
}

# plot of Cumulative Infections for all
# ctrl-WW and ctrl-trad
plot(days, model_out[[1]][,"Cumulative Infections"], type='l', lwd=1.5, col=plot_colors[2],
     ylim=c(0,4000), xlab="", ylab="", xaxs='i', yaxs='i')
lines(days, model_out[[2]][,"Cumulative Infections"], lwd=1.5, col=plot_colors[2], lty=2)
# mid-WW and mid-trad
lines(days, model_out[[3]][,"Cumulative Infections"], lwd=1.5, col=plot_colors[3])
lines(days, model_out[[4]][,"Cumulative Infections"], lwd=1.5, col=plot_colors[3], lty=2)
# high-WW and high-trad
lines(days, model_out[[5]][,"Cumulative Infections"], lwd=1.5, col=plot_colors[4])
lines(days, model_out[[6]][,"Cumulative Infections"], lwd=1.5, col=plot_colors[4], lty=2)
mtext(side=1, text="Day", line=2.2)
mtext(side=2, text="Total infections (count)", line=2.2)
mtext(side=3, text=panel_labels[1], line=0, cex=1, adj=0);
legend(1, 4000, legend=c("Control","Moderate","Severe","Wastewater","Traditional"),
       lty=c(NA,NA,NA,1,2), pch=c(15,15,15,NA,NA), col=plot_colors[c(2,3,4,1,1)], lwd=1.5, bty="n")
grid()

}

# save relevant parameters and configuration information
#save(model_out, parameters_ensemble, parameter_names, days, file = "../output/ensemble_results.RData")

##==============================================================================
## End
##==============================================================================
