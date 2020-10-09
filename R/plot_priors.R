##==============================================================================
## plot_priors.R
##
## Plotting histograms of samples from the parameters' prior distributions. This
## is done for the control case (presented in the main text) and the large
## (high risk) university case and small-college case from the Supplemental
## Material.
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================

rm(list=ls())

if (Sys.info()["user"]=="tony") {setwd("/Users/tony/codes/SEIR-WW/R")}
if (Sys.info()["user"]=="aewsma") {setwd("/Users/aewsma/codes/SEIR-WW/R")}

parameter_names <- c("initial_infected", "Rt", "days_to_incubation", "time_to_recovery"
                     ,"frequency_exogenous_shocks", "new_infections_per_shock"
                     ,"test_sensitivity", "test_specificity"
                     ,"pct_advancing_to_symptoms", "symptom_case_fatality_ratio"
                     ,"frequency_of_screening"
                     ,"frequency_wastewater_samples", "frac_wastewater", "frequency_of_screening_wastewater", "threshold_wastewater", "lag_screening_wastewater"
                     ,"building_size"
                     ,"frac_noncompliance"
                     )

source("scale_parameters.R")

##==================
## control scenario
##==================
source("priors_mid.R")
priors <- priors_func(parameter_names)

n_samples <- 100000
n_parameters <- length(parameter_names)
U <- mat.or.vec(nr=n_samples, nc=n_parameters)
colnames(U) <- parameter_names

# generate large sample from Unif(0,1)...
for (pp in 1:n_parameters) {U[,pp] <- runif(n_samples)}

# then scale it to the prior distribution for parameter pp
X <- scale_parameters(U, parameter_names, priors)

# labels for the parameters
parameter_labels <- c("Initial infections", "Rt", "Incubation period (d)", "Recovery time (d)"
                      , "Time-scale for\nexogenous shocks (d)", "New infections/shock", "Test sensitivity", "Test specificity"
                      , "% of infections that\nadvance to symptoms", "Symptomatic case fatality ratio", "Time-scale for traditional\nscreenings (d)", "Time-scale for\nwastewater sampling (d)"
                      , "Fraction of individuals\ncontributing to wastewater", "Time-scale for\nwastewater screenings (d)", "Wastewater detection\nthreshold (count)", "Time lag for wastewater\nscreenings (d)"
                      , "Building size (count)", "Fraction of individuals\nwho are noncompliant"
                      )

pdf('../figures/priors_mid.pdf',width=7,height=9.5,colormodel='cmyk', pointsize=11)
par(mfrow=c(5,4), mai=c(.6,.42,.1,.1))
for (pp in 1:length(parameter_names)) {
  if (parameter_names[pp]=="frequency_of_screening_wastewater" | parameter_names[pp]=="frequency_wastewater_samples") {
    breaks <- seq(from=0.5, to=8.5, by=1)
    hist(X[,pp], main="", xlab="", ylab="", freq=FALSE, breaks=breaks)
  } else {
    hist(X[,pp], main="", xlab="", ylab="", freq=FALSE)
  }
  mtext(side=1, text=parameter_labels[pp], line=3.4, cex=0.85)
  mtext(side=2, text="Density", line=2.3, cex=0.85)
}
dev.off()

##===============
## high scenario
##===============
source("priors_high.R")
priors <- priors_func(parameter_names)

U <- mat.or.vec(nr=n_samples, nc=n_parameters)
colnames(U) <- parameter_names

# generate large sample from Unif(0,1)...
for (pp in 1:n_parameters) {U[,pp] <- runif(n_samples)}

# then scale it to the prior distribution for parameter pp
X <- scale_parameters(U, parameter_names, priors)

pdf('../figures/priors_high.pdf',width=7,height=9.5,colormodel='cmyk', pointsize=11)
par(mfrow=c(5,4), mai=c(.6,.42,.1,.1))
for (pp in 1:length(parameter_names)) {
  if (parameter_names[pp]=="frequency_of_screening_wastewater" | parameter_names[pp]=="frequency_wastewater_samples") {
    breaks <- seq(from=0.5, to=8.5, by=1)
    hist(X[,pp], main="", xlab="", ylab="", freq=FALSE, breaks=breaks)
  } else {
    hist(X[,pp], main="", xlab="", ylab="", freq=FALSE)
  }
  mtext(side=1, text=parameter_labels[pp], line=3.4, cex=0.85)
  mtext(side=2, text="Density", line=2.3, cex=0.85)
}
dev.off()

##================
## small scenario
##================
source("priors_small.R")
priors <- priors_func(parameter_names)

U <- mat.or.vec(nr=n_samples, nc=n_parameters)
colnames(U) <- parameter_names

# generate large sample from Unif(0,1)...
for (pp in 1:n_parameters) {U[,pp] <- runif(n_samples)}

# then scale it to the prior distribution for parameter pp
X <- scale_parameters(U, parameter_names, priors)

pdf('../figures/priors_small.pdf',width=7,height=9.5,colormodel='cmyk', pointsize=11)
par(mfrow=c(5,4), mai=c(.6,.42,.1,.1))
for (pp in 1:length(parameter_names)) {
  if (parameter_names[pp]=="frequency_of_screening_wastewater" | parameter_names[pp]=="frequency_wastewater_samples") {
    breaks <- seq(from=0.5, to=8.5, by=1)
    hist(X[,pp], main="", xlab="", ylab="", freq=FALSE, breaks=breaks)
  } else {
    hist(X[,pp], main="", xlab="", ylab="", freq=FALSE)
  }
  mtext(side=1, text=parameter_labels[pp], line=3.4, cex=0.85)
  mtext(side=2, text="Density", line=2.3, cex=0.85)
}
dev.off()


##==============================================================================
## End
##==============================================================================
