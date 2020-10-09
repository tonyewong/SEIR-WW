##==============================================================================
## analysis.R
##
## This script calculates some numbers to report in the manuscript, based on the
## output files from the sensitivity analysis.
## If you are running your own analysis, the time-stamp on the `Sobol_file` set
## below will be different.
##
## Questions?  Tony Wong (aewsma@rit.edu)
##==============================================================================

if (Sys.info()["user"]=="tony") {setwd("/Users/tony/codes/SEIR-WW/R")}
if (Sys.info()["user"]=="aewsma") {setwd("/Users/aewsma/codes/SEIR-WW/R")}


## Sensitivity analysis numbers
Sobol_file <- "../output/sobol_infections_sensinfections_ns10000_nb1000_secondTRUE_priorsmid_24Sep2020.rds"
s.out <- readRDS(Sobol_file)

print(paste("First-order sensitivity to Rt =",round(s.out$S["Rt",1],3)))
print(paste("Total sensitivity to Rt =",round(s.out$T["Rt",1],3)))

source("sobol_functions.R")
s1st <- s.out$output.1st.tot
s1st1 <- stat_sig_s1st(s1st
                       ,method="congtr"
                       ,greater=0.01
                       ,sigCri='either')
idx_sig <- which(s1st1[,"sig"]==1)
s1_sig <- round(s.out$S[idx_sig,1],2); names(s1_sig) <- rownames(s.out$S)[idx_sig]
print(rev(sort(s1_sig)))


##==============================================================================
## End
##==============================================================================
