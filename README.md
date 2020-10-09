# README.md

## Purpose

This repository contains codes, input files, and output files for the SEIR modeling analysis describes in the manuscript [Evaluating the Sensitivity of SARS-CoV-2 Infection Rates on College Campuses to Wastewater Surveillance](preprint url here).

Questions? I'm happy to chat: Tony Wong (aewsma@rit.edu)

## Installation and compilation

This model is written in R, so no compilation is necessary.

To install, from here on Github, Clone or Download the repository.

Open an R terminal or RStudio. First, you must have the relevant packages installed, because the codes make heavy use of various R libraries for additional functionality. This can be done using the `install_packages.R` script in the `R` directory.

## Workflow

Here, we describe the workflow to reproduce the work described in the manuscript [Evaluating the Sensitivity of SARS-CoV-2 Infection Rates on College Campuses to Wastewater Surveillance](preprint url here).

### 1) Preliminaries

Note that many files have a `setwd()` command at the beginning to change the working directory to the `R` directory within this code repository. For example:
```
if (Sys.info()["user"]=="tony") {setwd("/Users/tony/Google Drive (aewsma@g.rit.edu)/research/covid/seir_ww/R")}
if (Sys.info()["user"]=="aewsma") {setwd("/Users/aewsma/Google Drive/research/covid/seir_ww/R")}
```
(The `if` statements are because I would run these codes on both my work and personal computer.) You will need to change this to match your own directory structure. If you aren't sure what this should look like, open R and navigate to the `R` directory from these codes. Then enter the command `getwd()` in R. Copy-paste the result into the quotes within the `setwd()` command. If you only have one computer that you're using because you have a healthier work-life balance than me, then you could get rid of the `if` statements and change the above to:
```
setwd([whatever your R directory is])
```

### 2) Making sure things work: conservation of people

The model compartments represent different states that individuals can be in at any given time. This includes various "sick" states (asymptomatic, true positive, and symptomatic), not-sick states (susceptible, false positive, and recovered), sort-of-sick state (exposed), and a state for those who have died. Through exposures, incubation, development of symptoms, recovery, and surveillance testing, the model has "fluxes" of people between these states. Akin to how an ecosystem model must conserve water and energy, the SEIR model must conserve people. That is, the overall number of individuals that you start the model with must equal the overall number at the end of a simulation.

As a step to make sure things work okay, it is a good idea to run the script `balance_check.R` before and after making any substantial modifications to the code. It will run a few model simulations with parameter values set from the input file `input/balance_check_set_parameters.csv`. Then this script will perform a check to make sure that the total census across all model compartments is equal to the prescribed overall campus population (which is an input forcing). The maximum imbalance for each simulation will be printed to the screen. If things go poorly, you'll lose an appreciable number of individuals. If things go well, you'll lose $O(10^{-11})$ individuals (just do to numerical accuracy).

### 3) Generating Figure 2 (comparison of surveillance strategies)

Run `seir_comparison.R`, which requires the input file `input/comparison_set_parameters.csv`. This CSV file sets the parameter values for the simulations that make up Figure 2 in the manuscript, which the `seir_comparison.R` script will automatically generate.

### 4) Generating Figure 3 (sensitivity analysis)

Run `seir_sensitivity_driver.R`. This will take a few hours on a modern desktop workstation. Figures 3, S6 and S7 will be generated.

### 5) Generating Figure 4 (risk analysis)

Run `seir_ensemble_risk.R`. This will take a few hours on a modern desktop workstation. Figure 4 will be generated. You can change the settings at the top of the script to include different combinations of wastewater/traditional screening times that you would like to test, but keep in mind that for each combination, there are 8 "failure" cases (scenarios for the parameters $R_t$, $N_{exo}$ and $f_{nc}$), and for each of those 8 scenarios, an ensemble of 3,000 simulations will be run.

### 6) Generating supplemental figures and/or analysis

`seir_ensemble_risk.R` with the `lag_screening_wastewater` parameter varying between 1 and 6 days.

To generate Figures S1-S3, which show the prior distributions from the sensitivity and risk analyses, run the script `plot_priors.R`.

To generate Figure S4, which plots the estimates of reliability as the number of simulations used increases, run the script `plot_sample_size.R`.

To obtain some of the numbers quoted in the manuscript, in particular, those from the sensitivity analysis, run the script `analysis.R`.

## Redistribution and use license

Please redistribute, please modify and please use this code. This code is part of SEIR-WW. SEIR-WW is free software: you can redistribute it and/or modify it under the terms of the MIT License (see `LICENSE.txt`). The original R version of the model was coded and shared by the University of Wisconsin-Madison Data Science Hub. Substantial modifications and additions have been made by Tony Wong at Rochester Institute of Technology. Therefore, both groups are listed on the copyright included here.

**MIT License**

Copyright (c) 2020 UW-Madison Data Science Hub
Copyright (c) 2020 Tony Wong

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

The software is provided "as is", without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose and noninfringement. In no event shall the authors or copyright holders be liable for any claim, damages or other liability, whether in an action of contract, tort or otherwise, arising from, out of or in connection with the software or the use or other dealings in the software.

## Questions?

There will always be questions. Whether you are trying to get the model to run, or wondering what a line of code does, we would love to talk to you. The corresponding author (Tony Wong) can be reached at aewsma@rit.edu.
