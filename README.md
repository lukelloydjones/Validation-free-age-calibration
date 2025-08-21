# Validation free age calibration

This repository contains R code and Template Model Builder (TMB) source files for fitting 
close-kin mark-recapture (CKMR) models that account for errors in chronological age 
measurements. Referred to as auto-calibration. Repository contains analytical scripts
used to generate the results for the manuscript "Validation-free estimation of chronological 
age via close-kin". The repo. also contains a test run to see how the method works and 
is summarised.

---

## Key Features

* **TMB-based Estimation:** Utilizes TMB for fast maximum likelihood estimation of model parameters 
and standard errors using algorithmic differentiation.
* **Aging uncertainty** Explicitly models mean and variance structures between noisy age and
chronological age leading to better estimates of population demographic parameters.
* **CKMR incorporated:** Provides a workflow for analyzing CKMR data under a half-sibling 
pair model to infer demographic quantities given this uncertain age model.

## Model


---

## Requirements

The code is built and tested recent version of the R statistical software. 
To run the 'test', you will need the following installed:

* **R Packages:**
    * `TMB`: The core package for compiling and running the statistical models.
    * `mgcv`: For setting up smooth matrices.
    * `dplyr`: For efficient data manipulation (recommended but not strictly required by the model).
    * `Rfast`: Helps with some faster matrix sampling and calculations.
    * `offarray` and `mvbutils`: From the MVB CKMR univeRse. 
       To install
        options(repos = unique( c(
                 mvb = 'https://markbravington.r-universe.dev',
                 getOption( 'repos')[ 'CRAN']
                )))
        install.packages( "offarray")


## Directory contents

* **rscripts/**
    * autocal/
        * Two scripts one for the normal model and another for the gamma. Runs the discrete
          deconvolution model for school shark for each of the models and produces summaries.
          Requires school shark data that is not stored publicly yet.
    * simulations_v2/
        * R scripts for each of 4 simulation scenarios. There is an 'all_simulation.R'
          script that produces all the simulated data using CKMRpop and then each of the
          model subdirectories performs the analyses and summaries. The suffix '_m4.R'
          runs the discrete deconvolution model on a high performance computer as it
          takes time to get through the 50*9 model runs.
* **cpp/**
    * Four TMB likelihoods are included here including the pairwise difference model, the
      normal linear model and the segmented regression models with normal and gamma error.

* **test/**
    * A small test run of the linear-only version of the model. The R script should be
      detailed enough to run through the method. 




