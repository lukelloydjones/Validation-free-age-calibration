# Validation-free-age-calibration

This repository contains R code and TMB (Template Model Builder) source files for fitting 
close-kin mark-recapture (CKMR) models that account for errors in chronological age 
measurements. Referred to as auto-calibration. Repository contains analytical scripts
used to generate the results for the manuscript "Validation-free estimation of chronological 
age via close-kin". The repo. also contains a test run to see how the method works,
used and summarised.

---

## Key Features

* **TMB-based Estimation:** Utilizes TMB for fast maximum likelihood estimation of model parameters 
and standard errors using algorithmic differentiation.
* **Aging uncertainty** Explicitly models mean and variance structures between noisy age and
chronological age leading to better estimates population demographic parameters.
* **CKMR incorporated:** Provides a workflow for analyzing CKMR data under a half-sibling 
pair model to infer demographic quantities given this uncertain age model.

## Model


---

## Requirements

The code is built and tested in R. To run the 'test', you will need the following installed:

* **R:** A recent version of the R statistical software.
* **R Packages:**
    * `TMB`: The core package for compiling and running the statistical models.
    * `mgcv`: For setting up smooth matrices.
    * `dplyr`: For efficient data manipulation (recommended but not strictly required by the model).
    * `Rfast`: Helps with some faster matrix sampling and calculations.
    * 'offarray' and 'mvbutils': From the MVB CKMR univeRse.
        options(repos = unique( c(
                 mvb = 'https://markbravington.r-universe.dev',
                 getOption( 'repos')[ 'CRAN']
                )))
        install.packages( "offarray")


## Directory contents

* **rscripts/**
    *

* **cpp/**
    *

* **test/**
    *