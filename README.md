# Validation-free-age-calibration
This repository contains R code and TMB (Template Model Builder) source files for fitting close-kin mark-recapture (CKMR) models that account for errors-in-variables. These models are designed to estimate demographic parameters, such as population size, from genetic kinship data while explicitly acknowledging and correcting for uncertainty in relatedness assignments.

---

## Key Features

* **TMB-based Estimation:** Utilizes TMB for fast and robust maximum likelihood estimation of model parameters.
* **Errors-in-Variables (EIV) Framework:** Explicitly models and corrects for potential errors in kinship assignments, leading to more accurate and less biased population estimates.
* **CKMR Analysis:** Provides a complete workflow for analyzing kinship data to infer demographic quantities.
* **Modular Code:** Organized R and C++ files for ease of understanding and modification.

---

## Requirements

The code is built and tested in R. To use this project, you will need:

* **R:** A recent version of the R statistical software.
* **R Packages:**
    * `TMB`: The core package for compiling and running the statistical models.
    * `Matrix`: Required by TMB.
    * `Rcpp`: For compiling the C++ code.
    * `data.table` or `dplyr`: For efficient data manipulation (recommended but not strictly required by the model).

You can install the required packages in R with the following command:

```R
install.packages(c("TMB", "Rcpp", "Matrix"))
Repository Structure
The repository is organized as follows:


.
├── R/                  # R functions for data preparation, model fitting, and plotting results.
├── src/                # TMB C++ source files (.cpp) containing the core model logic.
├── data/               # Placeholder for example data.
├── README.md           # This file.
├── LICENSE             # The repository license.
└── main.R              # An example script demonstrating how to run the model.
Getting Started
Follow these steps to set up the environment and run the model.

Clone the repository:

Bash

git clone [https://github.com/your-username/your-repository-name.git](https://github.com/your-username/your-repository-name.git)
cd your-repository-name
Compile the TMB Model:
Open R and navigate to the project directory. The TMB::compile function will create a dynamic library from the C++ source file.

R

library(TMB)
compile("src/eiv_ckmr.cpp") # Change filename if different
dyn.load(TMB::dynlib("src/eiv_ckmr"))
Run an Example Analysis:
The main.R script provides a full working example. It will load an example dataset, prepare it for the model, fit the model, and print a summary of the results.

R

source("R/model_functions.R") # Assumes a file with helper functions exists
source("main.R")

The Model

The statistical model is a TMB implementation of a capture-recapture likelihood, with the crucial addition of an errors-in-variables component. The model likelihood is based on the joint probability of observing a set of kinship pairs (k) given the true but unknown population size (N). The core of the model revolves around a likelihood function L that integrates over the uncertainty in the true relationship types.

$$L(N, \theta | \text{data}) = \sum_{r} P(\text{data} | r) \cdot P(r | N, \theta) 

$$where $r$ is the true relationship type (e.g., Parent-Offspring, Full-Sibling), $P(\\text{data} | r)$ is the probability of the observed data given a true relationship, and $P(r | N, \\theta)$ is the expected frequency of that relationship type in the population. This approach is particularly powerful because it uses the full genetic data to infer both the population size and the parameters of the EIV model, such as the probabilities of correctly identifying different kinship types. 

----- 

## Output and Interpretation Upon successful fitting, the model output will provide: 
* **Parameter Estimates:** Maximum likelihood estimates for population size ($N$) and other model parameters. 
* **Standard Errors:** Associated standard errors for all estimated parameters, derived from the Hessian matrix. 
* **AIC/BIC:** Information criteria for model comparison. 
* **Convergence Diagnostics:** Reports on whether the model converged successfully. The R functions in the `R/` folder will help you visualize these results, including confidence intervals and diagnostic plots. 
----- 



## License This project is licensed under the MIT License. See the `LICENSE` file for details. ``` ```$$