# ==============================================================================
# Script to try and simulate the bilinear mode FM + lucky litter for 
# a more SHS like population
# Author: Luke Lloyd-Jones and Jack Easton
# Date started: 17/03/2025
# Date updated: 06/06/2025
# ==============================================================================

library(CKMRpop)
library(arules)
library(mvbutils)  
library(offarray)  
library(atease)    
library(plyr)
library(dplyr)
library(DNAage)
library(deconvodisc)
library(mgcv)
library(TMB)


# Compilation
# -----------

compile("../epi_aging/auto_calibration/m4_simulation/cpp/hsp_auto_tmb.cpp")
dyn.load(dynlib("../epi_aging/auto_calibration/m4_simulation/cpp/hsp_auto_tmb"))

model1 <- 'linear_fac'
out.rs <- out.rs <- paste0("autocal_simulation_", model1, "/")
system(paste0("ls ", out.rs, "*/processed* > ", out.rs, "/prcs_list.txt"))

lst <- read.table(paste0(out.rs, "/prcs_list.txt"))
nreps <- 50
for (i in seq(1, length(lst$V1)))
{
  file.i <- lst$V1[i]
  sub1 <- gsub(out.rs, "", file.i)
  sub2 <- gsub("/processed_.*", "", sub1)
  episd <- as.numeric(gsub(".*\\_", "", sub2))
  outsub <- gsub("\\/", "", sub1)
  
  # Make the results out file
  out.fl <- paste0(out.rs, sub2)
  
  # ----------------------------------
  # Grab the processed population 
  # ----------------------------------
  
  prcs <- get(load(lst$V1[i]))
  
  
  cmat        <- prcs$comp.mat
  cmat$deltaA <- cmat$ages_2 - cmat$ages_1
  cmat$deltaG <- cmat$epi_ages_2 - cmat$epi_ages_1
  cmat$deltaY <- cmat$samp_years_2-cmat$samp_years_1
  HSP_mat     <- as.matrix(cmat$HSP)
  HSP_size    <- dim(HSP_mat)
  
  sig_dA <- sd(cmat$deltaA)
  
  # Here I remove the intracohorts from cmat.
  # Not sure how realistic that is?
  # Need to threshold it
  
  #cmat <- cmat[cmat$B2 != cmat$B1, ]
  
  
  # =====================
  # With normal HSP model
  # =====================
  
  par <- c(log(1000), 0.01, log(0.1))
  fitto_norm <- nlminb(par,
                       NEG(lglickie_norm),
                       prcs = prcs,
                       fishSim = 1)
  params <- fitto_norm$par
  save(params,  file = paste0(out.fl, "/norm_par_", i %% nreps + 1, ".Rdata"))
  
  
  # ================================================
  # Vanilla with standard - not even a waffle - cone
  # ================================================
  
  dim(prcs$samples.2)[1]
  sum(prcs$comp.mat$HSP)
  
  data <- list(deltaY = as.array(cmat$samp_years_2 - cmat$samp_years_1),
               deltaG = as.array(cmat$epi_ages_2 - cmat$epi_ages_1),
               hspsYN = as.array(cmat$HSP))
  
  
  epi_nllk1 <- MakeADFun(data,
                         parameters = list(theta = params[1],
                                           psi   = params[3],
                                           beta  = 2,
                                           logsi = log(1.5)),
                         DLL = "hsp_auto_tmb")
  
  fit1 = nlminb(epi_nllk1$par,
                epi_nllk1$fn,
                epi_nllk1$gr,
                control = list(trace = 1, xf.tol = 1e-14, 
                               rel.tol = 2e-14, sing.tol = 2e-14))
  
  rep <- sdreport(epi_nllk1)
  summary(rep)
  out.1 <- list(fit1, rep)
  save(out.1,  file = paste0(out.fl, "/model_1_rep", i %% nreps + 1, ".Rdata"))
  
  
}


