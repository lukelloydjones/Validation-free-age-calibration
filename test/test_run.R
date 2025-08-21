# ==============================================================================
# Script to run auto-calibration from a data set simulated using CKMRpop
#   - Use a gamma GLM model to simulate the noisy age|true age relationship
#   - Fit a normal HSP CKMR model
#   - Fit a pairwise difference model
#   - Fit normal discrete deconvolution model
#   - Fit gamma discrete deconvolution model
#   - all assumes you are in the /test folder
# Author: Luke Lloyd-Jones
# Date started: 21/08/2025
# Date updated: 21/08/2025
# ==============================================================================

library(mvbutils)
library(offarray)
library(dplyr)
library(mgcv)
library(TMB)
library(Rfast)

source("../rscripts/autocal/test_helpers.R")

# Compile the TMB models
# ----------------------

# Base pair-wise difference model

compile("../cpp/hsp_auto_tmb.cpp")
dyn.load(dynlib("../cpp/hsp_auto_tmb"))

# Normal discrete deconvolution model

compile("../cpp/epi_auto_deconv_nllk_lcky_known_age_v2.cpp")
dyn.load(dynlib("../cpp/epi_auto_deconv_nllk_lcky_known_age_v2"))


# Load the pre-fab data
# ---------------------

prcs <- get(load("data/0.0045_0.5/processed_13.Rdata"))
pop  <- get(load("data/0.0045_0.5/pop_13.Rdata"))

cmat        <- prcs$comp.mat
cmat$deltaA <- cmat$ages_2 - cmat$ages_1
cmat$deltaG <- cmat$epi_ages_2 - cmat$epi_ages_1
cmat$deltaY <- cmat$samp_years_2-cmat$samp_years_1
HSP_mat     <- as.matrix(cmat$HSP)
HSP_size    <- dim(HSP_mat)
sig.dA      <- sd(cmat$deltaA)


# ---------------------
# With normal HSP model
# ---------------------

par <- c(log(1000), 0.01, log(0.1))
fitto_norm <- nlminb(par,
                     NEG(lglickie_norm),
                     prcs = prcs,
                     fishSim = 1)
params <- fitto_norm$par
params # Log scale abundance, ROI and mortality


# -------------------
# Pairwise difference
# -------------------

dim(prcs$samples.2)[1]
sum(prcs$comp.mat$HSP)

data <- list(deltaY = as.array(cmat$samp_years_2 - cmat$samp_years_1),
             deltaG = as.array(cmat$epi_ages_2 - cmat$epi_ages_1),
             hspsYN = as.array(cmat$HSP))

epi_nllk1 <- MakeADFun(data,
                       parameters = list(theta = params[1] - 0.5,
                                         psi   = params[3],
                                         beta  = 0.1,
                                         logsi = log(1.5)),
                       DLL = "hsp_auto_tmb")

fit1 = nlminb(epi_nllk1$par,
              epi_nllk1$fn,
              epi_nllk1$gr,
              control = list(trace = 1, xf.tol = 1e-14, 
                             rel.tol = 2e-14, sing.tol = 2e-14))
fit1
epi_nllk1$gr(fit1$par)
rep <- sdreport(epi_nllk1)
summary(rep) 

# 1/theta is adult abundance
# beta is the rate of change for noisy age | true age - should be close to 1
# exp(psi) is the mortality estimate 
# From above

                # Normal model  Pairwise difference
# Adult abund - 878.1778        983.12646716
# ROI         - -0.004          0 - not estimated assumed 0
# Mort        - 0.296986        0.2042723
# Beta        - -               # 0.78 


# ===============================
# Discrete deconvolution
# ===============================

all_ps_kin <- prcs$all.ps.kin.gy
all_ps_mt  <- prcs$all.ps.mt.gy

# Some more data prep

g <- prcs$samples.2$epi_ages_bin
y <- prcs$samples.2$samp_years
g_range <- prcs$mid.pts
grid_sz <- prcs$mid.pts[2] - prcs$mid.pts[1]
y_range <- seq(min(y), max(y))

samp_epi_age <- matrix(0,
                       nrow = length(g_range),
                       ncol = length(y_range))
rownames(samp_epi_age) <-  as.character(g_range)
colnames(samp_epi_age) <-  as.character(y_range)

tab <- table(g, y)
samp_epi_age[rownames(tab), colnames(tab)] <- tab

mx_pt_age <- max(prcs$samples.2$ages)
a_range   <- seq(0, mx_pt_age)
y0        <- min(y_range) - max(a_range) 
ylast     <- max(y_range) - min(a_range)
b_range   <- y0 %upto% ylast

# Get the non-zero indices

g_set <- cbind(g_range, seq(1, length(g_range)))
rownames(g_set) <- g_set[, 1]

y_set <- cbind(y_range, seq(1, length(y_range)))
rownames(y_set) <- y_set[, 1]

nz_indcs <- which(all_ps_mt != 0, arr.ind = T)
nz_indcs_2 <- data.frame(nz_indcs,
                         g1_y1 = rownames(all_ps_mt)[nz_indcs[, 1]],
                         g2_y2 = colnames(all_ps_mt)[nz_indcs[, 2]])

nz_indcs_2$Y1    <- gsub("_.*", "", nz_indcs_2$g1_y1)
nz_indcs_2$Y2    <- gsub("_.*", "", nz_indcs_2$g2_y2)
nz_indcs_2$G1    <- gsub(".*_", "", nz_indcs_2$g1_y1)
nz_indcs_2$G2    <- gsub(".*_", "", nz_indcs_2$g2_y2)
nz_indcs_2$G1Idx <- g_set[nz_indcs_2$G1, 2]
nz_indcs_2$G2Idx <- g_set[nz_indcs_2$G2, 2]
nz_indcs_2$Y1Idx <- y_set[nz_indcs_2$Y1, 2]
nz_indcs_2$Y2Idx <- y_set[nz_indcs_2$Y2, 2]
nz_indcs <- nz_indcs_2

# Set for cpp

nz_indcs <- as.matrix(nz_indcs[, c(1, 2, 9, 10, 11, 12)])
nz_indcs_cpp <- nz_indcs - 1

# -------------------
# Deconvolution setup
# -------------------

dc_prep <- deconvPrep(samp_epi_age, length(a_range), length(g_range), k = -1)


# -------------------
# Data prep for TMB
# -------------------

episd <- 0.5 # Known from simulation

# Some known age individuals
smps.lw.ags <- prcs$samples.2[prcs$samples.2$ages < 2, ]
kwn.a <- smps.lw.ags$ages[1:30]
kwn.g <- smps.lw.ags$epi_ages[1:30]

data <- list(kwn_a       = as.array(kwn.a), #as.array(seq(1, 30)), ##as.array(seq(1, 30)), #as.array(kwn.a),
             knw_g       = as.array(kwn.g), #as.array(seq(1, 30)), #as.array(knw.g), #as.array(seq(1, 30)),
             comb_yg     = as.array(all_ps_mt),
             kin_yg      = as.array(all_ps_kin),
             g_range     = g_range - grid_sz/2,
             a_range     = a_range,
             y_range     = y_range,
             b_range     = b_range,
             fng_rate    = 1 - 0,
             sig_prior   = 1,
             Mu_log_sig_ga = log(episd),
             SD_log_sig_ga = log(1.1),
             nz_inds     = as.array(nz_indcs_cpp),
             count       = as.array(samp_epi_age),
             Desmat_fix  = as.array(dc_prep$Desmat_fix),
             Desmat_rand = as.array(dc_prep$Desmat_rand),
             chol_S      = as.array(dc_prep$chol_S))

# Map what you like
#map  <- list(log_sig_ga = factor(NA))
#map  <- list(sig_ga = factor(NA))
map  <- list()
sig_ga_str <- episd

# Setup the TMB object
gamma_str <- matrix(0, nrow = 1, ncol = length(dc_prep$beta_start))
epi_autocal_nllk_lcky_l <- MakeADFun(data,
                                     parameters = list(N0_ad  = log(2000),
                                                       Z1     = log(0.1),
                                                       R1     = 0,
                                                       alpha  = 0,
                                                       Beta   = 1,
                                                       u      = dc_prep$u_start,
                                                       gammas = gamma_str,
                                                       logsdu = log(2),
                                                       log_sig_ga = log(sig_ga_str),
                                                       log_c  = log(1.1)),
                                     random = "u",
                                     map = map,
                                     DLL = "epi_auto_deconv_nllk_lcky_known_age_v2")


# -----------
# Numeric fit
# -----------

fitlcky_l <- nlminb(epi_autocal_nllk_lcky_l$par,
                   epi_autocal_nllk_lcky_l$fn,
                   epi_autocal_nllk_lcky_l$gr, 
                   lower = c(3,    -5,  -0.5, -10, -0.5,
                             gamma_str - 5, log(0.0001), log(0.1), log(0.4)),
                   control = list(trace = 1, rel.tol = 1e-14, sing.tol = 1e-14))

repo_lcky_l <- sdreport(epi_autocal_nllk_lcky_l, ignore.parm.uncertainty = F)

# Get out the parameters plus their transformations, random predictions, 
# predictions for abundance at each t (log(N_ad)) and the predicted noisy age 
# given true age 
repo_lcky_ls <- summary(repo_lcky_l)

              # Normal model    Pairwise difference   # Discrete deconvo
# Adult abund - 878.1778        983.12646716          807.8259 at lowest interval
# ROI         - -0.004          0                     -0.004
# Mort        - 0.296986        0.2042723             0.291031
# Beta        - -               # 0.78                1.027
# Alpha                                               -0.1135014 
# sig_g|a                                             0.5735679


# ---------
# Summaries 
# ---------

# The age distributions for each sampling year and the prediction g|y matrix
u_best           <- summary(repo_lcky_l, 'random')[,1]
dim(u_best)      <- dim(dc_prep$u_start)
gammas_best      <- summary(repo_lcky_l)[,1] %such.that% (names(.) == 'gammas')
dim(gammas_best) <- dim(gammas_best)

logtruhat <- dc_prep$Desmat_fix %*% gammas_best + dc_prep$Desmat_rand %*% u_best

truhat    <- exp(logtruhat)
pr_truhat <- truhat / rep( colSums( truhat), each = nrow( truhat))
round(pr_truhat, 3)

g_a <- epi_autocal_nllk_lcky_l$report()$g_a_beta_theta

samp_epi_age_hat <- g_a %*% truhat
round(samp_epi_age_hat) 


# Estimated and observed
# ----------------------

png(filename = paste0("results/test_summaries.png"), 
    bg = "white", 
    width =  14, height = 8, 
    units = 'in', res = 300)
par(mfrow = c(2, 3))

matplot( unclass(samp_epi_age_hat), 
         type='l', 
         main='Predicted distribution of noisy age in each sampling year', 
         xlab='Age', ylab='Freq')

matplot( unclass(samp_epi_age), 
         type='l', 
         main='Observed distribution of true age in each sampling year', 
         xlab='Age', ylab='Freq')

# Estimate a|y
matplot( unclass(pr_truhat), 
         type='l', 
         main='Estimated smooth distributions of sample age in each year', 
         xlab='Age', ylab='Freq')


# Abundance
# ---------

repo_lcky_lbp_s <- summary(repo_lcky_l)
abund    <- repo_lcky_ls[grep("log\\(N_ad\\)", rownames(repo_lcky_ls)), 1]
abund_se <- repo_lcky_ls[grep("log\\(N_ad\\)", rownames(repo_lcky_ls)), 2]
plot(b_range, exp(abund) , ylim = c(0,3000), typ = 'line', lwd = 2,
     ylab = "Adult abundance", xlab = 'Year')
lwr <- exp(abund - 1.96*abund_se)
upr <- exp(abund + 1.96*abund_se)
lines(b_range, lwr)
lines(b_range, upr)

pop_sz <- data.frame(data.frame(pop$census_postkill) %>%
                       group_by(year) %>%
                       dplyr::summarise(total = sum(female[age > 11]) + 
                                                sum(male[age > 11])))
lines(pop_sz$year, pop_sz$total)

abund_df <- data.frame(Abund = exp(abund), Upper = upr, Lower = lwr)


# Relationship age and vertebral age
# ----------------------------------

age_rel <- repo_lcky_ls[grep("ghat2", rownames(repo_lcky_ls)), 1]
age_rel_se <- repo_lcky_ls[grep("ghat2", rownames(repo_lcky_ls)), 2]

plot(seq(0,length(age_rel)-1), age_rel, 
     ylim = c(0, max(age_rel)), xlim = c(0, length(age_rel)),
     xlab = "True age", ylab = "Vertebral age", typ = 'line')
lwr <- age_rel - 1.96*age_rel_se
upr <- age_rel + 1.96*age_rel_se
lines(seq(0, length(age_rel)-1), lwr)
lines(seq(0, length(age_rel)-1), upr)
abline(a=0, b=1, col = 'red') # True line


# Simulate the prediction intervals
# ---------------------------------

tot_sim <- 10000
lgsdu_no <- grep("logsdu", names(repo_lcky_l$par.fixed))
#lgsdu.no <- c(14, 15)
mvn_draws <- rmvnorm(tot_sim, repo_lcky_l$par.fixed[-lgsdu_no], 
                     repo_lcky_l$cov.fixed[-lgsdu_no, -lgsdu_no])
mvn_draws_sub <- mvn_draws[, c("alpha", "Beta", "log_sig_ga")]
mvn_draws_sub[mvn_draws_sub[, 1] <= 0, 1] <- 0.00001
mx_age <- 30
res <- matrix(0, nrow = tot_sim, ncol = mx_age + 1)
for (a in seq(0, mx_age))
{
  print(a)
  for (n in seq(1, tot_sim))
  {
    mu  <- mvn_draws_sub[n, 1] + mvn_draws_sub[n, 2] * a
    sig <- exp(mvn_draws_sub[n, 3])
    res[n, a+1] <- rnorm(1, mu, sig)
  }
}
pt025 <- apply(res, 2, function(x) quantile(x, prob = 0.025))
pt975 <- apply(res, 2, function(x) quantile(x, prob = 0.975))
pt5   <- apply(res, 2, function(x) quantile(x, prob = 0.5))
lines(seq(0, mx_age), pt5,   col = 'red')
lines(seq(0, mx_age), pt025, col = 'red', lty = "dashed")
lines(seq(0, mx_age), pt975, col = 'red', lty = "dashed")
arel_df <- data.frame(Abund = age_rel, Upper = upr, Lower = lwr)
dev.off()
