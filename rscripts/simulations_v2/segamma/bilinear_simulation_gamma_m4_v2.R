# ==============================================================================
# Script to try and simulate the bilinear mode FM + lucky litter for 
# a more SHS like population
# Author: Luke Lloyd-Jones and Jack Easton
# Date started: 17/03/2025
# Date updated: 17/03/2025
# ==============================================================================

#library(CKMRpop)
#library(arules)
library(mvbutils)  
#library(offarray)  
#library(atease)    
#library(plyr)
library(dplyr)
#library(DNAage)
#library(deconvodisc)
library(mgcv)
library(TMB)
library(Rfast)

# Compilation
# -----------

source("code/rscripts/autocal/decovoPrep.R")

model1 <- 'seggamma'
out.rs <- out.rs <- paste0("autocal_simulation_", model1, "/")

system(paste0("ls ", out.rs, "*/processed* > ", out.rs, "/prcs_list.txt"))

args = commandArgs(trailingOnly = TRUE) 
args <- 15

lst <- read.table(paste0(out.rs, "/prcs_list.txt"))
file.i <- lst$V1[args[1]]

prcs <- get(load(file.i))
pop <- get(load(gsub("processed", "pop", file.i)))

# Make the results out file
ot.st <- gsub("/.*", "", file.i)
out.fl <- paste0(out.rs, ot.st)
dir.create(out.fl)

sub1 <- gsub(out.rs, "", file.i)
sub2 <- gsub("/processed_.*", "", sub1)
#sub2 <- gsub("/prcs_.*", "", sub1)
episd <- as.numeric(gsub(".*\\_", "", sub2))
outsub <- gsub("\\/", "", sub1)


# ----------------------
# Pre-processing
# ----------------------

plot(prcs$samples.2$ages, prcs$samples.2$epi_ages)
abline(a=0, b=1)

cmat        <- prcs$comp.mat
cmat$deltaA <- cmat$ages_2 - cmat$ages_1
cmat$deltaG <- cmat$epi_ages_2 - cmat$epi_ages_1
cmat$deltaY <- cmat$samp_years_2-cmat$samp_years_1
HSP_mat     <- as.matrix(cmat$HSP)
HSP_size    <- dim(HSP_mat)

sig_dA <- sd(cmat$deltaA)


# ===============================
# All guns blazing - lucky litter
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

smps.lw.ags <- prcs$samples.2[prcs$samples.2$ages < 2, ]
kwn.a <- smps.lw.ags$ages[1:30]
kwn.g <- smps.lw.ags$epi_ages[1:30]

data.gma <- list(kwn_a       = as.array(kwn.a),
                 knw_g       = as.array(kwn.g),
                 comb_yg     = as.array(all_ps_mt),
                 kin_yg      = as.array(all_ps_kin),
                 g_range     = g_range - grid_sz/2,
                 a_range     = a_range,
                 y_range     = y_range,
                 b_range     = b_range,
                 fng_rate    = 1 - 0.08,
                 phi_prior   = 1,
                 Mu_phi      = log(episd),
                 SD_phi      = log(1.05),
                 nz_inds     = as.array(nz_indcs_cpp),
                 count       = as.array(samp_epi_age),
                 Desmat_fix  = as.array(dc_prep$Desmat_fix),
                 Desmat_rand = as.array(dc_prep$Desmat_rand),
                 chol_S      = as.array(dc_prep$chol_S))

#map  <- list(log_sig_ga = factor(NA))
#map  <- list(sig_ga = factor(NA))
map  <- list(theta = factor(NA))
sig_ga_str <- episd

#dyn.unload(dynlib("code/cpp/epi_auto_deconv_nllk_lcky_known_age_gamma_bilinear"))
compile("code/cpp/epi_auto_deconv_nllk_lcky_known_age_gamma_bilinear.cpp")
dyn.load(dynlib("code/cpp/epi_auto_deconv_nllk_lcky_known_age_gamma_bilinear"))

gamma_str <- matrix(0, nrow = 1, ncol = length(dc_prep$beta_start))
epi_autocal_nllk_lcky_l <- MakeADFun(data.gma,
                                     parameters = list(N0_ad  = log(2000),
                                                       Z1     = log(0.1),
                                                       R1     = 0,
                                                       alpha  = 0.01,
                                                       beta   = 1.2,
                                                       bp     = 0.5,
                                                       theta  = 10,
                                                       u      = dc_prep$u_start,
                                                       gammas = gamma_str,
                                                       logsdu = log(0.01),
                                                       log_phi = log(sig_ga_str),
                                                       log_c  = log(1.1)),
                                     random = "u",
                                     map = map,
                                     DLL = "epi_auto_deconv_nllk_lcky_known_age_gamma_bilinear")


# -----------
# Numeric fit
# -----------

fitlcky = nlminb(epi_autocal_nllk_lcky_l$par,
                 epi_autocal_nllk_lcky_l$fn,
                 epi_autocal_nllk_lcky_l$gr, #N0_ad, R1,  alpha, beta, bp, theta,   gammas(8)
                 lower = c(3,    -5.3, -5,  0.0001,    0.0001, 0.001,  
                           gamma_str - 5,   log(0.0001), log(0.001), log(0.6)), # log(0.001),
                 upper = c(15,    3,    5,   0.5,     2,     2,  
                           gamma_str + 5,   log(0.4),     log(0.15), log(3)), # log(0.001),
                 control = list(trace = 1))

repo.lcky.l  <- sdreport(epi_autocal_nllk_lcky_l, ignore.parm.uncertainty = F)
repo.lcky.ls <- summary(repo.lcky.l)

save(repo.lcky.l,  file = paste0(out.fl, "/model_5_rep", i, ".Rdata"))


# ---------
# Summaries 
# ---------

u_best           <- summary(repo.lcky.l, 'random')[,1]
dim(u_best)      <- dim(dc_prep$u_start)
gammas_best      <- summary(repo.lcky.l)[,1] %such.that% (names(.) == 'gammas')
dim(gammas_best) <- dim(gammas_best)

logtruhat <- dc_prep$Desmat_fix %*% gammas_best + dc_prep$Desmat_rand %*% u_best

truhat    <- exp(logtruhat)
pr_truhat <- truhat / rep( colSums( truhat), each = nrow( truhat))
round(pr_truhat, 3)

g_a <- epi_autocal_nllk_lcky_l$report()$g_a_beta_theta

samp_epi_age_hat <- g_a %*% truhat
round(samp_epi_age_hat) 

par(mfrow = c(2, 2))
matplot( unclass(samp_epi_age_hat), 
         type='l', 
         main='Observed *adult* age from P.E.', 
         xlab='Age', ylab='Freq')

matplot( unclass(samp_epi_age), 
         type='l', 
         main='Observed *adult* age from P.E.', 
         xlab='Age', ylab='Freq')

matplot( unclass(pr_truhat), 
         type='l', 
         main='Observed *adult* age from P.E.', 
         xlab='Age', ylab='Freq')


# Abundance
# ---------
repo.lcky.lbp.s <- summary(repo.lcky.l)
abund    <- repo.lcky.ls[grep("log\\(N_ad\\)", rownames(repo.lcky.ls)), 1]
abund.se <- repo.lcky.ls[grep("log\\(N_ad\\)", rownames(repo.lcky.ls)), 2]
plot(b_range, exp(abund) , ylim = c(0,3000), typ = 'line', lwd = 2)
lwr <- exp(abund - 1.96*abund.se)
upr <- exp(abund + 1.96*abund.se)
lines(b_range, lwr)
lines(b_range, upr)
abline(h = 100)

pop.sz <- data.frame(data.frame(pop$census_postkill) %>%
                       group_by(year) %>%
                       dplyr::summarise(total = sum(female[age > 11]) + sum(male[age > 11])))
lines(pop.sz$year, pop.sz$total)

abund.df <- data.frame(Abund = exp(abund), Upper = upr, Lower = lwr)

save(abund.df, file = paste0(out.res, "abund_logdu", logsdu_str, ".Rdata"))


# Relationship age and vertebral age
# ----------------------------------

age.rel <- repo.lcky.ls[grep("ghat2", rownames(repo.lcky.ls)), 1]
age.rel.se <- repo.lcky.ls[grep("ghat2", rownames(repo.lcky.ls)), 2]

plot(seq(0,length(age.rel)-1), age.rel, 
     ylim = c(0, max(age.rel)), xlim = c(0, length(age.rel)),
     xlab = "True age", ylab = "Vertebral age", typ = 'line')
lwr <- age.rel - 1.96*age.rel.se
upr <- age.rel + 1.96*age.rel.se
lines(seq(0, length(age.rel)-1), lwr)
lines(seq(0, length(age.rel)-1), upr)
a <- seq(0, 30)
y <- sapply(a, function(x) segment2(x, cp = 10, beta1=1,beta2=0.36))
lines(a, y, col = 'blue')

# Simulate the prediction intervals
# ---------------------------------
tot.sim <- 10000
lgsdu.no <- grep("logsdu", names(repo.lcky.l$par.fixed))
#lgsdu.no <- c(14, 15)
mvn.draws <- rmvnorm(tot.sim, repo.lcky.l$par.fixed[-lgsdu.no], 
                     repo.lcky.l$cov.fixed[-lgsdu.no, -lgsdu.no])
mvn.draws.sub <- mvn.draws[, c("alpha", "beta", "bp", "log_phi")]
mvn.draws.sub[mvn.draws.sub[, 1] <= 0, 1] <- 0.00001
mx.age <- 30
res <- matrix(0, nrow = tot.sim, ncol = mx.age + 1)
for (a in seq(0, mx.age))
{
  print(a)
  for (n in seq(1, tot.sim))
  {
    mu <- segment2(a, cp = 10, alpha = mvn.draws.sub[n, 1], 
                   beta1 = mvn.draws.sub[n, 2], 
                   beta2 = mvn.draws.sub[n, 3])
    shape <- 1 / exp(mvn.draws.sub[n, 4])
    scale <- exp(mvn.draws.sub[n, 4]) * mu
    res[n, a+1] <- rgamma(1, shape = shape, scale = scale)
  }
}
pt025 <- apply(res, 2, function(x) quantile(x, prob = 0.025))
pt975 <- apply(res, 2, function(x) quantile(x, prob = 0.975))
pt5   <- apply(res, 2, function(x) quantile(x, prob = 0.5))
lines(seq(0, mx.age), pt5,   col = 'red')
lines(seq(0, mx.age), pt025, col = 'red', lty = "dashed")
lines(seq(0, mx.age), pt975, col = 'red', lty = "dashed")
arel.df <- data.frame(Abund = age.rel, Upper = upr, Lower = lwr)

