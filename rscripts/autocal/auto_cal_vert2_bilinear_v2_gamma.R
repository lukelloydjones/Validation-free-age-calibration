# ==============================================================================
# Try and get the HSP model going with the vertebral ages - version 2
# Author: Luke Lloyd-Jones
# Date started: 17/04/2024
# Date updated: 02/04/2025
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
library(ggwebthemes)
library(knitr)
library(cowplot)
library(Rfast)
source("code/rscripts/simulations_v2/helper_funcs.R")

#dyn.unload(dynlib("code/cpp/epi_auto_deconv_nllk_lcky_known_age_bilinear"))
#compile("code/cpp/epi_auto_deconv_nllk_lcky_known_age_tweedie_bilinear.cpp")
#dyn.load(dynlib("code/cpp/epi_auto_deconv_nllk_lcky_known_age_tweedie_bilinear"))


# Load the data and prepare
# -------------------------

kin <- get(load("data/og_data/kp2_for_Luke.rda"))
all <- get(load("data/og_data/MetaData_for_Luke.rda"))

all.sub <- all[, c("RegID", "Year_Collected", "Rings_Jan")]
all(kin$RegID1 %in% all.sub$RegID)
all(kin$RegID2 %in% all.sub$RegID)

kin.hsp <- kin[kin$isHSP, ]

# Create the pairwise matrix
all.sub <- all[, c("RegID.1", "Year_Collected", "Rings_Jan")]
p <- all.sub %>%
  ggplot( aes(x=Rings_Jan)) + xlab("Vertebral age") + ylab("Count") +
  geom_histogram(color = "#e9ecef", alpha = 0.6, position = 'identity') +
  theme_web_classic() +
  labs(fill = "")

png(filename = paste0("results/figures/dist_vert_age.png"), 
    bg = "white", 
    width =  8, height = 6, 
    units = 'in', res = 300)
print(p)
dev.off()

# Only use ages up to 11

only11 <- 0
if (only11)
{
  all.sub <- all.sub[all.sub$Rings_Jan <= 10, ]
}


# ---------
# Order them by guessing the ages
# ---------

all.sub$By <- all.sub$Year_Collected - all.sub$Rings_Jan


# ---------
# Known age
# ---------

shs.mth <- read.csv("results/reorganised_data/shs_meth_plus_factor.csv", 
                    header = T)
zro.age <- shs.mth[grep("SS", shs.mth$X), "Vertebrate_Age"]
kwn.a   <- rep(0, length(zro.age))
knw.g   <- zro.age
rownames(knw.g) <- NULL; colnames(knw.g) <- NULL


# -------------------------------
# OK let's try the proper version
# -------------------------------

n <- dim(all.sub)[1]

g <- expand.grid(row = 1:n, col = 1:n) # grid
m <- matrix(1, nrow = n, ncol = n)
upr.inds <- g[upper.tri(m, diag = F), ]
head(upr.inds)


# ----------------------------------------------------------
# Process the samples for HSP analysis
# ----------------------------------------------------------

grid_sz <- 2
epi_ages <- all.sub$Rings_Jan
mn <- round_any(min(epi_ages), grid_sz, f = floor)
mx <- round_any(max(epi_ages), grid_sz, f = ceiling)
brks <- seq(mn, mx, grid_sz)
epi_ages_bin <- cut(epi_ages, breaks = brks)
mid.pts      <- brks + grid_sz/2
levels(epi_ages_bin) <- mid.pts
epi_ages_bin <- as.numeric(as.character(epi_ages_bin))
epi_ages_bin[is.na(epi_ages_bin)] <- min(epi_ages_bin, na.rm = T)
all.sub$epi_ages_bin <- epi_ages_bin


n <- dim(all.sub)[1]
g <- expand.grid(row = 1:n, col = 1:n)
m <- matrix(1, nrow = n, ncol = n)
upr.inds <- g[upper.tri(m, diag = F), ]
head(upr.inds)
samples.2 <- all.sub

comp.mat <- cbind(all.sub[upr.inds$row, ], all.sub[upr.inds$col, ])
colnames(comp.mat) <- c(paste0(colnames(all.sub), "_1"), 
                        paste0(colnames(all.sub), "_2"))
head(comp.mat)
head(comp.mat)


# Attempt at ordering the birth data by guess at model
# ------------------------------

set1 <- comp.mat[, grep("_1", colnames(comp.mat))]
set2 <- comp.mat[, grep("_2", colnames(comp.mat))]

wch.ls <- comp.mat$By_2 < comp.mat$By_1
set1.tmp <- set1
set1[wch.ls, ] <- set2[wch.ls, ]
set2[wch.ls, ] <- set1.tmp[wch.ls, ]
comp.mat <- cbind(set1, set2)

comp.mat$MATCH1 <- paste0(comp.mat$RegID.1_1, "_",  comp.mat$RegID.1_2)
kin.hsp$MATCH   <- paste0(kin.hsp$RegID1, "_", kin.hsp$RegID2)
comp.mat$MATCH2 <- paste0(comp.mat$RegID.1_2, "_",  comp.mat$RegID.1_1)

comp.mat$HSP <- 0
comp.mat$HSP[which(comp.mat$MATCH1 %in% kin.hsp$MATCH)] <- 1
comp.mat$HSP[which(comp.mat$MATCH2 %in% kin.hsp$MATCH)] <- 1

all.sub$Index <- seq(1, dim(all.sub)[1])
index.mat <- data.frame(Index1 = all.sub[comp.mat$ID_1, "Index"], 
                        Index2 = all.sub[comp.mat$ID_2, "Index"])
Y <- all.sub$Year_Collected
Y_range   <- seq(min(Y), max(Y))

mx.pt.age <- max(all.sub$Rings_Jan) + 10
A_range <- seq(0, mx.pt.age)
y0      <- min(Y_range) - max(A_range)
ylast   <- max(Y_range) - min(A_range)
B_range <- y0 %upto% ylast


# -----------------------
# Back to our model setup
# -----------------------

G_range <- mid.pts
comp.mat$Y1_G1 <- paste0(comp.mat$Year_Collected_1, "_", comp.mat$epi_ages_bin_1)
comp.mat$Y2_G2 <- paste0(comp.mat$Year_Collected_2, "_", comp.mat$epi_ages_bin_2)
all.poss <- paste0(rep(Y_range, times = length(G_range)),  "_", 
                   rep(G_range,  each = length(Y_range)))

all.ps.mt.gy <- matrix(0, nrow = length(all.poss), ncol = length(all.poss))
colnames(all.ps.mt.gy) <- all.poss
rownames(all.ps.mt.gy) <- all.poss

cmp.tab <- table(comp.mat$Y1_G1, comp.mat$Y2_G2)
all.ps.mt.gy[rownames(cmp.tab), colnames(cmp.tab)] <- cmp.tab
kin.tab <- table(comp.mat$Y1_G1[comp.mat$HSP == 1], comp.mat$Y2_G2[comp.mat$HSP == 1])
all.ps.kin.gy <- all.ps.mt.gy * 0
all.ps.kin.gy[rownames(kin.tab), colnames(kin.tab)] <- kin.tab


# --------
# Optimise
# --------

all_ps_kin <- all.ps.kin.gy
all_ps_mt  <- all.ps.mt.gy

# Some more data prep

g <- all.sub$epi_ages_bin
y <- all.sub$Year_Collected
g_range <- mid.pts
y_range <- seq(min(y), max(y))

samp_epi_age <- matrix(0,
                       nrow = length(g_range),
                       ncol = length(y_range))
rownames(samp_epi_age) <-  as.character(g_range)
colnames(samp_epi_age) <-  as.character(y_range)

tab <- table(g, y)
samp_epi_age[rownames(tab), colnames(tab)] <- tab

mx_pt_age <- max(all.sub$Rings_Jan) + 10
a_range   <- seq(0, mx_pt_age)
y0        <- min(y_range) - max(a_range) # SHOULDN'T really be data-driven
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

knw.g <- knw.g + 0.0001
episd <- 0.05
sdsd  <- 1.2
# Prep the data

data.gma <- list(kwn_a       = as.array(kwn.a),
                 knw_g       = as.array(knw.g),
                 comb_yg     = as.array(all.ps.mt.gy),
                 kin_yg      = as.array(all.ps.kin.gy),
                 g_range     = g_range - grid_sz/2,
                 a_range     = a_range,
                 y_range     = y_range,
                 b_range     = b_range,
                 fng_rate    = 1 - 0.08,
                 phi_prior   = 1,
                 Mu_phi      = log(episd),
                 SD_phi      = log(sdsd),
                 nz_inds     = as.array(nz_indcs_cpp),
                 count       = as.array(samp_epi_age),
                 Desmat_fix  = as.array(dc_prep$Desmat_fix),
                 Desmat_rand = as.array(dc_prep$Desmat_rand),
                 chol_S      = as.array(dc_prep$chol_S))

map  <- list(R1 = factor(NA), theta = factor(NA))#, alpha = factor(NA))
#map  <- list(R1 = factor(NA), logsdu = factor(NA), theta = factor(NA))
#map  <- list(sig_ga = factor(NA), R1 = factor(NA), Z1 = factor(NA))
#map  <- list(sig_ga = factor(NA))

gamma_str <- matrix(0, nrow = 1, ncol = length(dc_prep$beta_start))
hist(exp(rnorm(10000, log(episd), log(sdsd))))

# ----------------------------
# If we look at bi-linear model
# ----------------------------

#dyn.unload(dynlib("code/cpp/epi_auto_deconv_nllk_lcky_known_age_gamma_bilinear"))
compile("code/cpp/epi_auto_deconv_nllk_lcky_known_age_gamma_bilinear.cpp")
dyn.load(dynlib("code/cpp/epi_auto_deconv_nllk_lcky_known_age_gamma_bilinear"))

out.res <- "results/autocal_vert_bl_gamma/"
dir.create(out.res)
logsdu_str <- 0.001
theta <- 11
epi_autocal_nllk_lcky_bl_gma <- MakeADFun(data.gma,
                                          parameters = list(N0_ad  = log(50000),
                                                            Z1     = log(0.2),
                                                            R1     = 0,
                                                            log_alpha  = log(0.001),
                                                            beta   = 1,
                                                            bp     = 0.5,
                                                            theta  = theta,
                                                            u      = dc_prep$u_start,
                                                            gammas = gamma_str,
                                                            logsdu  = log(logsdu_str),
                                                            log_phi = log(episd),
                                                            log_c   = log(1.1)),
                                          random = "u",
                                          map = map,
                                          DLL = "epi_auto_deconv_nllk_lcky_known_age_gamma_bilinear")

fitlcky = nlminb(epi_autocal_nllk_lcky_bl_gma$par,
                 epi_autocal_nllk_lcky_bl_gma$fn,
                 epi_autocal_nllk_lcky_bl_gma$gr, #N0_ad, R1,  alpha, beta, bp, theta,   gammas(8)
                 lower = c(6,    -5.3,  -10, 0.0001, 0.001,  
                           gamma_str - 5,   log(0.001), log(0.01), log(0.6)), 
                 control = list(trace = 1, eval.max = 500))

repo.lcky.lbp   <- sdreport(epi_autocal_nllk_lcky_bl_gma)
repo.lcky.lbp.s <- summary(repo.lcky.lbp)

save(repo.lcky.lbp, file = paste0(out.res, "repo_", theta, ".Rdata"))
save(epi_autocal_nllk_lcky_bl_gma, file = paste0(out.res, "adfun_", theta, ".Rdata"))
save(repo.lcky.lbp.s, file = paste0(out.res, "repo_s_", theta, ".Rdata"))

# ---------
# Summaries 
# ---------

u_best           <- summary( repo.lcky.lbp, 'random')[,1]
dim(u_best)      <- dim(dc_prep$u_start)
gammas_best      <- summary( repo.lcky.lbp)[,1] %such.that% (names(.) == 'gammas')
dim(gammas_best) <- dim( gammas_best)

logtruhat <- dc_prep$Desmat_fix %*% gammas_best + dc_prep$Desmat_rand %*% u_best

truhat    <- exp(logtruhat)
pr_truhat <- truhat / rep( colSums( truhat), each = nrow( truhat))
round(pr_truhat, 3)

#g_a <- prep_g_a(a_range, g_range, data$alpha, fitlcky$par["Beta"], sig_ga_str)
g_a <- epi_autocal_nllk_lcky_bl_gma$report()$g_a_beta_theta

samp_epi_age_hat <- g_a %*% truhat
round(samp_epi_age_hat) 

# Create a data frame from the predicted values
# Create data frame for predicted values (samp_epi_age_hat)
ages <- seq(0, 27, by = 2)  # Ages corresponding to rows (0, 2, 4, ..., 26)
years <- 2010:2017  # Column names
colnames(samp_epi_age_hat) <- years
pred_data <- data.frame(
  Age = rep(ages, times = length(years)),
  Value = as.vector(samp_epi_age_hat),
  Year = as.factor(rep(years, each = length(ages))),
  Type = "Predicted"
)

# Create data frame for true values (samp_epi_age)
colnames(samp_epi_age) <- years
true_data <- data.frame(
  Age = rep(ages, times = length(years)),
  Value = as.vector(samp_epi_age),
  Year = as.factor(rep(years, each = length(ages))),
  Type = "Observed"
)

# Combine both data frames
combined_data <- rbind(pred_data, true_data)

# Plot using ggplot2
pred.gs <- ggplot(combined_data, aes(x = Age, y = Value, color = Year, linetype = Type, group = interaction(Year, Type))) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(title = "",
       x = "Vertebral age (Years)",
       y = "Count") +
  theme_web_bw() +
  theme(legend.position = c(0.8, 0.5)) +
  scale_color_brewer(palette = "Set3") +  # Distinct colors for years
  scale_linetype_manual(values = c("Predicted" = "solid", 
                                   "Observed" = "dotted"))  # Solid for predicted, dotted for true

ages <- seq(0, 36, by = 1)  # Ages corresponding to rows (0, 2, 4, ..., 26)
years <- 2010:2017  # Column names
colnames(samp_epi_age_hat) <- years
pred_ages <- data.frame(
  Age = rep(ages, times = length(years)),
  Value = as.vector(pr_truhat),
  Year = as.factor(rep(years, each = length(ages))),
  Type = "Predicted"
)

priors <- ggplot(pred_ages, aes(x = Age, y = Value, 
                                color = Year)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(title = "",
       x = "Age (Years)",
       y = "Density") +
  theme_web_bw() +
  theme(legend.position = c(0.8, 0.5)) +
  scale_color_brewer(palette = "Set3") 

# Abundance
# ---------

abund    <- repo.lcky.lbp.s[grep("log\\(N_ad\\)", rownames(repo.lcky.lbp.s)), 1]
abund.se <- repo.lcky.lbp.s[grep("log\\(N_ad\\)", rownames(repo.lcky.lbp.s)), 2]
plot(b_range, exp(abund) /1000, ylim = c(0, 220), typ = 'line', lwd = 2)
lwr <- exp(abund - 1.96*abund.se)/1000
upr <- exp(abund + 1.96*abund.se)/1000
abund2 <- exp(abund)/1000
lines(b_range, lwr/1000)
lines(b_range, upr/1000)
abline(h = 100)

abund.df <- data.frame(Abund = exp(abund)/1000, Upper = upr, Lower = lwr)
abund.df$Brange <- b_range

abund.g <- ggplot(abund.df, aes(x = Brange, y = Abund, ymin = Lower, ymax = Upper)) +
  geom_line(linewidth = 1) +
  geom_ribbon(alpha = 0.3) + 
  ylim(c(0, 200)) +
  theme_web_bw() +
  labs(title = "",
       x = "Year",
       y = "Number of adults (1000s)") +
  theme(legend.position = "") +
  scale_color_brewer(palette = "Set3") 

save(abund.df, file = paste0(out.res, "abund.Rdata"))


# ----------------------------------
# Relationship age and vertebral age
# ----------------------------------

age.rel <- repo.lcky.lbp.s[grep("ghat2", rownames(repo.lcky.lbp.s)), 1]
age.rel.se <- repo.lcky.lbp.s[grep("ghat2", rownames(repo.lcky.lbp.s)), 2]

plot(seq(0,length(age.rel)-1), age.rel, 
     ylim = c(0, length(age.rel)), xlim = c(0, length(age.rel)),
     xlab = "True age", ylab = "Vertebral age", typ = 'line')
lwr <- age.rel - 1.96*age.rel.se
upr <- age.rel + 1.96*age.rel.se
lines(seq(0, length(age.rel)-1), lwr)
lines(seq(0, length(age.rel)-1), upr)
abline(a = 0, b = 1)

age.rel.df <- data.frame(Rel = age.rel, Upper = upr, Lower = lwr)
age.rel.df$Age <- seq(0,length(age.rel)-1)

# Simulate the prediction intervals
tot.sim <- 10000
lgsdu.no <- grep("logsdu", names(repo.lcky.lbp$par.fixed))
#lgsdu.no <- c(14, 15)
mvn.draws <- rmvnorm(tot.sim, 
                     repo.lcky.lbp$par.fixed[-lgsdu.no], 
                     repo.lcky.lbp$cov.fixed[-lgsdu.no, -lgsdu.no])
mvn.draws.sub <- mvn.draws[, c("log_alpha", "beta", "bp", "log_phi")]
mx.age <- 30
res <- matrix(0, nrow = tot.sim, ncol = mx.age + 1)
for (a in seq(0, mx.age))
{
  print(a)
  for (n in seq(1, tot.sim))
  {
    mu <- segment2(a, cp = theta, alpha = exp(mvn.draws.sub[n, 1]), 
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

age.rel.df$pt5 <- pt5
age.rel.df$pt025 <- pt025
age.rel.df$pt975 <- pt975

rel.g <- ggplot(age.rel.df, aes(x = Age, y = Rel)) +
  # Prediction interval ribbon (wider)
  geom_ribbon(aes(ymin = pt025, ymax = pt975), fill = "blue", alpha = 0.2) +
  # Confidence interval ribbon (narrower, layered on top)
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "blue", alpha = 0.5) +
  theme_web_bw() +
  # Predicted line
  geom_line(color = "black") +
  labs(x = "True age", y = "Vertebral age") 

save(age.rel.df, file = paste0(out.res, "agerel.Rdata"))


# ------------------------------
# Put them all on same plot
# ------------------------------

png(filename = paste0(out.res, "/autocal_relationship_", theta, ".png"), 
    bg = "white", 
    width =  14, height = 12, 
    units = 'in', res = 300)
print(plot_grid(pred.gs, priors, abund.g, rel.g, labels = "auto", nrow = 2))
dev.off()


# ------------------------------------------------------------
# Summary tables
# ------------------------------------------------------------

sig_ga_str <- 0.05
# Table the key parameters
# ------------------------

res.s <- get(load(paste0(out.res, "repo_s_11.Rdata")))
res.s.key <- data.frame(res.s[c("N0_ad", "Z1", "log_alpha", "beta", "bp", "log_phi", "log_c"), ])
res.s.key$Lower <- 0
res.s.key$Upper <- 0
res.s.key$Sigma <- sig_ga_str
res.s.key["N0_ad", "Upper"] <- exp(res.s.key["N0_ad", "Estimate"] + 2 * res.s.key["N0_ad", "Std..Error"])
res.s.key["N0_ad", "Lower"] <- exp(res.s.key["N0_ad", "Estimate"] - 2 * res.s.key["N0_ad", "Std..Error"])
res.s.key["N0_ad", "Estimate"] <- exp(res.s.key["N0_ad", "Estimate"])

res.s.key["Z1", "Upper"] <- exp(res.s.key["Z1", "Estimate"] + 2 * res.s.key["Z1", "Std..Error"])
res.s.key["Z1", "Lower"] <- exp(res.s.key["Z1", "Estimate"] - 2 * res.s.key["Z1", "Std..Error"])
res.s.key["Z1", "Estimate"] <- exp(res.s.key["Z1", "Estimate"])

res.s.key["alpha", "Upper"] <- exp((res.s.key["log_alpha", "Estimate"] + 2 * res.s.key["log_alpha", "Std..Error"]))
res.s.key["alpha", "Lower"] <- exp((res.s.key["log_alpha", "Estimate"] - 2 * res.s.key["log_alpha", "Std..Error"]))
res.s.key["alpha", "Estimate"] <- exp((res.s.key["log_alpha", "Estimate"]))

res.s.key["beta", "Upper"] <- (res.s.key["beta", "Estimate"] + 2 * res.s.key["beta", "Std..Error"])
res.s.key["beta", "Lower"] <- (res.s.key["beta", "Estimate"] - 2 * res.s.key["beta", "Std..Error"])
res.s.key["beta", "Estimate"] <- (res.s.key["beta", "Estimate"])

res.s.key["bp", "Upper"] <- (res.s.key["bp", "Estimate"] + 2 * res.s.key["bp", "Std..Error"])
res.s.key["bp", "Lower"] <- (res.s.key["bp", "Estimate"] - 2 * res.s.key["bp", "Std..Error"])
res.s.key["bp", "Estimate"] <- (res.s.key["bp", "Estimate"])

res.s.key["log_phi", "Upper"] <- exp(res.s.key["log_phi", "Estimate"] + 2 * res.s.key["log_phi", "Std..Error"])
res.s.key["log_phi", "Lower"] <- exp(res.s.key["log_phi", "Estimate"] - 2 * res.s.key["log_phi", "Std..Error"])
res.s.key["log_phi", "Estimate"] <- exp(res.s.key["log_phi", "Estimate"])

res.s.key["log_c", "Upper"] <- exp(res.s.key["log_c", "Estimate"] + 2 * res.s.key["log_c", "Std..Error"])
res.s.key["log_c", "Lower"] <- exp(res.s.key["log_c", "Estimate"] - 2 * res.s.key["log_c", "Std..Error"])
res.s.key["log_c", "Estimate"] <- exp(res.s.key["log_c", "Estimate"])

res.s.key.r <- round(res.s.key, 4)
kable(t(res.s.key.r), format = "latex")
