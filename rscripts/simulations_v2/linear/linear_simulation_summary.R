# ==============================================================================
# All 4 methods simulation - V2
#   - M1 - Assumes no error in the measured g-scores
#   - M2 - Has more honest SEs on Beta
#   - M3 - Should be complete in terms of errors
#   - M5 - FM + lucky litter
# Author: Luke Lloyd-Jones 
# Date started: 31/03/2025
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
library(ggplot2)
library(ggwebthemes)
library(knitr)
library(matrixStats)
source("code/rscripts/simulations_v2/helper_funcs.R")

out.rs <- "autocal_simulation_linear//"
dir.create(paste0(out.rs, "results"))

year.pss <- paste0("Y", seq(64, 97))
N.est.dummy <- data.frame(matrix(NA, ncol = length(year.pss)))
colnames(N.est.dummy) <- year.pss
M4N.est.dummy <- data.frame(matrix(NA, ncol = length(year.pss)))
colnames(M4N.est.dummy) <- year.pss

fracs <- c(0.0025, 0.0035, 0.0045)
sds   <- c(0.5, 1, 2)
reps  <- 50

for (samp_frac in fracs)
{
  for (episd in sds)
  {
    # Make the results out file
    out.fl <- paste0(out.rs, samp_frac, "_", episd)
    
    for (i in seq(1, reps))
    {
      print(paste0(samp_frac, "_", episd, "_",  i ))
      # Pop size
      pop <- get(load(file = paste0(out.fl, "/pop_", i, ".Rdata")))
      
      adults <- pop$census_postkill[pop$census_postkill$age > 11, ]
      pop.sz.all <- data.frame(pop$census_postkill %>% 
                             group_by(year) %>% 
                             dplyr::summarise(total = sum(female) + sum(male)))
      
      pop.sz <- data.frame(adults %>% 
                             group_by(year) %>% 
                             dplyr::summarise(total = sum(female) + sum(male)))
      rownames(pop.sz) <- pop.sz[, 1]
      
      N.tru <- data.frame(t(data.frame(pop.sz[seq(65, 96), 2])))
      names(N.tru) <- paste0('Y', seq(65, 96))
      N.tru$SmpFrc <- samp_frac
      N.tru$EpiSd  <- episd
      N.tru$Iter   <- i
      
      # Processed data
      prcs <- get(load(file = paste0(out.fl, "/processed_", i, ".Rdata")))
      #N.tru$Alpha <- coefficients(lm(prcs$samples.2$epi_ages ~ 
      #                                 prcs$samples.2$ages))[1]
      
      # Standard HSP model
      params <- get(load(file = paste0(out.fl, "/norm_par_", i, ".Rdata")))
      
      N.est <- data.frame(t(data.frame(getNPred(prcs, params))[, 2]))
      names(N.est) <- paste0("Y", data.frame(getNPred(prcs, params))[, 1])
      N.est.dummy[, names(N.est)]  <- N.est
      N.est <- N.est.dummy
      
      N.est$SmpFrc <- samp_frac
      N.est$EpiSd  <- episd
      N.est$Iter   <- i
      N.est$Va <- var(prcs$samples.2$ages)
      N.est$Va_VG_A <- var(prcs$samples.2$ages) / episd^2
      N.est$Mort <- exp(params[3])
      N.est$Alpha <- coefficients(lm(prcs$samples.2$epi_ages ~ 
                                       prcs$samples.2$ages))[1]
      
      print(var(prcs$samples.2$ages))
      
      # Auto-cal HSP model 1
      # --------------------
      
      repo    <- get(load(file = paste0(out.fl, "/model_1_rep", i, ".Rdata")))
      
      repo.s  <- summary(repo[[2]])
      model1.ests <- data.frame(t(repo.s[, 1]))
      model1.ests$SmpFrc <- samp_frac
      model1.ests$EpiSd  <- episd
      model1.ests$Iter   <- i
      model1.ests$Converg <- repo[[1]]$convergence
      
      model1.ses  <- data.frame(t(repo.s[, 2]))
      model1.ses$SmpFrc <- samp_frac
      model1.ses$EpiSd  <- episd
      model1.ses$Iter   <- i
      
      
      # Auto-cal HSP model 4
      # --------------------

        if (file.exists(paste0(out.rs, 
                               out.rs, 
                               "/model_5_rep_", samp_frac, "_", 
                               episd, "_", i, ".Rdata")))
        {
          repo4   <- get(load(file = paste0(out.rs, 
                                            out.rs, 
                                            "/model_5_rep_", samp_frac, "_", 
                                            episd, "_", i, ".Rdata")))
          repo4.s <- summary(repo4)
          
          model4.ests <- summary(repo4)[, 1]
          
          model4.Ns <- exp(model4.ests[grep("N_ad", names(model4.ests)) ])
          n.nms <- length(model4.Ns)
          names(model4.Ns) <- year.pss[(length(year.pss)-n.nms+1):length(year.pss)]
          M4N.est.dummy[, names(model4.Ns)]  <- model4.Ns
          model4.Ns <- M4N.est.dummy
          
          model4.Ns$SmpFrc <- samp_frac
          model4.Ns$EpiSd  <- episd
          model4.Ns$Iter   <- i
          model4.Ns$Va_VG_A <- var(prcs$samples.2$ages) / episd^2
          
          model4.Gs <- data.frame(repo4.s[grep("ghat2", rownames(repo4.s)), 1])
          model4.Gs$SmpFrc <- samp_frac
          model4.Gs$EpiSd  <- episd
          model4.Gs$Iter   <- i
          colnames(model4.Gs) <- c("Est", "SmpFrc", "EpiSd", "Iter")
          
          model4.ests.df <- data.frame(t(model4.ests[c("N0_ad", "Z1", "R1", "alpha", "Beta", "logsdu", "log_sig_ga", "log_c", "Z1_t")]))
          model4.ests.df$SmpFrc <- samp_frac
          model4.ests.df$EpiSd  <- episd
          model4.ests.df$Iter   <- i
          model4.ests.df$Converg <- 0

          
        } else {
          next
        }
        # model4.ses  <- data.frame(t(repo4.s[, 2]))
        # model4.ses <- model4.ses[, c("N0_ad", "Z1", "R1", "alpha", "Beta", "logsdu", "Z1_t")]
        # model4.ses$SmpFrc <- samp_frac
        # model4.ses$EpiSd  <- episd
        # model4.ses$Iter   <- i
        
        # Few extras
        
        N.est$NSamp <- dim(prcs$samples.2)[1]
        N.est$HSP   <- sum(prcs$all.ps.kin)
        
        # Storage
        
        if (samp_frac == fracs[1] & episd == sds[1] & i == 1)
        {
          N.tru.all <- N.tru
          N.est.all <- N.est
          
          model1.ests.all <- model1.ests
          model1.ses.all  <- model1.ses
          
          model4.ests.all <- model4.ests.df
          model4.Ns.all   <- model4.Ns
          model4.Gs.all   <- model4.Gs
          
        } else {
          
          N.tru.all <- rbind(N.tru.all, N.tru)
          N.est.all <- rbind(N.est.all, N.est[, colnames(N.est.all)])
          
          model1.ests.all <- rbind(model1.ests.all, model1.ests)
          model1.ses.all  <- rbind(model1.ses.all, model1.ses)
          
          model4.ests.all <- rbind(model4.ests.all, model4.ests.df)
          model4.Ns.all   <- rbind(model4.Ns.all, model4.Ns[, colnames(model4.Ns.all)])
          model4.Gs.all   <- rbind(model4.Gs.all, model4.Gs)
          
        }
    }
  }
}


# Save
# ----

save(N.tru.all,       file = paste0(out.rs, 'results/N.tru.all.Rda'))
save(N.est.all,       file = paste0(out.rs, 'results/N.est.all.Rda'))
save(model1.ests.all, file = paste0(out.rs, 'results/model1.ests.all.Rda'))
save(model4.Ns.all,   file = paste0(out.rs, 'results/model4.Ns.all.Rda'))
save(model4.ests.all, file = paste0(out.rs, 'results/model4.ests.all.Rda'))
save(model4.ests.all, file = paste0(out.rs, 'results/model4.ests.all.Rda'))
save(model4.Gs.all, file = paste0(out.rs, 'results/model4.Gs.all.Rda'))

get(load(paste0(out.rs, 'results/N.tru.all.Rda')))
get(load(paste0(out.rs, 'results/N.est.all.Rda')))
get(load(paste0(out.rs, 'results/model1.ests.all.Rda')))
get(load(paste0(out.rs, 'results/model4.Ns.all.Rda')))
get(load(paste0(out.rs, 'results/model4.ests.all.Rda')))

N.est.all %>% group_by(SmpFrc, EpiSd) %>% 
              summarise(Mn1 = mean(Mort), Sd = sd(Mort))

N.est.all %>% group_by(SmpFrc, EpiSd) %>% summarise(Mn1 = mean(NSamp), Mn2 = mean(HSP))
N.df <- data.frame(N.est.all %>% group_by(SmpFrc, EpiSd) %>% summarise(Mn1 = mean(Y90), 
                                                    Mn2 = mean(Y91), 
                                                    Mn3 = mean(Y92), 
                                                    Mn4 = mean(Y93), 
                                                    Mn5 = mean(Y94), 
                                                    Mn6 = mean(Y95), 
                                                    Mn7 = mean(Y96)))

mdl1.df <- data.frame(model1.ests.all %>% 
                        group_by(SmpFrc, EpiSd) %>% 
                        summarise(Mn1 = mean(exp(psi)), Sd = sd(exp(psi))))
mean(model1.ests.all$X1.thetaT)
mean(exp(model1.ests.all$psi))

N.m4.df <- data.frame(model4.Ns.all %>% group_by(SmpFrc, EpiSd) %>% summarise(Mn1 = mean(Y90), 
                                                                       Mn2 = mean(Y91), 
                                                                       Mn3 = mean(Y92), 
                                                                       Mn4 = mean(Y93), 
                                                                       Mn5 = mean(Y94), 
                                                                       Mn6 = mean(Y95), 
                                                                       Mn7 = mean(Y96)))
mean(rowMeans(N.m4.df[, -c(1, 2)]))

N.hspmod.df <- data.frame(N.est.all %>% group_by(SmpFrc, EpiSd) %>% summarise(Mn1 = mean(Y90), 
                                                                              Mn2 = mean(Y91), 
                                                                              Mn3 = mean(Y92), 
                                                                              Mn4 = mean(Y93), 
                                                                              Mn5 = mean(Y94), 
                                                                              Mn6 = mean(Y95), 
                                                                              Mn7 = mean(Y96)))
mean(rowMeans(N.hspmod.df[, -c(1, 2)]))
model4.ests.all <- data.frame(model4.ests.all)
mdl4.df <- data.frame(model4.ests.all %>% group_by(SmpFrc, EpiSd) %>% summarise(Mn1 = mean(Z1_t),
                                                                                Sd = sd(Z1_t)))

# -----------------------------------------------
# Try and plot the trend estimates for each model
# -----------------------------------------------

bs.out <- paste0(out.rs, "results")

# True
mat.N    <- as.matrix(N.tru.all[, grep("Y", colnames(N.tru.all))])
tru.N.mn <- colMeans(mat.N)
tru.N.qs <- colQuantiles(mat.N, 
                         probs = c(0.025, 0.975))
tru.Ns.sum <- cbind(tru.N.mn, tru.N.qs)
df_avg <- data.frame(N.tru.all %>%
                       group_by(SmpFrc, EpiSd) %>%
                       dplyr::summarise(across(.cols = starts_with("Y"), 
                                        .fns = mean, 
                                        .names = "{.col}",
                                        na.rm = TRUE)) %>%
                       ungroup())

out.name <- paste0(bs.out, "/true_abund")
png(filename = paste0(out.name, ".png"), 
    bg = "white", width =  14, height = 14, units = 'in', res = 300)
plotAbund(N.tru.all, sds = sds, fracs = fracs)
dev.off()


# Normal HSP model
# ----------------
df_avg <- data.frame(N.est.all %>%
                       group_by(SmpFrc, EpiSd,) %>%
                       summarise(across(.cols = starts_with("Y"), 
                                        .fns = mean, 
                                        .names = "{.col}",
                                        na.rm = TRUE)) %>%
                       ungroup())

out.name <- paste0(bs.out, "/norm_hsp_abund")
png(filename = paste0(out.name, ".png"), 
    bg = "white", width =  14, height = 14, units = 'in', res = 300)
plotAbund(N.est.all, sds = sds, fracs = fracs, maxy = 2300)
dev.off()

# Model 1

out.name <- paste0(bs.out, "/model1_abund")
png(filename = paste0(out.name, ".png"), 
    bg = "white", width =  14, height = 14, units = 'in', res = 300)
plotAbund(model1.ests.all, MDL = 1, sds = sds, fracs = fracs, maxy = 2300)
dev.off()

# Model 4

df_avg.m4 <- data.frame(model4.Ns.all %>%
                          group_by(SmpFrc, EpiSd, Iter) %>%
                          summarise(across(.cols = starts_with("Y"), 
                                           .fns = mean, 
                                           .names = "{.col}",
                                           na.rm = TRUE)) %>%
                          ungroup())

out.name <- paste0(bs.out, "/model4_abund")
png(filename = paste0(out.name, ".png"), 
    bg = "white", width =  14, height = 14, units = 'in', res = 300)
plotAbund(model4.Ns.all, sds = sds, fracs = fracs, maxy = 2300)
dev.off()


# -----------------------------------------------
# What everyone wants is the tick rate
# -----------------------------------------------

Avgs <- data.frame(N.est.all %>% group_by(SmpFrc, EpiSd) %>% dplyr::summarise(MeanVarRat = mean(Va_VG_A),
                                                                       MeanNsamp = mean(NSamp),
                                                                       MeanHSP   = mean(HSP)))
scn.avgs <- c("0.0025, 0.5, Nsamp = 300, NHSP = 70, V(A):V(G|A) ~ 70:1",
              "0.0025, 1,   Nsamp = 300, NHSP = 70, V(A):V(G|A) ~ 17:1",
              "0.0025, 2,   Nsamp = 300, NHSP = 70, V(A):V(G|A) ~ 4:1",
              "0.0035, 0.5, Nsamp = 420, NHSP = 140, V(A):V(G|A) ~ 70:1",
              "0.0035, 1,   Nsamp = 420, NHSP = 140, V(A):V(G|A) ~ 17:1",
              "0.0035, 2,   Nsamp = 420, NHSP = 140, V(A):V(G|A) ~ 4:1",
              "0.0045, 0.5, Nsamp = 530, NHSP = 220, V(A):V(G|A) ~ 70:1",
              "0.0045, 1,   Nsamp = 530, NHSP = 220, V(A):V(G|A) ~ 17:1",
              "0.0045, 2,   Nsamp = 530, NHSP = 220, V(A):V(G|A) ~ 4:4")

model1.ests.all.bta <- model1.ests.all[, c("beta", "SmpFrc", "EpiSd", "Iter")]
colnames(model1.ests.all.bta)[1] <- "Beta"
model1.ests.all.bta$Model   <- "Model 1"
#model1.ests.all.bta$VA_VG_A <- "Model 1"

model4.ests.all.bta <- model4.ests.all[, c("Beta", "SmpFrc", "EpiSd", "Iter")]
model4.ests.all.bta$Model <- "Model 4"

head(model4.ests.all.bta)
# Bind 

all.ests.bta <- rbind(model1.ests.all.bta,
                      model4.ests.all.bta)

all.ests.bta <- rbind(model1.ests.all.bta,
                      model4.ests.all.bta)
all.ests.bta$Scen <- paste0(all.ests.bta$SmpFrc, " - ", all.ests.bta$EpiSd)
all.ests.bta$Scen <- factor(all.ests.bta$Scen)
levels(all.ests.bta$Scen) <- scn.avgs


all.ests.bta$Model <- as.factor(all.ests.bta$Mode)
levels(all.ests.bta$Model) <- c("Pairwise difference", "Discrete deconvolution")
kable(all.ests.bta, format = 'latex')

all.ests.bta %>% group_by(SmpFrc, EpiSd, Model) %>% dplyr::summarise(mn = mean(Beta),
                                                              med = median(Beta),
                                                              sd = sd(Beta))


# Plotto
g <- ggplot(aes(y = Beta, x = Model), data = all.ests.bta) +
  geom_boxplot() +
  facet_wrap(~Scen, scales = "free_y") + 
  theme_web_bw() +
  ylab("Beta") + 
  xlab("") + ylim(0.5, 1.5) +
  theme(text  = element_text(size = 20, face = "bold"), 
        legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x  = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 15),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1),"cm")) 

out.name <- paste0(bs.out, "/beta_estimates.png")
png(filename = paste0(out.name, ".png"), 
    bg = "white", width =  18, height = 14, units = 'in', res = 300)
print(g)
dev.off()


# -----------------------------------------------
# Intercept
# -----------------------------------------------

model4.ests.all.bta <- model4.ests.all[, c("alpha", "SmpFrc", "EpiSd", "Iter")]
model4.ests.all.bta$Model <- "Model 4"

# Bind 

all.ests.bta <- rbind(model4.ests.all.bta)
all.ests.bta$Scen <- paste0(all.ests.bta$SmpFrc, " - ", all.ests.bta$EpiSd)
all.ests.bta$Scen <- factor(all.ests.bta$Scen)
levels(all.ests.bta$Scen) <- scn.avgs

all.ests.bta$Model <- as.factor(all.ests.bta$Mode)
levels(all.ests.bta$Model) <- c("Discrete deconvolution")

all.ests.bta %>% group_by(SmpFrc, EpiSd, Model) %>% dplyr::summarise(mn = mean(alpha),
                                                                     med = median(alpha),
                                                                     sd = sd(alpha))

# Plotto
g <- ggplot(aes(y = alpha, x = Model), data = all.ests.bta) +
  geom_boxplot() +
  facet_wrap(~Scen, scales = "free_y") + 
  theme_web_bw() +
  ylab("Alpha") + 
  xlab("") + ylim(-1, 1) +
  geom_hline(yintercept = 10, col = 'red') +
  theme(text  = element_text(size = 20, face = "bold"), 
        legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x  = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 15),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1),"cm")) 

out.name <- paste0(bs.out, "/alpha_estimates")
png(filename = paste0(out.name, ".png"), 
    bg = "white", width =  18, height = 14, units = 'in', res = 300)
print(g)
dev.off()


# -----------------------------------------------
# Variance
# -----------------------------------------------

model4.ests.all.bta <- model4.ests.all[, c("log_sig_ga", "SmpFrc", "EpiSd", "Iter")]
model4.ests.all.bta$sig_ga <- exp(model4.ests.all.bta$log_sig_ga)
model4.ests.all.bta$Model <- "Model 4"

# Bind 

all.ests.bta <- rbind(model4.ests.all.bta)
all.ests.bta$Scen <- paste0(all.ests.bta$SmpFrc, " - ", all.ests.bta$EpiSd)
all.ests.bta$Scen <- factor(all.ests.bta$Scen)
levels(all.ests.bta$Scen) <- scn.avgs

all.ests.bta$Model <- as.factor(all.ests.bta$Mode)
levels(all.ests.bta$Model) <- c("Discrete deconvolution")

# Plotto
# ------

g <- ggplot(aes(y = sig_ga, x = Model), data = all.ests.bta) +
  geom_boxplot() +
  geom_hline(data = all.ests.bta, aes(yintercept = EpiSd), col = "red") +
  facet_wrap(~Scen, scales = "free_y") + 
  theme_web_bw() +
  ylab("Sigma(G|A)") + 
  xlab("") + ylim(0, 2.5) +
  theme(text  = element_text(size = 20, face = "bold"), 
        legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x  = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 15),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1),"cm")) 

out.name <- paste0(bs.out, "/sig_g_a_estimates")
png(filename = paste0(out.name, ".png"), 
    bg = "white", width =  18, height = 14, units = 'in', res = 300)
print(g)
dev.off()



# -----------------------------------------------
# Plot predicted curves and the true curve
# -----------------------------------------------

# Get the observed data
# ---------------------

fracs <- c(0.0025, 0.0035, 0.0045)
sds   <- c(0.5, 1, 2)
reps  <- 50
model4.Gs.all$Age <- 0
for (samp_frac in fracs)
{
  for (episd in sds)
  {
    # Make the results out file
    out.fl <- paste0(out.rs, samp_frac, "_", episd)
    print(out.fl)
    for (i in seq(1, reps))
    {
      gis <- which(model4.Gs.all$SmpFrc == samp_frac & 
            model4.Gs.all$EpiSd ==   episd &
              model4.Gs.all$Iter ==   i)
      model4.Gs.all$Age[gis] <- seq(0, length(gis))
      # Processed data
      prcs <- get(load(file = paste0(out.fl, "/processed_", i, ".Rdata")))
      head(prcs$samples.2)
      age.vage.i <- data.frame(Ages = prcs$samples.2$ages,
                               Vages = prcs$samples.2$epi_ages)
      if (i == 1)
      {
        age.vage <- age.vage.i
      } else {
        age.vage <- rbind(age.vage, age.vage.i)
      }
    }
    
    # Make the percentiles at 0.025, 0.975
    a <- seq(0, max(age.vage$Ages))
    res.a <- sapply(a, function(a) {quantile(age.vage$Vages[age.vage$Ages==a], 
                                             probs = c(0.025, 0.975))})
    res.a.df <- data.frame(Age = a, Lwr = res.a[1, ], Upr = res.a[2, ], 
                           SampFrac = samp_frac,
                           Episd    = episd)
    res.a.df$Lwr[res.a.df$Age == 0] <- quantile(rnorm(1000, 0, res.a.df$Episd[1]), 0.025)
    res.a.df$Upr[res.a.df$Age == 0] <- quantile(rnorm(1000, 0, res.a.df$Episd[1]), 0.975)
    
    if (samp_frac == fracs[1] & episd == sds[1])
    {
      res.quant <- res.a.df
    } else {
      res.quant <- rbind(res.quant, res.a.df)
    }
  }
}

model4.Gs.all$Scen <- paste0(model4.Gs.all$SmpFrc, " - ", model4.Gs.all$EpiSd)
model4.Gs.all$Scen <- factor(model4.Gs.all$Scen)

res.quant$Scen <- paste0(res.quant$SampFrac, " - ", res.quant$Episd)

plot_data <- merge(model4.Gs.all, res.quant[, c("Age", "Lwr", "Upr", "Scen")], 
                   by = c("Age", "Scen"), all.x = TRUE)
#   geom_line(data = model4.Gs.all, aes(y = True, x = Age), col = 'black',
#linewidth = 5) +
#geom_line() +

#plot_data$Lwr[is.na(plot_data$Lwr)] <- 0
#plot_data$Upr[is.na(plot_data$Upr)] <- 0

g <- ggplot(aes(y = Est, x = Age, col = factor(Iter)), data = plot_data) +
  geom_line() +
  geom_ribbon(aes(ymin = Lwr, ymax = Upr), fill = "grey80", 
              alpha = 0.7, col = NA) +
  facet_wrap(~Scen, scales = "free_y") + 
  theme_web_bw() +
  ylab("Estimated age") + 
  xlab("True age") + ylim(-5, 25) + xlim(0, 20) +
  theme(text  = element_text(size = 20, face = "bold"), 
        legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x  = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 15),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1),"cm")) 


out.name <- paste0(bs.out, "/linear_curve_ests")
png(filename = paste0(out.name, ".png"), 
    bg = "white", width =  12, height = 12, units = 'in', res = 300)
  print(g)
dev.off()



# -----------------------------------------------
# More slick summary tables
# -----------------------------------------------

mod1 <- data.frame(model1.ests.all[model1.ests.all$beta > -1, ] %>% # 4 weird ones.
                     group_by(SmpFrc, EpiSd) %>% 
                     dplyr::summarise(mn_beta = mean(beta),
                                      sd_beta = sd(beta),
                                      mn_Z1  = mean(exp(psi)),
                                      sd_Z1  = sd(exp(psi)),
                                      mn_c  = mean(si),
                                      sd_c  = sd(si)))

kable(round(mod1, 3), format = 'latex')

mod4 <- data.frame(model4.ests.all %>% group_by(SmpFrc, EpiSd) %>% 
  dplyr::summarise(mn_alpha = mean(alpha),
                   sd_alpha = sd(alpha),
                   mn_beta  = mean(Beta),
                   sd_beta  = sd(Beta),
                   mn_Z1  = mean(Z1_t),
                   sd_Z1  = sd(Z1_t),
                   mn_R1  = mean(R1),
                   sd_R1  = sd(R1),
                   mn_sigmaG  = mean(exp(log_sig_ga)),
                   sd_sigmaG  = sd(exp(log_sig_ga)),
                   mn_psi = mean(exp(log_c)),
                   sd_psi  = sd(exp(log_c))))

kable(round(mod4, 3), format = 'latex')
