# ==============================================================================
# Script that takes a common large sampled data set and subsamples to smaller 
# sets of samples for each of the simulation runs. Script generates the 
# the replicates for analysis and does it for linear, factor, bi-linear and 
# bi-linear gamma.
# Author: Luke Lloyd-Jones
# Date started: 04/06/2025
# Date updated: 04/06/2025
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
library(mgcv)
library(Rfast)

# Function for the segmented regression mean
segment2 <- function(a, cp = 11,
                     beta1 = 1.1,
                     beta2 = 0.35,
                     alpha = 0){
  if (a <= cp) {
    trunc_vage2 <- alpha + beta1 * a   
  } else {
    trunc_vage2 <- alpha + beta1 * cp + beta2 * (a - cp)   
  }
  trunc_vage2
}

# Modified function to work with pre-sampled data
processSamplesLinearSegLinNormalGamma <- function(slurpie, beta1 = 1, beta2 = 0.36,
                                        alpha = 0, theta = 10,
                                        phi_episd = 0.5,
                                        grid_sz = 1,
                                        model = 'normal') {
  samples <- as.data.frame(slurpie$samples)
  samples$ages <- lapply(seq_along(samples$samp_years_list_post),
                         function(x) {
                           y <- samples$samp_years_list_post[[x]] - samples$born_year[x]
                           return(y)
                         })
  if (model == 'gamma')
  {
    samples$epi_ages <- lapply(samples$ages, function(a) {
                               mu <- segment2(a, beta1 = beta1, beta2 = beta2,
                                              alpha = alpha, cp = theta)
                               shape   = 1 / phi_episd;
                               scale   = phi_episd * mu;
                               y  <- rgamma(1, shape = shape, scale = scale)
                               return(y)
                               })
  } else if (model == 'normal') {
    
    samples$epi_ages <- lapply(samples$ages, function(a) {
      y = segment2(a, beta1 = beta1, beta2 = beta2, alpha = alpha, cp = theta) + rnorm(1, 0, phi_episd)
      return(y)
    })
    
  }
  
  rownames(samples) <- samples$ID      
  crel <- as.data.frame(compile_related_pairs(samples))
  
  # Now process    
  samples.2 <- samples[, c("ID", "sex", "samp_years_list", "ages")]
  samples.2$samp_years <- unlist(samples$samp_years_list)
  samples.2$epi_ages   <- unlist(samples$epi_ages)
  samples.2$ages       <- unlist(samples$ages)
  
  mn <- round_any(min(samples.2$epi_ages), grid_sz, f = floor)
  mx <- round_any(max(samples.2$epi_ages), grid_sz, f = ceiling)
  brks <- seq(mn, mx, grid_sz)
  epi_ages_bin <- cut(samples.2$epi_ages, breaks = brks, include.lowest = T)
  
  mid.pts <- brks + grid_sz/2
  mid.pts <- mid.pts[1:length(levels(epi_ages_bin))]
  levels(epi_ages_bin) <- mid.pts
  epi_ages_bin <- as.numeric(as.character(epi_ages_bin))
  samples.2$epi_ages_bin <- epi_ages_bin
  crel.hsps <- crel[crel$dom_relat == "Si" & crel$max_hit == 1, ]
  
  n <- dim(samples.2)[1]
  g <- expand.grid(row = 1:n, col = 1:n)
  m <- matrix(1, nrow = n, ncol = n)
  upr.inds <- g[upper.tri(m, diag = F), ]
  
  comp.mat <- cbind(samples.2[upr.inds$row, ], samples.2[upr.inds$col, ])
  colnames(comp.mat) <- c(paste0(colnames(samples.2), "_1"),
                          paste0(colnames(samples.2), "_2"))
  comp.mat$MATCH1 <- paste0(comp.mat$ID_1, "_", comp.mat$ID_2)
  crel.hsps$MATCH <- paste0(crel.hsps$id_1, "_", crel.hsps$id_2)
  comp.mat$MATCH2 <- paste0(comp.mat$ID_2, "_", comp.mat$ID_1)
  comp.mat$HSP    <- 0
  comp.mat$HSP[which(comp.mat$MATCH1 %in% crel.hsps$MATCH)] <- 1
  comp.mat$HSP[which(comp.mat$MATCH2 %in% crel.hsps$MATCH)] <- 1
  
  samples.2$Index <- seq(1, dim(samples.2)[1])
  index.mat <- data.frame(Index1 = samples.2[comp.mat$ID_1, "Index"],
                          Index2 = samples.2[comp.mat$ID_2, "Index"])
  
  Y         <- samples.2$samp_years
  Y_range   <- seq(min(Y), max(Y))
  mx.pt.age <- max(samples.2$ages)
  A_range   <- seq(0, mx.pt.age)
  
  y0 <- min(Y_range) - max(A_range)
  ylast <- max(Y_range) - min(A_range)
  B_range <- y0 %upto% ylast
  
  comp.mat$B1 <- comp.mat$samp_years_1 - comp.mat$ages_1
  comp.mat$B2 <- comp.mat$samp_years_2 - comp.mat$ages_2
  
  all.ps.comp <- matrix(0, nrow = length(B_range), ncol = length(B_range))
  rownames(all.ps.comp) <- as.character(B_range)
  colnames(all.ps.comp) <- as.character(B_range)
  cmp.tab <- table(comp.mat$B1, comp.mat$B2)
  all.ps.comp[rownames(cmp.tab), colnames(cmp.tab)] <- cmp.tab
  
  kin.tab <- table(comp.mat$B1[comp.mat$HSP == 1],
                   comp.mat$B2[comp.mat$HSP == 1])
  all.ps.kin <- all.ps.comp * 0
  all.ps.kin[rownames(kin.tab), colnames(kin.tab)] <- kin.tab
  
  # GY version
  G_range <- mid.pts
  comp.mat$Y1_G1 <- paste0(comp.mat$samp_years_1, "_", comp.mat$epi_ages_bin_1)
  comp.mat$Y2_G2 <- paste0(comp.mat$samp_years_2, "_", comp.mat$epi_ages_bin_2)
  all.poss <- paste0(rep(Y_range, times = length(G_range)),
                     "_", rep(G_range, each = length(Y_range)))
  all.ps.mt.gy <- matrix(0, nrow = length(all.poss), ncol = length(all.poss))
  colnames(all.ps.mt.gy) <- all.poss
  rownames(all.ps.mt.gy) <- all.poss
  cmp.tab <- table(comp.mat$Y1_G1, comp.mat$Y2_G2)
  all.ps.mt.gy[rownames(cmp.tab), colnames(cmp.tab)] <- cmp.tab
  
  kin.tab <- table(comp.mat$Y1_G1[comp.mat$HSP == 1],
                   comp.mat$Y2_G2[comp.mat$HSP == 1])
  all.ps.kin.gy <- all.ps.mt.gy * 0
  all.ps.kin.gy[rownames(kin.tab), colnames(kin.tab)] <- kin.tab
  
  # AY version
  comp.mat$Y1_A1 <- paste0(comp.mat$samp_years_1, "_", comp.mat$ages_1)
  comp.mat$Y2_A2 <- paste0(comp.mat$samp_years_2, "_", comp.mat$ages_2)
  all.poss <- paste0(rep(Y_range, times = length(A_range)),
                     "_", rep(A_range, each = length(Y_range)))
  all.ps.mt.ay <- matrix(0, nrow = length(all.poss), ncol = length(all.poss))
  colnames(all.ps.mt.ay) <- all.poss
  rownames(all.ps.mt.ay) <- all.poss
  cmp.tab <- table(comp.mat$Y1_A1, comp.mat$Y2_A2)
  all.ps.mt.ay[rownames(cmp.tab), colnames(cmp.tab)] <- cmp.tab
  
  kin.tab <- table(comp.mat$Y1_A1[comp.mat$HSP == 1],
                   comp.mat$Y2_A2[comp.mat$HSP == 1])
  all.ps.kin.ay <- all.ps.mt.ay * 0
  all.ps.kin.ay[rownames(kin.tab), colnames(kin.tab)] <- kin.tab
  
  return(list(samples.2 = samples.2, comp.mat = comp.mat, index.mat = index.mat,
              mid.pts = mid.pts, all.ps.comp = all.ps.comp, all.ps.kin = all.ps.kin,
              all.ps.mt.gy = all.ps.mt.gy, all.ps.kin.gy = all.ps.kin.gy,
              all.ps.mt.ay = all.ps.mt.ay, all.ps.kin.ay = all.ps.kin.ay))
}

# Simulation parameters
n_replicates <- 50
target_samp_fracs <- c(0.0025, 0.0035, 0.0045)
target_samp_fracs <- c(0.0045)

# Population parameters
samp_start_year <- 90
samp_stop_year <- 98
cohort_size <- 4000
mate_fidel <- 0

cat("Master population generated. Starting subsampling simulations...\n")

# Simulation parameters for processing

for (model1 in c('linear', "linear_fac", "segnormal", "seggamma"))
for (model1 in c("linear_fac", "segnormal", "seggamma"))
for (model1 in c("seggamma"))
{
#model1 <- 'linear'
if (model1 == 'linear')
{
  grid_sz <- 2
  beta1 <- 1
  beta2 <- 1
  alpha <- 0
  theta <- 10
  model <- 'normal'
  phi_values <- c(0.5, 1, 2)
  
} else if (model1 == "linear_fac")
{
  
  grid_sz <- 20
  beta1   <- 10
  beta2   <- 10
  alpha   <- 20
  theta   <- 10
  model   <- 'normal'
  phi_values <- c(2, 4, 8)
  
} else if (model1 == "segnormal")
{
  grid_sz <- 2
  beta1 <- 1
  beta2 <- 0.36
  alpha <- 0
  theta <- 10
  model <- 'normal'
  phi_values <- c(0.5, 1, 2)
  
} else if (model1 == "seggamma")
{
  grid_sz <- 2
  beta1 <- 1
  beta2 <- 0.36
  alpha <- 0
  theta <- 10
  model <- 'gamma'
  phi_values <- c(0.005, 0.01, 0.02)
}

out.rs <- paste0("autocal_simulation_", model1, "/")
dir.create(out.rs, showWarnings = FALSE)

# Progress tracking
total_sims <- length(target_samp_fracs) * length(phi_values) * n_replicates
sim_count <- 0

# Main simulation loop
for (samp_frac in target_samp_fracs) {
  for (phi in phi_values) {
    
    # Create output directory for this parameter combination
    out.fl <- paste0(out.rs, samp_frac, "_", phi)
    dir.create(out.fl, showWarnings = FALSE)
    
    cat(sprintf("Processing samp_frac = %.3f, phi = %.3f\n", samp_frac, phi))
  
    
    # Run replicates for this parameter combination
    for (i in 1:n_replicates) {
      sim_count <- sim_count + 1
      
      if (sim_count %% 10 == 0) {
        cat(sprintf("  Completed %d/%d simulations (%.1f%%)\n", 
                    sim_count, total_sims, 100 * sim_count / total_sims))
      }
      
      # Create subsample from master population
      pop <- DNAage::getCKMRpop(
        samp_frac = samp_frac,
        samp_start_year = samp_start_year,
        mat_yr = 11,
        mx_age = 25,
        cohort_size = cohort_size,
        samp_stop_year = samp_stop_year,
        mate_fidel = mate_fidel
      )
      
      # Calculate population size over time
      pop.sz <- data.frame(data.frame(pop$census_postkill) %>%
                             group_by(year) %>%
                             dplyr::summarise(total = sum(female[age > 11]) + sum(male[age > 11])))
      plot(pop.sz$year, pop.sz$total)
      
      pop.sz <- data.frame(pop$census_postkill %>%
                             group_by(year) %>%
                             dplyr::summarise(total = sum(female[age > 0]) + sum(male[age > 0])))
      plot(pop.sz$year, pop.sz$total)
      
      # Process the subsample
      prcs <- processSamplesLinearSegLinNormalGamma(
                 pop,
                 grid_sz = grid_sz,
                 phi_episd = phi,
                 beta1 = beta1,
                 beta2 = beta2,
                 alpha = alpha,
                 theta = theta,
                 model = model)
      
      # Save results
      save(pop, file = paste0(out.fl, "/pop_", i, ".Rdata"))
      save(prcs, file = paste0(out.fl, "/processed_", i, ".Rdata"))
      
      print(sum(prcs$all.ps.kin))
      print(sum(prcs$all.ps.kin.gy))
      plot(prcs$samples.2$ages, prcs$samples.2$epi_ages)
      abline(a=0, b=1)
    }
  }
}
}
