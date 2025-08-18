# Segmented regression function
# -----------------------------

segment2 <- function(a, cp = 11, 
                     beta1 = 1.1,
                     beta2 = 0.35,
                     alpha = 0)
{
  if (a <= cp)
  {
    trunc_vage2 <- alpha + beta1 * a 
  } else {
    trunc_vage2 <-  alpha + beta1 * cp  + beta2 * (a - cp) 
  }
  trunc_vage2
}


# Abundance trend predictions
# ---------------------------

getNPred <- function(prcs, params, fishSim = 0)
{
  Y <- prcs$samples.2$samp_years
  Y_range <- seq(min(Y), max(Y))
  mx.pt.age <- max(prcs$samples.2$ages)
  A_range <- seq(1, mx.pt.age)
  y0 <- min(Y_range) - max(A_range)
  if (fishSim == 1) {
    ylast <- max(Y_range)
  } else {
    ylast <- max(Y_range) - min(A_range)
  }
  B_range <- y0 %upto% ylast
  Nad_y0 <- exp(params[1])
  RoI <- params[2]
  Z <- exp(params[3])
  Nad_Y <- offarray(0, dimseq = list(Y = B_range))
  Nad_Y <- autoloop(Y = B_range, {
    Nad_y0 * exp(RoI * (Y - y0))
  })
  return(Nad_Y)
}

# Plot abundance?
# ---------------

plotAbund <- function(est, MDL = 0, trunAbun = tru.Ns.sum, 
                      fracs = c(0.003, 0.006, 0.01),
                      sds = c(0.5, 1, 1.5),
                      max.iter = 30, 
                      maxy = NULL)
{
  par(mfrow = c(3, 3))
  set <- seq(85, 96)
  if (MDL != 1 & MDL != 2 & MDL != 3 )
  {
    min.y <- min(est[, paste0("Y", set)])
    max.y <- quantile(c(as.matrix(est[, paste0("Y", set)])), probs = 0.95) + 
      0.1 * quantile(c(as.matrix(est[, paste0("Y", set)])), probs = 0.95)
  } else if (MDL == 1 | MDL == 2) {
    min.y <- min(est[, "X1.thetaT"])
    max.y <- max(est[, "X1.thetaT"]) +  max(est[, "X1.thetaT"])*0.1
  } else {
    min.y <- min(est[, "N"])
    max.y <- max(est[, "N"]) +  max(est[, "N"])*0.1
  }
  
  if (!is.null(maxy))
  {
    max.y <- maxy
  }
  
  for (samp_frac in fracs)
  {
    for (episd in sds)
    {
      # Make the results out file
      out.fl <- paste0(out.rs, samp_frac, "_", episd)
      k <- 0
      ests.sub <- est[est$SmpFrc == samp_frac & 
                        est$EpiSd == episd ,]
      max.iter <- max(ests.sub$Iter)
      
      for (i in seq(1, max.iter))
      {
        if (sum(est$SmpFrc == samp_frac & 
                est$EpiSd == episd & 
                est$Iter == i) > 0)
        {
          k  <- k + 1
          if (MDL != 1 & MDL != 2 & MDL != 3)
          {
            std <- est[est$SmpFrc == samp_frac & 
                         est$EpiSd == episd & 
                         est$Iter == i, paste0("Y", set)]
          } else if (MDL == 1 | MDL == 2) {
            std <- mean(est[est$SmpFrc == samp_frac & 
                              est$EpiSd == episd & 
                              est$Iter == i, "X1.thetaT"])
            std <- rep(std, length(set))
            
          } else if (MDL == 3) {
            std <- est[est$SmpFrc == samp_frac & 
                         est$EpiSd == episd & 
                         est$Iter == i, "N"]
            std <- rep(std, length(set))
          }
        } else {
          next
        }
        if (k == 1)
        {
          plot(set, std, type = 'l', 
               ylab = "Adult abundance",
               xlim = c(85, 96),
               ylim = c(min.y, max.y),
               xlab = "Year", 
               main = paste0("Scenario ", samp_frac, " - ", episd))
        } else {
          lines(set, std, type = 'l', 
                ylab = "Adult abundance",
                xlab = "Year")
          if (i == max.iter)
          {
            if (MDL != 1 & MDL != 2 & MDL != 3)
            {
              mns <- colMeans(est[est$SmpFrc == samp_frac & 
                                    est$EpiSd == episd, paste0("Y", set)])
              lwr <- sapply(est[est$SmpFrc == samp_frac & 
                                 est$EpiSd == episd, paste0("Y", set)],
                            function(x) quantile(x, prob = 0.025))
              upr <-  sapply(est[est$SmpFrc == samp_frac & 
                                   est$EpiSd == episd, paste0("Y", set)],
                             function(x) quantile(x, prob = 0.975))
            } else if (MDL == 1 | MDL == 2) {
              mns <- rep(mean( est[est$SmpFrc == samp_frac & 
                                     est$EpiSd == episd, "X1.thetaT"]), 
                         length(set))
              lwr <- rep(quantile( est[est$SmpFrc == samp_frac & 
                               est$EpiSd == episd, "X1.thetaT"], prob = 0.025), 
                               length(set))
              upr <- rep(quantile( est[est$SmpFrc == samp_frac & 
                                         est$EpiSd == episd, "X1.thetaT"], prob = 0.975), 
                         length(set))
            } else if (MDL == 3) {
              mns <- rep(mean( est[est$SmpFrc == samp_frac & 
                                     est$EpiSd == episd, "N"]), 
                         length(set))
            }
            lines(set, mns, type = 'l', 
                  col = "lightgreen", lwd = 3,
                  ylab = "Adult abundance",
                  xlab = "Year")
            lines(set, tru.Ns.sum[paste0("Y", set), 1], 
                  col = "red", lwd = 3)
            lines(set, lwr, 
                   col = "lightgreen", lwd = 2, lty = 2)
            lines(set, upr, 
                   col ="lightgreen", lwd = 2, lty = 2)
          }
        }
        
      }
    }
  }
}