deconvPrep <- function(count, ndomain, nrange, k = 9)
{
  
  ## convmat acts on a _range_ of allowed inputs, and maps them...
  ## ... to a (presumably bigger) _range_ of outputs
  
  DROP <- !length( dim( count))
  if( DROP){
    count <- matrix( count, length( count), 1L)
  } else {
    stopifnot( length( dim( count)) == 2)
  }
  stopifnot( nrow( count) == nrange)
  
  ndatasets <- ncol( count)
  
  faketru <- rep.int( 0L, ndomain)
  xtru <- seq_len( ndomain)
  
  flunge <- gam( faketru ~ s(xtru, k = k), 
                 family = poisson(link=log), 
                 fit = FALSE)
  
  #flunge <- gam( faketru ~ s(xtru), 
  #               family = poisson(link=log), 
  #               fit = FALSE)
  
  S    <- flunge$S[[1]] # assume only 1
  X    <- flunge$X
  off  <- flunge$off
  ifix <- 1 %upto% (off-1)
  
  nfix <- length( ifix) # better be 1!
  stopifnot( nfix==1)
  
  nrand <- ncol( X) - nfix
  
  #nullcols_S <- which( rowSums(abs(S)) == 0)
  nullcols_S <- which( rowSums(abs(S)) < 1e-10)
  S_nonnull  <- if ( length( nullcols_S)) S[ -nullcols_S, -nullcols_S] else S
  chol_Snn   <- chol( S_nonnull)
  chol_S     <- 0 * S
  chol_S[ 1:ncol( S_nonnull), 1:ncol( S_nonnull)] <- chol_Snn
  
  irand_reord <- off-1 + c( 1:nrand %except% nullcols_S, nullcols_S)
  
  # Don't bother with sparse matrices: small and dense is likely...
  # ... and 'Matrix' is a bit of a PITA
  
  Desmat_fix <- X[, ifix, drop=FALSE]
  Desmat_rand <- X[, irand_reord, drop=FALSE]
  
  
  # Start with intercept (presumably first col) only
  # Single dataset was:
  # beta_start <- c( log( sum( count)), rep( 0, nfix-1))
  # u_start <- rep( 0, nrand)
  
  beta_start <- matrix( 0, nfix, ndatasets)
  beta_start[1,] <- log( colSums( count))
  u_start <- matrix( 0, nrand, ndatasets)
  
  return(list(beta_start. = beta_start, 
              u_start     = u_start, 
              Desmat_fix  = Desmat_fix, 
              Desmat_rand = Desmat_rand, 
              chol_S      = chol_S))
  
}