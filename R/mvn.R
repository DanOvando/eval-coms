#-------------------------------------------------------------
# Function to create multivariate-normal distribution for r-k
#-------------------------------------------------------------
mvn   <- function(n,mean.log.r,sd.log.r,mean.log.k,sd.log.k,cor.log.rk) {
  cov.log.rk <- cor.log.rk*sd.log.r*sd.log.k # covariance with empirical correlation and prior variances  covar.log.rk = matrix(NA, ncol=2,nrow=2)   # contract covariance matrix
  covar.log.rk      <- matrix(NA, ncol=2,nrow=2) # covariance matrix
  covar.log.rk[1,1] <- sd.log.r^2                # position [1,1] is variance of log.r
  covar.log.rk[2,2] <- sd.log.k^2               # position [2,2] is variance of log.k
  covar.log.rk[1,2] = covar.log.rk[2,1] = cov.log.rk     # positions [1,2] and [2,1] are correlations
  mu.log.rk  <- (c(mean.log.r,mean.log.k))      # vector of log.means
  mvn.log.rk <- rmvnorm(n,mean=mu.log.rk,sigma=covar.log.rk,method="svd")
  return(mvn.log.rk)
}
