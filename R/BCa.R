#' Calculate BCa (Bias-Corrected and Accelerated) Bootstrap Confidence Interval
#'
#' Computes the BCa confidence interval for a statistic using bootstrap and jackknife estimates.
#' The BCa method adjusts for both bias and skewness in the bootstrap distribution.
#'
#' @param obs Observed value of the statistic (scalar)
#' @param data Original dataset (data frame or matrix)
#' @param bootSample Vector of bootstrap replicates of the statistic (length = R)
#' @param bootjack Vector of jackknife replicates of the statistic (length = n)
#' @param level Confidence level (e.g., 0.95 for 95% confidence interval)
#'
#' @return A numeric vector of length 2 containing the lower and upper confidence limits.
#'
#' @details
#' The BCa method combines:
#'   - Bias correction (z0) based on the proportion of bootstrap replicates below the observed value
#'   - Acceleration (a) estimated from jackknife influence values
#'   - Normal quantiles adjusted for bias and acceleration
#'
#' @references
#' Efron, B. (1987). Better Bootstrap Confidence Intervals. *Journal of the American Statistical Association*, 82(397), 171-185.
#'
BCa <- function(obs, data, bootSample, bootjack, level){
R <- length(bootSample)
b <- qnorm((sum(bootSample > obs)+sum(bootSample==obs)/2)/R)
  
n <- nrow(data) 
n1 <- n-1 
obsn <- obs*n
pv <- i <- 0 
while(i < n){
  i = i+1 
pv[i] = obsn-n1*bootjack[i]
}
pv<-pv[!is.na(pv)]

je <- mean(pv)-pv
a <- sum(je^3)/(6*sum(je^2))^(3/2)

alpha <- (1-level)*2 
z <- qnorm(c(alpha/2,1-alpha/2)) # Std. norm. limits
p <- pnorm((z-b)/(1-a*(z-b))-b) # correct & convert to proportions

quantile(bootSample,p=p) # ABC percentile lims.      
}


