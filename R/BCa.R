#' Calculate BCa (Bias-Corrected and Accelerated) Bootstrap Confidence Interval
#'
#' @description
#' Computes the BCa confidence interval for a statistic using bootstrap and jackknife estimates.
#' The BCa method adjusts for both bias and skewness in the bootstrap distribution and generally
#' provides better coverage than standard bootstrap confidence intervals.
#'
#' @param obs numeric; The observed value of the statistic computed from the original sample
#' @param data data.frame or matrix; The original dataset from which bootstrap samples are drawn
#' @param bootSample numeric vector; Bootstrap replicates of the statistic. Must be of length R 
#'        where R is the number of bootstrap replicates
#' @param bootjack numeric vector; Jackknife replicates of the statistic. Must be of length n 
#'        where n is the number of observations in data
#' @param level numeric; Confidence level between 0 and 1 (e.g., 0.95 for 95% confidence interval)
#'
#' @return A numeric vector of length 2 containing:
#' \itemize{
#'   \item lower: The lower confidence limit
#'   \item upper: The upper confidence limit
#' }
#'
#' @details
#' The BCa method computes confidence intervals using three components:
#' 
#' 1. Bias correction (z0):
#' #$# z_0 = \Phi^{-1}(\frac{\#\{\theta^*(b) \leq \hat{\theta}\}}{B}) #$#
#' where B is the number of bootstrap replicates.
#'
#' 2. Acceleration (a):
#' #$# a = \frac{\sum_{i=1}^n (\bar{\theta}_{(\cdot)} - \theta_{(i)})^3}{6[\sum_{i=1}^n (\bar{\theta}_{(\cdot)} - \theta_{(i)})^2]^{3/2}} #$#
#' where θ_(i) are the jackknife values.
#'
#' 3. The confidence limits are then computed as:
#' #$# \alpha_{1,2} = \Phi(z_0 + \frac{z_0 \pm z_{\alpha}}{1 - a(z_0 \pm z_{\alpha})}) #$#
#'
#' @section Warning:
#' The function assumes that:
#' - The bootstrap samples are independent
#' - The statistic of interest is smooth (continuously differentiable)
#' - There are sufficient bootstrap replicates (typically >1000)
#' - The jackknife values are properly computed
#'
#' @examples
#' \dontrun{
#' # Example with mean as the statistic
#' data <- rnorm(100)
#' obs_mean <- mean(data)
#' 
#' # Generate bootstrap samples
#' R <- 1000
#' boot_means <- replicate(R, mean(sample(data, replace = TRUE)))
#' 
#' # Generate jackknife samples
#' jack_means <- sapply(1:length(data), 
#'                      function(i) mean(data[-i]))
#' 
#' # Calculate 95% CI
#' BCa(obs_mean, data, boot_means, jack_means, 0.95)
#' }
#'
#' @references
#' Efron, B. (1987). Better Bootstrap Confidence Intervals. 
#' *Journal of the American Statistical Association*, 82(397), 171-185.
#'
#' DiCiccio, T. J., & Efron, B. (1996). Bootstrap confidence intervals.
#' *Statistical Science*, 11(3), 189-212.
#'
#' @seealso 
#' \code{\link{quantile}} for the underlying quantile computation
#'
#' @export
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


