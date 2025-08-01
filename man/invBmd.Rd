% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/invBmd.R
\name{invBmd}
\alias{invBmd}
\title{Inverse Regression for Benchmark Dose Calculation}
\usage{
invBmd(
  object,
  bmr,
  level = 0.9,
  slope,
  backgType = "modelBased",
  backg = NA,
  catLev = NA,
  extFactor = 10,
  def,
  useSD = useSD,
  sandwich.vcov = FALSE
)
}
\arguments{
\item{object}{A fitted dose-response model object containing the model function,
derivatives, coefficients, and data information}

\item{bmr}{numeric; Benchmark response level for which to calculate the BMD}

\item{level}{numeric; Confidence level for BMD confidence intervals (default: 0.9)}

\item{slope}{character; Direction of dose-response relationship, either
"increasing" or "decreasing"}

\item{backgType}{character; Type of background response calculation. Options include:
"modelBased", "hybridSD", "absolute"}

\item{backg}{numeric; Background response value (used with certain backgType options).
Default: NA}

\item{catLev}{numeric; Category level for multi-level responses. Default: NA}

\item{extFactor}{numeric; Extension factor for maximum dose search range (default: 10)}

\item{def}{character; BMD definition type. Options include:
\itemize{
\item "additional"/"added": Additional risk over background
\item "excess": Excess risk over background
\item "relative": Relative change from background
\item "extra": Extra risk (between background and maximum)
\item "hybridAdd": Hybrid additional risk using statistical criteria
\item "hybridExc": Hybrid excess risk using statistical criteria
\item "point": Point estimate (not supported for inverse regression)
}}

\item{useSD}{numeric; Standard deviation value used in hybrid definitions}

\item{sandwich.vcov}{logical; Whether to use sandwich variance-covariance matrix
for robust standard errors (default: FALSE)}
}
\value{
A matrix with 1 row and 3 columns:
\itemize{
\item \strong{BMD}: Benchmark dose estimate
\item \strong{BMDL}: Lower confidence limit for BMD
\item \strong{BMDU}: Upper confidence limit for BMD
}

If root finding fails for confidence limits, BMDL is set to 0 and BMDU to Inf.
}
\description{
Helper function for \code{bmd} that calculates Benchmark Dose (BMD) and confidence
intervals using inverse regression. This function handles multiple BMD definitions
and background response types by dynamically modifying the dose-response function
and its derivatives to solve for the dose corresponding to a specified benchmark response.
}
\details{
The function implements inverse regression by solving:
\deqn{f(BMD) = BMR}{f(BMD) = BMR}
where f is the transformed dose-response function based on the specified definition.

\strong{BMD Definitions:}

\strong{Additional Risk:} \deqn{f(x) = |f_0(0) - f_0(x)|}{f(x) = |f_0(0) - f_0(x)|}

\strong{Excess Risk:} \deqn{f(x) = \frac{|f_0(0) - f_0(x)|}{f_0(0)}}{f(x) = |f_0(0) - f_0(x)|/f_0(0)} (decreasing)
\deqn{f(x) = \frac{f_0(x) - f_0(0)}{1 - f_0(0)}}{f(x) = (f_0(x) - f_0(0))/(1 - f_0(0))} (increasing)

\strong{Relative Risk:} \deqn{f(x) = \frac{|f_0(0) - f_0(x)|}{f_0(0)}}{f(x) = |f_0(0) - f_0(x)|/f_0(0)}

\strong{Extra Risk:} \deqn{f(x) = \frac{|f_0(0) - f_0(x)|}{|f_0(0) - f_0(\infty)|}}{f(x) = |f_0(0) - f_0(x)|/|f_0(0) - f_0(infinity)|}

\strong{Hybrid Definitions:} Use normal distribution probabilities:
\deqn{f(x) = \Phi\left(\frac{f_0(0) - k \cdot SD - f_0(x)}{SD}\right) - (1-\Phi(k))}{f(x) = Phi((f_0(0) - k*SD - f_0(x))/SD) - (1-Phi(k))}
where k is the number of standard deviations (default: 2).
}
\note{
This is an internal helper function for BMD calculation. It assumes the input
object has the standard structure from dose-response model fitting, including
components like \code{fct}, \code{dataList}, \code{type}, etc.

\strong{Important considerations:}
\itemize{
\item The search range is limited to \verb{[0, extFactor × max(dose)]}
\item Hybrid definitions require careful specification of useSD and backg
\item Function modification uses string manipulation and may be sensitive to model structure
}
}
\section{Function Modification}{

The function dynamically modifies the model's derivative function to account
for different BMD definitions by:
\itemize{
\item Parsing the derivative function body
\item Adding appropriate transformation terms
\item Handling parameter indexing for fixed vs. estimated parameters
\item Supporting both 4-parameter and 5-parameter models
}
}

\section{Confidence Intervals}{

Confidence intervals are calculated using the delta method:
\deqn{CI = BMD \pm t_{\alpha/2} \sqrt{h(BMD)^T \Sigma h(BMD)}}{CI = BMD +/- t * sqrt(h(BMD)' * Sigma * h(BMD))}
where h(BMD) is the gradient vector and Σ is the variance-covariance matrix.

For continuous data, t-distribution quantiles are used; for other data types,
normal distribution quantiles are used.
}

\section{Error Handling}{

The function includes several error checks:
\itemize{
\item Stops if "point" definition is used (not supported)
\item Stops if hybrid definitions are used without appropriate backgType
\item Stops if no background value is provided when required
\item Stops if no solution exists for the specified BMR
}
}

\examples{
\dontrun{
# Typically called internally by bmd(), but can be used directly:
# Need to be updated
# Fit a dose-response model
model <- drm(response ~ dose, data = dose_data, fct = LL.4())

# Calculate BMD for 10\% additional risk
bmd_result <- invBmd(object = model, 
                     bmr = 0.1, 
                     level = 0.9,
                     slope = "increasing",
                     def = "additional",
                     useSD = sd(model$data$response))

# Extract results
bmd_estimate <- bmd_result[1, "BMD"]
bmd_lower <- bmd_result[1, "BMDL"] 
bmd_upper <- bmd_result[1, "BMDU"]

# Using hybrid definition with 1.5 standard deviations
bmd_hybrid <- invBmd(object = model,
                     bmr = 0.05,
                     def = "hybridAdd",
                     backgType = "hybridSD",
                     backg = 1.5,
                     useSD = sd(model$data$response))
}
}
\seealso{
\code{\link{bmd}} for the main BMD calculation function,
\code{\link{uniroot}} for the root finding algorithm,
\code{\link{vcov}} for variance-covariance matrix extraction,
\code{\link[sandwich]{sandwich}} for robust variance estimation
}
