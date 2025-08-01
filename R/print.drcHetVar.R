#' Print Method for drcHetVar Objects
#'
#' Prints a summary of a drcHetVar model object.
#'
#' @param x A drcHetVar model object
#' @param ... Additional arguments (not used)
#' @param digits Number of significant digits to use for printing values
#'
#' @return Invisibly returns the x object
#' @export
print.drcHetVar <- function(x, ..., digits = max(3, getOption("digits") - 3)){
  object <- x
  classList <- class(object)
  cat(paste("\n", "A 'drcHetVar' model.", "\n", sep = ""))
  cat("\nCall:\n", deparse(object$call), "\n\n", sep = "")
  
  # Coeffecients
  cat("Curve Coefficients:\n")
  print.default(format(object$curvePar, digits = digits), 
                print.gap = 2, quote = FALSE)
  cat("\n")
  
  cat("Variance Coefficients:\n")
  print.default(format(object$sigmaPar, digits = digits), 
                print.gap = 2, quote = FALSE)
  cat("\n")
  
  # END
  invisible(object)
}
