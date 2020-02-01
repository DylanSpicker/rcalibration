#' The Regression Calibration Estimator
#'
#' This function implements the standard regression calibration estimator, based on replicated measurements.
#' 
#' @inherit generalizedRC params return
#' @seealso [rcalibration::generalizedRC()] which this function is just a wrapper for
#' @export
#' @examples
#' generalizedRC(W, weights="equal")
  
RC <- function(W, Z=NULL, return_var=FALSE) {
  # Proxy method for generalized RC w/ equal weights
  generalizedRC(W, Z, weights='e', return_var=return_var)
}