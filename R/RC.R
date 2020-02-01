#' The Regression Calibration Estimator
#'
#' This function implements the standard regression calibration estimator, based on replicated measurements.
#' 
#' @inheritsParams generalizedRC
#' @return Either a matrix of imputed values of size n x p (if return_var is FALSE), or a list which contains elements
#'    $X.hat (the aforementioned imputed matrix), $fitRC (a function which can be used to make the same correction), and
#'    $weights (the weights used in the correction).
#' @seealso [generalizedRC()] which this function is just a wrapper for
#' @export
#' @examples
#' generalizedRC(W, weights="equal")
  
RC <- function(W, Z=NULL, return_var=FALSE) {
  # Proxy method for generalized RC w/ equal weights
  generalizedRC(W, Z, weights='e', return_var=return_var)
}