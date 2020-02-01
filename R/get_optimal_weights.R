  #' An Estimator for "Optimal" Weight Selection
  #'
  #' This function returns the "optimal" weights for a set of 
  #' 
  #' @inheritsParams generalizedRC
  #' @return A list containing $weights, the optimally computed weights, as well as 
  #'    $M_j, the error-covariance structure matrix.
  #' @seealso [getMj()] which this function calls for the $Mj return, [generalizedRC()] which uses these weights
  #' @export
  #' @examples
  #' getOptimalWeights(W)

getOptimalWeights <- function(W) {
  #################
  # Run Error Checking for W
  #################
  if(! inherits(W, "list")) stop("W must be a list, each element of which is a matrix of error-prone covariates.")
  if(length(W) < 2) stop("W must have length >= 2.")
  
  # Convert each element to a matrix in W
  tryCatch({
    W <- lapply(W, as.matrix)
  }, warning = function(msg){
    warning(paste0("A warning was issued when trying to convert the elements of 'W' to be matrices. (", msg, ")"))
  }, error = function(msg){
    stop(paste0("An error occurred while trying to convert the elements of 'W' to be matrices. Please ensure that all elements are matrices. ERROR: ", msg))
  })
  if(length(unique(lapply(W, dim))) != 1) stop("Every element in W must have the same dimmensions.")

  n <- nrow(W[[1]])
  p <- ncol(W[[1]])
  k <- length(W)

  # Compute the weights
  M_j <- getMj(W)
  delta_j <- lapply(M_j, function(x){ return(sum(diag(x))) }) # Final Weights
  scale <- Reduce("+", delta_j)
  
  list(weights=unlist(lapply(delta_j, function(x){return(x/scale)})),
              M_j = M_j)
}