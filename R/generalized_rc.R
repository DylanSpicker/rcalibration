#' A Generalized Regression Calibration Estimator
#'
#' This function computes the imputed values based on a list of error-prone proxies of the true covariate.
#' 
#' @param W A list of length 'k' containing matrices of error-prone proxy measurements of the covariate. Matrices should all be n (observations) x p 
#'    (dimension of covariates).
#' @param Z A matrix containing all error-free covariates for use in estimation. Matrix should be n (observations) x q (dimension). 
#'    Use NULL if no such covariates exist. Defaults to NULL.
#' @param weights Either a string from {'numeric', 'optimal', 'equal'} (only required up to the point of unique identification) or a vector containing 'k' numbers, 
#'    summing to one, which serve as the convex combination of weights. Defaults to 'numeric'.
#' @param return_var A boolean represent whether the correction function and weights should be returned (TRUE) or only the imputed values. Defaults to FALSE.
#' @return Either a matrix of imputed values of size n x p (if return_var is FALSE), or a list which contains elements
#'    $X.hat (the aforementioned imputed matrix), $fitRC (a function which can be used to make the same correction), and
#'    $weights (the weights used in the correction).
#' @export
#' @examples
#' generalizedRC(W, weights="equal")

generalizedRC <- function(W, Z=NULL, weights="numeric", return_var=FALSE, ...) {
  #####################
  # Parameter Checking
  ######################
  # Run checks for 'W'
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

  ################ 
  # Compute Actual Parameters
  ###########################
  n <- nrow(W[[1]])
  p <- ncol(W[[1]])
  k <- length(W)

  # Run checks for 'Z'
  if (! is.null(Z)) {
    tryCatch({
      Z <- as.matrix(Z)
    }, warning = function(msg){
      warning(paste0("A warning was issued when trying to convert Z to be a matrix. (", msg, ")"))
    }, error = function(msg){
      stop(paste0("An error occurred while trying to convert Z to be a matrix. ERROR: ", msg))
    })

    if(nrow(Z) != n) stop(paste0("Non-conformable matrix sizes. W implies ", n, " observations, where Z implies ", nrow(Z), "."))

  }

  # Check on weights
  weight_techiques <- c('n', 'o','e')
  if (! (is.vector(weights) && length(weights) == k && is.numeric(weights) && (abs(sum(weights) - 1) <= 10e-16)) && 
      ! (inherits(weights, "character") && length(weights) == 1 && sum(startsWith(tolower(weights), weight_techiques)))) {
    stop(paste0("Weights must either (i) a vector of length k (", k, ") of numbers which sum to 1, or (ii) a string indication one of {[o]ptimal, [e]qual}."))
  }

  ################################
  # All Parameter Inputs are good!
  ###############################
  
  M_j <- NULL

  # Find weights
  if (length(weights) ==  1) {
    if (startsWith(weights, "n")){
      weights.obj <- solveWeights(W, ...)
      weights <- weights.obj$weights
      M_j <- weights.obj$M_j
    } else if (startsWith(weights, "o")) {
      weights.obj <- getOptimalWeights(W, ...)
      weights <- weights.obj$weights
      M_j <- weights.obj$M_j
    } else if (startsWith(weights, "e")) {
      weights <- rep(1/k, k)
    }
  } 

  if(is.null(M_j)) {
    # Compute the Mj estimator for pre-defined weights
    M_j <- getMj(W, ...)
  }

  # We have usable weights in 'weights', compute the estimator to use
  Xstar <- Reduce("+", lapply(1:k, function(ii){ as.matrix(weights[ii]*W[[ii]]) }))
  mu_x_hat <- colMeans(Xstar)
  SigmaXstar <- cov(Xstar)
  SigmaXX <- SigmaXstar - Reduce("+", lapply(1:k, function(ii){ (weights[ii]**2)*M_j[[ii]] }))

  # Set Z parameters to be null by default
  mu_z_hat <- NULL
  SigmaZ <- NULL
  SigmaXZ <- NULL

  correction <- function(Xstar, Z=NULL){
      if (! is.null(Z) ) warning("This model was fit without a 'Z', but you provided one for the correction. It is being ignored. ")
      return(t(mu_x_hat + SigmaXX %*% solve(SigmaXstar) %*% (t(Xstar) - mu_x_hat)))
  }

  if(! is.null(Z)) {
    mu_z_hat <- colMeans(Z)
    SigmaZ <- cov(Z)
    SigmaXZ <- cov(Xstar, Z)

    # Ensure that there is, in fact, variance present. 
    if (all(abs(SigmaXZ) <= 1e-6) || all(abs(SigmaZ <= 1e-6))) {
      warning("The provided Z has effectively no variance. The model is being fit without it.")
      Z <- NULL
    } else {
      correction <- function(Xstar, Z) {
        return(t(mu_x_hat + cbind(SigmaXX, SigmaXZ) %*% 
                            solve(rbind(cbind(SigmaXstar, SigmaXZ),
                                      cbind(t(SigmaXZ), SigmaZ))) %*% 
                            rbind(t(Xstar) - mu_x_hat, t(Z) - mu_z_hat)))
      }
    }

  }

  # Check if we need to construct list
  if (return_var) {
    return(list(X.hat = correction(Xstar, Z), 
                fitRC = correction,
                weights = weights))
  }

  correction(Xstar, Z)

}