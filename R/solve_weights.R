#' A numeric estimator for Optimal Weight Selection
#'
#' This function iteratively solves the equations to estimate for the weights which minimize the MSE of the 
#'  linear approximation to the truth. That is, this function returns the estimated weights, when these weights
#'  are considered parameters in the BLUP equation.
#' 
#' @inheritParams generalizedRC
#' @param maxit The maximum number of iterations to run. Defaults to 500.
#' @param epsilon The tolerance at which estimates are deemed to have converged (when all weights change by less than epsilon). Defaults to 1e-10.
#' @return A list containing $weights, the optimally computed weights, as well as 
#'    $M_j, the error-covariance structure matrix.
#' @seealso [rcalibration::getMj()] which this function calls for the $Mj return, [rcalibration::generalizedRC()] which uses these weights
#' @export
#' @examples
#' solveWeights(W)

solveWeights <- function(W, maxit=500, epsilon=1e-10) {
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

    grabWeights <- function(cur_beta, M_j) {
        inverse_traces <- lapply(M_j, function(mtx){
            1/sum(diag(t(cur_beta)%*%cur_beta%*%mtx))
        })

        unlist(inverse_traces)/(Reduce("+", inverse_traces))
    }

    n <- nrow(W[[1]])
    p <- ncol(W[[1]])
    k <- length(W)

    # Compute the weights
    tryCatch({
        M_j <- getMj(W, enforce.psd=TRUE)
        cur_weights <- rep(1/k, k)
    }, warning=function(w){
        message("When solving for Mj, the following warning was received: '")
        message(w)
        message("As a result, equal weighting will be used.")
        return(list(weights=rep(1/k,k), M_j=M_j))
    })
    
    for(ii in 1:maxit) {
        cur_w <- Reduce("+", lapply(1:k, function(idx){ cur_weights[idx]*W[[idx]] }))
        Sigma_WW <- cov(cur_w)
        Sigma_XX <- Sigma_WW - Reduce("+", lapply(1:k, function(idx){ (cur_weights[idx]**2)*M_j[[idx]] }))
        cur_beta <- Sigma_XX%*%solve(Sigma_WW)

        new_weights <- grabWeights(cur_beta, M_j)

        if (max(new_weights - cur_weights) <= epsilon) {
            return(list(weights=new_weights, M_j=M_j))
        }

        cur_weights <- new_weights
    }

    warning("The process failed to converge. Using most recently computed weights.")
    list(weights=new_weights, M_j=M_j)
}
