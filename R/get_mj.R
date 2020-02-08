#' An Estimator for error-covariance matrix.
#'
#' The generalized estimator works on error models with covariance of the form cov(X) + Mj; this function estimates the residual Mj term.
#' 
#' @inheritParams generalizedRC
#' @return A list M_j, the error-covariance structure matrix.
#' @seealso [rcalibration::getOptimalWeights()] which calls this function to compute weights [rcalibration::generalizedRC()] which uses this if non-optimal weights are selected
#' @export
#' @examples
#' getMj(W)

getMj <- function(W) {
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

    Xstar.bar <- Reduce("+", W)/k

    M <- Reduce("+", lapply(1:k, function(idx) { 
        Xstar.cen <- W[[idx]] - Xstar.bar

        M_p <- matrix(rep(0, p*p), nrow=p, ncol=p)
        for (ii in 1:n) {
            M_p <- M_p + as.matrix(Xstar.cen[ii,])%*%t(as.matrix(Xstar.cen[ii,]))
        }

        1/(k*(k-1)) * M
    }))
    
    SigmaXX <- cov(Xstar.bar) - M

    lapply(W, function(w){ w - SigmaXX})

}