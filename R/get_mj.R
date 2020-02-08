#' An Estimator for error-covariance matrix.
#'
#' The generalized estimator works on error models with covariance of the form cov(X) + Mj; this function estimates the residual Mj term.
#' 
#' @inheritParams generalizedRC
#' @param enforce.psd Should the returned matrix be projected to be positive semi-definite. Can be a list (of the same length as W) or a boolean. 
#'  If TRUE, the corresponding Mj will be computed as the eigen-decomposition of Mj, with eigenvalues cast to 0. Defaults to FALSE.
#' @return A list M_j, the error-covariance structure matrix.
#' @seealso [rcalibration::getOptimalWeights()] which calls this function to compute weights [rcalibration::generalizedRC()] which uses this if non-optimal weights are selected
#' @export
#' @examples
#' getMj(W)

getMj <- function(W, enforce.psd=FALSE, ...) {
    ########
    # Parameter Checking
    #####################
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

    if( (! inherits(enforce.psd, "list") && ! is.logical(enforce.psd) ) && 
        ( inherits(enforce.psd, "list") && (length(enforce.psd) != length(W) || 
                                            ! all(unlist(lapply(enforce.psd, is.logical)))))) stop("enforce.psd must either be a list of TRUE/FALSE values, or TRUE/FALSE itself.")


    if ( ! inherits(enforce.psd, "list") ) enforce.psd <- lapply(1:length(W), function(x){ enforce.psd }) 

    #######
    # Compute Estimator
    ####################
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

        1/(n*k*(k-1)) * M_p
    }))

    SigmaXX <- cov(Xstar.bar) - M

    lapply(1:k, function(idx) {
        Mj <- cov(W[[idx]]) - SigmaXX
        warn <- FALSE
        if (enforce.psd[[idx]]) {
            eigen.d <- eigen(Mj)
            neg_ev <- which(eigen.d$values < 0)
            k_tr <- sum(eigen.d$values)

            if (length(neg_ev) == length(eigen.d$values)) {
                warn <- TRUE
            } else if(length(neg_ev) > 0) {
                eigen.d$values[neg_ev] <- 0
                eigen.d$values <- eigen.d$values*k_tr/sum(eigen.d$values)
                
                Mj <- eigen.d$vectors %*% diag(eigen.d$values) %*% solve(eigen.d$vectors)
            }
        }
        if (warn) warning(paste0("The computed matrix for proxy #", idx, " is negative definite. In this situaiton it is recommended to use fixed weights."))
        Mj
    })
    
}