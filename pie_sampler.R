##' Implementation of PIE algorithm
##'
##' PIE algorithm for combining MCMC samples of one-dimensional
##'     functionals from $k$ subset posterior distributions. Samples
##'     from every subset posterior distribution is obtained after
##'     modifying the likelihood using stochastic approximation. This
##'     is achieved by raising the likelihood to the power of $k$ before
##'     sampling. See Algorithm 1 in Li et al. (2017) for details.
##' @title PIE algorithm for combining MCMC samples of one-dimensional functionals
##' @param samplesList A list of length $K$ with numeric vectors of
##'     subset posterior samples as its elements. The length of the
##'     numeric vectors can differ.
##' @param meshsize A scalar specifiying grid size of the quantiles
##'     used to approximate the CDF.
##' @return Quantiles of the PIE posterior distribution.
##' @author Cheng Li (stalic@nus.edu.sg) and Sanvesh Srivastava (sanvesh-srivastava@uiowa.edu)
fitPie <- function (samplesList, meshsize = 1e-4) {
    qseq <- seq(0, 1, by = meshsize)
    qlist <- lapply(samplesList, function(x) quantile(x, prob = qseq))
    piePost <- rowMeans(do.call(cbind, qlist))
    piePost
}

##' Implementation of multivariate PIE algorithm
##'
##' Multivariate PIE algorithm for combining MCMC samples of
##'     multivariate parameters from $k$ subset posterior
##'     distributions. Samples from every subset posterior
##'     distribution is obtained after modifying the likelihood using
##'     stochastic approximation. This is achieved by raising the
##'     likelihood to the power of $k$ before sampling. See Section
##'     4.4 in Li et al. (2017) for details.
##' @title PIE algorithm for combining MCMC samples of multivariate parameters
##' @param samplesList A list of length $k$ with matrices of
##'     subset posterior samples as its elements. The MCMC samples are
##'     arranged along the rows and the parameter dimensions vary
##'     along the columns. The number of rows can differ but the
##'     number of columns must be the same.
##' @param idxList A list containing the indices of parameters for
##'     which the joint distributions are to be computed.
##' @param meshsize A scalar specifiying grid size of the quantiles
##'     used to approximate the CDF as in \ref{fitPie}.
##' @return List of samples from the PIE posterior distribution. The
##'     length is equal to the length of idxList.
##' @author Cheng Li (stalic@nus.edu.sg) and Sanvesh Srivastava (sanvesh-srivastava@uiowa.edu)
fitMultVarPie <- function (samplesList, idxList, meshsize) {
    nsub <- length(samplesList)

    pieJt <- list()
    for (dd in seq_along(idxList)) {
        idx <- idxList[[dd]]
        ## list of inv-covariances of subset-posteriors
        infEst <- lapply(samplesList, function (x) solve(cov(x[ , idx])))
        ## calculate Cov = (Cov1^{-1} + ... + CovK^{-1})^{-1}
        overallCov <- solve(Reduce("+", infEst))
        ## this will be used for scaling subset posterior samples
        overallCovSqrt <- chol(overallCov)
        overallCholInv <- solve(t(overallCovSqrt))
        ## centering
        centSamp <- list()
        mn <- colMeans(do.call(rbind, samplesList)[ , idx])
        centSamp <- lapply(samplesList, function(x) sweep(x[ , idx], 2, mn, "-"))
        ## scaling by overallCholInv
        orthSamp <- list()
        for (kk in 1:nsub) {
            orthSamp[[kk]] <- t(overallCholInv %*% t(centSamp[[kk]]))
        }
        ## apply PIE algorithm to the marginals
        tmp <- list()
        for (kk in 1:nsub) {
            tmp[[kk]] <- do.call(cbind,
                                 lapply(split(orthSamp[[kk]], col(orthSamp[[kk]])),
                                        function (x) {
                                            quantile(x, prob = seq(0.0, 1.0, length = 1 / meshsize))
                                        }))
        }
        orthPie <- matrix(NA, nrow = 1 / meshsize, ncol = length(idx))
        for (iii in seq_along(idx)) {
            orthPie[ , iii] <- rowMeans(do.call(cbind, lapply(tmp, function(x) x[ , iii])))
        }
        ## sample from the estimated marginals
        tmp <- list()
        for (iii in seq_along(idx)) {
            sidx <- sample(1:(1 / meshsize), size = 1 / meshsize, replace = TRUE)
            tmp[[iii]] <- orthPie[sidx, iii]
        }
        ## reverse the centering and scaling
        pieJt[[dd]] <- t(crossprod(overallCovSqrt, t(do.call(cbind, tmp))) + mn)
    }

    pieJt
}

##' The PIE algorithm
##'
##' This function implements most general version of the PIE algorithm
##' for combining posterior samples of possibly multivariate
##' parameters. The default strategy is to use the PIE
##' algorithm for calculating the Wasserstein barycenter of the
##' univariate marginals. If 'jointId' is provided as a list of
##' indices, then we use the hueristic algorithm in Section 4.4 of Li
##' et al. (2017) to compute the Wasserstein barycenter multivariate
##' distributions.
##' @title PIE algorithm for combining MCMC samples in univariate and
##' multivariate parameter samples
##' @param samplesList A list of length $K$ consisting of vectors or
##' matrices. If the list contains matrix, then the rows index the
##' MCMC samples and columns index the parameter dimensions. The
##' matrices and vectors can differ in the number of rows and length,
##' respectively, but the number of columns of every matrix must be
##' equal.
##' @param meshsize grid size of the quantiles used to approximate the CDF.
##' @param marginal return the barycenter of marginals for
##'     multivariate parameters.
##' @param jointId A list of indices. If specified, then 'pie' returns
##' the Wasserstein barycenter of the $K$ joint distributions
##' corresponding to the every element of the list.
##' @return An object of class 'pie', with kernel density estimates of
##' joints and marginals if specified.
##' @author  Cheng Li (stalic@nus.edu.sg) and Sanvesh Srivastava (sanvesh-srivastava@uiowa.edu)
pie <- function (samplesList, meshsize = 1e-4, marginal = TRUE, jointId = NULL, ...) {
    library(KernSmooth)

    if (class(samplesList) != "list") {
        stop("subset posterior samples must be a list")
    }

    if (class(samplesList[[1]]) == "matrix") {
        if (!all(sapply(samplesList, function(x) ncol(x)) == ncol(samplesList[[1]]))) {
            stop("subset posterior samples must have the same number of columns")
        }
        if (is.null(jointId)) {
            cat("Joint indices are not provided. The estimated posterior distribution is the Wasserstein barycenter of the marginals.\n")
        }
    }

    nsub <- length(samplesList)
    if (length(samplesList) == 1) {
        cat("No. of subsets = 1. Nothing to compute.")
        return(samplesList)
    }

    jdens <- NULL
    dens <- NULL
    startTime <- proc.time()
    if (class(samplesList[[1]]) == "numeric") {
        ndim <- 1
        samps <- samplesList
        piePost <- fitPie(samps, meshsize = meshsize)
        rr <- range(unlist(samps))
        bw <- dpik(piePost, range.x = rr)
        dens <- list(bkde(piePost, bandwidth = bw, range.x = rr))
    }

    if (marginal) { ## apply PIE algorithm to the marginals
        ndim <- ncol(samplesList[[1]])
        dens <- list()
        for (dd in 1:ndim) {
            samps <- lapply(samplesList, function(x) x[ , dd])
            piePost <- fitPie(samps, meshsize = meshsize)
            rr <- range(unlist(samps))
            bw <- dpik(piePost, range.x = rr, ...)
            dens[[dd]] <- bkde(piePost, bandwidth = bw, range.x = rr, ...)
        }
    }

    if (!is.null(jointId)) {
        pieSamp <- fitMultVarPie(samplesList, jointId, meshsize)
        jdens <- list()
        for (jj in seq_along(pieSamp)) {
            rng <- list()
            bw <- numeric()
            for (idd in seq_len(ncol(pieSamp[[jj]]))) {
                rng[[idd]] <- range(pieSamp[[jj]][ , idd])
                bw <- append(bw, dpik(pieSamp[[jj]][ , idd], range.x = rng[[idd]], ...))
            }
            jdens[[jj]] <- bkde2D(pieSamp[[jj]], bandwidth = bw, range.x = rng)
        }
    }
    endTime <- proc.time()

    res <- list(marginals = dens,
                joints = jdens,
                time = endTime[3] - startTime[3]
                )

    class(res) <- "pie"

    res
}

## generic plot function for the 'pie' class
plot.pie <- function(x, ...) {
    marg <- x$marginals
    jts <- x$joints

    ndim <- length(marg)
    par(ask = TRUE)
    if (!is.null(marg)) {
        for (dd in 1:ndim) {
            plot(x = marg[[dd]], ...)
        }
    }
    if (!is.null(jts)) {
        for (jj in seq_along(jts)) {
            contour(jts[[jj]]$x1, jts[[jj]]$x2, jts[[jj]]$fhat, ...)
        }
    }
    par(ask = FALSE)
}

## generic summary function for the 'pie' class
summary.pie <- function(x, ...) {
    library(matrixStats)
    marg <- x$marginals
    jts <- x$joints

    if (!is.null(marg)) {
        ll <- lapply(marg, function(x) quantile(x$y, probs = c(0.025, 0.5, 0.975)))
        qmat <- do.call(rbind, ll)
        rownames(qmat) <- paste0("dim", 1:2)
        cat("Marginal credible intervals: \n ")
        print(format(qmat, justify = "centre", digit = 2))
    }

    if (!is.null(jts)) {
        for (jj in seq_along(jts)) {
            gdim <- dim(jts[[jj]]$fhat)
            cat("Joint", jj, " is calculated on ", gdim[1], "by", gdim[2], "grid.\n")
        }
    }
}
