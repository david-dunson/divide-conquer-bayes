##' Data augmentation for variable selection using GDP prior in subsets
##'
##' Sample from the posterior distribution of regression coefficients
##' and noise variance in linear regression using GDP prior after
##' stochastic approximation. This modifies an existing data
##' augmentation algorithm by raising the likelihood to the power of
##' 'nrep', which is $k$, and fixing $eta = 1$. See Argmagan et al.
##' (2013) for details.
##' @title Stochastic approximation based Gibbs sampler for variable
##'     selection
##' @param yvec response.
##' @param xmat the design matrix.
##' @param nrep the exponent in stochastic approx.
##' @param ngrid the grid size for sampling the $\alpha$ parameter in
##'     the GDP prior.
##' @param niter no. of sampling iterations.
##' @param nburn no. of samples to discard as burnins.
##' @param nthin no. of samples to discard during thinning.
##' @return list of posterior samples for regression coefficients and
##'     noise variance.
##' @author Cheng Li (stalic@nus.edu.sg) and Sanvesh Srivastava (sanvesh-srivastava@uiowa.edu)
sampleSubGdp <- function (yvec, xmat, nrep, ngrid = 100, niter, nburn, nthin) {
    library(mgcv)
    library(MCMCpack)
    library(matrixStats)

    ndim <- ncol(xmat)
    nobs <- nrow(xmat)

    cts <- 0
    sigmaSamp <- rep(0.0, (niter - nburn) / nthin)
    alphaSamp <- rep(0.0, (niter - nburn) / nthin)
    etaSamp <- rep(0.0, (niter - nburn) / nthin)
    betaSamp <- matrix(0.0, (niter - nburn) / nthin, ndim)

    eta <- 1
    alpha <- 1
    betas <- runif(ndim, 1, 10)
    sigma2 <- 1
    lambdas <- runif(ndim, 1, 10)
    taus <- runif(ndim, 1, 10)
    agrid <- seq(0.01, 0.99, length = ngrid)
    logWtAlpha <- rep(0, length(agrid))
    nuErr <- 2; aErr <- 1; aA <- 10^2

    gramMat <- crossprod(xmat, xmat)

    startTime <- proc.time()
    for (its in 1:niter) {
        covBeta <- sigma2 * chol2inv(chol(nrep * gramMat + diag(1 / taus)))
        muBeta <- nrep * covBeta %*% crossprod(xmat, yvec) / sigma2
        betas <- as.numeric(muBeta + crossprod(chol(covBeta), rnorm(ndim)))

        ## update post. for hyper parameters in half t
        aErr <- rinvgamma(1, shape = 0.5 * (nuErr + 1), scale = nuErr / sigma2 + 1 / aA^2)
        ## update post. for error var.
        resids <- yvec - as.numeric(xmat %*% betas)
        shp <- 0.5 * (nobs * nrep + nuErr)
        scl <- nrep * sum(resids^2) / 2 + nuErr / aErr
        sigma2 <- rinvgamma(1, shape = shp, scale = scl)

        for (pp in 1:ndim) {
            lambdas[pp] <- rgamma(1, shape = alpha + 1, scale = eta + abs(betas[pp]) / sqrt(sigma2))
            taus[pp] <- 1 / rig(1, mean = abs(lambdas[pp] * sqrt(sigma2) / betas[pp]), scale = 1 / lambdas[pp]^2)
        }

        logWtAlpha <- ndim * log((1 - agrid) / agrid) - sum(log(1 + abs(betas) / (eta * sqrt(sigma2)))) / agrid
        probsAlpha <- exp(logWtAlpha - max(logWtAlpha))
        idx <- sample(seq_along(agrid), 1, prob = probsAlpha)
        alpha <- 1 / agrid[idx] - 1

        if (its %% 1000 == 0) cat("iteration: ", its, "\n")

        if (its > nburn & its %% nthin == 0) {
            cts <- cts + 1
            sigmaSamp[cts] <- sigma2
            betaSamp[cts, ] <- betas
            alphaSamp[cts] <- alpha
        }
    }
    endTime <- proc.time()

    list(
        sigmaSamp = sigmaSamp,
        betaSamp = betaSamp,
        alphaSamp = alphaSamp,
        time = endTime - startTime
    )
}
