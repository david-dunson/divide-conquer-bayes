##' Data augmentation for posterior sampling in finite mixture of Gaussians
##'
##' Sample from the posterior distribution of parameters in a mixture
##'     of L multivariate Gaussian distributions, where L is fixed.
##'     The means and covariances of L Gaussian can be different. The
##'     parameters in the model are L mean vectors, L covariance
##'     matrices, and L mixture component probabilities. See
##'     Bishop (2006) for details.
##' @title Gibbs sampler for mixture of multivariate Gaussian distributions
##' @param dataMat matrix with rows containing the observations and
##'     columns indexing the dimensions.
##' @param ncomp number of mixture components (L).
##' @param niter no. of sampling iterations.
##' @param nburn no. of samples to discard as burnins.
##' @param nthin no. of samples to discard during thinning.
##' @return list of posterior samples for mean and covariance
##'     parameters of 'ncomp' Gaussian distributions and  'ncomp'
##'     mixture component probabilities.
##' @author Cheng Li (stalic@nus.edu.sg) and Sanvesh Srivastava (sanvesh-srivastava@uiowa.edu)
mvnMix <- function (dataMat, ncomp = 2, niter = 10000, nburn = 5000, nthin = 5) {
    library(matrixStats)
    library(mvtnorm)
    library(MCMCpack)

    nobs <- nrow(dataMat)
    ndim <- ncol(dataMat)
    alpha <- rep(1 / ncomp, ncomp)

    probs <- as.numeric(rdirichlet(1, alpha))
    zMat <- rmultinom(nobs, 1, probs)
    densMat <- matrix(0.0, nrow = nobs, ncol = ncomp)
    muMat <- matrix(0.0, nrow = ncomp, ncol = ndim)
    sigArr <- aperm(array(diag(1.0, ndim), c(ndim, ndim, ncomp)), perm = c(3, 1, 2))

    kap0 <- 0.01; m0 <- rep(0.0, ndim) # mu hyper-pars
    nu0 <- 2; s0 <- 2 * nu0 * diag(1.0, ndim) # sigma hyper-pars

    probsSamp <- matrix(0.0, nrow = (niter - nburn) / nthin, ncol = ncomp)
    muMatSamp <- array(0.0, dim = c(ncomp, ndim, (niter - nburn) / nthin))
    sigMatSamp <- array(0.0, dim = c(ncomp, ndim, ndim, (niter - nburn) / nthin))

    cts <- 0
    startTime <- proc.time()
    for (ii in 0:niter) {
        idxList <- lapply(split(zMat, row(zMat)),
                          function (x) {
                              which(x == 1)
                          })
        ns <- sapply(idxList, length)
        # sample probs
        probs <- (rdirichlet(1, ns + alpha))

        if (ii %% 100 == 0) {cat("gibbs: ", ii, "\n"); cat("ns: ", ns, "\n")}

        for (jj in seq_along(ns)) {
            datMean <- colMeans(dataMat[idxList[[jj]], , drop = FALSE])
            ## mean
            muCov <- sigArr[jj, , ] / (kap0 + ns[jj])
            muMean <- (kap0 * m0 + ns[jj] * datMean) / (kap0 + ns[jj])
            muMat[jj, ] <- rmvnorm(1, mean = muMean, sigma = muCov, method = "chol")
            ## cov
            mat1 <- (kap0 * ns[jj] /  (kap0 + ns[jj])) * tcrossprod(datMean - m0, datMean - m0)
            centMat <- dataMat[idxList[[jj]], , drop = FALSE] - matrix(datMean, nrow = length(idxList[[jj]]), ncol = ndim, byrow = TRUE)
            mat2 <- crossprod(centMat, centMat)
            covSclMat <- mat1 + mat2 + s0
            covDf <- ns[jj] + nu0 + 1
            sigArr[jj, , ]<- riwish(covDf, covSclMat)
            densMat[ , jj] <- dmvnorm(dataMat, muMat[jj, ], sigArr[jj, , ], log = TRUE)
        }

        lprobs <- densMat + log(matrix(probs, nrow = nobs, ncol = ncol(densMat), byrow = TRUE))
        eprobs <- exp(lprobs - rowMaxs(lprobs)) / rowSums(exp(lprobs - rowMaxs(lprobs)))

        for (kk in seq_len(nobs)) {
            ppp <- eprobs[kk, ]
            zMat[ , kk] <- rmultinom(1, 1, ppp)
        }

        if ((ii > nburn) && (ii %% nthin == 0)) {
            cts <- cts + 1
            probsSamp[cts, ] <- probs
            muMatSamp[ , , cts] <- muMat
            sigMatSamp[ , , , cts] <- sigArr
        }
    }
    endTime <- proc.time()


    list(
        'mu' = muMatSamp,
        'cov' = sigMatSamp,
        'prob' = probsSamp,
        'time' = endTime[3] - startTime[3]
    )
}

## rm(list=ls())
## setwd("~/Dropbox/code/")

## genData <- function (muList = list('1' = c(1, 2), '2' = c(4, 5)),
##                      sigMat = matrix(c(1, 0.5, 0.5, 2), 2, 2),
##                      probs = c(0.3, 0.7),
##                      nobs = 1000) {
##     library(mvtnorm)
##     zs <- rmultinom(nobs, 1, probs)

##     idxList <- lapply(split(zs, row(zs)),
##                       function (x) {
##                           which(x == 1)
##                       })



##     dataMat <- matrix(0.0, nrow = nobs, ncol = 2)
##     clusts <- numeric(nobs)
##     for (ii in seq_along(muList)) {
##         idx <- idxList[[ii]]
##         clusts[idx] <- ii
##         dataMat[idx, ] <- rmvnorm(length(idx), muList[[ii]], sigMat)
##     }
##     rownames(dataMat) <- paste("data", seq_len(nrow(dataMat)), clusts, sep = "_")

##     list(data = dataMat,
##          cluster = clusts,
##          zs = zs
##          )
## }

## set.seed(12345)
## res <- genData(nobs=1000)

## train <- res$data

## # consistency check
## full <- mvnMix(train, ncomp = 2, niter = 10000, nburn = 5000, nthin = 5)
## plot(ts(t(full$mu[1, , ])))
## plot(ts(t(full$mu[2, , ])))
## plot(rbind(t(full$mu[1, , ]), t(full$mu[2, , ])))

## fullCorr <- list(numeric(1000), numeric(1000))
## for (ss in 1:1000) {
##     fullCorr[[1]][ss] <- cov2cor(full$cov[1, , , ss])[1, 2]
##     fullCorr[[2]][ss] <- cov2cor(full$cov[2, , , ss])[1, 2]
## }

## plot(fullCorr[[1]], ylim = range(unlist(fullCorr)))
## points(fullCorr[[2]], col = "red")
