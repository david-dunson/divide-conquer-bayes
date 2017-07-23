##' MCMC sampler for subset posterior distributions for linear
##' mixed-effects model
##'
##' Sample from the posterior distribution of fixed effects and
##' covariance matrix of random effects in linear mixed-effects models
##' after stochastic approximation. This modifies an existing
##' Stan-based MCMC sampler by raising the likelihood to the power of
##' 'nrep', which is $k$. This is done using the 'stoc_approx_log' in
##' Stan; see the 'sub_lme.stan' file for details.
##' @title Stochastic approximation based MCMC sampler for linear mixed-effects model
##' @param yvec response.
##' @param xmat design matrix for the fixed effects.
##' @param zmat design matrix for the random effects.
##' @param group the no. of samples.
##' @param nrep the exponent in stochastic approx.
##' @param niter no. of sampling iterations.
##' @param nburn no. of samples to discard as burnins.
##' @param nthin no. of samples to discard during thinning.
##' @param id any natural number (used in set.seed)
##' @param subNo which data subset? (used in set.seed)
##' @return list of posterior samples for fixed effects and covariance matrix of random effects.
##' @author Cheng Li (stalic@nus.edu.sg) and Sanvesh Srivastava (sanvesh-srivastava@uiowa.edu)
sampleFromSubMixMdl <- function (yvec, xmat, zmat, group, nrep, niter, nburn, nthin, id, subNo) {
    library(inline)
    library(Rcpp)
    library(rstan)

    grpLbl <- sort(unique(group))
    ngroup <- length(grpLbl)

    ranefList <- list()
    grpIdx <- list()
    for (ii in 1:ngroup) {
        grpIdx[[ii]] <- which(group == grpLbl[ii])
        ranefList[[ii]] <- zmat[grpIdx[[ii]], , drop = FALSE]
    }
    ranefMat <- do.call(rbind, ranefList)
    fixefMat <- xmat[unlist(grpIdx), ]

    pos2 <- cumsum(sapply(grpIdx, length))
    pos1 <- c(1, pos2[-ngroup] + 1)

    ordY <- yvec[unlist(grpIdx)]
    ordGrp <- group[unlist(grpIdx)]

    idx <- seq_along(unlist(grpIdx))
    simList = list(nobs = length(yvec[idx]),
                   nfixef = ncol(xmat),
                   nranef = ncol(zmat),
                   ngroup = length(unique(ordGrp[idx])),
                   nrep = nrep,
                   xmat = fixefMat[idx, ],
                   zmat = ranefMat[idx, ],
                   group = ordGrp[idx],
                   yvec = ordY[idx],
                   pos1 = pos1,
                   pos2 = pos2)

    seeds <- (1:1000) * as.numeric(gsub(":", "", substr(Sys.time(), 12, 19)))

    stanCode <- readChar("sub_lme.stan", file.info("sub_lme.stan")$size)
    startTime <- proc.time()
    mdl <- stan(model_code = stanCode, data = simList, iter = niter, warmup = nburn, chains = 1, thin = nthin, seed = seeds[id * subNo])
    endTime <- proc.time()

    lst <- mdl@sim$samples[[1]]
    bs <- grep("fixef|covRanef", names(lst))
    sampdf <- do.call(cbind, lst[bs])

    list(samples = sampdf[(nrow(sampdf) - (niter - nburn) / nthin + 1):nrow(sampdf), ],
         time = endTime - startTime)
}
