##' MCMC sampler for linear mixed-effects model
##'
##' Sample from the posterior distribution of fixed effects and
##' covariance matrix of random effects in linear mixed-effects models
##'
##' MCMC sampler for linear mixed-effects model implemented using
##'     Stan. This sampler is based on the materials available online,
##'     including Stan manual. See the 'full_lme.stan' file for
##'     details.
##' @title MCMC sampler for linear mixed-effects model
##' @param yvec response.
##' @param xmat design matrix for the fixed effects.
##' @param zmat design matrix for the random effects.
##' @param group the no. of samples.
##' @param niter no. of sampling iterations.
##' @param nburn no. of samples to discard as burnins.
##' @param nthin no. of samples to discard during thinning.
##' @param id any natural number (used in set.seed).
##' @return list of posterior samples for fixed effects and covariance matrix of random effects.
##' @author Cheng Li (stalic@nus.edu.sg) and Sanvesh Srivastava (sanvesh-srivastava@uiowa.edu)
sampleFromMixMdl <- function (yvec, xmat, zmat, group, niter, nburn, nthin, id) {
    library(inline)
    library(Rcpp)
    library(rstan)

    gg <- ordered(as.character(group), levels = sort(unique(group)))
    group <- as.integer(gg)

    simList = list(
        nobs = length(yvec),
        nfixef = ncol(xmat),
        nranef = ncol(zmat),
        ngroup = length(unique(group)),
        xmat = xmat,
        zmat = zmat,
        group = group,
        yvec = yvec)

    seeds <- (1:1000) * as.numeric(gsub(":", "", substr(Sys.time(), 12, 19)))

    stanCode <- readChar("full_lme.stan", file.info("full_lme.stan")$size)
    startTime <- proc.time()
    mdl <- stan(model_code = stanCode, data = simList, iter = niter, warmup = nburn, chains = 1, thin = nthin, seed = seeds[id])
    endTime <- proc.time()

    lst <- mdl@sim$samples[[1]]
    bs <- grep("fixef|covRanef", names(lst))
    sampdf <- do.call(cbind, lst[bs])

    list(samples = sampdf[(nrow(sampdf) - (niter - nburn) / nthin + 1):nrow(sampdf), ],
         time = endTime - startTime)
}
