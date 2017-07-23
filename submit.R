cmdArgs <- commandArgs(trailingOnly = TRUE)

mtd <- as.numeric(cmdArgs[1])
id <- as.numeric(cmdArgs[2])

if (mtd == 1) {
    source("mcmc_lme_sampler.R")
    cvs <- rep(1:10, each = 2) # indexes the replications
    ndims <- rep(1:2, times = 10) # indexes the two dimension in every replication

    cid <- cvs[id]
    did <- ndims[id]

    cvtrain <- readRDS("data/mixed.rds") # change appropriately to import the data
    train <- cvtrain[[cid]][[did]] # extract data

    res <- sampleFromMixMdl(train$y, train$x, train$z, train$group, 10000, 5000, 5, id)
    fname <- paste0("result/mcmc/mcmc_lme_cv_", cid, "_p_", did, ".rds")
    saveRDS(res, fname)
} else if (mtd == 4) {
    source("sub_lme_sampler.R")
    cvs <- rep(1:10, each = 2)  # indexes the replications
    ndims <- rep(1:2, times = 10) # indexes the two dimension in every replication

    tmp <- cbind(cvs, ndims)
    # 10 indexes the number of subsets and 20 indexes the number of
    # replications for the two dimensions
    wids <- cbind(tmp[rep(1:nrow(tmp), each = 10), ], rep(1:10, times = 20))

    cid <- wids[id, 1]
    did <- wids[id, 2]
    sid <- wids[id, 3]

    cvtrain <- readRDS(paste0("data/sub_mixed_cv_", cid,"_k10", ".rds"))
    train <- cvtrain[[did]][[sid]]

    res <- sampleFromSubMixMdl(train$y, train$x, train$z, train$group, train$nrep, 10000, 5000, 5, id)
    fname <- paste0("result/sub/samp/sub_mixed_cv_", cid, "_p_", did, "_k_", sid, "_nsub10.rds")
    saveRDS(res, fname)
} else if (mtd == 5) {
    source("sub_lme_sampler.R")
    cvs <- rep(1:10, each = 2) # indexes the replications
    ndims <- rep(1:2, times = 10) # indexes the two dimension in every replication

    tmp <- cbind(cvs, ndims)
    # first 20 indexes the number of subsets and the next 20 indexes the number of
    # replications for the two dimensions
    wids <- cbind(tmp[rep(1:nrow(tmp), each = 20), ], rep(1:20, times = 20))

    cid <- wids[id, 1]
    did <- wids[id, 2]
    sid <- wids[id, 3]

    cvtrain <- readRDS(paste0("data/sub_mixed_cv_", cid,"_k20", ".rds"))
    train <- cvtrain[[did]][[sid]]

    res <- sampleFromSubMixMdl(train$y, train$x, train$z, train$group, train$nrep, 10000, 5000, 5, id)
    fname <- paste0("result/sub/samp/sub_mixed_cv_", cid, "_p_", did, "_k_", sid, "_nsub20.rds")
    saveRDS(res, fname)
} else {
    ## implement the combination using PIE algorithm after importing
    ## the sub*** files from result/sub/samp directory.
}
