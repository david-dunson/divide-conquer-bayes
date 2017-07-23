* Files:
  - Files with names mcmc_XXX_sampler.R, where XXX = {lme, variable_selection, mix}, represent the MCMC sampler that uses the full data for sampling.
  - Files with names sub_XXX_sampler.R, where XXX = {lme, variable_selection, mix}, represent the MCMC sampler modified using stochastic approximation.
  - Subset posterior samples of parameters are obtained using sub_XXX_sampler.R.
  - The two samplers with XXX = lme are for linear mixed-effects models.
  - The two samplers with XXX = variable_selection are for variable selection using the GDP prior in linear regression models.
  - The two samplers with XXX = mix are for fitting mixture of multivariate Gaussians, where the number of mixture components is known apriori.
  - File pie_sampler.R includes the 'pie' function that that implements the univariate and multivariate versions of the PIE algorithm for combining parameter samples obtained from the $k$ subset posterior distributions. 
  - File submit.R shows how to implement any distributed Bayesian sampling algorithm for linear mixed-effects model on an SGE cluster.
  
* Comments:
  - A generic approach to implement any distributed Bayesian sampling algorithm is as follows:
    1. Randomly split the samples into $k$ subsets.
    2. Use sub_XXX_sampler.R files to run $k$ subset posterior sampling algorithms in parallel (on a cluster or across different threads of a multicore processor).
    3. Store the posterior samples from every subset.
    4. Import the subset posterior samples and use the 'pie' function in pie_sampler.R file to combine the collection of $k$ subset posterior samples.
  - Our future work seeks to combine steps 3. and 4. into a single step.
  - The four steps outlined previously are easily implemented on a cluster. The submit.R file shows how to do this on a SGE cluster through a qsub files.  

* Citations:
  1. Li, C., Srivastava, S., & Dunson, D. B. (2017). Simple, scalable and accurate posterior interval estimation. Biometrika, online print asx033.
  2. Srivastava, S., Cevher, V., Tranh-Dinh, Q., and Dunson, D.B. (2015). WASP: Scalable Bayes via barycenters of subset posteriors. In Artificial Intelligence and Statistics (pp. 912-920).
  3. Srivastava, S., Li, C., & Dunson, D. B. (2015). Scalable Bayes via barycenter in Wasserstein space. arXiv preprint arXiv:1508.05880.

* Contact:
  - Please contact Cheng Li <stalic@nus.edu.sg> and/or Sanvesh Srivastava <sanvesh-srivastava@uiowa.edu> if you have any questions related to the code.
