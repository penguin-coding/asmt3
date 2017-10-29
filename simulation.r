# I confirm that the attached is my own work, except where clearly indicated
# in the text. 

source('BCaHelperFunctions.r') # Use Len's code to help with BCa Bootstrap

non.parametric.sample <- function(data, n){
  # purpose : produces n random samples of size 'size' from the supplied data
  #
  # inputs  : data - numeric vector of univariate observations
  #           n    - positive integer number of samples to be drawn
  #
  # output  : size*n dimension matrix. Each column is a generated sample.
  
  if (class(data)!='numeric') stop('input data must be numeric')
  if (n%%1!=0 | n<1) stop('n must be a positive integer')
  
  return(replicate(n, sample(x=data, size=length(data),replace=T)))
}

get.bca.alphas <- function(data, est, boot.est, alpha, func){
  # purpose : calculates values alpha1 and alpha2, to be taken as the quantiles
  #           used to compute a BCa confidence interval
  #
  # inputs  : data     - numeric vector of univariate observations 
  #           est      - estimate of statistic of interest from original sample
  #           boot.est - estimate of statistic of interest from each of the 
  #                      bootstrap samples
  #           alpha    - a value between 0 and 1 such that a (1-a)*100 %
  #                      confidence interval will be produced.
  #           func     - a function such that func(data) produces est. It is the
  #                      function which provides an estimate of the statistic of
  #                      interest given an input dataset. 
  #
  # output  : a named vector of alphas such that quantile(data, probs=alphas)
  #           returns the BCa confidence interval for the statistic of interest
  #
  # notes: - this method is treated as private, and so does not type check its
  #          inputs. It should only be called by functions which have already
  #          checked their inputs and produce consistent outputs of the correct
  #          format.
  
  alpha <-  alpha/2 # transform alpha s.t. we get a (1-2a)*100 % CI
  zhat0 <- get.zhat0(est, boot.est)
  ahat <- get.ahat(data, est, func)
  
  # Calculate the values of alpha1 and alpha2 according to the formula specified
  # in the assignment 3 outline pdf:
  alpha1 <- zhat0 + ( (zhat0 + qnorm(alpha)) / (1-ahat*(zhat0 + qnorm(alpha))) )
  alpha2 <- zhat0 + ((zhat0 + qnorm(1-alpha))/(1-ahat*(zhat0 + qnorm(1-alpha))))
  
  alpha1 <- pnorm(alpha1) ; alpha2 <- pnorm(alpha2)
  alphas <- c(alpha1,alpha2) ; names(alphas) <- c('alpha1','alpha2')
  
  return(alphas)
}

non.par.bootstrap <- function(data, n=999, alpha = 0.05, func = mean,
                           method = 'percentile'){
  # purpose : produces a 1 - alpha % confidence interval for the 'func' of the
  #           data, using a non-parametric bootstrap of n samples of size 'size'
  #
  # inputs  : data   - numeric vector of observations from a univariate 
  #                    distribution
  #           n      - the number of resamples to perform for the bootstrap
  #           alpha  - a number between 0 and 1, such that a 1-alpha %
  #                    confidence interval will be produced
  #           func   - function to be used to calculate the statistic of
  #                    interest for each sample. The default is 'mean'. 
  #           method - 'percentile' or 'bca'. Specifies the method to be used to 
  #                    obtain the quantiles to sample from to obtain the
  #                    interval
  #          
  # output  : named vector containing the lower and upper bounds of the interval
  
  # Input checks and setting default values: 
  if (class(data)!='numeric') stop('input data must be numeric')
  size <- length(data) 
  if (n%%1!=0 | n<1) stop('n must be a positive integer')
  if (alpha<0 | alpha>1) stop('alpha must be between 0 and 1')
  func <- match.fun(func) # to allow the user to pass in the func or func name
  if (!is.function(func)) stop('invalid function supplied as func argument')
  if (!method %in% c('percentile','BCa')) stop('invalid method')
  
  samples <- non.parametric.sample(data, n)       # generate the random samples
  samples <- cbind(samples,data)                  # add in the observed data
  samples <- apply(samples,2,func)                # calculate statistics
  
  lower <- alpha/2      # percentile method
  upper <- 1 - alpha/2  # intervals
  
  if (method=='BCa'){
    alphas <- get.bca.alphas(data, func(data), samples, alpha, func)
    lower <- alphas[1]
    upper <- alphas[2]
    }
  
  CI <- quantile(samples, probs=c(lower, upper))
  return(CI)
}

D <- rnorm(100)
non.par.bootstrap(D,method='BCa')