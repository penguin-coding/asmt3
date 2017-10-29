gy <- function(x,mu1,sd1,mu2,sd2,lambda1){
  # purpose : Evaluates the pdf of the fitted bimodal distribution, given 
  #           the estimated parameters. 
  # inputs  : x       - the scalar or vector at which the function must be
  #                     evaluated
  #           mu1     - the mean of the first component of the mixture
  #           sd1     - the standard deviation of the first component of the
  #                     mixture
  #           mu2     - the mean of the second component of the mixture
  #           sd2     - the standard deviation of the second component of the
  #                     mixture
  #           lambda1 - the proportion of the modelled population which arises
  #                     from the first component of the mixture
  # outputs : A scalar or vector quantity ; the density at x. 
  
  #Calculates the likelihood. 
  return(lambda1*dnorm(x,mu1,sd1) + (1-lambda1)*dnorm(x,mu2,sd2))
}