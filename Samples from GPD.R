# samples from GPD with shape parameter = k, 
#                       scale parameter = σ,
#                       no of samples = n.
# inverse-transform method.

GPD <- function(k, sigma) # single observation from GPD
{
  u <- runif(1,0,1)
  if(k == 0)
  {
    X <- -sigma * (log(1-u))
  }
  else
  {
    X <- (sigma/k) * (1- ((1-u)^k))
  }
  return(X)
}

rGPD <- function(n, k, sigma) # n number of observations
{
  vec <- numeric(length = n)
  for (i in 1:n) 
  {
    vec[i] <- GPD(k, sigma)
  }
  return(vec)
}