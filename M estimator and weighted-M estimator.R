# Tukey biweight function

rho_c <- function(u, c)
{
  if (abs(u) > c) 
  {
    rtn <- (c^2)/6
  }
  else
  {
    rtn <- ((u^2)/2) * ( 1 - ((u/c)^2) + (((u/c)^4)/3) )
  }
  return(rtn)
}

GPD_cdf <- function(x, k, sigma)
{
  if(k != 0)
  {
    cdf <- 1 - ((1-((k*x)/sigma))^(1/k))
  }
  else
  {
    cdf <- 1 - exp(-(x/sigma))
  }
  return(cdf)
}


# function needs to be minimized for M-estimators

M_est_func <- function(X, theta)
{
  n <- length(X)
  x <- sort(X)
  k <- theta[1]
  sigma <- theta[2]
  if (sigma > 0)
  {
    if (k <= 0 | (sigma/k) >= max(X))
    {
      r <- ((1:n-0.5)/n)- GPD_cdf(x, k, sigma)
      return(sum (sapply(r, rho_c, c=4.6851, simplify=T)))
    }
    else return(10^10)
  }
  else return(10^10)
}


# function needs to be minimized for weighted M-estimators

wei_M_est_func <- function(X, theta)
{
  n <- length(X)
  x <- sort(X)
  k <- theta[1]
  sigma <- theta[2]
  if (sigma > 0)
  {
    if (k <= 0 | (sigma/k) >= max(X))
    {
      cdf <- GPD_cdf(x, k, sigma)
      u <- (((1:n-0.5)/n)- cdf) / (sqrt(cdf * (1-cdf)))
      return(sum (sapply(u, rho_c, c=4.6851, simplify=T)))
    }
    else return(10^10)
  }
  else return(10^10)
}


# the M-estimators 
# Estimates by ZJ method are used as the initial values

M_est <- function(X)
{
  optim(theta_ZJ(X), M_est_func, data = X)$par
}

# the M-estimators 
# Estimates by M-estimator are used as the initial values

f.wm <- function(X)
{
  optim(M_est(X), wei_M_est_func, data = X)$par
}
