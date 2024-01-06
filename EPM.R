GPD_cdf <- function(x, k, sigma=1)
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

f <- function(X, k)
{
  fu <- numeric(length(X))
  for (i in length((X))) 
  {
    fu[i] <- (i/length(X)) - GPD_cdf(X[i], k)
  }
  return(sum(fu))
}


bisection <- function(X, start, end, tol = 1e-5, max_iter=100)
{
  iter <- 1
  while (iter <= max_iter) 
  {
    iter = iter + 1
    c <- (start+end)/2
    if(f(X, c) == 0 && ((end-start)/2)<tol )
    {
      return(c)
    }
    if(sign(f(X, c))==sign(f(X, start)))
    {
      start <- c
    }
    else
    {
      end <- c
    }
  }
  return(c)
}

# Elemental Percentile Method (EPM)

delta_ij <- function(X, x_i, x_j, i, j) # x_in, x_jn are order statistics
{
  n <- length(X)
  p_i <- i/(n + 1)
  c_i <- log(1 - p_i)
  p_j <- j/(n + 1)
  c_j <- log(1 - p_j)
  d <- (c_j*x_i) - (c_i*x_j)
  if(d == 0) return(print("Inf")) 
  else
  {
    delta_0 <- (x_i*x_j*(c_j - c_i)) / d
    if(delta_0 > 0)
    {
      rtn <- bisection(X, x_j, delta_0)
      return(rtn)
    }
    else
    {
      rtn <- bisection(X, delta_0, 0)
      return(rtn)
    }
  }
}

k_ij <- function(x, i, n, delta) # x_in, delta_ij_est, c_i
{
  if (delta == "Inf"){return(0)}
  else
  {
    p_i <- i/(n + 1)
    c_i <- log(1 - p_i)
    rtn <- log(abs((1 - x)/delta)) / c_i
    return(rtn)
  }
}

sigma_ij <- function(k, delta)
{
  rtn <- k * delta
  return( rtn )
}

theta_EPM <- function(X)
{
  X <- sort(X)
  n <- length(X)
  theta_ij <- matrix(NaN, n, n)
  kij <- matrix(NaN, n, n)
  sig_ij <- matrix(NaN, n, n)
  
  for (i in 1:(n-1)) 
  {
    for (j in (i+1):n) 
    {
      theta_ij[i,j] <- delta_ij(X, X[i], X[j], i, j)
      kij[i,j] <- k_ij(X[i], i, n, theta_ij[i,j])
      sig_ij[i,j] <- sigma_ij(kij[i,j], theta_ij[i,j])
    }
  }
  kij <- na.omit(as.vector(kij))
  sig_ij <- na.omit(as.vector(sig_ij))
  k_epm <- median(kij)
  sig_epm <- median(sig_ij)
  return(list("k" = k_epm,
              "sigma" = sig_epm))
}
theta_EPM(X)

# simulation study for EPM 

measure_epm <- function(n, k, sigma, times)
{
  percentiles <- seq(90,10,-10)
  b <- matrix(0, times, length(k))
  k_est <- matrix(0, times, length(k))
  abs_bias <- numeric(length(k))
  mse <- numeric(length(k))
  for (i in 1:length(k)) 
  {
    for (j in 1:times) 
    {
      X <- rgpd(n, loc = 0, scale = sigma, shape = k[i])
      k_est[j,i] <- theta_EPM(X)# median(quantile(sort(X), probs = percentiles / 100))
      b[j,i] <- abs(k_est[j,i]-k[i])
    }
    abs_bias[i] <- mean(b[,i])
    mse[i] <- n*mse(rep(k[i], times), k_est[,i])
  }
  return(list("Absolute Bias" = abs_bias, 
              "MSE" = mse))
}
dat_epm <- measure_epm(20,k,1,5)
