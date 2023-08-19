# Methods of Moments (MOM) 
# for k > -1

k_MOM <- function(X)
{
  x_mean <- mean(X) # sample mean
  x_var <- var(X) # sample variance
  k <- (( x_mean^2 )/( x_var-1 ))/2
  return(k)
}

sigma_MOM <- function(X)
{
  x_mean <- mean(X)
  x_var <- var(X)
  sigma <- (x_mean * (( x_mean^2 )/( x_var+1 )))/2
  return(sigma)
}

# Probability-Weighted Moments (PWM)
# k > -0.5

u_calc <- function(X)
{
  n <- length(X)
  X <- sort(X)
  rtn <- 0
  for (i in 1:n) 
  {
    val <- ((n-i)/(n-1))/X[i]
    rtn <- rtn + val
  }
  return(rtn)
}

k_PWM <- function(X)
{
  x_mean <- mean(X)
  u <- u_calc(X)
  k <- (x_mean/( x_mean - (2*u))) - 2
  return(k)
}

sigma <- function(X)
{
  x_mean <- mean(X)
  u <- u_calc(X)
  sigma <- (2 * x_mean * u)/(x_mean - ( 2*u ))
  return(sigma)
}