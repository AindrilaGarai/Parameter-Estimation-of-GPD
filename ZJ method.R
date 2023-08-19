# ZJ method (Bayesian perspective)
# works well when k < 0.5

# log- likelihood of GPD

log_likelihood_GPD <- function(X, n, theta_est) 
{
  new_k <- -(1/n)*( sum( log( 1- ( theta_est*X )))) 
  new_sigma <- new_k / theta_est
  if( new_k == 0 )
  {
    rtn <- -(n*log(new_sigma)) - ((sum(X))/new_sigma)
  }
  else
  {
    rtn <- -(n*log(new_sigma)) +(((1/new_k)-1)*sum(log(1-((new_k*X)/new_sigma))))
  }
}


# an intermediate step to calculate thetas.

theta_j <- function(X, j, n, m, k_star, sigma_star)
{
  term_1 <- (n-1) / ((n+1)*max(X))
  term_2 <- sigma_star / k_star
  term_3 <- 1-(((j - 0.5)/m)^k_star)
  
  rtn <- term_1 - (term_2 * term_3)
  return(rtn)
}


# an intermediate step to calculate weights.

weight_j <- function(X, theta__k, theta__j)
{
  n <- length(X)
  l_theta_k <- sapply(theta__k, function(y) log_likelihood_GPD(X, n, y) )
  l_theta_j <- log_likelihood_GPD(X, n, theta__j)
  den <- sum(exp(l_theta_k - l_theta_j))
  return( 1/den )
}


# final code to get estimated scale and shape parameter using ZJ method.

theta_ZJ <- function(X) 
{
  n <- length(X)
  m <- 20 + floor(sqrt(n))
  data <- sort(X)
  p <- (3:9)/10
  xp <- data[round(n*(1-p)+.5)]
  xq <- data[round(n*(1-p*p)+.5)]
  k <- log(xq/xp-1,p)
  a <- k*xp/(1-p^k)
  sigma_star <- 1/(2*median(a))
  k_star <- -1
  
  weight <- numeric(m)
  theta <- numeric(m)
  
  for (j in 1:m) 
  {
    theta[j] <- theta_j(X, j, n, m, k_star, sigma_star)
  }
  for (j in 1:m) 
  {
    theta__k <- theta # vector
    theta__j <- theta[j] # a single value
    weight[j] <- weight_j(X, theta__k, theta__j)
  }
  output <- sum(weight * theta)
  k <- -mean(log(1-output*X))
  return(list("shape parameter" = k,
              "scale parameter" = k/output))
}




