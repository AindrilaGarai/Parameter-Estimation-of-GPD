X <- rGPD(100,2,1)

S_data <- function(X)
{
  X <- sort(X)
  n <- length(X)
  j <- sample(1:n,1)
  sj <- X[j]
  rtn <- X/sj
  return(rtn)
}
S <- S_data(X)

neg_exp <- function(N)
{
  u <- runif(N,0,1)
  sam <- log(u)
  return(list(sam,
              exp(sam)))
}

intermediate <- function(S, u, k)
{
  rtn <- prod(1-(u*S))^((1/k)-1)
  return(rtn)
}

sampling <- function(S, n, k, N=100)
{
  if(k<0)
  {
    vec <- numeric(N)
    for (i in 1:N) 
    {
      sam <- neg_exp(1)[[1]]
      vec[i] <- intermediate(S, sam, k) * ((sam/k)^(n-1)) / (abs(k) * exp(sam))
    }
    return(sum(vec))
  }
  else
  {
    vec <- numeric(N)
    for (i in 1:N) 
    {
      sam <- runif(1,0,1/max(S))
      vec[i] <- intermediate(S, sam, k) * ((sam/k)^(n-1)) / (abs(k) * (1/sam))
    }
    return(sum(vec))
  }
  
}

likelihood <- function(k)
{
  n <- length(S)
  if(k == 0)
  {
    rtn <- (factorial(n) * factorial(n-1))/((sum(S))^n)
  }
  else
  {
    rtn <- factorial(n) * sampling(S, n, k)
  }
  rtn
}

k <- seq(-4,4,0.1)
vec <- numeric(81)
for (i in 1:81)
{
  vec[i] <- log(likelihood(k[i]))
}
vec

new_est <- function(S)
{
  maxLik(likelihood, start=1)#theta_ZJ(S)$`shape parameter`)
}
theta_ZJ(X)
new_est(S)


sigma <- function(X, k)
{
  if(k==0)
  {
    rtn <- mean(X)
  }
  if(k >= 1)
  {
    rtn <- k * max(X)
  }
  else
  {
    rtn <- mean(X)*(1/k)
  }
  return(rtn)
}
sigma(X, 1.131705)

sigma_bc <- function(sigma,k,n)
{
  if(k<0)
  {
    rtn <- sigma * (1-((3 - (5*k) - (4*(k^2)))/(n*(1-3*k))) )
  }
  if(k >= 1)
  {
    rtn <- sigma / (1- (factorial(n)*factorial(k)/factorial(n+k)))
  }
  else
  {
    rtn <- sigma
  }
  return(rtn)
}

library(Metrics)
library(maxLik)
library(ggplot2)

k <- seq(-4,4,0.5)

# simulation study for ZJ method by fixing scale parameter = 1

measure_new <- function(n, k, sigma, times)
{
  b <- matrix(0, times, length(k))
  k_est <- matrix(0, times, length(k))
  abs_bias <- numeric(length(k))
  rmse <- numeric(length(k))
  for (i in 1:length(k)) 
  {
    for (j in 1:times) 
    {
      X <- rGPD(n, k[i], sigma)
      S <- S_data(X)
      k_est[j,i] <- coef(new_est(S))
      b[j,i] <- abs(k_est[j,i]-k[i])
    }
    abs_bias[i] <- mean(b[,i])
    rmse[i] <- rmse(rep(k[i], times), k_est[,i])
  }
  return(list("Absolute Bias" = abs_bias, 
              "RMSE" = rmse))
}
dat_s_20 <- measure_new(20,k,1,50)
save(dat_s_20, file = "dat_s_20.Rdata")

dat_s_50 <- measure_new(50,k,1,100)
save(dat_s_50, file = "dat_s_50.Rdata")

dat_s_100 <- measure_new(100,k,1,100)
save(dat_s_100, file = "dat_s_100.Rdata")


dat_new$`Absolute Bias` <- (dat_new$`Absolute Bias` - min(dat_new$`Absolute Bias`))/(max(dat_new$`Absolute Bias`)-min(dat_new$`Absolute Bias`))

dat_n20_bias <- data.frame(k,dat_new$`Absolute Bias`) 
head(dat_n20_bias)

ggplot(dat_n20_bias, aes(k))+
  geom_line(aes(y=dat_new..Absolute.Bias.), color = "red")+
  ylab("Absolute Bias")+
  ylim(c(0,2.5))+
  theme_classic()

