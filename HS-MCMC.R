# MCMC for Bayesian fusion with Horseshoe prior
# Simulations with linear chain graphs
# Author: Sayantan Banerjee
# Date created: Dec 28, 2020
# Last modified on: Dec 28, 2020


##########################################################################
# Load required packages
require(igraph)
require(MASS)

# Previous steps:
# Step 1: Generate the linear chain graph and the true signal (gen-linear-chain.R)
# Step 2: Generate the data
# Step 3: MCMC (this code)

# Arguments required
# n -- number of nodes
# theta_0 -- true signal (for computing MSE)
# y -- data

denoise_MCMC = function(n = 100, theta_0, y){
  
  # Specifying hyperparameters in prior distribution
  
  a_s = b_s = 0.5 # sigma2 parameters
  lambda = rep(1,n+1)
  
  lambda[1] = sqrt(5) # fixed hyperparameter for theta[1]
  lambda[n+1] = Inf
  
  lambda_sq = lambda^2

  ## Posterior sampling
  
  # Initialize
  
  n.iter = 2000 # number of MCMC samples
  burnin = 200 # burn-in samples
  
  theta = c(y,0)
  mu = rep(0,n)
  zeta = rep(1,n)
  sigma2 = var(y)
  tau2 = 0.01
  nu = rep(1,n)
  xi = 1
 
  ## Store the results
  
  keep.theta = matrix(0, n.iter, n)
  keep.sigma2 = rep(0, n.iter)
  keep.lambda_sq = matrix(0, n.iter,n-1)
  keep.tau2 = rep(0,n.iter)
  keep.nu = matrix(0, n.iter, n-1)
  keep.xi = rep(0,n.iter)
  
  # MCMC-time!  
  
  for(iter in 1:n.iter){
    
  cat("Running iteration:",iter,"\n")
  
  ##### Sampling theta s ####
  
  # Sampling theta[1]
  
  zeta[1] = sigma2/(1 + 1/(lambda_sq[2]*tau2) + 1/lambda_sq[1])
  mu[1] = zeta[1]*(y[1] + theta[2]/(lambda_sq[2]*tau2))/sigma2
  
  theta[1] = rnorm(1, mean = mu[1], sd = sqrt(zeta[1]))
    
  # Sampling rest of the theta s
  
  for(i in 2:n){
    zeta[i] = sigma2/(1 + 1/(lambda_sq[i+1]*tau2) + 1/(lambda_sq[i]*tau2))
    
    mu[i] = zeta[i]*(y[i] + theta[i+1]/(lambda_sq[i+1]*tau2) +
                       theta[i-1]/(lambda_sq[i]*tau2))/sigma2
    
    theta[i] = rnorm(1, mean = mu[i], sd = sqrt(zeta[i]))
  }
  
  
  
  keep.theta[iter, ] = theta[1:n]
  
  ##### Sampling lambda squares and nu s####
  
  for(i in 2:n){
    lambda_sq[i] = 1/rgamma(1, 1, 1/nu[i] + (theta[i] - theta[i-1])^2/(2*sigma2*tau2))
    nu[i] = 1/rgamma(1, 1, 1 + 1/lambda_sq[i])
  }
  
  keep.lambda_sq[iter, ] = lambda_sq[2:n]
  keep.nu[iter, ] = nu[2:n]
  
  ##### Simulating sigma2 ####
  
  param1 = n + a_s
  param2 = b_s + (sum((y-theta[1:n])^2) + sum( (diff(theta[1:n]))^2/lambda_sq[2:n] )/tau2 + 
                                    theta[1]^2/lambda_sq[1])/2
  
  sigma2 = 1/rgamma(1, param1, param2)
  keep.sigma2[iter] = sigma2 
  
  ##### Simulating tau2 ####
  
  tau2 = 1/rgamma(1, n/2, 1/xi + sum( (diff(theta[1:n]))^2/lambda_sq[2:n])/(2*sigma2))
  #cat("NAN detected:", sum(is.nan(tau2)),"\n")
  
  keep.tau2[iter] = tau2
  
  ##### Simulating xi ####
  
  xi = 1/rgamma(1, 1, 1 + 1/tau2)
  
  keep.xi[iter] = xi
  
  }
  
  bayes.theta <- apply(keep.theta[(burnin+1):n.iter, ], 2, mean)
  bayes.theta.q <- apply(keep.theta[(burnin+1):n.iter, ], 2,
                         quantile, probs = c(0.025,0.975))
  sigma.theta <- mean(keep.sigma2[(burnin+1):n.iter])
  MSE <- sum((bayes.theta - theta_0)^2)/n
  adj.MSE <- sum((bayes.theta - theta_0)^2)/sum(theta_0^2)
  
  omega = omega.hat = blockmat = matrix(0,n,n)
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      omega[i,j] = sum(theta_0[i] == theta_0[j])
      omega.hat[i,j] = sum(abs(bayes.theta[i] - bayes.theta[j]) <
                             abs( 0.5*(y[i] - y[j]) )  )
      blockmat[i,j] = abs(bayes.theta[i] - bayes.theta[j])
    }
  }
  cover = rep(0,n)
  for(i in 1:n){
    cover[i] = sum(theta_0[i] > bayes.theta.q[1, i] & theta_0[i] < bayes.theta.q[2,i])
  }
  
  R = sum(abs(omega.hat - omega))
  W = sum(blockmat*omega)/sum(omega)
  B = min(blockmat[blockmat*(1 - omega) > 0])
  cover.miss = sum(cover == 0)/n
  
  out <- NULL
  out$bayes.theta <- bayes.theta
  out$bayes.theta.q <- bayes.theta.q
  out$R <- R
  out$W <- W
  out$B <- B
  out$MSE <- MSE
  out$adj.MSE <- adj.MSE
  out$cover.miss <- cover.miss
  
  return(out)
}


