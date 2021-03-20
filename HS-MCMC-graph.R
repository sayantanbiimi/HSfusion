# MCMC for Bayesian fusion with Horseshoe prior for graph de-noising
# Author: Sayantan Banerjee
# Date created: Mar 19, 2021
# Last modified on: Mar 19, 2021


##########################################################################
# Load required packages
require(igraph)
require(MASS)

# Previous steps:
# Step 1: Generate the linear chain graph and the true signal (gen-linear-chain.R)
# Step 2: Generate the data
# Step 3: Get edge incidence matrix
# Step 4: MCMC (this code)

# Arguments required
# n -- number of nodes
# theta_0 -- true signal (for computing MSE)
# y -- data

denoise_MCMC = function(n = 100, theta_0, y, g_chain, E_T, root){
  
  # Specifying hyperparameters in prior distribution
  
  a_s = b_s = 0.5 # sigma2 parameters
  lambda = rep(1,nrow(E_T))
  
  lambda_sq0 = 5 # fixed hyperparameter for theta[root]
  
  lambda_sq = lambda^2
  
  ## Posterior sampling
  
  # Initialize
  
  n.iter = 1000 # number of MCMC samples
  burnin = 100 # burn-in samples
  
  theta = y
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
    
    # Sampling theta[root]
    
    edge.set <- which(E_T[ , root] != 0)
    nei.set <- rep(0,length(edge.set))
    for(l in 1:length(edge.set)){
      nei.set[l] <- setdiff(which(E_T[edge.set[l], ] != 0), root)
    }
    
    zeta[root] = sigma2/(1 + sum(1/(lambda_sq[edge.set]*tau2)) +
                           1/lambda_sq0)
    mu[root] = zeta[root]*(y[root] + 
                             sum(theta[nei.set]/(lambda_sq[edge.set]*tau2)))/sigma2
    
    theta[root] = rnorm(1, mean = mu[root], sd = sqrt(zeta[root]))
    
    #cat("NAN detected:", sum(is.nan(zeta)),"\n")
    
    # Sampling rest of the theta s
    
    for(i in setdiff(1:n,root)){
      edge.set <- which(E_T[ , i] != 0)
      nei.set <- rep(0,length(edge.set))
      for(l in 1:length(edge.set)){
        nei.set[l] <- setdiff(which(E_T[edge.set[l], ] != 0), i)
      }
      
      zeta[i] = sigma2/(1 + sum(1/(lambda_sq[edge.set]*tau2)))
      
      # cat("NAN detected:", sum(is.nan(zeta)),"\n")
      mu[i] = zeta[i]*(y[i] + 
                         sum(theta[nei.set]/(lambda_sq[edge.set]*tau2)))/sigma2
      
      theta[i] = rnorm(1, mean = mu[i], sd = sqrt(zeta[i]))
    }
    
    
    
    keep.theta[iter, ] = theta
    
    ##### Sampling lambda squares and nu s####
    
    K <- nrow(E_T)
    for(k in 1:K){
      id <- ends(g_chain,k)
      lambda_sq[k] = 1/rgamma(1, 1, 1/nu[k] + 
                                (theta[id[1]] - theta[id[2]])^2/(2*sigma2*tau2))
      nu[k] = 1/rgamma(1, 1, 1 + 1/lambda_sq[k])
      
      
      #cat("NAN detected:", sum(is.nan(lambda_sq)),"\n")
      
      #cat("NAN detected:", sum(is.nan(nu)),"\n")
      keep.lambda_sq[iter, k] <- lambda_sq[k]
      keep.nu[iter, k] <- nu[k]
      
    }
    
    ##### Simulating sigma2 ####
    
    T <- (E_T%*%theta)*lambda_sq^(-1/2)
    param1 = n + a_s
    param2 = b_s + theta[root]^2/(2*lambda_sq0) + 0.5*sum((y-theta[1:n])^2) +
      0.5*c(t(T)%*%T)/tau2
    
    sigma2 = 1/rgamma(1, param1, param2)
    # cat("NAN detected:", sum(is.nan(sigma2)),"\n")
    
    
    keep.sigma2[iter] = sigma2 
    
    ##### Simulating tau2 ####
    
    tau2 = 1/rgamma(1, n/2, 1/xi + c(t(T)%*%T)/(2*sigma2))
    #cat("NAN detected:", sum(is.nan(tau2)),"\n")
    
    keep.tau2[iter] = tau2
    
    ##### Simulating xi ####
    
    xi = 1/rgamma(1, 1, 1 + 1/tau2)
    #cat("NAN detected:", sum(is.nan(xi)),"\n")
    
    keep.xi[iter] = xi
    
  }
  
  out <- NULL
  out$keep.theta <- keep.theta
  out$burnin <- burnin
  out$n.iter <- n.iter
  return(out)
}


