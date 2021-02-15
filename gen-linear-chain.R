# Function to generate signals on linear chain graphs
# Author: Sayantan Banerjee
# Date created: Jan 07, 2021
# Last modified on: Jan 07, 2021

##################################################################################
# Load required packages
require(igraph)
require(MASS)

##################################################################################
# Function: gen_linear_chain
# Arguments: signal type
# Output: g       -- a linear chain graph with 100 nodes and 10 components
#         theta_0 -- true signal over g
# Details: For 'type' = "even", each of the 10 components are evenly distributed
#          with component length 10
#          For 'type' = "uneven", the smaller components have length 5
#          For 'type' = "veryuneven", the smaller components have length 2
##################################################################################


gen_linear_chain = function(type = "even"){
  
  # Step 1: Construct a linear chain graph G
  n = 100 # number of nodes
  K = 10  # number of components
  g = make_graph(c(1,rep(2:(n-1),each=2),n), directed = F)
  
  comp = NULL
  
  if(type == "even"){
    comp_len = c(rep(n/K, K))
  } else if(type == "uneven"){
    comp_len = rep(c(5,15),5)
  } else if(type == "veryuneven"){
    comp_len = rep(c(2,18),5)
  }
  cusum_comp_len = c(0,cumsum(comp_len))
  for(i in 2:(K+1)){
    start = cusum_comp_len[i-1] + 1
    end = cusum_comp_len[i]
    comp[[i-1]] = c(start:end)
  }
  
  # Step 3: Signal generation
  theta_unique = c(1,0,-2,0,3,0,-4,0,4,0)
  theta_0 = rep(0,n)
  for(i in 1:K){
    theta_0[comp[[i]]] <- theta_unique[i]
  }
  out = NULL
  out$theta_0 = theta_0
  out$g = g
  return(out)
}

