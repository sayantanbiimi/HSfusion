# Function to generate DFS chain graph and the corresponding
# edge incidence matrix from an input graph

#####################################################################
# Function : get_DFS
# Arguments: g -- input graph
#            root -- root node for DFS
# Value: g_chain -- DFS chain graph
#        E_T -- edge incidence matrix of the DFS chain
#        root -- root node for DFS  
#####################################################################

require(igraph)

get_DFS = function(g, root){
  g_dfs = as.vector(dfs(g, root = root)$order)
  n = length(g_dfs)
  edgelist = c(g_dfs[1], rep(g_dfs[2:(n-1)], each = 2), g_dfs[n])
  g_chain = make_graph(edgelist, directed = F)

  K = length(E(g_chain))
  E_T = matrix(0,K,n)
  for(i in 1:K){
    E_T[i, ends(g_chain,i)[ , 1]] = -1
    E_T[i, ends(g_chain,i)[ , 2]] = 1
  }
  out = NULL
  out$g_chain = g_chain
  out$E_T = E_T
  out$root = root
  return(out)
}

