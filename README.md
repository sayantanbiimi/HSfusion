# HSfusion
This repository hosts the codes to implement Horseshoe shrinkage based methods for Bayesian fusion estimation. The original paper can be found in arXiv at https://arxiv.org/abs/2102.07378

The file 'gen-linear-chain.R' generates piecewise constant signals. Users can change the dimension of the signal along with their patterns and values. The 'HS-MCMC.R' file implements the Bayesian fusion estimation via Horseshoe shrinkage.

The 'HS-MCMC-graph.R' file implements the Bayesian fusion estimation via Horseshoe shrinkage for the graph signal de-noising problem. Users can start with an arbitrary graph G, and then find the DFS-chain graph along with the corresponding edge incidence matrix using the source code 'get-DFS.R'.
