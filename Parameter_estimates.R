# Script for generating processed data tables from a fitted HMSC-hpc model
# Make sure the libraries are loaded from the right directory

.libPaths(c(.libPaths(), "/cluster/projects/nn9725k/eivind/R_packages/4.5"))

# Load libraries
library(coda)
library(Hmsc)
library(jsonify)
set.seed(666) # Remember to set seed 


# If the fitted object has already been stored, simply load from here
m <- readRDS("fitted_m.rds")


# Genera post estimates for visualization, check that only estimates that are included in the model is called. 
# For example, Sigma and Alpha are only estimated in spatial models, Rho is estimated in models with a phylogenetic component
# Calling a parameter that is not estimated throws and error and exits the script

postBeta <- getPostEstimate(m, "Beta")
postGamma <- getPostEstimate(m, "Gamma")
etaPost <- getPostEstimate(m, "Eta")
lambdaPost <- getPostEstimate(m, "Lambda")
alphaPost <- getPostEstimate(m, "Alpha")
sigmaPost <- getPostEstimate(m, "Sigma")
#rhoPost <- getPostEstimate(m, "Rho")
omegaPost <- getPostEstimate(m, "Omega", thin = 5, start = 100) # Omega is often very computationally intensive, thin and start later in the chain

# Include the post estimates you want it this list
Post_estimates <- list( "postBeta" = postBeta,
                        "postGamma" = postGamma,
                        "postEta" = etaPost,
                        "postLambda" = lambdaPost,
                        "postSigma" = sigmaPost,
                        "postAlpha" = alphaPost
                        "postOmega" = omegaPost)

saveRDS(Post_estimates, "Post_estimates.rds")

# Omega correlations are very heavy to compute so estimate separately
# The 'thin = ' argument is telling the function to only evaluate every n-samples when computing species-to-species associations
# This function is the slowest and most computationally intensive so fina a tradeoff between low thin and capacity.
Omega.post <- computeAssociations(m, thin = 5) 

saveRDS(Omega.post, "Omega_correlations.rds")
