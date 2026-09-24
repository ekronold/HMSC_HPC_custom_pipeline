## Fit model from chains


.libPaths(c(.libPaths(), "/cluster/projects/nn9725k/eivind/R_packages/4.5"))

# Load libraries
library(coda)
library(Hmsc)
library(jsonify)
set.seed(666) # Remember to set seed

# Load the unfitted model object

m.unf <- readRDS("unfitted_m.preds.rds") # Load the unfitted model object

# Load the fitted chains
C0 <- from_json(readRDS("post_chain00_file.rds")[[1]])
C1 <- from_json(readRDS("post_chain01_file.rds")[[1]])
C2 <- from_json(readRDS("post_chain02_file.rds")[[1]])
C3 <- from_json(readRDS("post_chain03_file.rds")[[1]])
C4 <- from_json(readRDS("post_chain04_file.rds")[[1]])
chainlist <- list("Chain_1" = C0[[1]],
                  "Chain_2" = C1[[1]],
                  "Chain_3" = C2[[1]],
                  "Chain_4" = C3[[1]],
                   "Chain_5" = C4[[1]])

# Combine the chains into a fitted model object, make sure the samples, thin and transient match the values in the hpc run
fitted.m <- importPosteriorFromHPC(m.unf, chainlist, nSamples = 500, thin = 10, transient = 20000)

# Save the fitted object
saveRDS(fitted.m, "fitted_m.rds") # Save the fitted object
