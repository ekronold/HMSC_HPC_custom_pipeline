# Script for generating processed data tables from a fitted HMSC-hpc model
# Make sure the libraries are loaded from the right directory

.libPaths(c(.libPaths(), "/cluster/projects/nn9725k/eivind/R_packages/4.5"))

# Load libraries
library(coda)
library(Hmsc)
library(jsonify)
set.seed(666) # Remember to set seed 


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

# Load the fitted model
m <- readRDS("fitted_m.rds")

mpost <- convertToCodaObject(m, Omega = FALSE) # This may take a long time
# Omega is the species-to-species correlations and is computationally very heavy
# Estimated on its own later in the script

saveRDS(mpost, "mpost.rds")

#mpost <- readRDS("mpost.rds")

summary.beta <- summary(mpost$Beta)

saveRDS(summary.beta, "Summary_beta.rds")

summary.gamma <- summary(mpost$Gamma)

saveRDS(summary.gamma, "Summary_gamma.rds")

eff.sizeB <- effectiveSize(mpost$Beta)

summary.lambda <- summary(mpost$Lambda)

saveRDS(summary.lambda, "Summary_lambda.rds")

eff.sizeL <- effectiveSize(mpost$Lambda[[1]])

summary.Eta <- summary(mpost$Eta)

saveRDS(summary.Eta, "Summary_eta.rds")

eff.sizeE <- effectiveSize(mpost$Eta[[1]])

# Gelman diagnostics of parameters become too big for R with many speciesXparameters
# Subsample randomly 1000 and investigate
idx <- sample(1:nvar(mpost$Beta), 1000)

gel.diagB <- gelman.diag(
  mpost$Beta[, idx],
  multivariate = FALSE
)$psrf

write.csv(eff.sizeB, "Eff_size_beta.csv")
write.csv(gel.diagB, "Gelman_diag_beta.csv")

pdf("Effective_size_beta.pdf")
hist(eff.sizeB)
dev.off()

pdf("Gelman_diag_beta.pdf")
hist(gel.diagB[,1])
dev.off()

pdf("Chain_mixing_beta.pdf")
plot(mpost$Beta)
dev.off()

eff.sizeG <- effectiveSize(mpost$Gamma)
gel.diagG <- gelman.diag(mpost$Gamma, multivariate= FALSE)$psrf

write.csv(eff.sizeG, "Eff_size_gamma.csv")
write.csv(gel.diagG, "Gelman_diag_gamma.csv")

pdf("Effective_size_gamma.pdf")
hist(eff.sizeG)
dev.off()

pdf("Gelman_diag_gamma.pdf")
hist(gel.diagG)
dev.off()

pdf("Chain_mixing_gamma.pdf")
plot(mpost$Gamma)
dev.off()

gel.diagL <- gelman.diag(mpost$Lambda[[1]], multivariate=FALSE)$psrf

write.csv(eff.sizeL, "Eff_size_lambda.csv")
write.csv(gel.diagL, "Gelman_diag_lambda.csv")

pdf("Effective_size_labmda.pdf")
hist(eff.sizeL)
dev.off()

pdf("Gelman_diag_lambda.pdf")
hist(gel.diagL)
dev.off()

gel.diagE <- gelman.diag(mpost$Eta[[1]], multivariate=FALSE)$psrf

write.csv(eff.sizeE, "Eff_size_eta.csv")
write.csv(gel.diagE, "Gelman_diag_eta.csv")

pdf("Effective_size_eta.pdf")
hist(eff.sizeE)
dev.off()

pdf("Gelman_diag_eta.pdf")
hist(gel.diagE)
dev.off()

# Omega diagnostics are very heavy. Compute on a random small subset
mpost.Omega <- convertToCodaObject(m, Beta = FALSE,
  Gamma = FALSE,
  V = FALSE,
  Sigma = FALSE,
  Rho = FALSE,
  Eta = FALSE,
  Lambda = FALSE,
  Alpha = FALSE,
  Omega = TRUE,
  Psi = FALSE,
  Delta = FALSE,
  start = 200) # Set the starting point later along the sampling chain. Can be adjusted. Starting later in the chain speeds up computation but may reduce accuracy
saveRDS(mpost.Omega, "mpostOmega.rds")

tmp <- mpost.Omega$Omega
maxOmega = 1000 # Subset of 1000 OTUs can also be adjusted as needed
ma.omega <- list()
for (i in 1:length(tmp[[1]])){
  z <- dim(tmp[[1]][[i]])[2]
  sel = sample(1:z, size = maxOmega)
  
  for(j in 1:length(tmp)){
  tmp[[j]] = tmp[[j]][,sel]
  }
  psrf <- gelman.diag(tmp[[1]], multivariate = FALSE)$psrf
  ma.omega[[i]] <- psrf[,1]
}

saveRDS(ma.omega, "psrf_list_Omega.rds")

pdf("Omega_psrf.pdf")
for (f in 1:length(ma.omega)){
  hist(ma.omega[[f]], main = paste0("Omega psrf, random level ", f))
}
dev.off()


