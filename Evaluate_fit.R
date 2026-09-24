# evaluate fit and crossvalidate


.libPaths(c(.libPaths(), "/cluster/projects/nn9725k/eivind/R_packages/4.4"))

# Load libraries
library(coda)
library(Hmsc)
library(jsonify)
set.seed(666) # Remember to set seed

m <- readRDS("fitted_m.rds") # Load the fitted model

nfolds <- 2 # n-fold crossvalidation

preds <- computePredictedValues(m)

# Evaluate model fit on the whole prediction set
MF <- evaluateModelFit(hM=m, predY=preds)

# Initiate cross-validation by partitioning data into n folds
partition <- createPartition(m, nfolds = nfolds)

# Compute predictions on subsets, this is time-consuming as it fits the model again on each fold
# Consider the trade-off between nfolds and time. 2 is likely fine for large models

preds <- pcomputePredictedValues(m, partition=partition, thin = 5, start = 200)

# Evaluate the cross-validated model fit
MFCV <- evaluateModelFit(hM=m, predY=preds)

# Also compute the value of WAIC (Widely Applicable Information Criterion) 

WAIC <- computeWAIC(m)

print(WAIC)

# Visualise the model fit
# Compare the full model with the cross-validated models
# TjurR2 is an adjusted R2 measurement for the 0/1 model
pdf("TjurR2_evaluate_fit.pdf")
plot(MF$TjurR2,MFCV$TjurR2, 
     xlim=c(-1,1),ylim=c(-1,1),
     xlab = "explanatory power",
     ylab = "predictive power",
     main=paste0("Tjur R2.\n",
                 "mean(MF) = ",as.character(mean(MF$TjurR2,na.rm=TRUE)),
                 ", mean(MFCV) = ",as.character(mean(MFCV$TjurR2,na.rm=TRUE))))
abline(0,1)
abline(v=0)
abline(h=0)
dev.off()

# R2 is the explanatory power of the abundance model
pdf("R2_evaluate_fit.pdf")
plot(MF$R2,MFCV$R2,
     xlim=c(-1,1),ylim=c(-1,1),
     xlab = "explanatory power",
     ylab = "predictive power",
     main=paste0("R2.\n",
                 "mean(MF) = ",as.character(mean(MF$R2,na.rm=TRUE)),
                 ", mean(MFCV) = ",as.character(mean(MFCV$R2,na.rm=TRUE))))
abline(0,1)
abline(v=0)
abline(h=0)
dev.off()


# AUC, area-under-curve is another measure of model fit
pdf("AUC_evaluate_fit.pdf")
plot(MF$AUC,MFCV$AUC,xlim=c(0,1),ylim=c(0,1),
     xlab = "explanatory power",
     ylab = "predictive power",
     main=paste0("AUC. \n",
                 "mean(MF) = ",as.character(mean(MF$AUC,na.rm=TRUE)),
                 ", mean(MFCV) = ",as.character(mean(MFCV$AUC,na.rm=TRUE))))
abline(0,1)
abline(v=0.5)
abline(h=0.5)

dev.off()
