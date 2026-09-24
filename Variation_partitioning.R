# Set library path and load libraries

.libPaths(c(.libPaths(), "/cluster/projects/nn9725k/eivind/R_packages/4.5"))

# Load libraries
library(coda)
library(Hmsc)
library(jsonify)
set.seed(666) # Remember to set seed

# Variation partitioning

m <- readRDS("fitted_m.rds") # Load the fitted model

# Create the objects



preds = computePredictedValues(m) # predicted values
VP = computeVariancePartitioning(m) # partition variance
saveRDS(VP, "VarPart.rds")
# These two functions has many options for different grouping of variables and such.
# Check documentation for more settings

vals = VP$vals # Extract the values for variance partitioned based on parameter and species

MF = evaluateModelFit(hM=m, predY=preds)
saveRDS(MF, "ModelFit.rds")
R2 = MF$R2 # R2 on the abundance model
TjurR2 = MF$TjurR2 # R2 on the presence/absence model

# Combine the values
vals = rbind(vals,R2, TjurR2)

# Separate the two models to save one dataframe for each 
# This is not neccessary when running a single model with the same distribution on all data.
# Adapt as needed
#vals.pa <- vals %>% data.frame() %>% 
#  select(contains(".pa")) %>% 
#  rownames_to_column("Parameter") %>% 
#  filter(!Parameter == "R2")

#vals.abund <- vals %>% data.frame() %>% 
#  select(contains(".abund")) %>% 
#  rownames_to_column("Parameter") %>% 
#  filter(!Parameter == "TjurR2") 


# These dataframes will now contain the proportion of explained variation for each OTU and each predictor
# It is recommended to arrange by descending R2, so that you can sort the responses by effect size
# The Beta summary will contain posterior probability values that can also be used to filter OTUs
# OTUs inherit names from the model

#
write.csv(vals, "VarPart_model.csv", quote = F, row.names = F)
#write.csv(vals.pa, "VarPart_pa_model.csv", quote = F, row.names = F)
#write.csv(vals.abund, "VarPart_abund_model.csv", quote = F, row.names = F)
# 

