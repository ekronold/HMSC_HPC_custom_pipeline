# This script deviates some from the normal HMSC package and makes use of a custom function to extract predicted values.
# By following the script as is you will generate predicted values across your whole model for every OTU and every trait
# for each predictor and output them all as tidy dataframes that you can filter based on OTU names, predictor names and trait names


# Set library path and load libraries

.libPaths(c(.libPaths(), "/cluster/projects/nn9725k/eivind/R_packages/4.5"))

# Load libraries
library(coda)
library(Hmsc)
library(jsonify)
library(tidyverse)
set.seed(666) # Remember to set seed

# Variation partitioning

m <- readRDS("fitted_m.rds") # Load the fitted model


# Define the custom function:

# This function is modified from plotGradient() but instead of plotting it stores the values
gradientData <- function(hM,
                         Gradient,
                         predY,
                         measure = c("Y", "T", "S"),
                         index = 1,
                         q = c(0.025, 0.5, 0.975)) {

  measure <- match.arg(measure)

  ## x values
  switch(class(hM$X)[1L],
         matrix = {
           xx <- Gradient$XDataNew[,1]
         },
         list = {
           if (measure == "Y") {
             xx <- Gradient$XDataNew[[index]][,1]
           } else {
             xx <- Gradient$XDataNew[[1]][,1]
           }
         })

  ## Species responses
  if (measure == "Y") {

    tmp <- abind::abind(predY, along = 3)

    qpred <- apply(tmp, c(1,2), quantile,
                   probs = q, na.rm = TRUE)

    qpred <- qpred[,,index]

    out <- data.frame(
      response   = hM$spNames[index],
      response_id = index,
      response_type = "species",
      x = xx,
      lower = qpred[1,],
      median = qpred[2,],
      upper = qpred[3,]
    )

  }

  ## Trait responses
  if (measure == "T") {

    if (all(hM$distr[,1] == 1)) {
      predT <- lapply(predY, function(a)
        (exp(a) %*% hM$Tr) /
          matrix(rep(rowSums(exp(a)), hM$nt), ncol = hM$nt)
      )
    } else {
      predT <- lapply(predY, function(a)
        (a %*% hM$Tr) /
          matrix(rep(rowSums(a), hM$nt), ncol = hM$nt)
      )
    }

    predT <- abind::abind(predT, along = 3)

    qpred <- apply(predT, c(1,2), quantile,
                   probs = q, na.rm = TRUE)

    qpred <- qpred[,,index]

    out <- data.frame(
      response = hM$trNames[index],
      response_id = index,
      response_type = "trait",
      x = xx,
      lower = qpred[1,],
      median = qpred[2,],
      upper = qpred[3,]
    )

  }

  ## Summed response
  if (measure == "S") {

    predS <- abind::abind(lapply(predY, rowSums), along = 2)

    qpred <- apply(predS, 1, quantile,
                   probs = q, na.rm = TRUE)

    out <- data.frame(
      response = "Summed response",
      response_id = 1,
      response_type = "sum",
      x = xx,
      lower = qpred[1,],
      median = qpred[2,],
      upper = qpred[3,]
    )

  }

  out$model_type <- switch(
    measure,
    Y = if(hM$distr[1,1] == 2)
      "presence_absence"
    else
      "abundance",

    T = "community_trait_prediction",

    S = if(hM$distr[1,1] == 2)
      "species_richness_from_occurrence"
    else
      "summed_abundance"
  )

  out
}


# Assume model object is named 'm'

### 1 - Set up a vector of predictors
# Environmental predictors in the model
envs <- colnames(m$XData)
# Alternatively if the Xformula is fit on a subset of predictors:
#envs <- colnames(m$XFormula)



### 2 - Create gradients for each predictor
### NOTE: If there is no spatial random level structure, remove the 'coordinates' argument

Gradient.list <- list()
Gradient.list.marginal <- list()
for(g in 1:length(envs)){
  covariate <- envs[g]
  Gradient.list[[g]] <- constructGradient(m, focalVariable = covariate,
                                          coordinates = list(m$ranLevels$site))

  # Look into the documentation for various effect of changing the value for non.focalVariables
  # ?constructGradient
  Gradient.list.marginal[[g]] <- constructGradient(m,focalVariable = covariate,
                                                   non.focalVariables = 1, coordinates = list(m$ranLevels$site))
}
names(Gradient.list) <- envs
names(Gradient.list.marginal) <- envs

### 3- Make predictions along the gradient for each parameter

Pred.list <- list()
Pred.marginal.list <- list()
for(p in 1:length(Gradient.list)){
  Pred.list[[p]] = predict(m, Gradient=Gradient.list[[p]], expected = TRUE)
  Pred.marginal.list[[p]] = predict(m, Gradient=Gradient.list.marginal[[p]], expected = TRUE)
}
names(Pred.list) <- names(Gradient.list)
names(Pred.marginal.list) <- names(Gradient.list.marginal)


### 4 - Apply the custom prediction function across the variables

# Species response
Species.pred.list <- list()
Species.pred.marginal.list <- list()

for(s in 1:length(Gradient.list)){
  Species.pred.list[[s]] <- dplyr::bind_rows(
    lapply(seq_len(m$ns), function(i)
      gradientData(m, Gradient.list[[s]],
                   Pred.list[[s]],
                   measure = "Y",
                   index = i))
  )
  print(paste0(names(Gradient.list[s]), ": done"))
}
names(Species.pred.list) <- names(Gradient.list)
for(s in 1:length(Gradient.list.marginal)){
  Species.pred.marginal.list[[s]] <- dplyr::bind_rows(
    lapply(seq_len(m$ns), function(i)
      gradientData(m, Gradient.list.marginal[[s]],
                   Pred.marginal.list[[s]],
                   measure = "Y",
                   index = i))
  )
  print(paste0(names(Gradient.list.marginal[s]), ": done"))
}

# Trait response
Trait.pred.list <- list()
Trait.marginal.pred.list <- list()

for(t in 1:length(Gradient.list)){
  Trait.pred.list[[t]] <- dplyr::bind_rows(
    lapply(seq_len(m$nt), function(i)
      gradientData(m, Gradient.list[[t]],
                   Pred.list[[t]],
                   measure = "T",
                   index = i))
  )
  print(paste0(names(Gradient.list[t]), ": done"))
}
names(Trait.pred.list) <- names(Gradient.list)


for(t in 1:length(Gradient.list.marginal)){
  Trait.marginal.pred.list[[t]] <- dplyr::bind_rows(
    lapply(seq_len(m$nt), function(i)
      gradientData(m, Gradient.list.marginal[[t]],
                   Pred.marginal.list[[t]],
                   measure = "T",
                   index = i))
  )
  print(paste0(names(Gradient.list.marginal[t]), ": done"))
}
names(Trait.marginal.pred.list) <- names(Gradient.list.marginal)


# Summed response (richness model)
Summed.pred.list <- list()
Summed.marginal.pred.list <- list()

for(r in 1:length(Gradient.list)){
  Summed.pred.list[[r]] <- gradientData(
    hM = m,
    Gradient = Gradient.list[[r]],
    predY = predY,
    measure = "S")
  print(paste0(names(Gradient.list[r]), ": done"))
}
names(Summed.pred.list) <- names(Gradient.list)

for(r in 1:length(Gradient.list.marginal)){
  Summed.marginal.pred.list[[r]] <- gradientData(
    hM = m,
    Gradient = Gradient.list.marginal[[r]],
    predY = predY,
    measure = "S")
  print(paste0(names(Gradient.list.marginal[r]), ": done"))
}
names(Summed.marginal.pred.list) <- names(Gradient.list.marginal)

## Save all the list objects

saveRDS(Species.pred.list, "Species_pred_list.rds")
saveRDS(Species.pred.marginal.list, "Species_pred_marginal_list.rds")

saveRDS(Trait.pred.list, "Trait_pred_list.rds")
saveRDS(Trait.marginal.pred.list, "Trait_marginal_pred_list.rds")

saveRDS(Summed.pred.list, "Summed_pred_list.rds")
saveRDS(Summed.marginal.pred.list, "Summed_marginal_pred_list.rds")


### Done. Each list now contain the predicted response of each:
### Species, Trait and the community richness/ summed abundances
### As a single data.frame formatted for ggplot for each predictor


# # Example, assuming 'CN_ratio' and 'Stand_age' are predictors:
# Trait_CN_response <- Trait.pred.list[["CN_ratio"]]
# 
# ggplot(Trait_CN_response) +
#   geom_line(aes(x, median, color = response)) +
#   geom_line(aes(x, upper, color = response), linetype = 2) +
#   geom_line(aes(x, lower, color = response), linetype = 2) +
#   xlab("CN_ratio") +
#   theme_classic()
# 
# 
# Trait_Stand_age <- Trait.pred.list[["Stand_age"]]
# 
# ggplot(Trait_Stand_age) +
#   geom_line(aes(x, median, color = response)) +
#   geom_line(aes(x, upper, color = response), linetype = 2) +
#   geom_line(aes(x, lower, color = response), linetype = 2) +
#   xlab("CN_ratio") +
#   theme_classic()
