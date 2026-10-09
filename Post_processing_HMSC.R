### HMSC Model post-processing

library(tidyverse)
library(Hmsc)
library(vegan)
library(corrplot)
library(openxlsx)
library(reshape)

#### Load the model object
m <- readRDS("unfitted_m.preds.rds")
m.fitted <- readRDS("fitted_m.maximal.rds")

#### Load the post estimates
mPost <- readRDS("results/Post_estimates.rds")

## Extract the input data
Y <- data.frame(m$Y)
X <- data.frame(m$XData)
Tr <- m$TrData
Tax <- read.csv("data/Taxonomy.csv") %>% select(!"X")
Cords <- data.frame(m$ranLevels$site$s) %>% 
  rownames_to_column("Sites")
SpatialKnots <- m$ranLevels$site$sKnot

# Spatial random factor:
ggplot() +
  geom_point(data = Cords, aes(X, Y), fill = "purple", size = 2, shape = 21) +
  geom_point(data = SpatialKnots, aes(x, y), fill = "gold", size = 3, shape = 21) +
  theme_classic()

#### Evaluate model fit ###### 

MF <- readRDS("results/Evaluate_model_fit.rds")
# RMSE
summary(MF$RMSE)
hist(MF$RMSE)
mean(MF$RMSE)

# Check the two models, probit and abundance-given-presence
# probit:
summary(MF$RMSE[1:1835])
hist(MF$RMSE[1:1835])

# abund:
summary(MF$RMSE[1836:3670])
hist(MF$RMSE[1836:3670])

# AUC and TjurR2, only valid for probit model. In essence, this is the accuracy in predicting a presence
summary(na.omit(MF$AUC))
hist(na.omit(MF$AUC))

summary(na.omit(MF$TjurR2))
hist(na.omit(MF$TjurR2))

# How well does the models ability to predict if an OTU is present correlate with how well it can predict the abundance?
plot(na.omit(MF$AUC), na.omit(MF$TjurR2)) 
cor(na.omit(MF$AUC), na.omit(MF$TjurR2))

# R2, explanatory power in predicting abundances
summary(na.omit(MF$R2))
hist(na.omit(MF$R2))

# Save the model fit estimates for each OTU 
Fit_per_OTU <- data.frame("OTU" = Tax$OTUnr,
                          "AUC" = na.omit(MF$AUC),
                          "R2" = na.omit(MF$R2),
                          "TjurR2" = na.omit(MF$TjurR2))

# Explanation of the terms:
#
# AUC:
# Measures discrimination of presence vs absence.
# High AUC (>0.80) means the model is good at ranking sites, 
# assigning higher predicted occurrence probabilities to sites where 
# the OTU is observed as present than to sites where it is absent.
#
# TjurR2:
# Measures the separation between predicted occurrence probabilities 
# for observed presences and absences.
# Higher values indicate that the model predicts clearly different 
# probabilities for occupied vs unoccupied sites.
# TjurR2 can be low to moderate even when AUC is high because AUC 
# measures ranking, whereas TjurR2 depends on the magnitude of the 
# difference between predicted probabilities.
#
# R2:
# Measures the explanatory power of the abundance-given-presence model.
# Among sites where the OTU is present, it describes how well the model 
# explains variation in abundance based on the environmental predictors.
# High R2 indicates that abundance variation is strongly related to the 
# predictors included in the model.

# Look at the PSRF (potential scale reduction factor) and Effective Sample Size stats

EffSizeBeta <- read.csv("results/convergence/Eff_size_beta.csv", sep = ",")
EffSizeGamma <- read.csv("results/convergence/Eff_size_gamma.csv")
EffSizeLambda <- read.csv("results/convergence/Eff_size_lambda.csv")
EffSizeEta <- read.csv("results/convergence/Eff_size_eta.csv")
hist(EffSizeBeta$x)
head(EffSizeBeta, n = 22)
hist(EffSizeGamma$x)
hist(EffSizeLambda$x)
hist(EffSizeEta$x)
# Generate PSRF stats in the HPC scripts

# Parameter summaries contain detailed information on the estimates for each parameter. Can be investigated if necessary 
# They are structured in a slightly messy way so getting proper columns need a bit of cleaning
SummaryBeta <- readRDS("results/summary/Summary_beta.rds")
SummaryGamma <- readRDS("results/summary/Summary_gamma.rds")
head(SummaryBeta$quantiles, n = 22)
head(SummaryBeta$statistics, n = 22)

BetaQuantiles <- as.data.frame(SummaryBeta$quantiles) %>% 
  # Set the names which are in a very difficult notation
  rownames_to_column("HMSC_name")

# Set the colnames more neatly
colnames(BetaQuantiles)[2:6] <- c("Q2_5","Q25","Q50","Q75","Q97_5")

# Clean the names for the Parameters
BetaQuantiles$Parameter <- sub(
  "B\\[(.*) \\(C[0-9]+\\),.*",
  "\\1",
  BetaQuantiles$HMSC_name
)
# Clean the names for the OTUs
BetaQuantiles$OTU <- sub(
  ".*\\), (.*) \\(S[0-9]+\\)\\]",
  "\\1",
  BetaQuantiles$HMSC_name
)

# Remove the comples HMSC name
BetaQuantiles <- BetaQuantiles %>% 
  select(!HMSC_name)

# Do the same with the summary statistic
BetaStatistic <- as.data.frame(SummaryBeta$statistics) %>% 
  rownames_to_column("HMSC_name")

colnames(BetaStatistic)[2:5] <- c("Mean", "SD","Naive_SE","Time_series_SE")
BetaStatistic$Parameter <- sub(
  "B\\[(.*) \\(C[0-9]+\\),.*",
  "\\1",
  BetaStatistic$HMSC_name
)
BetaStatistic$OTU <- sub(
  ".*\\), (.*) \\(S[0-9]+\\)\\]",
  "\\1",
  BetaStatistic$HMSC_name
)

BetaStatistic <- BetaStatistic %>% select("OTU",
                                          "Parameter",
                                          "Mean",
                                          "SD",
                                          "Naive_SE",
                                          "Time_series_SE")

# Merge the quantiles and statistic data
BetaSummaryClean <- merge(
  BetaStatistic,
  BetaQuantiles,
  by = c("OTU", "Parameter")
)

# Filter responses with a positive estimate, i.e. the lower quartile is still above 0
BetaPositive <- subset(
  BetaSummaryClean,
  Q2_5 > 0
)
# Filter negative responses in the same way
BetaNegative <- subset(
  BetaSummaryClean,
  Q97_5 < 0
)

# Separate the .abund and .pa models

BetaPositive.pa <- BetaPositive %>% 
  separate(OTU, c("OTU", "Model"), sep = "\\.") %>% 
  filter(Model == "pa")
BetaPositive.abund <- BetaPositive %>% 
  separate(OTU, c("OTU", "Model"), sep = "\\.") %>% 
  filter(Model == "abund")

BetaNegative.pa <- BetaNegative %>% 
  separate(OTU, c("OTU", "Model"), sep = "\\.") %>% 
  filter(Model == "pa")
BetaNegative.abund <- BetaNegative %>% 
  separate(OTU, c("OTU", "Model"), sep = "\\.") %>% 
  filter(Model == "abund")

# Combine the pa and abund responses

Beta.pa <- rbind(BetaPositive.pa, BetaNegative.pa) %>% 
  filter(Parameter != "(Intercept)")
Beta.abund <- rbind(BetaPositive.abund, BetaNegative.abund) %>%
  filter(Parameter != "(Intercept)")

# Plot the responses
ggplot(Beta.pa, aes(x = OTU, y = Mean)) +
  geom_point() +
  coord_flip()+ 
  geom_errorbar(
    aes(
      ymin = Q2_5,
      ymax = Q97_5
    ),
    width = 0.2
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  facet_wrap(~ Parameter, scales = "free_y") +
  theme_classic()

ggplot(Beta.abund, aes(x = OTU, y = Mean)) +
  geom_point() +
  coord_flip()+ 
  geom_errorbar(
    aes(
      ymin = Q2_5,
      ymax = Q97_5
    ),
    width = 0.2
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  facet_wrap(~ Parameter, scales = "free_y") +
  theme_classic()
  

#### Investigate Beta parameters, OTU responses to predictor variables
PostBeta <- mPost$postBeta # Contains the parameter estimates as well as the posterior probabilities

#plotBeta(m, mPost$postBeta)

PostBetaEst <- data.frame(PostBeta$mean) # contains the mean estimate for each OTU both for the abundance and 0/1 model

PostBetaEst[1:nrow(PostBetaEst),1:5]  # The basic output does not contain the predictor names by default
# Set predictor names
m$XFormula
colnames(X)
rownames(PostBetaEst) <- c("Intercept", colnames(X)) # First column contains the intercept values, the rest are predictors in the order they are formulated in the XFormula. Should be in the same order as they appear in the XData, print both to control

PostBetaEst.pa <- PostBetaEst %>% 
  # Extract the 0/1 part of the model
  select(contains(".pa")) %>%  
  # Transpose for ggplotting
  t() %>% data.frame() %>% 
  # The intercept is not really informative, so remove it
  select(!Intercept)

PostBetaEst.abund <- PostBetaEst %>% 
  # Extract the abundance model
  select(contains(".abund")) %>% 
  # Transpose for ggplotting
  t() %>% data.frame() %>% 
  # The intercept is not really informative, so remove it
  select(!Intercept)

# Control the estimates by posterior probability
# The basic output contains two posterior probability matrixes
# One is the probability of a positive response, one is the probability of a negative response
# These are inverse of each other, i.e. the support + supportNeg for each OTU sums to 1
sum(PostBeta$support[100], PostBeta$supportNeg[100])

PostBetaSup <- data.frame(PostBeta$support)
# This is in the same format as before, so set names, filter and transpose
rownames(PostBetaSup) <- c("Intercept", colnames(X))

PostBetaSup.pa <- PostBetaSup %>% 
  select(contains(".pa")) %>%  
  t() %>% data.frame() %>% 
  select(!Intercept)

PostBetaSup.abund <- PostBetaSup %>% 
  select(contains(".abund")) %>%  
  t() %>% data.frame() %>% 
  select(!Intercept)

# Following a p = 0.05 logic, apply a two-sided 95% probability to the estimates
PostBetaSup.pa[PostBetaSup.pa > 0.025 & PostBetaSup.pa < 0.975] = NA 
PostBetaSup.abund[PostBetaSup.abund > 0.025 & PostBetaSup.abund < 0.975] = NA

# Set the estimates without strong evidence of the parameter estimates to 0
PostBetaEst.pa.sig <- PostBetaEst.pa
PostBetaEst.abund.sig <- PostBetaEst.abund

# Set the estimates below the probability threshold to NA
PostBetaEst.pa.sig[is.na(PostBetaSup.pa)] = 0
PostBetaEst.abund.sig[is.na(PostBetaSup.abund)] = 0

# Set to long format for ggplot and add the taxonomy
# All estimates
PostBetaEst.pa.long <- PostBetaEst.pa %>%
  rownames_to_column("OTUnr") %>% 
  melt() %>% 
  set_names(c("OTUnr","Parameter","Response")) %>% 
  mutate("OTUnr" = str_replace(.$OTUnr, ".pa", "")) %>% 
  left_join(Tax, by = "OTUnr")

PostBetaEst.abund.long <- PostBetaEst.abund %>%
  rownames_to_column("OTUnr") %>% 
  melt() %>% 
  set_names(c("OTUnr","Parameter","Response")) %>% 
  mutate("OTUnr" = str_replace(.$OTUnr, ".abund", "")) %>% 
  left_join(Tax, by = "OTUnr")

# Only estimates wihtin the two-sided 95% probability threshold
PostBetaEst.pa.sig.long <- PostBetaEst.pa.sig %>%
  rownames_to_column("OTUnr") %>% 
  melt() %>% 
  set_names(c("OTUnr","Parameter","Response")) %>% 
  filter(Response != 0) %>% 
  mutate("OTUnr" = str_replace(.$OTUnr, ".pa", "")) %>% 
  left_join(Tax, by = "OTUnr")

PostBetaEst.abund.sig.long <- PostBetaEst.abund.sig %>%
  rownames_to_column("OTUnr") %>% 
  melt() %>% 
  set_names(c("OTUnr","Parameter","Response")) %>% 
  filter(Response != 0) %>% 
  mutate("OTUnr" = str_replace(.$OTUnr, ".abund", "")) %>% 
  left_join(Tax, by = "OTUnr")

# Plot the parameter estimates in ggplot:
ggplot(PostBetaEst.pa.long, aes(x = Parameter, y = OTUnr, fill = Response)) +
  labs(x = "Environmental Niche", y = "OTU", fill = "Response") +
  geom_raster()+
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 5))


ggplot(PostBetaEst.abund.long, aes(x = Parameter, y = OTUnr, fill = Response)) +
  labs(x = "Environmental Niche", y = "OTU", fill = "Response") +
  geom_raster()+
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 5))


ggplot(PostBetaEst.pa.sig.long, aes(x = Parameter, y = OTUnr, fill = Response)) +
  labs(x = "Environmental Niche", y = "OTU", fill = "Response") +
  geom_raster()+
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 5))


ggplot(PostBetaEst.abund.sig.long, aes(x = Parameter, y = OTUnr, fill = Response)) +
  labs(x = "Environmental Niche", y = "OTU", fill = "Response") +
  geom_raster()+
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 5))

#### Investigate the Gamma parameters, trait responses to environment
# In the hurdle model, these are estimated both as a richness response (0/1) and a change in abundances in the same output

# Check the coverage of traits, i.e. how many OTUs the traits consist of
summary(Tr)

PostGamma <- mPost$postGamma

PostGammaEst <- PostGamma$mean
rownames(PostGammaEst) <- c("Intercept", colnames(X))
# These by default does not have the trait identities either, fill in
colnames(PostGammaEst) <- c("Intercept", colnames(Tr))

PostGammaEst <- PostGammaEst %>% 
  t() %>% data.frame() %>% 
  select(!Intercept) %>% 
  rownames_to_column("Guild") %>% 
  filter(Guild != "Intercept") 
         
PostGammaSup <- data.frame(PostGamma$support)
# This is in the same format as before, so set names, filter and transpose
rownames(PostGammaSup) <- c("Intercept", colnames(X))
colnames(PostGammaSup) <- c("Intercept", colnames(Tr))

PostGammaSup <- PostGammaSup %>% 
  t() %>% data.frame() %>% 
  select(!Intercept) %>% 
  rownames_to_column("Guild") %>% 
  filter(Guild != "Intercept")

# Following a p = 0.05 logic, apply a two-sided 95% probability to the estimates
PostGammaSup[PostGammaSup > 0.025 & PostGammaSup < 0.975] = NA 

PostGammaEst.sig <- PostGammaEst
PostGammaEst.sig[is.na(PostGammaSup)] = 0

PostGammaEst.long <- PostGammaEst %>%
  melt() %>% 
  set_names(c("Guild","Parameter","Response"))

PostGammaEst.sig.long <- PostGammaEst.sig %>%
  melt() %>% 
  set_names(c("Guild","Parameter","Response"))

ggplot(PostGammaEst.long, aes(x = Parameter, y = Guild, fill = Response)) +
  labs(x = "Environmental Niche", y = "Guild", fill = "Response") +
  geom_raster()+
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(PostGammaEst.sig.long, aes(x = Parameter, y = Guild, fill = Response)) +
  labs(x = "Environmental Niche", y = "Guild", fill = "Response") +
  geom_raster()+
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#### Investigate the Eta parameters, residual ordination of the sites 
# This is the potentially structured residual effects from the spatial model and the posterior probabilities are not very meaningful
# Instead, investigate the residual ordination with the spatial model in mind:

PostEtaScores <- data.frame(mPost$postEta$mean)
# Four axis as the model was specified with four latent factors
# This again by default does not come with the site names, so find them in the model
Sites <- levels(m$rL$site$pi)

PostEtaScores$Sites <- Sites

# Add the coordinates
PostEta <- PostEtaScores %>% 
  left_join(Cords, by = "Sites")

# Make a plot and see how well the spatial structure might be captured in the residuals

ggplot(PostEta) +
  geom_point(aes(X, Y, color = X1), size = 5) +
  geom_point(data = SpatialKnots, aes(x, y), color = "red") +
  scale_color_viridis_c() +
  theme_classic()

ggplot(PostEta, aes(X, Y, color = X2)) +
  geom_point(size = 5) +
  scale_color_viridis_c() +
  theme_classic()

ggplot(PostEta, aes(X, Y, color = X3)) +
  geom_point(size = 5) +
  scale_color_viridis_c() +
  theme_classic()

ggplot(PostEta, aes(X, Y, color = X4)) +
  geom_point(size = 5) +
  scale_color_viridis_c() +
  theme_classic()

# A structured spatial effect would show up on the map as clusters or gradients

# Can also be plotted as a typical ordination
# Add the environmental predictors so color the sites by categories
Env <- X %>% rownames_to_column("Sites")

Categories <- Env %>% select(Sites)
Categories$Type <- "Forest"
Categories$Type[Env$High_impact_non_forest == 1] = "High_impact_non_forest"
Categories$Type[Env$Wetland == 1] = "Wetland"
Categories$Management <- "Managed"
Categories$Management[Env$selectively_logged == 1] = "Selectively_logged"
Categories$Management[Env$near_natural == 1] = "Near_natural"

PostEtaEnv <- PostEta %>% left_join(Categories, by = "Sites")

ggplot(PostEtaEnv) +
  geom_point(aes(X1, X2, color = Management, shape = Type), size = 4) +
  theme_classic()

ggplot(PostEtaEnv) +
  geom_point(aes(X1, X3, color = Management, shape = Type), size = 4) +
  theme_classic()

ggplot(PostEtaEnv) +
  geom_point(aes(X1, X4, color = Management, shape = Type), size = 4) +
  theme_classic()

# Alpha is an estimate of how the spatial correlation decays with distance

PostAlpha <- data.frame("Alpha" = mPost$postAlpha$mean,
                        "Support" = mPost$postAlpha$support)
PostAlpha
# Large value with strong support for the first axis, might indicate a spatial correlation across the whole system
# Smaller value for the fourth axis, might indictae an intermediate distance structure independent of the first axis
# Low support for axis 2 and 3, might not be very important for the whole system
# Interpet a bit with care


# Lambda loading can also be investigated, these represent species loadings

PostLambda <- mPost$postLambda

PostLambdaScore <- data.frame(PostLambda$mean)

# Again the OTUs need to be given names
colnames(PostLambdaScore) <- colnames(Y)

# transpose for ggplot, separate the two models and add taxonomy

PostLambdaScore.pa <- PostLambdaScore %>% 
  select(contains(".pa")) %>%
  t() %>% data.frame() %>% 
  rownames_to_column("OTUnr") %>% 
  mutate("OTUnr" = str_replace(.$OTUnr, ".pa", "")) %>% 
  left_join(Tax, by = "OTUnr")

PostLambdaScore.abund <- PostLambdaScore %>% 
  select(contains(".abund")) %>%
  t() %>% data.frame() %>% 
  rownames_to_column("OTUnr") %>% 
  mutate("OTUnr" = str_replace(.$OTUnr, ".abund", "")) %>% 
  left_join(Tax, by = "OTUnr")  

# This can be correlated to the Eta scores and look into species with potential loading following the Eta structure
# Or plotted as a residual species-ordination

# Clear patterns in Lambda may indicate responses to environmental factors that are not included in the model, or biotic associations

ggplot(PostLambdaScore.pa, aes(X1, X2)) +
  geom_point()
ggplot(PostLambdaScore.pa, aes(X1, X3)) +
  geom_point()
ggplot(PostLambdaScore.pa, aes(X1, X4)) +
  geom_point()

# Both Eta and Lambda can be compared with a null-model without predictors so investigate the change in loadings/scores by adding predictors

##### Invertigate the Omega component, potential co-occurrences
# Currently output in two objects from the pipeline, both elements has a support value

# Post estimates
PostOmega <- mPost$postOmega
# Omega correlations derived from the post estimates by the function computeAssociations() which can be set to compute
# correlations between species for each random level
# In the spatial model, there is only one random level
OmegaCor <- readRDS("results/Omega_correlations.rds")

# The correlations can be visualised simply with coppplot

OmegaCorMean <- OmegaCor[[1]]$mean
OmegaCorSupport <- OmegaCor[[1]]$support

hist(OmegaCorMean) 
hist(OmegaCorSupport)

# Apply a 0.95 probability threshold to the correlation
OmegaCorMean.sig <- OmegaCorMean
OmegaCorMean.sig[OmegaCorSupport < 0.95] = 0

# Plotting the entire correlation matrix is heavy
# corrplot(OmegaCorMean.sig, # The correlation matrix to plot
#          method = "square", # Plot filled squares
#          #type = "lower", # As the matrix is symmetric plot only the lower section
#          addgrid.col = NA, # As we are plotting very many squared, remove the small grid around each cell
#          tl.cex = 0.1, # Size of the labels on the axis
#          cl.cex = 0.5, # Size of the text on the legend
#          order = "AOE", # How to cluster the correlations
#          diag = FALSE
#          )
# 
# # It is most informative to split the correlation into three matrices:
# # .pa against .pa, .abund against .abund
# 
# OmegaCor1 <- OmegaCorMean.sig[1:1835, 1:1835]
# OmegaCor2 <- OmegaCorMean.sig[1836:3670, 1836:3670]
# 
# 
# corrplot(OmegaCor1, # The correlation matrix to plot
#          method = "square", # Plot filled squares
#          type = "lower", # As the matrix is symmetric plot only the lower section
#          addgrid.col = NA, # As we are plotting very many squared, remove the small grid around each cell
#          tl.cex = 0.1, # Size of the labels on the axis
#          cl.cex = 0.5, # Size of the text on the legend
#          order = "AOE", # How to cluster the correlations
#          diag = FALSE # Wether or not to plot the diagonal
# ) 
# 
# corrplot(OmegaCor2, # The correlation matrix to plot
#          method = "square", # Plot filled squares
#          type = "lower", # As the matrix is symmetric plot only the lower section
#          addgrid.col = NA, # As we are plotting very many squared, remove the small grid around each cell
#          tl.cex = 0.1, # Size of the labels on the axis
#          cl.cex = 0.5, # Size of the text on the legend
#          order = "AOE", # How to cluster the correlations
#          diag = FALSE # Wether or not to plot the diagonal
# ) 

# In addition to the correlations, it is possible to use the PostOmega object to extract positive/negative associations 
# and filter on species groups

PostOmegaPos <- PostOmega$support
PostOmegaNeg <- PostOmega$supportNeg

PostOmegaMean.pos <- data.frame(PostOmega$mean)
PostOmegaMean.neg <- data.frame(PostOmega$mean)

PostOmegaMean.pos[PostOmegaPos < 0.95] = 0
PostOmegaMean.neg[PostOmegaNeg < 0.95] = 0

PostOmegaMean.pos.pa <- PostOmegaMean.pos[1:1835, 1:1835]
PostOmegaMean.neg.pa <- PostOmegaMean.neg[1:1835, 1:1835]

rownames(PostOmegaMean.pos.pa) <- Tax$OTUnr
colnames(PostOmegaMean.pos.pa) <- Tax$OTUnr
rownames(PostOmegaMean.neg.pa) <- Tax$OTUnr
colnames(PostOmegaMean.neg.pa) <- Tax$OTUnr

# Are there positive co-occurrence patterns between ECM fungi and litter saptrotrophs?

OmegaEst.ECM.SAP.pos <- PostOmegaMean.pos.pa[rownames(PostOmegaMean.pos.pa) %in% Tax$OTUnr[Tax$Guild == "ectomycorrhizal"], colnames(PostOmegaMean.pos.pa) %in% Tax$OTUnr[Tax$Guild == "litter_saprotroph"]] %>%
  rownames_to_column("OTUnr") %>%
  melt() %>% 
  set_names("ECM", "LitSAP", "Estimate") %>% 
  filter(Estimate > 0)

ggplot(OmegaEst.ECM.SAP.pos, aes(x = LitSAP, y = ECM, fill = Estimate)) +
  labs(x = "Litter saprotrophs", y = "Ectomycorrhizal", fill = "Estimate") +
  geom_raster()+
  scale_fill_gradient(low = "gold",
                       high = "magenta")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Same sub-sample logic can be applied to the negative association and to the OmegaCor object

OmegaEst.ECM.SAP.neg <- PostOmegaMean.neg.pa[rownames(PostOmegaMean.neg.pa) %in% Tax$OTUnr[Tax$Guild == "ectomycorrhizal"], colnames(PostOmegaMean.neg.pa) %in% Tax$OTUnr[Tax$Guild == "litter_saprotroph"]] %>%
  rownames_to_column("OTUnr") %>%
  melt() %>% 
  set_names("ECM", "LitSAP", "Estimate") %>% 
  filter(Estimate < 0)

ggplot(OmegaEst.ECM.SAP.neg, aes(x = LitSAP, y = ECM, fill = Estimate)) +
  labs(x = "Litter saprotrophs", y = "Ectomycorrhizal", fill = "Estimate") +
  geom_raster()+
  scale_fill_gradient(low = "purple",
                      high = "white")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##### Plot the variation partitioning

VarPart <- read.csv("results/VarPart_model.csv")

# Set the rownames again
rownames(VarPart) <- c(colnames(X),"Spatial", "R2", "TjurR2")

# This table contain the VarPart for both parts of the model. Separate them before plotting

VarPart.pa <- VarPart %>% select(contains(".pa")) %>% 
  t() %>% data.frame() %>% arrange(desc(Spatial)) %>% 
  rownames_to_column("OTUnr") %>% 
  mutate("OTUnr" = str_replace(.$OTUnr, ".pa", ""))
OTU_order <- factor(VarPart.pa$OTUnr, levels = VarPart.pa$OTUnr)

# Format long while keeping the R2 separate
VarPart.pa.R2 <- VarPart.pa %>% 
  select(OTUnr, TjurR2)

VarPart.pa.long <- VarPart.pa %>% 
  select(!TjurR2) %>% 
  rownames_to_column("OTUnr") %>% 
  melt() %>% 
  set_names("OTUnr", "Parameter", "VarPart") %>% 
  left_join(VarPart.pa.R2, by = "OTUnr") %>% 
  mutate("OTUnr" = str_replace(.$OTUnr, ".pa", "")) %>% 
  left_join(Tax, by = "OTUnr") %>% 
  filter(Parameter != "R2") %>% 
  group_by(Parameter) %>% 
  arrange(desc(VarPart))
VarPart.pa.long$OTUnr <- factor(VarPart.pa.long$OTUnr, levels = OTU_order)

# Plot the full table

ggplot(VarPart.pa.long,
       aes(x = OTUnr,
           y = VarPart,
           fill = Parameter)) +
  geom_col() +
  theme_bw() +
  #coord_flip() +
  theme(
    axis.text.x = NULL)


