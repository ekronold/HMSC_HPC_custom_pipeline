# Simple spatial model fitting HMSC 
setwd("example_model/")
Y <- read.csv("YData.csv") %>% 
  column_to_rownames("X")
XData <- read.csv("XData.csv") %>% 
  column_to_rownames("X")
TrData <- read.csv("TrData.csv") %>% 
  column_to_rownames("X") %>% 
  as.matrix()
XY <- read.csv("XY_cords.csv") %>% 
  column_to_rownames("X.1")

# Set up the spatial component
coords <- as.matrix(XY)

# convex hull
hull_ind <- chull(coords)
hull_coords <- rbind(coords[hull_ind, ],
                     coords[hull_ind[1], ])

# spatial polygon object
hull_sp <- SpatialPolygons(list(
  Polygons(list(Polygon(hull_coords)), ID=1)
))

# choose reasonable number of knots
nKnots <- 8  # This represent equally spaced points within the hull for the gaussian predictive model

kn <- as.data.frame(
  spsample(hull_sp, nKnots, type="hexagonal")
)

# Check spatial points and hull
plot(XY)
points(hull_coords, col = "blue", pch = 19)
lines(hull_coords, col = "orange")
points(kn, col = "red", pch = 17)

studyDesign <- data.frame(site=as.factor(rownames(XData)))

rLSite <- HmscRandomLevel(
  sData = XY,
  sMethod = "GPP",
  sKnot = kn
)

nf = 2 # number of latent factors 
rLSite <- setPriors(rLSite, nfMin = nf, nfMax = nf)

# Check the data
Y[1:10, 1:10] # P/A
Y[1:10, 120:130] # scaled abundance given presence


names(XData) 
head(XData)
names(TrData)
head(TrData)
XFormula = ~ CN_ratio + Stand_age + Elevation
TrFormula = ~ .


#Create the model
m <- Hmsc(
  Y = Y,
  XData = XData, XFormula = XFormula,
  Tr = TrData, TrFormula = TrFormula,
  # Set the distribution to probit for the 0/1 data and normal for the abundance part
  distr = c(rep("probit", ncol(Y)/2), rep("normal", ncol(Y)/2)),
  studyDesign = studyDesign,
  ranLevels = list(site = rLSite)
)
# Check that model runs
sampleMcmc(m, samples = 2)
# 

# Run the small local model to start playing with downstream functions 
# m_local <- sampleMcmc(
#   m,
#   samples   = 100,
#   thin      = 5,
#   transient = 5000,
#   nChains   = 2,
#   engine    = "R"
# )


# Save a init object to run on HPC
# init.obj.test <- sampleMcmc(
#   m,
#   samples = ,
#   thin = ,
#   transient = ,
#   nChains = ,
#   engine = "HPC"
# )
