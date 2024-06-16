###################################################
#
# Example Code predMOB tree and predMOB forest, 
#    including calculation of variable importance 
#    measures (permutation importance and mean minimal depth)
# 
####################################################


require("partykit")
require("model4you") # only required for objfun function to calculate permutation importance 
#source("helpfunctions.R")
# NW: path adjusted
source("C:/Users/magis/OneDrive - Boehringer Ingelheim/Documents/Masterarbeit/Programming/predMOB/helpfunctions.R")
# source("C:/Dissertation/Paper/Paper 1 - SiM - Identification of predictive biomarkers/R Code/helpfunctions.R")

require(data.table) # only for rbindlist function
require(dplyr) # used to determine mean minimal depth
require(randomForestExplainer) # used to determine mean minimal depth

#########################################################

### generate data for a simple example with 3 binary biomarkers x1-x3 and a binary treatment variable T
n <- 200    # number of observations
data <- data.frame(ID=1:n,
                   x1=rbinom(n, 1,0.5),
                   x2=rbinom(n, 1,0.5),
                   x3=rbinom(n, 1,0.5),
                   T=rbinom(n, 1,0.5))

# generate normally distributed outcome Y such that x2 is prognostic and x3 is predictive
predictor <- 0.6*data$T+0.5*data$x2+0.8*data$x3*data$T-0.75
data$Y <- rnorm(n, mean=predictor, sd=0.25)

# define effect-coded treatment variable (experimental group: +0.5, control group: -0.5)
data$T.effect <- ifelse(data$T==1, 1, -1)/2   

# check simulation by using a standard linear model with interactions
summary(lm(Y~T*(x1+x2+x3),data=data))

#########################################################

### single predMOB tree
# first part of the formula defines base model, 
# potential predictive factors follow separated by a pipe symbol (cf. partykit:::mob)
# CAVE: the formula for the base model must not include an intercept
predMOB.tree <- lmtree(Y ~ T.effect - 1 | x1 + x2 + x3, data = data, mtry=3, alpha=0.05, bonferroni=FALSE)
print(predMOB.tree)
plot(predMOB.tree)

#########################################################

### predMOB forest
ntree <- 100      # number of trees in the forest
nvarspersplit <- 2 # number of candidate variables to select for each split

# prepare manual subsampling
# output: list of 0/1 vectors of length n (1: used in subsample for construction of the tree, 0: out-of-bag sample)
sampleweights <- lapply(1:ntree, FUN=function(x, data){
  i <- sample(1:nrow(data), size=ceiling(0.632*nrow(data)), replace=FALSE)
  boot.wt <- as.numeric(1:nrow(data) %in% i)
  
  return(boot.wt)
}, data=data)

# grow forest 
# grow one single tree per subsample defined by one element of list 'sampleweights'
# and calculate permutation importance for this tree based on out-of-bag sample
grow.predMOB.forest <- lapply(sampleweights, FUN=function(weights, data, nvarspersplit){
  data.train <- data[which(weights==1),]
  data.oob <- data[which(weights==0),]
  
  fit.tree <- lmtree(Y ~ T.effect - 1 | x1 + x2 + x3, data = data.train, mtry=nvarspersplit, alpha=0.05, bonferroni=FALSE)
  
  permimp <- VI.lmtree(fit.tree, data.oob)
  
  return(list('tree'=fit.tree, 'permimp'=permimp))
}, data=data, nvarspersplit=nvarspersplit)

# extract forest as list of trees
predMOB.forest <- lapply(grow.predMOB.forest, FUN=function(x) x[['tree']])

# extract permutation importance and average over all trees
permimp <- colMeans(rbindlist(lapply(grow.predMOB.forest, FUN=function(x) data.frame(t(x[['permimp']])))))

# determine mean minimal depth over all trees
mindepth <- mindepth.mob.forest(predMOB.forest, mean_sample="all_trees")


permimp

#####################################################

### R Version
# sessionInfo()
# R version 4.0.2 (2020-06-22)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 18363)
# 
# Matrix products: default
# 
# Random number generation:
#   RNG:     Mersenne-Twister 
# Normal:  Inversion 
# Sample:  Rounding 
# 
# locale:
#   [1] LC_COLLATE=German_Germany.1252  LC_CTYPE=German_Germany.1252    LC_MONETARY=German_Germany.1252
# [4] LC_NUMERIC=C                    LC_TIME=German_Germany.1252    
# 
# attached base packages:
#   [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] randomForestExplainer_0.10.1 dplyr_1.0.2                  data.table_1.13.0            model4you_0.9-5             
# [5] partykit_1.2-9               mvtnorm_1.1-1                libcoin_1.0-6                TH.data_1.0-10              
# [9] MASS_7.3-51.6                survival_3.1-12             
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.5         plyr_1.8.6         RColorBrewer_1.1-2 pillar_1.4.6       compiler_4.0.2     tools_4.0.2       
# [7] digest_0.6.25      rpart_4.1-15       lifecycle_0.2.0    tibble_3.0.3       gtable_0.3.0       lattice_0.20-41   
# [13] pkgconfig_2.0.3    rlang_0.4.7        Matrix_1.2-18      GGally_2.0.0       rstudioapi_0.11    ggrepel_0.8.2     
# [19] gridExtra_2.3      htmlwidgets_1.5.1  generics_0.0.2     vctrs_0.3.2        DT_0.15            tidyselect_1.1.0  
# [25] reshape_0.8.8      glue_1.4.1         R6_2.4.1           Formula_1.2-3      ggplot2_3.3.2      purrr_0.3.4       
# [31] magrittr_1.5       htmltools_0.5.0    scales_1.1.1       ellipsis_0.3.1     splines_4.0.2      colorspace_1.4-1  
# [37] sandwich_2.5-1     munsell_0.5.0      inum_1.0-1         crayon_1.3.4       zoo_1.8-8         
 
