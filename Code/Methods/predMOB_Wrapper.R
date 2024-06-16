######################################
######################################
# predMOB
######################################
######################################
# NOTE:
# the help functions of predMOB as well as the application of predMOB itself was 
# provided by the author of predMOB, J. Krzykalla;
# the wrapper function was developed during this master thesis

################################
# Load libraries and code
################################
require("partykit")
require("model4you") # only required for objfun function to calculate permutation importance 
require(data.table) # only for rbindlist function
require(dplyr) # used to determine mean minimal depth
require(randomForestExplainer) # used to determine mean minimal depth
recent_path <-paste0(getwd(),"/Code")
source(paste0(recent_path,'/Methods/Help_Functions.R'))
source(paste0(recent_path,'/Methods/Code_inherited_adjusted/helpfunctions_predMOB.R'))

#####################################
# Wrapper Function for predMOB
#####################################

PredMOB_Wrapper <- function(X, trt, y, ntree=100, nvarspersplit=2){
  # X:              biomarker data, of type data frame with biomarkers represented as columns; binary biomarker are of type factor, continuous of type numeric
  # trt:            treatment variable as a vector, of type integer
  # y:              response as a vector, of type factor if binary, of type numeric if continuous
  # ntree:          number of trees in the forest
  # nvarspersplit:  number of candidate variables to select for each split
  
  
  # error in case provided data is not supported by the original function
  error = check_data_types(X=X, y=y, method="predMOB")
  if (error==TRUE){
    return(NULL) 
  }
  
  # define effect-coded treatment variable (experimental group: +0.5, control group: -0.5)
  trt_eff <- ifelse(trt==1, 1, -1)/2  
  
  # combine data
  data <- cbind(X, trt_eff, y)
  colnames(data) <- c(colnames(X), "trt_eff", "y") # trt_eff

  ##########################################
  # based on code provided by  J. Krzykalla
  
  # ### single predMOB tree
  # first part of the formula defines base model, 
  # potential predictive factors follow separated by a pipe symbol (cf. partykit:::mob)
  # CAVE: the formula for the base model must not include an intercept
  fo <- as.formula(paste("y", paste("trt_eff - 1 |", paste(colnames(X), collapse=" + ")), sep=" ~ "))

  ### predMOB forest
  # parameters can be provided when calling the function
  # ntree <- 100 # number of trees in the forest
  # nvarspersplit <- 2 # number of candidate variables to select for each split
  
  # capture start time
  start = Sys.time()
  
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
    
    fit.tree <- lmtree(fo, data = data.train, mtry=nvarspersplit, alpha=0.05, bonferroni=FALSE)
    # example tree:
    # fit.tree <- lmtree(Y ~ T.effect - 1 | x1 + x2 + x3, data = data.train, mtry=nvarspersplit, alpha=0.05, bonferroni=FALSE)
    
    permimp <- VI.lmtree(fit.tree, data.oob)
    
    return(list('tree'=fit.tree, 'permimp'=permimp))
  }, data=data, nvarspersplit=nvarspersplit)
  
  # extract forest as list of trees
  predMOB.forest <- lapply(grow.predMOB.forest, FUN=function(x) x[['tree']])
  
  # extract permutation importance and average over all trees
  permimp <- colMeans(rbindlist(lapply(grow.predMOB.forest, FUN=function(x) data.frame(t(x[['permimp']])))))
  
  # determine mean minimal depth over all trees
  # not required in this case
  # mindepth <- mindepth.mob.forest(predMOB.forest, mean_sample="all_trees")
  
  # capture end time and calculate run time
  end = Sys.time()
  runtime = as.numeric(difftime(end, start, units='mins')) # only minutes are returned
  
  return(list("scores"=permimp, "runtime"=runtime)) # vector with bm names as names
}

