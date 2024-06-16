##################################
##################################
# Code for INFO+
##################################
##################################

################################
# Load libraries and code
################################
library(MASS) # To generate synthetic data by sampling a Multivariate Normal
library(infotheo) # Information theoretic library  
recent_path <-paste0(getwd(),"/Code")
# source(paste0(recent_path,"/2018-Bioinformatics-Predictive-Biomarker-Discovery-master/Functions-GenerateData.R")) # Function to generate synthetic data
# source(paste0(recent_path,"/2018-Bioinformatics-Predictive-Biomarker-Discovery-master/InformationTheory-PredictiveRankings.R")) # Functions to derive predictive rankings
source(paste0(recent_path,"/Methods/Code_inherited_Adjusted/InformationTheory-PredictiveRankings.R")) # Function to generate synthetic data
source(paste0(recent_path,'/Methods/Help_Functions.R'))

##################################
# Wrapper Function for INFO+
##################################

Info_Plus <- function(X,trt, y, top_k=NULL){ 
  # X:      biomarker data, of type data frame with biomarkers represented as columns; binary biomarker are of type factor, continuous of type numeric
  # trt:    treatment variable as a vector, of type integer
  # y:      response as a vector, of type factor if binary, of type numeric if continuous
  # top_k:  the number of biomarkers, the score should be computed for; only for the top-k parameters the score is returned
  
  # error in case provided data is not supported by the original function
  error = check_data_types(X=X, y=y, method="INFO_plus")
  if (error==TRUE){
    return(NULL)
  }
  
  # response and binary biomarkers are converted to an integer
  data_prep <- convert_variables(y=y,X=X) 
  y=data_prep$y
  X=data_prep$X
  # same results in case of using
  # INFOplus.Output_Categorical.Covariates_Continuous (binary X as integer, factor not allowed)
  # INFOplus.Output_Categorical.Covariates_Categorical (binary X as factor or integer)
  # for continuous biomarkers

  # define that function of infotheo should be used and not of the recipes package
  # otherwise an error occurs when running
  discretize <- infotheo::discretize 

  # since scores for all biomarkers should be calculated, top_k is set accordingly
  if (is.null(top_k)){
    top_k = ncol(X)
  }
  
  # capture start time
  start = Sys.time()
  # calling of original function INFO+, which captures second order interactions (returns the top_k = 5 biomarkers)
  infoplus <- INFOplus.Output_Categorical.Covariates_Continuous(X,y,trt,top_k) #$ranking # this function returns the ranking
  
  # capture end time 
  end = Sys.time()
  
  # preparation of scores
  ranking <- infoplus$ranking_scores
  scores <- infoplus$scores
  names(ranking) <- colnames(X)
  names(scores) <- colnames(X)
  # $scores: score of each feature in the sequence as they are in the dataset 
  # $ranking: the ranking of the scores from lower to higher e.g. 4 5 1 = 4th BM has the lowest score, 1 the highest
  # $ranking_scores: ranking of each feature in the sequence as they are in the dataset
  
  # calculate run time
  runtime = as.numeric(difftime(end, start, units='mins')) # only minutes are returned
  
  return(list("ranking"= ranking, "scores"= scores, "runtime"=runtime))

}
