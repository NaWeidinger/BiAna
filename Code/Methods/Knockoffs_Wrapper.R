##################################
##################################
# Code for Knockoffs
##################################
##################################

################################
# Load libraries and code
################################
library(glmnet) # is needed for the linear interaction KO filter 
library(grf) # is needed for the causal forest (CF) variable importance KO filter 
library(knockoff) # is needed for generating the knockoffs and calculating the threshold
library(pracma) # for check isempty
recent_path <-paste0(getwd(),"/Code")
source(paste0(recent_path,'/Methods/Help_Functions.R'))

# for predictive knockoffs
source(paste0(recent_path,'/Methods/Code_inherited_adjusted/predictive_knockoff_filters.R'))
source(paste0(recent_path,'/Methods/Code_inherited_adjusted/seqknockoff/knockoff_filters.R'))

# for prognostic knockoffs
source(paste0(recent_path,'/Methods/Code_inherited_adjusted/seqknockoff/internal.R'))
source(paste0(recent_path,'/Methods/Code_inherited_adjusted/seqknockoff/performance.R'))
source(paste0(recent_path,'/Methods/Code_inherited_adjusted/seqknockoff/plot.R'))
source(paste0(recent_path,'/Methods/Code_inherited_adjusted/seqknockoff/knockoff_filters.R'))
source(paste0(recent_path,'/Methods/Code_inherited_adjusted/seqknockoff/simdata.R'))
source(paste0(recent_path,'/Methods/Code_inherited_adjusted/seqknockoff/simknockoffs.R'))


#############################################
# Wrapper Function for predictive knockoffs
#############################################
# Original R-Code: https://github.com/sechidis/2021-SiM-Predictive-Knockoffs
# R-Code was adjusted to also return the W-statistic and the threshold

Run_Knockoffs_predictive <- function(X, trt, y, response_type, filter="LI", fdr_nominal=0.2){
  # X: biomarker data, of type data frame with biomarkers represented as columns; binary biomarker are of type factor, continuous of type numeric
  # trt: treatment variable as a vector, of type integer
  # y: response as a vector, of type factor if binary, of type numeric if continuous
  # filter: filter that should be applied, either "LI" or "CF"
  # fdr_nominal: false discovery rate
  
  # error in case provided data is not supported by the original function, also depending on the filter
  if (filter=="CF"){
    error = check_data_types(X=X, y=y, method="KO_pred_CF")
    if (error==TRUE){
      # CF_filter==FALSE
      return(NULL) 
    }
  } else {
    error = check_data_types(X=X, y=y, method="KO_pred_LI")
    if (error==TRUE){
      # LI_filter==FALSE
      return(NULL) 
    }
  }
  
  # conversion of X from dataframe to matrix required for X tilde calculation
  X_tilde = data.frame(create.second_order(as.matrix(X)))  
  
  # a binary response have to be a factor according to documentation
  # set famlily_val according to the response scale type
  if (response_type=="continuous"){
    family_val ="gaussian"
  } else { 
    y <- as.factor(y)
    family_val ="binomial"
  }
  
  # capture start time 
  start = Sys.time()
  
  # compute results, call original function
  W_Statistic = NULL
  selected = NULL
  if (filter=="LI"){
    results = linear_model_predictive_filter(X, X_tilde, y, trt, fdr_nominal = fdr_nominal, family = family_val)
  } else { 
    family_val=NULL
    results = causal_forest_predictive_filter(X, X_tilde, y, trt, fdr_nominal = fdr_nominal, family = family_val)
  }
  
  # capture end time and calculate run time
  end = Sys.time()
  runtime = as.numeric(difftime(end, start, units='mins')) # only minutes are returned
  
  # save W-statistic
  W_Statistic <- results$W_Statistic
  
  # for selected attribute
  selected = rep(0, length(W_Statistic))
  # in case some biomarkers are identified as predictive, the selected attribute has to be set 
  if (isempty(results$selected_features_indices)==FALSE){
    selected[results$selected_features_indices] <- 1
  }
  
  # naming of the variables
  names <- colnames(X)
  names(selected) <- names
  names(W_Statistic) <- names
  
  # score - threshold can be used as "new" score
  scores_thres <- W_Statistic - results$threshold

  return(list("selected"=selected, "scores"=W_Statistic, "runtime"=runtime, "threshold"=results$threshold,"scores_thres"=scores_thres))
}

#############################################
# Wrapper Function for prognostic knockoffs
#############################################
# Original R-Code: https://github.com/kormama1/seqknockoff
# R-Code was adjusted to also return the W-statistic and the threshold

Run_Knockoffs_prognostic <- function(X, trt, y, y_type, fdr=0.2){
  # X:      biomarker data, of type data frame with biomarkers represented as columns; binary biomarker are of type factor, continuous of type numeric
  # trt:    treatment variable as a vector, of type integer
  # y:      response as a vector, of type factor if binary, of type numeric if continuous
  # y_type: scale type of the response; allowed for this wrapper function are "continuous", "binary"
  # fdr:    false discovery rate
  
  pbo_data <- Return_pbo_data(X, trt, y)
  X= pbo_data$X
  y=pbo_data$y

  # binary response have to be a factor according to documentation
  if (y_type=="continuous"){
    family_val ="gaussian"
  } else { 
    family_val ="binomial"
  }
  
  # capture start time 
  start = Sys.time()
  # compute results, call original function
  results <- knockoff_filter(X, y, fdr=fdr , family=family_val)
  
  # capture end time and calculate run time
  end = Sys.time()
  runtime = as.numeric(difftime(end, start, units='mins')) # only minutes are returned
  
  # save W-statistic
  W_Statistic <- results$W_Statistic
  
  # for selected attribute
  selected = rep(0, length(W_Statistic))
  # in case some biomarkers are identified as predictive, the selected attribute has to be set 
  if (isempty(results$selected_features_indices)==FALSE){
    selected[results$selected_features_indices] <- 1
  }
  
  # naming of the variables
  names <- colnames(X)
  names(selected) <- names
  names(W_Statistic) <- names
  
  # score - threshold can be used as "new" score
  scores_thres <- W_Statistic - results$threshold

  return(list("selected"=selected, "scores"=W_Statistic, "runtime"=runtime, "threshold"=results$threshold,"scores_thres"=scores_thres))
}
