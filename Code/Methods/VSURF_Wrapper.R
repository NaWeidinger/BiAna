##################################
##################################
# Code for VSURF
##################################
##################################

################################
# Load libraries and code
################################
library(VSURF)
library(stringr) # for str_remove
library(doParallel)
library(dplyr) # for editing data frames e.g. mutate, select
recent_path <-paste0(getwd(),"/Code")
source(paste0(recent_path,'/Methods/Help_Functions.R'))

#############################################
# Wrapper Function for VSURF
#############################################

VSURF_Wrapper <- function(X, trt, y, bm_type="predictive", cores=NULL){
  # X:        biomarker data, of type data frame with biomarkers represented as columns; binary biomarker are of type factor, continuous of type numeric
  # trt:      treatment variable as a vector, of type integer
  # y:        response as a vector, of type factor if binary, of type numeric if continuous
  # bm_type:  determines if predictive or prognostic biomarkers should be identified, either "prognostic" or "predictive"
  # cores:    the number of cores the original function is executed on; can be specified; otherwise it is set to number of available cores - 1
  
  # conversion of trt to factor
  data_prep <- convert_variables(trt=trt)
  trt=data_prep$trt
  
  if (bm_type=="predictive")  { 
    # in case  predictive biomarkers should be identified, interaction terms of the biomarkers with trt have to be added to X
    X <- Add_interact_dataset(X, trt)
  } else {
    # otherwise - in case prognostic biomarkers should be identified, the data has to be restricted to placebo subjects only
    pbo_data <- Return_pbo_data(X, trt, y)
    X= pbo_data$X
    y=pbo_data$y
  }
  
  # if VSURF should be run in parallel
  if (is.null(cores)){
    cores = detectCores()-1
  }
  
  # capture start time
  start = Sys.time()
  
  # compute results, call original function
  vsurf <- VSURF(X, y, parallel = TRUE, ncores=cores) #A logical indicating if you want VSURF to run in parallel on multiple cores (default to FALSE)
  
  # capture end time and calculate run time
  end = Sys.time()
  runtime = as.numeric(difftime(end, start, units='mins')) # only minutes are returned
  
  
  # construct selected attribute
  sel_ind <-vsurf$varselect.pred
  selected <- rep(0,ncol(X))
  selected[sel_ind] <- 1
  names(selected) <- colnames(X)
  selected_interpretation <- rep(0,ncol(X))
  selected_interpretation[vsurf$varselect.interp] <- 1
  names(selected_interpretation) <- colnames(X)
  
  # construct scores
  all_scores <- rep(0, ncol(X))
  all_scores[vsurf$imp.mean.dec.ind] <- vsurf$imp.mean.dec
  names(all_scores) <- colnames(X)

  if (bm_type=="predictive"){
    # for predictive effects, only the interaction effects have to be regarded; thus filtering is required
    id_cols <- grep(":trt", colnames(X))
    all_scores <- all_scores[id_cols]
    all_selected <- selected[id_cols]
    selected_interpretation <- selected_interpretation[id_cols]
    names(all_scores) <- str_remove(names(all_scores), ":trt")
    names(all_selected) <- str_remove(names(all_selected), ":trt")
    names(selected_interpretation) <- str_remove(names(selected_interpretation), ":trt")
    selected = all_selected
  } else {
    names(all_scores) = colnames(X)
    names(selected) = colnames(X)
    names(selected_interpretation) = colnames(X)
  }
  
  # score - threshold can be used as "new" score
  scores_thres <- all_scores - vsurf$min.thres
  
  return(list("vsurf_element"=vsurf, "selected" = selected, "scores" = all_scores,"threshold"=vsurf$min.thres, "scores_thres"=scores_thres,"runtime"=runtime, "selected_interpretation"=selected_interpretation)) #
}

