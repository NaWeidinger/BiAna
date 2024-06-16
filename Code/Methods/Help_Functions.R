################################
# Load libraries and code
################################
library(dplyr) # for editing data frames e.g. mutate, select

################################################
# Function for conversion of X, trt, y to 
#respective variable types required by a method
################################################

convert_variables<- function(X=NULL, trt=NULL, y=NULL){
  # X: biomarkers as data frame
  # trt: treatment as vector
  # y: response as vector
  # only if the respective parameters are provided to the function, a conversion is performed
  
  if(!is.null(trt)){
    # convert trt to factor
    trt = as.factor(trt)
  }  
  
  if(!is.null(y)){
    if (is.factor(y)){
      # convert y to integer only if binary (factor)
      y = as.integer(as.character(y))
    }
  }
  
  if(!is.null(X)){
    # convert X
    # select biomarkers that are factors and convert them to integer
    X <- X %>% mutate_if(is.factor, as.character)
    X <- X %>% mutate_if(is.character, as.integer)
    
  }
  
  return(list("X"=X, "trt"=trt, "y"=y))
}




######################################################
# Function to check if dataset / simulation 
# setting is supported by the method
# (if not, a warning is shown and return value 
# of the method is NULL)
#######################################################
check_data_types <- function(X, y, method){
  # X:      biomarkers as data frame
  # y:      response as vector
  # method: method for which the combination of scale type of X and y should be verified
  
  error = FALSE
  
  if (method=="VT" && !is.factor(y)){
    warning("Continuous response (y) not allowed/implemented for Virtual Twins.")
    error= TRUE
  } else if (method=="INFO_plus"){
    if (!is.factor(y)){
      warning("Continuous response (y) not allowed/implemented for INFO+.")
      error= TRUE
    } 
  }  else if (method=="KO_pred_CF" || method=="KO_pred_LI" ){
      if (any(sapply(X, is.factor))){
        # no factors allowed in X
        warning("Binary covariates not allowed/implemented for predictive knockoffs.")
        error= TRUE
      }
      if (method=="KO_pred_CF" && is.factor(y)){
        warning("Binary response not allowed/implemented for predictive knockoffs with CF filter.")
        error= TRUE
      }
 
  } else if (method=="predMOB" && is.factor(y)){
    warning("Binary response not allowed/implemented for predMOB.")
    error= TRUE
  } else if (method=="PPLasso"){
    if (is.factor(y)){
      warning("Binary response not allowed/implemented for PPLasso.")
      error= TRUE
    } 
    if (any(sapply(X, is.factor))){
      warning("Binary covariates not allowed/implemented for PPLasso.")
      error= TRUE
    }
  }
  return(error)
}

#########################################################
# Function for adding interaction effects to the data
# (BM*TRT)
#########################################################
Add_interact_dataset <- function(X, trt){
  # X:    covariates with continuous as numeric and binary as factor variable type
  # trt:  treatment variable as integer or factor variable type
  # trt is a integer when provided to the wrapper function of VSURF, but is converted for VSURF to a factor.
  # If this function is used for within another wrapper function, trt should always be provided as 
  # a factor or integer variable type
  
  # trt should be integer or factor
  conv_trt = FALSE
  # convert trt to factor if not already a factor (if integer)
  if (!is.factor(trt)){
    conv_trt = TRUE # indicator for later back-transformation
    trt = as.factor(trt)
  }

  # combine X and trt for adding interaction covariates
  data <- cbind(X, trt) 

  # since at least trt is a factor variable, 
  # factor columns are converted to integer
  col_ind_fact = which(sapply(data, is.factor))
  data <- data %>% mutate_at(names(col_ind_fact), as.character)
  data <- data %>% mutate_at(names(col_ind_fact), as.integer)

  # apply model.matrix in order to add interaction variables with trt
  data =as.data.frame(model.matrix(~ .*trt+0, data))

  # convert previous factor variables to factor again
  data <- data %>% mutate_at(names(col_ind_fact), as.factor)
  
  # interaction covariates with  trt and another factor covariate, should have the 
  # variable type factor again. Thus the respective interaction covariates need to be converted.
  for (i in 1:length(col_ind_fact)){
    if (!identical(names(col_ind_fact[i]),"trt")){
      col_ind_inter <- paste0(names(col_ind_fact[i]),":trt")
      data <- data %>% mutate_at(col_ind_inter, as.factor)
    }
  }

  # if trt was previously an integer, a conversion should take place
  if (conv_trt==TRUE){
    data$trt = as.integer(as.character(data$trt))
  }
  
  # provide X combined with trt and the interaction covariates as data frame
  X_interact <- as.data.frame(data)
  # only trt is converted to integer, the other factor (interaction) covariates remain as factors
  # and have to be processed according to functions they are used in
  
  # return data frame with X, trt and interaction covariates
  return(X_interact)
}

######################################################
# Function for filterin to placebo subjects 
######################################################

Return_pbo_data <- function(X, trt, y){
  # X: biomarkers as data frame
  # trt: treatment as vector
  # y: response as vector

  df= data.frame(X, trt, y)
  df <-df[(df$trt==0),]
  rownames(df) <- NULL
  
  trt <- df$trt 
  X <- df %>% dplyr::select(-contains(c("trt", "y")))
  y <- df$y
 
  return(list("X"=X, "trt"=trt, "y"=y))
}


