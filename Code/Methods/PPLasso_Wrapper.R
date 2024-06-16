##################################
##################################
# Code for PPLasso
##################################
##################################

################################
# Load libraries and code
################################
library(PPLasso)
library(cvCovEst)
library(ggplot2)
require(gridExtra)
recent_path <-paste0(getwd(),"/Code")
source(paste0(recent_path,'/Methods/Help_Functions.R'))

#############################################
# Wrapper Function for PPLasso
#############################################

execute_PPLasso <- function(X, trt, y){
  # X:    biomarker data, of type data frame with biomarkers represented as columns; binary biomarker are of type factor, continuous of type numeric
  # trt:  treatment variable as a vector, of type integer
  # y:    response as a vector, of type factor if binary, of type numeric if continuous
  
  # error in case provided data is not supported by the original function
  error = check_data_types(X=X, y=y, method="PPLasso")
  if (error==TRUE){
    return(NULL) 
  }
  
  dataset <- data.frame(X, trt, y)
  p=ncol(X)
  
  # sort dataset with regards to treatment (first group should be placebo), as required
  # for the original function
  dataset <- dataset[order(dataset$trt),] #ascending sorting (placebos = 0 first)
 
  # create bm dataset without y and trt
  X_bm=dataset[,!names(dataset) %in% c("y", "trt")]
  trt=dataset[,names(dataset) %in% c("trt")]
  y=dataset[,names(dataset) %in% c("y")]
  
  # create X1 and X2 and TRT1 and TRT2 separated for the two treatment groups
  X1 = filter(dataset, trt == 0) # receiving standard treatment / placebo
  X1=X1[,!names(X1) %in% c("y", "trt")]
  
  X2 =  filter(dataset, trt == 1) 
  X2=X2[,!names(X2) %in% c("y", "trt")]
  
  n1 = nrow(X1)
  n2 = nrow(X2)

  # conversions to matrices
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  y <- as.matrix(y)
  X_bm <- as.matrix(X_bm)
  
  # Help function for estimation of Sigma
  cv_cov_est_out <- cvCovEst(
    dat = X_bm,
    estimators = c(
      linearShrinkLWEst, denseLinearShrinkEst,
      thresholdingEst, poetEst, sampleCovEst
    ),
    estimator_params = list(
      thresholdingEst = list(gamma = c(0.2, 0.4)),
      poetEst = list(lambda = c(0.1, 0.2), k = c(1L, 2L))
    ),
    cv_loss = cvMatrixFrobeniusLoss,
    cv_scheme = "v_fold",
    v_folds = 5, 
    parallel= TRUE 
  )
  
  Sigma_est <- cov2cor(cv_cov_est_out$estimate)
  
  # capture start time 
  start = Sys.time()
  
  # compute results, call original function
  mod <- ProgPredLasso(X1 = X1, X2 = X2, Y = y, cor_matrix = Sigma_est)
  # X1 rows related to treatment t0
  # X2 rows related to treatment t1
  # When t0 stands for the standard treatment (placebo), prognostic biomarkers are defined as those having non-zero coefficients in beta1
  # and non prognostic (resp. non predictive) biomarkers correspond to the indices having null coefficients in beta1
  # When t0 stands for the standard treatment (placebo), predictive biomarkers are defined as those having non-zero coefficients in 
  # resp. in beta1-beta2 and non predictive biomarkers correspond to the indices having null coefficients in beta1-beta2
  
  # capture end time and calculate run time
  end = Sys.time()
  runtime = as.numeric(difftime(end, start, units='mins')) # only minutes are returned
  
  # identified as prognostic:
  prog_coeff <- mod$beta.min[1:p]
  prog_ind <- which(mod$beta.min[1:p]!=0)
  # identified as predictive
  pred_coeff <- mod$beta.min[(p+1):(2*p)]
  pred_ind <- which(mod$beta.min[(p+1):(2*p)]!=0)

  # build selected attribute
  selected_prog <- rep(0, ncol(X_bm))
  selected_pred <- rep(0, ncol(X_bm))
  selected_prog[prog_ind] <- 1
  selected_pred[pred_ind] <- 1

  # naming of the variables
  names(selected_prog) = colnames(X)
  names(selected_pred) = colnames(X)
  names(prog_coeff) = colnames(X)
  names(pred_coeff) = colnames(X)

  return(list("selected_prog" = selected_prog, "selected_pred" = selected_pred, "scores_prog"=prog_coeff, "scores_pred"=pred_coeff, "runtime"=runtime))
}

