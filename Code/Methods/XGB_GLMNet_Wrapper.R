######################################
######################################
# Code for XGB and GLMNet approaches
######################################
######################################

################################
# Load libraries and code
################################
library(doParallel)
library(xgboost)
library(flashlight)
library(MachineShop)
library(stringr) # for str_remove
registerDoParallel(cores = 4) # recommended: 32 cores
recent_path <-paste0(getwd(),"/Code")
source(paste0(recent_path,'/Methods/Help_Functions.R'))


##########################################
# Function for defining tuning 
# parameter for XGB and GLMNet approaches
##########################################

XGB_GLM_Tune_params <- function(method="XGB", rand_grid_size=2, folds=1, repeats=1){ #50, 5,3
  # method:         either "XGB" or "GLMNet"
  # rand_grid_size: size of the random grid
  # folds:          number of folds for the cross validation procedure used within the hyperparamter tuning
  # repeats:        number of repeats for the cross validation procedure used within the hyperparamter tuning
  
  # control is cross-validation and tuning set is automatically created with TuningGrid
  control <- MachineShop::CVControl(folds = folds, repeats = repeats) 
  rand_grid =  MachineShop::TuningGrid(size = 500, random = rand_grid_size)
  
  # for XGB XGBTreeModel is selected and the respective model is defined
  if (method=="XGB"){
    tune_model<-  MachineShop::TunedModel( 
      XGBTreeModel, 
      grid =  rand_grid,
      control = control
    )
  } else {
    # for GLMNet GLMNetModel is selected and the respective model is defined
    tune_model<-  MachineShop::TunedModel(
      GLMNetModel,
      grid =  rand_grid, 
      control = control
    )
  }
  
  return(list("tune_model"=tune_model,"control"=control)) 
}

###############################################
# Function for execution of the hyperparameter 
# tuning and construction of the final model
###############################################

HP_search <- function(dataset, tune_model) {
  # dataset: data frame containing of the biomarkers (X), the treatment (trt) and the response (y)
  # tune_model: model used for tuning, usually from XGB_GLM_Tune_params
  
  model_fit  <- MachineShop::fit(y ~ ., data = dataset,
                                 model = tune_model)
  
  # only required if parameters are desired to be provided
  ml_mod <-as.MLModel(model_fit) #tuned model
  ml_mod_sum <- summary(ml_mod)
  model_possibilities <- ml_mod_sum$TrainingStep1
  selected_model <- model_possibilities[model_possibilities$selected == TRUE,]

  return(list("params" = selected_model$params, "model_fit"=model_fit)) 
}

#################################################
# Wrapper Function for XGB and GLMNet approaches
#################################################

execute_GLM_XGB <- function(X, trt, y, type_y, method="XGB", approach="all", rand_grid_size=3, folds=2, repeats=1){ #NWe 50, 5,3
  # X:              biomarker data, of type data frame with biomarkers represented as columns; binary biomarker are of type factor, continuous of type numeric
  # trt:            treatment variable as a vector, of type integer
  # y:              response as a vector, of type factor if binary, of type numeric if continuous
  # rand_grid_size: size of the random grid
  # folds:          number of folds for the cross validation procedure used within the hyperparamter tuning
  # repeats:        number of repeats for the cross validation procedure used within the hyperparamter tuning
  
  # convert binary biomarkers to integers
  data_prep <- convert_variables(X=X)
  X=data_prep$X
  
  if (method=="GLMNet" && (approach=="Hstat" || approach=="Shap_inter")){
    stop("GLMNet cannot be combined with H-statistic or Shapley interaction values. No framework available. Please choose another method/approach.")
  }
  
  # construct data frame and define X as biomarkers with trt variable
  dataset <- data.frame(X,trt,y) 
  X=dataset[,!names(dataset) %in% c("y")]
  n <- nrow(X)
  m <- ncol(X)
  cols <- colnames(X)
  
  # construct tuning parameters
  tuning_params <- XGB_GLM_Tune_params(method=method, rand_grid_size=rand_grid_size, folds=folds, repeats=repeats)
  
  ###############################################
  # approaches for predictive biomarkers
  ###############################################
  
  # in case XGB with H-Statistic and/or Shap interaction values should be computed, the dataset can be taken for modelling as it is
  if (method=="XGB" && (approach=="all" ||approach=="Hstat"|| approach=="Shap_inter")){
    
    # capture start time
    start_pred = Sys.time()
    
    # compute results, call original function
    selected_model <- HP_search(dataset, tuning_params$tune_model)
    
    # capture end time and calculate run time
    end_pred = Sys.time()
    runtime_general_pred = as.numeric(difftime(end_pred, start_pred, units='mins')) 
    
    ##############################
    # XGB: Shapl interaction
    
    # if Shap interaction values should be computed
    if (approach=="all" || approach=="Shap_inter"){
      
      # capture start time for Shap interaction values
      start_shapinter = Sys.time()

      # computation of the Shap interaction values
      pred_inter <- xgboost:::predict.xgb.Booster(selected_model$model_fit, as.matrix(X), predinteraction=TRUE)
      
      # capture end time and calculate run time
      end_shapinter = Sys.time()
      runtime_shap_inter = as.numeric(difftime(end_shapinter, start_shapinter, units='mins'))
      
      # in case of a binary response, the results are given per class
      # consequently, the results have to be filtered for the second class
      if (type_y=="binary"){
        pred_inter <- pred_inter[[2]]
      }
      
      # calculate the mean shap interaction values
      mean_shap_pre <- apply(abs(pred_inter), c(2,3), mean)
      # set off-diagonals * 2 and remove bias
      mean_shap <- mean_shap_pre*2
      diag(mean_shap) <- diag(mean_shap_pre)
      bias = mean_shap[nrow(mean_shap),ncol(mean_shap)]
      mean_shap <- mean_shap[-nrow(mean_shap),-ncol(mean_shap)] #remove bias
      # filter for relevant variables: interaction of a biomarker with the treatment
      shap_interaction <- as.data.frame(mean_shap["trt",])
      # naming of the vector and further preparations
      colnames(shap_interaction) <- "shap_inter"
      shap_interaction$bm <- rownames(shap_interaction)
      rownames(shap_interaction) <- 1:nrow(shap_interaction)
      shap_interaction <- shap_interaction[,c(2,1)] 
      shap_interaction <- shap_interaction[-nrow(shap_interaction),] 
      names <- shap_interaction$bm
      shap_interaction <- shap_interaction$shap_inter 
      names(shap_interaction) <- names
      print(shap_interaction)
      
    # in case Shap interaction values should not be computed, the respective variables are set to NULL
    } else {
      shap_interaction=NULL
      runtime_shap_inter=NULL
    }

    ##############################
    # XGB: H-Statistic

    # if H-Statistic values should be computed
    if (approach=="all" ||approach=="Hstat"){

      # define custom predict function, depending on type of the response, applicable for H-Stastic and Shap values
      if (type_y=="binary"){
        custom_pred_fct <- function(m,X) {
          predict(m, data.matrix(X[, cols]), type="prob")
        }
      } else {# assumed to be continuous then
        custom_pred_fct <- function(m,X) {
          predict(m, data.matrix(X[, cols]))
        }
      }
      
      # capture start time for H-Statistic
      start_hstat= Sys.time()
      
      # compute H-Statistic
      fl <- flashlight(model = selected_model$model_fit, data = data.frame(y, X), y = "y", label = method,
                       predict_fun = custom_pred_fct)
      
      imp <- light_interaction(fl, v = cols, pairwise = TRUE)
      
      # capture end time and calculate run time
      end_hstat = Sys.time()
      runtime_Hstat = as.numeric(difftime(end_hstat, start_hstat, units='mins'))
      
      # filter for interactions of biomarkers with treatment
      h_interact_trt <- filter(imp$data, grepl("trt", variable, fixed = TRUE)) 
      scores <- as.data.frame(h_interact_trt)
      scores$variable <- str_remove(scores$variable, ":trt")
      scores=scores[,!names(scores) %in% c("label", "error")]
      # preparation and naming of the vector
      Hstat <- data.frame(bm=scores$variable, Hstat=scores$value)
      order <- colnames(X[,!names(X) %in% c("trt")])
      Hstat <- Hstat[match(order, Hstat$bm),]
      names <- Hstat$bm
      Hstat <- Hstat$Hstat 
      names(Hstat) <- names
    # in case the H-Statistic should not be computed, the respective variables are set to NULL
    } else {
      Hstat =NULL
      runtime_Hstat =NULL
    }
  # in case GLMNet should be constructed, the respective variables are set to NULL
  } else {
    shap_interaction=NULL
    runtime_shap_inter=NULL
    Hstat =NULL
    runtime_Hstat =NULL
    runtime_general_pred = NULL
  }
  
  ###############################################
  # approaches for prognostic biomarkers
  ###############################################
  
  # in case Shap values and / or variable importance for GLMNet or XGB should be calculated (prognostic approaches)
  if (approach=="all" || approach=="Shap" || approach=="VI"){ # independent of method

    # restrict dataset to placebo subjects only
    X=X[,!names(X) %in% c("trt")]
    pbo_data <- Return_pbo_data(X, trt, y)
    X= pbo_data$X
    trt=pbo_data$trt
    y=pbo_data$y
    dataset <- data.frame(X,y) # dataset
    
    # capture start time for prognostic approaches
    start_prog = Sys.time()
    
    # run hyperparameter tuning and final model construction on restricted dataset
    selected_model <- HP_search(dataset, tuning_params$tune_model)
    
    # capture end time and calculate run time
    end_prog = Sys.time()
    runtime_general_prog = as.numeric(difftime(end_prog, start_prog, units='mins')) 
    
  } else {
    # if no prognostic approaches should be applied, set run time to NULL
    runtime_general_prog = NULL
  }
  
  n <- nrow(X) # probably not required
  m <- ncol(X)
  cols <- colnames(X)
  
  #####################
  # Shap mean

  # if Shap values should be calculated
  if (approach=="all" || approach=="Shap"){ 
    
    # define custom predict function, depending on type of the response, applicable for H-Statistic and Shap values
    if (type_y=="binary"){
      custom_pred_fct <- function(m,X) {
        predict(m, data.matrix(X[, cols]), type="prob")
      }
    } else { # assumed to be continuous then
      custom_pred_fct <- function(m,X) {
        predict(m, data.matrix(X[, cols]))
      }
    }
    # capture start time for Shap values
    start_shapmean = Sys.time()
    
    # compute Shap values
    fl <- flashlight(model = selected_model$model_fit, data = data.frame(y, X), y = "y", label = method,
                     predict_fun = custom_pred_fct)
    fl <- add_shap(fl)
    shap_mean <- light_importance(fl, type="shap")
    
    # capture end time and calculate run time for Shap values
    end_shapmean = Sys.time()
    runtime_shap_mean = as.numeric(difftime(end_shapmean, start_shapmean, units='mins')) 
    
    # preparation of the vector, e.g. order and naming
    shap <- data.frame(bm=shap_mean$data$variable, shap_mean= shap_mean$data$value)
    order <- colnames(X)
    shap <- shap[match(order, shap$bm),]
    rownames(shap) <- 1:nrow(shap)    # Assign sequence to row names
    names <- shap$bm
    shap <- shap$shap_mean 
    names(shap) <- names
    
  } else {
    # if no Shap values should be calculated, set variables to NULL
    shap=NULL
    runtime_shap_mean=NULL
  }
  
  #######################
  # variable importance

  # if Shap values should be calculated
  if (approach=="all" ||approach=="VI"){ # independent of method

    # capture start time for variable importance
    start_vi=Sys.time()
    
    # compute variable importance
    imp <- MachineShop::varimp(object=selected_model$model_fit)  #scale = FALSE?
    
    # capture end time and calculate run time for variable importance
    end_vi = Sys.time()
    runtime_VI = as.numeric(difftime(end_vi, start_vi, units='mins')) 
    
    # preparation
    imp_orig <- imp
    bm_names <- rownames(imp) #[-ind_trt]
    
    # select the respective column, depending on type of the response (different metrics)
    if (type_y=="binary"){
      imp <- imp$Permute.mean.brier
    } else {
      imp <- imp$Permute.mean.rmse
    }
    
    # preparation
    names(imp) <- bm_names
    order <- colnames(X)
    imp <- imp[match(order, names(imp))]

  } else {
    # if no variable importance should be calculated, set variables to NULL
    imp= NULL
    runtime_VI=NULL
  }

  return(list("varimp"=imp, "shap_mean"=shap, "shap_interaction"=shap_interaction, "Hstat"=Hstat, "runtime_general_pred"=runtime_general_pred,"runtime_general_prog"=runtime_general_prog, "runtime_VI"=runtime_VI,"runtime_shap_mean"=runtime_shap_mean,"runtime_shap_inter"=runtime_shap_inter,"runtime_Hstat"=runtime_Hstat))
  
}
