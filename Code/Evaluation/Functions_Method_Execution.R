################################
# Load libraries and code
################################
library(dplyr) # for editing data frames e.g. mutate, select

################################################################
# Function for preprocessing of the data and method execution
################################################################
Execute_Methods <- function(dataset, methods=c("KO_prog","KO_pred","PPLasso","VT","INFO_plus","SIDES","predMOB","VSURF_inter","VSURF_prog", "XGB"), simulation=TRUE, file_name_results="results.csv"){
  # dataset:            prepared dataset as in Method_Execution.R
  # methods:            methods to be executed, see default values for possible values
  # simulation:         TRUE if the function is applied to simulations
  # file_name_results:  file name for saving the results
  
  if ("y_cont" %in% colnames(dataset)){
    dataset <- dataset %>% dplyr:::rename(y = y_cont)
  } else if ("y_binary" %in% colnames(dataset)){
    dataset <- dataset %>% dplyr:::rename(y = y_binary)
  }
  
  if (simulation == TRUE){
    sim_info <- paste(dataset$Complexity[1], dataset$Dataset_No[1], dataset$Input_Data_Type[1], sep = ",", collapse=NULL)
    # determine response type
    if (is.factor(dataset$y)){
      sim_info <- paste(sim_info, "binary", sep = ",", collapse=NULL)
    } else { # otherwise continuous
      sim_info <- paste(sim_info, "continuous", sep = ",", collapse=NULL)
    }
    dataset=dataset[,!names(dataset) %in% c("Complexity", "Dataset_No", "Input_Data_Type")]
  }

  # split of dataset required for the methods 
  X=dataset[,!names(dataset) %in% c("y", "trt")]
  trt=dataset[,names(dataset) %in% c("trt")]
  y=dataset[,names(dataset) %in% c("y")]
  
  # for each method, the function Call_each_method is executed (enables implementation of parallelization)
  for (method in methods){
    rows <- Call_each_method(X,trt,y, method, simulation=simulation, file_name_results=file_name_results, sim_info=sim_info)
    
    # results are written to the results file
    write.table(rows, file = file_name_results, sep = ",",
                append = TRUE, quote = FALSE,
                col.names = FALSE, row.names = FALSE)
    
  }
}
###########################################
# Function for execution of each method
###########################################
Call_each_method <- function(X, trt, y, method, simulation, file_name_results, sim_info=NULL){
  # X:                  biomarker data
  # trt:                treatment variable
  # y:                  response
  # method:             method to be executed
  # simulation:         TRUE if the function is applied to simulations
  # file_name_results:  file name for saving the results
  # sim_info:           information regarding the simulation scenario (will be written to the results file)

  if (simulation==TRUE && is.null(sim_info)){
    stop("If simulation is run, the information with regard to the scenario have to be provided to the function")
  } 

  if (is.factor(y)==TRUE){
    y_type="binary"
  } else {
    y_type= "continuous"
  }
  
  #######################
  # Prognostic Knockoffs 
  #######################
  if (method=="KO_prog"){
    ko_prog = Run_Knockoffs_prognostic(X=X, trt=trt, y=y, y_type=y_type, fdr=0.2) #simulation=simulation
    
    # in case of a threshold set to infinity, the scores and selected attribute are set to 0
    if (is.infinite(ko_prog$threshold)){
      ko_prog$selected <- rep(0,length(ko_prog$selected)) # this is probably redundant,
      # since comparison with Inf works in R
      ko_prog$scores <- rep(0, length(ko_prog$scores))
    } else {
      # otherwise scores - threshold is the new score
      ko_prog$scores <- ko_prog$scores - ko_prog$threshold
    }
    
    rows <- paste(paste(sim_info, collapse=""),paste(c("prognostic-Knockoffs", "prognostic",ko_prog$scores,ko_prog$selected,ko_prog$runtime), collapse=","), sep=",")
    return(rows)
    
    ######################
    # Knockoffs predictive
    ######################
  } else if (method=="KO_pred"){
    
    ko_pred_LI = Run_Knockoffs_predictive(X=X, trt=trt, y=y, response_type=y_type, filter="LI", fdr_nominal=0.2)
    ko_pred_CF = Run_Knockoffs_predictive(X=X, trt=trt, y=y, response_type=y_type, filter="CF", fdr_nominal=0.2)
    
    # in case of no results
    if (is.null(ko_pred_LI)){
      rowLI=NULL
    } else {
      # in case of a threshold set to infinity, the scores and selected attribute are set to 0
      if (is.infinite(ko_pred_LI$threshold)){
        ko_pred_LI$selected <- rep(0,length(ko_pred_LI$selected)) 
        ko_pred_LI$scores <- rep(0, length(ko_pred_LI$scores))
      } else {
        # otherwise scores - threshold is the new score
        ko_pred_LI$scores <- ko_pred_LI$scores - ko_pred_LI$threshold
      }
      rowLI<- paste(paste(sim_info, collapse=""),paste(c("predictive-Knockoffs (LI)", "predictive", ko_pred_LI$scores,ko_pred_LI$selected, ko_pred_LI$runtime), collapse=","), sep=",")
    }
    
    # in case of no results
    if (is.null(ko_pred_CF)){
      rowCF=NULL
    } else {
      # in case of a threshold set to infinity, the scores and selected attribute are set to 0
      if (is.infinite(ko_pred_CF$threshold)){
        # if threshold is infinity, set all selected to 0
        ko_pred_CF$selected <- rep(0,length(ko_pred_CF$selected)) # this is probably redundant,
        ko_pred_CF$scores <- rep(0, length(ko_pred_CF$scores))
      } else {
        # otherwise scores - threshold is the new score
        ko_pred_CF$scores <- ko_pred_CF$scores - ko_pred_CF$threshold
      }
      rowCF <- paste(paste(sim_info, collapse=""),paste(c("predictive-Knockoffs (CF)", "predictive", ko_pred_CF$scores,ko_pred_CF$selected, ko_pred_CF$runtime), collapse=","), sep=",")
    }    
    return(rbind(rowLI,rowCF))
    
    ######################
    # PPLasso
    ###################### 
  } else if (method=="PPLasso"){
    pplasso = execute_PPLasso(X=X, trt=trt, y=y)
    
    # in case of no results
    if (is.null(pplasso)){
      return(NULL)
    } else {
      row_pred <- paste(paste(sim_info, collapse=""),paste(c("PPLasso", "predictive", pplasso$scores_pred,pplasso$selected_pred, pplasso$runtime), collapse=","), sep=",")
      row_prog <- paste(paste(sim_info, collapse=""),paste(c("PPLasso", "prognostic", pplasso$scores_prog,pplasso$selected_prog, pplasso$runtime), collapse=","), sep=",")
      return(rbind(row_pred,row_prog))
    }
    
    ######################
    # Virtual Twins
    ###################### 
  } else if (method=="VT"){
    vt = VT_Wrapper(X=X, trt=trt, y=y) # NWe tree_type="one"
    
    # in case of no results
    if (is.null(vt)){
      return(NULL)
    } else {
      rows <- paste(paste(sim_info, collapse=""),paste(c("VirtualTwins", "predictive", vt$scores,rep(NA,ncol(X)), vt$runtime), collapse=","), sep=",")
      return(rows)
    }

    ######################
    # INFO+
    ###################### 
  } else if (method=="INFO_plus"){
    Info_Plus = Info_Plus(X=X,trt=trt, y=y, top_k=NULL) #y_type=y_type, 
    
    # in case of no results
    if (is.null(Info_Plus)){
      return(NULL)
    } else {
      rows <- paste(paste(sim_info, collapse=""),paste(c("INFOplus", "predictive", Info_Plus$scores,rep(NA,ncol(X)), Info_Plus$runtime), collapse=","), sep=",")
      return(rows)
    }
   
    ######################
    # SIDES
    ###################### 
  } else if (method=="SIDES"){
    
    # define biomarker scale types for SIDES execution
    type_X = rep("continuous", ncol(X))
    type_X[which(sapply(X,is.factor)==TRUE)] <- "nominal"

    sides = execute_SIDES(X=X, trt=trt, y=y, type_X=type_X, type_y=y_type)
    rows <- paste(paste(sim_info, collapse=""),paste(c("SIDES", "predictive", sides$scores, sides$selected, sides$runtime), collapse=","), sep=",")

    return(rows)
    
    ######################
    # PredMOB
    ###################### 
  } else if (method=="predMOB"){
    
    pred_mob = PredMOB_Wrapper(X=X, trt=trt, y=y, ntree=100, nvarspersplit=2)
    
    # in case of no results
    if (is.null(pred_mob)){
      return(NULL)
    } else {
      rows <- paste(paste(sim_info, collapse=""),paste(c("PredMOB", "predictive", pred_mob$scores,rep(NA,ncol(X)), pred_mob$runtime), collapse=","), sep=",")
      return(rows)
    }
    
    ######################
    # VSURF interactions
    ######################  
  } else if (method=="VSURF_inter"){
    vsurf_inter = VSURF_Wrapper(X, trt, y, bm_type="predictive", cores=32) #, type="predictive"
    rows <- paste(paste(sim_info, collapse=""),paste(c("VSURF-interaction", "predictive", (vsurf_inter$scores-vsurf_inter$threshold),vsurf_inter$selected, vsurf_inter$runtime), collapse=","), sep=",")
    return(rows)

    ######################
    # VSURF
    ###################### 
  } else if (method=="VSURF_prog"){
    vsurf = VSURF_Wrapper(X, trt, y, bm_type="prognostic", cores=32) #, type="predictive"
    rows <- paste(paste(sim_info, collapse=""),paste(c("VSURF", "prognostic", (vsurf$scores-vsurf$threshold),vsurf$selected,vsurf$runtime), collapse=","), sep=",")
    return(rows)

    ######################
    # XGB
    ###################### 
  } else if (method=="XGB"){
    xgb = execute_GLM_XGB(X, trt, y, type_y=y_type, method="XGB", approach="all")#,approach="Shap_inter" #NWe calc_metrics=FALSE,
   
    row1 <- paste(paste(sim_info, collapse=""),paste(c("XGB-Varimp", "prognostic", xgb$varimp,rep(NA,ncol(X)),(xgb$runtime_general_prog+xgb$runtime_VI)), collapse=","), sep=",")
    row2 <- paste(paste(sim_info, collapse=""),paste(c("XGB-Shap", "prognostic", xgb$shap_mean,rep(NA,ncol(X)),(xgb$runtime_general_prog+xgb$runtime_shap_mean)), collapse=","), sep=",")
    row3 <- paste(paste(sim_info, collapse=""),paste(c("XGB-Shap-inter", "predictive", xgb$shap_interaction,rep(NA,ncol(X)),(xgb$runtime_general_pred+xgb$runtime_shap_inter)), collapse=","), sep=",")
    row4 <- paste(paste(sim_info, collapse=""),paste(c("XGB-HStat", "predictive", xgb$Hstat,rep(NA,ncol(X)),(xgb$runtime_general_pred+xgb$runtime_Hstat)), collapse=","), sep=",")
    
    return(rbind(row1,row2,row3, row4))
    
    ######################
    # GLMNet
    ###################### 
  } else if (method=="GLMNet"){
    print("GLMNet will be executed")
    glmnet = execute_GLM_XGB(X, trt, y, type_y=y_type, method="GLMNet", approach="all")#,approach="Shap_inter"
    row1 <- paste(paste(sim_info, collapse=""),paste(c("GLMNet-Varimp", "prognostic", glmnet$varimp,rep(NA,ncol(X)),(glmnet$runtime_general_prog+glmnet$runtime_VI)), collapse=","), sep=",")
    row2 <- paste(paste(sim_info, collapse=""),paste(c("GLMNet-Shap", "prognostic", glmnet$shap_mean,rep(NA,ncol(X)),(glmnet$runtime_general_prog+glmnet$runtime_shap_mean)), collapse=","), sep=",")
    
    return(rbind(row1,row2))
  }
  
}
