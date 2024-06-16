###########################################
# Load libraries and code
###########################################
library(reshape2)  
library(data.table)
library(dplyr) 
library(ggpubr)
library(readr)
library(dplyr)
library(hablar)
library(stringr)
library(ggplot2)
library(ggtext)
library(scales)

#######################################
# Function for combining result files
#######################################
# function for read in of files containing results
combine_results <- function(subfolder="Results", prefix){
  # subfolder:  subfolder name for loading the results
  # prefix:     prefix of the files; only files with these prefix are loaded and combined
  pattern <- paste0("^", prefix, ".*\\.csv$")
  files <- list.files(subfolder, pattern=pattern, full.names = TRUE)
  i=1
  for (f in files){
    # for (f in files[-1]){
    if (i==1){
      DF <-  read.csv(files[1])
      i=i+1
    } else {
      df <- read.csv(f)      # read the file
      DF <- rbind(DF, df)    # append the current file
    }
  }
  return(DF)
}

############################################
# FunctionScaling of the scores 
# by preserving the 0
############################################
Score_scaling <- function(data, range_max=100){
  # data:         as provided by visualize_scores_selected
  # range_max:    maximum value which should be present after scaling; the maximum absolute value in the data
  #               is scaled to this given value
  df_temp = data 
  if (max(abs(df_temp))!= 0){
    df_temp = df_temp/max(abs(df_temp))*range_max
  }
  # return value is the scaled data frame
  return(df_temp)
}

############################################################
# Visualisation of scores and selected attribute
############################################################

visualize_scores_selected <- function(data, vis_type="BM", type_BM="predictive", type_eval="Scores", complexity="E1",Input_Data_Type="continuous", Response_Type="continuous", BM_ind=NULL, scaled=TRUE, add_info=NULL){
  # data:             evaluation data frame with the results from the simulation settings
  # vis_type:         "BM" or "Method"; indicatior if the visualisations should be grouped by BMs or by the methods
  # type_BM:          "predictive" or "prognostic" to be shown in the visualization
  # type_eval:        "Score" or "Selected" to be displayed in the visualization
  # complexity:         complexity of the simulation setting for which the vizualization should be generated
  # Input_Data_Type:  scale type of the biomarkers, either "continuous","binary", "mixed1" or "mixed2"
  # Response_Type:    scale type of the response, either "continuous" or "binary"
  # BM_ind:           has to be provided if only results of specific biomarkers should be included in the visualization
  #                   especially useful for datasets with many biomarkers
  # scaled:           true or false, only used for vis_type = "Method", for vis_type = "BM"
  #                   scores are always scaled since they are otherwise not comparable
  # add_info:         additional information that should be added to the title of the visualization

  info_columns = c(1,2,3,4,5,6)
  n_bm = (ncol(data)-length(info_columns)-1)/2 # since scores and selected are present, -1 for run_time

  # exclude columns, not required for visualization
  if (type_eval=="Selected"){
    df <- data %>% dplyr::select(-contains("Score"))
    # rename columns in order to show only the BM names
    colnames(df)<-gsub("Selected_","",colnames(df))
  } else {#if (type_eval=="Scores"){
    df <- data %>% dplyr::select(-contains("Selected"))
    # rename columns in order to show only the BM names
    colnames(df)<-gsub("Score_","",colnames(df))
  }
  # remove column run_time for this evaluation
  df <- df %>% dplyr::select(-contains("Run_time"))
  # remove NAs if present
  df <- na.omit(df)

  # filter data according to complexity, biomarker scale type and response scale type
  df <-df[(df$Complexity ==complexity),]
  df <-df[(df$Input_Data_Type ==Input_Data_Type),]
  df <-df[(df$Response_Type ==Response_Type),]

  print(head(df))

  # remove predictive or prognostic methods depending on what should be visualized and set titles of the visualization
  if (type_BM == "predictive"){
    df <-df[(df$Type=="predictive"),]
    lab_BM = "Predictive Biomarkers"
  } else {
    df <-df[(df$Type=="prognostic"),]
    lab_BM = "Prognostic Biomarkers"
  }
  lab_BM = paste0(lab_BM, " - ",  type_eval, " (grouped by ", vis_type, ")")

  # filter BMs according to BM_ind (if provided)
  if (!is.null(BM_ind)){
    relevant_cols = c(info_columns, BM_ind+length(info_columns))
    df = df[, relevant_cols]
  }

  # in case of scaling or BM visualisation, the scores have to be scaled  
  if (type_eval=="Scores" && (vis_type=="BM" || scaled==TRUE)){
    df_scaled = apply(X=df[-info_columns], MARGIN=1, FUN=Score_scaling, range_max=100)
    df[-info_columns] = t(df_scaled)
  } 
  
  # prepare data for visualization
  df_melt <- as.data.frame(melt(data.table(df), info_columns))
  
  # shorten labels
  df_melt["Method"][df_melt["Method"] == "predictive-Knockoffs (LI)"] <- "pred. KO (LI)"
  df_melt["Method"][df_melt["Method"] == "predictive-Knockoffs (CF)"] <- "pred. KO (CF)"
  df_melt["Method"][df_melt["Method"] == "VSURF-interaction"] <- "VSURF-inter"

  #assuming that all methods are run on the same count of datasets, the number of datasets is determined (for the selection rate)
  no_datasets <- length(unique(df_melt$Dataset_No))
  
  # for the selected attribute, the data has to be aggregated (sum of selected attribute)
  df_adj <- setNames(aggregate(df_melt$value, by=list(df_melt$Method, df_melt$Type, df_melt$variable), FUN=sum),  # setNames function
                     c("Method", "Type", "variable", "count"))      # df_melt_filtered = df_melt[(df_melt$value==1),]

  # define coloring for the biomarkers, depending on the truth in the dataset
  man_color <- rep("coral1", n_bm)
  # for E1 and E2, M1, M2, H: 1,4,5,8 prognostic, 9,11,14,15 predictive
  # for M1 and H additional BM6 is prognostic and predictive
  if (type_BM=="predictive"){
    if (complexity != "N"){
      man_color[c(9,11,14,15)] <- "palegreen3"
    }
    if (complexity == "M2" || complexity =="H"){
      man_color[6] <- "palegreen3"
    }
  } else {
    if (complexity != "N"){
      man_color[c(1,4,5,8)] <- "palegreen3"
    }
    if (complexity == "M2" || complexity =="H"){
      man_color[6] <- "palegreen3"
    }
  } 
  if (!is.null(BM_ind)){
    # in case only specific BMs should be displayed, the colors need to be restricted to these
    man_color = man_color[BM_ind]
  }

  lab_sub = paste0("(Covariate Type: ",Input_Data_Type,", Response Type: ", Response_Type, ", Complexity: ",complexity ,  ")")
  
  # add additional information to the visualization titles, if provided
  if (!is.null(add_info)){
    lab_sub = paste0(lab_sub,", ", add_info,")")
  } else {
    lab_sub = paste0(lab_sub,")")
  }
    
  if (type_eval=="Scores"){
    # in case of scores should be visualized
    if (vis_type=="BM"){
      p <- ggplot(df_melt, aes(x=Method, y=value, fill=variable)) + # fill=Method
        geom_boxplot()+ ylim(0,100) + theme_minimal() +geom_vline(xintercept=0)+
        facet_wrap(~variable) + xlab("Method") + ylab("Scores") +  theme_minimal()+
        labs(title = lab_BM, subtitle = lab_sub)+ 
        coord_flip() +  scale_fill_manual(values=man_color) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + 
        theme(legend.position = "none")
    } else {
      p <- ggplot(df_melt, aes(x=variable, y=value, fill=variable)) + #fill=Method
        geom_boxplot() +
        facet_wrap(~Method, scale="free", nrow = 2) + xlab("Biomarker") + ylab("Scores")+ 
        labs(title = lab_BM, subtitle = lab_sub)+ 
        coord_flip()+  scale_fill_manual(values=man_color) + theme_minimal() +
        theme(axis.text.x = element_text(angle = 90)) +
        theme(axis.text.y = ggtext::element_markdown(colour = man_color)) + 
        theme(legend.position = "none")
    }
  } else { 
    # in case of the selected attribute should be visualized
    # sum over all selected attribute, per complexity, method, biomarker type (predictive prognostic), biomarker scale type, response scale type
      df_adj$percentage <- df_adj$count/no_datasets*100
      if (vis_type=="BM"){ 
        p <- ggplot(df_adj, aes(fill=variable,x=Method, y=percentage)) +geom_bar(position="dodge",stat="identity") +
          facet_wrap(~variable) +
          labs(title = lab_BM, subtitle = lab_sub)+ 
          ylab("Selection rate [%]") +  
          coord_flip() + theme_minimal() +
          scale_fill_manual(values=man_color) + 
          theme(axis.text.x = element_text(angle = 90)) +theme(legend.position = "none") 
      } else {
        p <-ggplot(df_adj, aes(fill=variable,x=variable, y=percentage)) +geom_bar(position="dodge",stat="identity") +
          facet_wrap(~Method) + 
          labs(title = lab_BM, subtitle = lab_sub)+ 
          ylab("Selection rate [%]") + xlab("Biomarker") + theme(legend.position = "none") + 
          coord_flip() + theme_minimal() +
          scale_fill_manual(values=man_color)  + 
          theme(axis.text.x = element_text(angle = 90)) +
          theme(axis.text.y = ggtext::element_markdown(colour = man_color)) + theme(legend.position = "none")
      }
  }
  # return the created plot
  return(p)
}


################################################
# (Helper) Function for Calculation of Metrics:
# TPR, TNR, PPV, NPV
################################################

Metrics_calculation <- function(x, true_ind, selected=TRUE){
  # x:          row of data frame provided by Calculate_Visualize_Metrics
  # true_ind:   indices of true prognostic/predictive biomarkers 
  # selected:   TRUE if calculation should be done for the selected attribute
  
  # count of true prognostic/predictive biomarkers
  count_true_rel = length(true_ind)

  
  if (selected==TRUE){
    # for the seelected attribute, the prognostic/predictive biomarker indices have to be determined
    sel_ind <- which(x==1)
  } else { 
    # in case of scores: the q highest scores with q being the true count of predictive/prognostic biomarkers in the dataset
    max_values <- x[order(x, decreasing = TRUE)[1:count_true_rel]] # get q max values from the scores
    if (all(max_values==0) || var(x) == 0){ # in case of all scores =0 or all same values
      sel_ind=NULL
    } else {
      sel_ind <- which(x %in% max_values) 
    }
  }
  names(sel_ind)=NULL

  # metrics can be calculated accordingly
  bm_ind = seq(1,length(x), 1) 
  
  TP = intersect(true_ind, sel_ind)
  FN = setdiff(true_ind, sel_ind)
  TN = intersect(setdiff(bm_ind, true_ind), setdiff(bm_ind, sel_ind))
  FP = setdiff(sel_ind, true_ind)
  TP_count=length(TP)
  FN_count=length(FN)
  TN_count=length(TN)
  FP_count=length(FP)
  
  PPV=TP_count/(FP_count+TP_count)
  TPR=TP_count/(FN_count+TP_count) # sensitivity, recall
  NPV=TN_count/(TN_count+FN_count)
  TNR=TN_count/(TN_count+FP_count) # specificity

  # return values are PPV, NPV, TNR, TPR
  c(PPV,NPV,TNR,TPR)
  
}

################################################
# Function for Visualization of Metrics:
# TPR, TNR, PPV, NPV
################################################

Calculate_Visualize_Metrics <- function(data, type_BM="predictive", complexity="E1", type_eval="Scores", Input_Data_Type="continuous", Response_Type="continuous", method="all",add_info=NULL){
  # data:             evaluation data frame with the results from the simulation settings
  # type_BM:          "predictive" or "prognostic" to be shown in the visualization
  # complexity:         complexity of the simulation setting for which the vizualization should be generated
  # type_eval:        "Score" or "Selected" to be displayed in the visualization
  # method:           default "all" to show results for all methods, but can be used to filter results 
  # Input_Data_Type:  scale type of the biomarkers, either "continuous","binary", "mixed1" or "mixed2"
  # Response_Type:    scale type of the response, either "continuous" or "binary"
  # add_info:         additional information that should be added to the title of the visualization
  
  info_columns = c(1,2,3,4,5,6)
  n_bm = (ncol(data)-length(info_columns)-1)/2 # since scores and selected are present, -1 for run_time
  
  # shorten labels
  data["Method"][data["Method"] == "predictive-Knockoffs (LI)"] <- "pred. KO (LI)"
  data["Method"][data["Method"] == "predictive-Knockoffs (CF)"] <- "pred. KO (CF)"
  data["Method"][data["Method"] == "VSURF-interaction"] <- "VSURF-inter"
  
  # determine true prognostic / predictive biomarkers based on the complexity
  true_pred = rep(0, n_bm)
  true_prog = rep(0, n_bm)
  if (type_BM=="predictive" ){ 
    if (complexity != "N"){
      true_pred[c(9,11,14,15)] <- 1
    }
    if (complexity == "M2" || complexity =="H"){
      true_pred[6] <- 1
    }
  }
  # if prognostic
  if (type_BM=="prognostic" ){ 
    if (complexity != "N"){
      true_prog[c(1,4,5,8)] <- 1
    }
    if (complexity == "M2" || complexity =="H"){
      true_prog[6] <- 1
    }
  } 

  true_prog_ind <- which(true_prog==1)
  true_pred_ind <- which(true_pred==1)
  
  # exclude columns, not required for visualization
  if (type_eval=="Selected"){
    df <- data %>% dplyr::select(-contains("Score"))
    # rename columns in order to show only the BM names
    colnames(df)<-gsub("Selected_","",colnames(df))
  } else {#if (type_eval=="Scores"){
    df <- data %>% dplyr::select(-contains("Selected"))
    # rename columns in order to show only the BM names
    colnames(df)<-gsub("Score_","",colnames(df))
  }
  # remove column run_time for this evaluation
  df <- df %>% dplyr::select(-contains("Run_time"))
  # remove NAs if present
  df <- na.omit(df)
  
  # filter data according to complexity, biomarker scale type and response scale type
  if (complexity != "all"){
    df <-df[(df$Complexity ==complexity),]
  }
    df <-df[(df$Input_Data_Type ==Input_Data_Type),]
  df <-df[(df$Response_Type ==Response_Type),]
  
  # remove predictive or prognostic methods depending on what should be visualized and set titles of the visualization
  if (type_BM == "predictive"){
    df <-df[(df$Type=="predictive"),]
    lab_BM = "Predictive Biomarkers"
  } else {# prognostic
    df <-df[(df$Type=="prognostic"),]
    lab_BM = "Prognostic Biomarkers"
  }
  lab_BM = paste0(lab_BM , " - ", type_eval)
  
  # calculate metrics
  df_short <- df[,-info_columns]
  if (type_BM=="predictive" ){ 
    if (type_eval=="Selected"){
      metrics <- apply(X=df_short,MARGIN=1,FUN=Metrics_calculation, true_ind=true_pred_ind, selected=TRUE)
      # metrics <- FDR_FNR_calculation(x=df_short, true_ind=true_pred, selected=TRUE)
    } else {
      metrics <- apply(X=df_short,MARGIN=1,FUN=Metrics_calculation, true_ind=true_pred_ind, selected=FALSE) 
    }
  
  } else { # prognostic
    if (type_eval=="Selected"){
      metrics <- apply(X=df_short,MARGIN=1,FUN=Metrics_calculation, true_ind=true_prog_ind, selected=TRUE) 
    } else {
      metrics <- apply(X=df_short,MARGIN=1,FUN=Metrics_calculation, true_ind=true_prog_ind, selected=FALSE) 
    }
  }
  df$PPV <- metrics[1,]*100 # PPV
  df$NPV <- metrics[2,]*100 # NPV
  df$TNR <- metrics[3,]*100 # TNR
  df$TPR <- metrics[4,]*100 # TPR
  
  # prepare data frame for visualisation
  melt_cols <-  seq(1,(length(info_columns)+n_bm),1) #c(info_columns, ((length(info_columns)+1):(length(info_columns)+n_bm)))
  df_melt <- as.data.frame(melt(data.table(df), c(melt_cols)))
  
  # add additional information to the visualization titles, if provided
  lab_sub = paste0("(Complexity: ", complexity, ", Response Type: ",Response_Type, ", Covariate Type: ",Input_Data_Type) #, ")"
  if (!is.null(add_info)){
    lab_sub = paste0(lab_sub,", ", add_info,")")
  } else {
    lab_sub = paste0(lab_sub,")")
  }
  
  # visualization of the results
  if (complexity !="all"){
    p <- ggplot(df_melt, aes(x=Method, y=value, fill=Method)) + 
      geom_boxplot()+ geom_jitter(color="black", size=0.2, alpha=0.3)+ ylim(0,100)+
      facet_wrap(~variable) + 
      xlab("Method") + ylab("Metric value") +
      labs(title = lab_BM, subtitle = lab_sub)+ 
      coord_flip() +  theme_minimal() + theme(legend.position = "none")
  } else {
    colnames(df_melt)[colnames(df_melt) == "variable"] ="Metric"
    lab_BM = paste0(lab_BM , " - ", method)
    df_melt <-df_melt[(df_melt$Method == method),]
    p <- ggplot(df_melt, aes(x=Metric, y=value, fill=Metric)) + 
      geom_boxplot()+ 
      geom_jitter(color="black", size=0.4, alpha=0.3) + ylim(0,100) +
      facet_wrap(~factor(Complexity, levels=c("N","E1","E2","M1","M2","H"))) +
      xlab("Metric") + ylab("Metric value") +
      labs(title = lab_BM, subtitle = lab_sub)+ 
      coord_flip() +  theme_minimal() +theme(legend.position = "none")
  }
  
  # return the created plot
  return(p)
}

dataset <- read.csv(file_name_results)
visualize_scores_selected(data=dataset, vis_type="Method", type_BM="predictive", type_eval="Scores", complexity="E1",Input_Data_Type="continuous", Response_Type="continuous", BM_ind=NULL, scaled=TRUE, add_info=NULL)
  

                    