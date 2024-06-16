################################
# Load libraries and code
################################
library(dplyr)
recent_path <-paste0(getwd(),"/Code")
source(paste0(recent_path,'/Evaluation/Functions_Method_Execution.R'))
source(paste0(recent_path,'/Methods/predMOB_Wrapper.R'))
source(paste0(recent_path,'/Methods/VSURF_Wrapper.R'))
source(paste0(recent_path,'/Methods/PPLasso_Wrapper.R'))
source(paste0(recent_path,'/Methods/XGB_GLMNet_Wrapper.R'))
source(paste0(recent_path,'/Methods/VirtualTwins_Wrapper.R'))
source(paste0(recent_path,'/Methods/Knockoffs_Wrapper.R'))
source(paste0(recent_path,'/Methods/INFO_plus_Wrapper.R'))
source(paste0(recent_path,'/Methods/SIDES_Wrapper.R'))

source(paste0(recent_path,'/Simulations/Simulation_Settings.R'))
# definition of simulation,attributes in file
# n_sim = 100 #100 # count of datasets per simulation scenario
# n = 200 #200 # count of subjects/rows
# n_bm= 20 # count of biomarkers, less than 15 not allowed

# source(paste0(recent_path,'/logging_function.R'))
# logit=NULL

# # define file name of simulated data and of the results file
# subfolder_files = "/Data/"
# file_name_simulations=paste0(getwd(),subfolder_files,"simulated_data.csv")
# file_name_results=paste0(getwd(),subfolder_files,"results.csv")

# # define parameters for results generation
# selected_complexities=c("N","E1","E2","M1","M2","H") # 6 different complexities to choose from
# selected_methods = c("KO_prog","KO_pred","PPLasso","VT","INFO_plus","SIDES","predMOB","VSURF_inter","VSURF_prog","XGB","GLMNet")
# selected_datasets = seq(1,n_sim,1)
# selected_data_input_types = c("continuous","binary", "mixed1", "mixed2")#, "binary", "mixed1", "mixed2")# c("continuous", "binary")
# selected_response_types = c("y_cont","y_binary") #"y_binary",

# read in data
sim_datasets <- read.csv(file_name_simulations, colClasses=c("character", "double", "character",rep("double",n_bm+3)))#for only specific columns colClasses=c("time"="character")
View(sim_datasets)

########################################
# Execution of methods
########################################
# create empty file with header for the results file
header1 <- c("Complexity","Dataset_No","Input_Data_Type","Response_Type", "Method", "Type") #BM_Type instead of Type? xxx
header2 <- paste0("Score_bm",seq(1,n_bm,1))
header3 <- paste0("Selected_bm",seq(1,n_bm,1))
header4 <- c("Run_time")
header= c(header1,header2,header3,header4)
header = matrix(header,ncol=length(header))
if (!file.exists(file_name_results)){ #"sim_data.R"
  # if file not existent write header
  write.table(header, file = file_name_results, sep = ",",
              append = TRUE, quote = FALSE,
              col.names = FALSE, row.names = FALSE)
} 

#selected_methods = c("KO_prog","KO_pred","PPLasso","VT","INFO_plus","SIDES","predMOB","VSURF_inter","VSURF_prog","XGB","GLMNet")
# selected_methods = c("KO_pred") #,"KO_pred"
# selected_methods = c("PPLasso") # OK
# selected_methods = c("KO_prog") # OK
# selected_methods = c("VT") # run wo error
# selected_methods = c("INFO_plus") # run wo error
selected_methods = c("SIDES") # probably OK
selected_methods = c("predMOB") # probably OK
selected_methods = c("VSURF_inter") #OK
selected_methods = c("VSURF_prog") #OK
selected_methods = c("XGB") #OK
selected_methods = c("GLMNet") #OK


# loop through all complexities
for (sel_complexity in 1:length(selected_complexities)){
  # filter data according to selected complexity
  curr_datasets_scen <-sim_datasets[(sim_datasets$Complexity==selected_complexities[sel_complexity]),]
  # for each data type that is selected for execution
  for (sel_type in 1:length(selected_data_input_types)){
    curr_dataset_dt <-curr_datasets_scen[(curr_datasets_scen$Input_Data_Type==selected_data_input_types[sel_type]),]
    for (sel_dataset in 1:length(selected_datasets)){
      # filter data according to dataset (no.)
      curr_dataset_no <-curr_dataset_dt[(curr_dataset_dt$Dataset_No==selected_datasets[sel_dataset]),]
      # in case of mixed biomarkers, the data type has to be set respectively
      if (curr_dataset_no$Input_Data_Type[1]=="mixed1" || curr_dataset_no$Input_Data_Type[1]=="mixed2"){
        mixed_type_binary <- apply(curr_dataset_no,2,function(x) {all(x %in% 0:1)})
        mixed_type_binary_ind <-which(mixed_type_binary)
        names(mixed_type_binary_ind) =NULL
        mixed_type_binary_ind <- mixed_type_binary_ind[! mixed_type_binary_ind %in% c(1, 2, 3)]
        mixed_type_binary_ind <- mixed_type_binary_ind[! mixed_type_binary_ind %in% seq(ncol(curr_dataset_no)-2,ncol(curr_dataset_no),1)]
        curr_dataset_no[,(mixed_type_binary_ind)] <- lapply(curr_dataset_no[,(mixed_type_binary_ind)], factor)
      } else if (curr_dataset_no$Input_Data_Type[1]=="binary"){
        # in case of binary biomarkers, all biomarkers have to be converted to a factor
        curr_dataset_no[,seq(4,(ncol(curr_dataset_no)-3),1)] <- lapply(curr_dataset_no[,seq(4,(ncol(curr_dataset_no)-3),1)], factor)
      }
      # for each response type
      for (sel_resp_type in 1:length(selected_response_types)){
        # conversion of a binary response to factor
        if (selected_response_types[sel_resp_type]== "y_binary"){
          # remove y_cont column
          curr_dataset <- curr_dataset_no %>% dplyr::select(-c("y_cont"))
          # convert type of y_binary to factor
          #NWe
          curr_dataset[,"y_binary"] <- as.factor(curr_dataset[,"y_binary"])
          #curr_dataset <- curr_dataset %>% dplyr::rename(y = y_binary)

        } else {
          # in case of a continuous response, the binary response column can be removed
          # remove y_binary column
          curr_dataset <- curr_dataset_no %>% dplyr::select(-c("y_binary"))
          # NWe
          #curr_dataset <- curr_dataset %>% dplyr::rename(y = y_cont)
        }
        print(curr_dataset)
        rownames(curr_dataset) <- NULL        
        Execute_Methods(curr_dataset, selected_methods, simulation=TRUE, file_name_results=file_name_results) 
        print("after Methods")
      }
    }
  }
}

