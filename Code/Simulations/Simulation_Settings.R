# set parameters for simulation settings

n_sim=100 #count of datasets per simulation scenario
n=200    #sample size / count of subjects (rows)
n_bm=20 #count of biomarkers (covariates); less than 15 not allowed

selected_complexities=c("E1") #c("N","E1","E2","M1","M2","H") # 6 different scenarios to choose from
selected_methods = c("KO_prog","KO_pred","PPLasso","VT","INFO_plus","SIDES","predMOB","VSURF_inter","VSURF_prog","XGB","GLMNet")
selected_datasets = seq(1,n_sim,1)
selected_data_input_types = c("continuous","binary", "mixed1", "mixed2")#, "binary", "mixed1", "mixed2")# c("continuous", "binary")
selected_response_types = c("y_cont","y_binary") #"y_binary",

subfolder_files = "/Data/"
file_name_simulations=paste0(getwd(),subfolder_files,"simulated_data.csv")
file_name_results=paste0(getwd(),subfolder_files,"results.csv")
