################################
# Load libraries and code
################################
library(dplyr) 
recent_path <-paste0(getwd(),"/Code/Simulations/")
source(paste0(recent_path,'Simulation_Scenarios.R'))

#############################
# parameters for simulations
#############################
# read in from simulation settings file:
source(paste0(recent_path,'Simulation_Settings.R'))

# subfolder_files = paste0(getwd(),"/Data/")
# file_name_simulations=paste0(subfolder_files,"simulated_data.csv")

# default parameter values are already included in the function Create_Datasets
# definition of simulation scenarios
# n_sim = 100 # count of datasets per simulation scenario
# n = 200 # count of subjects/rows
# n_bm= 20 # count of biomarkers (covariates); less than 15 not allowed
# corr_non_rel=0.00
# corr_rel=0.5
# data_input_types = c("continuous", "binary", "mixed1", "mixed2")
# all_complexities=c("N","E1","E2","M1","M2","H")
# eff_type="linear"
# mean=0 # mean of the normal distributed biomarkers
# var=1 # variance of the normal distributed biomarkers
# prop_bin=0.5 # probability p for binary biomarkers (related to Bernoulli distribution)
#############################################################

Create_Datasets <- function(complexities=c("N","E1","E2","M1","M2","H"), n_sim=100, n=200, n_bm=20, data_input_types=c("continuous", "binary", "mixed1", "mixed2"),corr_non_rel=0.00,corr_rel=0.5, seed=1, rho=0.5, eff_type="linear"){
  # n_sim:            count of datasets per simulation scenario
  # n:                sample size / count of subjects (rows)
  # n_bm:             count of biomarkers (covariates); less than 15 not allowed
  # corr_non_rel:     correlation of the non-related biomarkers 0.00
  # corr_rel:         correlation value of biomarkers that do not have a (true) correlation
  # data_input_types: scale types of the biomarkers for which data should be generated, subset of c("continuous", "binary", "mixed1", "mixed2")
  # complexities:     complexities for which data should be generated, defines the correlation structure between the biomarkers, different complexity levels are available: c("N","E1","E2","M1","M2","H")
  # seed:             seed for data generation for reproducability
  # rho:              correlation for the AR(1) setting
  # eff_type:         either "linear" or "square" depending on which effect type should be used for predictive, prognostic and interaction effects
  
  # create empty file with header
  header1 <- c("Complexity","Dataset_No","Input_Data_Type")
  header2 <- paste0("bm",seq(1,n_bm,1)) 
  header3 <- c("trt", "y_binary","y_cont")
  header= c(header1,header2,header3)
  header = matrix(header,ncol=length(header))
  print(file_name_simulations)
  write.table(header, file = file_name_simulations, sep = ",",
              append = FALSE, quote = FALSE,
              col.names = FALSE, row.names = FALSE)
  
  # set seed for reproducability
  set.seed(seed)  
  
  # for each simulation complexity
  for (complexity in 1:length(complexities)){
    # for each scale type of the biomarkers
    for (data_input_type in 1:length(data_input_types)){
      cm <- Generate_Corr_Matrix(complexity=complexities[complexity],data_types=data_input_types[data_input_type],corr_non_rel=corr_non_rel,corr_rel=corr_rel, n_bm=n_bm, rho=rho)
      # n_sim datasets should be generated
      for (simulation_no in 1:n_sim){
        # create dataset with respective complexity, dataset no. and  correlation matrix 
        # length(all_complexities) * n_sim * n * length(data_input_types) = count of rows of created file (without header)
        # --> 6 * 100 * 200 * 4 = 480.000 rows + Header 
        # generate a dataset with the specified parameters
        gen_data <- Create_Base_Scenario_Data(n_subj=n, sigma_star=cm$sigma.star, bin_count=cm$bin_count, cont_count=cm$cont_count, bm_names=cm$bm_names, mean=0, var=1, prop_bin=0.5,complexity = cm$complexity , reorder=cm$reorder,  eff_type=eff_type)
        # save respective data frame to simulation file
        df <- data.frame(complexities[complexity], simulation_no, data_input_types[data_input_type],gen_data$bm, gen_data$trt, gen_data$y_bin, gen_data$y_cont)
        print(file_name_simulations)
        write.table(df, file = file_name_simulations, sep = ",",
                    append = TRUE, quote = FALSE,
                    col.names = FALSE, row.names = FALSE)
      }
    }
  }
}

########################################
# Creation of datasets 
########################################
set.seed(123)
if (!file.exists(file_name_simulations)){ 
  # create datasets for simulation settings
  #Create_Datasets(complexities=all_complexities, n_sim=n_sim, n=n, n_bm=n_bm, data_input_types=data_input_types, corr_non_rel=corr_non_rel,corr_rel=corr_rel, seed = 1)
  #Create_Datasets() # this would be sufficient to create all datasets of all scenarios
  Create_Datasets(complexities=c("E1"), n_sim=n_sim, n=n, n_bm=n_bm, data_input_types=c("continuous", "binary", "mixed1", "mixed2"),corr_non_rel=0.00,corr_rel=0.5, seed=1, rho=0.5, eff_type="linear")
} else {
  print("Simulated data already exists.")
}



