################################
# Load libraries and code
################################
library(BinNor)
library(ggplot2)
library(corrgram)
library(corrplot)
library(lessR)
library(Matrix) 
#############################################################
# fixed attributes for all base scenarios
# n = 200 # sample size / count of subjects
# n_sim = 100 # count of simulated datasets
# n_bm = 20 # count of BM
# corr_rel = 0.5 
# corr_non_rel = 0.00 

##################################################
# function for generation the AR(1) correlations 
create_ar1_corrmatrix <- function(n, rho) {
  # n:    dimension of the matrix
  # rho: correlation for the AR(1) setting
  exp <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) -(1:n - 1))
  return(rho^exp)
}

######################################################################
# Function for creation of a correlation matrix
######################################################################
Construct_base_Corr_Matrix <- function(n_bm, complexity="E1",corr_rel=0.5,corr_non_rel=0.0, rho=0.5){
  # n_bm:           count of biomarkers, less than 15 not allowed due to construction of the correlation matrix
  # complexity:     simulation complexity according to the master thesis; possible values: E1 (easy 1), E2 (easy 2), M1 (medium 1), M2 (medium 2), H(hard)  
  # corr_rel:       correlation of biomarkers that are correlated with each other
  # corr_non_rel:   correlation of biomarkers that are correlated with each other
  # rho:            correlation for the AR(1) setting; only relevant to be defined if AR(1) complexity is selected
  
  # if base complexity
  if (complexity=="N" || complexity=="E1" || complexity=="E2" ||complexity=="M1" ||complexity=="M2" ||complexity=="H"){
    # generate empty correlation matrix
    cm = matrix(rep(corr_non_rel,n_bm*n_bm), nrow = n_bm, ncol = n_bm)
    diag(cm) <- 1
    
    # add correlations pertinent to the complexity
    if (complexity=="E2" |complexity=="M1" |complexity=="M2" |complexity=="H"){
      cm[2,c(3,4,5,9)] = corr_rel
      cm[3,c(2,4,5,14)] = corr_rel
      cm[4,c(2,3,7)] = corr_rel
      cm[5,c(2,3)] = corr_rel
      cm[7,c(4,10)] = corr_rel
      cm[8,10] = corr_rel
      cm[9,c(2,10)] = corr_rel
      cm[10,c(7,8,9)] = corr_rel
      cm[12,c(13,14,15)] = corr_rel
      cm[13,c(12,14,15)] = corr_rel
      cm[14,c(3,12,13)] = corr_rel
      cm[15,c(12,13)] = corr_rel
      
    } 
    if (complexity=="M1" | complexity=="M2" | complexity=="H"){
      cm[4,c(5,8)] = corr_rel
      cm[5,4] = corr_rel
      cm[8,4] = corr_rel
      cm[9,15] = corr_rel
      cm[14,15] = corr_rel
      cm[15,c(9,14)] = corr_rel
    }  
    
    if (complexity=="M2"| complexity=="H"){
      cm[4,c(9,14)] = corr_rel
      cm[5,15] = corr_rel
      cm[7, c(8,9)] = corr_rel
      cm[8, c(7,9)] = corr_rel
      cm[9, c(4,7,8)] = corr_rel
      cm[14,4] = corr_rel
      cm[15,5] = corr_rel
      
    }
  } else if (complexity=="AR1"){
    # if AR(1) complexity, create_ar1_corrmatrix constructs the correlation matrix
    cm <- create_ar1_corrmatrix(n=n_bm, rho=rho)
    
  }
  
  # set row and column names (biomarker names)
  colnames(cm) = paste0("bm",seq(1,n_bm,1))
  rownames(cm) = paste0("bm",seq(1,n_bm,1))

  return(cm)
}
######################################################################
# Function for sigma.star generation
######################################################################
Generate_Corr_Matrix <- function(complexity="E1",data_types="continuous",corr_non_rel=0.05,corr_rel=0.5, n_bm=20, rho=0.5){
  # complexity:     simulation complexity according to the master thesis; possible values: E1 (easy 1), E2 (easy 2), M1 (medium 1), M2 (medium 2), H(hard)  
  # data_types:     scale types of the biomarkers for which data should be generated, can be one of the following "continuous", "binary", "mixed1", "mixed2"
  # corr_rel:       correlation of biomarkers that are correlated with each other
  # corr_non_rel:   correlation of biomarkers that are correlated with each other
  # n_bm:           count of biomarkers, less than 15 not allowed due to construction of the correlation matrix
  # rho:            correlation for the AR(1) setting; only relevant to be defined if AR(1) complexity is selected
  #
  # Generate_Corr_Matrix: Generation of sigma.star matrix for the data generation with the BinNor package in accordance with the 
  #                       simulation complexity and other settings
  # Steps:
  # - creates the correlation matrix of the respective scneario
  # - depending on the data types, this correlation matrix is reordered for the needs of compute.sigma.star of the BinNor package
  # - before execution of compute.sigma.star, the nearest positive definite matrix is created (otherwise sigma.star is not provided)
  # - compute.sigma.star finally computes the correlation matrix that can be used for data creation with the BinNor package
  # 
  # Output parameters (among others):
  # - sigma.star: correlation matrix used for data generation with BinNor package
  # - bin_count:  count of binary biomarkers that are created
  # - cont_count: count of continuous biomarkers that are created
  # - n_bm:       count of binary and continuous biomarkers
  # - bm_names:   names of the columns (bm1, bm2, etc.), required because sequence of biomarkers has to be changed for sigma.star creation in mixed type settings
  # - reorder:    for reordering of the biomarkers in case of mixed binary and continuous

  # correlation matrix is constructed according to the input parameters
  cm <-Construct_base_Corr_Matrix(n_bm=n_bm, complexity=complexity,corr_rel=corr_rel,corr_non_rel=corr_non_rel, rho=rho)
  
  # for creation of the data with BinNor, in case of mixed binary and continuous covariates,
  # the binary parameters have to be the first in the correlation matrix (for sigma star generation)
  # thus, the order of correlation matrix for mixed data type scenarios has to be changed
  if (data_types=="mixed1"){
    # Mixed type case 1:
    # Binary: 1,4,6,8,15 (prog/pred BM) 	- 2,7,12 (non-pred/non-prog BM)   --> 8 binary
    # Cont.: 5,9,11,14 (prog/pred BM)	    - 3,10,13 (non-pred/non-prog BM)  --> 7 normal/continuous
    binary_BM_ind =c(1,2,4,6,7,8,12,15)
    cont_BM_ind = c(3,5,9,10,11,13,14)
    # fill until n_bm roughly equally sized with binary and continuous biomarkers
    remain_bm_length <- n_bm-length(binary_BM_ind)-length(cont_BM_ind)
    ratio <- remain_bm_length/2
    # slightly less binary bm than continuous
    binary_BM_ind_remain <- seq(max(cont_BM_ind,binary_BM_ind)+1,max(cont_BM_ind,binary_BM_ind)+floor(ratio),1)
    cont_BM_ind_remain <- seq(max(cont_BM_ind,binary_BM_ind)+floor(ratio)+1,n_bm,1)
    binary_BM_ind <- c(binary_BM_ind, binary_BM_ind_remain)
    cont_BM_ind <- c(cont_BM_ind, cont_BM_ind_remain)
    col_order<- c(binary_BM_ind, cont_BM_ind)
    cm<- corReorder(R=cm, vars=col_order)
    # for reordering of the columns later
    re_col_order <-match(colnames(cm),colnames(cm)[col_order])
    prob_bin <- rep(0.5,length(binary_BM_ind))
    
  } else if (data_types=="mixed2"){
    # Mixed type case 2:
    # Binary: 5,9,11,14 (prog/pred BM)	- 3,10,13 (non-pred/prog BM)  --> 7 binary
    # Cont.: 1,4,6,8,15 (prog/pred BM) 	- 2,7,12 (non-pred/prog BM)   --> 8 normal/continuous
    # roughly the other way around as in mixed 1
    binary_BM_ind = c(3,5,9,10,11,13,14)
    cont_BM_ind =c(1,2,4,6,7,8,12,15)
    # fill until n_bm roughly equally sized with binary and continuous biomarkers
    remain_bm_length <- n_bm-length(binary_BM_ind)-length(cont_BM_ind)
    ratio <- remain_bm_length/2
    # slightly less binary bm than continuous
    binary_BM_ind_remain <- seq(max(cont_BM_ind,binary_BM_ind)+1,max(cont_BM_ind,binary_BM_ind)+floor(ratio),1)
    cont_BM_ind_remain <- seq(max(cont_BM_ind,binary_BM_ind)+floor(ratio)+1,n_bm,1)
    binary_BM_ind <- c(binary_BM_ind, binary_BM_ind_remain)
    cont_BM_ind <- c(cont_BM_ind, cont_BM_ind_remain)
    col_order<- c(binary_BM_ind, cont_BM_ind)
    cm<- corReorder(R=cm, vars=col_order)
    # for reordering of the columns later
    re_col_order <-match(colnames(cm),colnames(cm)[col_order])
    prob_bin <- rep(0.5,length(binary_BM_ind))
  } else if (data_types=="continuous"){
    # for only continuous biomarkers
    binary_BM_ind = NULL
    cont_BM_ind = seq(1,n_bm,1)
    prob_bin = NULL
    re_col_order=NULL
  } else if (data_types=="binary"){
    # for only binary biomarkers
    binary_BM_ind = seq(1,n_bm,1)
    cont_BM_ind = NULL 
    prob_bin <- rep(0.5,length(binary_BM_ind))
    re_col_order=NULL
  }
  # correlation matrix has to stay in this sequence until data creation is finished
  
  # create the nearest positive definite matrix to cm in order to be able to create sigma.star
  if (is.positive.definite(cm)==FALSE){
    cm_adj = nearPD(cm, corr = TRUE, keepDiag = TRUE, ensureSymmetry = TRUE)
    cm_adj = cm_adj$mat
    
  } else {
    cm_adj = cm
  }
  
  # generate sigma.star with the BinNor package as a preparation for data generation
  sigma.star=compute.sigma.star(no.bin=length(binary_BM_ind), no.nor=length(cont_BM_ind), prop.vec.bin=prob_bin,
                                corr.mat=cm_adj)
  
  sigma.star_matrix <- as.matrix(sigma.star$sigma_star) 

  return(list("sigma.star"= sigma.star_matrix,"bin_count"=length(binary_BM_ind), "cont_count"=length(cont_BM_ind), "n_bm"=n_bm, "bm_names"=colnames(cm), "cm_adj"=cm_adj, "cm"=cm, "complexity"=complexity, "reorder"=re_col_order))
  
}


####################################################
# Function for generation of biomarker, treatment
# and response data
####################################################
Create_Base_Scenario_Data <- function(n_subj, sigma_star, bin_count, cont_count, bm_names, mean=0, var=1, prop_bin=0.5, complexity="N", reorder=NULL, trt_eff =0.5, prog_eff_size = 0.5, pred_eff_size = 0.5, inter_act_eff_size=0.3, eff_type="linear"){
  # n_subj:             sample size / count of subjects (rows)
  # sigma_star:         correlation matrix used for data generation
  # bin_count:          count of binary biomarkers
  # cont_count:         count of continuous biomarkers
  # bm_names:           names of the biomarkers
  # mean:               mean of the normal distributed biomarkers
  # var=1:              variance of the normal distributed biomarkers
  # prop_bin:           probability p for binary biomarkers (related to Bernoulli distribution)
  # complexity:         simulation complexity according to the master thesis; possible values: 
  #                     E1 (easy 1), E2 (easy 2), M1 (medium 1), M2 (medium 2), H(hard)  
  # reorder:            for reordering of the biomarkers in case of mixed binary and continuous
  # trt_eff:            strength of the treatment effect
  # prog_eff_size:      strength of the prognostic biomarkers
  # pred_eff_size:      strength of the predictive biomarkers
  # inter_act_eff_size: strength of the interactions effects between biomarkers (only relevant with H complexity)
  # eff_type:           either "linear" or "square" depending on which effect type should be used 
  #                     for predictive, prognostic and interaction effects
  #
  # Steps:
  # - error handling
  # - generation of biomarker data with regard to sigma.star and simulation setting (e.g. in terms of the count of binary and cont. biomarkers)
  # - creation of treatment variable
  # - creation of outcome variable in accordance with the outcome type
  #
  # Output parameters
  # - bm:   created biomarker data according to sigma.star
  # - trt:  generated treatment variable
  # - y:    generated outcome variable
  
  # determination of bm_names
  bm_names <- colnames(sigma_star)

  # depending on the complexity, the true predictive and prognostic biomarkers can be determined
  if (!complexity=="N"){
    # in case of the null complexity, no biomarker is predictive nor prognostic
    # otherwise:
    prog_bm_ind <- c(1,4,5,8)
    pred_bm_ind <- c(9,11,14,15)
    if (complexity=="M2" | complexity =="H"){
      # for complexity M2 and H, we have to include BM6 which is prognostic and predictive
      prog_bm_ind <- sort(c(prog_bm_ind, 6))
      pred_bm_ind <- sort(c(pred_bm_ind, 6))
    }
  } else {
    # for the null sceanario
    prog_bm_ind = NULL
    pred_bm_ind =NULL
  }
  
  if (cont_count==0){
    mean=NULL
    var=NULL
  } else {
    mean <- rep(mean,cont_count)
    var <-  rep(var,cont_count)
  }
  if (bin_count==0){
    prob_bin = NULL
  } else {
    prob_bin <- rep(prop_bin,bin_count) # !! is also one value OK instead of a vector?
  }
  
  # generate biomarker data with specified parameters
  data=jointly.generate.binary.normal(n_subj,no.bin=bin_count,no.nor=cont_count,prop.vec.bin=prob_bin,
                                      mean.vec.nor=mean,var.nor=var, sigma_star=sigma_star,
                                      continue.with.warning=TRUE)
  colnames(data) <- bm_names

  
  if (bin_count>0 && cont_count>0) {
    # in case of a mixed setting, reordering of the columns is required
    data <- data[,reorder] #order(colnames(data))
  } 
  
  # generate binary treatment variable
  bin_n_trt=1
  bin_p_trt=0.5
  trt = rbinom(n=n_subj, size=bin_n_trt, prob=bin_p_trt)
  
  pred_effect=0
  prog_effect=0
  
  # predictive and prognostic effects are generated based on eff_type and the prognostic and predictive indices
  if (!complexity=="N"){
    if (eff_type=="linear"){
      for (i in 1:length(prog_bm_ind)){
        prog_effect= prog_effect + data[,prog_bm_ind[i]]
      }
      for (i in 1:length(pred_bm_ind)){
        pred_effect= pred_effect + data[,pred_bm_ind[i]]*trt
        
      }
    } else if (eff_type=="square"){
      for (i in 1:(length(prog_bm_ind))){
        # only the first half with squared effects (rounded down), for the rest with linear effects         
        if (i<= (floor(length(prog_bm_ind)/2))){
          prog_effect= prog_effect + data[,prog_bm_ind[i]]*data[,prog_bm_ind[i]]
        } else {
          prog_effect= prog_effect +  data[,prog_bm_ind[i]]
        }
      }
      for (i in 1:(length(pred_bm_ind))){
        # only the first half with squared effects (rounded down), for the rest with linear effects       
        if (i<= (floor(length(pred_bm_ind)/2))){
          pred_effect= pred_effect + data[,pred_bm_ind[i]]*data[,pred_bm_ind[i]]*trt
        } else {
          pred_effect= pred_effect + data[,pred_bm_ind[i]]*trt
        }
      }
    }
  }

  # for the hard complexity, interaction terms have to be added for the response
  # BM4*BM5
  # BM8*BM9
  # BM14*BM15
  # depending on eff_type (linear or square effects)
  inter_eff =0 
  if (complexity =="H"){
    inter_bm_ind <- c(4,5,8,9,14,15) # always the next two of the array...
    if (eff_type=="linear"){
      for (i in 1:(length(inter_bm_ind)/2)){
        inter_eff = inter_eff + inter_act_eff_size*data[,inter_bm_ind[(2*i)-1]]*data[,inter_bm_ind[(2*i)]]
      }
    } else if (eff_type=="square"){
      # only the first half with squared effects (rounded down), for the rest with linear effects         
      for (i in 1:(length(inter_bm_ind)/2)){
        if (i<= (floor(length(inter_bm_ind)/4))){
          inter_eff = inter_eff + inter_act_eff_size*data[,inter_bm_ind[(2*i)-1]]*data[,inter_bm_ind[(2*i)-1]]*data[,inter_bm_ind[(2*i)]]*data[,inter_bm_ind[(2*i)]]
        } else {
          inter_eff = inter_eff + inter_act_eff_size*data[,inter_bm_ind[(2*i)]]*data[,inter_bm_ind[(2*i)]]
        }
      }
    }
  }
  
  theta_pred = pred_eff_size
  theta_prog = prog_eff_size
  
  # for binary response
  logit_p_y <- theta_prog * prog_effect  * pred_eff_size * pred_effect +trt_eff*trt + inter_eff
  p_y <- 1/(1+exp(-logit_p_y))
  y_bin <-  rbinom(n_subj,1,p_y)

  # for continuous response
  mu_overall = pred_eff_size*pred_effect+prog_eff_size*prog_effect+ trt_eff*trt + inter_eff 

  # mu and logit_p_y is the same
  sigma_y2 =0.5 
  y_cont = rnorm(n=n_subj, mean = mu_overall, sd=sqrt(sigma_y2))
  
  return(list("bm"=data, "trt"=trt, "y_cont"=y_cont, "y_bin"=y_bin))
  
}



