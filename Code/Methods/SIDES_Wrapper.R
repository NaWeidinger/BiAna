##################################
##################################
# Code for SIDES
##################################
##################################

################################
# Load libraries and code
################################
library(SIDES)
recent_path <-paste0(getwd(),"/Code")
source(paste0(recent_path,'/Methods/Help_Functions.R'))

#############################################
# Wrapper Function for SIDES
#############################################

execute_SIDES <- function(X, trt, y, type_X, type_y, D=0, L=3, S=20, M=5, gamma=c(1,1,1), H=2, num_crit=1, alpha=0.10, nsim=100, ord.bin=3, upper_best=TRUE, prop_gpe=c(0.7,0.3)){
  # X:      biomarker data, of type data frame with biomarkers represented as columns; binary biomarker are of type factor, continuous of type numeric
  # trt:    treatment variable as a vector, of type integer
  # y:      response as a vector, of type factor if binary, of type numeric if continuous
  # type_X: vector providing the scale types of the biomarkers in X; allowed values for this wrapper functions are "continuous", "ordinal" or "nominal"
  # type_y: scale type of the response; allowed for this wrapper function are "continuous", "binary"
  # all other parameters as detailed in the master thesis and as required for the original function of SIDES
  
  # preparation of biomarkers and response as required for original function of SIDES
  data_prep <- convert_variables(y=y, X=X) 
  y=data_prep$y
  X=data_prep$X
  data=cbind(y,trt,X) 
  
  # capture start time 
  start = Sys.time()
  
  # compute results, call original function
  sides = SIDES(all_set=data,
                type_var=type_X, type_outcome=type_y,
                level_control=0, D=D, L=L, S=S, M=M, gamma=gamma, H=H, num_crit=num_crit,
                alpha=alpha, nsim=nsim, ord.bin=ord.bin, upper_best=upper_best, prop_gpe=prop_gpe)   #, seed=42  
  
  # capture end time and calculate run time
  end = Sys.time()
  runtime = as.numeric(difftime(end, start, units='mins')) # only minutes are returned
  
  confirmed = sides$confirmed 

  # preparation of results
  if (length(confirmed[[1]])>0){
    # if biomarkers are identified
    confirmed_all = c()
    for (i in 1:length(confirmed[[1]])){
      confirmed_all <- c(confirmed_all, confirmed[[1]][[i]][[1]])
    }
    tab <- as.data.frame(table(confirmed_all))
    scores <- rep(0, ncol(data))
    as.integer(as.character(tab$confirmed_all))
    scores[as.integer(as.character(tab$confirmed_all))] <-tab$Freq
    names(scores) <- colnames(data)
    scores=scores[!names(scores)== 'trt']
    scores=scores[!names(scores)== 'y']
    # does not work...
    selected <- rep(0, ncol(data))
    selected[as.integer(as.character(tab$confirmed_all))] <- 1
    names(selected) <- colnames(data)
    selected=selected[!names(selected)== 'trt']
    selected=selected[!names(selected)== 'y']
  } else {
    # if no biomarkers are identified
    scores <- rep(0, ncol(X))
    selected <- rep(0, ncol(X))
    names(scores) <- colnames(X)
    names(selected) <- colnames(X)
  }

  return(list( "scores"=scores, "selected"=selected, "runtime"=runtime))
}
