######################################
######################################
# Virtual Twins
######################################
######################################

################################
# Load libraries and code
################################
recent_path <-paste0(getwd(),"/Code")
source(paste0(recent_path,'/Methods/Help_Functions.R'))
library(varImp)
library(randomForest)
library(aVirtualTwins)
library(caret)

#####################################
# Wrapper Function for Virtual Twins
#####################################

VT_Wrapper <- function(X, trt, y){ 
  # X:    biomarker data, of type data frame with biomarkers represented as columns; binary biomarker are of type factor, continuous of type numeric
  # trt:  treatment variable as a vector, of type integer
  # y:    response as a vector, of type factor if binary, of type numeric if continuous
 
  # error in case provided data is not supported by the original function,
  error = check_data_types(X=X, y=y, method="VT")

  if (error==TRUE){
    return(NULL) 
  }
  
  # preparation of dataset
  dataset <- data.frame(X, trt, y)
  vt.o <- vt.data(dataset=dataset, outcome.field="y", treatment.field="trt", interactions=TRUE) 
  
  # capture start time
  start = Sys.time()
  
  model.rf <- randomForest(x = vt.o$getX(interactions = T),
                           y = vt.o$getY(),
                           ntree = 500)
  vt.f.rf <- vt.forest("one", vt.data = vt.o, model = model.rf, interactions = T)
  
  # prepare data for random forest application for predicting zi with xi (xi, zi)
  # with zi = vt.f.rf$difft
  zi = vt.f.rf$difft
  xi <- X 
  new_data <- cbind(xi, zi)
  
  # compute Variable Importance 
  rf_fin <- randomForest(zi ~ ., data= new_data, importance=TRUE)
  varimp <- caret::varImp(rf_fin) # , conditional=TRUE does not work with importance = TRUE?
  # "1" will be the favorable outcome 
  
  # capture end time and calculate run time
  end = Sys.time()
  runtime = as.numeric(difftime(end, start, units='mins')) # only minutes are returned
  
  # preparation of vector
  names <- rownames(varimp)
  varimp <- varimp$Overall
  names(varimp) <- names

  return(list("scores"=varimp, "runtime"=runtime))

}
