#' The two predictive knockoff filters we presented in Section 2.3
#'
#'
#' @param X data.frame (or tibble) the input X, nrow(X) are the number of examples and ncol(X) the number of variables (features)
#' @param X_tilde data.frame (or tibble) the knockoffs of X.
#' @param y response vector with \code{length(y) = nrow(X)}. Accepts "numeric" (family="gaussian") or binary "factor" (family="binomial").
#' @param t treatment allocation binary vector with \code{length(t) = nrow(X)}. 
#' @param fdr_nominal target false discovery rate. Can be a vector of multiple thresholds.
#' @param family should be "gaussian" if y is numeric, but "binomial" if y is a binary factor variable.
#'
#' @return a vector of selected indices

# linear interaction KO filter (presented in Section 2.3.1)
linear_model_predictive_filter = function(X, X_tilde, y, t, fdr_nominal, family){
  
  p = dim(X)[2]
  
  # Combine the original Xs and their interactions T:X
  X_original_interactions = matrix(0, dim(X)[1], dim(X)[2])
  for (col in 1:p){
    X_original_interactions[,col]=  X[,col]*t
  }
  
  # Generate interactions with knockoffs the knockoffs Xs and their interactions T:X
  X_knockoffs_interactions = matrix(0, dim(X)[1], dim(X)[2])
  for (col in 1:p){
    X_knockoffs_interactions[,col] =  X_tilde[,col]*t
  }
  
  # Combine data
  X_combined = cbind(X, X_tilde, X_original_interactions, X_knockoffs_interactions, t)
  
  
  # Derive the variable importances
  X_combined = as.matrix(X_combined)
  y = as.matrix(y)
  models_glmnet_cv = glmnet::cv.glmnet(X_combined, y, family = family, nfolds=10, alpha = 1.0,  standardize=TRUE)
 
  importance_scores = coef(models_glmnet_cv, s = "lambda.1se")[2:(dim(X_combined)[2]+1)] # Ignore intercept
  
  # Derive W statistic only for the interaction terms
  W = abs(importance_scores[(2*p+1):(3*p)]) - abs(importance_scores[(3*p+1):(4*p)])
  
  # Calculate the threshold
  threshold = knockoff::knockoff.threshold(W, fdr=fdr_nominal)
  # Find the selected features
  selected_features = which(W >= threshold)

  # NW:
  # return indices of selected features, but also return W-statistic of all features and the threshold
  return(list("selected_features_indices" = selected_features, "W_Statistic" = W,"threshold"=threshold))
  
  
}


# CF variable importance KO filter (presented in Section 2.3.2)
# The family argument is left NULL since the same function can be used for both continuous or binary outcomes
causal_forest_predictive_filter = function(X, X_tilde, y, t, fdr_nominal, family = NULL){

  p = dim(X)[2]

  
  # Build causal forests using initial features and knockoffs
  c.forest = grf::causal_forest(cbind(X, X_tilde), y, t, mtry = 2*p)
  
  # Derive the variable importance scores for initial featurs and knockoffs
  # We are using the default parameters for decay.exponent and max.depth
  importance_scores = grf::variable_importance(c.forest, decay.exponent = 2, max.depth = 5)
  
  # Calculate the W statistic as the difference
  W = abs(importance_scores[1:p]) - abs(importance_scores[(p+1):(2*p)])
  
  # Calculate the threshold
  threshold = knockoff::knockoff.threshold(W, fdr=fdr_nominal)
  # Find the selected features
  selected_features = which(W >= threshold)
  
  # NW:
  # return indices of selected features, but also return W-statistic of all features and the threshold
  return(list("selected_features_indices" = selected_features, "W_Statistic" = W, "threshold"=threshold))
}
