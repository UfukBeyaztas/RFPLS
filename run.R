source("bsq_comps.R") # for variable selection procedure
  source("bsq_funs.R") # for the RFPLS method
  source("rsm.R") # functions to be used in the RSM procedure
  
  # Load data
  load("ash_content.RData") # Y
  load("spectra.RData") # X(t)
  
  # Split data
  index = sample(1:length(Y), 168, replace = FALSE)
  
  Y_train = Y[index]
  Y_test = Y[-index]
  X_train = list()
  X_test = list()
  for(i in 1:length(X_all)){
    X_train[[i]] = X_all[[i]][index,]
    X_test[[i]] = X_all[[i]][-index,]
  }
  
  # Variable selection
  npls = 2
  vs_model = bsq_base(Y = Y_train, X = X_train, npls = npls)
  projections = do.call(cbind, vs_model)
  rsm_model = regRSM(y = Y_train, x = projections, B = 1000)
  empt_mat = rep(NA, npls*length(X_train))
  empt_mat[sort(rsm_model$model)] = sort(rsm_model$model)
  empt_mat = matrix(empt_mat, nrow = length(X_train), ncol = npls, byrow = TRUE)
  variable_indices = which(rowSums(empt_mat, na.rm = T) > 0)
  
  # Selected functional predictors
  Xs_train = list()
  Xs_test = list()
  for(i in 1:length(variable_indices)){
    Xs_train[[i]] = X_train[[variable_indices[i]]]
    Xs_test[[i]] = X_test[[variable_indices[i]]]
  }
  
  # Final model
  fmodel = bsq_fun(Y = Y_train, X = Xs_train, X_test = Xs_test, npls = 5)
  
  # MSPE
  mean((Y_test - fmodel$preds)^2)
