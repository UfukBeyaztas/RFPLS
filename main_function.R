main_fun = function(Y, X_tr, X_te, t_index, nfold, npls_max, nbf, method = c("classical", "robust"),
                    model = c("full", "true", "selected"), CV = TRUE, npls_use){
  # Y is a vector of scalar response
  # X_tr is a list containing functional predictors in the training sample
  # X_te is a list containing functional predictors in the test sample
  # t_index is a vector containing the index of true functional predictors
  # nfold: the value k in the k-fold cross-validation
  # npls_max is the maximum number of PLS components
  # nbf is the number of basis functions to estimate scalar-on-function regression model
  # method is a selection to indicate which method to use 
  # model is a selection to indicate which model to use 
  # CV: if true, the k-fold cross validation is used to determine the optimum number of PLS components
  # otherwise, npls_use is used
  
  
  method = match.arg(method)
  model = match.arg(model)

  if(model == "full"){
    if(CV){
      npls_optim = nfoldCV(Y = Y, X = X_tr, npls_max = npls_max, nfold = nfold, method = method)
    }else{
      npls_optim = npls_use
    }
    
    Amats_fpls = getAmat(X_tr = X_tr, X_te = X_te, nbf = nbf)
    fpls_model = fpls_fun(Y = Y, X_tr = Amats_fpls$Amat_tr, X_te = Amats_fpls$Amat_te,
                          npls = npls_optim, method = method)
    
    fits = fpls_model$fits_opt
    preds = fpls_model$preds_opt
    est_coefs = getCoef(object = Amats_fpls, plsComps = fpls_model$plscomps, coefs = fpls_model$coefs)$Bhat
    var_used = "All the variables were used"
  }else if(model == "true"){
    X_tr = X_tr[t_index]
    X_te = X_te[t_index]
    
    if(CV){
      npls_optim = nfoldCV(Y = Y, X = X_tr, npls_max = npls_max, nfold = nfold, method = method)
    }else{
      npls_optim = npls_use
    }
    
    Amats_fpls = getAmat(X_tr = X_tr, X_te = X_te, nbf = nbf)
    fpls_model = fpls_fun(Y = Y, X_tr = Amats_fpls$Amat_tr, X_te = Amats_fpls$Amat_te,
                          npls = npls_optim, method = method)
    
    fits = fpls_model$fits_opt
    preds = fpls_model$preds_opt
    est_coefs = getCoef(object = Amats_fpls, plsComps = fpls_model$plscomps, coefs = fpls_model$coefs)$Bhat
    var_used = "Only the significant variables were used"
  }else if(model == "selected"){
    var_indx = var_select(Y = Y, X = X_tr)
    
    X_tr = X_tr[var_indx]
    X_te = X_te[var_indx]
    
    if(CV){
      npls_optim = nfoldCV(Y = Y, X = X_tr, npls_max = npls_max, nfold = nfold, method = method)
    }else{
      npls_optim = npls_use
    }
    
    Amats_fpls = getAmat(X_tr = X_tr, X_te = X_te, nbf = nbf)
    fpls_model = fpls_fun(Y = Y, X_tr = Amats_fpls$Amat_tr, X_te = Amats_fpls$Amat_te,
                          npls = npls_optim, method = method)
    
    fits = fpls_model$fits_opt
    preds = fpls_model$preds_opt
    est_coefs = getCoef(object = Amats_fpls, plsComps = fpls_model$plscomps, coefs = fpls_model$coefs)$Bhat
    var_used = c("Selected variables are:", var_indx)
  }
  
  return(list(fitvals = fits, predictions = preds, coefficients = est_coefs, var_used = var_used))
}
