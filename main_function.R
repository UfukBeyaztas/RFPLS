main_fun = function(Y, X_tr, X_te, t_index, npls_max, npca_max, nbf, nfold, ncp_vs, rangeval,
                    method = c("classical", "robust"), model = c("full", "true", "selected"),
                    fmethod = c("fpls", "fpca")){
  
  # Y is a vector of scalar response
  # X_tr is a list containing functional predictors in the training sample
  # X_te is a list containing functional predictors in the test sample
  # t_index is a vector containing the index of true functional predictors
  # npls_max is the maximum number of PLS components
  # npca_max is the maximum number of PCA components
  # nbf is the number of basis functions to estimate scalar-on-function regression model
  # nfold is the value of "k" in the k-fold cross-validation
  # ncp_vs is a vector containing number of components for each functional predictor used in the variable selection procedure
  # rangeval is a list containing the interval values where the functional dataset is evaluated
  # method is a selection to indicate which method to use 
  # model is a selection to indicate which model to use 
  # fmethod is a selection to indicate whether to use  FPLS or FPCA
  
  
  method = match.arg(method)
  model = match.arg(model)
  fmethod = match.arg(fmethod)
  
  if(model == "full"){
    if(fmethod == "fpls"){
      ncomp_optim = nfoldCV(Y = Y, X = X_tr, ncomp = npls_max, nfold = nfold, nbf = nbf,
                            rangeval = rangeval, method = method, fmethod = fmethod)
      fpls_model = fpls_fun(Y = Y, X_tr = X_tr, X_te = X_te, nbf = nbf, npls = ncomp_optim, 
                            rangeval = rangeval, method = method)
      
      fits = fpls_model$fits
      preds = fpls_model$preds
      est_coefs = getCoef(object = fpls_model)$Bhat
    }else if(fmethod == "fpca"){
      ncomp_optim = nfoldCV(Y = Y, X = X_tr, ncomp = npca_max, nfold = nfold, nbf = nbf,
                            rangeval = rangeval, fmethod = fmethod)
      fpca_model = fpca(Y = Y, X_tr = X_tr, X_te = X_te, npca = ncomp_optim, nbf = nbf, rangeval = rangeval)
      fits = fpca_model$fits
      preds = fpca_model$preds
      est_coefs = getCoef_fpca(object = fpca_model, npca = ncomp_optim)$coefs
    }
  }else if(model == "true"){
    X_tr = X_tr[t_index]
    X_te = X_te[t_index]
    
    if(fmethod == "fpls"){
      ncomp_optim = nfoldCV(Y = Y, X = X_tr, ncomp = npls_max, nfold = nfold, nbf = nbf,
                            rangeval = rangeval, method = method, fmethod = fmethod)
      fpls_model = fpls_fun(Y = Y, X_tr = X_tr, X_te = X_te, nbf = nbf, npls = ncomp_optim, 
                            rangeval = rangeval, method = method)
      
      fits = fpls_model$fits
      preds = fpls_model$preds
      est_coefs = getCoef(object = fpls_model)$Bhat
    }else if(fmethod == "fpca"){
      ncomp_optim = nfoldCV(Y = Y, X = X_tr, ncomp = npca_max, nfold = nfold, nbf = nbf,
                            rangeval = rangeval, fmethod = fmethod)
      fpca_model = fpca(Y = Y, X_tr = X_tr, X_te = X_te, npca = ncomp_optim, nbf = nbf, rangeval = rangeval)
      fits = fpca_model$fits
      preds = fpca_model$preds
      est_coefs = getCoef_fpca(object = fpca_model, npca = ncomp_optim)$coefs
    }
  }else if(model == "selected"){
    var_indx = var_select(Y = Y, X = X_tr, nbf = nbf, ncp = ncp_vs,
                          rangeval = rangeval, method = method, fmethod = fmethod)
    
    X_tr = X_tr[var_indx]
    X_te = X_te[var_indx]
    
    rangeval2 = rangeval[var_indx]
    nbf2 = nbf[var_indx]
    
    if(fmethod == "fpls"){
      ncomp_optim = nfoldCV(Y = Y, X = X_tr, ncomp = npls_max, nfold = nfold, nbf = nbf2,
                            rangeval = rangeval2, method = method, fmethod = fmethod)
      fpls_model = fpls_fun(Y = Y, X_tr = X_tr, X_te = X_te, nbf = nbf2, npls = ncomp_optim, 
                            rangeval = rangeval2, method = method)
      
      fits = fpls_model$fits
      preds = fpls_model$preds
      est_coefs = getCoef(object = fpls_model)$Bhat
    }else if(fmethod == "fpca"){
      ncomp_optim = nfoldCV(Y = Y, X = X_tr, ncomp = npca_max, nfold = nfold, nbf = nbf2,
                            rangeval = rangeval2, fmethod = fmethod)
      fpca_model = fpca(Y = Y, X_tr = X_tr, X_te = X_te, npca = ncomp_optim, nbf = nbf2, rangeval = rangeval2)
      fits = fpca_model$fits
      preds = fpca_model$preds
      est_coefs = getCoef_fpca(object = fpca_model, npca = ncomp_optim)$coefs
    }
  }
  
  return(list(fitvals = fits, predictions = preds, coefficients = est_coefs, ncomp_optim = ncomp_optim))
}
