main_fun = function(Y, X_tr, X_te, t_index, npls_max, npca_max, nbf, nbf_vs, rangeval,
                    method = c("classical", "robust"), model = c("full", "true", "selected"),
                    fmethod = c("fpls", "fpca")){
  
  # Y is a vector of scalar response
  # X_tr is a list containing functional predictors in the training sample
  # X_te is a list containing functional predictors in the test sample
  # t_index is a vector containing the index of true functional predictors
  # npls_max is the maximum number of PLS components
  # npca_max is the maximum number of PCA components
  # nbf is the number of basis functions to estimate scalar-on-function regression model
  # nbf_vs is a vector containing number of basis functions for each functional predictor used in the variable selection procedure
  # rangeval is a list containing the interval values where the functional dataset is evaluated
  # method is a selection to indicate which method to use 
  # model is a selection to indicate which model to use 
  # fmethod is a selection to indicate whether to use  FPLS or FPCA
  
  
  method = match.arg(method)
  model = match.arg(model)
  fmethod = match.arg(fmethod)
  
  if(model == "full"){
    if(fmethod == "fpls"){
      BIC_res = BIC_nc(Y = Y, X = X_tr, hmax = npls_max, method = method)
      ncomp_optim = BIC_res$hopt
      BIC_val = BIC_res$BIC_val
      Amats_fpls = getAmat(X_tr = X_tr, X_te = X_te, nbf = nbf, rangeval = rangeval)
      fpls_model = fpls_fun(Y = Y, X_tr = Amats_fpls$Amat_tr, X_te = Amats_fpls$Amat_te,
                            npls = ncomp_optim, method = method)
      
      fits = fpls_model$fits_opt
      preds = fpls_model$preds_opt
      est_coefs = getCoef(object = Amats_fpls, plsComps = fpls_model$plscomps, coefs = fpls_model$coefs)$Bhat
    }else if(fmethod == "fpca"){
      BIC_res = BIC_nc_pca(Y = Y, X = X_tr, npca_max = npca_max, nbf = nbf, rangeval = rangeval)
      ncomp_optim = BIC_res$hopt
      BIC_val = BIC_res$BIC_val
      fpca_model = fpca(Y = Y, X_tr = X_tr, X_te = X_te, npca = ncomp_optim, nbf = nbf, rangeval = rangeval)
      fits = fpca_model$fits
      preds = fpca_model$preds
      est_coefs = getCoef_fpca(object = fpca_model, npca = ncomp_optim)$coefs
    }
  }else if(model == "true"){
    X_tr = X_tr[t_index]
    X_te = X_te[t_index]
    
    if(fmethod == "fpls"){
      BIC_res = BIC_nc(Y = Y, X = X_tr, hmax = npls_max, method = method)
      ncomp_optim = BIC_res$hopt
      BIC_val = BIC_res$BIC_val
      Amats_fpls = getAmat(X_tr = X_tr, X_te = X_te, nbf = nbf, rangeval = rangeval)
      fpls_model = fpls_fun(Y = Y, X_tr = Amats_fpls$Amat_tr, X_te = Amats_fpls$Amat_te,
                            npls = ncomp_optim, method = method)
      
      fits = fpls_model$fits_opt
      preds = fpls_model$preds_opt
      est_coefs = getCoef(object = Amats_fpls, plsComps = fpls_model$plscomps, coefs = fpls_model$coefs)$Bhat
    }else if(fmethod == "fpca"){
      BIC_res = BIC_nc_pca(Y = Y, X = X_tr, npca_max = npca_max, nbf = nbf, rangeval = rangeval)
      ncomp_optim = BIC_res$hopt
      BIC_val = BIC_res$BIC_val
      fpca_model = fpca(Y = Y, X_tr = X_tr, X_te = X_te, npca = ncomp_optim, nbf = nbf, rangeval = rangeval)
      fits = fpca_model$fits
      preds = fpca_model$preds
      est_coefs = getCoef_fpca(object = fpca_model, npca = ncomp_optim)$coefs
    }
  }else if(model == "selected"){
    var_indx = var_select(Y = Y, X = X_tr, nbf = nbf_vs, rangeval = rangeval)
    
    X_tr = X_tr[var_indx]
    X_te = X_te[var_indx]
    
    rangeval2 = rangeval[var_indx]
    
    if(fmethod == "fpls"){
      BIC_res = BIC_nc(Y = Y, X = X_tr, hmax = npls_max, method = method)
      ncomp_optim = BIC_res$hopt
      BIC_val = BIC_res$BIC_val
      Amats_fpls = getAmat(X_tr = X_tr, X_te = X_te, nbf = nbf, rangeval = rangeval2)
      fpls_model = fpls_fun(Y = Y, X_tr = Amats_fpls$Amat_tr, X_te = Amats_fpls$Amat_te,
                            npls = ncomp_optim, method = method)
      
      fits = fpls_model$fits_opt
      preds = fpls_model$preds_opt
      est_coefs = getCoef(object = Amats_fpls, plsComps = fpls_model$plscomps, coefs = fpls_model$coefs)$Bhat
    }else if(fmethod == "fpca"){
      BIC_res = BIC_nc_pca(Y = Y, X = X_tr, npca_max = npca_max, nbf = nbf, rangeval = rangeval2)
      ncomp_optim = BIC_res$hopt
      BIC_val = BIC_res$BIC_val
      fpca_model = fpca(Y = Y, X_tr = X_tr, X_te = X_te, npca = ncomp_optim, nbf = nbf, rangeval = rangeval2)
      fits = fpca_model$fits
      preds = fpca_model$preds
      est_coefs = getCoef_fpca(object = fpca_model, npca = ncomp_optim)$coefs
    }
  }
  
  return(list(fitvals = fits, predictions = preds, coefficients = est_coefs, BIC_val = BIC_val, ncomp_optim = ncomp_optim))
}
