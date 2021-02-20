library(rlmDataDriven)
library(MASS)
options(warn = -1)

bsq_fun = function(Y, X, X_test, npls){
  # Y is a vector of scalar response
  # X is a list of functional predictors in the training sample
  # X_test is a list of functional predictors in the test sample
  # npls is the number of PLS iterations

  s_fun = function(data){
    return((data - mean(data)) / sd(data))
  }
  
  np = length(X)
  N = length(Y)
  M = ncol(X[[1]])
  
  N_train = round(0.5*N)
  Y_train = Y[1:N_train]
  Y_test = Y[-(1:N_train)]
  
  ########## Determination of number of PLS basis
  
  Xs = list()
  Xs_test = list()

  
  for(i in 1:np){
    Xs[[i]] = scale(X[[i]][1:N_train,], center = TRUE, scale = FALSE)
    Xs_test[[i]] = scale(X[[i]][-(1:N_train),], center = TRUE, scale = FALSE)
  }

  pred_vals = list()

  for(ii in 1:npls){
    
    it_coefs = matrix(NA, nrow = M, ncol = np)

    for(j in 1:M){
      s_X = matrix(NA, nrow = N_train, ncol = np)

      for(jj in 1:np)
        s_X[,jj] = s_fun(Xs[[jj]][,j])

      it_model = lm(Y_train~s_X)

      it_model_rb = rlm(Y_train~s_X, psi = psi.bisquare, c = 4.685)

      it_model_rlmdd = rlmDD(Y_train, s_X, it_model$coef, it_model_rb$coef, method = "Bisquare",plot = FALSE)

      it_coefs[j,] = it_model_rlmdd$esti$coefficients[2:(np+1)]
    }
    
    for(j in 1:np)
      it_coefs[,j] = it_coefs[,j] / sqrt(sum(it_coefs[,j]^2))

    if(ii == 1){
      comps = list()
      comps_test = list()
      base = list()

      for(j in 1:np){
        comps[[j]] = Xs[[j]] %*% it_coefs[,j]
        comps_test[[j]] = Xs_test[[j]] %*% it_coefs[,j]
        base[[j]] = it_coefs[,j]
      }
    }else{
      for(j in 1:np){
        comps[[j]] = cbind(comps[[j]], Xs[[j]] %*% it_coefs[,j])
        comps_test[[j]] = cbind(comps_test[[j]], Xs_test[[j]] %*% it_coefs[,j])
        base[[j]] = cbind(base[[j]], it_coefs[,j])
      }
    }
    
    for(j in 1:np){
      for(jj in 1:M){
        Xs[[j]][,jj] = Xs[[j]][,jj] - mean(Xs[[j]][,jj]) - as.numeric(cov(Xs[[j]][,jj], comps[[j]][,ii]) / 
                                                                        var(comps[[j]][,ii])) * (comps[[j]][,ii] - mean(comps[[j]][,ii]))
        
        Xs_test[[j]][,jj] = Xs_test[[j]][,jj] - mean(Xs_test[[j]][,jj]) - as.numeric(cov(Xs_test[[j]][,jj], comps_test[[j]][,ii]) / 
                                                                        var(comps_test[[j]][,ii])) * (comps_test[[j]][,ii] - mean(comps_test[[j]][,ii]))
      }
    }

    fit_model = lm(Y_train~do.call(cbind, comps))

    fit_model_rb = rlm(Y_train~do.call(cbind, comps), psi = psi.bisquare, c = 4.685)

    fit_model_rlmdd = rlmDD(Y_train, do.call(cbind, comps), fit_model$coef, fit_model_rb$coef, method = "Bisquare",plot = FALSE)

    fit_coefs = fit_model_rlmdd$esti$coefficients

    pred_vals[[ii]] = cbind(1, do.call(cbind, comps_test)) %*% fit_coefs
  }
  
  mspe = numeric()
  for(i in 1:npls)
    mspe[i] = mean((Y_test - pred_vals[[i]])^2)
  
  npls_opt = which.min(mspe)

  
  ########## Model based on optimum number of PLS
  
  
  Xs = list()
  Xs_test = list()
  
  
  for(i in 1:np){
    Xs[[i]] = scale(X[[i]], center = TRUE, scale = FALSE)
    Xs_test[[i]] = scale(X_test[[i]], center = TRUE, scale = FALSE)
  }
  
  
  pred_vals = list()
  
  for(ii in 1:npls_opt){
    
    it_coefs = matrix(NA, nrow = M, ncol = np)
    
    for(j in 1:M){
      s_X = matrix(NA, nrow = N, ncol = np)
      
      for(jj in 1:np)
        s_X[,jj] = s_fun(Xs[[jj]][,j])
      
      it_model = lm(Y~s_X)
      it_model_rb = rlm(Y~s_X, psi = psi.bisquare, c = 4.685)
      it_model_rlmdd = rlmDD(Y, s_X, it_model$coef, it_model_rb$coef, method = "Bisquare",plot = FALSE)
      
      it_coefs[j,] = it_model_rlmdd$esti$coefficients[2:(np+1)]
    }
    
    for(j in 1:np)
      it_coefs[,j] = it_coefs[,j] / sqrt(sum(it_coefs[,j]^2))
    
    if(ii == 1){
      comps = list()
      comps_test = list()
      base = list()
      
      for(j in 1:np){
        comps[[j]] = Xs[[j]] %*% it_coefs[,j]
        comps_test[[j]] = Xs_test[[j]] %*% it_coefs[,j]
        base[[j]] = it_coefs[,j]
      }
    }else{
      for(j in 1:np){
        comps[[j]] = cbind(comps[[j]], Xs[[j]] %*% it_coefs[,j])
        comps_test[[j]] = cbind(comps_test[[j]], Xs_test[[j]] %*% it_coefs[,j])
        base[[j]] = cbind(base[[j]], it_coefs[,j])
      }
    }
    
    for(j in 1:np){
      for(jj in 1:M){
        Xs[[j]][,jj] = Xs[[j]][,jj] - mean(Xs[[j]][,jj]) - as.numeric(cov(Xs[[j]][,jj], comps[[j]][,ii]) / 
                                                                        var(comps[[j]][,ii])) * (comps[[j]][,ii] - mean(comps[[j]][,ii]))
        
        Xs_test[[j]][,jj] = Xs_test[[j]][,jj] - mean(Xs_test[[j]][,jj]) - as.numeric(cov(Xs_test[[j]][,jj], comps_test[[j]][,ii]) / 
                                                                        var(comps_test[[j]][,ii])) * (comps_test[[j]][,ii] - mean(comps_test[[j]][,ii]))
      }
    }
    
    fit_model = lm(Y~do.call(cbind, comps))
    fit_model_rb = rlm(Y~do.call(cbind, comps), psi = psi.bisquare, c = 4.685)
    fit_model_rlmdd = rlmDD(Y, do.call(cbind, comps), fit_model$coef, fit_model_rb$coef, method = "Bisquare",plot = FALSE)
    
    fit_coefs = fit_model_rlmdd$esti$coefficients
    
    pred_vals[[ii]] = cbind(1, do.call(cbind, comps_test)) %*% fit_coefs
  }
  
  pred_opt = pred_vals[[npls_opt]]

  return(list("basis" = base, "comps" = comps, "comps_test" = comps_test, "preds" = pred_opt))
}
