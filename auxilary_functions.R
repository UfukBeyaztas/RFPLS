#_____________________
#______Packages_______
#_____________________
library(rlmDataDriven)
library(MASS)
library(fda)
library(expm)
library(regRSM)
#_____________________
#_____________________

options(warn = -1)

#_____________________________________________________________________________________________________________________________
#___________________________________________________________Design_matrix_____________________________________________________
#_____________________________________________________________________________________________________________________________

getAmat = function(X_tr, X_te, nbf){
  np = length(X_tr)
  
  dtp = list()
  for(i in 1:np)
    dtp[[i]] = seq(0, 1, length=ncol(X_tr[[i]]))
  
  Bsp_b = list()
  Bsp_f = list()
  for(i in 1:np){
    Bsp_b[[i]] = create.bspline.basis(c(0,1), nbasis = nbf[i])
    Bsp_f[[i]] = eval.basis(dtp[[i]], Bsp_b[[i]])
  }
  
  Inn_prod = list()
  for(i in 1:np)
    Inn_prod[[i]] = inprod(Bsp_b[[i]], Bsp_b[[i]])
  
  Inn_prod_sqrt = list()
  
  for(i in 1:np)
    Inn_prod_sqrt[[i]] = sqrtm(Inn_prod[[i]])
  
  w_tr = list()
  w_te = list()
  for(i in 1:np){
    w_tr[[i]] = matrix(dtp[[i]], nrow = nrow(X_tr[[i]]), ncol = ncol(X_tr[[i]]), byrow=T)
    w_te[[i]] = matrix(dtp[[i]], nrow = nrow(X_te[[i]]), ncol = ncol(X_te[[i]]), byrow=T)
  }
  
  Wx_tr = list()
  Wx_te = list()
  for(i in 1:np){
    Wx_tr[[i]] = t(smooth.basis(argvals=t(w_tr[[i]]), y=t(X_tr[[i]]), fdParobj=Bsp_b[[i]])$fd$coefs)
    Wx_te[[i]] = t(smooth.basis(argvals=t(w_te[[i]]), y=t(X_te[[i]]), fdParobj=Bsp_b[[i]])$fd$coefs)
  }
  
  Amat_tr = list()
  Amat_te = list()
  for(i in 1:np){
    Amat_tr[[i]] = Wx_tr[[i]] %*% Inn_prod_sqrt[[i]]
    Amat_te[[i]] = Wx_te[[i]] %*% Inn_prod_sqrt[[i]]
  }
  
  return(list(Amat_tr = Amat_tr, Amat_te = Amat_te,
              Bsp_funs = Bsp_f, InnProds = Inn_prod_sqrt))
}

#_____________________________________________________________________________________________________________________________
#_____________________________________________________________________________________________________________________________


#_______________________________________________________________________________________________________________________________________________________________
#_____________________________________________________________________FPLS_function_____________________________________________________________________________
#_______________________________________________________________________________________________________________________________________________________________
fpls_fun = function(Y, X_tr, X_te, npls, method = c("classical", "robust")){
  np = length(X_tr)
  N = dim(X_tr[[1]])[1]
  M = numeric()
  for(im in 1:np)
    M[im] = dim(X_tr[[im]])[2]
  
  pred_vals = list()
  fit_vals = list()
  coef_mat = list()
  comps = list()
  comps_test = list()
  base = list()
  
  for(ii in 1:npls){
    for(jm in 1:np){
      it_coefs = matrix(NA, nrow = M[jm], ncol = 1)
      
      for(j in 1:M[jm]){
        s_X = matrix(X_tr[[jm]][,j], nrow = N, ncol = 1)
        
        it_model = lm(Y~s_X)
        
        if(method == "classical"){
          it_coefs[j,] = it_model$coefficients[2]
        }else if(method == "robust"){
          it_model_rb = rlm(Y~s_X, psi = psi.bisquare, c = 4.685, maxit = 1000)
          it_model_rlmdd = rlmDD(Y, s_X, it_model$coef, it_model_rb$coef, method = "Bisquare",plot = FALSE)
          it_coefs[j,] = it_model_rlmdd$esti$coefficients[2]
        }
      }
      
      it_coefs = it_coefs / sqrt(sum(it_coefs^2))
      
      if(ii == 1){
        comps[[jm]] = X_tr[[jm]] %*% it_coefs
        comps_test[[jm]] = X_te[[jm]] %*% it_coefs
        base[[jm]] = it_coefs
      }else{
        comps[[jm]] = cbind(comps[[jm]], X_tr[[jm]] %*% it_coefs)
        comps_test[[jm]] = cbind(comps_test[[jm]], X_te[[jm]] %*% it_coefs)
        base[[jm]] = cbind(base[[jm]], it_coefs)
      }
      
      for(jj in 1:M[jm]){
        X_tr[[jm]][,jj] = X_tr[[jm]][,jj] - mean(X_tr[[jm]][,jj]) - as.numeric(cov(X_tr[[jm]][,jj], comps[[jm]][,ii]) / 
                                                                                 var(comps[[jm]][,ii])) * (comps[[jm]][,ii] - mean(comps[[jm]][,ii]))
        
        X_te[[jm]][,jj] = X_te[[jm]][,jj] - mean(X_te[[jm]][,jj]) - as.numeric(cov(X_te[[jm]][,jj], comps_test[[jm]][,ii]) / 
                                                                                 var(comps_test[[jm]][,ii])) * (comps_test[[jm]][,ii] - mean(comps_test[[jm]][,ii]))
      }
    }
    
    fit_model = lm(Y~do.call(cbind, comps))
    
    if(method == "classical"){
      fit_coefs = fit_model$coefficients
    }else if(method == "robust"){
      fit_model_rb = rlm(Y~do.call(cbind, comps), psi = psi.bisquare, c = 4.685, maxit = 1000)
      fit_model_rlmdd = rlmDD(Y, do.call(cbind, comps), fit_model$coef, fit_model_rb$coef, method = "Bisquare",plot = FALSE)
      fit_coefs = fit_model_rlmdd$esti$coefficients
    }
    
    coef_mat[[ii]] = fit_coefs
    
    fit_vals[[ii]] = cbind(1, do.call(cbind, comps)) %*% fit_coefs
    pred_vals[[ii]] = cbind(1, do.call(cbind, comps_test)) %*% fit_coefs
  }
  
  opt_fits = fit_vals[[npls]]
  opt_preds = pred_vals[[npls]]
  opt_coefs = coef_mat[[npls]]
  
  return(list(all_preds = pred_vals, fits_opt = opt_fits, preds_opt = opt_preds, coefs = opt_coefs, plscomps = comps))
}

#_______________________________________________________________________________________________________________________________________________________________
#_______________________________________________________________________________________________________________________________________________________________


#___________________________________________________________________________________________________
#__________________________________Estimated_coefficient_functions__________________________________
#___________________________________________________________________________________________________

getCoef = function(object, plsComps, coefs){
  
  Amats = object$Amat_tr
  Bs_funs = object$Bsp_funs
  Inn_prd_sqrts = object$InnProds
  np = length(Amats)
  npls = dim(plsComps[[1]])[2]
  
  for(ij in 1:np)
    Amats[[ij]] = scale(Amats[[ij]])
  
  B0 = coefs[1]
  Brest = matrix(coefs[-1], ncol = npls, byrow = TRUE)
  
  V = list()
  for(i in 1:np)
    V[[i]] = ginv(t(Amats[[i]]) %*% Amats[[i]]) %*% t(Amats[[i]]) %*% plsComps[[i]]
  
  Cf = list()
  for(i in 1:np)
    Cf[[i]] = t(solve(Inn_prd_sqrts[[i]])) %*% V[[i]] %*% Brest[i,]
  
  Bhat = list()
  for(i in 1:np)
    Bhat[[i]] = t(Cf[[i]]) %*% t(Bs_funs[[i]])
  
  return(list(Bhat = Bhat))
}

#___________________________________________________________________________________________________
#___________________________________________________________________________________________________


#__________________________________________________________________________________________
#______________________________________Variable_selection__________________________________
#__________________________________________________________________________________________
var_select = function(Y, X, nbf=5){
  np = length(X)
  
  dtp = list()
  for(i in 1:np)
    dtp[[i]] = seq(0, 1, length=ncol(X[[i]]))
  
  Bsp_b = list()
  Bsp_f = list()
  for(i in 1:np){
    Bsp_b[[i]] = create.bspline.basis(c(0,1), nbasis = nbf)
    Bsp_f[[i]] = eval.basis(dtp[[i]], Bsp_b[[i]])
  }
  
  Inn_prod = list()
  for(i in 1:np)
    Inn_prod[[i]] = inprod(Bsp_b[[i]], Bsp_b[[i]])
  
  Inn_prod_sqrt = list()
  for(i in 1:np)
    Inn_prod_sqrt[[i]] = sqrtm(Inn_prod[[i]])
  
  w = list()
  for(i in 1:np)
    w[[i]] = matrix(dtp[[i]], nrow = nrow(X[[i]]), ncol = ncol(X[[i]]), byrow=T)
  
  Wx = list()
  for(i in 1:np)
    Wx[[i]] = t(smooth.basis(argvals=t(w[[i]]), y=t(X[[i]]), fdParobj=Bsp_b[[i]])$fd$coefs)
  
  Amat = list()
  for(i in 1:np)
    Amat[[i]] = Wx[[i]] %*% Inn_prod_sqrt[[i]]
  
  Dmat = do.call(cbind, Amat)
  rsm_model = regRSM(y = Y, x = Dmat, B = 1000)
  empt_mat = rep(NA, nbf*length(X))
  empt_mat[sort(rsm_model$model)] = sort(rsm_model$model)
  empt_mat = matrix(empt_mat, nrow = length(X), ncol = nbf, byrow = TRUE)
  variable_indices = which(rowSums(empt_mat, na.rm = T) > 0)
  
  return(variable_indices)
}

#__________________________________________________________________________________________
#__________________________________________________________________________________________


#___________________________________________________________________________________________
#________________________________K-fold_cross-validation____________________________________
#___________________________________________________________________________________________

nfoldCV = function(Y, X, npls_max, nfold, method = c("classical", "robust")){
  n = dim(X[[1]])[1]
  fp = length(X)
  oIndx = 1:n
  samp_dat = oIndx[sample(n)]
  folds = cut(seq(1,n), breaks = nfold, labels=FALSE)
  
  nfMspe = numeric()
  
  for(i in 1:nfold){
    indx_test = which(folds == i, arr.ind = TRUE)
    
    nfY_tr = Y[-indx_test]
    nfY_te = Y[indx_test]
    
    nfX_tr = list()
    nfX_te = list()
    
    for(j in 1:fp){
      nfX_tr[[j]] = X[[j]][-indx_test,]
      nfX_te[[j]] = X[[j]][indx_test,]
    }
    
    fold_fpls = fpls_fun(Y = nfY_tr, X_tr = nfX_tr, X_te = nfX_te,
                         npls = npls_max, method = method)
    
    mspe_i = numeric()
    for(k in 1: npls_max)
      mspe_i[k] = mean((nfY_te - fold_fpls$all_preds[[k]])^2)
    
    nfMspe[i] = which.min(mspe_i)
    
  }
  return(round(mean(nfMspe)))
}

#___________________________________________________________________________________________
#___________________________________________________________________________________________

#___________________________________________________________________________________________
#__________________________________CPD_and_interval_score___________________________________
#___________________________________________________________________________________________

interval_score <- function(holdout, lb, ub, alpha){
  lb_ind = ifelse(holdout < lb, 1, 0)
  ub_ind = ifelse(holdout > ub, 1, 0)
  score = (ub - lb) + 2/alpha * ((lb - holdout) * lb_ind + (holdout - ub) * ub_ind)
  cover = 1 - (length(which(lb_ind == 1)) + length(which(ub_ind == 1)))/length(holdout)
  cpd = abs(cover - (1 - alpha))
  return(c(mean(score), cpd))
}
#___________________________________________________________________________________________
#___________________________________________________________________________________________
