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
#_______________________________________________________Design_matrix_(FPLS)__________________________________________________
#_____________________________________________________________________________________________________________________________

getAmat = function(X_tr, X_te, nbf, rangeval){
  np = length(X_tr)
  
  dtp = list()
  for(i in 1:np)
    dtp[[i]] = seq(rangeval[[i]][1], rangeval[[i]][2], length=ncol(X_tr[[i]]))
  
  Bsp_b = list()
  Bsp_f = list()
  for(i in 1:np){
    Bsp_b[[i]] = create.bspline.basis(c(rangeval[[i]][1], rangeval[[i]][2]), nbasis = nbf[i])
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


#_________________________________________________________________________________________________________________________________________________________________________
#______________________________________________________________________________FPLS_function______________________________________________________________________________
#_________________________________________________________________________________________________________________________________________________________________________
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
  
  return(list(all_fits = fit_vals, fits_opt = opt_fits, preds_opt = opt_preds, coefs = opt_coefs, plscomps = comps))
}

#_________________________________________________________________________________________________________________________________________________________________________
#_________________________________________________________________________________________________________________________________________________________________________


#___________________________________________________________________________________________________
#_______________________________Estimated_coefficient_functions_(PLS)_______________________________
#___________________________________________________________________________________________________

getCoef = function(object, plsComps, coefs){
  
  Amats = object$Amat_tr
  Bs_funs = object$Bsp_funs
  Inn_prd_sqrts = object$InnProds
  np = length(Amats)
  npls = dim(plsComps[[1]])[2]
  
  
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
#_________________________FPCA_components_and_scores_for_regression________________________
#__________________________________________________________________________________________

getPCA = function(X_tr, X_te, nbasis, ncomp, rangeval){
  data = rbind(X_tr, X_te)
  n_train = dim(X_tr)[1]
  n_test = dim(X_te)[1]
  n = dim(data)[1]
  p = dim(data)[2]
  dimnames(data)=list(as.character(1:n), as.character(1:p))
  grid_points = seq(rangeval[1], rangeval[2], length.out = p)
  bs_basis = create.bspline.basis(rangeval, nbasis = nbasis)
  evalbase = eval.basis(grid_points, bs_basis)
  fdobj = fdPar(bs_basis, int2Lfd(2), lambda=0)
  pcaobj = smooth.basisPar(grid_points, t(data), bs_basis, Lfdobj=NULL, lambda=0)$fd
  dpca = pca.fd(pcaobj[1:n_train,], nharm = ncomp, fdobj)
  PCAscore = dpca$scores
  PCAcoef = dpca$harmonics$coefs
  fdXnew = pcaobj[-(1:n_train),]
  fdXnew$coefs = t(scale(t(fdXnew$coefs), scale = F))
  PCAscore_test = inprod(fdXnew, dpca$harmonics)
  return(list(PCAcoef = PCAcoef, PCAscore = PCAscore, PCAscore_test = PCAscore_test, evalbase = evalbase))
}

#__________________________________________________________________________________________
#__________________________________________________________________________________________


#__________________________________________________________________________________________
#_________________________________FPCA_parameter_estimation________________________________
#__________________________________________________________________________________________
fpca_est = function(Y, sco_X){
  sco_X = cbind(1, sco_X)
  Bhat = ginv(t(sco_X) %*% sco_X) %*% t(sco_X) %*% Y
  return(Bhat)
}
#__________________________________________________________________________________________
#__________________________________________________________________________________________


#__________________________________________________________________________________________
#______________________________________FPCA_function_______________________________________
#__________________________________________________________________________________________

fpca = function(Y, X_tr, X_te, npca, nbf, rangeval){
  np = length(X_tr)
  
  fpca_X = list()
  for(fij in 1:np){
    fpca_X[[fij]] = getPCA(X_tr = X_tr[[fij]], X_te = X_te[[fij]],
                           nbasis = nbf[fij],  ncomp = npca, rangeval = rangeval[[fij]])
  }
  
  sco_tr = list()
  sco_te = list()
  for(i in 1:np){
    sco_tr[[i]] = fpca_X[[i]]$PCAscore
    sco_te[[i]] = fpca_X[[i]]$PCAscore_test
  }
  
  
  Bhat = fpca_est(Y = Y, sco_X = do.call(cbind, sco_tr))
  fits = cbind(1, do.call(cbind, sco_tr)) %*% Bhat
  preds =  cbind(1, do.call(cbind, sco_te)) %*% Bhat
  
  return(list(fits = fits, preds = preds, obj = fpca_X, Bhat = Bhat))
}

#__________________________________________________________________________________________
#__________________________________________________________________________________________


#__________________________________________________________________________________________
#__________________________________FPCA_coefficient_estimation_____________________________
#__________________________________________________________________________________________

getCoef_fpca = function(object, npca){
  np = length(object$obj)
  comp_X = list()
  evb = list()
  for(i in 1:np){
    comp_X[[i]] = object$obj[[i]]$PCAcoef
    evb[[i]] = object$obj[[i]]$evalbase
  }
  
  Bhat = as.matrix(object$Bhat[-1])
  
  
  V = list()
  for(i in 1:np)
    V[[i]] = comp_X[[i]]
  
  Cf = list()
  km = 1
  for(i in 1:np){
    Cf[[i]] = V[[i]] %*% as.matrix(Bhat[km: (km+npca-1),])
    km = i*npca+1
  }
  
  Bhat = list()
  for(i in 1:np)
    Bhat[[i]] = t(Cf[[i]]) %*% t(evb[[i]])
  
  return(list(coefs = Bhat))
}

#__________________________________________________________________________________________
#__________________________________________________________________________________________


#__________________________________________________________________________________________
#______________________________________Variable_selection__________________________________
#__________________________________________________________________________________________
var_select = function(Y, X, nbf=5, rangeval){
  np = length(X)
  
  dtp = list()
  for(i in 1:np)
    dtp[[i]] = seq(rangeval[[i]][1], rangeval[[i]][2], length=ncol(X[[i]]))
  
  Bsp_b = list()
  Bsp_f = list()
  for(i in 1:np){
    Bsp_b[[i]] = create.bspline.basis(c(rangeval[[i]][1], rangeval[[i]][2]), nbasis = nbf[[i]])
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
  empt_mat = rep(NA, sum(nbf))
  empt_mat[sort(rsm_model$model)] = sort(rsm_model$model)
  mat_list = list()
  mat_list[[1]] = empt_mat[1:nbf[1]]
  ik = nbf[1]
  for(i in 2:np){
    mat_list[[i]] = empt_mat[(ik+1):(ik+nbf[i])]
    ik = ik + nbf[i]
  }
  mat_ind = matrix(NA, ncol = 1, nrow = np)
  for(i in 1:np)
    mat_ind[i,] = sum(mat_list[[i]], na.rm = TRUE)
  variable_indices = which(rowSums(mat_ind, na.rm = T) > 0)
  
  return(variable_indices)
}

#__________________________________________________________________________________________
#__________________________________________________________________________________________


#__________________________________________________________________________________________
#_________________________________________BIC_function_____________________________________
#__________________________________________________________________________________________

BIC_fun = function(Y, Yfit, h, method = c("classical", "robust")){
  n = length(Y)
  arg_bic = (Y - Yfit)^2
  
  bic_index = sort.int(arg_bic, decreasing = FALSE,
                       index.return = TRUE)
  ntrim = round(0.9 * n)
  index_trunc = bic_index$ix[1:ntrim]
  
  if(method == "classical"){
    BIC_val = n * log(sum(arg_bic) / n) +
      h * log(n)
  }else if(method == "robust"){
    BIC_val = ntrim * log(sum(arg_bic[index_trunc]) / ntrim) +
      h * log(ntrim)
  }
  return(BIC_val)
}

#__________________________________________________________________________________________
#__________________________________________________________________________________________


#__________________________________________________________________________________________
#________________________________BIC_function_for_h_(PLS)__________________________________
#__________________________________________________________________________________________

BIC_nc = function(Y, X, hmax, method = c("classical", "robust")){
  
  method = match.arg(method)
  BIC_model = fpls_fun(Y = Y, X_tr = X, X_te = X,
                       npls = hmax, method = method)
  BIC_vec = numeric()
  for(i in 1:hmax){
    BIC_vec[i] = BIC_fun(Y = Y, Yfit = BIC_model$all_fits[[i]], h = i, method = method)
  }
  
  BIC_val = min(BIC_vec)
  hopt = which(BIC_vec == BIC_val)
  
  return(list(BIC_val = BIC_val, hopt = hopt))
}

#__________________________________________________________________________________________
#__________________________________________________________________________________________


#__________________________________________________________________________________________
#________________________________BIC_function_for_h_(FPCA)_________________________________
#__________________________________________________________________________________________

BIC_nc_pca = function(Y, X, npca_max, nbf = nbf,
                      rangeval, method = "classical"){
  
  fits_k = list()
  for(k in 1:npca_max){
    model_fpca = fpca(Y = Y, X_tr = X, X_te = X,
                      npca=k, nbf = nbf, rangeval = rangeval)
    fits_k[[k]] = model_fpca$fits
  }
  
  BIC_vec = numeric()
  for(i in 1:npca_max){
    BIC_vec[i] = BIC_fun(Y = Y, Yfit = fits_k[[i]], h = i, method = method)
  }
  
  BIC_val = min(BIC_vec)
  hopt = which(BIC_vec == BIC_val)
  
  return(list(BIC_val = BIC_val, hopt = hopt))
}

#__________________________________________________________________________________________
#__________________________________________________________________________________________


#_____________________________________________________
#_______________________Trimmed_MSE___________________
#_____________________________________________________
tmse = function(Y, fit){
  
  ersq = (Y - fit)^2
  trim_index = sort.int(ersq, decreasing = FALSE,
                        index.return = TRUE)
  ntrim = round(0.9 * length(Y))
  index_trunc = trim_index$ix[1:ntrim]
  newerr = ersq[index_trunc]
  
  trimmed_mse = mean(newerr)
  return(trimmed_mse)
}

#_____________________________________________________
#_____________________________________________________


#_____________________________________________________
#_____________________Trimmed_R^2_____________________
#_____________________________________________________

trsq = function(Y, fit){
  ersq = (Y - fit)^2
  trim_index = sort.int(ersq, decreasing = FALSE,
                        index.return = TRUE)
  ntrim = round(0.9 * length(Y))
  index_trunc = trim_index$ix[1:ntrim]
  
  newY = Y[index_trunc]
  newfir = fit[index_trunc]
  
  trimmed_rsq = cor(newY, newfir)^2
  return(trimmed_rsq)
}

#_____________________________________________________
#_____________________________________________________
