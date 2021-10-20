#_____________________
#______Packages_______
#_____________________
library(MASS)
library(fda)
library(expm)
library(regRSM)
library(robustX)  
library(matrixStats)
library(sprm)
library(rlmDataDriven)
#_____________________
#_____________________

options(warn = -1)

#_____________________________________________________________________________________________________________________________
#_______________________________________________________Design_matrix_(FPLS)__________________________________________________
#_____________________________________________________________________________________________________________________________

getAmat = function(X_tr, X_te, nbf, rangeval){
  np = length(X_tr)
  
  Amat_tr = list()
  Amat_te = list()
  Bsp_funs = list()
  InnProds = list()
  
  for(i in 1:np){
    data_tr = X_tr[[i]]
    data_te = X_te[[i]]
    n1 = dim(data_tr)[1]
    n2 = dim(data_te)[1]
    p = dim(data_tr)[2]
    dimnames(data_tr)=list(as.character(1:n1), as.character(1:p))
    dimnames(data_te)=list(as.character(1:n2), as.character(1:p))
    grid_points = seq(rangeval[[i]][1], rangeval[[i]][2], length.out = p)
    bs_basis = create.bspline.basis(rangeval[[i]], nbasis = nbf[i])
    Bsp_funs[[i]] = eval.basis(grid_points, bs_basis)
    InnProds[[i]] = sqrtm(inprod(bs_basis, bs_basis))
    fdobj = fdPar(bs_basis, int2Lfd(2), lambda=0)
    fdobj_tr = smooth.basisPar(grid_points, t(data_tr), bs_basis, Lfdobj=NULL, lambda=0)$fd
    fdobj_te = smooth.basisPar(grid_points, t(data_te), bs_basis, Lfdobj=NULL, lambda=0)$fd
    Amat_tr[[i]] = t(fdobj_tr$coefs) %*% InnProds[[i]]
    Amat_te[[i]] = t(fdobj_te$coefs) %*% InnProds[[i]]
  }
  
  return(list(Amat_tr = Amat_tr, Amat_te = Amat_te,
              Bsp_funs = Bsp_funs, InnProds = InnProds))
}

#_____________________________________________________________________________________________________________________________
#_____________________________________________________________________________________________________________________________


#_____________________________________________________________________________________________
#_____________________________________________SIMPS___________________________________________
#_____________________________________________________________________________________________

simpls = function(X, Xt, Y, a){
  Y = as.matrix(Y)
  n = nrow(X)
  k = ncol(X)
  m = ncol(Y)
  
  Ps = matrix(0, k, a)
  Cs = matrix(0, m, a)
  Rs = matrix(0, k, a)
  Ts = matrix(0, n, a)
  Tst = matrix(0, n, a)
  
  mx = apply(X, 2, mean)
  mxt = apply(Xt, 2, mean)
  sdx = apply(X, 2, sd)
  sdxt = apply(Xt, 2, sd)
  X = sapply(1:k, function(i) (X[,i]-mx[i]))
  Xt = sapply(1:k, function(i) (Xt[,i]-mxt[i]))
  my = apply(Y, 2, mean)
  sdy = apply(Y, 2, sd)
  Y = sapply(1:m, function(i) (Y[,i]-my[i]))
  S = t(X)%*%Y
  
  Snew = S
  for (i in 1:a){    
    rs = svd(Snew)$u[,1,drop=FALSE]
    rs = rs/norm(rs,type="F")
    ts = X %*% rs
    tst = Xt %*% rs
    tsn = ts/norm(ts,type="F")
    ps = t(X) %*% tsn
    cs = t(Y) %*% tsn
    Rs[,i] = rs
    Ts[,i] = ts
    Tst[,i] = tst
    Ps[,i] = ps
    Cs[,i] = cs
    Snew = Snew-Ps[,1:i] %*% solve(t(Ps[,1:i]) %*% Ps[,1:i]) %*% t(Ps[,1:i]) %*% Snew
  }
  beta = Rs %*% solve(t(Ps) %*% Rs) %*% t(Cs)
  
  return(list(beta = beta, Ts = Ts, Tst = Tst, Rs = Rs))
}
#_____________________________________________________________________________________________
#_____________________________________________________________________________________________


#_________________________________________________________________________________________________________
#_____________________________________________daprpr_function_____________________________________________
#_________________________________________________________________________________________________________
# This function was taken from the sprm R package: https://github.com/cran/sprm/blob/master/R/daprpr.R
daprpr = function(Data, center.type, scale.type){
  
  if(any(is.na(Data))==TRUE){stop("Routine cannot process NAs")}
  if(missing(center.type)){center.type="median"}
  if(missing(scale.type)){scale.type="qn"}
  if(center.type ==" l1median"){
    Data.center = l1median(Data)
  } else {
    Data.center = apply(Data, 2, center.type)
  }
  Data.scaled = (Data - matrix(Data.center, nrow=dim(Data)[1], ncol=dim(Data)[2], byrow=TRUE))
  if(!(scale.type=="no")){
    Data.scale = apply(Data,2,scale.type)
    if (any(1/Data.scale>1e19) | any(is.nan(Data.scale))){
      Data.scale = apply(Data,2,sd)
      if (any(1/Data.scale>1e19) | any(is.nan(Data.scale))){
        Data.scale = rep(1,ncol(Data))
        warning("Routine used scale.type='no' to avoide division by zero or infinity.")
      } else {
        warning("Routine used scale.type='sd' to avoide division by zero or infinity.")
      }
    }
    Data.scaled = Data.scaled/matrix(Data.scale,nrow=dim(Data)[1],ncol=dim(Data)[2],byrow=TRUE)
  } else {
    Data.scale = rep(1,ncol(Data))
  }
  attr(Data.scaled,"Center") = Data.center
  attr(Data.scaled,"Scale") = Data.scale
  attr(Data.scaled,"Type") = c(center.type,scale.type)
  return(Data.scaled)
}
#_________________________________________________________________________________________________________
#_________________________________________________________________________________________________________


#______________________________________________________________________________________________________________________________
#___________________________________________________________nipls_function_____________________________________________________
#______________________________________________________________________________________________________________________________
# This function was taken from the sprm R package: https://github.com/cran/sprm/blob/master/R/nipls.R
nipls = function(data, a){
  Z = data
  Xh = scale(as.matrix(Z[,2:ncol(Z)]),center=TRUE,scale=FALSE)
  X0 = Xh
  yh = scale(as.vector(Z[,1]),center=TRUE,scale=FALSE)
  my = attr(yh,"scaled:center")
  y0 = yh
  Tpls = NULL 
  W = NULL 
  P = NULL 
  C = NULL 
  B = NULL 
  bh = 0
  Xev = matrix(0,nrow=1,ncol=a)
  Yev = matrix(0,nrow=1,ncol=a)
  for(i in 1:a){
    wh = t(Xh)%*%yh
    wh = wh/norm(wh,"F")
    th = Xh%*%wh
    nth = norm(th,"F")
    ch = t(yh)%*%th/(nth^2)
    ph = t(Xh)%*%th/(nth^2)
    yh = yh - th * as.numeric(ch) 
    Xh = Xh - th%*%t(ph)
    W = cbind(W,wh)
    P = cbind(P,ph) 
    C = rbind(C,ch) 
    Tpls = cbind(Tpls,th)
    Xev[i] = (nth^2*norm(ph,"F")^2)/sum(X0^2)*100
    Yev[i] = sum(nth^2*as.numeric(ch^2))/sum(y0^2)*100
  }
  R = W %*% solve(t(P)%*%W)
  B = R%*%C
  yp = X0%*%B + my
  if (any(is.nan(Tpls))){
    stop("NaN generated in Tpls")
  }
  if(dim(Tpls)[1]==0){
    stop("No variables have been retained in Sparse PRM model!")
  }
  return(list(W=W,loadings=P, C=C, scores=Tpls, coefficients=B, Xev=Xev, Yev=Yev, fitted.values=yp, R=R, T2=X0%*%R))
}
#______________________________________________________________________________________________________________________________
#______________________________________________________________________________________________________________________________


#___________________________________________________________________________________________________________________________________________________-
#___________________________________________________________________prms_function___________________________________________________________________
#___________________________________________________________________________________________________________________________________________________
# This is a slightly modified version of the original prms function in: https://github.com/cran/sprm/blob/master/R/prms.R
prms2 = function (formula, data, Xt, a, fun="Hampel", probp1 = .95, hampelp2 = .975, hampelp3 = .999,
                  center = "median", scale = "qn", numit=100, prec=0.01){
  
  if(!class(formula)=="formula"){formula = formula(formula)}
  if(is.data.frame(data) | is.list(data)){
    mt = terms(formula, data=data)
    yname = dimnames(attr(mt,"factors"))[[1]][1]
    if(is.list(data)){
      datnames = names(data)
    } else {
      datnames = colnames(data)
    }
    ic = attr(mt, "intercept")
    if (ic==0){
      data = tryCatch({data = cbind(data[[which(datnames==yname)]], model.matrix(mt, data))},
                      error=function(err){
                        error = TRUE
                        return(error)
                      }) 
    } else{
      data = tryCatch({data = cbind(data[[which(datnames==yname)]],model.matrix(mt, data)[,-1])},
                      error=function(err){
                        error = TRUE
                        return(error)
                      }) 
    }
    if (is.logical(data)){
      stop("Data cannot be matched with formula.")
    } else {
      colnames(data)[1] = dimnames(attr(mt,"factors"))[[1]][1]
    }    
  } else {
    stop("Wrong data fromat.")
  }
  
  data = as.matrix(data)
  n = nrow(data)
  q = ncol(data)
  rnames = rownames(data) # restore original rownames in the end
  rownames(data) = 1:n  # 1:n are indices and names of w etc.
  p = q - 1 
  
  if(length(a)>1){
    warning("Only the first element of a is used.")
    a = a[1]
  }
  if(a>n|a>p){
    stop("The number of components a is too large.")
  }
  if (a<=0){
    stop("The number of components a has to be positive.")
  }
  if(!any(fun == c("Hampel", "Huber", "Fair"))){
    stop("Invalid weighting function. Choose Hampel, Huber or Fair for parameter fun.")
  }
  if(probp1>1|probp1<=0){
    stop("Parameter probp1 is a probability. Choose a value between 0 and 1")
  }
  if(fun=="Hampel"){
    if (!(probp1<hampelp2 & hampelp2<hampelp3 & hampelp3<=1)){
      stop("Wrong choise of parameters for Hampel function. Use 0<probp1<hampelp2<hampelp3<=1")
    }
  }
  
  dimensions = 0
  datam = data
  datamt = Xt
  
  datamc = daprpr(datam,center,scale)
  datamct = daprpr(Xt,center,scale)
  datac = attr(datamc,"Center")
  datact = attr(datamct,"Center")
  datas = attr(datamc,"Scale")
  datast = attr(datamct,"Scale")
  attr(datac,"Type") = center
  attr(datact,"Type") = center
  y0 = datam[,1]
  ys = datamc[,1]
  ns =nrow(datamc)
  qs = ncol(datamc)
  ps = qs - 1
  zerows = vector(length=0)
  wx = sqrt(apply(datamc[,2:qs]^2, 1, sum))
  wx = wx/median(wx)
  wy = abs(datamc[,1])
  
  if (length(wy)/2>sum(wy==0)){ # not too many zeros
    wy = wy/(1.4826*median(wy))
  } else{
    wy = wy/(1.4826*median(wy[wy!=0]))
  }
  
  probct = qnorm(probp1)
  
  hampelb = qnorm(hampelp2)
  hampelr = qnorm(hampelp3)
  wx[which(wx <= probct)] = 1 
  wx[which(wx > probct & wx <= hampelb)] = probct/abs(wx[which(wx > probct & wx <= hampelb)])
  wx[which(wx > hampelb & wx <= hampelr)] = probct*(hampelr-abs(wx[which(wx > hampelb & wx <= hampelr)]))/
    (hampelr -hampelb)*1/abs(wx[which(wx > hampelb & wx <= hampelr)])
  wx[which(wx > hampelr)] = 0
  wy[which(wy <= probct)] = 1 
  wy[which(wy > probct & wy <= hampelb)] = probct/abs(wy[which(wy > probct & wy <= hampelb)])
  wy[which(wy > hampelb & wy <= hampelr)] = probct*(hampelr-abs(wy[which(wy > hampelb & wy <= hampelr)]))/
    (hampelr -hampelb)*1/abs(wy[which(wy > hampelb & wy <= hampelr)])
  wy[which(wy > hampelr)] = 0 
  
  w = wx * wy
  if(any(w<1e-6)){
    w0 = which(w<1e-6)
    w = replace(w,list=w0,values=1e-6)
    we = w
  } else {
    wxe = wx
    wye = wy
    we = w
  }
  dataw = as.data.frame(datamc * sqrt(we))
  loops = 1
  rold = 10^-5
  difference = 1
  while ((difference > prec) && loops < numit) {    
    res.nipls = nipls(data=dataw,a=a)
    yp = fitted(res.nipls)
    r = datamc[,1] - yp
    b = coef(res.nipls)
    Tpls = res.nipls$scores/sqrt(we)
    if (length(r)/2>sum(r==0)){ 
      r = abs(r)/(1.4826*median(abs(r))) 
    } else{
      r = abs(r)/(1.4826*median(abs(r[r!=0])))
    }
    
    scalet = scale
    if(scale=="no"){scalet="qn"}
    dt = daprpr(Tpls,center,scalet)
    wtn = sqrt(apply(dt^2, 1, sum))
    wtn = wtn/median(wtn)
    
    probct = qnorm(probp1)
    hampelb = qnorm(hampelp2)
    hampelr = qnorm(hampelp3)
    wye = r
    wye[which(r <= probct)] = 1
    wye[which(r > probct & r <= hampelb)] = probct/abs(r[which(r > probct & r <= hampelb)])
    wye[which(r > hampelb & r <= hampelr)] = probct*(hampelr-abs(r[which(r > hampelb & r <= hampelr)]))/
      (hampelr -hampelb)*1/abs(r[which(r > hampelb & r <= hampelr)])
    wye[which(r > hampelr)] = 0
    wye = as.numeric(wye)
    
    probct = qchisq(probp1,a)
    hampelb = qchisq(hampelp2, a)
    hampelr = qchisq(hampelp3, a)
    wte = wtn
    wte[which(wtn <= probct)] = 1 
    wte[which(wtn > probct & wtn <= hampelb)] = probct/abs(wtn[which(wtn > probct & wtn <= hampelb)])
    wte[which(wtn > hampelb & wtn <= hampelr)] = probct*(hampelr-abs(wtn[which(wtn > hampelb & wtn <= hampelr)]))/
      (hampelr -hampelb)*1/abs(wtn[which(wtn > hampelb & wtn <= hampelr)])
    wte[which(wtn > hampelr)] = 0
    
    
    difference = abs(sum(b^2) - rold)/rold
    rold = sum(b^2)
    we = wye * wte
    if(any(we<1e-6)){
      w0 = which(we<1e-6)
      we = replace(we,list=w0,values=1e-6)
      zerows = unique(c(zerows,as.numeric(names(w0))))
    }
    
    if(length(zerows)>=(n/2)){
      break
    }
    dataw = as.data.frame(datamc * sqrt(we))
    loops = loops + 1
  }
  if (difference > prec){
    warning(paste("Method did not converge. The scaled difference between norms of the coefficient vectors is ", round(difference, digits=4)))
  }
  
  w = we
  w[zerows] = 0 
  wt = wte
  wt[zerows] = 0 
  wy = wye
  wy[zerows] = 0 
  
  P = res.nipls$loadings
  W = res.nipls$W
  R = res.nipls$R
  Tpls = scale(datam[,2:qs],center=datac[2:qs],scale=datas[2:qs]) %*% R 
  Tplst = scale(datamt, center = datact, scale = datast) %*% R
  
  output = list(scores = Tpls, scoresT = Tplst, R=R, loadings = P)
  return(output)
}


#_________________________________________________________________________________________________________________________________________________________________________
#______________________________________________________________________________FPLS_function______________________________________________________________________________
#_________________________________________________________________________________________________________________________________________________________________________
fpls_fun = function(Y, X_tr, X_te, nbf, npls, rangeval, method = c("classical", "robust")){
  
  getMat = getAmat(X_tr = X_tr, X_te = X_te, nbf = nbf, rangeval = rangeval)
  X_tr = getMat$Amat_tr
  X_te = getMat$Amat_te
  Bsp_funs = getMat$Bsp_funs
  InnProds = getMat$InnProds
  
  np = length(X_tr)
  N = dim(X_tr[[1]])[1]
  M = numeric()
  for(im in 1:np){
    M[im] = dim(X_tr[[im]])[2]
  }
  
  comps = list()
  comps_test = list()
  base = list()
  
  for(jm in 1:np){
    X = X_tr[[jm]]
    Xt = X_te[[jm]]
    data = as.data.frame(X)
    data$y = Y
    
    if(method == "classical"){
      it_model = simpls(X,Xt,Y,npls)
      comps[[jm]] = it_model$Ts
      comps_test[[jm]] = it_model$Tst
      base[[jm]] = it_model$R
    }else if(method == "robust"){
      it_model = prms2(y~., data = data, Xt, a = npls, fun = "Hampel")
      comps[[jm]] = it_model$scores
      comps_test[[jm]] = it_model$scoresT
      base[[jm]] = it_model$R
    }
  }
  
  if(method == "classical"){
    fmodel = lm(Y~do.call(cbind, comps))
    intercept = fmodel$coefficients[1]
    coefs = fmodel$coefficients[-1]
  }else if(method == "robust"){
    fit_model = lm(Y~do.call(cbind, comps))
    fit_model_rb = rlm(Y~do.call(cbind, comps), psi = psi.bisquare, c = 4.685, maxit = 1000)
    fit_model_rlmdd = rlmDD(Y, do.call(cbind, comps), fit_model$coef, fit_model_rb$coef, method = "Bisquare",plot = FALSE)
    fit_coefs = fit_model_rlmdd$esti$coefficients
    intercept = fit_coefs[1]
    coefs = fit_coefs[-1]
  }
  
  fitted_values = do.call(cbind, comps) %*% coefs + intercept
  pred_vals = do.call(cbind, comps_test) %*% coefs + intercept
  
  return(list(fits = fitted_values, preds = pred_vals, coefs = coefs, plsbase = base, plsscore = comps,
              Bsp_funs = Bsp_funs, InnProds = InnProds))
}

#_________________________________________________________________________________________________________________________________________________________________________
#_________________________________________________________________________________________________________________________________________________________________________

#___________________________________________________________________________________________________
#_______________________________Estimated_coefficient_functions_(PLS)_______________________________
#___________________________________________________________________________________________________

getCoef = function(object){
  
  Bsp_funs = object$Bsp_funs
  Inn_prd_sqrts = object$InnProds
  plsBase = object$plsbase
  coefs = object$coefs
  np = length(Bsp_funs)
  npls = dim(plsBase[[1]])[2]
  
  Brest = matrix(coefs, ncol = npls, byrow = TRUE)
  
  Cf = list()
  for(i in 1:np)
    Cf[[i]] = t(solve(Inn_prd_sqrts[[i]])) %*% plsBase[[i]] %*% Brest[i,]
  
  Bhat = list()
  for(i in 1:np)
    Bhat[[i]] = t(Cf[[i]]) %*% t(Bsp_funs[[i]])
  
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
  return(list(PCAcoef = PCAcoef, PCAscore = PCAscore,
              PCAscore_test = PCAscore_test, evalbase = evalbase))
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
var_select = function(Y, X, nbf, ncp, rangeval, method, fmethod){
  np = length(X)
  
  if(fmethod == "fpca"){
    Dmat = list()
    for(i in 1:np)
      Dmat[[i]] = getPCA(X_tr = X[[i]], X_te = X[[i]], nbasis = nbf[[i]], ncomp = ncp[[i]], rangeval = rangeval[[i]])$PCAscore
    Dmat = do.call(cbind, Dmat)
  }
  if(fmethod == "fpls"){
    Dmat = do.call(cbind,
                   fpls_fun(Y = Y, X_tr = X, X_te = X, nbf = nbf, npls = ncp[1], rangeval = rangeval, method = method)$plsscore)
  }
  rsm_model = regRSM(y = Y, x = Dmat, B = 1000)
  empt_mat = rep(NA, sum(ncp))
  empt_mat[sort(rsm_model$model)] = sort(rsm_model$model)
  mat_list = list()
  mat_list[[1]] = empt_mat[1:ncp[1]]
  ik = ncp[1]
  for(i in 2:np){
    mat_list[[i]] = empt_mat[(ik+1):(ik+ncp[i])]
    ik = ik + ncp[i]
  }
  mat_ind = matrix(NA, ncol = 1, nrow = np)
  for(i in 1:np)
    mat_ind[i,] = sum(mat_list[[i]], na.rm = TRUE)
  variable_indices = which(rowSums(mat_ind, na.rm = T) > 0)
  
  return(variable_indices)
}
#__________________________________________________________________________________________
#__________________________________________________________________________________________


#_____________________________________________________
#_______________________Trimmed_MSE___________________
#_____________________________________________________
tmse = function(Y, fit, tp){
  
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

#___________________________________________________________________________________________
#___________________________________K-fold_cross-validation_________________________________
#___________________________________________________________________________________________
nfoldCV = function(Y, X, ncomp, nfold, nbf, rangeval, method = c("classical", "robust"),
                   fmethod = c("fpca", "fpls")){
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
    
    preds_k = list()
    for(k in 1:ncomp){
      if(fmethod == "fpca"){
        fold_res = fpca(Y = nfY_tr, X_tr = nfX_tr, X_te = nfX_te,
                        npca=k, nbf = nbf, rangeval = rangeval)
      }else if(fmethod == "fpls"){
        fold_res = fpls_fun(Y = nfY_tr, X_tr = nfX_tr, X_te = nfX_te, nbf = nbf,
                            npls = k, rangeval = rangeval, method = method)
        
      }
      preds_k[[k]] = fold_res$preds
    }
    
    
    
    
    
    mspe_i = numeric()
    for(k in 1: ncomp)
      mspe_i[k] = tmse(nfY_te, preds_k[[k]])
    
    nfMspe[i] = which.min(mspe_i)
    
  }
  return(round(mean(nfMspe)))
}
#___________________________________________________________________________________________
#___________________________________________________________________________________________


#___________________________________________________________________________________________
#________________________________________REISEE_function____________________________________
#___________________________________________________________________________________________
reisee = function(beta, beta_hat, domain){
  f = (beta - beta_hat)^2
  ncol.f = dim(f)[2]
  gap.mat = diff(domain)
  r1 = (sum(f[-1]*gap.mat) + sum(f[-ncol.f]*gap.mat))/2
  
  f1 = beta^2
  ncol.f1 = dim(f1)[2]
  r2 = (sum(f1[-1]*gap.mat) + sum(f1[-ncol.f1]*gap.mat))/2
  
  s1 = sum(r1*diff(domain))/2
  s2 = sum(r2*diff(domain))/2
  
  return(s1/s2)
}
#___________________________________________________________________________________________
#___________________________________________________________________________________________
