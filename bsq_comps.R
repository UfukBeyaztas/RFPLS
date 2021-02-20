library(rlmDataDriven)
library(MASS)
options(warn = -1)

bsq_base = function(Y, X, npls){
  
  s_fun = function(data){
    return((data - mean(data)) / sd(data))
  }
  
  np = length(X)
  N = length(Y)
  M = ncol(X[[1]])
  
  Xs = list()
  for(i in 1:np)
    Xs[[i]] = scale(X[[i]], center = TRUE, scale = FALSE)
  
  for(ii in 1:npls){
    
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
      for(j in 1:np)
        comps[[j]] = Xs[[j]] %*% it_coefs[,j]
    }else{
      for(j in 1:np)
        comps[[j]] = cbind(comps[[j]], Xs[[j]] %*% it_coefs[,j])
    }
    
    for(j in 1:np){
      for(jj in 1:M){
        Xs[[j]][,jj] = Xs[[j]][,jj] - mean(Xs[[j]][,jj]) - as.numeric(cov(Xs[[j]][,jj], comps[[j]][,ii]) / 
                                                                        var(comps[[j]][,ii])) * (comps[[j]][,ii] - mean(comps[[j]][,ii]))
      }
    }
  }
  
  return(comps)
}
