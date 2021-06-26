library(fda.usc)

data_generation2 = function(n, j, index, ntrain){
  s = seq(0, 1, length.out = j)
  
  gammaexpcov = function(n.point, t){ 
    Sigma = matrix(0, n.point,n.point)
    
    for(i in 1:(n.point-1)){
      for(j in (i+1):n.point){
        Sigma[i,j] = exp(-(10*abs(t[j]-t[i]))^2)
      }
    }
    
    Sigma = (Sigma+t(Sigma))
    diag(Sigma) = 1
    return(Sigma)
  }
  
  Sigma.X = gammaexpcov(j, s)
  
  eig = eigen(Sigma.X)
  eigval = eig$values
  b = cumsum(eigval)/sum(eigval)
  n.comp.sigma = which(b>0.999999999999)[1]
  eig0 = eigval[1:n.comp.sigma]
  sigma.5 =  eig$vectors[,1:n.comp.sigma] %*% sqrt(diag(eig0))
  range(abs(sigma.5 %*% t(sigma.5) - Sigma.X))
  
  V=list()
  
  for(i in 1:2){
    z = matrix(rnorm(n.comp.sigma*n, 0,1), n.comp.sigma, n)
    V[[i]] = t(sigma.5 %*% z )
  }
  
  
  fX = list()
  for(ij in 1:3){
    
    ksi = list()
    for(ik in 1:5){
      ksi[[ik]] = rnorm(n, 0, sd = (4*ik^(-3/2)))
    }
    
    phi = list()
    for(ik in 1:5){
      phi[[ik]] = sin(ik * pi * s) - cos(ik * pi * s)
    }
    
    fX[[ij]] = Reduce("+", lapply(1:5, function(k){ksi[[k]] %*% t(phi[[k]])}))
  }
  
  fX[[4]] = V[[1]]+10
  fX[[5]] = V[[2]]+10
  
  vBeta = list()
  vBeta[[1]] = 4*sqrt(s)
  vBeta[[2]] = 2*exp(-(s - 0.5)^2)
  vBeta[[3]] = 3*cos(pi * s)
  vBeta[[4]] = exp(-(s^2))
  vBeta[[5]] = 2 * sqrt(s)
  
  for(ij in 1:5){
    fX[[ij]] = fdata(fX[[ij]], argvals = s)
    vBeta[[ij]] = fdata(vBeta[[ij]], argvals = s)
  }
  
  err = rnorm(n, mean=0, sd=1)
  
  fY = Reduce("+", lapply(1:length(index), function(k){inprod.fdata(fX[[index[k]]], vBeta[[index[k]]])})) 
  fYe = fY + err
  
  X_train = list()
  X_test = list()
  
  for(ij in 1:5){
    X_train[[ij]] = fX[[ij]]$data[1:ntrain,]
    X_test[[ij]] = fX[[ij]]$data[-(1:ntrain),]
  }
  
  return(list("Y_tr" = fYe[1:ntrain,], "Y_te" = fYe[-(1:ntrain),], "X_tr" = X_train, "X_te" = X_test, tcoefs = vBeta))
  
}
