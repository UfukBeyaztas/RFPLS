library(fda.usc)

data_generation2 = function(n, j, ntrain){
  s = seq(0, 1, length.out = j)
  
  fX = list()
  for(ij in 1:3){
    
    ksi = list()
    for(ik in 1:5){
      ksi[[ik]] = rnorm(n, 0, sd = (4*ik^(-3/2)))
    }
    
    phi = list()
    for(ik in 1:5){
      phi[[ik]] = 2*sin(ik * pi * s) - cos(ik * pi * s)
    }
    
    fX[[ij]] = Reduce("+", lapply(1:5, function(k){ksi[[k]] %*% t(phi[[k]])}))
  }

  vBeta = list()
  vBeta[[1]] = sin(2*pi* s)
  vBeta[[2]] = sin(2*pi* s)
  vBeta[[3]] = cos(2*pi * s)

  for(ij in 1:3){
    fX[[ij]] = fdata(fX[[ij]], argvals = s)
    vBeta[[ij]] = fdata(vBeta[[ij]], argvals = s)
  }
  
  err = rnorm(n, mean=10, sd=1)
  
  fY = Reduce("+", lapply(1:length(fX), function(k){inprod.fdata(fX[[k]], vBeta[[k]])})) 
  fYe = fY + err
  
  X_train = list()
  X_test = list()
  
  for(ij in 1:3){
    X_train[[ij]] = fX[[ij]]$data[1:ntrain,]
    X_test[[ij]] = fX[[ij]]$data[-(1:ntrain),]
  }
  
  return(list("Y_tr" = fYe[1:ntrain,], "Y_te" = fYe[-(1:ntrain),], "X_tr" = X_train, "X_te" = X_test, tcoefs = vBeta))
  
}
