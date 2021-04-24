rm(list=ls())
source("auxilary_functions.R")
source("main_function.R")
source("data_generation.R")
source("data_generation2.R")

# Generate the data
sim_data = data_generation(n=400, j=200, index=c(1,2,3), ntrain=200)
sim_data2 = data_generation2(n=400, j=200, index=c(1,2,3), ntrain=200)

Y_train = sim_data$Y_tr # generated scalar response (training sample)
Y_test = sim_data$Y_te # generated scalar response (test sample)
X_train = sim_data$X_tr # generated functional predictors (training sample)
X_test = sim_data$X_te # generated functional predictors (test sample)

# Outliers
nout = 0.2 * length(Y_train)
out_indx = sample(1:length(Y_train), nout, replace=FALSE)

for(iij in 1:3){
  X_train[[iij]][out_indx,] = sim_data2$X_tr[[iij]][out_indx,] 
}

Y_train[out_indx] = sim_data2$Y_tr[out_indx]

nbf = rep(20, 5)

# Note that k-fold cross-validation increases the computational cost
# A specific number of pls component can be used with npls_use
# In the followings, use 2-fold cross-validation

# FPLS (Full model)
pls_full = main_fun(Y = Y_train, X_tr = X_train, X_te = X_test, nfold = 2, npls_max = 5, nbf = nbf,
                    method = "classical", model = "full", CV = TRUE)

# RFPLS (Full model)
rpls_full = main_fun(Y = Y_train, X_tr = X_train, X_te = X_test, nfold = 2, npls_max = 5, nbf = nbf,
                    method = "robust", model = "full", CV = TRUE)

# FPLS (True model)
pls_true = main_fun(Y = Y_train, X_tr = X_train, X_te = X_test, nfold = 2, npls_max = 5, nbf = nbf,
                    method = "classical", model = "true", CV = TRUE)

# RFPLS (True model)
rpls_true = main_fun(Y = Y_train, X_tr = X_train, X_te = X_test, nfold = 2, npls_max = 5, nbf = nbf,
                     method = "robust", model = "true", CV = TRUE)

# FPLS (Selected model)
pls_selected = main_fun(Y = Y_train, X_tr = X_train, X_te = X_test, nfold = 2, npls_max = 5, nbf = nbf,
                     method = "classical", model = "selected", CV = TRUE)

# RFPLS (Selected model)
rpls_selected = main_fun(Y = Y_train, X_tr = X_train, X_te = X_test, nfold = 2, npls_max = 5, nbf = nbf,
                     method = "robust", model = "selected", CV = TRUE)

# Selected variables
pls_selected$var_used
rpls_selected$var_used


# MSPE
# Full model
mean((Y_test - pls_full$predictions)^2)
mean((Y_test - rpls_full$predictions)^2)
# True model
mean((Y_test - pls_true$predictions)^2)
mean((Y_test - rpls_true$predictions)^2)
# Selected model
mean((Y_test - pls_selected$predictions)^2)
mean((Y_test - rpls_selected$predictions)^2)


# R^2
# Full model
cor(Y_test, pls_full$predictions)^2
cor(Y_test, rpls_full$predictions)^2

# True model
cor(Y_test, pls_true$predictions)^2
cor(Y_test, rpls_true$predictions)^2

# Selected model
cor(Y_test, pls_selected$predictions)^2
cor(Y_test, rpls_selected$predictions)^2

# IMSE
imse_pls1 = numeric()
imse_rpls1 = numeric()

for(ims in 1:3){
  imse_pls1[ims] = mean((c(sim_data$tcoefs[[ims]]$data) - c(pls_true$coefficients[[ims]]))^2)
  imse_rpls1[ims] = mean((c(sim_data$tcoefs[[ims]]$data) - c(rpls_true$coefficients[[ims]]))^2)
}

mean(imse_pls1)
mean(imse_rpls1)
