rm(list=ls())
source("auxilary_functions.R")
source("main_function.R")
source("data_generation.R")
source("data_generation2.R")

nout = 20 # number of outliers (20% percentage)

rangeval = list(c(0,1), # rangeval values for each of the functional predictors
                c(0,1),
                c(0,1),
                c(0,1),
                c(0,1))

# here only two fold cross-validation is used
nfold = 2


set.seed(1234) #
# Generate the data
sim_data = data_generation(n=400, j=200, ntrain=200)
sim_data2 = data_generation2(n=400, j=200, ntrain=200)

Y_train = sim_data$Y_tr # generated scalar response (training sample)
Y_test = sim_data$Y_te # generated scalar response (test sample)
X_train = sim_data$X_tr # generated functional predictors (training sample)
X_test = sim_data$X_te # generated functional predictors (test sample)

out_indx = sample(1:length(Y_train), nout, replace=FALSE) # randomly select the outlier indices

for(iij in 1:3){
  X_train[[iij]][out_indx,] = sim_data2$X_tr[[iij]][out_indx,] # contaminate predictors with the leverage points
}

Y_train[out_indx] = sim_data2$Y_tr[out_indx] # conraminate response with the vertical outliers

# FPCA
pca_full = main_fun(Y = Y_train, X_tr = X_train, X_te = X_test, npca_max = 5, nbf = rep(20,5), nfold = nfold,
                    rangeval = rangeval, method = c("classical"), fmethod = c("fpca"))

# FPLS
pls_full = main_fun(Y = Y_train, X_tr = X_train, X_te = X_test, npls_max = 5, nbf = rep(20,5), nfold = nfold,
                    rangeval = rangeval, method = c("classical"), fmethod = c("fpls"))

# RFPLS
rpls_full = main_fun(Y = Y_train, X_tr = X_train, X_te = X_test, npls_max = 5, nbf = rep(20,5), nfold = nfold,
                     rangeval = rangeval, method = c("robust"), fmethod = c("fpls"))


# MSPE
tmse(Y_test, pca_full$predictions) # 1.709957
tmse(Y_test, pls_full$predictions) # 1.436411
tmse(Y_test, rpls_full$predictions) # 0.6061528


# R^2
trsq(Y_test, pca_full$predictions) # 0.8728854
trsq(Y_test, pls_full$predictions) # 0.8972966
trsq(Y_test, rpls_full$predictions) # 0.9338817

# RISEE
# FPCA
# \beta_1(t)
reisee(sim_data$tcoefs[[1]]$data, pca_full$coefficients[[1]], seq(0,1, length.out=200)) # 1.870509
# \beta_2(t)
reisee(sim_data$tcoefs[[2]]$data, pca_full$coefficients[[2]], seq(0,1, length.out=200)) # 1.991387
# \beta_3(t)
reisee(sim_data$tcoefs[[3]]$data, pca_full$coefficients[[3]], seq(0,1, length.out=200)) # 2.147665

# FPLS
# \beta_1(t)
reisee(sim_data$tcoefs[[1]]$data, pls_full$coefficients[[1]], seq(0,1, length.out=200)) # 1.386727
# \beta_2(t)
reisee(sim_data$tcoefs[[2]]$data, pls_full$coefficients[[2]], seq(0,1, length.out=200)) # 4.330007
# \beta_3(t)
reisee(sim_data$tcoefs[[3]]$data, pls_full$coefficients[[3]], seq(0,1, length.out=200)) # 4.937152

# RFPLS
# \beta_1(t)
reisee(sim_data$tcoefs[[1]]$data, rpls_full$coefficients[[1]], seq(0,1, length.out=200)) # 1.437968
# \beta_2(t)
reisee(sim_data$tcoefs[[2]]$data, rpls_full$coefficients[[2]], seq(0,1, length.out=200)) # 0.2259788
# \beta_3(t)
reisee(sim_data$tcoefs[[3]]$data, rpls_full$coefficients[[3]], seq(0,1, length.out=200)) # 0.4031647
