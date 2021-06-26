rm(list=ls())
source("auxilary_functions.R")
source("main_function.R")
source("data_generation.R")
source("data_generation2.R")

rangeval = list(c(0,1),
                c(0,1),
                c(0,1),
                c(0,1),
                c(0,1))

# Generate the data
set.seed(123456)
sim_data = data_generation(n=400, j=200, index=c(1,2,3), ntrain=200)
sim_data2 = data_generation2(n=400, j=200, index=c(1,2,3), ntrain=200)

Y_train = sim_data$Y_tr # generated scalar response (training sample)
Y_test = sim_data$Y_te # generated scalar response (test sample)
X_train = sim_data$X_tr # generated functional predictors (training sample)
X_test = sim_data$X_te # generated functional predictors (test sample)

# In the following, the generated data are contaminated by outliers.
# To evaluate the methods when no outliers present in the data, ignore lines 20-27
# Outliers
nout = 20 # 10%
out_indx = sample(1:length(Y_train), nout, replace=FALSE)

for(iij in 1:3){
  X_train[[iij]][out_indx,] = sim_data2$X_tr[[iij]][out_indx,] 
}

Y_train[out_indx] = sim_data2$Y_tr[out_indx]


# FPLS (Full model)
pls_full = main_fun(Y = Y_train, X_tr = X_train, X_te = X_test, npls_max = 5, nbf = rep(20,5),
                    rangeval = rangeval, method = c("classical"), model = c("full"), fmethod = c("fpls"))

# RFPLS (Full model)
rpls_full = main_fun(Y = Y_train, X_tr = X_train, X_te = X_test, npls_max = 5, nbf = rep(20,5), 
                     rangeval = rangeval, method = c("robust"), model = c("full"), fmethod = c("fpls"))

# FPCA (Full model)
pca_full = main_fun(Y = Y_train, X_tr = X_train, X_te = X_test, npca_max = 5, nbf = rep(20,5), 
                    rangeval = rangeval, method = c("classical"), model = c("full"), fmethod = c("fpca"))

# FPLS (True model)
pls_true = main_fun(Y = Y_train, X_tr = X_train, X_te = X_test, t_index = c(1,2,3), npls_max = 5,
                    rangeval = rangeval[c(1,2,3)], nbf = rep(20,3), method = c("classical"), model = c("true"), fmethod = c("fpls"))

# RFPLS (True model)
rpls_true = main_fun(Y = Y_train, X_tr = X_train, X_te = X_test, t_index = c(1,2,3), npls_max = 5,
                     rangeval = rangeval[c(1,2,3)], nbf = rep(20,3), method = c("robust"), model = c("true"), fmethod = c("fpls"))

# FPCA (True model)
pca_true = main_fun(Y = Y_train, X_tr = X_train, X_te = X_test, t_index = c(1,2,3), npca_max = 5, nbf = rep(20,3),
                    rangeval = rangeval[c(1,2,3)], method = c("classical"), model = c("true"), fmethod = c("fpca"))

# FPLS (Selected model)
pls_selected = main_fun(Y = Y_train, X_tr = X_train, X_te = X_test, npls_max = 5, nbf = rep(20,5), 
                        rangeval = rangeval, nbf_vs = rep(5,5), method = c("classical"), model = c("selected"), fmethod = c("fpls"))

# RFPLS (Selected model)
rpls_selected= main_fun(Y = Y_train, X_tr = X_train, X_te = X_test, npls_max = 5, nbf = rep(20,5),
                        rangeval = rangeval, nbf_vs = rep(5,5), method = c("robust"), model = c("selected"), fmethod = c("fpls"))

# FPCA (Selected model)
pca_selected = main_fun(Y = Y_train, X_tr = X_train, X_te = X_test, npca_max = 5, nbf = rep(20,5),
                        rangeval = rangeval, nbf_vs = rep(5,5), method = c("classical"), model = c("selected"), fmethod = c("fpca"))



# MSE
# Full model
tmse(Y_train, pca_full$fitvals) # 2.360425
tmse(Y_train, pls_full$fitvals) # 2.285565
tmse(Y_train, rpls_full$fitvals) # 0.8439575

# True model
tmse(Y_train, pca_true$fitvals) # 2.214008
tmse(Y_train, pls_true$fitvals) # 2.231887
tmse(Y_train, rpls_true$fitvals) # 0.8714299

# Selected model
tmse(Y_train, pca_selected$fitvals) # 2.214008
tmse(Y_train, pls_selected$fitvals) # 2.231887
tmse(Y_train, rpls_selected$fitvals) # 0.8714299


# MSPE
# Full model
tmse(Y_test, pca_full$predictions) # 2.279117
tmse(Y_test, pls_full$predictions) # 2.251953
tmse(Y_test, rpls_full$predictions) # 1.026346

# True model
tmse(Y_test, pca_true$predictions) # 1.635951
tmse(Y_test, pls_true$predictions) # 2.069308
tmse(Y_test, rpls_true$predictions) # 0.9144857

# Selected model
tmse(Y_test, pca_selected$predictions) # 1.635951
tmse(Y_test, pls_selected$predictions) # 2.069308
tmse(Y_test, rpls_selected$predictions) # 0.9144857


# R^2
# Full model
trsq(Y_train, pca_full$fitvals) # 0.9407093
trsq(Y_train, pls_full$fitvals) # 0.9462906
trsq(Y_train, rpls_full$fitvals) # 0.9734357

# True model
trsq(Y_train, pca_true$fitvals) # 0.9445404
trsq(Y_train, pls_true$fitvals) # 0.9460039
trsq(Y_train, rpls_true$fitvals) # 0.9728872

# Selected model
trsq(Y_train, pca_selected$fitvals) # 0.9445404
trsq(Y_train, pls_selected$fitvals) # 0.9460039
trsq(Y_train, rpls_selected$fitvals) # 0.9728872


# R_p^2
# Full model
trsq(Y_test, pca_full$predictions) # 0.9440397
trsq(Y_test, pls_full$predictions) # 0.9478304
trsq(Y_test, rpls_full$predictions) # 0.9732456

# True model
trsq(Y_test, pca_true$predictions) # 0.9597614
trsq(Y_test, pls_true$predictions) # 0.9516234
trsq(Y_test, rpls_true$predictions) # 0.9748124

# Selected model
trsq(Y_test, pca_selected$predictions) # 0.9597614
trsq(Y_test, pls_selected$predictions) # 0.9516234
trsq(Y_test, rpls_selected$predictions) # 0.9748124



# IMSE
imse_pca = numeric()
imse_pls = numeric()
imse_rpls = numeric()

for(ims in 1:3){
  imse_pca[ims] = mean((c(sim_data$tcoefs[[ims]]$data) - c(pca_true$coefficients[[ims]]))^2)
  imse_pls[ims] = mean((c(sim_data$tcoefs[[ims]]$data) - c(pls_true$coefficients[[ims]]))^2)
  imse_rpls[ims] = mean((c(sim_data$tcoefs[[ims]]$data) - c(rpls_true$coefficients[[ims]]))^2)
}

imse_pca # 1.1614318 0.4575989 0.6843156
imse_pls # 1.7204460 0.4402417 0.7041253
imse_rpls # 0.1985453 0.5832299 0.5827725

