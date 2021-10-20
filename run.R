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
    sim_data = data_generation(n=400, j=200, index=c(1,2,3), ntrain=200)
    sim_data2 = data_generation2(n=400, j=200, index=c(1,2,3), ntrain=200)
    
    Y_train = sim_data$Y_tr # generated scalar response (training sample)
    Y_test = sim_data$Y_te # generated scalar response (test sample)
    X_train = sim_data$X_tr # generated functional predictors (training sample)
    X_test = sim_data$X_te # generated functional predictors (test sample)
    
    out_indx = sample(1:length(Y_train), nout, replace=FALSE) # randomly select the outlier indices
    
    for(iij in 1:3){
      X_train[[iij]][out_indx,] = sim_data2$X_tr[[iij]][out_indx,] # contaminate predictors with the leverage points
    }
    
    Y_train[out_indx] = sim_data2$Y_tr[out_indx] # conraminate response with the vertical outliers
    
    # FPCA (Full model)
    pca_full = main_fun(Y = Y_train, X_tr = X_train, X_te = X_test, npca_max = 5, nbf = rep(20,5), nfold = nfold,
                        rangeval = rangeval, method = c("classical"), model = c("full"), fmethod = c("fpca"))
    
    # FPLS (Full model)
    pls_full = main_fun(Y = Y_train, X_tr = X_train, X_te = X_test, npls_max = 5, nbf = rep(20,5), nfold = nfold,
                        rangeval = rangeval, method = c("classical"), model = c("full"), fmethod = c("fpls"))
    
    # RFPLS (Full model)
    rpls_full = main_fun(Y = Y_train, X_tr = X_train, X_te = X_test, npls_max = 5, nbf = rep(20,5), nfold = nfold,
                         rangeval = rangeval, method = c("robust"), model = c("full"), fmethod = c("fpls"))
    
    
    # FPCA (True model)
    pca_true = main_fun(Y = Y_train, X_tr = X_train, X_te = X_test, t_index = c(1,2,3), npca_max = 5, nbf = rep(20,5), nfold = nfold,
                        rangeval = rangeval[c(1,2,3)], method = c("classical"), model = c("true"), fmethod = c("fpca"))
    
    # FPLS (True model)
    pls_true = main_fun(Y = Y_train, X_tr = X_train, X_te = X_test, t_index = c(1,2,3), npls_max = 5, nfold = nfold,
                        rangeval = rangeval[c(1,2,3)], nbf = rep(20,3), method = c("classical"), model = c("true"), fmethod = c("fpls"))
    
    # RFPLS (True model)
    rpls_true = main_fun(Y = Y_train, X_tr = X_train, X_te = X_test, t_index = c(1,2,3), npls_max = 5, nfold = nfold,
                         rangeval = rangeval[c(1,2,3)], nbf = rep(20,3), method = c("robust"), model = c("true"), fmethod = c("fpls"))
    
    
    # FPCA (Selected model)
    pca_selected = main_fun(Y = Y_train, X_tr = X_train, X_te = X_test, npca_max = 5, nbf = rep(20,5), nfold = nfold,
                            rangeval = rangeval, ncp_vs = rep(5,5), method = c("classical"), model = c("selected"), fmethod = c("fpca"))
    
    # FPLS (Selected model)
    pls_selected = main_fun(Y = Y_train, X_tr = X_train, X_te = X_test, npls_max = 5, nbf = rep(20,5), nfold = nfold,
                            rangeval = rangeval, ncp_vs = rep(5,5), method = c("classical"), model = c("selected"), fmethod = c("fpls"))
    
    # RFPLS (Selected model)
    rpls_selected= main_fun(Y = Y_train, X_tr = X_train, X_te = X_test, npls_max = 5, nbf = rep(20,5), nfold = nfold,
                            rangeval = rangeval, ncp_vs = rep(5,5), method = c("robust"), model = c("selected"), fmethod = c("fpls"))
    
    
    # MSPE
    # Full model
    tmse(Y_test, pca_full$predictions) # 3.430534
    tmse(Y_test, pls_full$predictions) # 2.125717
    tmse(Y_test, rpls_full$predictions) # 0.7990948
    
    # True model
    tmse(Y_test, pca_true$predictions) # 2.090089
    tmse(Y_test, pls_true$predictions) # 1.703035
    tmse(Y_test, rpls_true$predictions) # 0.7014038
    
    # Selected model
    tmse(Y_test, pca_selected$predictions) # 2.090089
    tmse(Y_test, pls_selected$predictions) # 1.974743
    tmse(Y_test, rpls_selected$predictions) # 0.7014038
    
    
    # R^2
    # Full model
    trsq(Y_test, pca_full$predictions) # 0.06581048
    trsq(Y_test, pls_full$predictions) # 0.3615189
    trsq(Y_test, rpls_full$predictions) # 0.6249167
    
    # True model
    trsq(Y_test, pca_true$predictions) # 0.3190196
    trsq(Y_test, pls_true$predictions) # 0.4898033
    trsq(Y_test, rpls_true$predictions) # 0.6510442
    
    # Selected model
    trsq(Y_test, pca_selected$predictions) # 0.3190196
    trsq(Y_test, pls_selected$predictions) # 0.4188317
    trsq(Y_test, rpls_selected$predictions) # 0.6510442
    
    

    # RISEE
    # FPCA
    # \beta_1(t)
    reisee(sim_data$tcoefs[[1]]$data, pca_true$coefficients[[1]], seq(0,1, length.out=200)) # 12.23619
    # \beta_2(t)
    reisee(sim_data$tcoefs[[2]]$data, pca_true$coefficients[[2]], seq(0,1, length.out=200)) # 0.5466141
    # \beta_3(t)
    reisee(sim_data$tcoefs[[3]]$data, pca_true$coefficients[[3]], seq(0,1, length.out=200)) # 0.8148587
    
    # FPLS
    # \beta_1(t)
    reisee(sim_data$tcoefs[[1]]$data, pls_true$coefficients[[1]], seq(0,1, length.out=200)) # 15.8621
    # \beta_2(t)
    reisee(sim_data$tcoefs[[2]]$data, pls_true$coefficients[[2]], seq(0,1, length.out=200)) # 17.98122
    # \beta_3(t)
    reisee(sim_data$tcoefs[[3]]$data, pls_true$coefficients[[3]], seq(0,1, length.out=200)) # 19.66167
    
    # RFPLS
    # \beta_1(t)
    reisee(sim_data$tcoefs[[1]]$data, rpls_true$coefficients[[1]], seq(0,1, length.out=200)) # 0.3455245
    # \beta_2(t)
    reisee(sim_data$tcoefs[[2]]$data, rpls_true$coefficients[[2]], seq(0,1, length.out=200)) # 0.373456
    # \beta_3(t)
    reisee(sim_data$tcoefs[[3]]$data, rpls_true$coefficients[[3]], seq(0,1, length.out=200)) # 0.5194473
    
    
