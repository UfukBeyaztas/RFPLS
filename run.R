source("bsq_comps.R") # for variable selection procedure
source("bsq_funs.R") # for the RFPLS method
source("rsm.R") # functions to be used in the RSM procedure
source("data_generation.R") # daga generation

# Generate the data
sim_data = data_generation(n=300, j=100, index=c(1,2,3), ntrain=100)

Y_train = sim_data$Y_tr # generated scalar response (training sample)
Y_test = sim_data$Y_te # generated scalar response (test sample)
X_train = sim_data$X_tr # generated functional predictors (training sample)
X_test = sim_data$X_te # generated functional predictors (test sample)

# Variable selection
npls = 2
vs_model = bsq_base(Y = Y_train, X = X_train, npls = npls)
projections = do.call(cbind, vs_model)
rsm_model = regRSM(y = Y_train, x = projections, B = 1000)
empt_mat = rep(NA, npls*length(X_train))
empt_mat[sort(rsm_model$model)] = sort(rsm_model$model)
empt_mat = matrix(empt_mat, nrow = length(X_train), ncol = npls, byrow = TRUE)
variable_indices = which(rowSums(empt_mat, na.rm = T) > 0) # selected variables

# Selected functional predictors
Xs_train = list()
Xs_test = list()
for(i in 1:length(variable_indices)){
  Xs_train[[i]] = X_train[[variable_indices[i]]]
  Xs_test[[i]] = X_test[[variable_indices[i]]]
}

# Final model
fmodel = bsq_fun(Y = Y_train, X = Xs_train, X_test = Xs_test, npls = 5)

# MSPE
mean((Y_test - fmodel$preds)^2)

