# Triple step algorithm to find optimal model
#
# Inputs: 
# 'Ys' : list with i entries of matrices containing the Yx and Yr
# 'YsEval' : list with i entries of matrices containing the YxEval and YrEval
#
# Output: 
# All_store: list with 10 entries 
# All_store[[1]] : list containing all optimized model for each combination of lambdas of the 2 sequences
# All_store[[2]] : matrix containing the iterations needed to optimize each model. If an optimized model reached the maximum of iterations, one may consider rerunning the estimation with different starting values
# All_store[[3]] : matrix storing the column indices of the columns found active in step 2 of the optimization 
# All_store[[4]] : the regularization parameters used for estimation
# All_store[[5]] : data used for model estimation
# All_store[[6]] : data used for model evaluation
# All_store[[7]] : indicator of optimal model
# All_store[[8]] : optimal lambda combination
# All_store[[9]] : optimal model
# All_store[[10]] : mean squared error of the models on the evaluation data

# Estimation procedure
# The code iterates over the lambda combinations and prints to the console the 
# current lambda combination it is deriving and the iteration it is at 
# in the following form: 'lambda combination' | 'iteration'
# This information is separately provided for each of the 3 steps of the algorithm

rm(list=ls(all=TRUE))
source("BetaUpdate_SGNAR.R")

# choose number of groups per continent. Can be either 3 or 10
Groups = 3 

if (Groups == 3) {
  # for estimation of 3 groups per continent
  load("BlockchainData3Groups.RData")
} else if (Groups == 10) {
  # load for estimation of 10 groups per continent
  load("BlockchainData10Groups.RData")
}

# year to estimate the adjacency matrix for
# Can be either '2012', '2013', '2014', '2015', '2016'
which_year = "2012"
ite = switch(which_year, "2012" = 1, "2013" = 2, "2014" = 3, "2015" = 4, "2016" = 5)

# length of dataset for model evaluation
eval_day_length = 196


### derivation starts from here
Y = Ys[[ite]]
Yeval = YsEval[[ite]]

Yx = Y[,-dim(Y)[2]]
Yr = Y[,-1]
YxEval = Yeval[,-dim(Yeval)[2]]
YrEval = Yeval[,-1]
N = dim(Y)[1]
alpha = 1/N

########################################################################################
########################################################################################
## Derivation of the lambda sequence 
btest = matrix(0, N, N)
diag(btest) = 0 
test = Yr - btest %*% Yx
test2 = Yx
test3 = test2 %*% t(test)/dim(test2)[2]
test3 = test3 * 100000
lambdatest = function(x) {abs(sum((sign(test3) * max(mean(abs(test3)) - (x*alpha),0))^2) - (x*alpha)^2)}
lambdamax2 = c()
# derive the largest lambda which penalizes each group out
for (lam in 1:N) {
  btest = matrix(0, N, N)
  btest[,lam] = 0
  test = (Yr - btest %*% Yx)
  test2 = Yx[lam,]
  test3 = t(test2 %x% diag(N)) %*% c(test)
  test3 = test3
  lambdatest = function(x) {abs(sum((sign(test3) * pmax(abs(test3) - (x*alpha),0))^2) - (x*alpha)^2)}
  lambdamax2[lam] = optimize(lambdatest, interval = c(1e-10,100000))$minimum
  lambdamax2[lam] = lambdamax2[lam]
}
k=1
lambda2 = max(lambdamax2)
while (lambda2[k] > 0.0000001) {k=k+1; lambda2[k] = lambda2[k-1]*0.5} 
# no regularization on the diagonal
lambda1 = 0
lambdas = expand.grid(lambda1,lambda2)
# END LAMBDA SEQUENCE
#####################################################################################
#####################################################################################
#
# 1. part of code which identifies the active columns. For details, see the 3 step algorithm in the SGNAR paper
############################################
max.iterations = 100
eps = 1e-6

store = list()
iterate = matrix(0, length(lambda1), length(lambda2))
for (j in 1:length(lambda1)) {
  erg = list()
  b = matrix(0, N, N)
  b_old = matrix(0, N, N)
  for (co in 1:length(lambda2)) {
    k=0
    bfull = matrix(1, N, N)
    while (any(abs(b - bfull) > eps) & k<=max.iterations) {
      k = k + 1
      print(paste(co, "|", k))
      bfull = b 
      b = BetaUpdate_SGNAR(0:N, b_old, lambda1[j], lambda2[co], Yx, Yr, alpha) 
      b_old = b
    }
    erg[[co]] = list(b)
    iterate[j,co]=k
  }
  store[[j]] = erg
}

# fix the columns for optimization of parameters 
which_columns = matrix(NA, length(store[[1]]), N)
for (i in 1:length(store[[1]])) {
  curr_val = store[[1]][[i]][[1]]
  diag(curr_val) = NA
  which_columns[i,] = apply(round(curr_val, 4), 2, function(x) !all(x == 0, na.rm = TRUE))
}
#  which_max_col = max(unlist(apply(diff(t(apply(which_columns, 1, as.numeric))),2, function(x) which(x == -1)))) + 1

# 2. part of code which computes the parameter matrix unpenalized. The result serves as starting value for the SGNAR estimation

eps2 = 1e-8
max.iterations = 1000
store_initial = list()
iterate = matrix(0, length(lambda1), length(lambda2))
unique_which_columns = unique(which_columns)
for(j in 1:length(lambda1)) {
  erg = list()
  b_old = matrix(0, N, N)
  for (co in 1:dim(unique_which_columns)[1]) {
    k=0
    b = matrix(0, N, N)
    bfull = matrix(1, N, N)
    iset = c(0, which(unique_which_columns[co,]==TRUE))
    while (any(abs(b - bfull) > eps2) & k<=max.iterations) {
      k = k + 1
      print(paste(co, "|", k))
      bfull = b 
      b = BetaUpdate_SGNAR(iset, b_old, 0, 0, Yx, Yr, alpha) 
      b_old = b
    }
    erg[[co]] = list(b)
    iterate[j,co]=k
  }
  store_initial[[j]] = erg
}

# 3. part of code which derives the parameter matrix with SGNAR, using parts 1. and 2. to improve runtime
which_duplicated = which(!duplicated(which_columns))
eps2 = 1e-8
max.iterations = 1000000
store_new = list()
iterate = matrix(0, length(lambda1), length(lambda2))
select_initial = 0
for(j in 1:length(lambda1)) {
  erg = list()
  b_old = matrix(0, N, N)
  mse3 = c()
  for (co in 1:length(lambda2)) {
    k=0
    b = matrix(0, N, N)
    if (is.element(co, which_duplicated)) {
      select_initial = select_initial + 1
      b_old = store_initial[[1]][[select_initial]][[1]]
    }
    bfull = matrix(1, N, N)
    iset = c(0, which(which_columns[co,]==TRUE))
    while (any(abs(b - bfull) > eps2) & k<=max.iterations) {
      k = k + 1
      print(paste(co, "|", k))
      bfull = b 
      b = BetaUpdate_SGNAR(iset, b_old, lambda1[j], lambda2[co], Yx, Yr, alpha) 
      b_old = b
    }
    erg[[co]] = list(b)
    iterate[j,co]=k
  }
  store_new[[j]] = erg
}

################################
# Model evaluation
################################
cut_day = eval_day_length
mse = c()
counter = 0
b_store = list()
for (i in 1:length(store_new)) {
  for (j in 1:length(store_new[[i]])) {
    counter = counter + 1
    mse[counter] = mean((YrEval[,1:cut_day] - store_new[[i]][[j]][[1]] %*% YxEval[,1:cut_day])^2)
    b_store[[counter]] = store_new[[i]][[j]][[1]]
  }
}

opt_mod = which.min(mse)
lambdas_opt = lambdas[opt_mod,]
b_opt = b_store[[opt_mod]]

All_store = list(store_new, iterate, which_columns, lambdas, Y, Yeval, opt_mod, lambdas_opt, b_opt, mse)

