# This function updates the parameters in the specified columns of the parameter matrix. 
#
# Note: Only run the file 'main_ParaMatEst.R'. This code is an internal file needed for the algorithm. 
#
# Parameters: 
# 'iset' : vector of numbers which indicate the columns to iterate over. '0' refers to the parameters on the diagonal which are a separate group
# 'b' : the current parameter matrix which will be updated 
# 'lambda11' : vector of regularization parameters for the parameters on the diagonal of the parameter matrix
# 'lambda22' : vector of regularization parameters for the parameters off-diagonal of the parameter matrix
# 'Yx' : matrix of past observations to estimate parameters (regressors)
# 'Yr' : matrix of current observations which are the target (dependent) variables
# 'alpha' : mixing parameter for the impact of the regularization parameter on the individual parameters and each entire column of parameters

BetaUpdate_SGNAR = function(iset, b, lambda11, lambda22, Yx, Yr, alpha) {
  # reordering the iset so that the diagonal is estimated first and then in decreasing order
  # the columns which have likely the most explanatory effect in the model. The ordering 
  # improvest the runtime of the algoritm 
  if (is.element(0, iset)) {
    new_iset = 0
    iset = iset[-which(0 == iset)]
  } 
  store_testval = c()
  if (length(iset) == 0) {
    iset = new_iset
  } else {
    for (i in iset) {
      a = b
      a[-i,i] = 0
      YR = (a %*% Yx - Yr)[-i,] 
      M = Yx[i,]
      testval = M %*% t(YR)/length(M)
      store_testval = c(store_testval, sqrt(sum((sign(testval) * pmax(abs(testval) - lambda22*alpha,0))^2)))
    }
    iset = c(new_iset, iset[order(store_testval, decreasing = TRUE)])
  }
  
  for (i in iset) {
    # active set estimation over iset
    if (i == 0) {
      a = b
      diag(a) = 0
      # unexplained information (error) by the current parameters
      YR = a %*% Yx - Yr 
      M = matrix(0, ncol = dim(Yx)[1], nrow = dim(Yx)[1] * dim(Yx)[2])
      for ( k in 1:dim(Yx)[2]) {M[(1 + dim(Yr)[1]*(k-1)):(dim(Yr)[1]*k),]=diag(Yx[,k])} 
      testval = t(M) %*% matrix(YR,ncol=1)/length(YR)
    } else {
      a = b
      a[-i,i] = 0
      YR = (a %*% Yx - Yr)[-i,] 
      M = Yx[i,] 
      testval = M %*% t(YR)/length(M)
    }
    
    # set all parameters on the diagonal or in a column to 0 if SGNAR penalty dictates so
    if (i == 0 & sqrt(sum((sign(testval) * pmax(abs(testval) - lambda11*alpha,0))^2)) <= lambda11*(1-alpha)) {
      diag(b) = 0
    } else if (i > 0 & sqrt(sum((sign(testval) * pmax(abs(testval) - lambda22*alpha,0))^2)) <= lambda22*(1-alpha)) {
      b[-i,i] = 0
    } else {
      l = 1
      if (i == 0) {
        beta = diag(b)
        lambda = lambda11
      } else {
        beta = b[-i,i]
        lambda = lambda22
      }
      
      # parameter update for the diagonal
      if (i == 0) {
        const2 = crossprod(M)
        const3 = crossprod(M,c(YR))
        bupdate = beta + 1
        while (any(abs(beta - bupdate) > 1e-6)) { 
          # for details on the parameter update, see appendix in SGNAR paper
          beta = bupdate
          tval = 1
          grad = const2 %*% beta + const3  
          const1 = 0.5 * sum((M %*% beta + c(YR))^2)
          tgrad = t(grad)
          optStepSizeLeft = 2
          optStepSizeRight = 1
          while (optStepSizeLeft > optStepSizeRight) { 
            z = beta - tval * grad
            beta1 = beta[1]
            Sfunc = sign(mean(z))*max(mean(abs(z)) - lambda*alpha,0)
            Ufunc = max(1 - ((tval*lambda*(1-alpha)) / sqrt(sum(Sfunc^2))), 0) * Sfunc
            optStepSizeRight = const1 + sum(tgrad * (Ufunc - beta1)) + 1/(2*tval) * sum((Ufunc - beta)^2)
            optStepSizeLeft = 0.5 * sum((c(Yx) * Ufunc + c(YR))^2)
            tval = 0.8 * tval
          }
          bupdate = beta + l/(l + 3) * (Ufunc - beta)
          l = l + 1
        }
        if (i == 0) {
          diag(b) = bupdate
        } else {
          b[-i,i] = bupdate
        }
        # parameter update for the off-diagonal parameters (columns)
      } else if (i > 0) {
        const2 = sum(M^2) 
        const3 = YR %*% M 
        bupdate = beta + 1
        while (any(abs(beta - bupdate) > 1e-06)) { 
          tval = 1
          beta = bupdate
          grad = const2 * beta + const3 
          const1 = 0.5 * sum((M %x% beta + c(YR))^2) 
          tgrad = t(grad)
          optStepSizeLeft = 2
          optStepSizeRight = 1
          while (optStepSizeLeft > optStepSizeRight) { 
            z = beta - tval * grad
            Sfunc = sign(z)*pmax(abs(z) - lambda*alpha,0)
            Ufunc = max(1 - ((tval*lambda*(1-alpha)) / sqrt(sum(Sfunc^2))), 0) * Sfunc
            optStepSizeRight = const1 + tgrad %*% (Ufunc - beta) + 1/(2*tval) * sum((Ufunc - beta)^2)
            optStepSizeLeft = 0.5 * sum((M %x% Ufunc + c(YR))^2)
            tval = 0.8 * tval
          }
          bupdate = beta + l/(l + 3) * (Ufunc - beta)
          l = l + 1
        }
        if (i == 0) {
          diag(b) = bupdate
        } else {
          b[-i,i] = bupdate
        }
      }
    }
  }
  return(b)
}
