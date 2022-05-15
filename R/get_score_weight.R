get_score_weight <- function(X1, Y1, gmm,
                             type, side,
                             quantiles,
                             outfun, outparams,
                             wtfun,
                             trainprop, trainid1, nested=FALSE){
  #############
  # type~{CQR, mean} should match outfun~{quantRF,RF}
  # X1, Y1 are observed data
  # trainid1 are indices of observed data in the preliminary set
  #############

  Xtrain <- X1[trainid1, ,drop=FALSE]
  Ytrain <- Y1[trainid1]

  Xval <- X1[-trainid1, ,drop=FALSE]
  Yval <- Y1[-trainid1]
  nval <- dim(Xval)[1]


  # Learn output function on preliminary set with observed Y
  outparams <- c(list(Y = Ytrain, X = Xtrain, quantiles = quantiles), outparams)
  Ymodel <- function(X){
    do.call(outfun, c(outparams, list(Xtest = X)))
  }

  # score and weight of the validation set
  Yhat <- Ymodel(Xval)

  Yscore <- conformalScore(Yval, Yhat, type, side)
  if(nested){
    wt <- NA
    wtfun <- NA
  }
  else{
    wt <- wtfun(Xval,gmm)
  }

  object <- list(Yscore = Yscore, wt = wt,
                 Ymodel = Ymodel, wtfun = wtfun,
                 gmm = gmm,
                 type = type,
                 side = side,
                 quantiles = quantiles,
                 trainprop = trainprop,
                 Xtrain=Xtrain,
                 Xval = Xval,
                 trainid = trainid1)
  class(object) <- "conformalSensitivity"
  return(object)
}
