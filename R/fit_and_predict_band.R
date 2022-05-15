fit_and_predict_band <- function(object, X_test, cfun, quantiles=c(0.4,0.6),outparams=list()){
  ########################################################################################
  ##  fit the regression for predicting interval, and predict the interval given data.  ##
  ########################################################################################

  chat <- object$chat
  upper <- chat$upper
  lower <- chat$lower


  #print(mean(upper-lower))


  Xval <- object$X_val
  Yval <- object$Y_val
  Tval <- object$T_val
  ite <- upper-lower

  print(paste0("t = 1 interval length ", mean(ite[which(Tval==1)])))
  print(paste0("t = 0 interval length ", mean(ite[which(Tval==0)])))
  #ite_val <- object$ite_val


  #mz: handled infty in predict.nested()
  # ##################################################################
  # ##          replace all the idinfty with the following          ##
  # ##################################################################
  # idinfty_upper <- which(is.infinite(upper))
  # idinfty_lower <-which(is.infinite(lower))
  # idinfty_index <- union(idinfty_upper, idinfty_lower)
  # print(paste0('percent infty  ',mean(is.infinite(upper))))
  #
  # lower[idinfty_index] <- min(Yval[which(Tval==1)]) - max(Yval[which(Tval==0)])
  # upper[idinfty_index] <- max(Yval[which(Tval==1)]) - min(Yval[which(Tval==0)])
  #
  # if(scaled){
  #   lower <- lower* attr(s.x, 'scaled:scale') + attr(s.x, 'scaled:center')
  # }

  ##################################################################
  ##           train the interval using validation data           ##
  ##################################################################


  low_outparams <- c(list(Y = lower, X = Xval, quantiles = quantiles[1]), outparams)
  up_outparams <- c(list(Y = upper, X = Xval, quantiles = quantiles[2]), outparams)


  # low_outparams <- c(list(Y = lower, X = Xval), outparams)
  # up_outparams <- c(list(Y = upper, X = Xval), outparams)

  low_Cmodel <-  function(X){
    do.call(cfun, c(low_outparams, list(Xtest=X)))
  }
  up_Cmodel <-function(X){
    do.call(cfun, c(up_outparams, list(Xtest=X)))
  }


  #coverage <- mean((ite_val >= low_Cmodel(Xval)) & (ite_val <= up_Cmodel(Xval)),na.rm = TRUE)
  #print(paste0("using the predictor, the coverage is ", mean(coverage)))

  #################################################################
  ##                    predict the interval                     ##
  #################################################################


  interval <- data.frame(lower=low_Cmodel(X_test), upper= up_Cmodel(X_test), y1_mean= mean(ite[which(Tval==1)]), y0_mean = mean(ite[which(Tval==0)]))
  #debug <- data.frame(y1_mean= mean(ite[which(Tval==1)]), y0_mean = mean(ite[which(Tval==0)]))

  #print(paste0("after on the testing site original ", mean(interval[,2] - interval[,1])))

  return(interval)
}
