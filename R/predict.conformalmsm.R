predict.conformalmsm <- function(object, Xtest, alpha = 0.2, wthigh = 20, wtlow = 0.05){
  type <- object$type
  Yhat_test <- object$Ymodel(Xtest)
  wt_test <- object$wtfun(Xtest,object$gmm)

  if(type == 'mean'){
    # mean prediction
    Yslack <- cutoff_SA(object$Yscore,object$wt,wt_test,alpha)
    Ylo <- Yhat_test - Yslack
    Yup <- Yhat_test + Yslack
    interval <- data.frame(lower = Ylo, upper = Yup)}
  else if(type == 'CQR'){
    # quantile regression
    Yslack <- cutoff_SA(object$Yscore,object$wt,wt_test,alpha)
    Ylo<- Yhat_test[, 1] - Yslack
    Yup <- Yhat_test[, 2] + Yslack
    interval <- data.frame(lower = Ylo, upper = Yup)
  }
  return(interval)
}
