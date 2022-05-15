predict.nested <- function(object, X, Y, T, alpha = 0.2, wthigh = 20, wtlow = 0.05){
  ##################################################################
  ##    use group 2 data to create a bunch of prediction bands    ##
  ##################################################################
  PSmodel <- object$PSmodel
  wtfun <- object$wtfun
  group2id <- object$group2id
  ite <- object$ite
  Xtest <- X[group2id, , drop=FALSE]
  Ytest <- Y[group2id]
  Ttest <- T[group2id]
  #ite_test <- ite[group2id]

  inds1 <- which(Ttest == 1)
  inds0 <- which(Ttest == 0)

  Xtest_treated <- Xtest[inds1, , drop=FALSE]
  Ytest_treated <- Ytest[inds1]
  Xtest_control <- Xtest[inds0, , drop=FALSE]
  Ytest_control <- Ytest[inds0]

  #################################################################
  ##               create updated prediction bands               ##
  #################################################################
  c1 <- predict.conformalmsm(object$treated, Xtest_control,alpha = alpha)
  c0 <- predict.conformalmsm(object$control, Xtest_treated,alpha = alpha)

  copy <- c0
  c0$lower<- Ytest_treated- copy$upper
  c0$upper <-Ytest_treated- copy$lower

  c1$lower <- c1$lower - Ytest_control
  c1$upper <- c1$upper - Ytest_control

  chat <- rbind(c1,c0)


  ##################################################################
  ##            checking the coverage for group 2 data            ##
  ##################################################################


  #handle infty
  inf_idx <- union(which(is.infinite(chat[,1])), which(is.infinite(chat[,2])))
  chat[inf_idx,1] <- min(Y[which(T==1)]) - max(Y[which(T==0)])
  chat[inf_idx,2] <- max(Y[which(T==1)]) - min(Y[which(T==0)])
  
  #print(paste0("how many infinity ", length(inf_idx), inf_idx))
  #print("i am looking at the prediction rate for the second part")
  # ite_true <- append(ite_test[inds0], ite_test[inds1])
  # coverage <- mean((ite_true >= chat[, 1]) & (ite_true <= chat[, 2]),na.rm = TRUE)
  # 
  # print(paste0("at part 2, the average coverage is ", coverage))




  ##################################################################
  ##          all the data used for the prediction bands          ##
  ##################################################################
  X_val = rbind(Xtest_control,Xtest_treated )
  Y_val = append(Ytest_control,Ytest_treated)
  T_val = append(rep(0, length(Ytest_control)), rep(1, length(Ytest_treated)))
  #ite_val = append(ite_test[inds0], ite_test[inds1])

  
  return(list(chat = chat,
              X_val = X_val,
              Y_val = Y_val,
              T_val = T_val
              #ite_val = ite_val
              ))
}
