nested_conformalSA <- function(X, Y1, Y0, T, gmm,
                             type=type, side="two",
                             quantiles=quantiles,
                             outfun=outfun, outparams=list(),
                             psfun=Boosting, psparams=list(),
                             group1_prop=0.5){

  # # setting the unobserved value to na
  # Y1[which(T==0)] <- NA
  # Y0[which(T==1)]<- NA

  #################################################################
  ##             splitting the data into two groups.             ##
  #################################################################
  n = dim(X)[1]
  group1id <- sample(n,floor(group1_prop*n) )
  id <- seq(1, n)
  group2id<- id[!(id %in% group1id)]

  ##################################################################
  ##  USE GROUP 1 to fit the outcome predictors and get Yscores.  ##
  ##################################################################

  obj_treated <- conformal_SA(X[group1id, , drop=FALSE], Y1[group1id], gmm, type = type,  quantiles=quantiles, outfun=outfun, nested=TRUE)
  obj_control <- conformal_SA(X[group1id, , drop=FALSE], Y0[group1id], gmm, type = type, quantiles=quantiles, outfun=outfun, nested=TRUE)

  ########################################################################
  ##  fitting propensity score using all the training data from group1  ##
  ########################################################################

  X_train <- rbind(obj_treated$Xtrain, obj_control$Xtrain)
  T_train_treated <- rep(1, dim(obj_treated$Xtrain)[1])
  T_train_control <- rep(0, dim(obj_control$Xtrain)[1])
  T_train <- append(T_train_treated, T_train_control)

  psparams0 <- psparams
  psparams <- c(list(Y = T_train, X = X_train),psparams0)

  PSmodel <- function(X){
    do.call(psfun, c(psparams, list(Xtest = X)))
  }

  ##################################################################
  ##       calculate the weights on validation of Group 1         ##
  ##################################################################

  wtfun1 <- function(X,gmm){
    ps <- PSmodel(X)
    res <- list(low = (ps)/((1-ps)*gmm), high=ps*gmm/(1-ps))
    return(res)
  }

  wtfun0 <- function(X,gmm){
    ps <- PSmodel(X)
    res<- list(low = (1-ps)/(ps*gmm), high=(1-ps)*gmm/ps)
    return(res)
  }

  obj_treated$wtfun <- wtfun1
  obj_control$wtfun <- wtfun0
  
  obj_treated$wt <- wtfun1(obj_treated$Xval, gmm)
  obj_control$wt <- wtfun0(obj_control$Xval, gmm)

  return(list(
    treated = obj_treated,
    control = obj_control,
    PSmodel = PSmodel,
    group2id = group2id
  ))
}
