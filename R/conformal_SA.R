conformal_SA <- function(X, Y, gmm,
                         type,
                         side="two",
                         quantiles=NULL,
                         outfun=outfun, outparams=list(),
                         psfun=Boosting, psparams=list(),
                         trainprop=0.75, conf_wt=NULL,
                         nested=FALSE){

  T <- as.numeric(!is.na(Y))
  inds1 <- which(T == 1)
  inds0 <- which(T == 0)
  n1 <- length(inds1)
  n0 <- length(inds0)
  trainid1 <- sample(n1, floor(n1 * trainprop))
  trainid0 <- sample(n0, floor(n0 * trainprop))
  trainid <- c(inds1[trainid1], inds0[trainid0])

  # convert string to function.
  # (str_outfun is in utils.R, the regressors are in conformal_learners.R)
  if (is.null(outfun)){
    outfun <- switch(type,
                     CQR = quantRF,
                     mean = RF)
  } else if (is.character(outfun)){
    outfun <- str_outfun(outfun[1])
  } else if (is.function(outfun)){
    check_outfun(outfun, type)
  } else {
    stop("outfun must be NULL or a string or a function")
  }
  # print(type)
  # print(outfun)

  if (is.null(psfun)){
    psfun <- Boosting
  } else if (is.character(psfun)){
    psfun <- str_psfun(psfun[1])
  } else if (is.function(psfun)){
    check_psfun(psfun)
  } else {
    stop("psfun must be NULL or a string or a function")
  }

  # Step 1
  # Preliminary set
  Xtrain <- X[trainid, , drop=FALSE]
  Ttrain <- T[trainid]

  psparams0 <- psparams
  psparams <- c(list(Y = Ttrain, X = Xtrain),psparams0)

  # Learn weight function on preliminary set

  wtfun <- function(X,gmm){
    ps <- do.call(psfun, c(list(Xtest = X), psparams))
    if(is.null(conf_wt)){
      return(list(low = 1+(1-ps)/(ps*gmm), high=1+(1-ps)*gmm/ps, pscore=ps))
    } else{
      return(conf_wt(ps,gmm))
    }
  }

  X1 <- X[inds1, ,drop=FALSE]
  Y1 <- Y[inds1]


  # Step 2
  # X1 Y1 are all data in the treated group
  object <- get_score_weight(X1, Y1, gmm,
                             type, side,
                             quantiles,
                             outfun, outparams,
                             wtfun,
                             trainprop, trainid1, nested = nested)
  return(object)
}
