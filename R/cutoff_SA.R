cutoff_SA <- function(score, wt, wt_test, alpha){
  #alpha is miscoverage
  wtlow <- wt$low
  wthigh <- wt$high
  ord <- order(score)
  score <- score[ord]
  wtlow <- wtlow[ord]
  wthigh <- wthigh[ord]
  ntest <- length(wt_test$low)
  wt_combine <- t(replicate(ntest, wtlow))
  wt_combine <- cbind(wt_combine,wt_test$high)
  c(ntest,n_combine) := dim(wt_combine)

  findn <- rep(1, ntest) #index of cutoff
  for (i in n_combine:2) {
    probmass <- rowSums(wt_combine[,i:n_combine,drop=FALSE])/rowSums(wt_combine)
    findbool <- (probmass>=(alpha+1e-12)) #True if the cutoff is found
    if (sum(findbool)>0){
      findn[findbool] <- pmax(findn[findbool],i)
    }
    if (sum(findbool) == ntest){
      break  #all cutoffs are found
    }
    wt_combine[,i-1] <- wthigh[i-1]
  }
  
  id_unbounded <- (findn==n_combine) 
  cutoff <- rep(0,ntest)

  if (sum(id_unbounded)>0){
    cutoff[id_unbounded] <- Inf
    cutoff[!id_unbounded] <- score[findn[!id_unbounded]]
  }else{
    cutoff <- score[findn]
  }
  return(cutoff)
}
