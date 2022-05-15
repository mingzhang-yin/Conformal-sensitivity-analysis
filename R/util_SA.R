samplecf <- function(X,taufun, sdfun, case, gmm, area=NULL, errdist=rnorm){
  Y <- rep(NA, dim(X)[1])
  while(1){
    idnotaccept <- which(is.na(Y))
    if(length(idnotaccept)==0){
      break
    }
    Xsub <- X[idnotaccept,,drop=FALSE]
    n <- dim(Xsub)[1]
    Yproposal <- taufun(Xsub) + sdfun(Xsub) * errdist(n)
    
    if(case==1){
      if(is.null(area)){
        stop("need area for this type of confounding")
      }
      area_sub <- area[idnotaccept,,drop=FALSE]
      m <- taufun(Xsub)
      p_A <- pnorm((area_sub[,2]-m)/sdfun(Xsub)) - pnorm((area_sub[,1]-m)/sdfun(Xsub))
      const <- p_A/gmm + gmm*(1-p_A)
      U <- rep(NA,length(const))
      idx <- (const>=1)
      U[idx] <- (gmm-p_A[idx])/(gmm*(1-p_A[idx]))
      U[!idx] <- p_A[!idx]/(1-gmm*(1-p_A[!idx]))
      tilt <- as.numeric((Yproposal<area_sub[,1]) | (Yproposal>area_sub[,2]))
      tilt[which(tilt==0)] <- 1/(gmm*U[which(tilt==0)])
      probaccept <- tilt
    } 
    
    if(case==2){
      n_sd <- qnorm((2*gmm+1)/(2*gmm+2))
      area <- cbind(taufun(Xsub)-n_sd*sdfun(Xsub), taufun(Xsub)+n_sd*sdfun(Xsub))
      tilt <- as.numeric((Yproposal<area[,1]) | (Yproposal>area[,2]))
      tilt[which(tilt==0)] <- 1/gmm**2
      probaccept <- tilt
    }
    
    if(case==3){ 
      tilt <- (sin(3*(Yproposal-0.2))*(gmm-1)/2)+gmm/2+1/2
      probaccept <- tilt/gmm
    }
    
    if(case==4){ 
      tilt <- (1/(1+exp(-3*(Yproposal-2))))*(gmm-1)+1
      probaccept <- tilt/gmm
    }
    
    idaccept_now <- as.logical(runif(n)<probaccept)
    idaccept <- idnotaccept[idaccept_now]
    Y[idaccept] <- Yproposal[idaccept_now]
    
  }
  return(Y)
}

summary_CI <- function(target,cf_sample,CI){
  
  diff <- CI[, 2] - CI[, 1]
  len <- mean(diff[is.finite(diff)])
  n_inf <- sum(is.infinite(diff))
  
  id_miss <- which(is.na(target))
  id_obs <- which(!is.na(target))  
  cover_obs <- sum((target[id_obs] >= CI[id_obs, 1]) & (target[id_obs] <= CI[id_obs, 2]))
  cover_miss <- sum((cf_sample>=CI[id_miss,1])*(cf_sample<=CI[id_miss,2]))
  cr <- (cover_obs+cover_miss)/length(target)

  return(list(cr = cr, len = len, n_inf=n_inf))
}

# samplecf <- function(X,taufun, sdfun, errdist, tilting, prob_down=NULL){
#   Y <- rep(NA, dim(X)[1])
#   while(1){
#     idnotaccept <- which(is.na(Y))
#     if(length(idnotaccept)==0){
#       break
#     }
#     Xsub <- X[idnotaccept,,drop=FALSE]
#     n <- dim(Xsub)[1]
#     Yproposal <- taufun(Xsub) + sdfun(Xsub) * errdist(n)
#     if(!is.null(prob_down)) {
#       n_sd <- qnorm((1+prob_down)/2)
#       area <- cbind(taufun(Xsub)-n_sd*sdfun(Xsub), taufun(Xsub)+n_sd*sdfun(Xsub))
#       probaccept <- tilting(Yproposal,rg=area)$prob
#     }else{
#       probaccept <- tilting(Yproposal)$prob
#     }
#     idaccept_now <- as.logical(runif(n)<probaccept)
#     idaccept <- idnotaccept[idaccept_now]
#     Y[idaccept] <- Yproposal[idaccept_now]
#   }
#   return(Y)
# }

# samplecf <- function(X,taufun, sdfun, errdist, tilting, gmm=NULL){
#   Y <- rep(NA, dim(X)[1])
#   while(1){
#     idnotaccept <- which(is.na(Y))
#     if(length(idnotaccept)==0){
#       break
#     }
#     Xsub <- X[idnotaccept,,drop=FALSE]
#     n <- dim(Xsub)[1]
#     Yproposal <- taufun(Xsub) + sdfun(Xsub) * errdist(n)
#     if(!is.null(gmm)) {
#       n_sd <- qnorm((2*gmm+1)/(2*gmm+2))
#       area <- cbind(taufun(Xsub)-n_sd*sdfun(Xsub), taufun(Xsub)+n_sd*sdfun(Xsub))
#       probaccept <- tilting(Yproposal,rg=area)$prob
#     }else{
#       probaccept <- tilting(Yproposal)$prob
#     }
#     idaccept_now <- as.logical(runif(n)<probaccept)
#     idaccept <- idnotaccept[idaccept_now]
#     Y[idaccept] <- Yproposal[idaccept_now]
#   }
#   return(Y)
# }


# tilting <- function(y,rg,bound=gmm_star){
#   tilt <- as.numeric((y>rg[,1]) & (y<rg[,2]))
#   tilt[which(tilt==0)] <- bound**2
#   prob <- tilt/(bound**2)
#   return(list(tilt=tilt,prob=prob))
# }


# summary_CI <- function(target,cf_sample,CI){
#   id_miss <- which(is.na(target))
#   id_obs <- which(!is.na(target))
#   diff <- CI[, 2] - CI[, 1]
#   len <- mean(diff[is.finite(diff)])
#   n_inf <- sum(is.infinite(diff))
#
#   cover_obs <- sum((target[id_obs] >= CI[id_obs, 1]) & (target[id_obs] <= CI[id_obs, 2]))
#   cover_miss <- apply((cf_sample>=CI[id_miss,1])*(cf_sample<=CI[id_miss,2]),2,sum)
#   cr <- (cover_obs+cover_miss)/length(target)
#   return(list(cr = cr, len = len, n_inf=n_inf))
# }

# MIS <- function(x,CI,alpha){
#   mis <- CI[, 2]-CI[, 1] + 2*(x-CI[, 2])*as.numeric(x >= CI[, 2])/alpha +
#                              2*(CI[, 1]-x)*as.numeric(x <= CI[, 1])/alpha
#   return(mis)
# }


# summary_CI <- function(target,cf_sample,CI,alpha){
#   id_miss <- which(is.na(target))
#   id_obs <- which(!is.na(target))
# 
#   diff <- CI[, 2] - CI[, 1]
#   len <- mean(diff[is.finite(diff)])
#   n_inf <- sum(is.infinite(diff))
# 
#   cover_obs <- sum((target[id_obs] >= CI[id_obs, 1]) & (target[id_obs] <= CI[id_obs, 2]))
#   cover_miss <- apply((cf_sample>=CI[id_miss,1])*(cf_sample<=CI[id_miss,2]),2,sum)
#   cr <- (cover_obs+cover_miss)/length(target)
# 
#   #cr <- (cover_obs)/length(id_obs)
#   #cr <- (cover_miss)/length(id_miss)
# 
#   MIS_obs <- sum(MIS(target[id_obs],CI[id_obs, ],alpha),na.rm=TRUE)
#   MIS_MISS <- apply(MIS(cf_sample,CI[id_miss,],alpha),2,sum,na.rm=TRUE)
#   mis <- (MIS_obs+MIS_MISS)/length(target)
#   return(list(cr = cr, len = len, n_inf=n_inf, MIS=mis))
# }



