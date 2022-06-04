##################################################################
##                Sensitivity of counterfactuals                ##
##################################################################
# setwd("~/Conformal-sensitivity-analysis")
library("cfcausal")
library("devtools")
library("dplyr")
library("ggplot2")
library("nloptr")
library("parallel")
devtools::load_all(".")
if(exists("cfcausal::summary_CI")){
  rm(list = c("summary_CI"))
}
rm(list = ls())
#### Get parameters
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("--n", type = "integer", default = 3000, help = "Sample size")
parser$add_argument("--d", type = "integer", default = 20, help = "Dimension")

parser$add_argument("--gmm_star", type = "double", default = 1.5, help = "SA parameter, >=1")
parser$add_argument("--alpha", type="double", default=0.2, help="mis-coverage")
parser$add_argument("--dtype", type = "character", default = 'hete', help = 'data type')
parser$add_argument("--cftype", type = "integer", default = 2, help = 'confounding type {1,2,3}')

parser$add_argument("--seed", type = "double", default = 1, help = "random seed")
parser$add_argument("--trial", type = "integer", default = 1, help = "id of trial")
# parser$add_argument("--ntrial", type = "integer", default = 10, help = "number of trials")
parser$add_argument("--path", type = "character", default = '/proj/sml_netapp/projects/conformal/exp1')
parser$add_argument("--fct", type = "double", default = 1, help = 'shrink factor of band')

args <- parser$parse_args()
n <- args$n
d <- args$d
alpha <- args$alpha
gmm_star <- args$gmm_star

dtype <- args$dtype
cftype <- args$cftype
path <- args$path
fct <- args$fct
path <- paste0(path,'/loop_gmm/', dtype, '/')
# ntrial<- args$ntrial
seed <- args$seed
q<- c(alpha/2, 1- (alpha/2))
ntest <- 500

#create a new path for files
folder<- paste0(path, gmm_star, "/")
print(folder)
dir.create(folder, recursive=TRUE, showWarnings = FALSE)

##---------------------------------------------------------------
##                       Data generation                       --
##---------------------------------------------------------------

if(dtype=='homo'){
    Xfun <- function(n, d){
      matrix(runif(n * d), nrow = n, ncol = d)
    }
    sdfun <- function(X){
      rep(1, nrow(X))
    }
} else{
    rho <- 0.9
    Xfun <- function(n, d){
      X <- matrix(rnorm(n * d), nrow = n, ncol = d)
      fac <- rnorm(n)
      X <- X * sqrt(1 - rho) + fac * sqrt(rho)
      pnorm(X)
    }
    sdfun <- function(X){
      runif(nrow(X),0.5,1.5)
    }
}

taufun <- function(X){
  # browser()
  2 / (1 + exp(-5 * (X[, 1] - 0.5))) * 2 / (1 + exp(-5 * (X[, 2] - 0.5)))
}
pscorefun <- function(X){
  (1 + pbeta(1-X[, 1], 2, 4)) / 4
}
errdist <- rnorm

get_Yobs <- function(X){
  return(taufun(X) + sdfun(X) * errdist(dim(X)[1]))
}

shrink <- function(set,fc){
  newset <- set
  idx <- is.finite(set[,1])
  center <- (set[idx,2] + set[idx,1])/2
  halflen <- (set[idx,2] - set[idx,1])/2
  newset[idx,] <- cbind(center-halflen*fc, center+halflen*fc)
  return(newset)
}

##----------------------------------------------------------------
##                          Estimation                          --
##----------------------------------------------------------------


# record_mean <- matrix(0,nrow=ntrial,ncol=3)
# record_cqr <- matrix(0,nrow=ntrial,ncol=3)

# for(sss in 1:5){
  ##----------------------------------------------------------------
  ##                           Training                           --
  ##----------------------------------------------------------------

  X <- Xfun(n,d)
  Y <- get_Yobs(X)
  ps_true <- pscorefun(X)
  T <- as.numeric(runif(n)<ps_true)
  Y[!T] <- NA
  
  obj_mean <- conformal_SA(X, Y, gmm_star, type = "mean", outfun='RF')
  # obj_cqr <- conformal_SA(X, Y, gmm_star, type = "CQR", quantiles=q, outfun='quantRF')
  
  ##---------------------------------------------------------------
  ##                           Testing                           --
  ##---------------------------------------------------------------
  
  # interval estimate
  
  Xtest <- Xfun(ntest,d)
  pstest_true <- pscorefun(Xtest)
  Ttest <- as.numeric(runif(ntest)<pstest_true)
  
  object <- obj_mean
  
  type <- object$type
  side <- object$side
  Yhat_test <- object$Ymodel(Xtest)
  wt_test <- object$wtfun(Xtest,object$gmm) #test weights
  
  wt <- object$wt  #calibration weights
  wtlow <- wt$low
  wthigh <- wt$high
  ps_val <- wt$pscore
  
  score <- object$Yscore
  ord <- order(score)
  score <- score[ord]
  wtlow <- wtlow[ord]
  wthigh <- wthigh[ord]
  ps_val <- ps_val[ord]
  
  
  objective <- function(w, K0) {
      return(-1 * sum(tail(w, K0)) / sum(w))
  }
  
  # cutoff <- rep(1, ntest)
  # start_time <- Sys.time()
  # for(j in 1:1){
  #
  #   wt_combine_low <- c(wtlow, wt_test$low[j])
  #   wt_combine_high <- c(wthigh, wt_test$high[j])
  #   ps_combine <- c(ps_val, wt_test$pscore[j])
  #
  #   n_combine <- length(wt_combine_low)
  #   con <- function(w, coef0=ps_combine) {
  #       return(mean(w * coef0) - 1)
  #   }
  #
  #   for (k in n_combine:2) {
  #     w0 <- c(wt_combine_low[1:k-1], wt_combine_high[k:n_combine])
  #     alpha_hat <- sum(wt_combine_high[k:n_combine])/sum(w0)
  #     if(alpha_hat>=(alpha+1e-12)){
  #       break
  #     }
  #   }
  #
  #   if(abs(con(w0))<1e-6){
  #     cut_position <- k
  #   }else{
  #     # binary search
  #     L = 1; R = k;
  #     print(c(k, n_combine))
  #     m <- ifelse(R>50, R-30, ceiling(0.5*L + 0.5*R))
  #     while(R-L>1){
  #       res <- slsqp(w0, objective, heq = con, lower = wt_combine_low, upper = wt_combine_high, K0 = (n_combine-m+1),
  #           control = list(maxeval = 1000, xtol_rel = 1e-5, ftol_abs = 1e-4))
  #       alpha_hat <- (-res$value)
  #       if(alpha_hat>=(alpha+1e-12)){
  #         L <- m
  #       }else{
  #         R <- m
  #       }
  #       w0 <- res$par
  #       print(abs(con(w0)))
  #       m <- ceiling(0.5*L + 0.5*R)
  #     }
  #
  #     cut_position <- R
  #     print(cut_position)
  #   }
  #
  #   cutoff[j] <- ifelse(cut_position==n_combine, Inf, score[cut_position])
  #
  # }
  # end_time <- Sys.time()
  # print(paste0("Time=", end_time-start_time))



  sim <- function(j, wtlow, wthigh, wt_test, ps_val, alpha){
    wt_combine_low <- c(wtlow, wt_test$low[j])
    wt_combine_high <- c(wthigh, wt_test$high[j])
    ps_combine <- c(ps_val, wt_test$pscore[j])
  
    n_combine <- length(wt_combine_low)
    con <- function(w, coef0=ps_combine) {
        return(mean(w[1:length(w)-1] * coef0[1:length(coef0)-1]) - 1)
    }
  
    for (k in n_combine:2) {
      w0 <- c(wt_combine_low[1:k-1], wt_combine_high[k:n_combine])
      alpha_hat <- sum(wt_combine_high[k:n_combine])/sum(w0)
      if(alpha_hat>=(alpha+1e-12)){
        break
      }
    }
    cut_position <- k
  
    if((abs(con(w0))>1e-6)){
      # binary search
      L = 1; R = k;
      m <- ifelse(R>50, R-10*gmm_star, ceiling(0.5*L + 0.5*R))
      while(R-L>1){
        res <- slsqp(w0, objective, heq = con, lower = wt_combine_low, upper = wt_combine_high, K0 = (n_combine-m+1),
            control = list(maxeval = 1000, xtol_rel = 1e-5, ftol_abs = 1e-4))
        alpha_hat <- (-res$value)
        if(alpha_hat>=(alpha+1e-12)){
          L <- m
        }else{
          R <- m
        }
        w0 <- res$par
        m <- ceiling(0.5*L + 0.5*R)
      }
  
      cut_position <- R
      print(c(cut_position, k, n_combine))
    }
  
    cutoff_j <- ifelse(cut_position==n_combine, Inf, score[cut_position])
    return(c(j, cutoff_j))
  }
  
  start_time <- Sys.time()
  numWorkers <- detectCores()-1
  cl <- makeCluster(numWorkers,type="FORK")
  res <- parLapply(cl, 1:ntest, sim, wtlow=wtlow, wthigh=wthigh, wt_test=wt_test, ps_val=ps_val, alpha=alpha)
  stopCluster(cl)
  end_time <- Sys.time()
  
  
  results = do.call(rbind, res)
  results = results[order(results[,1],decreasing=FALSE),]
  
  cutoff <- results[,2]
  Ylo <- Yhat_test - cutoff
  Yup <- Yhat_test + cutoff
  ci_mean <- data.frame(lower = Ylo, upper = Yup)
  
  
  # Ylo<- Yhat_test[, 1] - Yslack
  # Yup <- Yhat_test[, 2] + Yslack
  # interval <- data.frame(lower = Ylo, upper = Yup)
  
  
  # test outcome & evaluation
  id1 <- which(Ttest==1)
  id0 <- which(Ttest==0)
  Ytest <- rep(NA,ntest)
  #print("before get yobs")
  Ytest[id1] <- get_Yobs(Xtest[id1,, drop=F])
  #print("after get y obs")
  Ytest_cf <- samplecf(Xtest[id0,, drop=F],taufun, sdfun, case=cftype, gmm=gmm_star)
  
  out_mean <- cfcausal::summary_CI(Ytest,Ytest_cf,ci_mean)
  
  # record_mean[trial,] <- c(out_mean$cr, out_mean$len, out_mean$n_inf)
  
  ##----------------------------------------------------------------
  ##                         Save results                         --
  ##----------------------------------------------------------------
  
  
  
  # #coverage data
  # data_cov <- data.frame(Coverage=c(as.vector(record_mean[,1]),as.vector(record_cqr[,1])),
  #             group=rep(c("CSA-M","CSA-Q"),each=ntrial))
  # data_len <- data.frame(Interval_length=c(record_mean[,2], record_cqr[,2]),
  #                    group=rep(c("CSA-M","CSA-Q"),each=ntrial))
  
  
  # write.csv(data_cov, paste0(folder,'coverage', '.csv'), row.names = FALSE)
  # write.csv(data_len, paste0(folder,'len','.csv'), row.names = FALSE)
  
  
  #coverage data
  data_cov <- data.frame(Coverage=min(out_mean$cr), group="CSSA-M")
  data_len <- data.frame(Interval_length=out_mean$len, group="CSSA-M")
  
  write.csv(data_cov, paste0(folder,'coverage_', args$trial, '.csv'), row.names = FALSE)
  write.csv(data_len, paste0(folder,'len_', args$trial, '.csv'), row.names = FALSE)
  # print("hiiiiiiiiiiiiiiii")
  print(paste(out_mean$cr,out_mean$len, out_mean$n_inf))
  print(end_time-start_time)

# }