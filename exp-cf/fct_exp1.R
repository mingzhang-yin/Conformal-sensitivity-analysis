##################################################################
##                Sensitivity of counterfactuals                ##
##################################################################
library("cfcausal")
library("devtools")
if(exists("cfcausal::summary_CI")){
  rm(list = c("summary_CI"))
}
library("dplyr")
library("ggplot2")
devtools::load_all(".")
rm(list = ls())
#### Get parameters

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("--n", type = "integer", default = 3000, help = "Sample size")
parser$add_argument("--d", type = "integer", default = 20, help = "Dimension")

parser$add_argument("--gmm_star", type = "double", default = 3, help = "SA parameter, >=1")
parser$add_argument("--alpha", type="double", default=0.2, help="mis-coverage")
parser$add_argument("--dtype", type = "character", default = 'hete', help = 'data type, homo or hete')
parser$add_argument("--cftype", type = "integer", default = 2, help = 'confounding type {1,2,3}')

parser$add_argument("--seed", type = "double", default = 1, help = "random seed")
parser$add_argument("--ntrial", type = "integer", default = 5, help = "number of trials")
parser$add_argument("--path", type = "character", default = './exp1/')
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
path <- paste0(path,'loop_fct/', dtype,'/gmm_', gmm_star,'/')

ntrial<- args$ntrial
seed <- args$seed
save <- TRUE
q<- c(alpha/2, 1- (alpha/2))


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

record_ite <- matrix(0,nrow=ntrial,ncol=3)
record_mean <- matrix(0,nrow=ntrial,ncol=3)
record_cqr <- matrix(0,nrow=ntrial,ncol=3)

for (trial in 1:ntrial){
  ##----------------------------------------------------------------
  ##                           Training                           --
  ##----------------------------------------------------------------

  X <- Xfun(n,d)
  Y <- get_Yobs(X)
  ps <- pscorefun(X)
  T <- as.numeric(runif(n)<ps)
  Y[!T] <- NA

  obj_mean <- conformal_SA(X, Y, gmm_star, type = "mean", outfun='RF')
  obj_cqr <- conformal_SA(X, Y, gmm_star, type = "CQR", quantiles=q, outfun='quantRF')
  obj_ite <- conformalCf(X, Y, type = 'mean', outfun ='RF', useCV = FALSE)

  ##---------------------------------------------------------------
  ##                           Testing                           --
  ##---------------------------------------------------------------

  # interval estimate
  ntest <- 10000
  Xtest <- Xfun(ntest,d)
  pstest <- pscorefun(Xtest)
  Ttest <- as.numeric(runif(ntest)<pstest)

  ci_mean <- predict.conformalmsm(obj_mean, Xtest,alpha = alpha)
  ci_cqr <- predict.conformalmsm(obj_cqr, Xtest,alpha = alpha)
  ci_ite <-predict(obj_ite, Xtest, alpha = alpha)



  # test outcome & evaluation
  id1 <- which(Ttest==1)
  id0 <- which(Ttest==0)
  Ytest <- rep(NA,ntest)
  Ytest[id1] <- get_Yobs(Xtest[id1,])
  Ytest_cf <- samplecf(Xtest[id0,],taufun, sdfun, case=cftype, gmm=gmm_star)

  get_fct<- function(set, Ytest, Ytest_cf){
    fcts <- seq(1, 0.5, by=-0.01)
    for(fct in fcts){
      ci <-shrink(set, fc=fct)
      out <-  summary_CI(Ytest,Ytest_cf,ci)
      coverage <- out$cr
        if(coverage < 0.8){
          record <- fct + 0.01
          break
        }}
    return(record)}


  fct_mean <- get_fct(ci_mean, Ytest, Ytest_cf)
  ci_mean <- shrink(ci_mean, fct_mean)
  out_mean <- summary_CI(Ytest,Ytest_cf,ci_mean)
  print("now it's the cqr")

  fct_cqr <- get_fct(ci_cqr, Ytest,Ytest_cf)
  ci_cqr <- shrink(ci_cqr, fct_cqr)
  out_cqr <- summary_CI(Ytest,Ytest_cf,ci_cqr)

  print(paste(min(out_mean$cr),out_mean$len, fct_mean))
  print(paste(min(out_cqr$cr),out_cqr$len,fct_cqr))


  print('##############')

  record_mean[trial,] <- c(out_mean$cr, out_mean$len, fct_mean)
  record_cqr[trial,] <- c(out_cqr$cr, out_cqr$len, fct_cqr)
}

##----------------------------------------------------------------
##                         Save results                         --
##----------------------------------------------------------------

#create a new path for files
folder <- path
dir.create(folder, recursive=TRUE, showWarnings = FALSE)

#fct data
data <- data.frame(fct_min=c(record_mean[,3], record_cqr[,3]),
                   group=rep(c("CSA-M","CSA-Q"),each=ntrial))
write.csv(data, paste0(folder,'fct','.csv'), row.names = FALSE)

