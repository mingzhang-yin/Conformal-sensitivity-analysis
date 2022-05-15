######## Sensitivity of ITE (synthetic)########
library("devtools")
if(exists("cfcausal::summary_CI")){
  rm(list = c("summary_CI"))
}
devtools::load_all(".")
library("cfcausal")
library("dplyr")
library("ggplot2")
library("bannerCommenter")
options(scipen=999)


#### Get parameters
if (!interactive()){
  print("i am here")
  suppressPackageStartupMessages(library("argparse"))
  parser <- ArgumentParser()
  parser$add_argument("--n", type = "integer", default = 3000, help = "Sample size")
  parser$add_argument("--d", type = "integer", default = 20, help = "Dimension")
  parser$add_argument("--gmm_star", type = "double", default = 3, help = "SA parameter, >=1")
  parser$add_argument("--alpha", type="double", default=0.2, help="miscoverage")
  parser$add_argument("--cftype", type="integer", default=2, help="confounding type")
  parser$add_argument("--dtype", type="character", default='homo', help="data type, homo or het")
  parser$add_argument("--fct", type="double", default=1, help="shrinkage, <=1")
  parser$add_argument("--save", type="logical", default=TRUE, help="save")
  parser$add_argument("--seed", type = "double", default = 1, help = "random seed")
  parser$add_argument("--ntrial", type = "integer", default =50, help = "number of trials")
  parser$add_argument("--path", type = "character", default = './loop_fct/', help = "save location")
  args <- parser$parse_args()
  n <- args$n
  d <- args$d
  alpha <- args$alpha
  ntrial<- args$ntrial
  dtype <- args$dtype
  cftype<- args$cftype
  fct <- args$fct
  seed <- args$seed
  gmm_star <- args$gmm_star
  save <- args$save
  path = args$path
} else {
  n <- 1000
  d <- 10
  ntrial <- 2
  gmm_star <- 2#SA parameter
  alpha <- 0.2   #mis-coverage
  dtype <- 'hete' #data type
  cftype <- 2  #confounding type
  fct <- 1 #shrink length
  save<- TRUE
  path = '../../../res/conformal/dep/'
}

q <- c(alpha/2, 1-alpha/2)
#########

######## experiment ##########
##homoscedastic errors + independent covariates

#Data generation
if(fct <=0){
  fct <- exp(fct)
}

Xfun <- function(n, d){
  matrix(runif(n * d), nrow = n, ncol = d)
}

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


get_Y1obs <- function(X){
  return(taufun(X) + sdfun(X) * errdist(dim(X)[1]))
}

taufun0 <- function(X){
  taufun(X) + 10*sin(X[, 3])*(1/(1+exp(-5*X[, 3])))
}

get_Y0obs <- function(X){
  return(taufun0(X) + sdfun(X) * errdist(dim(X)[1]))
}

shrink <- function(set,fc){
  newset <- set
  idx <- is.finite(set[,1])
  center <- (set[idx,2] + set[idx,1])/2
  halflen <- (set[idx,2] - set[idx,1])/2
  newset[idx,] <- cbind(center-halflen*fc, center+halflen*fc)
  return(newset)
}


print_list <- list("sa_mean", "sa_cqr", "ite_nuc", "sa_naive")
record <- replicate(length(print_list),matrix(0,nrow=ntrial,ncol=3), simplify=FALSE)

for (trial in 1:ntrial){
  ##---------------------------------------------------------------
  ##                        Generate observed data                         -
  ##---------------------------------------------------------------
  X <- Xfun(n,d)
  ps <- pscorefun(X)
  T <- as.numeric(runif(n)<ps)

  Y0 <- get_Y0obs(X)
  Y1 <- get_Y1obs(X)

  Y_obs <- Y1*T + Y0*(1-T)

  Y1[which(T==0)] <- NA
  Y0[which(T==1)] <- NA

  ##----------------------------------------------------------------
  ##                          bonferroni                           -
  ##----------------------------------------------------------------

  obj1_ite <- conformal_SA(X, Y1, gmm_star, type = "mean", outfun='RF')
  obj0_ite <- conformal_SA(X, Y0, gmm_star, type = "mean", outfun='RF')

  ##----------------------------------------------------------------
  ##      inexact ite method assuming no unobserved confounder     -
  ##----------------------------------------------------------------


  CIfun_inexact <- conformalIte(X, Y_obs, T, alpha = alpha,
                                algo = "nest", exact=FALSE, type = "CQR",
                                #lofun = 'RF', upfun = 'RF', citype = "mean",
                                quantiles = c(alpha/2, 1- (alpha/2)), outfun = "quantRF",  useCV = FALSE)


  ##----------------------------------------------------------------
  ##                            Train on Group1
  ##----------------------------------------------------------------
  obj_mean <- nested_conformalSA(X, Y1, Y0, T, gmm_star, type = "mean",quantiles=list(), outfun='RF')
  obj_cqr <- nested_conformalSA(X, Y1, Y0, T, gmm_star, type = "CQR",quantiles=q, outfun='quantRF')

  ##----------------------------------------------------------------
  ##                  getting prediction bands on Group2
  ##----------------------------------------------------------------

  obj_bands_mean <- predict.nested(obj_mean, X, Y_obs, T, alpha = alpha)
  obj_bands_cqr <- predict.nested(obj_cqr, X, Y_obs, T, alpha = alpha)


  ##----------------------------------------------------------------
  ##                        generate testing data                    -
  ##----------------------------------------------------------------

  ##Testing
  ntest <- 10000
  Xtest <- Xfun(ntest,d)
  pstest <- pscorefun(Xtest)
  Ttest <- as.numeric(runif(ntest)<pstest)
  id1 <- which(Ttest==1)
  id0 <- which(Ttest==0)

  Y0test <- rep(NA,ntest)
  Y0test[id0] <- get_Y0obs(Xtest[id0,])
  Y0test[id1] <- samplecf(Xtest[id1,], taufun0, sdfun, case=cftype, gmm=gmm_star)

  Y1test <- rep(NA,ntest)
  Y1test[id1] <- get_Y1obs(Xtest[id1,])
  Y1test[id0] <- samplecf(Xtest[id0,],taufun, sdfun, case=cftype, gmm=gmm_star)

  ##----------------------------------------------------------------
  ##                       ITE & evaluation                    -
  ##----------------------------------------------------------------

  ite <- Y1test - Y0test
  ci_mean_copy <-fit_and_predict_band(obj_bands_mean,Xtest, 'quantRF')
  ci_cqr_copy <- fit_and_predict_band(obj_bands_cqr,Xtest, 'quantRF')

  ci_mean <-shrink(ci_mean_copy, fc=fct)
  ci_cqr <- shrink(ci_cqr_copy, fc=fct)
  ci_mean[, 3] <- ci_mean_copy[,3]
  ci_mean[, 4] <- ci_mean_copy[, 4]

  ci_cqr[, 3] <- ci_cqr_copy[,3]
  ci_cqr[, 4] <- ci_cqr_copy[, 4]

  #bonferroni

  ci0_ite <- predict.conformalmsm(obj0_ite, Xtest,alpha = alpha/2)
  ci1_ite <- predict.conformalmsm(obj1_ite, Xtest,alpha = alpha/2)
  ci_ite <- cbind(ci1_ite[,1] - ci0_ite[,2], ci1_ite[,2] - ci0_ite[,1])

  #ite-nuc
  ci_inexact <- CIfun_inexact(Xtest)

  ci_list <- list(ci_mean, ci_cqr, ci_inexact, ci_ite)



  for(i in 1:length(ci_list)){
    ci <- ci_list[[i]]
    coverage <- mean((ite >= ci[, 1]) & (ite <= ci[, 2]),na.rm = TRUE)
    diff <- ci[, 2] - ci[, 1]
    len <- mean(diff[is.finite(diff)])
    n_inf <- sum(is.infinite(diff))

    print(paste0(print_list[i], " coverage, ",coverage, ', lens ', len))
    record[[i]][trial,] <- c(coverage,len,n_inf)
  }
  print("#################")
}


##----------------------------------------------------------------
##                         Save results                         --
##----------------------------------------------------------------

#create a new path for files

folder<- paste0(path, dtype)
if(fct <1){
  folder <- paste0(folder,"/", "fct_", fct,"/")
}else{
  folder<- paste0(folder, "/","gmm_tmp",gmm_star,"/")
}
print(folder)
dir.create(folder, recursive=TRUE, showWarnings = FALSE)

#coverage data
coverage <-c()
for (i in 1:length(print_list)){coverage[[i]]<- as.vector(record[[i]][,1])}

data <- data.frame(Coverage=unlist(coverage),
                   group=rep(c("CSA-M","CSA-Q","ITE-NUS", "CSA-B"),
                             each=ntrial))

if(save){
  write.csv(data, paste0(folder,'coverage', '.csv'), row.names = FALSE)
}

#length data
Interval_length <-c()
for (i in 1:length(print_list)){Interval_length[[i]]<- as.vector(record[[i]][,2])}
data <- data.frame(Interval_length= unlist(Interval_length),
                   group=rep(c("CSA-M","CSA-Q","ITE-NUS", "CSA-B"),each=ntrial))

if(save){
  write.csv(data, paste0(folder,'len','.csv'), row.names = FALSE)
}

