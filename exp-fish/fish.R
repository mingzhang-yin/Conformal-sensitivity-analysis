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

suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()
parser$add_argument("--gmm_star", type = "double", default = 3, help = "SA parameter, >=1")
parser$add_argument("--alpha", type="double", default=0.2, help="miscoverage")

parser$add_argument("--save", type="logical", default=TRUE, help="save")
parser$add_argument("--seed", type = "double", default = 1, help = "random seed")
parser$add_argument("--ntrial", type = "integer", default = 50, help = "number of trials")
parser$add_argument("--path", type = "character", default = '/proj/sml_netapp/projects/claudia/conformal/fish/', help = "save location")

args <- parser$parse_args()
alpha <- args$alpha
gmm_star <- args$gmm_star
ntrial<- args$ntrial
seed <- args$seed
save <- args$save
path = args$path
q<- c(alpha/2, 1- (alpha/2))


##---------------------------------------------------------------
##      loading observed data and pre-processing               -
##---------------------------------------------------------------

load('./data/fish.Rda')
A <- as.numeric(nhanes.fish$fish.level == "high")
X <- nhanes.fish[, c("gender", "age", "income", "income.missing", "race", "education", "smoking.ever", "smoking.now")]
X$race <- factor(X$race)
X1 <- model.matrix(~ . - 1, X)
Y_all <- log2(nhanes.fish$o.LBXTHG)

record <- replicate(2,matrix(0,nrow=ntrial,ncol=3), simplify=FALSE)


folder<- paste0(path, 'alpha_',alpha,'/','gmm_',gmm_star, '/')
dir.create(folder, recursive=TRUE, showWarnings = FALSE)

for (iter in 1:ntrial){
  n<- length(Y_all)
  trainprop <- 0.8
  set.seed(123)
  trainid <- sample(n, floor(n * trainprop))
  set.seed(NULL)
  print(paste0("alpha is ",alpha))
  ##---------------------------------------------------------------
  ##                     splitting                    -
  ##---------------------------------------------------------------


  Y_obs <- Y_all[trainid]
  X <- X1[trainid,]
  T_obs <- A[trainid]
  Y1 <- Y_obs
  Y1[which(T_obs==0)] <- NA
  Y0 <- Y_obs
  Y0[which(T_obs==1)] <- NA

  id <- seq(1, n)
  testid<- id[!(id %in% trainid)]
  Xtest <- X1[testid,]


  ##----------------------------------------------------------------
  ##            getting prediction bands            --
  ##----------------------------------------------------------------



  obj_mean <- nested_conformalSA(X, Y1, Y0, T_obs, gmm_star, type = "mean",quantiles=list(), outfun='RF')
  obj_cqr <- nested_conformalSA(X, Y1, Y0, T_obs, gmm_star, type = "CQR",quantiles=q, outfun='quantRF')

  obj_bands_mean <- predict.nested(obj_mean, X, Y_obs, T_obs, alpha = alpha)
  obj_bands_cqr <- predict.nested(obj_cqr, X, Y_obs, T_obs, alpha = alpha)


  ##----------------------------------------------------------------
  ##                       ITE & evaluation                    -
  ##----------------------------------------------------------------


  # Inexact-SA
  ci_mean <- fit_and_predict_band(obj_bands_mean, Xtest,'RF')
  ci_cqr <- fit_and_predict_band(obj_bands_cqr, Xtest,'RF')

  ci_list <- list(ci_mean,ci_cqr)
  print_list <- list("ci_mean", "ci_cqr")
  data <- cbind(ci_mean, ci_cqr)
  colnames(data) <- c("mean_low", "mean_high", "cqr_low", "cqr_high")
  df <- as.data.frame(t(data))
  write.csv(data, file=paste0(folder, 'ntrial_', iter, '.csv'))


  for(i in 1:length(ci_list)){
    ci <- ci_list[[i]]
    negative <- sum(ci[, 2] <=0)
    positive <- sum(ci[, 1] >=0)
    neg_percentage <- negative/dim(ci)[1]
    pos_percentage <- positive/dim(ci)[1]
    diff <- ci[, 2] - ci[, 1]
    len <-mean(diff[is.finite(diff)])
    print(print_list[i])
    print(paste0( 'neg_percentage: ', neg_percentage,
                  ' pos_percentage: ', pos_percentage, "  len: ", len))
    print("###############")
    record[[i]][iter,] <- c(neg_percentage,pos_percentage,len)
  }

}






