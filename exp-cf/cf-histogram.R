library("ggplot2")
setwd('~/Desktop/Working/conformal-sensitivity/src')
library("devtools")
devtools::load_all(".")
rm(list = ls())
gmm <- 4

##sample test data
Xfun <- function(n, d){
  matrix(runif(n * d), nrow = n, ncol = d)
}
taufun <- function(X){
  2 / (1 + exp(-5 * (X[, 1] - 0.5))) * 2 / (1 + exp(-5 * (X[, 2] - 0.5)))
}
sdfun <- function(X){
  rep(1, nrow(X))
}
pscorefun <- function(X){
  1-(1 + pbeta(1-X[, 1], 2, 4)) / 4
}
errdist <- rnorm

ntest <- 100000
d <- 20
set.seed(1)
Xtest <- Xfun(1,d)
Xtest <- Xtest[rep(1, times = ntest), ]
#plot
Yobs <- taufun(Xtest) + sdfun(Xtest) * errdist(dim(Xtest)[1])

for(case in 1:3){
  if(case==1){
    rg <- t(as.matrix(c(0,2)))
    rg <- rg[rep(1, times = ntest), ]
    Ycf <- samplecf(Xtest, taufun, sdfun, case=case, area=rg, gmm=gmm)
  }else{
    Ycf <- samplecf(Xtest, taufun, sdfun, case=case, gmm=gmm)
  }

  data<-data.frame(y=c(Ycf,Yobs),group=rep(c("Cf","Obs"),each=length(Yobs)))
  p <- ggplot(data,aes(x=y,fill=group))+geom_density(alpha=0.2) +
       theme(legend.position=c(0.9,0.8),legend.key.size = unit(0.85, 'cm'),
             legend.text = element_text(size=14),
             #legend.title = element_text(size=14),
             legend.title = element_blank(),
             axis.text = element_text(size=16),
             axis.title=element_text(size=16))+
       coord_cartesian(xlim=c(-4, 4.5),ylim=c(0, 0.7))+
       labs(x="Y(1)")
  print(p)
}



