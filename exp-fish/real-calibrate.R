Aall <- as.numeric(nhanes.fish$fish.level == "high")
Xall <- nhanes.fish[, c("gender", "age", "income", "income.missing",
"race", "education", "smoking.ever", "smoking.now")]
Xall$race <- factor(Xall$race)
X1 <- model.matrix(~ . - 1, Xall)
Yall <- log2(nhanes.fish$o.LBXTHG)

psfun <- Boosting

n<- length(Yall)
trainprop <- 0.8
set.seed(123)
trainid <- sample(n, floor(n * trainprop))
set.seed(NULL)

X <- X1[trainid,]
Y <- Yall[trainid]
A <- Aall[trainid]
Xtest <- X1[-trainid,]

ps <- do.call(psfun, list(Xtest = Xtest, Y = A , X = X))

record <- c()
for(j in 1:dim(X)[2]){
  ps_j <- do.call(psfun, list(Xtest = Xtest[,-c(j)], Y = A , X = X[,-c(j)]))
  odds_j <- (ps/(1-ps))/(ps_j/(1-ps_j))
  record <- c(record,odds_j)
}

idx <- which(record>=1)
R <- record
R[idx] <- record[idx]
R[-idx] <- 1/record[-idx]

quantile(R, seq(0, 1, by=.05))
summary(R)
boxplot(R)

pdf("../doc/fig/fish/calib.pdf",width=5.5,height=6)
hist(R, breaks = 10, prob=TRUE, xlab='(Adjusted) odds ratio',main=NULL,
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
dens <- density(R, from=1, width=1.5)
lines(dens)
par(mar = c(4, 4, 0.1, 0.1))
dev.off()
