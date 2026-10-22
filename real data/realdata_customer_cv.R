
rm(list=ls())
gc()

setwd("D:\\GAO Junzhuo\\ADMM_Huber_decen_vert")

library(caret)

thres <- function(a, t) {
  sign(a) * pmax(abs(a) - t, 0)
}

trunc <- function(a, tau) {
  pmin(pmax(a, -tau), tau)
}

loss <- function(x, tau) {
  (abs(x) <= tau) * x ^ 2 / 2 + (abs(x) > tau) * (tau * abs(x) - tau ^ 2 / 2)
}

mydata <- read.csv("marketing_campaign.csv")
mydata <- na.omit(mydata)
mydata <- mydata[, -c(1, 2, 3, 4, 8, 27, 28)]
mylist <- c("Income", "Recency", "MntWines", "MntFruits", "MntMeatProducts", "MntFishProducts",
            "MntSweetProducts", "MntGoldProds", "NumWebPurchases", "NumCatalogPurchases", "NumStorePurchases",
            "NumWebVisitsMonth")
for (i in 1:length(mylist)) {
  mydata[, mylist[i]] <- scale(mydata[, mylist[i]], center=T, scale=T)
}
mydata[, "NumDealsPurchases"] <- ifelse(mydata[, "NumDealsPurchases"] > 0,
                                        log(mydata[, "NumDealsPurchases"]),
                                        mydata[, "NumDealsPurchases"])
# mydata[, "NumDealsPurchases"] <- sqrt(mydata[, "NumDealsPurchases"])
# mydata[, "NumDealsPurchases"] <- scale(mydata[, "NumDealsPurchases"], center=T, scale=T)

n.train <- 1500
n.test <- dim(mydata)[1] - n.train
d <- dim(mydata)[2]

ga <- 0.5
m <- 4
A <- matrix(c(1, 1, 0, 1,
              1, 1, 1, 0,
              0, 1, 1, 1,
              1, 0, 1, 1), ncol=4)
dj <- c(5 + 1, 6, 4, 6)
djk <- c(0, 5 + 1, 11 + 1, 15 + 1, 21 + 1)

set.seed(1)
id.train <- sample(1:dim(mydata)[1], n.train, replace=F)
x.train <- mydata[id.train, c("Income", "Kidhome", "Teenhome", "Recency", "Complain",
                              "MntWines", "MntFruits", "MntMeatProducts", "MntFishProducts", "MntSweetProducts", "MntGoldProds",
                              "NumWebPurchases", "NumCatalogPurchases", "NumStorePurchases", "NumWebVisitsMonth",
                              "AcceptedCmp1", "AcceptedCmp2", "AcceptedCmp3", "AcceptedCmp4", "AcceptedCmp5", "Response")]
x.train <- cbind(1, as.matrix(x.train))
x.test <- mydata[-c(id.train), c("Income", "Kidhome", "Teenhome", "Recency", "Complain",
                                 "MntWines", "MntFruits", "MntMeatProducts", "MntFishProducts", "MntSweetProducts", "MntGoldProds",
                                 "NumWebPurchases", "NumCatalogPurchases", "NumStorePurchases", "NumWebVisitsMonth",
                                 "AcceptedCmp1", "AcceptedCmp2", "AcceptedCmp3", "AcceptedCmp4", "AcceptedCmp5", "Response")]
x.test <- cbind(1, as.matrix(x.test))
y.train <- mydata[id.train, "NumDealsPurchases"]

########## contamination ##########
set.seed(1)
id.cont <- sample(1:n.train, n.train * 0.1, replace=F)
y.train[id.cont] <- y.train[id.cont] - 1.5
###################################

y.test <- mydata[-c(id.train), "NumDealsPurchases"]

xm.train <- array(NA, c(n.train, m, max(dj)))
for (j in 1:m) {
  xm.train[, j, 1:dj[j]] <- x.train[, (djk[j] + 1):djk[j + 1]]
}

# VFLHuber
## cross validation
tau <- 1 # sqrt(n.train / log(d)) / 25
k.fold <- 5
itr.max <- 500
cc <- 20
p1.set <- c(0.05, 0.10, 0.15)
p2.set <- c(0.05, 0.10, 0.15)
err <- array(NA, c(k.fold, length(p1.set), length(p2.set)))
err.mean <- matrix(0, length(p1.set), length(p2.set))
set.seed(1)
fold <- createFolds(1:n.train, k=k.fold, list=T)
for (k in 1:k.fold) {
  cat(k)
  for (p1.idx in 1:length(p1.set)) {
    p1 <- p1.set[p1.idx]
    for (p2.idx in 1:length(p2.set)) {
      p2 <- p2.set[p2.idx]
      b <- matrix(0, itr.max + 1, d)
      u <- v <- a <- array(0, c(itr.max + 1, m, n.train * (1 - 1 / k.fold)))
      for (itr in 1:itr.max) {
        for (j in 1:m) {
          b[itr + 1, (djk[j] + 1):djk[j + 1]] <-
            thres(-t(xm.train[-as.vector(unlist(fold[k])), j, 1:dj[j]]) %*% u[itr, j, ] / (2 * n.train * (1 - 1 / k.fold) * p2), p1 / (2 * p2))
        }
        for (j in 1:m) {
          su <- sum(A[j, ]) * u[itr, j, ] - apply(u[itr, A[j, ]!=0, ], 2, sum)
          v[itr + 1, j, ] <- u[itr, j, ] - (-xm.train[-as.vector(unlist(fold[k])), j, 1:dj[j]] %*% b[itr + 1, (djk[j] + 1):djk[j + 1]] + ga * su + a[itr, j, ]) / cc
        }
        for (j in 1:m) {
          sv <- sum(A[j, ]) * v[itr + 1, j, ] - apply(v[itr + 1, A[j, ]!=0, ], 2, sum)
          a[itr + 1, j, ] <- a[itr, j, ] + ga * sv
          u[itr + 1, j, ] <- trunc((v[itr + 1, j, ] - y.train[-as.vector(unlist(fold[k]))] / (m * cc)) / (1 + 1 / (m * cc)), tau)
        }
      }
      x.train.test <- c()
      for (j in 1:m) {
        x.train.test <- cbind(x.train.test, xm.train[as.vector(unlist(fold[k])), j, 1:dj[j]])
      }
      # err[k, p1.idx, p2.idx] <- mean((y.train[as.vector(unlist(fold[k]))] - x.train.test %*% b[itr + 1, ]) ^ 2)
      err[k, p1.idx, p2.idx] <- mean(loss(y.train[as.vector(unlist(fold[k]))] - x.train.test %*% b[itr + 1, ], tau))
    }
  }
  err.mean <- err.mean + err[k, , ] / k.fold
}
p1 <- p1.set[which(err.mean==min(err.mean), arr.ind=T)[1]]
p2 <- p2.set[which(err.mean==min(err.mean), arr.ind=T)[2]]
p1
p2

# tau <- 1
# p1 <- 0.05
# p2 <- 0.05
# itr.max <- 500
# cc <- 20
b.hub <- matrix(0, itr.max + 1, d)
mse.hub <- med.hub <- loss.hub <- rep(NA, itr.max + 1)
mse.hub[1] <- mean((x.test %*% c(b.hub[1, ]) - y.test) ^ 2)
med.hub[1] <- mean(abs(x.test %*% c(b.hub[1, ]) - y.test))
loss.hub[1] <- mean(loss(y.train - x.train %*% b.hub[1, ], tau))
u <- v <- a <- array(0, c(itr.max + 1, m, n.train))
for (itr in 1:itr.max) {
  if (itr %% 100==0) {
    cat(itr)
  }
  for (j in 1:m) {
    b.hub[itr + 1, (djk[j] + 1):djk[j + 1]] <- 
      thres(-t(xm.train[, j, 1:dj[j]]) %*% u[itr, j, ] / (2 * n.train * p2), p1 / (2 * p2))
  }
  mse.hub[itr + 1] <- mean((x.test %*% b.hub[itr + 1, ] - y.test) ^ 2)
  med.hub[itr + 1] <- mean(abs(x.test %*% b.hub[itr + 1, ] - y.test))
  loss.hub[itr + 1] <- mean(loss(y.train - x.train %*% b.hub[itr + 1, ], tau))
  for (j in 1:m) {
    su <- sum(A[j, ]) * u[itr, j, ] - apply(u[itr, A[j, ]!=0, ], 2, sum)
    v[itr + 1, j, ] <- u[itr, j, ] - 
      (-xm.train[, j, 1:dj[j]] %*% b.hub[itr + 1, (djk[j] + 1):djk[j + 1]] + ga * su + a[itr, j, ]) / cc
  }
  for (j in 1:m) {
    sv <- sum(A[j, ]) * v[itr + 1, j, ] - apply(v[itr + 1, A[j, ]!=0, ], 2, sum)
    a[itr + 1, j, ] <- a[itr, j, ] + ga * sv
    u[itr + 1, j, ] <- trunc((v[itr + 1, j, ] - y.train / (m * cc)) / (1 + 1 / (m * cc)), tau)
  }
}
plot(loss.hub)
plot(mse.hub)
plot(med.hub)
mse.hub[itr.max + 1]
med.hub[itr.max + 1]
b.hub[itr.max + 1, ]

# CenHuber
cc <- 20
b.cen <- matrix(0, itr.max + 1, d)
mse.cen <- med.cen <- loss.cen <- rep(NA, itr.max + 1)
mse.cen[1] <- mean((x.test %*% b.cen[1, ] - y.test) ^ 2)
med.cen[1] <- mean(abs(x.test %*% b.cen[1, ] - y.test))
loss.cen[1] <- mean(loss(y.train - x.train %*% b.cen[1, ], tau))
x.cen <- z.cen <- v.cen <- matrix(0, itr.max + 1, n.train)
for (itr in 1:itr.max) {
  if (itr %% 100==0) {
    cat(itr)
  }
  b.cen[itr + 1, ] <- thres(-t(x.train) %*% x.cen[itr, ] / (2 * n.train * p2), p1 / (2 * p2))
  mse.cen[itr + 1] <- mean((x.test %*% b.cen[itr + 1, ] - y.test) ^ 2)
  med.cen[itr + 1] <- mean(abs(x.test %*% b.cen[itr + 1, ] - y.test))
  loss.cen[itr + 1] <- mean(loss(y.train - x.train %*% b.cen[itr + 1, ], tau))
  x.cen[itr + 1, ] <- (1 - ga / cc) * x.cen[itr, ] + x.train %*% b.cen[itr + 1, ] / cc -
    v.cen[itr, ] / cc + ga * z.cen[itr, ] / cc
  z.cen[itr + 1, ] <- trunc((x.cen[itr + 1, ] + v.cen[itr, ] / ga - y.train / ga) / (1 + 1 / ga), tau)
  v.cen[itr + 1, ] <- v.cen[itr, ] + ga * (x.cen[itr + 1, ] - z.cen[itr + 1, ])
}
plot(loss.cen)
plot(mse.cen)
plot(med.cen)
mse.cen[itr.max + 1]
med.cen[itr.max + 1]
b.cen[itr.max + 1, ]

# qqnorm(y.train - x.train %*% b.cen[itr.max + 1, ], main="QQ Plot: Comparison with Normal Distribution")
# abline(0, 1, col="red")
# hist(y.train - x.train %*% b.cen[itr.max + 1, ])
# seq1 <- y.train - x.train %*% b.cen[itr.max + 1, ]
# seq2 <- y.train - x.train %*% b.cen[itr.max + 1, ]
# dat <- data.frame(Category=rep(c("CenHuber", "CenHuber (10% contamination)"), each=length(seq1)),
#                   Value=c(seq1, seq2))
# boxplot(Value~Category, data=dat,
#         xlab="",
#         ylab="Residuals",
#         col=c("lightblue", "lightcoral"),
#         border="black",
#         cex.axis=2,
#         cex.lab=1.5)

# LocHuber
cc <- 20
b.loc <- matrix(0, itr.max + 1, d)
for (j in 1:m) {
  cat(j)
  x.loc <- z.loc <- v.loc <- matrix(0, itr.max + 1, n.train)
  for (itr in 1:itr.max) {
    b.loc[itr + 1, (djk[j] + 1):djk[j + 1]] <- thres(-t(xm.train[, j, 1:dj[j]]) %*% x.loc[itr, ] / (2 * n.train * p2), p1 / (2 * p2))
    x.loc[itr + 1, ] <- (1 - ga / cc) * x.loc[itr, ] + xm.train[, j, 1:dj[j]] %*% b.loc[itr + 1, (djk[j] + 1):djk[j + 1]] / cc -
      v.loc[itr, ] / cc + ga * z.loc[itr, ] / cc
    z.loc[itr + 1, ] <- trunc((x.loc[itr + 1, ] + v.loc[itr, ] / ga - y.train / ga) / (1 + 1 / ga), tau)
    v.loc[itr + 1, ] <- v.loc[itr, ] + ga * (x.loc[itr + 1, ] - z.loc[itr + 1, ])
  }
}
mse.loc <- med.loc <- loss.loc <- rep(NA, itr.max + 1)
for (itr in 0:itr.max) {
  mse.loc[itr + 1] <- mean((x.test %*% c(b.loc[itr + 1, ]) - y.test) ^ 2)
  med.loc[itr + 1] <- mean(abs(x.test %*% c(b.loc[itr + 1, ]) - y.test))
  loss.loc[itr + 1] <- mean(loss(y.train - x.train %*% b.loc[itr + 1, ], tau))
}
plot(loss.loc)
plot(mse.loc)
plot(med.loc)
mse.loc[itr.max + 1]
med.loc[itr.max + 1]
b.loc[itr.max + 1, ]

c(mse.cen[length(mse.cen)],
  mse.hub[length(mse.hub)],
  mse.loc[length(mse.loc)])

c(med.cen[length(med.cen)],
  med.hub[length(med.hub)],
  med.loc[length(med.loc)])

# VFLMed
## cross validation
# k.fold <- 5
# itr.max <- 1500
# cc <- 20
# p1.set <- c(0.05, 0.10, 0.15)
# p2.set <- c(0.05, 0.10, 0.15)
# err <- array(NA, c(k.fold, length(p1.set), length(p2.set)))
# err.mean <- matrix(0, length(p1.set), length(p2.set))
# for (k in 1:k.fold) {
#   cat(k)
#   for (p1.idx in 1:length(p1.set)) {
#     p1 <- p1.set[p1.idx]
#     for (p2.idx in 1:length(p2.set)) {
#       p2 <- p2.set[p2.idx]
#       b <- matrix(0, itr.max + 1, d)
#       u <- v <- a <- array(0, c(itr.max + 1, m, n.train * (1 - 1 / k.fold)))
#       for (itr in 1:itr.max) {
#         for (j in 1:m) {
#           b[itr + 1, (djk[j] + 1):djk[j + 1]] <-
#             thres(-t(xm.train[-as.vector(unlist(fold[k])), j, 1:dj[j]]) %*% u[itr, j, ] / (2 * n.train * (1 - 1 / k.fold) * p2), p1 / (2 * p2))
#         }
#         for (j in 1:m) {
#           su <- sum(A[j, ]) * u[itr, j, ] - apply(u[itr, A[j, ]!=0, ], 2, sum)
#           v[itr + 1, j, ] <- u[itr, j, ] - (-xm.train[-as.vector(unlist(fold[k])), j, 1:dj[j]] %*% b[itr + 1, (djk[j] + 1):djk[j + 1]] + ga * su + a[itr, j, ]) / cc
#         }
#         for (j in 1:m) {
#           sv <- sum(A[j, ]) * v[itr + 1, j, ] - apply(v[itr + 1, A[j, ]!=0, ], 2, sum)
#           a[itr + 1, j, ] <- a[itr, j, ] + ga * sv
#           u[itr + 1, j, ] <- trunc(v[itr + 1, j, ] - y.train[-as.vector(unlist(fold[k]))] / (m * cc), 1)
#         }
#       }
#       # err[k, p1.idx, p2.idx] <- mean((y.train[as.vector(unlist(fold[k]))] - x.train.test %*% b[itr + 1, ]) ^ 2)
#       err[k, p1.idx, p2.idx] <- mean(abs(y.train[as.vector(unlist(fold[k]))] - x.train.test %*% b[itr + 1, ]))
#     }
#   }
#   err.mean <- err.mean + err[k, , ] / k.fold
# }
# p1 <- p1.set[which(err.mean==min(err.mean), arr.ind=T)[1]]
# p2 <- p2.set[which(err.mean==min(err.mean), arr.ind=T)[2]]
# p1
# p2

# p1 <- 0.10
# p2 <- 0.10
# itr.max <- 1500
# cc <- 20
# b.med <- matrix(0, itr.max + 1, d)
# mse.med <- med.med <- loss.med <- rep(NA, itr.max + 1)
# mse.med[1] <- mean((x.test %*% c(b.med[1, ]) - y.test) ^ 2)
# med.med[1] <- mean(abs(x.test %*% c(b.med[1, ]) - y.test))
# loss.med[1] <- mean(abs(y.train - x.train %*% b.med[1, ]))
# u <- v <- a <- array(0, c(itr.max + 1, m, n.train))
# for (itr in 1:itr.max) {
#   if (itr %% 100==0) {
#     cat(itr)
#   }
#   for (j in 1:m) {
#     b.med[itr + 1, (djk[j] + 1):djk[j + 1]] <- thres(-t(xm.train[, j, 1:dj[j]]) %*% u[itr, j, ] / (2 * n.train * p2), p1 / (2 * p2))
#   }
#   mse.med[itr + 1] <- mean((x.test %*% b.med[itr + 1, ] - y.test) ^ 2)
#   med.med[itr + 1] <- mean(abs(x.test %*% b.med[itr + 1, ] - y.test))
#   loss.med[itr + 1] <- mean(abs(y.train - x.train %*% b.med[itr + 1, ]))
#   for (j in 1:m) {
#     su <- sum(A[j, ]) * u[itr, j, ] - apply(u[itr, A[j, ]!=0, ], 2, sum)
#     v[itr + 1, j, ] <- u[itr, j, ] - (-xm.train[, j, 1:dj[j]] %*% b.med[itr + 1, (djk[j] + 1):djk[j + 1]] + ga * su + a[itr, j, ]) / cc
#   }
#   for (j in 1:m) {
#     sv <- sum(A[j, ]) * v[itr + 1, j, ] - apply(v[itr + 1, A[j, ]!=0, ], 2, sum)
#     a[itr + 1, j, ] <- a[itr, j, ] + ga * sv
#     u[itr + 1, j, ] <- trunc(v[itr + 1, j, ] - y.train / (m * cc), 1)
#   }
# }
# plot(loss.med)
# plot(mse.med)
# plot(med.med)
# mse.med[itr.max + 1]
# med.med[itr.max + 1]
# b.med[itr.max + 1, ]

# VFLLS
## cross validation
k.fold <- 5
itr.max <- 500
cc <- 20
p1.set <- c(0.05, 0.10, 0.15)
p2.set <- c(0.05, 0.10, 0.15)
err <- array(NA, c(k.fold, length(p1.set), length(p2.set)))
err.mean <- matrix(0, length(p1.set), length(p2.set))
for (k in 1:k.fold) {
  cat(k)
  for (p1.idx in 1:length(p1.set)) {
    p1 <- p1.set[p1.idx]
    for (p2.idx in 1:length(p2.set)) {
      p2 <- p2.set[p2.idx]
      b <- matrix(0, itr.max + 1, d)
      u <- v <- a <- array(0, c(itr.max + 1, m, n.train * (1 - 1 / k.fold)))
      for (itr in 1:itr.max) {
        for (j in 1:m) {
          b[itr + 1, (djk[j] + 1):djk[j + 1]] <-
            thres(-t(xm.train[-as.vector(unlist(fold[k])), j, 1:dj[j]]) %*% u[itr, j, ] / (2 * n.train * (1 - 1 / k.fold) * p2), p1 / (2 * p2))
        }
        for (j in 1:m) {
          su <- sum(A[j, ]) * u[itr, j, ] - apply(u[itr, A[j, ]!=0, ], 2, sum)
          v[itr + 1, j, ] <- u[itr, j, ] - (-xm.train[-as.vector(unlist(fold[k])), j, 1:dj[j]] %*% b[itr + 1, (djk[j] + 1):djk[j + 1]] + ga * su + a[itr, j, ]) / cc
        }
        for (j in 1:m) {
          sv <- sum(A[j, ]) * v[itr + 1, j, ] - apply(v[itr + 1, A[j, ]!=0, ], 2, sum)
          a[itr + 1, j, ] <- a[itr, j, ] + ga * sv
          u[itr + 1, j, ] <- (v[itr + 1, j, ] - y.train[-as.vector(unlist(fold[k]))] / (m * cc)) / (1 + 1 / (m * cc))
        }
      }
      err[k, p1.idx, p2.idx] <- mean((y.train[as.vector(unlist(fold[k]))] - x.train.test %*% b[itr + 1, ]) ^ 2)
    }
  }
  err.mean <- err.mean + err[k, , ] / k.fold
}
p1 <- p1.set[which(err.mean==min(err.mean), arr.ind=T)[1]]
p2 <- p2.set[which(err.mean==min(err.mean), arr.ind=T)[2]]
p1
p2

# p1 <- 0.10
# p2 <- 0.10
# itr.max <- 500
# cc <- 20
b.ols <- matrix(0, itr.max + 1, d)
mse.ols <- med.ols <- loss.ols <- rep(NA, itr.max + 1)
mse.ols[1] <- mean((x.test %*% c(b.ols[1, ]) - y.test) ^ 2)
med.ols[1] <- mean(abs(x.test %*% c(b.ols[1, ]) - y.test))
loss.ols[1] <- mean(loss(y.train - x.train %*% b.ols[1, ], tau))
u <- v <- a <- array(0, c(itr.max + 1, m, n.train))
for (itr in 1:itr.max) {
  if (itr %% 100==0) {
    cat(itr)
  }
  for (j in 1:m) {
    b.ols[itr + 1, (djk[j] + 1):djk[j + 1]] <- thres(-t(xm.train[, j, 1:dj[j]]) %*% u[itr, j, ] / (2 * n.train * p2), p1 / (2 * p2))
  }
  mse.ols[itr + 1] <- mean((x.test %*% b.ols[itr + 1, ] - y.test) ^ 2)
  med.ols[itr + 1] <- mean(abs(x.test %*% b.ols[itr + 1, ] - y.test))
  loss.ols[itr + 1] <- mean(loss(y.train - x.train %*% b.ols[itr + 1, ], tau))
  for (j in 1:m) {
    su <- sum(A[j, ]) * u[itr, j, ] - apply(u[itr, A[j, ]!=0, ], 2, sum)
    v[itr + 1, j, ] <- u[itr, j, ] - (-xm.train[, j, 1:dj[j]] %*% b.ols[itr + 1, (djk[j] + 1):djk[j + 1]] + ga * su + a[itr, j, ]) / cc
  }
  for (j in 1:m) {
    sv <- sum(A[j, ]) * v[itr + 1, j, ] - apply(v[itr + 1, A[j, ]!=0, ], 2, sum)
    a[itr + 1, j, ] <- a[itr, j, ] + ga * sv
    u[itr + 1, j, ] <- (v[itr + 1, j, ] - y.train / (m * cc)) / (1 + 1 / (m * cc))
  }
}
plot(loss.ols)
plot(mse.ols)
plot(med.ols)
mse.ols[itr.max + 1]
med.ols[itr.max + 1]
b.ols[itr.max + 1, ]

c(mse.cen[length(mse.cen)],
  mse.hub[length(mse.hub)],
  mse.loc[length(mse.loc)],
  #mse.med[length(mse.med)],
  mse.ols[length(mse.ols)])

c(med.cen[length(med.cen)],
  med.hub[length(med.hub)],
  med.loc[length(med.loc)],
  #med.med[length(med.med)],
  med.ols[length(med.ols)])

round(b.cen[length(mse.cen), ], 3)
round(b.hub[length(mse.hub), ], 3)
round(b.loc[length(mse.loc), ], 3)
round(b.med[length(mse.med), ], 3)
round(b.ols[length(mse.ols), ], 3)

# library(ggplot2)
# data <- data.frame(value=mydata[, "NumDealsPurchases"])
# ggplot(data, aes(x=value)) +
#   geom_histogram(aes(y=..count..), bins=6, fill="lightblue", color="black", alpha=0.7) +
#   theme_minimal() +
#   labs(title="Histogram with Frequency", x="NumDealsPurchases", y="Frequency") +
#   theme(plot.title=element_text(hjust=0.5, size=30, face="bold"),
#         axis.title=element_text(size=30),
#         axis.text=element_text(size=30))

