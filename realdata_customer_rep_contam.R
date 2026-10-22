
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

rep.max <- 100
itr.max <- 500
b.cen.rep <- b.hub.rep <- b.loc.rep <- b.ols.rep <- matrix(NA, rep.max, 22)
mse.cen.rep <- mse.hub.rep <- mse.loc.rep <- mse.ols.rep <- rep(NA, rep.max)
mae.cen.rep <- mae.hub.rep <- mae.loc.rep <- mae.ols.rep <- rep(NA, rep.max)

for (rep.idx in 1:rep.max) {
  cat(rep.idx)
  set.seed(rep.idx)
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
  set.seed(rep.idx)
  id.cont <- sample(1:n.train, n.train * 0.1, replace=F)
  y.train[id.cont] <- y.train[id.cont] - 5
  ###################################
  
  y.test <- mydata[-c(id.train), "NumDealsPurchases"]

  xm.train <- array(NA, c(n.train, m, max(dj)))
  for (j in 1:m) {
    xm.train[, j, 1:dj[j]] <- x.train[, (djk[j] + 1):djk[j + 1]]
  }

  # VFLHuber
  tau <- 1 # sqrt(n.train / log(d)) / 25
  cc <- 20
  p1 <- 0.05
  p2 <- 0.05
  b.hub <- matrix(0, itr.max + 1, d)
  mse.hub <- med.hub <- loss.hub <- rep(NA, itr.max + 1)
  mse.hub[1] <- mean((x.test %*% c(b.hub[1, ]) - y.test) ^ 2)
  med.hub[1] <- mean(abs(x.test %*% c(b.hub[1, ]) - y.test))
  loss.hub[1] <- mean(loss(y.train - x.train %*% b.hub[1, ], tau))
  u <- v <- a <- array(0, c(itr.max + 1, m, n.train))
  for (itr in 1:itr.max) {
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
  mse.hub.rep[rep.idx] <- mse.hub[itr.max + 1]
  mae.hub.rep[rep.idx] <- med.hub[itr.max + 1]
  b.hub.rep[rep.idx, ] <- b.hub[itr.max + 1, ]
  
  # CenHuber
  cc <- 20
  b.cen <- matrix(0, itr.max + 1, d)
  mse.cen <- med.cen <- loss.cen <- rep(NA, itr.max + 1)
  mse.cen[1] <- mean((x.test %*% b.cen[1, ] - y.test) ^ 2)
  med.cen[1] <- mean(abs(x.test %*% b.cen[1, ] - y.test))
  loss.cen[1] <- mean(loss(y.train - x.train %*% b.cen[1, ], tau))
  x.cen <- z.cen <- v.cen <- matrix(0, itr.max + 1, n.train)
  for (itr in 1:itr.max) {
    b.cen[itr + 1, ] <- thres(-t(x.train) %*% x.cen[itr, ] / (2 * n.train * p2), p1 / (2 * p2))
    mse.cen[itr + 1] <- mean((x.test %*% b.cen[itr + 1, ] - y.test) ^ 2)
    med.cen[itr + 1] <- mean(abs(x.test %*% b.cen[itr + 1, ] - y.test))
    loss.cen[itr + 1] <- mean(loss(y.train - x.train %*% b.cen[itr + 1, ], tau))
    x.cen[itr + 1, ] <- (1 - ga / cc) * x.cen[itr, ] + x.train %*% b.cen[itr + 1, ] / cc -
      v.cen[itr, ] / cc + ga * z.cen[itr, ] / cc
    z.cen[itr + 1, ] <- trunc((x.cen[itr + 1, ] + v.cen[itr, ] / ga - y.train / ga) / (1 + 1 / ga), tau)
    v.cen[itr + 1, ] <- v.cen[itr, ] + ga * (x.cen[itr + 1, ] - z.cen[itr + 1, ])
  }
  mse.cen.rep[rep.idx] <- mse.cen[itr.max + 1]
  mae.cen.rep[rep.idx] <- med.cen[itr.max + 1]
  b.cen.rep[rep.idx, ] <- b.cen[itr.max + 1, ]
  
  # LocHuber
  cc <- 20
  b.loc <- matrix(0, itr.max + 1, d)
  for (j in 1:m) {
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
  mse.loc.rep[rep.idx] <- mse.loc[itr.max + 1]
  mae.loc.rep[rep.idx] <- med.loc[itr.max + 1]
  b.loc.rep[rep.idx, ] <- b.loc[itr.max + 1, ]

  # VFLLS
  p1 <- 0.10
  p2 <- 0.10
  cc <- 20
  b.ols <- matrix(0, itr.max + 1, d)
  mse.ols <- med.ols <- loss.ols <- rep(NA, itr.max + 1)
  mse.ols[1] <- mean((x.test %*% c(b.ols[1, ]) - y.test) ^ 2)
  med.ols[1] <- mean(abs(x.test %*% c(b.ols[1, ]) - y.test))
  loss.ols[1] <- mean(loss(y.train - x.train %*% b.ols[1, ], tau))
  u <- v <- a <- array(0, c(itr.max + 1, m, n.train))
  for (itr in 1:itr.max) {
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
  mse.ols.rep[rep.idx] <- mse.ols[itr.max + 1]
  mae.ols.rep[rep.idx] <- med.ols[itr.max + 1]
  b.ols.rep[rep.idx, ] <- b.ols[itr.max + 1, ]
}

save(b.cen.rep, b.hub.rep, b.loc.rep, b.ols.rep,
     mse.cen.rep, mse.hub.rep, mse.loc.rep, mse.ols.rep,
     mae.cen.rep, mae.hub.rep, mae.loc.rep, mae.ols.rep,
     file = "realdata_customer_rep_contam.RData")

# load("realdata_customer_rep_contam.RData")
round(c(mean(mse.cen.rep), mean(mse.hub.rep), mean(mse.loc.rep), mean(mse.ols.rep)), 3)
round(c(sd(mse.cen.rep), sd(mse.hub.rep), sd(mse.loc.rep), sd(mse.ols.rep)), 3)
round(c(mean(mae.cen.rep), mean(mae.hub.rep), mean(mae.loc.rep), mean(mae.ols.rep)), 3)
round(c(sd(mae.cen.rep), sd(mae.hub.rep), sd(mae.loc.rep), sd(mae.ols.rep)), 3)

round(apply(b.cen.rep, 2, mean), 3)
round(apply(b.cen.rep, 2, sd), 3)
round(apply(b.hub.rep, 2, mean), 3)
round(apply(b.hub.rep, 2, sd), 3)
round(apply(b.loc.rep, 2, mean), 3)
round(apply(b.loc.rep, 2, sd), 3)
round(apply(b.ols.rep, 2, mean), 3)
round(apply(b.ols.rep, 2, sd), 3)

