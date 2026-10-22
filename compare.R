
rm(list=ls())
gc()

setwd("D:\\GAO Junzhuo\\ADMM_Huber_decen_vert")

library(MASS)

thres <- function(a, t) {
  sign(a) * pmax(abs(a) - t, 0)
}

trunc <- function(a, tau) {
  pmin(pmax(a, -tau), tau)
}

loss <- function(x, tau) {
  (abs(x) <= tau) * x ^ 2 / 2 + (abs(x) > tau) * (tau * abs(x) - tau ^ 2 / 2)
}

n <- 200
d <- 500
mux <- rep(0, d)
sigmax <- matrix(NA, d, d)
rhox <- 0.5
for (i in 1:d) {
  for (j in 1:d) {
    sigmax[i, j] <- rhox ^ abs(i - j)
  }
}
s <- 5
bt <- c(rep(1, s), rep(0, d - s))
ga <- 0.5
p1 <- 0.1
p2 <- 0.1

m <- 10
A <- matrix(0, m, m)
set.seed(1)
index <- rbinom(m * (m - 1) / 2, size=1, prob=0.5)
t <- 1
for(i in 1:(m - 1)) {
  for(j in (i + 1):m) {
    A[i, j] <- A[j, i] <- index[t]
    t <- t + 1
  }
}
diag(A) <- rep(1, m)
dj <- d / m

case <- 6
sim.max <- 100
itr.max <- 1000
StatErr.cen <- StatErr.hub <- StatErr.med <- StatErr.lsq <- matrix(NA, sim.max, itr.max + 1)
StatErr.loc <- matrix(0, sim.max, itr.max + 1)
for (sim.idx in 1:sim.max) {
  cat(sim.idx)
  set.seed(sim.idx)
  x <- mvrnorm(n, mux, sigmax)
  if (case==1) {
    e <- rnorm(n, 0, 1)
    y <- x %*% bt + e
  } else if (case==2) {
    e <- rgamma(n, 2, 1)
    y <- x %*% bt + e
  } else if (case==3) {
    e <- (0.5 + 0.5 * abs(x[, 1])) * rnorm(n, 0, 2)
    y <- x %*% bt + e
  } else if (case==4) {
    e <- rt(n, df=3)
    y <- x %*% bt + e
  } else if (case==5) {
    e <- rt(n, 3)
    y <- x %*% bt + e
    y[1:(0.1 * n)] <- y[1:(0.1 * n)] + 5
  } else {
    e <- rt(n, 3)
    y <- x %*% bt + e
    y[1:(0.2 * n)] <- y[1:(0.2 * n)] + 5
  }
  
  xm <- array(NA, c(n, m, dj))
  for (j in 1:m) {
    xm[, j, ] <- x[, ((j - 1) * dj + 1):(j * dj)]
  }
  
  tau <- sqrt(n / log(d)) / 5
  
  # VFLHuber
  cc <- 10
  b <- array(0, c(itr.max + 1, m, dj))
  StatErr.hub[sim.idx, 1] <- sqrt(sum((as.vector(t(b[1, , ])) - bt) ^ 2))
  u <- v <- a <- array(0, c(itr.max + 1, m, n))
  for (itr in 1:itr.max) {
    for (j in 1:m) {
      b[itr + 1, j, ] <- thres(-t(xm[, j, ]) %*% u[itr, j, ] / (2 * n * p2), p1 / (2 * p2))
    }
    StatErr.hub[sim.idx, itr + 1] <- sqrt(sum((as.vector(t(b[itr + 1, , ])) - bt) ^ 2))
    for (j in 1:m) {
      su <- sum(A[j, ]) * u[itr, j, ] - apply(u[itr, A[j, ]!=0, ], 2, sum)
      v[itr + 1, j, ] <- u[itr, j, ] - (-xm[, j, ] %*% b[itr + 1, j, ] + ga * su + a[itr, j, ]) / cc
    }
    for (j in 1:m) {
      sv <- sum(A[j, ]) * v[itr + 1, j, ] - apply(v[itr + 1, A[j, ]!=0, ], 2, sum)
      a[itr + 1, j, ] <- a[itr, j, ] + ga * sv
      u[itr + 1, j, ] <- trunc((v[itr + 1, j, ] - y / (m * cc)) / (1 + 1 / (m * cc)), tau)
    }
  }
  
  # CenHuber
  itr.max <- 1000
  cc <- 10
  b.cen <- matrix(0, itr.max + 1, d)
  StatErr.cen[sim.idx, 1] <- sqrt(sum((b.cen[1, ] - bt) ^ 2))
  x.cen <- z.cen <- v.cen <- matrix(0, itr.max + 1, n)
  for (itr in 1:itr.max) {
    b.cen[itr + 1, ] <- thres(-t(x) %*% x.cen[itr, ] / (2 * n * p2), p1 / (2 * p2))
    StatErr.cen[sim.idx, itr + 1] <- sqrt(sum((b.cen[itr + 1, ] - bt) ^ 2))
    x.cen[itr + 1, ] <- (1 - ga / cc) * x.cen[itr, ] + x %*% b.cen[itr + 1, ] / cc -
      v.cen[itr, ] / cc + ga * z.cen[itr, ] / cc
    z.cen[itr + 1, ] <- trunc((x.cen[itr + 1, ] + v.cen[itr, ] / ga - y / ga) / (1 + 1 / ga), tau)
    v.cen[itr + 1, ] <- v.cen[itr, ] + ga * (x.cen[itr + 1, ] - z.cen[itr + 1, ])
  }
  
  # LocHuber
  itr.max <- 1000
  cc <- 10
  b.loc <- matrix(0, itr.max + 1, d)
  StatErr.loc[sim.idx, 1] <- sum((b.loc[1, ] - bt) ^ 2)
  for (j in 1:m) {
    x.loc <- z.loc <- v.loc <- matrix(0, itr.max + 1, n)
    for (itr in 1:itr.max) {
      b.loc[itr + 1, ((j - 1) * dj + 1):(j * dj)] <- thres(-t(xm[, j, ]) %*% x.loc[itr, ] / (2 * n * p2), p1 / (2 * p2))
      StatErr.loc[sim.idx, itr + 1] <- StatErr.loc[sim.idx, itr + 1] + 
        sum((b.loc[itr + 1, ((j - 1) * dj + 1):(j * dj)] - bt[((j - 1) * dj + 1):(j * dj)]) ^ 2)
      x.loc[itr + 1, ] <- (1 - ga / cc) * x.loc[itr, ] + xm[, j, ] %*% b.loc[itr + 1, ((j - 1) * dj + 1):(j * dj)] / cc -
        v.loc[itr, ] / cc + ga * z.loc[itr, ] / cc
      z.loc[itr + 1, ] <- trunc((x.loc[itr + 1, ] + v.loc[itr, ] / ga - y / ga) / (1 + 1 / ga), tau)
      v.loc[itr + 1, ] <- v.loc[itr, ] + ga * (x.loc[itr + 1, ] - z.loc[itr + 1, ])
    }
  }
  StatErr.loc[sim.idx, ] <- sqrt(StatErr.loc[sim.idx, ])
  
  # VFLMed
  itr.max <- 1000
  cc <- 10
  b <- array(0, c(itr.max + 1, m, dj))
  StatErr.med[sim.idx, 1] <- sqrt(sum((as.vector(t(b[1, , ])) - bt) ^ 2))
  u <- v <- a <- array(0, c(itr.max + 1, m, n))
  for (itr in 1:itr.max) {
    for (j in 1:m) {
      b[itr + 1, j, ] <- thres(-t(xm[, j, ]) %*% u[itr, j, ] / (2 * n * p2), p1 / (2 * p2))
    }
    StatErr.med[sim.idx, itr + 1] <- sqrt(sum((as.vector(t(b[itr + 1, , ])) - bt) ^ 2))
    for (j in 1:m) {
      su <- sum(A[j, ]) * u[itr, j, ] - apply(u[itr, A[j, ]!=0, ], 2, sum)
      v[itr + 1, j, ] <- u[itr, j, ] - (-xm[, j, ] %*% b[itr + 1, j, ] + ga * su + a[itr, j, ]) / cc
    }
    for (j in 1:m) {
      sv <- sum(A[j, ]) * v[itr + 1, j, ] - apply(v[itr + 1, A[j, ]!=0, ], 2, sum)
      a[itr + 1, j, ] <- a[itr, j, ] + ga * sv
      u[itr + 1, j, ] <- trunc(v[itr + 1, j, ] - y / (m * cc), 1)
    }
  }
  
  # VFLLS
  itr.max <- 1000
  cc <- 10
  b <- array(0, c(itr.max + 1, m, dj))
  StatErr.lsq[sim.idx, 1] <- sqrt(sum((as.vector(t(b[1, , ])) - bt) ^ 2))
  u <- v <- a <- array(0, c(itr.max + 1, m, n))
  for (itr in 1:itr.max) {
    for (j in 1:m) {
      b[itr + 1, j, ] <- thres(-t(xm[, j, ]) %*% u[itr, j, ] / (2 * n * p2), p1 / (2 * p2))
    }
    StatErr.lsq[sim.idx, itr + 1] <- sqrt(sum((as.vector(t(b[itr + 1, , ])) - bt) ^ 2))
    for (j in 1:m) {
      su <- sum(A[j, ]) * u[itr, j, ] - apply(u[itr, A[j, ]!=0, ], 2, sum)
      v[itr + 1, j, ] <- u[itr, j, ] - (-xm[, j, ] %*% b[itr + 1, j, ] + ga * su + a[itr, j, ]) / cc
    }
    for (j in 1:m) {
      sv <- sum(A[j, ]) * v[itr + 1, j, ] - apply(v[itr + 1, A[j, ]!=0, ], 2, sum)
      a[itr + 1, j, ] <- a[itr, j, ] + ga * sv
      u[itr + 1, j, ] <- (v[itr + 1, j, ] - y / (m * cc)) / (1 + 1 / (m * cc))
    }
  }
}

plot(apply(StatErr.cen, 2, mean))
plot(apply(StatErr.loc, 2, mean))
plot(apply(StatErr.hub, 2, mean))
plot(apply(StatErr.med, 2, mean))
plot(apply(StatErr.lsq, 2, mean))

# save(StatErr.cen, StatErr.loc, StatErr.hub, StatErr.med, StatErr.lsq, sim.max, itr.max, file="compare6.RData")

library(ggplot2)
library(ggpubr)
library(latex2exp)

# load("compare6.RData")

df <- data.frame(group=rep(c("CenHuber", "LocHuber", "VFLHuber", "VFLMed", "VFLLS"), each=sim.max),
                 value=c(StatErr.cen[, itr.max + 1], StatErr.loc[, itr.max + 1], StatErr.hub[, itr.max + 1], 
                         StatErr.med[, itr.max + 1], StatErr.lsq[, itr.max + 1]))

p <- ggplot(df, aes(x=group, y=value, fill=group)) +
  geom_boxplot() +
  labs(x="", y="StatErr", title=TeX("$\\epsilon\\sim t(3)$, 20% contamination")) +
  scale_fill_manual(values=c("CenHuber"="skyblue", "LocHuber"="lightcoral", "VFLHuber"="orange", 
                             "VFLMed"="forestgreen", "VFLLS"="purple")) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=40),
        legend.position="none",
        axis.text.y=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=35, hjust=0.5))
p
# \\epsilon\\sim N(0,1)
# \\epsilon\\sim Gamma(2,1)
# $\\epsilon\\sim (0.5+0.5|x_1|)N(0,4)$
# $\\epsilon\\sim t(3)$
# $\\epsilon\\sim t(3)$, 10% contamination
# $\\epsilon\\sim t(3)$, 20% contamination

# ggsave("compare6.eps", plot=p, device="eps", width=8, height=8)

