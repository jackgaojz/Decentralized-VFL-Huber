
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

set.seed(1)
x <- mvrnorm(n, mux, sigmax)
e <- rnorm(n, 0, 1)
y <- x %*% bt + e

cc <- 10
p1 <- 0.1
p2 <- 0.1
tau <- 0.5

# CenHuber
ga <- 0.5
itr.max <- 1000
StatErr.cen <- Obj.cen <- ReObj.cen <- rep(NA, itr.max + 1)
b.cen <- matrix(0, itr.max + 1, d)
x.cen <- z.cen <- v.cen <- matrix(0, itr.max + 1, n)
StatErr.cen[1] <- sqrt(sum((b.cen[1, ] - bt) ^ 2))
Obj.cen[1] <- mean(loss(y - x %*% b.cen[1, ], tau)) +
  p2 * sum(b.cen[1, ] ^ 2) + p1 * sum(abs(b.cen[1, ]))
for (itr in 1:itr.max) {
  b.cen[itr + 1, ] <- thres(-t(x) %*% x.cen[itr, ] / (2 * n * p2), p1 / (2 * p2))
  StatErr.cen[itr + 1] <- sqrt(sum((b.cen[itr + 1, ] - bt) ^ 2))
  Obj.cen[itr + 1] <- mean(loss(y - x %*% b.cen[itr + 1, ], tau)) +
    p2 * sum(b.cen[itr + 1, ] ^ 2) + p1 * sum(abs(b.cen[itr + 1, ]))
  x.cen[itr + 1, ] <- (1 - ga / cc) * x.cen[itr, ] + x %*% b.cen[itr + 1, ] / cc -
    v.cen[itr, ] / cc + ga * z.cen[itr, ] / cc
  z.cen[itr + 1, ] <- trunc((x.cen[itr + 1, ] + v.cen[itr, ] / ga - y / ga) / (1 + 1 / ga), tau)
  v.cen[itr + 1, ] <- v.cen[itr, ] + ga * (x.cen[itr + 1, ] - z.cen[itr + 1, ])
}
plot(log(Obj.cen))
plot(log(StatErr.cen))
b.cen[itr.max + 1, ]

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
xm <- array(NA, c(n, m, dj))
for (j in 1:m) {
  xm[, j, ] <- x[, ((j - 1) * dj + 1):(j * dj)]
}

# VFLHuber
ga.set <- c(0.05, 0.1, 0.5, 2)
itr.max <- 800
StatErr <- OptErr <- Obj <- ReObj <- matrix(NA, itr.max + 1, length(ga.set))
b <- array(0, c(itr.max + 1, m, dj))
u <- v <- a <- array(0, c(itr.max + 1, m, n))
StatErr[1, ] <- sqrt(sum((as.vector(t(b[1, , ])) - bt) ^ 2))
OptErr[1, ] <- sqrt(sum((as.vector(t(b[1, , ])) - b.cen[itr.max + 1, ]) ^ 2))
Obj[1, ] <- mean(loss(y - x %*% as.vector(t(b[1, , ])), tau)) + 
  p2 * sum(b[1, , ] ^ 2) + p1 * sum(abs(b[1, , ]))
ReObj[1, ] <- mean(loss(y - x %*% as.vector(t(b[1, , ])), tau)) +
  p2 * sum(b[1, , ] ^ 2) + p1 * sum(abs(b[1, , ])) -
  mean(loss(y - x %*% b.cen[itr.max + 1, ], tau)) - 
  p2 * sum(b.cen[itr.max + 1, ] ^ 2) - 
  p1 * sum(abs(b.cen[itr.max + 1, ]))
for (ga.itr in 1:length(ga.set)) {
  cat(ga.itr)
  ga <- ga.set[ga.itr]
  for (itr in 1:itr.max) {
    if (ga.itr==4 & itr>=15) {
      break
    } else {
      for (j in 1:m) {
        b[itr + 1, j, ] <- thres(-t(xm[, j, ]) %*% u[itr, j, ] / (2 * n * p2), p1 / (2 * p2))
      }
      StatErr[itr + 1, ga.itr] <- sqrt(sum((as.vector(t(b[itr + 1, , ])) - bt) ^ 2))
      OptErr[itr + 1, ga.itr] <- sqrt(sum((as.vector(t(b[itr + 1, , ])) - b.cen[itr.max + 1, ]) ^ 2))
      Obj[itr + 1, ga.itr] <- mean(loss(y - x %*% as.vector(t(b[itr + 1, , ])), tau)) + 
        p2 * sum(b[itr + 1, , ] ^ 2) + p1 * sum(abs(b[itr + 1, , ]))
      ReObj[itr + 1, ga.itr] <- mean(loss(y - x %*% as.vector(t(b[itr + 1, , ])), tau)) +
        p2 * sum(b[itr + 1, , ] ^ 2) + p1 * sum(abs(b[itr + 1, , ])) -
        mean(loss(y - x %*% b.cen[itr.max + 1, ], tau)) - p2 * sum(b.cen[itr.max + 1, ] ^ 2) - p1 * sum(abs(b.cen[itr.max + 1, ]))
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
  }
}

library(ggplot2)
library(ggpubr)

tab1 <- data.frame(index=rep(0:itr.max, times=length(ga.set)),
                   StatErr=log(c(as.vector(StatErr[, 1]), 
                                 as.vector(StatErr[, 2]), 
                                 as.vector(StatErr[, 3]),
                                 as.vector(StatErr[, 4]))),
                   method=rep(c("0.05", "0.1", "0.5", "2"), each=itr.max + 1))

pic1 <- ggplot(tab1, aes(x=index, y=StatErr, color=method, linetype=method, shape=method)) +
  geom_line(linewidth=3, data=subset(tab1, (!(method=="2" & index > 14)))) +
  geom_point(size=5, data=subset(tab1, (!(method=="2" & index > 14)) & index %% 50==0)) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5)) +
  xlab(label=expression(t)) +
  ylab(label=expression(log(StatErr))) +
  labs(color=expression(gamma), shape=expression(gamma), linetype=expression(gamma)) +
  theme(panel.border=element_rect(color="black", fill=NA),
        axis.title=element_text(size=25),
        axis.text.x=element_text(size=25),
        axis.text.y=element_text(size=25),
        legend.text=element_text(size=25),
        legend.title=element_text(size=25),
        plot.title=element_text(size=25))

tab2 <- data.frame(index=rep(0:itr.max, times=length(ga.set)),
                   OptErr=log(c(as.vector(OptErr[, 1]), 
                                as.vector(OptErr[, 2]), 
                                as.vector(OptErr[, 3]),
                                as.vector(OptErr[, 4]))),
                   method=rep(c("0.05", "0.1", "0.5", "2"), each=itr.max + 1))

pic2 <- ggplot(tab2, aes(x=index, y=OptErr, color=method, linetype=method, shape=method)) +
  geom_line(linewidth=3, data=subset(tab2, (!(method=="2" & index > 14)))) +
  geom_point(size=5, data=subset(tab2, (!(method=="2" & index > 14)) & index %% 50==0)) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5)) +
  xlab(label=expression(t)) +
  ylab(label=expression(log(OptErr))) +
  labs(color="c", shape="c", linetype="c") +
  theme(panel.border=element_rect(color="black", fill=NA),
        axis.title=element_text(size=25),
        axis.text.x=element_text(size=25),
        axis.text.y=element_text(size=25),
        legend.text=element_text(size=25),
        legend.title=element_text(size=25),
        plot.title=element_text(size=25))

tab3 <- data.frame(index=rep(0:itr.max, times=length(ga.set)),
                   Obj=log(c(as.vector(Obj[, 1]), 
                             as.vector(Obj[, 2]), 
                             as.vector(Obj[, 3]),
                             as.vector(Obj[, 4]))),
                   method=rep(c("0.05", "0.1", "0.5", "2"), each=itr.max + 1))

pic3 <- ggplot(tab3, aes(x=index, y=Obj, color=method, linetype=method, shape=method)) +
  geom_line(linewidth=3, data=subset(tab3, (!(method=="2" & index > 14)))) +
  geom_point(size=5, data=subset(tab3, (!(method=="2" & index > 14)) & index %% 50==0)) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5)) +
  xlab(label=expression(t)) +
  ylab(label=expression(log(Obj))) +
  labs(color="c", shape="c", linetype="c") +
  theme(panel.border=element_rect(color="black", fill=NA),
        axis.title=element_text(size=25),
        axis.text.x=element_text(size=25),
        axis.text.y=element_text(size=25),
        legend.text=element_text(size=25),
        legend.title=element_text(size=25),
        plot.title=element_text(size=25))

tab4 <- data.frame(index=rep(0:itr.max, times=length(ga.set)),
                   ReObj=log(c(as.vector(ReObj[, 1]), 
                               as.vector(ReObj[, 2]), 
                               as.vector(ReObj[, 3]),
                               as.vector(ReObj[, 4]))),
                   method=rep(c("0.05", "0.1", "0.5", "2"), each=itr.max + 1))

pic4 <- ggplot(tab4, aes(x=index, y=ReObj, color=method, linetype=method, shape=method)) +
  geom_line(linewidth=3, data=subset(tab4, (!(method=="2" & index > 14)))) +
  geom_point(size=5, data=subset(tab4, (!(method=="2" & index > 14)) & index %% 50==0)) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5)) +
  xlab(label=expression(t)) +
  ylab(label=expression(log(ReObj))) +
  labs(color="c", shape="c", linetype="c") +
  theme(panel.border=element_rect(color="black", fill=NA),
        axis.title=element_text(size=25),
        axis.text.x=element_text(size=25),
        axis.text.y=element_text(size=25),
        legend.text=element_text(size=25),
        legend.title=element_text(size=25),
        plot.title=element_text(size=25))

p <- ggarrange(pic1, pic2, pic4, nrow=1, ncol=3, common.legend=T)
p

# ggsave("gamma.eps", plot=p, device="eps", width=18, height=6)

