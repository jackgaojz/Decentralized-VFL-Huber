
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

n.set <- c(200, 350, 500)
n0 <- n.set[length(n.set)]
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
x0 <- mvrnorm(n0, mux, sigmax)
e0 <- rnorm(n0, 0, 1)

ga <- 0.5
p1 <- 0.1
p2 <- 0.1
tau <- 0.5

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
# ring_graph <- function(m, deg) {
#   W <- 0 * diag(m)
#   for (i in 1:deg) {
#     diag(W[-(1:i), -((m - i + 1):m)]) <- 1  
#     W[i, (m - deg + i):m] <- 1  
#   }
#   return(W + t(W) + diag(m))
# }
# A <- ring_graph(m, deg=1)

itr.max <- 500
itr.max.cen <- 1000
cc <- 25
StatErr.cen <- Obj.cen <- ReObj.cen <- matrix(NA, itr.max.cen + 1, length(n.set))
StatErr <- OptErr <- Obj <- ReObj <- matrix(NA, itr.max + 1, length(n.set))
for (n.itr in 1:length(n.set)) {
  cat(n.itr)
  n <- n.set[n.itr]
  x <- x0[1:n, ]
  e <- e0[1:n]
  y <- x %*% bt + e
  
  dj <- d / m
  xm <- array(NA, c(n, m, dj))
  for (j in 1:m) {
    xm[, j, ] <- x[, ((j - 1) * dj + 1):(j * dj)]
  }
  
  # CenHuber
  b.cen <- matrix(0, itr.max.cen + 1, d)
  x.cen <- z.cen <- v.cen <- matrix(0, itr.max.cen + 1, n)
  StatErr.cen[1, n.itr] <- sqrt(sum((b.cen[1, ] - bt) ^ 2))
  Obj.cen[1, n.itr] <- mean(loss(y - x %*% b.cen[1, ], tau)) +
    p2 * sum(b.cen[1, ] ^ 2) + p1 * sum(abs(b.cen[1, ]))
  for (itr in 1:itr.max.cen) {
    b.cen[itr + 1, ] <- thres(-t(x) %*% x.cen[itr, ] / (2 * n * p2), p1 / (2 * p2))
    StatErr.cen[itr + 1, n.itr] <- sqrt(sum((b.cen[itr + 1, ] - bt) ^ 2))
    Obj.cen[itr + 1, n.itr] <- mean(loss(y - x %*% b.cen[itr + 1, ], tau)) +
      p2 * sum(b.cen[itr + 1, ] ^ 2) + p1 * sum(abs(b.cen[itr + 1, ]))
    x.cen[itr + 1, ] <- (1 - ga / cc) * x.cen[itr, ] + x %*% b.cen[itr + 1, ] / cc -
      v.cen[itr, ] / cc + ga * z.cen[itr, ] / cc
    z.cen[itr + 1, ] <- trunc((x.cen[itr + 1, ] + v.cen[itr, ] / ga - y / ga) / (1 + 1 / ga), tau)
    v.cen[itr + 1, ] <- v.cen[itr, ] + ga * (x.cen[itr + 1, ] - z.cen[itr + 1, ])
  }
  # plot(log(Obj.cen[, n.itr]))
  # plot(log(StatErr.cen[, n.itr]))
  # b.cen[itr.max.cen + 1, ]
  
  b <- array(0, c(itr.max + 1, m, dj))
  u <- v <- a <- array(0, c(itr.max + 1, m, n))
  StatErr[1, n.itr] <- sqrt(sum((as.vector(t(b[1, , ])) - bt) ^ 2))
  OptErr[1, n.itr] <- sqrt(sum((as.vector(t(b[1, , ])) - b.cen[itr.max.cen + 1, ]) ^ 2))
  Obj[1, n.itr] <- mean(loss(y - x %*% as.vector(t(b[1, , ])), tau)) + 
    p2 * sum(b[1, , ] ^ 2) + p1 * sum(abs(b[1, , ]))
  ReObj[1, n.itr] <- mean(loss(y - x %*% as.vector(t(b[1, , ])), tau)) +
    p2 * sum(b[1, , ] ^ 2) + p1 * sum(abs(b[1, , ])) -
    mean(loss(y - x %*% b.cen[itr.max.cen + 1, ], tau)) - 
    p2 * sum(b.cen[itr.max.cen + 1, ] ^ 2) - 
    p1 * sum(abs(b.cen[itr.max.cen + 1, ]))
  for (itr in 1:itr.max) {
    for (j in 1:m) {
      b[itr + 1, j, ] <- thres(-t(xm[, j, ]) %*% u[itr, j, ] / (2 * n * p2), p1 / (2 * p2))
    }
    StatErr[itr + 1, n.itr] <- sqrt(sum((as.vector(t(b[itr + 1, , ])) - bt) ^ 2))
    OptErr[itr + 1, n.itr] <- sqrt(sum((as.vector(t(b[itr + 1, , ])) - b.cen[itr.max.cen + 1, ]) ^ 2))
    Obj[itr + 1, n.itr] <- mean(loss(y - x %*% as.vector(t(b[itr + 1, , ])), tau)) + 
      p2 * sum(b[itr + 1, , ] ^ 2) + p1 * sum(abs(b[itr + 1, , ]))
    ReObj[itr + 1, n.itr] <- mean(loss(y - x %*% as.vector(t(b[itr + 1, , ])), tau)) +
      p2 * sum(b[itr + 1, , ] ^ 2) + p1 * sum(abs(b[itr + 1, , ])) -
      mean(loss(y - x %*% b.cen[itr.max.cen + 1, ], tau)) - p2 * sum(b.cen[itr.max.cen + 1, ] ^ 2) - p1 * sum(abs(b.cen[itr.max.cen + 1, ]))
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

library(ggplot2)
library(ggpubr)
library(latex2exp)

tab1 <- data.frame(index=rep(0:itr.max, times=length(n.set)),
                   StatErr=log(c(as.vector(StatErr))),
                   method=rep(c("200", "350", "500"), each=itr.max + 1))

pic1 <- ggplot(tab1, aes(x=index, y=StatErr, color=method, linetype=method, shape=method)) +
  geom_line(linewidth=3) +
  geom_point(size=5, data=subset(tab1, index %% 50==0)) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5)) +
  xlab(label=expression(t)) +
  ylab(label=expression(log(StatErr))) +
  labs(color=TeX("$n$"), shape=TeX("$n$"), linetype=TeX("$n$")) +
  theme(panel.border=element_rect(color="black", fill=NA),
        axis.title=element_text(size=25),
        axis.text.x=element_text(size=25),
        axis.text.y=element_text(size=25),
        legend.text=element_text(size=25),
        legend.title=element_text(size=25),
        plot.title=element_text(size=25))

tab2 <- data.frame(index=rep(0:itr.max, times=length(n.set)),
                   OptErr=log(c(as.vector(OptErr))),
                   method=rep(c("200", "350", "500"), each=itr.max + 1))

pic2 <- ggplot(tab2, aes(x=index, y=OptErr, color=method, linetype=method, shape=method)) +
  geom_line(linewidth=3) +
  geom_point(size=5, data=subset(tab2, index %% 50==0)) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5)) +
  xlab(label=expression(t)) +
  ylab(label=expression(log(OptErr))) +
  labs(color=TeX("n"), shape=TeX("n"), linetype=TeX("n")) +
  theme(panel.border=element_rect(color="black", fill=NA),
        axis.title=element_text(size=25),
        axis.text.x=element_text(size=25),
        axis.text.y=element_text(size=25),
        legend.text=element_text(size=25),
        legend.title=element_text(size=25),
        plot.title=element_text(size=25))

p <- ggarrange(pic1, pic2, nrow=2, ncol=1, common.legend=T)
p

# ggsave("n.eps", plot=p, device="eps", width=6, height=12)
