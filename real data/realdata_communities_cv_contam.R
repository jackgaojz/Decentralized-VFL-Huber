rm(list = ls())
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
  (abs(x) <= tau) * x^2 / 2 + (abs(x) > tau) * (tau * abs(x) - tau^2 / 2)
}

huber_train_loss <- function(x, y, beta, tau) {
  mean(loss(as.vector(y - x %*% beta), tau))
}

read_communities_data <- function() {
  name_lines <- readLines("communities.names", warn = FALSE)
  attr_lines <- grep("^@attribute ", name_lines, value = TRUE)
  var_names <- sub("^@attribute ([^ ]+) .*$", "\\1", attr_lines)

  dat <- read.csv("communities.data", header = FALSE, na.strings = "?",
                  stringsAsFactors = FALSE)
  colnames(dat) <- var_names

  drop_cols <- c("state", "county", "community", "communityname", "fold")
  y <- as.numeric(dat$ViolentCrimesPerPop)
  xdat <- dat[, setdiff(names(dat), c(drop_cols, "ViolentCrimesPerPop"))]
  xdat <- as.data.frame(lapply(xdat, function(z) as.numeric(z)))

  # Median imputation keeps the law-enforcement block in the VFL experiment.
  for (nm in names(xdat)) {
    med <- median(xdat[[nm]], na.rm = TRUE)
    xdat[[nm]][is.na(xdat[[nm]])] <- med
  }

  list(xdat = xdat, y = y)
}

make_semantic_partition <- function(xdat) {
  vars <- names(xdat)

  demographic_vars <- c(
    "population", "householdsize", "racepctblack", "racePctWhite",
    "racePctAsian", "racePctHisp", "agePct12t21", "agePct12t29",
    "agePct16t24", "agePct65up", "numbUrban", "pctUrban",
    "NumImmig", "PctImmigRecent", "PctImmigRec5", "PctImmigRec8",
    "PctImmigRec10", "PctRecentImmig", "PctRecImmig5",
    "PctRecImmig8", "PctRecImmig10", "PctSpeakEnglOnly",
    "PctNotSpeakEnglWell", "PctForeignBorn"
  )

  socioeconomic_vars <- c(
    "medIncome", "pctWWage", "pctWFarmSelf", "pctWInvInc",
    "pctWSocSec", "pctWPubAsst", "pctWRetire", "medFamInc",
    "perCapInc", "whitePerCap", "blackPerCap", "indianPerCap",
    "AsianPerCap", "OtherPerCap", "HispPerCap", "NumUnderPov",
    "PctPopUnderPov", "PctLess9thGrade", "PctNotHSGrad",
    "PctBSorMore", "PctUnemployed", "PctEmploy", "PctEmplManu",
    "PctEmplProfServ", "PctOccupManu", "PctOccupMgmtProf"
  )

  law_vars <- vars[grepl("^(Lemas|Polic|PctPolic|RacialMatchCommPol|OfficAssgnDrugUnits|NumKindsDrugsSeiz)", vars)]

  group1 <- c("Intercept", intersect(demographic_vars, vars))
  group2 <- intersect(socioeconomic_vars, vars)
  group4 <- intersect(law_vars, vars)
  assigned <- c(setdiff(group1, "Intercept"), group2, group4)
  group3 <- setdiff(vars, assigned)

  list(
    CensusDemographic = group1,
    Socioeconomic = group2,
    HousingFamilyMobility = group3,
    LawEnforcement = group4
  )
}

build_design <- function(xdat, groups) {
  xfull <- cbind(Intercept = 1, as.matrix(xdat))
  ordered_vars <- unlist(groups, use.names = FALSE)
  x <- xfull[, ordered_vars, drop = FALSE]
  dj <- as.numeric(sapply(groups, length))
  djk <- c(0, cumsum(dj))
  list(x = x, dj = dj, djk = djk, ordered_vars = ordered_vars)
}

make_xm <- function(x, dj) {
  n <- nrow(x)
  m <- length(dj)
  xm <- array(NA, c(n, m, max(dj)))
  start <- 1
  for (j in 1:m) {
    idx <- start:(start + dj[j] - 1)
    xm[, j, 1:dj[j]] <- x[, idx, drop = FALSE]
    start <- start + dj[j]
  }
  xm
}

fit_vfl <- function(x, y, dj, djk, A, p1, p2, tau = 1, ga = 0.5,
                    cc = 20, itr.max = 500, type = "huber",
                    return_loss = FALSE, early_stop_loss = Inf) {
  n <- nrow(x)
  d <- ncol(x)
  m <- length(dj)
  xm <- make_xm(x, dj)
  b <- matrix(0, itr.max + 1, d)
  u <- v <- a <- array(0, c(itr.max + 1, m, n))
  train_huber_loss <- rep(NA, itr.max + 1)
  train_huber_loss[1] <- huber_train_loss(x, y, b[1, ], tau)

  for (itr in 1:itr.max) {
    for (j in 1:m) {
      idx <- (djk[j] + 1):djk[j + 1]
      b[itr + 1, idx] <- thres(-t(xm[, j, 1:dj[j]]) %*% u[itr, j, ] /
                                (2 * n * p2), p1 / (2 * p2))
    }
    train_huber_loss[itr + 1] <- huber_train_loss(x, y, b[itr + 1, ], tau)
    if (!is.finite(train_huber_loss[itr + 1]) ||
        train_huber_loss[itr + 1] > early_stop_loss) {
      attr(b, "early_stop") <- TRUE
      if (return_loss) {
        return(list(beta = b, train_huber_loss = train_huber_loss,
                    early_stop = TRUE))
      }
      return(b)
    }

    for (j in 1:m) {
      idx <- (djk[j] + 1):djk[j + 1]
      su <- sum(A[j, ]) * u[itr, j, ] - apply(u[itr, A[j, ] != 0, ], 2, sum)
      v[itr + 1, j, ] <- u[itr, j, ] -
        (-xm[, j, 1:dj[j]] %*% b[itr + 1, idx] + ga * su + a[itr, j, ]) / cc
    }
    for (j in 1:m) {
      sv <- sum(A[j, ]) * v[itr + 1, j, ] - apply(v[itr + 1, A[j, ] != 0, ], 2, sum)
      a[itr + 1, j, ] <- a[itr, j, ] + ga * sv
      if (type == "huber") {
        u[itr + 1, j, ] <- trunc((v[itr + 1, j, ] - y / (m * cc)) /
                                   (1 + 1 / (m * cc)), tau)
      } else if (type == "median") {
        u[itr + 1, j, ] <- trunc(v[itr + 1, j, ] - y / (m * cc), 1)
      } else if (type == "ls") {
        u[itr + 1, j, ] <- (v[itr + 1, j, ] - y / (m * cc)) /
          (1 + 1 / (m * cc))
      }
    }
  }

  attr(b, "early_stop") <- FALSE
  if (return_loss) {
    return(list(beta = b, train_huber_loss = train_huber_loss,
                early_stop = FALSE))
  }
  b
}

fit_central_huber <- function(x, y, p1, p2, tau = 1, ga = 0.5,
                              cc = 20, itr.max = 500,
                              return_loss = FALSE, early_stop_loss = Inf) {
  n <- nrow(x)
  d <- ncol(x)
  b <- matrix(0, itr.max + 1, d)
  u <- z <- v <- matrix(0, itr.max + 1, n)
  train_huber_loss <- rep(NA, itr.max + 1)
  train_huber_loss[1] <- huber_train_loss(x, y, b[1, ], tau)

  for (itr in 1:itr.max) {
    b[itr + 1, ] <- thres(-t(x) %*% u[itr, ] / (2 * n * p2), p1 / (2 * p2))
    train_huber_loss[itr + 1] <- huber_train_loss(x, y, b[itr + 1, ], tau)
    if (!is.finite(train_huber_loss[itr + 1]) ||
        train_huber_loss[itr + 1] > early_stop_loss) {
      attr(b, "early_stop") <- TRUE
      if (return_loss) {
        return(list(beta = b, train_huber_loss = train_huber_loss,
                    early_stop = TRUE))
      }
      return(b)
    }
    u[itr + 1, ] <- (1 - ga / cc) * u[itr, ] + x %*% b[itr + 1, ] / cc -
      v[itr, ] / cc + ga * z[itr, ] / cc
    z[itr + 1, ] <- trunc((u[itr + 1, ] + v[itr, ] / ga - y / ga) /
                            (1 + 1 / ga), tau)
    v[itr + 1, ] <- v[itr, ] + ga * (u[itr + 1, ] - z[itr + 1, ])
  }

  attr(b, "early_stop") <- FALSE
  if (return_loss) {
    return(list(beta = b, train_huber_loss = train_huber_loss,
                early_stop = FALSE))
  }
  b
}

fit_local_huber <- function(x, y, dj, djk, p1, p2, tau = 1, ga = 0.5,
                            cc = 20, itr.max = 500,
                            return_loss = FALSE) {
  n <- nrow(x)
  d <- ncol(x)
  m <- length(dj)
  xm <- make_xm(x, dj)
  b <- matrix(0, itr.max + 1, d)
  train_huber_loss <- rep(NA, itr.max + 1)
  train_huber_loss[1] <- huber_train_loss(x, y, b[1, ], tau)

  for (j in 1:m) {
    x.loc <- z.loc <- v.loc <- matrix(0, itr.max + 1, n)
    idx <- (djk[j] + 1):djk[j + 1]
    for (itr in 1:itr.max) {
      b[itr + 1, idx] <- thres(-t(xm[, j, 1:dj[j]]) %*% x.loc[itr, ] /
                                (2 * n * p2), p1 / (2 * p2))
      x.loc[itr + 1, ] <- (1 - ga / cc) * x.loc[itr, ] +
        xm[, j, 1:dj[j]] %*% b[itr + 1, idx] / cc - v.loc[itr, ] / cc +
        ga * z.loc[itr, ] / cc
      z.loc[itr + 1, ] <- trunc((x.loc[itr + 1, ] + v.loc[itr, ] / ga - y / ga) /
                                  (1 + 1 / ga), tau)
      v.loc[itr + 1, ] <- v.loc[itr, ] + ga * (x.loc[itr + 1, ] - z.loc[itr + 1, ])
      train_huber_loss[itr + 1] <- huber_train_loss(x, y, b[itr + 1, ], tau)
    }
  }

  if (return_loss) {
    return(list(beta = b, train_huber_loss = train_huber_loss))
  }
  b
}

fit_vfl_nonhuber <- function(x, y, dj, djk, A, p1, p2, ga = 0.5,
                             cc = 20, itr.max = 500,
                             type = c("median", "ls"),
                             return_loss = FALSE) {
  type <- match.arg(type)
  n <- nrow(x)
  d <- ncol(x)
  m <- length(dj)
  xm <- make_xm(x, dj)
  b <- matrix(0, itr.max + 1, d)
  u <- v <- a <- array(0, c(itr.max + 1, m, n))
  train_loss <- rep(NA, itr.max + 1)
  res0 <- as.vector(y - x %*% b[1, ])
  train_loss[1] <- if (type == "median") mean(abs(res0)) else mean(res0^2)

  for (itr in 1:itr.max) {
    for (j in 1:m) {
      idx <- (djk[j] + 1):djk[j + 1]
      b[itr + 1, idx] <- thres(-t(xm[, j, 1:dj[j]]) %*% u[itr, j, ] /
                                (2 * n * p2), p1 / (2 * p2))
    }
    res <- as.vector(y - x %*% b[itr + 1, ])
    train_loss[itr + 1] <- if (type == "median") mean(abs(res)) else mean(res^2)

    for (j in 1:m) {
      idx <- (djk[j] + 1):djk[j + 1]
      su <- sum(A[j, ]) * u[itr, j, ] - apply(u[itr, A[j, ] != 0, ], 2, sum)
      v[itr + 1, j, ] <- u[itr, j, ] -
        (-xm[, j, 1:dj[j]] %*% b[itr + 1, idx] + ga * su + a[itr, j, ]) / cc
    }
    for (j in 1:m) {
      sv <- sum(A[j, ]) * v[itr + 1, j, ] - apply(v[itr + 1, A[j, ] != 0, ], 2, sum)
      a[itr + 1, j, ] <- a[itr, j, ] + ga * sv
      if (type == "median") {
        u[itr + 1, j, ] <- trunc(v[itr + 1, j, ] - y / (m * cc), 1)
      } else {
        u[itr + 1, j, ] <- (v[itr + 1, j, ] - y / (m * cc)) /
          (1 + 1 / (m * cc))
      }
    }
  }

  if (return_loss) {
    return(list(beta = b, train_loss = train_loss))
  }
  b
}
cv_vfl_huber <- function(x, y, dj, djk, A, p1.set, p2.set, tau.set, ga, cc,
                         itr.max, k.fold = 5, early_stop_loss = Inf) {
  fold <- createFolds(1:nrow(x), k = k.fold, list = TRUE)
  err <- array(NA, c(k.fold, length(p1.set), length(p2.set), length(tau.set)))

  for (k in 1:k.fold) {
    cat("CV fold", k, "for huber\n")
    val_id <- as.vector(unlist(fold[k]))
    tr_id <- setdiff(seq_len(nrow(x)), val_id)
    for (p1.idx in seq_along(p1.set)) {
      for (p2.idx in seq_along(p2.set)) {
        for (tau.idx in seq_along(tau.set)) {
          tau.use <- tau.set[tau.idx]
          b <- fit_vfl(x[tr_id, , drop = FALSE], y[tr_id], dj, djk, A,
                       p1 = p1.set[p1.idx], p2 = p2.set[p2.idx], tau = tau.use,
                       ga = ga, cc = cc, itr.max = itr.max, type = "huber",
                       early_stop_loss = early_stop_loss)
          if (isTRUE(attr(b, "early_stop"))) {
            err[k, p1.idx, p2.idx, tau.idx] <- Inf
            next
          }
          pred <- x[val_id, , drop = FALSE] %*% b[itr.max + 1, ]
          res <- y[val_id] - pred
          err[k, p1.idx, p2.idx, tau.idx] <- mean(res^2)
        }
      }
    }
  }

  err.mean <- apply(err, c(2, 3, 4), mean, na.rm = TRUE)
  err.mean[!is.finite(err.mean)] <- Inf
  if (all(is.infinite(err.mean))) {
    stop("All VFL Huber cross-validation errors are non-finite. Please increase cc or use a larger p2.")
  }
  best <- which(err.mean == min(err.mean), arr.ind = TRUE)[1, ]
  list(p1 = p1.set[best[1]], p2 = p2.set[best[2]], tau = tau.set[best[3]],
       err.mean = err.mean)
}

cv_vfl_nonhuber <- function(x, y, dj, djk, A, p1.set, p2.set, ga, cc,
                            itr.max, type = c("median", "ls"), k.fold = 5) {
  type <- match.arg(type)
  fold <- createFolds(1:nrow(x), k = k.fold, list = TRUE)
  err <- array(NA, c(k.fold, length(p1.set), length(p2.set)))

  for (k in 1:k.fold) {
    cat("CV fold", k, "for", type, "\n")
    val_id <- as.vector(unlist(fold[k]))
    tr_id <- setdiff(seq_len(nrow(x)), val_id)
    for (p1.idx in seq_along(p1.set)) {
      for (p2.idx in seq_along(p2.set)) {
        b <- fit_vfl_nonhuber(x[tr_id, , drop = FALSE], y[tr_id], dj, djk, A,
                            p1 = p1.set[p1.idx], p2 = p2.set[p2.idx],
                            ga = ga, cc = cc, itr.max = itr.max, type = type)
        pred <- x[val_id, , drop = FALSE] %*% b[itr.max + 1, ]
        res <- y[val_id] - pred
        if (type == "median") err[k, p1.idx, p2.idx] <- mean(abs(res))
        if (type == "ls") err[k, p1.idx, p2.idx] <- mean(res^2)
      }
    }
  }

  err.mean <- apply(err, c(2, 3), mean, na.rm = TRUE)
  err.mean[!is.finite(err.mean)] <- Inf
  if (all(is.infinite(err.mean))) {
    stop(paste("All", type, "cross-validation errors are non-finite. Please increase cc or use a larger p2."))
  }
  best <- which(err.mean == min(err.mean), arr.ind = TRUE)[1, ]
  list(p1 = p1.set[best[1]], p2 = p2.set[best[2]], err.mean = err.mean)
}

cv_central_huber <- function(x, y, p1.set, p2.set, tau.set = 1, ga, cc,
                             itr.max, k.fold = 5, early_stop_loss = Inf) {
  fold <- createFolds(1:nrow(x), k = k.fold, list = TRUE)
  err <- array(NA, c(k.fold, length(p1.set), length(p2.set), length(tau.set)))

  for (k in 1:k.fold) {
    cat("CV fold", k, "for central huber\n")
    val_id <- as.vector(unlist(fold[k]))
    tr_id <- setdiff(seq_len(nrow(x)), val_id)
    for (p1.idx in seq_along(p1.set)) {
      for (p2.idx in seq_along(p2.set)) {
        for (tau.idx in seq_along(tau.set)) {
          tau.use <- tau.set[tau.idx]
          b <- fit_central_huber(x[tr_id, , drop = FALSE], y[tr_id],
                                 p1 = p1.set[p1.idx], p2 = p2.set[p2.idx],
                                 tau = tau.use, ga = ga, cc = cc,
                                 itr.max = itr.max,
                                 early_stop_loss = early_stop_loss)
          if (isTRUE(attr(b, "early_stop"))) {
            err[k, p1.idx, p2.idx, tau.idx] <- Inf
            next
          }
          pred <- x[val_id, , drop = FALSE] %*% b[itr.max + 1, ]
          res <- y[val_id] - pred
          err[k, p1.idx, p2.idx, tau.idx] <- mean(res^2)
        }
      }
    }
  }

  err.mean <- apply(err, c(2, 3, 4), mean, na.rm = TRUE)
  err.mean[!is.finite(err.mean)] <- Inf
  if (all(is.infinite(err.mean))) {
    stop("All central Huber cross-validation errors are non-finite. Please increase cc or use a larger p2.")
  }
  best <- which(err.mean == min(err.mean), arr.ind = TRUE)[1, ]
  list(p1 = p1.set[best[1]], p2 = p2.set[best[2]], tau = tau.set[best[3]],
       err.mean = err.mean)
}
## Data preparation
communities <- read_communities_data()
groups <- make_semantic_partition(communities$xdat)
des <- build_design(communities$xdat, groups)
x.all <- scale(des$x, center = TRUE, scale = TRUE)
x.all[, "Intercept"] <- 1
y.all <- as.numeric(scale(communities$y, center = TRUE, scale = TRUE))

set.seed(1)
n.train <- 1500
id.train <- sample(seq_len(nrow(x.all)), n.train, replace = FALSE)
x.train <- x.all[id.train, , drop = FALSE]
x.test <- x.all[-id.train, , drop = FALSE]
y.train <- y.all[id.train]
y.test <- y.all[-id.train]

## Contaminate a fraction of the training responses only.
contam.rate <- 0.1
contam.shift <- 2
set.seed(1)
contam.id <- sample(seq_along(y.train), floor(contam.rate * length(y.train)),
                    replace = FALSE)
y.train.clean <- y.train
y.train[contam.id] <- y.train[contam.id] + contam.shift
cat("Contaminated", length(contam.id), "training responses by adding",
    contam.shift, "\n")

d <- ncol(x.train)
m <- 4
dj <- des$dj
djk <- des$djk

cat("Agent dimensions:\n")
print(data.frame(agent = names(groups), d_j = dj))

# A connected 4-agent communication graph.
A <- matrix(c(1, 1, 0, 1,
              1, 1, 1, 0,
              0, 1, 1, 1,
              1, 0, 1, 1), ncol = 4)

ga <- 0.5
cc <- 200
itr.huber <- 1000
itr.med <- 1000
itr.ls <- 1000
p1.set <- c(0.0008, 0.001, 0.003)
p2.set <- c(0.03, 0.05, 0.08)
tau.set <- c(0.75, 1, 1.25)
cv.stop.loss <- 5

## Central Huber with cross validation over (rho1, rho2, tau)
cv.cen <- cv_central_huber(x.train, y.train, p1.set, p2.set,
                           tau.set = tau.set, ga = ga, cc = cc,
                           itr.max = itr.huber,
                           early_stop_loss = cv.stop.loss)
p1.cen <- cv.cen$p1
p2.cen <- cv.cen$p2
tau.cen <- cv.cen$tau
cat("Selected CenHuber: p1 =", p1.cen, "p2 =", p2.cen,
    "tau =", tau.cen, "\n")

fit.cen <- fit_central_huber(x.train, y.train, p1.cen, p2.cen,
                             tau = tau.cen, ga = ga, cc = cc,
                             itr.max = itr.huber, return_loss = TRUE)
b.cen <- fit.cen$beta
loss.cen <- fit.cen$train_huber_loss
plot(loss.cen)

## VFLHuber with cross validation over (rho1, rho2, tau)
cv.huber <- cv_vfl_huber(x.train, y.train, dj, djk, A, p1.set, p2.set,
                   tau.set = tau.set, ga = ga, cc = cc, itr.max = itr.huber,
                   early_stop_loss = cv.stop.loss)
p1.huber <- cv.huber$p1
p2.huber <- cv.huber$p2
tau.huber <- cv.huber$tau
cat("Selected VFLHuber: p1 =", p1.huber, "p2 =", p2.huber,
    "tau =", tau.huber, "\n")

# p1.huber <- p1.cen
# p2.huber <- p2.cen
# tau.huber <- tau.cen
fit.hub <- fit_vfl(x.train, y.train, dj, djk, A, p1.huber, p2.huber,
                   tau = tau.huber, ga = ga, cc = cc, itr.max = itr.huber,
                   type = "huber", return_loss = TRUE)
b.hub <- fit.hub$beta
loss.hub <- fit.hub$train_huber_loss
plot(loss.hub)

fit.loc <- fit_local_huber(x.train, y.train, dj, djk, p1.huber, p2.huber,
                           tau = tau.huber, ga = ga, cc = 500,
                           itr.max = itr.huber, return_loss = TRUE)
b.loc <- fit.loc$beta
loss.loc <- fit.loc$train_huber_loss
plot(loss.loc)

## VFLMed with cross validation over (rho1, rho2)
cv.med <- cv_vfl_nonhuber(x.train, y.train, dj, djk, A, p1.set, p2.set,
                          ga = ga, cc = cc, itr.max = itr.med,
                          type = "median")
p1.med <- cv.med$p1
p2.med <- cv.med$p2
cat("Selected VFLMed: p1 =", p1.med, "p2 =", p2.med, "\n")
fit.med <- fit_vfl_nonhuber(x.train, y.train, dj, djk, A, p1.med, p2.med,
                            ga = ga, cc = cc, itr.max = itr.med,
                            type = "median", return_loss = TRUE)
b.med <- fit.med$beta
loss.med <- fit.med$train_loss
plot(loss.med)

## VFLLS with cross validation over (rho1, rho2)
cv.ls <- cv_vfl_nonhuber(x.train, y.train, dj, djk, A, p1.set, p2.set,
                         ga = ga, cc = cc, itr.max = itr.ls,
                         type = "ls")
p1.ls <- cv.ls$p1
p2.ls <- cv.ls$p2
cat("Selected VFLLS: p1 =", p1.ls, "p2 =", p2.ls, "\n")
fit.ols <- fit_vfl_nonhuber(x.train, y.train, dj, djk, A, p1.ls, p2.ls,
                            ga = ga, cc = cc, itr.max = itr.ls,
                            type = "ls", return_loss = TRUE)
b.ols <- fit.ols$beta
loss.ols <- fit.ols$train_loss
plot(loss.ols)

metric <- function(beta) {
  pred <- as.vector(x.test %*% beta)
  c(MSE = mean((pred - y.test)^2), MAE = mean(abs(pred - y.test)))
}

results <- rbind(
  CenHuber = metric(b.cen[itr.huber + 1, ]),
  VFLHuber = metric(b.hub[itr.huber + 1, ]),
  LocHuber = metric(b.loc[itr.huber + 1, ]),
  VFLMed = metric(b.med[itr.med + 1, ]),
  VFLLS = metric(b.ols[itr.ls + 1, ])
)

loss.path <- data.frame(
  iter = 0:itr.huber,
  CenHuber = loss.cen,
  VFLHuber = loss.hub,
  LocHuber = loss.loc,
  VFLMed = loss.med,
  VFLLS = loss.ols
)

print(round(results, 4))

