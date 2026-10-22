rm(list = ls())
gc()

setwd("C:\\Users\\JUNZHUO GAO\\Desktop\\ADMM_Huber_decen_vert")

library(MASS)
library(ggplot2)
library(ggpubr)

thres <- function(a, t) {
  sign(a) * pmax(abs(a) - t, 0)
}

clip_huber_score <- function(a, tau) {
  pmin(pmax(a, -tau), tau)
}

huber_loss <- function(x, tau) {
  (abs(x) <= tau) * x^2 / 2 + (abs(x) > tau) * (tau * abs(x) - tau^2 / 2)
}

make_ring_graph <- function(m, degree_each_side = 1) {
  A <- matrix(0, m, m)
  for (i in 1:m) {
    for (h in 1:degree_each_side) {
      j1 <- ((i - 1 + h) %% m) + 1
      j2 <- ((i - 1 - h) %% m) + 1
      A[i, j1] <- A[j1, i] <- 1
      A[i, j2] <- A[j2, i] <- 1
    }
  }
  diag(A) <- 1
  A
}

make_star_graph <- function(m, hub = 1, leaf_degree_each_side = 0) {
  A <- matrix(0, m, m)
  leaves <- setdiff(1:m, hub)
  A[hub, leaves] <- 1
  A[leaves, hub] <- 1
  
  if (leaf_degree_each_side > 0) {
    n_leaf <- length(leaves)
    max_degree_each_side <- floor(n_leaf / 2)
    if (leaf_degree_each_side > max_degree_each_side) {
      stop("leaf_degree_each_side is too large for the number of leaves.")
    }
    for (idx in seq_along(leaves)) {
      for (h in 1:leaf_degree_each_side) {
        j1 <- leaves[((idx - 1 + h) %% n_leaf) + 1]
        j2 <- leaves[((idx - 1 - h) %% n_leaf) + 1]
        A[leaves[idx], j1] <- A[j1, leaves[idx]] <- 1
        A[leaves[idx], j2] <- A[j2, leaves[idx]] <- 1
      }
    }
  }
  
  diag(A) <- 1
  A
}

graph_is_connected <- function(A) {
  m <- nrow(A)
  seen <- rep(FALSE, m)
  queue <- 1
  seen[1] <- TRUE
  while (length(queue) > 0) {
    cur <- queue[1]
    queue <- queue[-1]
    nb <- setdiff(which(A[cur, ] != 0), cur)
    new_nb <- nb[!seen[nb]]
    if (length(new_nb) > 0) {
      seen[new_nb] <- TRUE
      queue <- c(queue, new_nb)
    }
  }
  all(seen)
}

average_degree <- function(A) {
  mean(rowSums(A) - diag(A))
}

generate_data <- function(n = 200, d = 500, s = 5, rhox = 0.5,
                          noise_sd = 1, seed = 1) {
  mux <- rep(0, d)
  sigmax <- outer(1:d, 1:d, function(i, j) rhox^abs(i - j))
  bt <- c(rep(1, s), rep(0, d - s))
  set.seed(seed)
  x <- mvrnorm(n, mux, sigmax)
  e <- rnorm(n, 0, noise_sd)
  y <- as.vector(x %*% bt + e)
  list(x = x, y = y, bt = bt, sigmax = sigmax)
}

run_central_huber <- function(x, y, p1, p2, tau, gamma = 0.5,
                              c = 10, itr_max = 1000) {
  n <- nrow(x)
  d <- ncol(x)
  b <- matrix(0, itr_max + 1, d)
  u <- z <- v <- matrix(0, itr_max + 1, n)
  obj <- rep(NA, itr_max + 1)
  obj[1] <- mean(huber_loss(y - x %*% b[1, ], tau)) +
    p2 * sum(b[1, ]^2) + p1 * sum(abs(b[1, ]))

  for (itr in 1:itr_max) {
    b[itr + 1, ] <- thres(-t(x) %*% u[itr, ] / (2 * n * p2),
                          p1 / (2 * p2))
    obj[itr + 1] <- mean(huber_loss(y - x %*% b[itr + 1, ], tau)) +
      p2 * sum(b[itr + 1, ]^2) + p1 * sum(abs(b[itr + 1, ]))
    u[itr + 1, ] <- (1 - gamma / c) * u[itr, ] +
      x %*% b[itr + 1, ] / c - v[itr, ] / c + gamma * z[itr, ] / c
    z[itr + 1, ] <- clip_huber_score(
      (u[itr + 1, ] + v[itr, ] / gamma - y / gamma) / (1 + 1 / gamma),
      tau
    )
    v[itr + 1, ] <- v[itr, ] + gamma * (u[itr + 1, ] - z[itr + 1, ])
  }
  list(beta = b[itr_max + 1, ], obj = obj)
}

run_vfl_huber <- function(x, y, bt, beta_ref, A, p1, p2, tau,
                          gamma = 0.5, c = 10, itr_max = 800) {
  if (!graph_is_connected(A)) {
    stop("The supplied graph is not connected.")
  }
  n <- nrow(x)
  d <- ncol(x)
  m <- nrow(A)
  dj <- d / m
  if (dj != floor(dj)) {
    stop("d must be divisible by m.")
  }

  xm <- array(NA, c(n, m, dj))
  for (j in 1:m) {
    xm[, j, ] <- x[, ((j - 1) * dj + 1):(j * dj)]
  }

  stat_err <- opt_err <- obj <- re_obj <- rep(NA, itr_max + 1)
  b <- array(0, c(itr_max + 1, m, dj))
  u <- v <- alpha <- array(0, c(itr_max + 1, m, n))

  beta0 <- as.vector(t(b[1, , ]))
  ref_obj <- mean(huber_loss(y - x %*% beta_ref, tau)) +
    p2 * sum(beta_ref^2) + p1 * sum(abs(beta_ref))

  stat_err[1] <- sqrt(sum((beta0 - bt)^2))
  opt_err[1] <- sqrt(sum((beta0 - beta_ref)^2))
  obj[1] <- mean(huber_loss(y - x %*% beta0, tau)) +
    p2 * sum(beta0^2) + p1 * sum(abs(beta0))
  re_obj[1] <- obj[1] - ref_obj

  for (itr in 1:itr_max) {
    for (j in 1:m) {
      b[itr + 1, j, ] <- thres(
        -t(xm[, j, ]) %*% u[itr, j, ] / (2 * n * p2),
        p1 / (2 * p2)
      )
    }

    beta_now <- as.vector(t(b[itr + 1, , ]))
    stat_err[itr + 1] <- sqrt(sum((beta_now - bt)^2))
    opt_err[itr + 1] <- sqrt(sum((beta_now - beta_ref)^2))
    obj[itr + 1] <- mean(huber_loss(y - x %*% beta_now, tau)) +
      p2 * sum(beta_now^2) + p1 * sum(abs(beta_now))
    re_obj[itr + 1] <- obj[itr + 1] - ref_obj

    for (j in 1:m) {
      nb <- which(A[j, ] != 0)
      neigh_u <- matrix(u[itr, nb, ], nrow = length(nb))
      su <- length(nb) * u[itr, j, ] - colSums(neigh_u)
      v[itr + 1, j, ] <- u[itr, j, ] -
        (-xm[, j, ] %*% b[itr + 1, j, ] + gamma * su + alpha[itr, j, ]) / c
    }

    for (j in 1:m) {
      nb <- which(A[j, ] != 0)
      neigh_v <- matrix(v[itr + 1, nb, ], nrow = length(nb))
      sv <- length(nb) * v[itr + 1, j, ] - colSums(neigh_v)
      alpha[itr + 1, j, ] <- alpha[itr, j, ] + gamma * sv
      u[itr + 1, j, ] <- clip_huber_score(
        (v[itr + 1, j, ] - y / (m * c)) / (1 + 1 / (m * c)),
        tau
      )
    }
  }

  data.frame(
    iter = 0:itr_max,
    StatErr = stat_err,
    OptErr = opt_err,
    ReObj = re_obj
  )
}

plot_convergence <- function(results, legend_title = "graph") {
  tab <- do.call(rbind, results)
  tab$logStatErr <- log(tab$StatErr)
  tab$logOptErr <- log(tab$OptErr)
  tab$logReObj <- log(pmax(tab$ReObj, .Machine$double.eps))

  base_theme <- theme_bw() +
    theme(
      panel.border = element_rect(color = "black", fill = NA),
      axis.title = element_text(size = 25),
      axis.text.x = element_text(size = 25),
      axis.text.y = element_text(size = 25),
      legend.text = element_text(size = 25),
      legend.title = element_text(size = 25),
      plot.title = element_text(size = 25, hjust = 0.5)
    )

  p1 <- ggplot(tab, aes(iter, logStatErr, color = label, linetype = label, shape = label)) +
    geom_line(linewidth = 3) +
    geom_point(data = function(x) x[x$iter %% 50 == 0, ], size = 5) +
    xlab(expression(t)) +
    ylab(expression(log(StatErr))) +
    labs(color = legend_title, linetype = legend_title, shape = legend_title) +
    base_theme

  p2 <- ggplot(tab, aes(iter, logOptErr, color = label, linetype = label, shape = label)) +
    geom_line(linewidth = 3) +
    geom_point(data = function(x) x[x$iter %% 50 == 0, ], size = 5) +
    xlab(expression(t)) +
    ylab(expression(log(OptErr))) +
    labs(color = legend_title, linetype = legend_title, shape = legend_title) +
    base_theme

  p3 <- ggplot(tab, aes(iter, logReObj, color = label, linetype = label, shape = label)) +
    geom_line(linewidth = 3) +
    geom_point(data = function(x) x[x$iter %% 50 == 0, ], size = 5) +
    xlab(expression(t)) +
    ylab(expression(log(ReObj))) +
    labs(color = legend_title, linetype = legend_title, shape = legend_title) +
    base_theme

  ggarrange(p1, p2, p3, nrow = 1, ncol = 3, common.legend = TRUE)
}

n <- 200
d <- 500
m <- 25
p1 <- 0.1
p2 <- 0.1
tau <- 0.5
gamma <- 0.5
itr_max <- 800

dat <- generate_data(n = n, d = d, seed = 1)
x <- dat$x
y <- dat$y
bt <- dat$bt

cen <- run_central_huber(x, y, p1 = p1, p2 = p2, tau = tau,
                         gamma = gamma, c = 20, itr_max = 1000)

# Part 1: compare several step sizes under the star-like graph.
c_set <- c(2, 10, 25, 50)
c_diverge <- 2
c_diverge_itr <- 15
c_diverge_plot_max <- 10
r_fixed <- 0
graph_star <- make_star_graph(m, hub = 2, leaf_degree_each_side = r_fixed)

step_results <- list()
step_plot_results <- list()
for (cc in c_set) {
  itr_use <- if (cc == c_diverge) c_diverge_itr else itr_max
  tmp <- run_vfl_huber(x, y, bt, cen$beta, graph_star,
                       p1 = p1, p2 = p2, tau = tau,
                       gamma = gamma, c = cc, itr_max = itr_use)
  tmp$label <- if (cc == c_diverge) " 2" else paste0("", cc)
  step_results[[tmp$label[1]]] <- tmp
  if (cc == c_diverge) {
    step_plot_results[[tmp$label[1]]] <- tmp[tmp$iter <= c_diverge_plot_max, ]
  } else {
    step_plot_results[[tmp$label[1]]] <- tmp
  }
}

p_step <- plot_convergence(step_plot_results, legend_title = "c")
p_step
ggsave("star_stepsize.eps", plot = p_step, device = "eps",
       width = 18, height = 6)

# Part 2: compare star-like graphs with different leaf connectivities.
# Tune c_degree separately for each topology to keep every graph stable
# while making the convergence-speed difference visible.
r_set <- c(0, 5, 9)
graphs_degree <- list()
for (rr in r_set) {
  graphs_degree[[paste0(rr)]] <- make_star_graph(m, hub = 2,
                                                 leaf_degree_each_side = rr)
}

c_degree <- c(
  "0" = 25,
  "5" = 20,
  "9" = 15
)

degree_results <- list()
for (gname in names(graphs_degree)) {
  c_use <- c_degree[[gname]]
  tmp <- run_vfl_huber(x, y, bt, cen$beta, graphs_degree[[gname]],
                       p1 = p1, p2 = p2, tau = tau,
                       gamma = gamma, c = c_use, itr_max = itr_max)
  tmp$label <- gname
  tmp$c <- c_use
  degree_results[[gname]] <- tmp
}

p_degree <- plot_convergence(degree_results, legend_title = "r")
p_degree


ggsave("star_degree.eps", plot = p_degree, device = "eps",
       width = 18, height = 6)

