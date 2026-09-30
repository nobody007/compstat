## MCMC lecture - R code for all examples (Computational Statistics)
## Givens & Hoeting, Computational Statistics 2nd ed., Ch. 1.7 & 7

## ================= markov.R =================
## 3-state Markov chain: 1 = Sunny, 2 = Cloudy, 3 = Rainy
P <- matrix(c(0.7, 0.2, 0.1,
              0.3, 0.4, 0.3,
              0.2, 0.3, 0.5), 3, 3, byrow = TRUE)
## stationary distribution: left eigenvector of P with eigenvalue 1
e <- eigen(t(P))
pi.st <- Re(e$vectors[, 1]); pi.st <- pi.st / sum(pi.st)
round(pi.st, 4)
## pi_t = pi_0 P^t  converges to pi.st from any start
p0 <- c(0, 0, 1)
for (t in 1:10) p0 <- p0 %*% P
round(p0, 4)
## ergodic theorem: long-run frequencies of ONE simulated path
set.seed(1)
n <- 100000; x <- numeric(n); x[1] <- 3
for (t in 2:n) x[t] <- sample(1:3, 1, prob = P[x[t - 1], ])
round(table(x) / n, 4)

## ================= mh_indep.R =================
## Data: y_i ~ delta N(7, 0.5^2) + (1 - delta) N(10, 0.5^2),  n = 100, delta = 0.7
set.seed(2)
y <- ifelse(runif(100) < 0.7, rnorm(100, 7, 0.5), rnorm(100, 10, 0.5))
## log posterior of delta with Unif(0,1) prior (up to a constant)
logpost <- function(d) sum(log(d * dnorm(y, 7, 0.5) + (1 - d) * dnorm(y, 10, 0.5)))

## Independence chain Metropolis-Hastings, proposal g = Beta(a, b)
mh_indep <- function(N, a, b, x0 = 0.5) {
  x <- numeric(N); x[1] <- x0; acc <- 0
  for (t in 2:N) {
    xs <- rbeta(1, a, b)                                   # candidate
    logR <- logpost(xs) - logpost(x[t - 1]) +
            dbeta(x[t - 1], a, b, log = TRUE) - dbeta(xs, a, b, log = TRUE)
    if (log(runif(1)) < logR) { x[t] <- xs; acc <- acc + 1 }
    else x[t] <- x[t - 1]
  }
  list(x = x, acc.rate = acc / (N - 1))
}
good <- mh_indep(10000, 1, 1)      # Beta(1,1) = Unif(0,1)
bad  <- mh_indep(10000, 2, 10)     # Beta(2,10): mass near 0.17, far from target
c(good = good$acc.rate, bad = bad$acc.rate)
c(mean(good$x[-(1:500)]), mean(bad$x[-(1:500)]))

## ================= rw.R =================
## Random walk chain on u = logit(delta)  (same data / posterior as mh_indep.R)
set.seed(2)
y <- ifelse(runif(100) < 0.7, rnorm(100, 7, 0.5), rnorm(100, 10, 0.5))
logpost <- function(d) sum(log(d * dnorm(y, 7, 0.5) + (1 - d) * dnorm(y, 10, 0.5)))
## target in u:  f_U(u) = f_delta(expit(u)) * |d delta / du|   (Jacobian)
logpost_u <- function(u) { d <- plogis(u); logpost(d) + log(d * (1 - d)) }

rw_chain <- function(N, s, u0 = 0) {
  u <- numeric(N); u[1] <- u0; acc <- 0
  for (t in 2:N) {
    us <- u[t - 1] + rnorm(1, 0, s)                        # symmetric proposal
    if (log(runif(1)) < logpost_u(us) - logpost_u(u[t - 1])) {  # Metropolis ratio
      u[t] <- us; acc <- acc + 1
    } else u[t] <- u[t - 1]
  }
  list(delta = plogis(u), acc.rate = acc / (N - 1))
}
set.seed(3)
small <- rw_chain(10000, 0.02); mid <- rw_chain(10000, 0.5); big <- rw_chain(10000, 10)
c(small = small$acc.rate, mid = mid$acc.rate, big = big$acc.rate)

## ================= diagnostics.R =================
## Convergence diagnostics for the random walk chain (rw.R must be run first)
ess <- function(x) {                              # effective sample size
  r <- acf(x, lag.max = 200, plot = FALSE)$acf[-1]
  k <- which(r < 0.05)[1]; if (is.na(k)) k <- length(r)
  length(x) / (1 + 2 * sum(r[1:k]))
}
rhat <- function(chains) {                        # Gelman-Rubin (columns = chains)
  n <- nrow(chains); m <- ncol(chains)
  B <- n * var(colMeans(chains)); W <- mean(apply(chains, 2, var))
  sqrt(((n - 1) / n * W + B / n) / W)
}
set.seed(6)
starts <- c(-4, -1, 2, 5)                         # overdispersed starts (logit scale)
ch <- sapply(starts, function(u0) rw_chain(5000, 0.5, u0)$delta)
burn <- 1:500
c(Rhat = rhat(ch[-burn, ]), ESS.small = ess(small$delta), ESS.mid = ess(mid$delta),
  ESS.big = ess(big$delta))

## ================= gibbs_bvn.R =================
## Gibbs sampler for (X1, X2) ~ N2(0, [[1, rho], [rho, 1]])
## full conditionals:  X1 | X2 ~ N(rho X2, 1 - rho^2),  X2 | X1 ~ N(rho X1, 1 - rho^2)
gibbs_bvn <- function(N, rho, x0 = c(-3, 3)) {
  x <- matrix(NA, N, 2); x[1, ] <- x0
  for (t in 2:N) {
    x[t, 1] <- rnorm(1, rho * x[t - 1, 2], sqrt(1 - rho^2))
    x[t, 2] <- rnorm(1, rho * x[t, 1],     sqrt(1 - rho^2))   # uses NEW x1
  }
  x
}
set.seed(4)
g5  <- gibbs_bvn(5000, 0.5)
g99 <- gibbs_bvn(5000, 0.99)
c(cor(g5)[1, 2], cor(g99)[1, 2])
c(acf(g5[, 1], plot = FALSE)$acf[2], acf(g99[, 1], plot = FALSE)$acf[2])

## ================= gibbs_normal.R =================
## y_i ~ N(mu, sigma2);  priors mu ~ N(mu0, tau0^2),  sigma2 ~ InvGamma(a, b)
set.seed(5)
y <- rnorm(30, mean = 10, sd = 2); n <- length(y); ybar <- mean(y)
mu0 <- 0; tau02 <- 100; a <- 0.01; b <- 0.01
N <- 10000; mu <- numeric(N); s2 <- numeric(N); mu[1] <- 0; s2[1] <- 1
for (t in 2:N) {
  ## mu | sigma2, y  ~  Normal
  prec  <- 1 / tau02 + n / s2[t - 1]
  mu[t] <- rnorm(1, (mu0 / tau02 + n * ybar / s2[t - 1]) / prec, sqrt(1 / prec))
  ## sigma2 | mu, y  ~  InvGamma(a + n/2, b + sum((y - mu)^2)/2)
  s2[t] <- 1 / rgamma(1, a + n / 2, b + sum((y - mu[t])^2) / 2)
}
keep <- 1001:N                                   # burn-in 1000
c(mu = mean(mu[keep]), sigma = mean(sqrt(s2[keep])))
quantile(mu[keep], c(0.025, 0.975)); c(ybar = ybar, sd = sd(y))
