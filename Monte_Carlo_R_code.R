## Monte Carlo Simulation lecture - R code for all examples (Computational Statistics)
## Givens & Hoeting, Computational Statistics 2nd ed., Chapter 6

## ================= mc_pi.R =================
## pi = 4 * P(U1^2 + U2^2 <= 1),  (U1, U2) ~ Unif(0,1)^2
set.seed(1)
n <- 10000
x <- runif(n); y <- runif(n)
inside <- (x^2 + y^2 <= 1)            # indicator h(X)
pi.hat <- 4 * mean(inside)
se     <- 4 * sd(inside) / sqrt(n)    # Monte Carlo standard error
c(pi.hat = pi.hat, se = se)

## ================= mc_integral.R =================
## mu = int_0^1 exp(x^2) dx = E[exp(U^2)],  U ~ Unif(0,1)
set.seed(2)
n <- 100000
h <- exp(runif(n)^2)                  # h(U_i)
mu.hat <- mean(h)
se     <- sd(h) / sqrt(n)
c(mu.hat = mu.hat, se = se,
  lower = mu.hat - 1.96 * se, upper = mu.hat + 1.96 * se)
integrate(function(x) exp(x^2), 0, 1)$value   # 1.462652 (check)

## running mean: law of large numbers
run <- cumsum(h) / seq_along(h)
plot(run, type = "l", log = "x", xlab = "n", ylab = "running mean")
abline(h = 1.462652, lty = 2)

## LLN fails without E|X| < Inf : Cauchy
z <- rcauchy(n)
plot(cumsum(z) / seq_along(z), type = "l", log = "x")

## ================= inv_cdf.R =================
## Exponential(rate = lambda):  F^{-1}(u) = -log(1 - u) / lambda
rexp_inv <- function(n, lambda) -log(1 - runif(n)) / lambda

## Cauchy(0, 1):  F^{-1}(u) = tan(pi * (u - 1/2))
rcauchy_inv <- function(n) tan(pi * (runif(n) - 0.5))

## Discrete: generalized inverse  F^-(u) = min{x : F(x) >= u}
rdiscrete <- function(n, vals, prob) {
  cp <- cumsum(prob)
  vals[findInterval(runif(n), cp) + 1]
}

## Truncated normal on (a, b):  F^{-1}(F(a) + U (F(b) - F(a)))
rtnorm <- function(n, mu, sd, a, b) {
  Fa <- pnorm(a, mu, sd); Fb <- pnorm(b, mu, sd)
  qnorm(Fa + runif(n) * (Fb - Fa), mu, sd)
}

set.seed(4)
x <- rexp_inv(10000, lambda = 2)
c(mean = mean(x), var = var(x))                 # theory: 0.5, 0.25
d <- rdiscrete(10000, vals = 0:3, prob = c(0.1, 0.3, 0.4, 0.2))
table(d) / 10000
t <- rtnorm(10000, 0, 1, a = 1, b = 3)
c(mean = mean(t), min = min(t), max = max(t))

## ================= rej_beta.R =================
## Target f: Beta(2.7, 6.3),  envelope e(x) = M * g(x),  g = Unif(0,1)
f <- function(x) dbeta(x, 2.7, 6.3)
M <- optimize(f, c(0, 1), maximum = TRUE)$objective   # sup f(x)

rej_beta <- function(n) {
  out <- numeric(0); tries <- 0
  while (length(out) < n) {
    y <- runif(n)                      # 1. Y ~ g
    u <- runif(n)                      # 2. U ~ Unif(0,1)
    out <- c(out, y[u <= f(y) / M])    # 3. accept if U <= f(Y)/e(Y)
    tries <- tries + n
  }
  list(x = out[1:n], acc.rate = length(out) / tries)
}
set.seed(5)
r <- rej_beta(10000)
c(M = M, theory = 1 / M, acc.rate = r$acc.rate,
  mean = mean(r$x), true.mean = 2.7 / 9)

## ================= rej_post.R =================
## Poisson data, lognormal prior: log(lambda) ~ N(log 4, 0.5^2)
x <- c(8, 3, 4, 3, 1, 7, 2, 6, 2, 7)
logL <- function(lam) sum(x) * log(lam) - length(x) * lam   # log L + const
lam.hat <- mean(x)                         # MLE maximizes L(lambda | x)

## q(lam) = prior(lam) * L(lam)  <=  prior(lam) * L(lam.hat) = e(lam)
rej_post <- function(n) {
  out <- numeric(0); tries <- 0
  while (length(out) < n) {
    y <- exp(rnorm(n, log(4), 0.5))        # Y ~ prior
    u <- runif(n)
    keep <- u <= exp(logL(y) - logL(lam.hat))  # q(Y)/e(Y) = L(Y)/L(lam.hat)
    out <- c(out, y[keep]); tries <- tries + n
  }
  list(lambda = out[1:n], acc.rate = length(out) / tries)
}
set.seed(7)
r <- rej_post(10000)
c(acc.rate = r$acc.rate, post.mean = mean(r$lambda),
  quantile(r$lambda, c(0.025, 0.975)))

## ================= sir.R =================
## Sampling Importance Resampling (Rubin, 1987)
sir <- function(n, m, rg, dg, df) {
  y <- rg(m)                                  # 1. Y_1..Y_m ~ g
  w <- df(y) / dg(y)                          # 2. importance weights
  w <- w / sum(w)                             #    standardized
  x <- sample(y, n, replace = TRUE, prob = w) # 3. resample n values
  list(x = x, ess = 1 / sum(w^2))             # effective sample size
}

## Slash distribution: Y = Z / U,  Z ~ N(0,1), U ~ Unif(0,1)
dslash <- function(y) ifelse(y == 0, 1 / (2 * sqrt(2 * pi)),
                             (1 - exp(-y^2 / 2)) / (y^2 * sqrt(2 * pi)))
rslash <- function(n) rnorm(n) / runif(n)

set.seed(8)
bad  <- sir(5000, 100000, rnorm,  dnorm,  dslash)  # target slash, g = N(0,1)
good <- sir(5000, 100000, rslash, dslash, dnorm)   # target N(0,1), g = slash
c(ess.bad = bad$ess, ess.good = good$ess)
quantile(bad$x,  c(0.01, 0.99))   # slash: very heavy tails (true ~ +-40)
quantile(good$x, c(0.01, 0.99))   # N(0,1): +-2.33

## ================= sir_post.R =================
## SIR for the same posterior: g = prior, weights = likelihood
x <- c(8, 3, 4, 3, 1, 7, 2, 6, 2, 7)
logL <- function(lam) sum(x) * log(lam) - length(x) * lam
set.seed(9)
m <- 100000; n <- 5000
y  <- exp(rnorm(m, log(4), 0.5))           # Y_i ~ prior
lw <- logL(y)                              # log weights  (prior cancels)
w  <- exp(lw - max(lw)); w <- w / sum(w)   # stable standardization
lam <- sample(y, n, replace = TRUE, prob = w)
c(ess = 1 / sum(w^2), post.mean = mean(lam),
  quantile(lam, c(0.025, 0.975)))
