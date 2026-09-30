## Variance Reduction lecture - R code for all examples (Computational Statistics)
## Givens & Hoeting, Computational Statistics 2nd ed., Section 6.4

## ================= is_tail.R =================
## mu = P(Z > 4.5),  Z ~ N(0,1)      true: 1 - pnorm(4.5) = 3.398e-06
set.seed(1)
n <- 100000
## (a) naive Monte Carlo
z <- rnorm(n)
h <- as.numeric(z > 4.5)
c(naive = mean(h), se = sd(h) / sqrt(n))

## (b) importance sampling with g = N(4.5, 1)
x <- rnorm(n, mean = 4.5)
w <- dnorm(x) / dnorm(x, mean = 4.5)      # w = f/g = exp(-4.5 x + 4.5^2/2)
hw <- (x > 4.5) * w
c(IS = mean(hw), se = sd(hw) / sqrt(n), true = pnorm(4.5, lower.tail = FALSE))

## ================= network.R =================
## Network: A - a1..a4 - b1..b4 - B, 20 edges, each edge fails w.p. p = 0.05
## nodes: 1 = A, 2:5 = a1..a4, 6:9 = b1..b4, 10 = B
edges <- rbind(cbind(1, 2:5),                 # 1-4   : A - a_i
               cbind(2:4, 3:5),               # 5-7   : a_i - a_(i+1)
               cbind(2:5, 6:9),               # 8-11  : a_i - b_i
               c(2, 7), c(5, 8),              # 12-13 : a1 - b2, a4 - b3
               cbind(6:8, 7:9),               # 14-16 : b_i - b_(i+1)
               cbind(6:9, 10))                # 17-20 : b_i - B

## h(x) = 1 if A and B are NOT connected;  X: n x 20 matrix, 1 = edge works
h_net <- function(X) {
  reach <- matrix(FALSE, nrow(X), 10); reach[, 1] <- TRUE
  repeat {
    old <- reach
    for (k in 1:20) {
      u <- edges[k, 1]; v <- edges[k, 2]
      ok <- X[, k] == 1 & (reach[, u] | reach[, v])
      reach[ok, u] <- TRUE; reach[ok, v] <- TRUE
    }
    if (all(reach == old)) break
  }
  as.numeric(!reach[, 10])
}

## Importance sampling: simulate with failure prob p.star, reweight
net_is <- function(n, p = 0.05, p.star = 0.25) {
  B <- matrix(rbinom(n * 20, 1, p.star), n, 20)       # 1 = broken
  b <- rowSums(B)                                     # number of broken edges
  w <- (p / p.star)^b * ((1 - p) / (1 - p.star))^(20 - b)
  hw <- h_net(1 - B) * w
  c(est = mean(hw), se = sd(hw) / sqrt(n))
}
set.seed(2)
net_is(100000, p.star = 0.05)     # = naive Monte Carlo
net_is(100000, p.star = 0.25)     # importance sampling

## ================= antithetic.R =================
## Antithetic pairs via inverse CDF:  X = F^{-1}(U),  X' = F^{-1}(1 - U)
anti <- function(n, hq) {               # hq(u) = h(F^{-1}(u))
  u  <- runif(n)
  a  <- hq(u); b <- hq(1 - u)
  rho <- cor(a, b)
  est <- mean((a + b) / 2)
  c(est = est, se = sd((a + b) / 2) / sqrt(n), rho = rho,
    efficiency = 1 / (1 + rho))           # vs. naive MC with 2n draws
}
set.seed(3)
n <- 50000
## (1) monotone h, uniform:  E[exp(U^2)]
anti(n, function(u) exp(u^2))
## (2) monotone h, skewed Gamma(2,1) (not symmetric!):  E[log(1 + X)]
anti(n, function(u) log(1 + qgamma(u, shape = 2)))
## (3) non-monotone h:  E[Z^2] with Z = qnorm(U)  ->  Z' = -Z, same h!
anti(n, function(u) qnorm(u)^2)
## (4) non-monotone h:  E[cos(2 pi U)]
anti(n, function(u) cos(2 * pi * u))

## cost of the inverse CDF for non-symmetric distributions
system.time(rgamma(1e6, shape = 2))            # fast (rejection-type method)
system.time(qgamma(runif(1e6), shape = 2))     # slow (numerical inversion)

## ================= cv.R =================
## (1) mu = E[exp(U^2)],  control variate c(U) = 1 + U^2,  E[c(U)] = 4/3
set.seed(4)
n  <- 10000
u  <- runif(n)
h  <- exp(u^2); cc <- 1 + u^2
lam <- -cov(h, cc) / var(cc)                   # estimated optimal lambda
cv  <- h + lam * (cc - 4/3)
c(naive = mean(h), se.naive = sd(h) / sqrt(n),
  CV = mean(cv), se.CV = sd(cv) / sqrt(n), lambda = lam, rho = cor(h, cc))
## same thing as regression:  coef(lm(h ~ I(cc - 4/3)))[1]

## (2) Arithmetic Asian call option; control = geometric Asian call
S0 <- 100; K <- 100; r <- 0.05; sig <- 0.3; T <- 1; m <- 12
dt <- T / m
geo_price <- function() {                      # closed form (lognormal)
  mu.G <- log(S0) + (r - sig^2 / 2) * T * (m + 1) / (2 * m)
  s.G  <- sig * sqrt(T * (m + 1) * (2 * m + 1) / (6 * m^2))
  d2 <- (mu.G - log(K)) / s.G; d1 <- d2 + s.G
  exp(-r * T) * (exp(mu.G + s.G^2 / 2) * pnorm(d1) - K * pnorm(d2))
}
Z  <- matrix(rnorm(n * m), n, m)
logS <- log(S0) + t(apply((r - sig^2 / 2) * dt + sig * sqrt(dt) * Z, 1, cumsum))
S  <- exp(logS)                                 # n x m price paths
ari <- exp(-r * T) * pmax(rowMeans(S) - K, 0)   # target payoff
geo <- exp(-r * T) * pmax(exp(rowMeans(logS)) - K, 0)   # control
lam <- -cov(ari, geo) / var(geo)
cv  <- ari + lam * (geo - geo_price())
c(naive = mean(ari), se.naive = sd(ari) / sqrt(n),
  CV = mean(cv), se.CV = sd(cv) / sqrt(n), lambda = lam, rho = cor(ari, geo))

## ================= rb.R =================
## (1) t tail probability  P(T > 3),  T = Z / sqrt(V / nu),  V ~ chi^2_nu
## Rao-Blackwell:  E[ I(T > 3) | V ] = 1 - Phi(3 * sqrt(V / nu))
set.seed(5)
n <- 10000; nu <- 5
V <- rchisq(n, df = nu); Z <- rnorm(n)
h.naive <- as.numeric(Z / sqrt(V / nu) > 3)
h.rb    <- pnorm(3 * sqrt(V / nu), lower.tail = FALSE)
c(naive = mean(h.naive), se = sd(h.naive) / sqrt(n),
  RB = mean(h.rb), se.RB = sd(h.rb) / sqrt(n),
  true = pt(3, nu, lower.tail = FALSE))

## (2) Compound Poisson (aggregate claims): S = X_1 + ... + X_N,
##     N ~ Poisson(5), X_j ~ Exp(1).  Estimate P(S > 12).
## Rao-Blackwell:  S | N ~ Gamma(N, 1)  ->  E[I(S > 12) | N] = P(Gamma(N,1) > 12)
N  <- rpois(n, 5)
S  <- sapply(N, function(k) sum(rexp(k)))
h.naive <- as.numeric(S > 12)
h.rb    <- ifelse(N == 0, 0, pgamma(12, shape = pmax(N, 1), lower.tail = FALSE))
true <- sum(dpois(1:100, 5) * pgamma(12, 1:100, lower.tail = FALSE))
c(naive = mean(h.naive), se = sd(h.naive) / sqrt(n),
  RB = mean(h.rb), se.RB = sd(h.rb) / sqrt(n), true = true)
