## EM Algorithm lecture - R code for all examples (Computational Statistics)

## ================= gmm.R =================
em_gmm <- function(y, p = 0.5, mu = c(-1, 1), s = c(2, 2),
                   tol = 1e-8, maxit = 500) {
  ll_old <- -Inf
  for (t in 1:maxit) {
    ## E-step: responsibilities gamma_i
    d1 <- p * dnorm(y, mu[1], s[1])
    d2 <- (1 - p) * dnorm(y, mu[2], s[2])
    g  <- d1 / (d1 + d2)
    ## M-step: weighted MLEs
    p     <- mean(g)
    mu[1] <- sum(g * y) / sum(g)
    mu[2] <- sum((1 - g) * y) / sum(1 - g)
    s[1]  <- sqrt(sum(g * (y - mu[1])^2) / sum(g))
    s[2]  <- sqrt(sum((1 - g) * (y - mu[2])^2) / sum(1 - g))
    ll <- sum(log(p * dnorm(y, mu[1], s[1]) +
                  (1 - p) * dnorm(y, mu[2], s[2])))
    if (abs(ll - ll_old) < tol) break
    ll_old <- ll
  }
  list(pi = p, mu = mu, sigma = s, loglik = ll, iter = t)
}
set.seed(2026)
z <- rbinom(300, 1, 0.4)
y <- ifelse(z == 1, rnorm(300, 0, 1), rnorm(300, 4, 1.5))
str(em_gmm(y))

## ================= bvn.R =================
em_bvn <- function(Y, tol = 1e-8) {         # Y: n x 2, NA = missing
  m1 <- is.na(Y[, 1]); m2 <- is.na(Y[, 2]); n <- nrow(Y)
  mu <- colMeans(Y, na.rm = TRUE); S <- var(Y, use = "complete.obs")
  repeat {
    E <- Y; V <- matrix(0, n, 2)             # V: conditional variances
    ## E-step
    E[m2, 2] <- mu[2] + S[1,2]/S[1,1] * (Y[m2, 1] - mu[1])
    V[m2, 2] <- S[2,2] - S[1,2]^2 / S[1,1]
    E[m1, 1] <- mu[1] + S[1,2]/S[2,2] * (Y[m1, 2] - mu[2])
    V[m1, 1] <- S[1,1] - S[1,2]^2 / S[2,2]
    ## M-step
    mu.new <- colMeans(E)
    S.new  <- crossprod(E)/n + diag(colMeans(V)) - tcrossprod(mu.new)
    if (max(abs(c(mu.new - mu, S.new - S))) < tol) break
    mu <- mu.new; S <- S.new
  }
  list(mu = mu.new, Sigma = S.new)
}
set.seed(7); library(stats)
Z <- matrix(rnorm(400),200,2) %*% chol(matrix(c(4,3,3,9),2)); Z[,1]<-Z[,1]+5; Z[,2]<-Z[,2]+10
Z[sample(200,40),2] <- NA; Z[sample(which(!is.na(Z[,2])),15),1] <- NA
print(em_bvn(Z))

## ================= abo.R =================
abo_em <- function(nA = 186, nB = 38, nAB = 13, nO = 284,
                   p = 1/3, q = 1/3, r = 1/3, tol = 1e-10) {
  n <- nA + nB + nAB + nO
  repeat {
    ## E-step: split phenotype counts into genotypes
    nAA <- nA * p^2 / (p^2 + 2*p*r);  nAO <- nA - nAA
    nBB <- nB * q^2 / (q^2 + 2*q*r);  nBO <- nB - nBB
    ## M-step: gene counting
    new <- c(2*nAA + nAO + nAB,
             2*nBB + nBO + nAB,
             2*nO  + nAO + nBO) / (2*n)
    if (max(abs(new - c(p, q, r))) < tol) break
    p <- new[1]; q <- new[2]; r <- new[3]
  }
  c(p = new[1], q = new[2], r = new[3])
}
abo_em()

## ================= mlr.R =================
em_mlr <- function(X, y, tol = 1e-10) {    # y: NA = missing response
  o <- !is.na(y); n <- length(y); m <- sum(o)
  b <- rep(0, ncol(X)); s2 <- 1
  XtXi <- solve(crossprod(X))
  repeat {
    yh <- ifelse(o, y, X %*% b)                     # E-step
    b.new  <- drop(XtXi %*% crossprod(X, yh))       # M-step
    s2.new <- (sum((yh - X %*% b.new)^2) + (n - m) * s2) / n
    if (max(abs(c(b.new - b, s2.new - s2))) < tol) break
    b <- b.new; s2 <- s2.new
  }
  list(beta = b.new, sigma2 = s2.new)
}
set.seed(1); X <- cbind(1, rnorm(100), runif(100,0,2)); y <- drop(X%*%c(1,2,-1)) + rnorm(100)
y[71:100] <- NA
print(em_mlr(X,y)); f<-lm(y[1:70]~X[1:70,-1]); print(coef(f)); print(sum(resid(f)^2)/70)

## ================= mcem.R =================
mcem_cens <- function(x, cens, c, mu = 0, s = 2,
                      T = 40, m0 = 10, inc = 1.25) {
  n <- length(x); k <- sum(cens); path <- mu
  for (t in 1:T) {
    m <- ceiling(m0 * inc^t)                  # MC sample size grows
    ## MC E-step: Z | Z > c  by inverse-cdf sampling
    a <- pnorm(c, mu, s)
    Z <- matrix(qnorm(runif(k * m, a, 1), mu, s), k, m)
    S1 <- sum(x[!cens])   + sum(rowMeans(Z))
    S2 <- sum(x[!cens]^2) + sum(rowMeans(Z^2))
    ## M-step (closed form)
    mu <- S1 / n; s <- sqrt(S2 / n - mu^2); path <- c(path, mu)
  }
  list(mu = mu, sigma = s, path = path)
}
set.seed(3); x0 <- rnorm(100, 1, 1); cens <- x0 > 1.5; x <- pmin(x0, 1.5)
r <- mcem_cens(x, cens, 1.5); print(c(r$mu, r$sigma))

## ================= ect.R =================
ecm_t <- function(y, mu = median(y), s2 = var(y), nu = 10,
                  tol = 1e-9, maxit = 1000) {
  ll <- function(mu, s2, nu)
    sum(dt((y - mu)/sqrt(s2), nu, log = TRUE) - 0.5*log(s2))
  old <- ll(mu, s2, nu)
  for (t in 1:maxit) {
    d2 <- (y - mu)^2 / s2                            # E-step
    w  <- (nu + 1) / (nu + d2)
    El <- digamma((nu + 1)/2) - log((nu + d2)/2)
    mu <- sum(w * y) / sum(w)                        # CM-step 1
    s2 <- mean(w * (y - mu)^2)
    f  <- function(v) -digamma(v/2) + log(v/2) + 1 + mean(El - w)
    nu <- uniroot(f, c(0.05, 1000))$root             # CM-step 2
    new <- ll(mu, s2, nu)
    if (abs(new - old) < tol) break
    old <- new
  }
  c(mu = mu, sigma = sqrt(s2), nu = nu, iter = t)
}
set.seed(5); y <- 2 + 1.5 * rt(300, df = 4); print(ecm_t(y))
