expit <- function(a) 1 / (1 + exp(-a))
logit <- function(a) log(a) - log(1-a)

lambda <- function(w, beta, kappa) c(kappa + 2 * (1 - kappa) * expit(w %*% beta - beta[3]))

G0 <- function(u, lambda) c(lambda * u + (1 - lambda) * u^2)

G0inv <- function(u, lambda) {
  if(length(lambda) == 1) {
    if(lambda == 1) return(u)
    (-lambda + sqrt(lambda^2 + 4 * u * (1 - lambda)))/(2 * (1 -lambda))
  } else {
    ret <- rep(NA, length(u))
    ones <- lambda == 1
    ret[ones] <- u[ones]
    ret[!ones] <- (-lambda[!ones] + sqrt(lambda[!ones]^2 + 4 * u[!ones] * (1 - lambda[!ones])))/(2 * (1 -lambda[!ones]))
    return(ret)
  }
}

f0 <- function(a, p.0=.2, p.5=.2, p.1=.2) {
  ifelse(a == 0, p.0,
         ifelse(a == .5, p.5,
                ifelse(a == 1, p.1,
                       ifelse(0 < a & a < 1, (1-p.0-p.5-p.1) * dbeta(a, 2,2),
                              0))))
}

F0 <- function(a, p.0=.2, p.5=.2, p.1=.2) {
  (a >= 0) * p.0 + (a >= .5) * p.5 + (a >= 1) * p.1 + (1-p.0-p.5-p.1) * pbeta(a, 2,2)
}

F0inv <- function(u, p.0=.2, p.5=.2, p.1=.2) {
  tmp <- rep(-Inf, length(u))
  tmp[0 <= u & u <= p.0] <- 0
  vals <- p.0 < u & u < p.0 + (1-p.0-p.5-p.1) * pbeta(.5, 2, 2)
  tmp[vals] <- qbeta((u[vals] - p.0) / (1-p.0-p.5-p.1), 2,2)
  vals <- p.0 + (1-p.0-p.5-p.1) * pbeta(.5, 2, 2) <= u & u <= p.0 + p.5 + (1-p.0-p.5-p.1) * pbeta(.5, 2, 2)
  tmp[vals] <- .5
  vals <- p.0 + p.5 + (1-p.0-p.5-p.1) * pbeta(.5, 2, 2) < u & u < 1 - p.1
  tmp[vals] <- qbeta((u[vals] - p.0 - p.5) / (1-p.0-p.5-p.1), 2, 2)
  vals <- 1 - p.1 <= u
  tmp[vals] <- 1
  return(tmp)
}

g0 <- function(a, lambda, p.0=.2, p.5=.2, p.1=.2) {
  ifelse(a == 0, (G0(F0(0), lambda) - G0(0, lambda)) / p.0,
         ifelse(a == .5, (G0(F0(.5), lambda) - G0(F0(.5) - p.5, lambda)) / p.5,
                ifelse(a == 1, (G0(F0(1), lambda) - G0(F0(1) - p.1, lambda)) / p.1,
                       c(lambda + 2 * (1 - lambda) * a))))
}

mu0 <- function(a, w, gamma1, gamma2, gamma3) {
  c(gamma1[1] + w %*% gamma1[-1] + a * gamma2[1] + a * (w %*% gamma2[-1]) + (a - .5)^2 * gamma3)
}


theta0 <- function(a, gamma1, gamma2, gamma3) {
  gamma1[1] + gamma1[4] + a * (gamma2[1] + gamma2[4]) + (a - .5)^2 * gamma3
}

kappa <- .1
beta <- c(-1,1,-1)
gamma1 <- c(0, 2, 2, -2)
gamma2 <- c(1, 1, -1, -.5)
gamma3 <- 0

n <- 500

W <- matrix(c(rnorm(n), rnorm(n), rnorm(n, 1, 1)), ncol=3)
lambda.vals <- lambda(W, beta=beta, kappa=kappa)
U0 <- G0inv(runif(n), lambda.vals)
A <- F0inv(U0, p.0 = .2, p.5=.2, p.1=.2)
g0s <- g0(A, lambda.vals, p.0 = .2, p.5 = .2, p.1 = .2)
mu0s <- mu0(A, W, gamma1, gamma2, gamma3)
Y <- rnorm(n, mean = mu0s, sd = sqrt(1 + abs(mu0s) ))

V <- 5
folds <- sample(rep(1:V, length.out=n), n, replace=FALSE)
cv.folds <- lapply(1:V, function(j) which(folds == j))

g.hats.np.cv <- lapply(1:V, function(j) {
  g.fit <- con.mixed.dens.SL(A=A[folds != j], W=W[folds != j,], SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"), n.bins = 2:10, verbose = TRUE)
  c(predict.con.mixed.dens.SL(g.fit, new.A=A[folds == j], new.W=W[folds == j,]))
})
g.hats.np <- rep(NA, n)
for(j in 1:V) g.hats.np[folds == j] <- g.hats.np.cv[[j]]

plot(g0s, g.hats.np)
abline(0,1, col='red')

library(SuperLearner)

df <- data.frame(A, W)
names(df) <- c("A", "W1", "W2", "W3")

screen.noA <- function (Y, X, ...) {
  whichVariable <- rep(TRUE, ncol(X))
  whichVariable[names(X) == "A"] <- FALSE
  return(whichVariable)
}

library <- lapply(c("SL.mean", "SL.glm", "SL.gam", "SL.earth", "SL.glm.interaction"), function(alg) c(alg, "All", "screen.noA"))

mu.fits.np.cv <- lapply(1:V, function(j) {
  mu.fit <- SuperLearner(Y[folds != j], df[folds != j,], family = 'gaussian', SL.library = library, method = 'method.NNLS', verbose = FALSE)
  function(a, w) {
    newdf <- data.frame(A=a, w)
    names(newdf) <- c("A", "W1", "W2", "W3")
    c(predict(mu.fit, newdata=newdf)$pred)
  }
})

test.np <- causal.null.test(Y = Y, A = A, W = W, g.hats = g.hats.np.cv, mu.hat = mu.fits.np.cv, p=c(1,2,Inf), cv.folds = cv.folds)
