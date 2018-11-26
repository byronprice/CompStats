# [ MA 589: MCMC demo ]

# [ p ~ Beta(a, b) and X | p ~ Bin(n, p) ]
x <- 8; n <- 10 # #heads, #flips

comparebeta <- function (f, x, n, a, b, ...) {
  p <- f(x, n, a, b, ...)
  op <- par(mfrow = c(2, 1))
  plot(p, type = "l")
  hist(p, freq = FALSE)
  t <- seq(0, 1, length = 100)
  lines(t, dbeta(t, x + a, n - x + b))
  par(op)
}

# Independent chain sampler
mhcoin <- function (x, n, a, b, ns = 1000) {
  p <- numeric(ns)
  p[1] <- a / (a + b)  # prior mean
  for (i in 2:ns) {
    pc <- rbeta(1, a, b) # candidate
    logR <- x * log(pc) + (n - x) * log(1 - pc) -
      (x * log(p[i - 1]) + (n - x) * log(1 - p[i - 1]))
    if (logR > 0 || log(runif(1)) < logR) # accept?
      p[i] <- pc
    else # reject
      p[i] <- p[i - 1]
  }
  p
}

a <-  1; b <-  1; comparebeta(mhcoin, x, n, a, b) # ok
a <- 10; b <- 10; comparebeta(mhcoin, x, n, a, b) # reasonably ok
a <-  1; b <- 10; comparebeta(mhcoin, x, n, a, b) # gets stuck: prior far from lhood

# prior far from lhood: visualization
t <- seq(0, 1, length = 100)
plot(t, dbeta(t, a, b), type = "l", lty = 2,
     xlab = "p", ylab = "density") # prior
lines(t, dbeta(t, x, n - x)) # normalized likelihood
lines(t, dbeta(t, a + x, n - x + b), lwd = 2) # posterior

# Random walk sampler
rwcoin <- function (x, n, a, b, s = .01, ns = 1000) {
  p <- numeric(ns)
  p[1] <- a / (a + b)
  for (i in 2:ns) {
    pc <- p[i - 1] + runif(1, -s, s) # candidate
    if (pc < 0 || pc > 1) { # out of bounds?
      p[i] <- p[i - 1] # reject
      next
    }
    logR <- dbeta(pc, a + x, b + n - x, log = TRUE) -
      dbeta(p[i - 1], a + x, b + n - x, log = TRUE)
    if (logR > 0 || log(runif(1)) < logR) # accept?
      p[i] <- pc
    else
      p[i] <- p[i - 1]
  }
  p
}

comparebeta(rwcoin, x, n, a, b, .01) # random step is too small: slow drift
comparebeta(rwcoin, x, n, a, b, .50) # large step: gets stuck once in a while
comparebeta(rwcoin, x, n, a, b, .10) # "just right"



# [ X ~ N(mu, sigma) ]
mu <- c(5, -5)
sigma <- matrix(c(100, 40, 40, 25), nrow = 2)

# plot (x - mu)' * sigma^(-1) * (x - mu) = qchisq(alpha, 2)
ellipse <- function (mu, sigma, alpha = .95, ns = 100) {
  p <- matrix(nrow = 2, ncol = ns) # points
  t <- seq(0, 2 * pi, length = ns) # param (angle)
  e <- eigen(sigma)
  s <- sqrt(e$values * qchisq(alpha, 2))
  # scale
  p[1,] <- s[1] * cos(t)
  p[2,] <- s[2] * sin(t)
  e$vectors %*% p + mu
}

comparegaussian <- function (f, mu, sigma, ...) {
  x <- f(mu, sigma, ...)
  op <- par(mfrow = c(2, 2))
  plot(x[,1], type = "l", xlab = "sample", ylab = expression(X[1]))
  plot(x[,2], type = "l", xlab = "sample", ylab = expression(X[2]))
  plot(x[,1], x[,2], type = "l",
       xlab = expression(X[1]), ylab = expression(X[2]))
  plot(x[,1], x[,2],
       xlab = expression(X[1]), ylab = expression(X[2]),
       pch = 20, col = "gray50")
  p <- ellipse(mu, sigma)
  lines(p[1,], p[2,], lwd=2)
  par(op)
}


# Gibbs sampler
gibbsgaussian <- function (mu, sigma, ns = 1000) {
  x <- matrix(nrow = ns, ncol = 2)
  x[1,] <- mu
  for (i in 2:ns) {
    # x1 | x2
    m <- mu[1] + sigma[1, 2] / sigma[2, 2] * (x[i - 1, 2] - mu[2])
    s <- sqrt(sigma[1, 1] - sigma[1, 2] ^ 2 / sigma[2, 2])
    x[i, 1] <- rnorm(1, m, s)
    # x2 | x1
    m <- mu[2] + sigma[1, 2] / sigma[1, 1] * (x[i, 1] - mu[1])
    s <- sqrt(sigma[2, 2] - sigma[1, 2] ^ 2 / sigma[1, 1])
    x[i, 2] <- rnorm(1, m, s)
  }
  x
}

# Random walk
rwgaussian <- function (mu, sigma, srw = 1, ns = 1000) {
  C <- chol(sigma)
  x <- matrix(nrow = ns, ncol = 2)
  x[1,] <- mu
  for (i in 2:ns) {
    xc <- rnorm(2, x[i - 1,], srw) # xc ~ N(x[i-1,], srw^2 * I2)
    logR <- -.5 * crossprod(backsolve(C, xc - mu, transpose = TRUE)) +
      .5 * crossprod(backsolve(C, x[i - 1,] - mu, transpose = TRUE))
    if (logR > 0 || log(runif(1)) < logR) # accept?
      x[i,] <- xc
    else # reject
      x[i,] <- x[i - 1,]
  }
  x
}

comparegaussian(gibbsgaussian, mu, sigma) # nice mixing
comparegaussian(rwgaussian, mu, sigma, 1) # slow mixing
comparegaussian(rwgaussian, mu, sigma, 10) # large step: gets stuck occasionally
comparegaussian(rwgaussian, mu, sigma, 5) # ok
