# [ Hamiltonian MCMC: simple example ]

# Euclidean-Gaussian kinetics: mass matrix `M` is fixed
leapfrog <- function (grad_V, M, epsilon, T, q, p) {
  qt <- q; pt <- p
  pt <- pt - epsilon / 2 * grad_V(qt) # first half-step
  for (t in 1:(T - 1)) {
    qt <- qt + epsilon * pt / M
    pt <- pt - epsilon * grad_V(qt)
  }
  qt <- qt + epsilon * pt / M
  pt <- pt - epsilon / 2 * grad_V(qt) # last half-step
  list(q = qt, p = pt)
}

# Example
neg_log_gaussian <- function (mu, sigma2)
  function (x) sum((x - mu) ^ 2 / sigma2) / 2
grad_neg_log_gaussian <- function (mu, sigma2)
  function (x) (x - mu) / sigma2

sigma2 <- seq(.01, 1, by = .01) ^ 2
n <- length(sigma2)
mu <- rep(0, n)
V <- neg_log_gaussian(mu, sigma2)
grad_V <- grad_neg_log_gaussian(mu, sigma2)
M <- rep(1, n)
K <- neg_log_gaussian(0, M)
ns <- 1000

# RW
sigma_rw <- .01
q <- matrix(0, nrow = ns, ncol = n)
for (is in 2:ns) {
  q[is,] <- q[is - 1,]
  qc <- rnorm(n, q[is,], sigma_rw)
  logR <- -V(qc) + V(q[is,])
  if (logR >= 0 || log(runif(1)) < logR) q[is,] <- qc
}
qrw <- q
plot(qrw[, n], type = 'l') # trace
plot(apply(qrw, 2, sd), colMeans(qrw))
plot(sqrt(sigma2), apply(qrw, 2, sd)); abline(0, 1)

# Hamiltonian MC
epsilon <- .01; T <- 150
q <- matrix(0, nrow = ns, ncol = n)
for (is in 2:ns) {
  q[is,] <- q[is - 1,]
  pc <- rnorm(n, 0, sqrt(M))
  lf <- leapfrog(grad_V, M, epsilon, T, q[is,], pc)
  logR <- -(V(lf$q) + K(lf$p)) + (V(q[is,]) + K(pc))
  if (logR >= 0 || log(runif(1)) < logR) q[is,] <- lf$q
}
qhmc <- q
plot(qhmc[, n], type = 'l') # trace
plot(apply(qhmc, 2, sd), colMeans(qhmc), ylim = range(colMeans(qrw)))
plot(sqrt(sigma2), apply(qhmc, 2, sd)); abline(0, 1)

