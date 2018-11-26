# Bayesian simple logistic regression using Laplace trick as proposal
# y_i ~ind Bern[logit^(-1)(x_i * beta + offset_i)]
# beta ~ N(beta0, omega^(-1))

logit <- function (p) log(p / (1 - p))
inv_logit <- function (x) 1 / (1 + exp(-x))
LOGEPS <- log(.Machine$double.eps / 2)
log1pe <- function (x) {
  l <- ifelse(x > 0, x, 0)
  x <- ifelse(x > 0, -x, x)
  ifelse(x < LOGEPS, l, l + log1p(exp(x)))
}

log_posterior <- function (beta, y, x, offset, beta0 = 0, omega = 0) {
  eta <- offset + x * beta
  sum(y * eta - log1pe(eta)) - .5 * omega * (beta - beta0) ^ 2
}

step_laplace_mh <- function (beta, y, x, offset, beta0 = 0, omega = 0) {
  mu <- inv_logit(offset + x * beta)
  s <- sqrt(1 / (sum(mu * (1 - mu) * x ^ 2) + omega))
  betac <- rnorm(1, beta, s)
  muc <- inv_logit(offset + x * betac)
  sc <- sqrt(1 / (sum(muc * (1 - muc) * x ^ 2) + omega))
  log_R <- log_posterior(betac, y, x, offset, beta0, omega) +
    dnorm(beta, betac, sc, log = TRUE) -
    (log_posterior(beta, y, x, offset, beta0, omega) +
     dnorm(betac, beta, s, log = TRUE))
  ifelse(log_R >= 0 || log(runif(1)) < log_R, betac, beta)
}

sample_lpmh <- function (ns, y, x, offset, beta0 = 0, omega = 0,
                         start = (logit(mean(y)) - offset) / mean(x)) {
  beta <- numeric(ns)
  beta[1] <- start
  for (is in 2:ns)
    beta[is] <- step_laplace_mh(beta[is - 1], y, x, offset, beta0, omega)
  beta
}


# [ Example ]
# Urn with `W` white balls, `B` = exp(beta) blue balls:
# sample `n` balls, y_i = I(ball is blue)
W <- 90
n <- 100
yc <- 10 # #observed blue balls out of `n` draws
y <- integer(n); y[sample.int(n, yc)] <- 1
x <- rep(1, n)
offset <- -log(W)

# Laplace approximation around mode:
b <- seq(0, 4, length = 100)
plot(b, sapply(b, log_posterior, y, x, offset), type = 'l')
bmode <- logit(mean(y)) - offset 
lmode <- log_posterior(bmode, y, x, offset)
mu <- inv_logit(offset + x * bmode)
s <- sqrt(1 / (sum(mu * (1 - mu) * x ^ 2)))
lshift <- -dnorm(bmode, bmode, s, log = TRUE) + lmode
lines(b, dnorm(b, bmode, s, log = TRUE) + lshift, lty = 2)


mcmc_array <- function (ns, nchains, params) {
  nparams <- length(params)
  array(dim = c(ns, nchains, nparams),
        dimnames = list(iterations = NULL,
                        chains = paste0("chain", 1:nchains),
                        parameters = params))
}
ns <- 1000
nchains <- 4
sims <- mcmc_array(ns, nchains, "beta")
for (ic in 1:nchains)
  sims[, ic, ] <- sample_lpmh(ns, y, x, offset)

library(bayesplot)
rhat <- function (sims, ...)
  rstan::monitor(sims, print = FALSE, ...)[, "Rhat"]

mcmc_trace(sims)
mcmc_acf(sims)
rhat(sims)

mcmc_dens_overlay(sims)
mcmc_areas(sims)

