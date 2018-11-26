# [ MA 589: Gibbs sampler example ]
# Capture-recapture model

library(bayesplot)
rhat <- function (sims, ...)
  rstan::monitor(sims, print = FALSE, ...)[, "Rhat"]
neff_ratio <- function (sims, ...)
  rstan::monitor(sims, print = FALSE, ...)[, "n_eff"]
# format expected by `bayesplot`: iterations, chains, params
mcmc_array <- function (ns, nchains, params) {
  nparams <- length(params)
  array(dim = c(ns, nchains, nparams),
        dimnames = list(iterations = NULL,
                        chains = paste0("chain", 1:nchains),
                        parameters = params))
}

# [ Gibbs sampler ]
log_posterior <- function (N, alpha, captured, r) {
  lgamma(N + 1) - lgamma(N - r + 1) +
    sum(captured * log(alpha) + (N - captured) * log(1 - alpha))
}

sample_chain <- function (ns, captured, r, N0, alpha0) {
  k <- length(captured)
  N <- integer(ns)
  alpha <- matrix(nrow = ns, ncol = k)
  lp <- numeric(ns)
  N[1] <- N0; alpha[1,] <- alpha0
  lp[1] <- log_posterior(N[1], alpha[1], captured, r)
  for (is in 2:ns) {
    N[is] <- r + rnbinom(1, r + 1, 1 - prod(1 - alpha[is - 1,]))
    alpha[is,] <- rbeta(k, captured + 1, N[is] - captured + 1)
    lp[is] <- log_posterior(N[is], alpha[is,], captured, r)
  }
  list(N = N, alpha = alpha, lp = lp)
}
  

# [ Example: Fur seal pup study, example 7.6, pg 212 ]
captured <- c(30, 22, 29, 26, 31, 32, 35)
marked_new <- c(30, 8, 17, 7, 9, 8, 5)
r <- sum(marked_new)

ns <- 1000; nchains <- 4
params <- c("N", paste0("alpha", seq_along(captured)), "lp")
sims <- mcmc_array(ns, nchains, params)
N0 <- r; alpha0 <- captured / N0
for (ic in 1:nchains) {
  ch <- sample_chain(ns, captured, r, N0, alpha0)
  sims[, ic, ] <- cbind(ch$N, ch$alpha, ch$lp)
  #sims[, ic, 1] <- ch$N
  #sims[, ic, 1 + seq_along(captured)] <- ch$alpha
  #sims[, ic, 2 + length(captured)] <- ch$lp
}

color_scheme_set("mix-blue-red")
mcmc_areas(sims, pars = paste0("alpha", seq_along(captured)))
mcmc_violin(sims)
mcmc_dens_overlay(sims)
mcmc_pairs(sims, pars = c("N", "alpha1", "alpha2"),
           diag_fun = "hist", off_diag_fun = "hex")

mcmc_trace(sims)
mcmc_rhat(rhat(sims))
mcmc_neff(neff_ratio(sims))
mcmc_acf(sims)


# [ extra: MAP estimates ]
N <- N0; alpha <- alpha0; tol <- 1e-6
k <- length(captured)
g <- numeric(k + 1); H <- matrix(0, nrow = k + 1, ncol = k + 1)
l <- log_posterior(N, alpha, captured, r)
repeat {
  g[1] <- digamma(N + 1) - digamma(N - r + 1) + sum(log(1 - alpha))
  g[-1] <- captured / alpha - (N - captured) / (1 - alpha)
  H[1, 1] <- trigamma(N + 1) - trigamma(N - r + 1)
  H[1, -1] <- H[-1, 1] <- -1 / (1 - alpha)
  diag(H[-1, -1]) <- -captured / alpha ^ 2 - (N - captured) / (1 - alpha) ^ 2
  C <- chol(-H)
  delta <- backsolve(C, backsolve(C, g, transpose = TRUE))
  N <- N + delta[1]; alpha <- alpha + delta[-1]
  l.new <- log_posterior(N, alpha, captured, r)
  if (abs((l.new - l) / l) < tol) break
  l <- l.new
}

