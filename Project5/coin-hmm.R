# [ MA 589: Hidden Markov models ]

# soft max (log-sum-exp) of vector `x`
LOGEPS <- log(.Machine$double.eps / 2)
lse <- function (x) {
  m <- max(x); x <- x - m
  m + log(sum(exp(x[x > LOGEPS])))
}


# [ HMM functions ]
# compute 'forward' probabilities from sequence `s`, with initial
# probabilities `logI` (vector of size `ns` = |S|, S is the state space),
# emission probs `logE` (S by emissions matrix) and
# transition probs `logP` (S to S matrix), all in log scale
forward <- function (s, logI, logP, logE) {
  n <- length(s); ns <- length(logI) # = nrow(logE) = nrow(logP)
  f <- matrix(nrow = n, ncol = ns)
  f[1,] <- logI + logE[, s[1]]
  for (i in 2:n) # for each position in sequence
    for (j in 1:ns) # for each X_i
      f[i, j] <- lse(f[i - 1,] + logP[,j]) + logE[j, s[i]]
  f
}

# compute most likely assignment of states given data (MAP)
viterbi <- function (s, logI, logP, logE) {
  n <- length(s); ns <- length(logI) # = nrow(logE) = nrow(logP)
  m <- matrix(nrow = n, ncol = ns) # maxima
  b <- matrix(nrow = n, ncol = ns) # backtrack pointers
  # recurse
  m[1,] <- logI + logE[, s[1]]
  for (i in 2:n) { # for each position in sequence
    for (j in 1:ns) { # for each X_i
      u <- m[i - 1,] + logP[,j]
      m[i, j] <- max(u) + logE[j, s[i]]
      b[i - 1, j] <- which.max(u)
    }
  }
  # backtrack
  v <- numeric(n)
  v[n] <- which.max(m[n,]) # b[n]
  for (i in (n - 1):1)
    v[i] <- b[i, v[i + 1]]
  list(m = m, b = b, seq = v)
}

# Example:
# We have two coins: a "fair" coin with probability of heads .5, and a
# "loaded" coin, with prob. of heads pe=.9. Given that we're flipping the fair
# coin, there's a probability pf=.9 of staying at the fair coin for the next
# flip; similarly, if we're at the loaded coin, there's a pl=.9 probability of
# flipping the loaded coin in the next flip.
# We start with the fair coin and keep transitioning coin "states" (fair or
# loaded) and emitting heads (outcome 1) or tails (outcome 2) to generate a
# sequence of n flips. Our task is to compute the most likely assignment of
# states given the recorded outcomes, which we can do using Viterbi's
# algorithm.

# Parameters
pf <- .9 # self-transition prob for "fair" state
pl <- .9 # self-transition prob for "loaded" state
pe <- .9 # (heads) emission for "loaded" state

ns <- 2 # #states
P <- matrix(c(pf, 1 - pf, 1 - pl, pl), nrow = ns, byrow = T)
E <- matrix(c(.5, .5, pe, 1 - pe), nrow = ns, byrow = T)
logI <- log(c(1, 0)); logP <- log(P); logE <- log(E) # cache logs

sample_coin <- function (n) {
  # generate sequence (Markov chain transitions)
  t <- integer(n)
  t[1] <- 1 # always start at fair coin
  for (i in 2:n)
    t[i] <- sample.int(ns, 1, prob = P[t[i - 1],])
  # generate data (emissions)
  s <- integer(n)
  for (i in 1:n)
    s[i] <- sample.int(ns, 1, prob = E[t[i],])
  list(data = s, seq = t)
}

# test:
n <- 200
s <- sample_coin(n)
f <- forward(s$data, logI, logP, logE)
v <- viterbi(s$data, logI, logP, logE)
