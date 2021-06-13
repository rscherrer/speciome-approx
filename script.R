rm(list = ls())

library(tidyverse)

get_lambda <- function(x, xhat) {

  e1 <- exp(-s * (x + psi)^2 / psi^2)
  e2 <- exp(-s * (x - psi)^2 / psi^2)
  e1hat <- exp(-s * (xhat + psi)^2 / psi^2)
  e2hat <- exp(-s * (xhat - psi)^2 / psi^2)
  h11 <- h22 <- 1
  h12 <- h21 <- h

  N1 == f(N1, N2)
  N2 == g(N1, N2)

  R11 <- iota * h11 / (omicron + N1 * e1hat)
  R12 <- iota * h12 / (omicron + N2 * e1hat)
  R21 <- iota * h21 / (omicron + N1 * e2hat)
  R22 <- iota * h22 / (omicron + N2 * e2hat)
  W1 <- e1 * R11 + e2 * R21
  W2 <- e1 * R12 + e2 * R22
  r1 <- 1 - d + 0.5 * W1 * (1 + A)
  r2 <- 1 - d + 0.5 * W2 * (1 + A)
  lambda <- 0.5 * ((1 - m) * (r1 + r2) + sqrt((1 - m)^2 * (r1^2 - 2 * r1 * r2 + r2^2) + 4 * m^2 * r1 * r2))

}

FUN <- function() {

  pip <- expand_grid(xhat = seq(-5, 5, 1), x = seq(-5, 5, 1))

  pip <- pip %>% mutate(lambda = get_lambda(x, xhat))

  return(pip)

}

FUN()
