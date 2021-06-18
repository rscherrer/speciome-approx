library(testthat)

source("functions.R")

test_that("Model is loaded correctly", {

  model <- get_model()
  pars <- get_default_pars()
  x <- xres <- 1
  for (i in seq(pars)) eval(pars[[i]])
  N <- c(1000, 1000)
  for (i in seq(model)) eval(model[[i]])
  expect_true(exists("Lambda"))

})

test_that("Demographic equilibrium", {

  init <- c(1000, 1000)
  N <- find_equilibrium(0, get_default_pars(), init)
  expect_true(all(N > 0))

})

test_that("Invasion fitness of the resident is one", {

  init <- c(1000, 1000)
  expect_equal(get_lambda(0, 0, get_default_pars(), init), 1)
  expect_equal(get_lambda(-1, -1, get_default_pars(), init), 1)

})

test_that("Case where a mutant should invade", {

  init <- c(1000, 1000)
  lambda <- get_lambda(0.1, 0, get_default_pars(), init)
  expect_true(lambda > 1)

})

test_that("Pairwise invasibility plot", {

  xvalues <- seq(-1, 1, 1)
  pars <- get_default_pars()
  init <- c(1000, 1000)
  expect_true("ggplot" %in% class(plot_pip(xvalues, pars, init)))
  expect_true("ggplot" %in% class(plot_pip(xvalues, pars, init, binary = FALSE)))
  expect_true("tbl" %in% class(plot_pip(xvalues, pars, init, plotit = FALSE)))

})

test_that("Gradient is zero", {

  G <- get_gradient(0, get_default_pars(), c(1000, 1000))
  expect_equal(G, 0)

})

test_that("Uniroot does not error when ends are on the same side", {

  # A function with two roots between -1 and 1
  f <- function(x) cos(x)^4 - 4*cos(x)^3 + 8*cos(x)^2 - 5*cos(x) + 1/2

  # f(x) is negative at the two ends, regular uniroot should error
  expect_error(uniroot(f, c(-1, 1)))
  expect_true(is.null(uniroot_noerr(f, c(-1, 1))))

})

test_that("Sliding windows", {

  windows <- get_sliding_windows(0, 1, 10, 2)
  expect_equal(windows[[1]], c(0, 0.2))
  expect_equal(windows[[length(windows)]], c(0.8, 1))

})

test_that("Finding multiple roots", {

  # A function with two roots between -1 and 1
  f <- function(x) cos(x)^4 - 4*cos(x)^3 + 8*cos(x)^2 - 5*cos(x) + 1/2

  # Three roots between -1 and 2
  expect_equal(length(find_roots(f, -1, 2)), 3)

})

test_that("Zero should be a singularity", {

  pars <- get_default_pars()
  xeq <- find_singularities(from = -1, to = 1, pars, init = c(1000, 1000))
  expect_equal(xeq, 0)

})

test_that("Zero should be evolutionarily unstable", {

  pars <- get_default_pars()
  expect_false(is_stable(0, pars, init = c(1000, 1000)))

})

test_that("Zero should be convergent stable", {

  pars <- get_default_pars()
  expect_true(is_convergent(0, pars, init = c(1000, 1000)))

})

test_that("Dimorphic model is loaded correctly", {

  model <- get_model_di()
  pars <- get_default_pars()
  x1 <- x2 <- xres1 <- xres2 <- 1
  for (i in seq(pars)) eval(pars[[i]])
  N1 <- N2 <- c(1000, 1000)
  for (i in seq(model)) eval(model[[i]])
  expect_true(exists("Lambda1"))
  expect_true(exists("Lambda2"))

})

test_that("Dimorphic demographic equilibrium", {

  init <- c(1000, 1000, 1000, 1000)
  N <- find_equilibrium_di(0, 0, get_default_pars(), init)
  expect_true(all(N > 0))

})

test_that("Dimorphic gradients should be opposites", {

  pars <- get_default_pars()
  init <- rep(1000, 4)
  G <- get_gradient_di(-1, 1, pars, init)
  expect_equal(G[1], -G[2])

})

test_that("Monomorphic simulation", {

  pars <- get_default_pars()
  init <- rep(1000, 2)
  data <- simulate_mono(-1, 10, pars, init)
  expect_true(is_tibble(data))
  expect_equal(colnames(data), c("time", "x"))

})

test_that("Dimorphic simulation", {

  pars <- get_default_pars()
  init <- rep(1000, 4)
  data <- simulate_di(0, 10, pars, init)
  expect_true(is_tibble(data))
  expect_equal(colnames(data), c("time", "x1", "x2"))

})

test_that("Monomorphic simulation stops when branching point is reached", {

  pars <- get_default_pars()
  init <- rep(1000, 2)
  data1 <- simulate_mono(-0.1, 20, pars, init, branch = FALSE)
  data2 <- simulate_mono(-0.1, 20, pars, init, tol = 0.01)
  expect_true(nrow(data2) < nrow(data1))

})

test_that("Evolution from monomorphic to dimorphic", {

  pars <- get_default_pars()
  init <- rep(1000, 4)
  data <- simulate(-0.1, 50, pars, init)
  expect_equal(nrow(data), 51)
  expect_equal(colnames(data), c("time", "x1", "x2"))
  expect_equal(data$x1[1], data$x2[1]) # should be one trait at the beginning
  expect_false(last(data$x1) == last(data$x2)) # two traits at the end

})

test_that("Burn-in equilibrium density", {

  pars <- get_default_pars()
  N <- find_equilibrium_burnin(0, pars)
  expect_true(N > 0)

})

test_that("Burn-in selection gradient should be negative", {

  pars <- get_default_pars()
  G <- get_gradient_burnin(0, pars)
  expect_true(G < 0)

})

test_that("Simulate the burn-in", {

  pars <- get_default_pars()
  data <- simulate_burnin(0, 10, pars)
  expect_true(is_tibble(data))
  expect_equal(colnames(data), c("time", "x"))

})

test_that("Full simulation with burn-in", {

  pars <- get_default_pars()
  init <- rep(1000, 4)
  data <- simulate(0, 10, pars, init, burnin = 10)
  expect_equal(nrow(data), 21)
  expect_true(data$time[1] < 0)

})
