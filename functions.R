
# Set up the model, not evaluated
get_model <- function() {

  alist(

    w1 <- w0 * exp(-s * (x + psi)^2 / psi^2),
    w2 <- w0 * exp(-s * (x - psi)^2 / psi^2),
    w1res <- w0 * exp(-s * (xres + psi)^2 / psi^2),
    w2res <- w0 * exp(-s * (xres - psi)^2 / psi^2),

    R11 <- iota / (omicron + N[1] * w1res),
    R12 <- h * iota / (omicron + N[2] * w1res),
    R21 <- h * iota / (omicron + N[1] * w2res),
    R22 <- iota / (omicron + N[2] * w2res),

    W1 <- w1 * R11 + w2  * R21,
    W2 <- w1 * R12 + w2  * R22,

    A <- exp(-a * (x - xres)^2),

    r1 <- 1 - d + 0.5 * W1 * (1 + A),
    r2 <- 1 - d + 0.5 * W2 * (1 + A),

    M <- matrix(c(1 - m, m, m, 1 - m), 2, 2, byrow = TRUE),
    Q <- matrix(c(r1, 0, 0, r2), 2, 2, byrow = TRUE),
    Lambda <- M %*% Q

  )

}

# Generate default parameter values
get_default_pars <- function() {

  alist(
    m <- 0.1,
    d <- 0.2,
    s <- 1,
    psi <- 1,
    w0 <- 1,
    iota <- 400,
    omicron <- 100,
    a <- 1,
    h <- 1
  )

}

# Function to find the demographic equilibrium
find_equilibrium <- function(xres, pars, init) {

  # Of the resident population
  x <- xres

  # Unpack the parameters
  for (i in seq(pars)) eval(pars[[i]])

  # Model setup
  model <- get_model()

  # System of equations to solve
  f <- function(N) {

    for (i in seq(model)) eval(model[[i]])
    return(N - Lambda %*% N)

  }

  # Solve the demographic equilibrium
  N <- fsolve(f, init)$x

  return(N)

}

# Function to compute the invasion fitness
get_lambda <- function(x, xres, pars, init, fast = FALSE) {

  # Unpack the parameters
  for (i in seq(pars)) eval(pars[[i]])

  # Find the demographic equilibrium
  if (fast) N <- init else N <- find_equilibrium(xres, pars, init)

  # Evaluate the model at equilibrium
  model <- get_model()
  for (i in seq(model)) eval(model[[i]])

  # Compute invasion fitness
  0.5 * ((1 - m) * (r1 + r2) + sqrt((1 - m)^2 * (r1^2 - 2 * r1 * r2 + r2^2) + 4 * m^2 * r1 * r2))

}

# Function to plot a pairwise invasibility plot
plot_pip <- function(xvalues, pars, init, binary = TRUE, plotit = TRUE) {

  # Compute demographic equilibria of the different residents
  N <- map(xvalues, find_equilibrium, pars, init)

  # Create a matrix of trait values and population sizes
  data <- matrix(unlist(N), ncol = 2, byrow = TRUE)
  colnames(data) <- c("N1", "N2")
  data <- as_tibble(data)
  data <- data %>% add_column(xres = xvalues, .before = "N1")

  # Expand to all combinations of mutants and residents
  data <- data %>% expand(x = xvalues, nesting(xres, N1, N2))

  # Compute invasion fitness
  data <- data %>%
    mutate(
      lambda = pmap_dbl(
        list(x, xres, N1, N2),
        ~ get_lambda(..1, ..2, pars, c(..3, ..4), fast = TRUE)
      ),
      invades = lambda > 1
    )

  # Early exit
  if (!plotit) return(data)

  # Pairwise invasibility plot
  if (binary) {

    # Invades or not
    plot <- data %>%
      ggplot(aes(x = xres, y = x, fill = invades)) +
      geom_tile() +
      xlab(parse(text = "'Resident trait value'~hat(x)")) +
      ylab(parse(text = "'Mutant trait value'~x")) +
      labs(fill = parse(text = "lambda(x,hat(x))>1")) +
      scale_fill_manual(values = c("gray20", "gray80"))

  } else {

    # Actual value of the invasion fitness
    plot <- data %>%
      ggplot(aes(x = xres, y = x, fill = lambda)) +
      geom_tile() +
      xlab(parse(text = "'Resident trait value'~hat(x)")) +
      ylab(parse(text = "'Mutant trait value'~x")) +
      labs(fill = parse(text = "lambda(x,hat(x))")) +
      scale_fill_continuous(type = "viridis", option = "magma")

  }

  return(plot)

}

# Function to compute the selection gradient
get_gradient <- function(xres, pars, init) {

  # Evaluation is at the resident
  x <- xres

  # Unpack the parameters
  for (i in seq(pars)) eval(pars[[i]])

  # Find the demographic equilibrium
  N <- find_equilibrium(xres, pars, init = init)

  # Evaluate the model at equilibrium
  model <- get_model()
  for (i in seq(model)) eval(model[[i]])

  # Derivatives of the attack rates
  dw1 <- -2 * s * w0 / psi * (x + psi) * w1
  dw2 <- -2 * s * w0 / psi * (x - psi) * w2

  # Derivatives of the reproductive success
  dW1 <- R11 * dw1 + R21 * dw2
  dW2 <- R12 * dw1 + R22 * dw2

  # Return the selection gradient
  denom <- (1 - r2 + m * r2)^2 + m^2 * r1 * r2
  num <- (1 - m - (2 - 4 * m + m^2) * r2 + (1 - 3 * m + 2 * m^2) * r2^2) * dW1 + m^2 * r1 * dW2
  num / denom

}

# Uniroot without errors
uniroot_noerr <- function(f, interval, ...) {

  tryCatch(

    uniroot(f, interval, ...),
    error = function(err) return(NULL)

  )

}

# Function to get starts and ends of sliding windows
get_sliding_windows <- function(from, to, n, size, intervals = TRUE) {

  # Sequence of breaks
  breaks <- seq(from, to, length.out = n + 1)

  # Pick the starts and ends of each window
  starts <- breaks[1:(length(breaks) - size)]
  ends <- breaks[(1 + size):length(breaks)]

  # Return a list of intervals
  if (intervals) return(map2(starts, ends, ~ c(.x, .y)))

  # Or return start and end points separated
  return(list(starts = starts, ends = ends))

}

# Function to find the roots of a function with one argument
find_roots <- function(f, from, to, n = 10, size = 1, precis = NULL, ...) {

  if (size > n) stop("size must be smaller or equal to n")

  # Sliding windows
  windows <- get_sliding_windows(from, to, n, size)

  # Look for roots in multiple windows
  roots <- map(windows, function(win) uniroot_noerr(f, interval = win, ...))

  # Remove nulls
  roots <- roots[!map_lgl(roots, is.null)]

  # Simplify
  roots <- map_dbl(roots, ~ .x[["root"]])

  # Remove duplicates (at a certain level of precision)
  if (!is.null(precis)) roots <- round(roots, precis)
  roots[!duplicated(roots)]

}

# Function to find the roots of the selection gradient
find_singularities <- function(from, to, pars, init, ...) {

  # Equation for which to find the root(s)
  f <- function(x) get_gradient(x, pars, init)

  # Find the roots of the selection gradient
  find_roots(f, from, to, ...)

}

# Evolutionary stability criterion
is_stable <- function(xeq, pars, init) {

  # Evaluation is at equilibrium
  xres <- xeq
  x <- xres

  # Unpack the parameters
  for (i in seq(pars)) eval(pars[[i]])

  # Find the demographic equilibrium
  N <- find_equilibrium(xeq, pars, init)

  # Evaluate the model at equilibrium
  model <- get_model()
  for (i in seq(model)) eval(model[[i]])

  # Derivatives of the attack rates
  dw1 <- -2 * s * w0 / psi * (x + psi) * w1
  dw2 <- -2 * s * w0 / psi * (x - psi) * w2

  # Derivatives of the reproductive success
  dW1 <- R11 * dw1 + R21 * dw2
  dW2 <- R12 * dw1 + R22 * dw2

  # Second derivatives of the attack rates
  ddw1 <- -2 * s * w0 / psi * (w1 + (x + psi) * dw1)
  ddw2 <- -2 * s * w0 / psi * (w2 + (x - psi) * dw2)

  # Second derivatives of the reproductive success
  ddW1 <- R11 * ddw1 + R21 * ddw2
  ddW2 <- R12 * ddw1 + R22 * ddw2

  # Compute the second derivative of the fitness function
  denom <- 1 + m^2 * r1 * r2 + (m - 1) * (2 + (m - 1) * r2^2)
  num <- m^2 * r1 * ddW2 + 2 * (2 * m - 1) * (1 + (m - 1) * r2) * dW1 * dW2 + (1 + (m - 1) * r2) * (1 - m + (2 * m - 1) * r2) * ddW1 - a * (m^2 * r1 * W2 - (m - 1 + (2 - 4 * m + m^2) * r2 - (1 - 3 * m + 2 * m^2) * r2^2) * W1)

  # Is the equilibrium stable?
  num / denom < 0

}

# Convergence stability criterion
is_convergent <- function(xeq, pars, init, step = 0.0001) {

  # Compute the selection gradient around the equilibrium
  gdt_before <- get_gradient(xeq - step, pars, init)
  gdt_after <- get_gradient(xeq + step, pars, init)

  # Compute the difference
  delta <- gdt_after - gdt_before

  # Is the equilibrium convergence-stable?
  delta < 0

}

# Dimorphic model (after branching)
get_model_di <- function() {

  alist(

    w11 <- w0 * exp(-s * (x1 + psi)^2 / psi^2),
    w21 <- w0 * exp(-s * (x1 - psi)^2 / psi^2),
    w12 <- w0 * exp(-s * (x2 + psi)^2 / psi^2),
    w22 <- w0 * exp(-s * (x2 - psi)^2 / psi^2),
    w1res1 <- w0 * exp(-s * (xres1 + psi)^2 / psi^2),
    w2res1 <- w0 * exp(-s * (xres1 - psi)^2 / psi^2),
    w1res2 <- w0 * exp(-s * (xres2 + psi)^2 / psi^2),
    w2res2 <- w0 * exp(-s * (xres2 - psi)^2 / psi^2),

    R11 <- iota / (omicron + N1[1] * w1res1 + N2[1] * w1res2),
    R12 <- h * iota / (omicron + N1[2] * w1res1 + N2[2] * w1res2),
    R21 <- h * iota / (omicron + N1[1] * w2res1 + N2[1] * w2res2),
    R22 <- iota / (omicron + N1[2] * w2res1 + N2[2] * w2res2),

    W11 <- w11 * R11 + w21 * R21,
    W21 <- w11 * R12 + w21 * R22,
    W12 <- w12 * R11 + w22 * R21,
    W22 <- w12 * R12 + w22 * R22,

    A1 <- exp(-a * (x1 - xres1)^2),
    A2 <- exp(-a * (x2 - xres2)^2),

    r11 <- 1 - d + 0.5 * W11 * (1 + A1),
    r21 <- 1 - d + 0.5 * W21 * (1 + A1),
    r12 <- 1 - d + 0.5 * W12 * (1 + A2),
    r22 <- 1 - d + 0.5 * W22 * (1 + A2),

    M <- matrix(c(1 - m, m, m, 1 - m), 2, 2, byrow = TRUE),
    Q1 <- matrix(c(r11, 0, 0, r21), 2, 2, byrow = TRUE),
    Q2 <- matrix(c(r12, 0, 0, r22), 2, 2, byrow = TRUE),

    Lambda1 <- M %*% Q1,
    Lambda2 <- M %*% Q2

  )

}

# Function to find the dimorphic demographic equilibrium
find_equilibrium_di <- function(xres1, xres2, pars, init) {

  # Of the resident population
  x1 <- xres1
  x2 <- xres2

  # Unpack the parameters
  for (i in seq(pars)) eval(pars[[i]])

  # Model setup
  model <- get_model_di()

  # System of equations to solve
  f <- function(N) {

    # Decompose the vector of population sizes into the two species
    N1 <- N[1:2]
    N2 <- N[3:4]

    # Evaluate the model
    for (i in seq(model)) eval(model[[i]])

    # Assemble the two transition matrices into one
    Zeros <- matrix(0, 2, 2)
    Lambda <- rbind(cbind(Lambda1, Zeros), cbind(Zeros, Lambda2))

    return(N - Lambda %*% N)

  }

  # Solve the demographic equilibrium
  N <- fsolve(f, init)$x

  return(N)

}

# Function to compute the dimorphic selection gradient
get_gradient_di <- function(xres1, xres2, pars, init) {

  # Evaluation is at the resident
  x1 <- xres1
  x2 <- xres2

  # Unpack the parameters
  for (i in seq(pars)) eval(pars[[i]])

  # Find the demographic equilibrium
  N <- find_equilibrium_di(xres1, xres2, pars, init = init)

  # Separate population sizes for the two species
  N1 <- N[1:2]
  N2 <- N[3:4]

  # Evaluate the model at equilibrium
  model <- get_model_di()
  for (i in seq(model)) eval(model[[i]])

  # List of species-level values
  l <- list(
    x = c(x1, x2), # trait value
    w1 = c(w11, w12), # attack rate on resource 1
    w2 = c(w21, w22), # attack rate on resource 2
    r1 = c(r11, r12), # growth rate in habitat 1
    r2 = c(r21, r22) # growth rate in habitat 2
  )

  # For each species...
  G <- pmap_dbl(l, function(x, w1, w2, r1, r2) {

    # Derivatives of the attack rates
    dw1 <- -2 * s * w0 / psi * (x + psi) * w1
    dw2 <- -2 * s * w0 / psi * (x - psi) * w2

    # Derivatives of the reproductive success
    dW1 <- R11 * dw1 + R21 * dw2
    dW2 <- R12 * dw1 + R22 * dw2

    # Return the selection gradient
    denom <- (1 - r2 + m * r2)^2 + m^2 * r1 * r2
    num <- (1 - m - (2 - 4 * m + m^2) * r2 + (1 - 3 * m + 2 * m^2) * r2^2) * dW1 + m^2 * r1 * dW2
    num / denom

  })

  return(G)

}

# Function to simulate evolution
simulate_mono <- function(
  xstart, ntimes, pars, init, mu = 0.01, sigma = 0.1, branch = TRUE,
  tol = 0.0001
) {

  # Initialize
  x <- xstart
  xvalues <- rep(NA, ntimes + 1)
  xvalues[1] <- x

  # At each time point...
  for (t in 1:ntimes) {

    # Compute the selection gradient and demographic equilibrium
    G <- get_gradient(x, pars, init)
    N <- find_equilibrium(x, pars, init)

    # Break if we have reached a branching point, if needed
    if (branch & G < tol & G > -tol & !is_stable(x, pars, init)) break

    # Evolve
    kappa <- 0.5 * sum(N) * sigma^2 * mu
    dx <- kappa * G
    x <- x + dx

    # Record
    xvalues[t + 1] <- x

  }

  xvalues <- xvalues[!is.na(xvalues)]

  # Assemble into a table
  return(tibble(time = seq(xvalues) - 1, x = xvalues))

}

# Function to simulate evolution after branching
simulate_di <- function(xstart, ntimes, pars, init, mu = 0.01, sigma = 0.1, dodge = 0.001) {

  # Dodge the starting value by a little amount to start two lineages
  x1 <- xstart - dodge
  x2 <- xstart + dodge

  # Initialize
  x1values <- x2values <- rep(NA, ntimes + 1)
  x1values[1] <- x1
  x2values[1] <- x2

  # At each time point...
  for (t in 1:ntimes) {

    # Compute the selection gradient and demographic equilibrium
    G <- get_gradient_di(x1, x2, pars, init)
    N <- find_equilibrium_di(x1, x2, pars, init)

    # Evolve the first species
    kappa1 <- 0.5 * sum(N[1:2]) * sigma^2 * mu
    dx1 <- kappa1 * G[1]
    x1 <- x1 + dx1

    # Evolve the second species
    kappa2 <- 0.5 * sum(N[3:4]) * sigma^2 * mu
    dx2 <- kappa2 * G[2]
    x2 <- x2 + dx2

    # Record
    x1values[t + 1] <- x1
    x2values[t + 1] <- x2

  }

  # Assemble into a table
  return(tibble(time = 0:ntimes, x1 = x1values, x2 = x2values))

}

# Function to find the (non-extinct) burn-in demographic equilibrium
find_equilibrium_burnin <- function(xres, pars) {

  # Unpack the parameters
  for (i in seq(pars)) eval(pars[[i]])

  # Attack rate of the resident
  w <- w0 * exp(-s * (xres + psi)^2 / psi^2)

  # Compute equilibrium population size
  num <- w * iota - d * omicron
  denom <- d * w
  num / denom

}

# Function to compute the burn-in selection gradient
get_gradient_burnin <- function(xres, pars) {

  # Evaluation is at the resident
  x <- xres

  # Unpack the parameters
  for (i in seq(pars)) eval(pars[[i]])

  # Find the demographic equilibrium
  N <- find_equilibrium_burnin(xres, pars)

  # Attack rate of the resident
  w <- w0 * exp(-s * (xres + psi)^2 / psi^2)

  # Equilibrium resource concentration
  R <- iota / (omicron + N * w)

  # Selection gradient
  G <- -2 * s * w0 * (xres + psi) * w * R

  return(G)

}

# Function to simulate the burn-in
simulate_burnin <- function(xstart, ntimes, pars, mu = 0.01, sigma = 0.1) {

  # Initialize
  x <- xstart
  xvalues <- rep(NA, ntimes + 1)
  xvalues[1] <- x

  # At each time point...
  for (t in 1:ntimes) {

    # Compute the selection gradient and demographic equilibrium
    G <- get_gradient_burnin(x, pars)
    N <- find_equilibrium_burnin(x, pars)

    # Evolve
    kappa <- 0.5 * N * sigma^2 * mu
    dx <- kappa * G
    x <- x + dx

    # Record
    xvalues[t + 1] <- x

  }

  # Assemble into a table
  return(tibble(time = seq(xvalues) - 1, x = xvalues))

}

# Function to simulate evolution
simulate <- function(
  xstart, ntimes, pars, init, mu = 0.01, sigma = 0.1, tol = 0.0001,
  dodge = 0.001, burnin = 0
) {

  # If the burn-in must be simulated...
  if (burnin > 0) {

    # Simulate the burn-in
    data_burnin <- simulate_burnin(xstart, burnin, pars, mu, sigma)

    # Give burn-in time points negative values
    data_burnin <- data_burnin %>% mutate(time = time - n() + 1)

    # Record the end point of the burn-in
    xstart <- last(data_burnin$x)

    # Remove the end point of the burn-in (it is the start of the next phase)
    data_burnin <- data_burnin[-nrow(data_burnin),]

    # Reformat the burn-in so it can be appended later
    data_burnin <- data_burnin %>% mutate(x2 = x) %>% rename(x1 = "x")

  }

  # Monomorphic evolution until potential branching point
  data <- simulate_mono(xstart, ntimes, pars, init[1:2], mu, sigma, tol = tol)

  # Reformat the dataset
  data <- data %>% mutate(x2 = x) %>% rename(x1 = "x")

  # Number of time steps left
  ntimes <- ntimes - nrow(data) + 1

  # If there is simulation time left...
  if (ntimes > 0) {

    # The endpoint of monomorphic evolution is the new starting point
    xstart <- last(data$x1)

    # Simulate dimorphic evolution after the branching point
    data_di <- simulate_di(xstart, ntimes, pars, init, mu, sigma, dodge)

    # The first time step of dimorphic evolution is the last one from before
    data_di <- data_di[-1,]

    # Assemble data before and after the branching point
    data <- map_dfr(list(data, data_di), ~ .x) %>% mutate(time = seq(n()) - 1)

  }

  # Add the burn-in data if needed
  if (burnin > 0) data <- map_dfr(list(data_burnin, data), ~ .x)

  return(data)

}

# Example parameter set
get_example_pars <- function() {

  alist(
    m <- 0.01,
    d <- 0.2,
    s <- 1.4,
    psi <- 5,
    w0 <- 1,
    iota <- 400,
    omicron <- 100,
    a <- 0,
    h <- 1
  )

}
