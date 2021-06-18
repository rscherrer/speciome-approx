rm(list = ls())

library(pracma)
library(tidyverse)

theme_set(theme_classic())

source("functions.R")
source("tests.R")

#### Analysis ####

# Parameter values
pars <- get_example_pars()

# Approximate rates of evolution
mu <- 0.01 # mutation rate per locus * number of ecological loci
sigma <- 0.01 # 1 / number of ecological loci (variance in normalized effect sizes)

tic <- Sys.time()

# Simulate evolution
data <- simulate(0, ntimes = 100000, pars, init = rep(1000, 4), burnin = 20000, mu = mu, sigma = sigma)

toc <- Sys.time()

toc - tic

# Plot
data %>%
  pivot_longer(c(x1, x2), names_to = "species", values_to = "x") %>%
  ggplot(aes(x = time, y = x, group = species)) +
  geom_line() +
  xlab("Time (generations)") +
  ylab("Trait value")


# Initial conditions for solving the demographic equilibrium
init <- c(1000, 1000)

# Possible resident values
xvalues <- seq(-1, 1, 0.01)

# Pairwise invasibility plot (might take a while)
pip <- plot_pip(xvalues, pars, init)
pip2 <- plot_pip(xvalues, pars, init, binary = FALSE)

# Plot the selection gradient across trait values
tibble(
  x = xvalues,
  G = map_dbl(xvalues, get_gradient, pars, init)
) %>%
  ggplot(aes(x = x, y = G)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab(parse(text = "'Resident trait value'~hat(x)")) +
  ylab(parse(text = "'Selection gradient'~G(hat(x))"))

# Find evolutionary equilibria
xeq <- find_singularities(from = -1, to = 1, pars, init)

# Is the equilibrium stable?
is_stable(xeq, pars, init)

# Is the equilibrium attainable?
is_convergent(xeq, pars, init)

# Simulate evolution
data <- simulate(0, 300, pars, init, burnin = 100, mu = mu, sigma = mu)

data %>%
  pivot_longer(c(x1, x2), names_to = "species", values_to = "x") %>%
  ggplot(aes(x = time, y = x, group = species)) +
  geom_line() +
  xlab("Time (generations)") +
  ylab("Trait value")




