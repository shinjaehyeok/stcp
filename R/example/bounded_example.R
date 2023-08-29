# If stcp is not installed then run below commands
# install.packages("devtools")
# devtools::install_github("shinjaehyeok/stcp")

library(stcp)

# Bounded case ----
# Pre-change : 5 * Beta(1, 4) - 1 (m = 0, sig^2 = 5^2 * 4 / (5^2 * 6) = 2/3)
# Post-change: 5 * Beta(2, 3) - 1 (m = 1, sig^2 = 5^2 * 6 / (5^2 * 6) = 1)
# Change-point: v = 0 (immediate change), 200, 500, 1000 (no change)

max_sample <- 1000L
m_pre <- 0
m_post <- 1
bound_lower <- -1
bound_upper <- 5
v <- 500
alpha <- 1e-4

# Compute optimal delta star
delta_star <- m_post - m_pre
var_star <- 1
delta_upper <- bound_upper - m_pre
delta_lower <- 0.1

k_max <- 1e+3


# Generate sample
x_vec <- generator(
  max_sample,
  v,
  pre_sampler = function(n) {
    5 * rbeta(n, 1, 4) - 1
  },
  post_sampler = function(n) {
    5 * rbeta(n, 2, 3) - 1
  }
)

plot(
  1:max_sample,
  x_vec,
  pch = 20,
  xlab = "n",
  ylab = "X_n",
  main = "Simulated Data"
)
graphics::lines(
  x = c(0, v),
  y = c(m_pre, m_pre),
  col = 2,
  lwd = 2
)
graphics::lines(
  x = c(v, max_sample),
  y = c(m_post, m_post),
  col = 2,
  lwd = 2
)

# Build CP detectors
# When delta_lower = delta_upper = delta_star
stcp_star <- build_stcp_bounded(
  alpha,
  m_pre,
  delta_star,
  delta_star,
  is_test = FALSE,
  is_flipped = FALSE,
  bound_lower,
  bound_upper,
  k_max = 1000,
  tol = 1e-6,
  var_lower = var_star,
  var_upper = var_star
)

# When delta_lower < delta_star < delta_upper
stcp_mix <- build_stcp_bounded(
  alpha,
  m_pre,
  delta_lower,
  delta_upper,
  is_test = FALSE,
  is_flipped = FALSE,
  bound_lower,
  bound_upper,
  k_max = 1000,
  tol = 1e-6
)

stopped_time <- run_quick_simulation(x_vec,
                                     stcp_mix,
                                     v,
                                     stcp_star)
print(stopped_time)
