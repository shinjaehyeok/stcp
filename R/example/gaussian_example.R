# If EDCP is not installed then run below commands
# install.packages("devtools")
# devtools::install_github("shinjaehyeok/EDCP")

library(EDCP)

# Gaussian case ----
# Pre-change : N(0,1)
# Post-change: N(mu, 1)
# Change-point: v = 0 (immediate change), 200, 500, 1000 (no change)

max_sample <- 1000L
mu <- .5
sig <- 1
v <- 500
alpha <- 1e-4

# Compute optimal delta star
delta_star <- mu
# delta_star <- 2
delta_upper <- 10
delta_lower <- 0.001

psi_fn_list <- generate_sub_G_fn()
v_min <- 1
k_max <- 1e+3


# Generate sample
m_pre <- 0
x_vec <- generator(
  max_sample,
  v,
  pre_sampler = function(n) {
    rnorm(n)
  },
  post_sampler = function(n) {
    rnorm(n, mu, sig)
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
  y = c(0, 0),
  col = 2,
  lwd = 2
)
graphics::lines(
  x = c(v, max_sample),
  y = c(mu, mu),
  col = 2,
  lwd = 2
)

# Build CP detectors
# When delta_lower = delta_upper = delta_star
edcp_star <- build_edcp_exp(
  alpha,
  m_pre,
  delta_star,
  delta_star,
  is_test = FALSE,
  psi_fn_list,
  s_fn = function(x) {
    x - m_pre
  },
  v_fn = function(x) {
    1
  },
  v_min,
  k_max,
  tol = 1e-6
)

# When delta_lower < delta_star < delta_upper
edcp_mix <- build_edcp_exp(
  alpha,
  m_pre,
  delta_lower,
  delta_upper,
  is_test = FALSE,
  psi_fn_list,
  s_fn = function(x) {
    x - m_pre
  },
  v_fn = function(x) {
    1
  },
  v_min,
  k_max,
  tol = 1e-6
)


stopped_time <- run_quick_simulation(x_vec,
                                     edcp_mix,
                                     v,
                                     edcp_star)
print(stopped_time)




