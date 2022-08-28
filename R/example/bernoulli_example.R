# If stcp is not installed then run below commands
# install.packages("devtools")
# devtools::install_github("shinjaehyeok/stcp")

library(stcp)

# Bernoulli case ----
# Pre-change : Ber(0.2)
# Post-change: Ber(0.4)
# Change-point: v = 0 (immediate change), 200, 500, 1000 (no change)

max_sample <- 1000L
m_pre <- 0.2
m_post <- 0.4
v <- 500
alpha <- 1e-4

# Compute optimal delta star
delta_star <- m_post - m_pre
delta_upper <- 0.4
delta_lower <- 0.05

psi_fn_list <- generate_sub_B_fn(p = m_pre)
v_min <- 1
k_max <- 1e+3


# Generate sample
x_vec <- generator(
  max_sample,
  v,
  pre_sampler = function(n) {
    rbinom(n, 1, m_pre)
  },
  post_sampler = function(n) {
    rbinom(n, 1, m_post)
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
stcp_star <- build_stcp_exp(
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
stcp_mix <- build_stcp_exp(
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
                                     stcp_mix,
                                     v,
                                     stcp_star)
print(stopped_time)
