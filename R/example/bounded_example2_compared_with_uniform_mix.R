# If STCP is not installed then run below commands
# install.packages("devtools")
# devtools::install_github("shinjaehyeok/STCP")

library(STCP)

# Bounded case ----

max_sample <- 5000L
sig_scale <- 1
alpha_pre <- 1 / sig_scale
beta_pre <- 4 / sig_scale
alpha_post <- 1 / sig_scale
beta_post <- 3 / sig_scale
m_pre <- alpha_pre / (alpha_pre + beta_pre)
m_post <- alpha_post / (alpha_post + beta_post)
delta_lower <- 0.05
bound_lower <- 0
bound_upper <- 1
v <- max_sample / 2
alpha <- 1e-4
k_max <- 1e+3


# Generate sample
x_vec <- generator(
  max_sample,
  v,
  pre_sampler = function(n) {
    rbeta(n, alpha_pre, beta_pre)
  },
  post_sampler = function(n) {
    rbeta(n, alpha_post, beta_post)
  }
)

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
# Uniform mixture
# Compute parameters
by_gap <- 0.02
lambda <- seq(by_gap, 1, by = by_gap)
omega <- rep(1/length(lambda), length(lambda))

stcp_uniform <- build_stcp(alpha,
                           m_pre,
                           is_test = FALSE,
                           omega = omega,
                           lambda = lambda,
                           log_base_fn_generator = generate_log_bounded_base_fn,
                           m = m_pre,
                           bound_lower = bound_lower)

# When delta_lower < delta_star < delta_upper
stcp_mix <- build_stcp_bounded(
  alpha,
  m_pre,
  delta_lower,
  is_test = FALSE,
  k_max = 1000,
  tol = 1e-6
)

stopped_time <- run_quick_simulation(x_vec,
                                     stcp_mix,
                                     v,
                                     stcp_uniform)
print(stopped_time)


out <- run_stcp(x_vec, stcp_uniform)
print(out$stcp_obj$lambda[which.max(out$stcp_obj$log_e_vec)])

out2 <- run_stcp(x_vec, stcp_mix)
print(out2$stcp_obj$lambda[which.max(out2$stcp_obj$log_e_vec)])
