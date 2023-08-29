# If stcp is not installed then run below commands
# install.packages("devtools")
# devtools::install_github("shinjaehyeok/stcp")

library(stcp)

# Bernoulli case ----
# Null : Ber(0.5)
# Alternative: Ber(0.52)

max_sample <- 1e+4L
m_pre <- 0.5
m_post <- 0.52
alpha <- 1e-3


# Compute optimal delta star
delta_star <- m_post - m_pre
delta_upper <- 0.2
delta_lower <- 0.005

psi_fn_list <- generate_sub_B_fn(p = m_pre)
v_min <- 1
k_max <- 1e+3


# Generate sample
x_vec <- rbinom(max_sample, 1, m_pre)

# Build CP detectors
# When delta_lower = delta_upper = delta_star
stcp_star <- build_stcp_exp(
  alpha,
  m_pre,
  delta_star,
  delta_star,
  is_test = TRUE,
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
  is_test = TRUE,
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

run_star <- run_stcp(x_vec, stcp_star)
run_mix <- run_stcp(x_vec, stcp_mix)

plot(run_star)
lines(1:max_sample, run_mix$log_mix_e_vec, col = 2)
print(run_star)
print(run_mix)
