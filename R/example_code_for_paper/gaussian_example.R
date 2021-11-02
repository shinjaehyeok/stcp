# If EDCP is not installed then run below commands
# install.packages("devtools")
# devtools::install_github("shinjaehyeok/EDCP")

library(EDCP)

# Gaussian case ----
# Pre-change : N(0,1)
# Post-change: N(mu, sig)
# Change-point: v = 0 (immediate chane), 200, 500, 1000 (no change)

max_sample <- 1000L
mu <- 1
sig <- sqrt(2)
v <- 500
alpha <- 0.1/max_sample

# Compute optimal delta star
delta_star <- mu / sig^2
# delta_star <- 2
delta_upper <- 10
delta_lower <- 0.001

psi_fn_list <- generate_sub_G_fn()
v_min <- 1
k_max <- 1e+3


# Generate sample
x_vec <- generator(
  max_sample,
  v,
  pre_sampler = function(n){rnorm(n)},
  post_sampler = function(n){rnorm(n, mu, sig)})

plot(1:max_sample, x_vec, pch=20,
     xlab = "n", ylab = "X_n", main = "Simulated Data")
graphics::lines(x = c(0,v), y = c(0,0), col = 2, lwd = 2)
graphics::lines(x = c(v,max_sample), y = c(mu,mu), col = 2, lwd =2 )

stopped_time <- run_quick_simulation(
  x_vec,
  v,
  alpha,
  delta_star,
  delta_upper,
  delta_lower,
  psi_fn_list = psi_fn_list,
  v_min,
  k_max
)
print(stopped_time)
