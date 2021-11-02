# If EDCP is not installed then run below commands
# install.packages("devtools")
# devtools::install_github("shinjaehyeok/EDCP")

library(EDCP)

# Bounded case ----
# Pre-change : Beta(0.2)
# Post-change: B(0.4)
# Change-point: v = 0 (immediate chane), 200, 500, 1000 (no change)

max_sample <- 1000L
p_pre <- 0.2
p_post <- 0.4
v <- 500
alpha <- 0.1/max_sample

# Compute optimal delta star
delta_star <- p_post - p_pre
# delta_star <- 2
delta_upper <- 0.4
delta_lower <- 0.05

psi_fn_list <- generate_sub_B_fn(p = p_pre)
v_min <- 1
k_max <- 1e+3


# Generate sample
x_vec <- generator(
  max_sample,
  v,
  pre_sampler = function(n){rbinom(n, 1, p_pre)},
  post_sampler = function(n){rbinom(n, 1, p_post)})

plot(1:max_sample, x_vec, pch=20,
     xlab = "n", ylab = "X_n", main = "Simulated Data")
graphics::lines(x = c(0,v), y = c(p_pre, p_pre), col = 2, lwd = 2)
graphics::lines(x = c(v,max_sample), y = c(p_post, p_post), col = 2, lwd =2 )

stopped_time <- run_quick_simulation(
  x_vec - p_pre,
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
