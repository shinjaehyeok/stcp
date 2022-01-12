# If EDCP is not installed then run below commands
# install.packages("devtools")
# devtools::install_github("shinjaehyeok/EDCP")

library(EDCP)

# Ratio case ----
# log(X): N(1, 1) + log(unit)
# log(Y): log(r) + log(X) + N(0,sig^2)
# Pre-change : r = 0.68 (with daily fluctuation)
# Post-change: r = 0.7 (with daily fluctuation)
# Change-point: v = 0 (immediate change), 200, 500, 1000 (no change)

max_sample <- 1000L
unit <- 10
sig <- 0.01
#1 - 0.68 / 0.8 = 0.15
m_pre <- 0.68
m_pre_model <- 0.685
#1 - 0.7 / 0.8 = 0.125
m_post <- 0.7
bound_lower <- 0.63
bound_upper <- 0.73
v <- 500
alpha <- 1 / (24 * 365)
k_max <- 1e+3


# Generate sample
x_vec <- unit * exp(rnorm(max_sample, 0,0.1))
r_pre <- filter(rnorm(v, 0, 0.001), filter=rep(1,24), circular=TRUE) +  m_pre
r_post <- filter(rnorm(max_sample-v, 0, 0.001), filter=rep(1,24), circular=TRUE) +  m_post
# r_post <- filter(rnorm(max_sample-v, 0, 0.001), filter=rep(1,24), circular=TRUE) +
#   seq(m_pre, 2 * m_post - m_pre, length.out = max_sample - v)
y_pre_vec <- log(r_pre) + log(x_vec[1:v]) + rnorm(v,0,sig)
y_post_vec <- log(r_post) + log(x_vec[seq(v+1, max_sample)]) + rnorm(max_sample-v,0,sig)
y_vec <- exp(c(y_pre_vec, y_post_vec) - sig^2 / 2)


ma_day <- filter(y_vec / x_vec, rep(1 / 24, 24), sides = 1)
ma_week <- filter(y_vec / x_vec, rep(1 / 24 / 7, 24 * 7), sides = 1)

plot(
  1:max_sample,
  y_vec / x_vec,
  pch = 20,
  xlab = "n",
  ylab = "Y/X",
  main = "Simulated Data",
  type = "l"
)
plot(
  1:max_sample,
  y_vec / x_vec,
  pch = 20,
  xlab = "n",
  ylab = "Y/X",
  main = "Simulated Data"
)
graphics::lines(
  x = 1:max_sample,
  y = c(r_pre,r_post),
  lty = 2
)
graphics::lines(
  x = 1:max_sample,
  y = ma_day,
  col = 2,
  lty = 2
)
graphics::lines(
  x = 1:max_sample,
  y = ma_week,
  col = 2
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
  # y = c(m_pre, 2 * m_post - m_pre),
  col = 2,
  lwd = 2
)
graphics::lines(
  x = c(1, max_sample),
  y = c((m_pre + m_post) / 2, (m_pre + m_post) / 2),
  col = 2,
  lwd = 2,
  lty = 2
)




# Compute optimal delta star
# Note delta_star for bounded case is not the same scale of delta_lower and delta_upper
# In the future, it may be better to use different terminology to avoid confusion.
# delta_star <-
#   (m_pre - bound_lower) * (m_post - m_pre) / ((m_post - m_pre) ^ 2 + 1)
delta_star <- m_post - m_pre
var_star <- var(y_vec / x_vec)
delta_upper <- bound_upper - m_pre_model
# delta_upper <- 0.1
delta_lower <- 0.01


# Build CP detectors
# When delta_lower = delta_upper = delta_star
edcp_star <- build_edcp_bounded(
  alpha,
  m_pre_model,
  delta_star,
  delta_star,
  is_test = FALSE,
  bound_lower,
  bound_upper,
  k_max = 1000,
  tol = 1e-6,
  var_lower = var_star,
  var_upper = var_star
)

# When delta_lower < delta_star < delta_upper
edcp_mix <- build_edcp_bounded(
  alpha,
  m_pre_model,
  delta_lower,
  delta_upper,
  is_test = FALSE,
  bound_lower,
  bound_upper,
  k_max = 1000,
  tol = 1e-6
)

ratio_vec <- y_vec / x_vec
ratio_vec[ratio_vec > bound_upper] <- bound_upper
ratio_vec[ratio_vec < bound_lower] <- bound_lower

stopped_time <- run_quick_simulation(ratio_vec,
                                     edcp_mix,
                                     v,
                                     edcp_star)
print(stopped_time)
print("daily_average")
print(min(which(ma_day > (m_pre + m_post) / 2)))
print("weekly_average")
print(min(which(ma_week > (m_pre + m_post) / 2)))
