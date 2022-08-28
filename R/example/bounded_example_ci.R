# If stcp is not installed then run below commands
# install.packages("devtools")
# devtools::install_github("shinjaehyeok/stcp")

library(stcp)
set.seed(1)
# Bernoulli case ----
# 5 * Beta(2, 3) - 1 (m = 1, sig^2 = 5^2 * 6 / (5^2 * 6) = 1)

max_sample <- 1000L
a <- 2
b <- 3
const <- 1
shift <- 0
m <- const * a / (a + b) + shift
sig <- const * sqrt(a * b / (a + b)^2 / (a + b + 1))
bound_lower <- shift
bound_upper <- const + shift
alpha <- 0.025


# Compute target interval
n_lower <- max_sample / 10
n_upper <- max_sample

# Generate data
x_vec <-  const * rbeta(max_sample, a, b) + shift
x_bar <- cumsum(x_vec) / seq_along(x_vec)

# Build CI
# 0. Baseline: z-score based CI
ci_fixed <-
  sig * qnorm(alpha, lower.tail = FALSE) / sqrt(1:max_sample)

# 1. Use sub-G and sub-B. They are sub-optimal but can be computed in the online way.
sub_G_sig <- (bound_upper - bound_lower) / 2
psi_fn_list_generator <- function() {
  generate_sub_G_fn(sub_G_sig)
}
ci_model_G <- build_ci_exp(alpha,
                           n_upper,
                           n_lower,
                           psi_fn_list_generator)

ci_G <- compute_ci(x_vec, ci_model_G, width_upper = sub_G_sig * 100, ci_lower_trivial = bound_lower)

ci_model_Ber <- build_ci_exp(alpha,
                             n_upper,
                             n_lower,
                             generate_sub_B_fn)

ci_Ber <- compute_ci((x_vec - shift) / const, ci_model_Ber, ci_lower_trivial = 0)
ci_Ber$x_bar <- ci_Ber$x_bar * const + shift
ci_Ber$ci_lower <- ci_Ber$ci_lower * const + shift

#2. Brute-force based on stcp_bounded model. Slower but tighter in general.
bf_ci_model <- build_bf_ci_bounded (alpha,
                                    n_upper,
                                    n_lower,
                                    bound_lower)

bf_ci <- compute_bf_ci(x_vec,
                       bf_ci_model,
                       ci_lower_trivial = bound_lower,
                       max_num_ci = 100)


# Plot CI
plot(
  1:max_sample,
  x_bar,
  type = "l",
  xlab = "n",
  ylab = "X_bar",
  main = "Running average",
  ylim = c(bound_lower, bound_upper),
)
graphics::abline(h = m,
                 col = 2,
                 lwd = 2)

graphics::lines(1:max_sample,
                x_bar - ci_fixed,
                lty = 2,
                col = 1)
graphics::lines(1:max_sample,
                ci_G$ci_lower,
                lty = 2,
                col = 2)
graphics::lines(1:max_sample,
                ci_Ber$ci_lower,
                lty = 2,
                col = 3)
graphics::lines(bf_ci$n, bf_ci$ci_lower, lty = 2, col = 4)

# Compare width
plot(bf_ci$n,
     (bf_ci$x_bar - bf_ci$ci_lower) / ci_fixed[bf_ci$n],
     type = "l",
     ylim = c(1, 2),
     col = 4)
graphics::lines(1:max_sample, (x_bar - ci_Ber$ci_lower) / ci_fixed, lty = 2, col = 3)
graphics::lines(1:max_sample, (x_bar - ci_G$ci_lower) / ci_fixed, lty = 2, col = 2)

bf_ci_model
(bf_ci$x_bar - bf_ci$ci_lower) / ci_fixed[bf_ci$n]
