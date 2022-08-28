# If stcp is not installed then run below commands
# install.packages("devtools")
# devtools::install_github("shinjaehyeok/stcp")
library(stcp)
# In this analysis, we use tidyverse to handle NBA data
library(tidyverse)
# Load NBA regular season log from 2011 to 2022
dat <- read.csv("R/example_code_for_paper/NBA_2010_2020.csv")
summary(dat$ptsTeam[dat$ptsTeam > 0])
summary(dat$plusminusTeam)
hist(dat$plusminusTeam)

# Get Cleveland Cavaliers Stats
CLE_dat <- dat %>% dplyr::filter(slugTeam == "CLE") %>%
  select(
    yearSeason,
    slugSeason,
    typeSeason,
    dateGame,
    nameTeam,
    slugTeam,
    isWin,
    ptsTeam,
    plusminusTeam
  ) %>%
  mutate(monthGame = format(as.Date(dateGame), "%Y-%m")) %>%
  mutate(relative_pm = plusminusTeam / (ptsTeam - plusminusTeam / 2))

CLE_dat_2010 <- CLE_dat %>% filter(yearSeason == 2010)
CLE_dat <-
  CLE_dat %>% filter(yearSeason > 2010 & yearSeason <= 2018)

hist(CLE_dat$plusminusTeam)


year_summ <- CLE_dat %>% group_by(yearSeason) %>%
  summarise(
    win_rate_year = mean(isWin),
    pts_year = mean(ptsTeam),
    pm_year = mean(plusminusTeam),
  )
month_summ <- CLE_dat %>% group_by(monthGame) %>%
  summarise(
    win_rate_month = mean(isWin),
    pts_month = mean(ptsTeam),
    pm_month = mean(plusminusTeam),
  )

CLE_dat <-
  CLE_dat %>% left_join(year_summ, by = "yearSeason") %>% left_join(month_summ, by = "monthGame")

regular_season_start_end_date <- CLE_dat %>%
  group_by(yearSeason) %>%
  summarise(start_date = min(dateGame),
            end_date = max(dateGame))

year_summ <-
  year_summ %>% left_join(regular_season_start_end_date, by = "yearSeason")

# NBA data ----
# 1. Win rate ----
# Pre-change : win_rate <= 0.49
# post_change: win_rate >= 0.51

plot(
  as.Date(CLE_dat$dateGame),
  CLE_dat$win_rate_month,
  pch = 20,
  xlab = "Date",
  ylab = "win rate",
  main = "Monthly win rates with seasonal averages"
)
for (i in 1:nrow(year_summ)) {
  year_date_range <-
    year_summ[i, c("start_date", "end_date")] %>% as.character() %>% as.Date()
  win_rate_year <- year_summ[i, "win_rate_year"] %>% as.numeric
  lines(
    x = year_date_range,
    y = rep(win_rate_year, 2),
    col = 2,
    lwd = 2
  )
}

# Build model
alpha <- 1e-3 # Inverse of target ARL
p_pre <- 0.49
delta_lower <- 0.02

ber_model <- build_stcp_exp(
  alpha,
  p_pre,
  delta_lower,
  delta_upper = 1 - (p_pre + delta_lower),
  psi_fn_list = generate_sub_B_fn(p_pre)
  )

# Compute mixtures of SR and CUSUM e-detectors.
mix_SR_ber <- run_stcp(CLE_dat$isWin, ber_model)
mix_CS_ber <- run_stcp(CLE_dat$isWin, ber_model, is_SR_type = FALSE)


# Stopping times of mixtures of e-SR and e-CUSUM procedures
mix_SR_stop_ber <- CLE_dat$dateGame[mix_SR_ber$stopped_ind] %>% as.Date()
mix_CS_stop_ber <- CLE_dat$dateGame[mix_CS_ber$stopped_ind] %>% as.Date()

# Plot log e-detector
# Size 600 * 450 for the paper
# plot(mix_SR_ber,
#      n = as.Date(CLE_dat$dateGame),
#      xlab = "Date",
#      main = paste0("CP detected at ", mix_SR_stop_ber, " (SR) / ", mix_CS_stop_ber, " (CUSUM)"),
#      draw_detect_line = FALSE,
#      col = 2) # unable default detection line

plot(as.Date(CLE_dat$dateGame),
     mix_SR_ber$log_mix_e_vec,
     xlab = "n (Date scale)",
     ylab = expression('log(M'['n']*')'),
     main = paste0("CP detected at ", mix_SR_stop_ber, " (SR) / ", mix_CS_stop_ber, " (CUSUM)"),
     type = "l",
     col = 2)
lines(as.Date(CLE_dat$dateGame),
      mix_CS_ber$log_mix_e_vec,
      col = 3)

# Draw customized detection line for the paper
abline(h = mix_SR_ber$stcp_obj$log_one_over_alpha)
abline(v = mix_SR_stop_ber, col = 2, lty = 2)
abline(v = mix_CS_stop_ber, col = 3, lty = 2)



# Draw detected line onto the original winning rate plot.
plot(
  as.Date(CLE_dat$dateGame),
  CLE_dat$win_rate_month,
  pch = 20,
  xlab = "Date",
  ylab = "win rate",
  main = "Monthly win rates with seasonal averages"
)
for (i in 1:nrow(year_summ)) {
  year_date_range <-
    year_summ[i, c("start_date", "end_date")] %>% as.character() %>% as.Date()
  win_rate_year <- year_summ[i, "win_rate_year"] %>% as.numeric
  lines(
    x = year_date_range,
    y = rep(win_rate_year, 2),
    col = 2,
    lwd = 2
  )
}
abline(v = mix_SR_stop_ber,
       col = 2,
       lty = 2,
       lwd = 2)




# 2. +/-  ----
# Pre-change : +/- =< -1
# Post-change: +/- > 1
# Assume +/- for each game is always between -80 and 80

plot(
  as.Date(CLE_dat$dateGame),
  CLE_dat$plusminusTeam,
  pch = 20,
  xlab = "Game Date",
  ylab = "+/-",
  main = "Plus-Minus of the Cavaliers"
)
for (i in 1:nrow(year_summ)) {
  year_date_range <-
    year_summ[i, c("start_date", "end_date")] %>% as.character() %>% as.Date()
  pm_year <- year_summ[i, "pm_year"] %>% as.numeric
  lines(
    x = year_date_range,
    y = rep(pm_year, 2),
    col = 2,
    lwd = 3
  )
}

# Build model
alpha <- 1e-3 # Inverse of target ARL
bound_lower <- -80
bound_upper <- 80
m_pre <- -1
delta_lower <- 2


bounded_model <- build_stcp_bounded(
  alpha,
  m_pre,
  delta_lower,
  bound_lower = bound_lower,
  bound_upper = bound_upper
)

# Compute mixtures of SR and CUSUM e-detectors.
mix_SR_bounded <- run_stcp(CLE_dat$plusminusTeam, bounded_model)
mix_CS_bounded <- run_stcp(CLE_dat$plusminusTeam, bounded_model, is_SR_type = FALSE)


# Stopping time of the mixture of SR procedure.
mix_SR_stop_bounded <- CLE_dat$dateGame[mix_SR_bounded$stopped_ind] %>% as.Date()
mix_CS_stop_bounded <- CLE_dat$dateGame[mix_CS_bounded$stopped_ind] %>% as.Date()


# Plot log e-detector
# Size 600 * 450 for the paper
# plot(mix_SR_bounded,
#      n = as.Date(CLE_dat$dateGame),
#      xlab = "Date",
#      main = paste0("CP detected at ", mix_SR_stop_bounded, " (SR) / ", mix_CS_stop_bounded, " (CUSUM)"),
#      draw_detect_line = FALSE, # unable default detection line
#      col = 2)
# plot(mix_CS_bounded,
#      n = as.Date(CLE_dat$dateGame),
#      draw_detect_line = FALSE,
#      add = TRUE,
#      col = 3)

plot(as.Date(CLE_dat$dateGame),
     mix_SR_bounded$log_mix_e_vec,
     xlab = "n (Date scale)",
     ylab = expression('log(M'['n']*')'),
     main = paste0("CP detected at ", mix_SR_stop_bounded, " (SR) / ", mix_CS_stop_bounded, " (CUSUM)"),
     type = "l",
     col = 2)
lines(as.Date(CLE_dat$dateGame),
      mix_CS_bounded$log_mix_e_vec,
      col = 3)

# Draw customized detection line for the paper
abline(h = mix_SR_bounded$stcp_obj$log_one_over_alpha)
abline(v = mix_SR_stop_bounded, col = 2, lty = 2)
abline(v = mix_CS_stop_bounded, col = 3, lty = 2)

# M_n scale
plot(seq_along(mix_SR_bounded$log_mix_e_vec),
     exp(mix_SR_bounded$log_mix_e_vec),
     xlab = "n",
     ylab = expression('M'['n']),
     main = paste0("CP detected at ", mix_SR_bounded$stopped_ind, "-th game (SR) / ", mix_CS_bounded$stopped_ind, "-th game (CUSUM)"),
     type = "l",
     col = 2,
     ylim = c(0, 4 / alpha))
lines(seq_along(mix_CS_bounded$log_mix_e_vec),
      exp(mix_CS_bounded$log_mix_e_vec),
      col = 3)

# Draw customized detection line for the paper
abline(h = exp(mix_SR_bounded$stcp_obj$log_one_over_alpha))
abline(v = mix_SR_bounded$stopped_ind, col = 2, lty = 2)
abline(v = mix_CS_bounded$stopped_ind, col = 3, lty = 2)



# Draw detected line onto the original +/- plot.
plot(
  as.Date(CLE_dat$dateGame),
  CLE_dat$plusminusTeam,
  pch = 20,
  xlab = "Date",
  ylab = "+/-",
  main = "+/- with detected CP"
)
for (i in 1:nrow(year_summ)) {
  year_date_range <-
    year_summ[i, c("start_date", "end_date")] %>% as.character() %>% as.Date()
  pm_year <- year_summ[i, "pm_year"] %>% as.numeric
  lines(
    x = year_date_range,
    y = rep(pm_year, 2),
    col = 2,
    lwd = 2
  )
}

abline(v = mix_SR_stop_bounded,
       col = 2,
       lty = 2,
       lwd = 2)


# Uniform mixture example
# Normalize observations
alpha <- 1e-3 # Inverse of target ARL
bound_lower <- -80
bound_upper <- 80
m_pre <- -1
m_pre_normalized <-(m_pre - bound_lower) / (bound_upper - bound_lower)
pm_normalized <- (CLE_dat$plusminusTeam - bound_lower) / (bound_upper - bound_lower)

# Build model
# Compute parameters
by_gap <- 0.02
lambda <- seq(by_gap, 1, by = by_gap)
omega <- rep(1/length(lambda), length(lambda))

stcp_uniform <- build_stcp(alpha,
                           m_pre_normalized,
                           is_test = FALSE,
                           omega = omega,
                           lambda = lambda,
                           log_base_fn_generator = generate_log_bounded_base_fn,
                           m = m_pre_normalized,
                           bound_lower = 0)

# Compute mixtures of SR and CUSUM e-detectors.
uniform_SR <- run_stcp(pm_normalized, stcp_uniform)
uniform_CS <- run_stcp(pm_normalized, stcp_uniform, is_SR_type = FALSE)


# Stopping time of the mixture of SR procedure.
uniform_SR_stop <- CLE_dat$dateGame[uniform_SR$stopped_ind] %>% as.Date()
uniform_CS_stop <- CLE_dat$dateGame[uniform_CS$stopped_ind] %>% as.Date()


plot(as.Date(CLE_dat$dateGame),
     uniform_SR$log_mix_e_vec,
     xlab = "n (Date scale)",
     ylab = expression('log(M'['n']*')'),
     main = paste0("CP detected at ", uniform_SR_stop, " (SR) / ", uniform_CS_stop, " (CUSUM)"),
     type = "l",
     col = 2)
lines(as.Date(CLE_dat$dateGame),
      uniform_CS$log_mix_e_vec,
      col = 3)
lines(as.Date(CLE_dat$dateGame),
      mix_SR_bounded$log_mix_e_vec,
      col = 4, lty = 2)
lines(as.Date(CLE_dat$dateGame),
      mix_CS_bounded$log_mix_e_vec,
      col = 5, lty = 2)

# Draw customized detection line for the paper
abline(h = uniform_SR$stcp_obj$log_one_over_alpha)
abline(v = uniform_SR_stop, col = 2, lty = 2)
abline(v = uniform_CS_stop, col = 3, lty = 2)
abline(v = mix_SR_stop_bounded, col = 4, lty = 2)
abline(v = mix_CS_stop_bounded, col = 5, lty = 2)
