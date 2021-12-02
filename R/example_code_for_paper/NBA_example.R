# If EDCP is not installed then run below commands
# install.packages("devtools")
# devtools::install_github("shinjaehyeok/EDCP")
library(EDCP)
# In this analysis, we use tidyverse to handle NBA data
library(tidyverse)
# Load NBA regular season log from 2011 to 2022
dat <- read.csv("R/example_code_for_paper/NBA_2010_2020.csv")
summary(dat$ptsTeam[dat$ptsTeam > 0 ])
summary(dat$plusminusTeam)
hist(dat$plusminusTeam)

# Get Cleveland Cavaliers Stats
CLE_dat <- dat %>% dplyr::filter(slugTeam == "CLE") %>%
  select(yearSeason, slugSeason, typeSeason, dateGame, nameTeam, slugTeam, isWin, ptsTeam, plusminusTeam) %>%
  mutate(monthGame = format(as.Date(dateGame), "%Y-%m")) %>%
  mutate(
  relative_pm = plusminusTeam / (ptsTeam - plusminusTeam/2)
)

CLE_dat_2010 <- CLE_dat %>% filter(yearSeason == 2010)
CLE_dat <- CLE_dat %>% filter(yearSeason > 2010 & yearSeason <= 2018)

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

CLE_dat <- CLE_dat %>% left_join(year_summ, by = "yearSeason") %>% left_join(month_summ, by = "monthGame")

regular_season_start_end_date <- CLE_dat %>%
  group_by(yearSeason) %>%
  summarise(start_date = min(dateGame), end_date = max(dateGame))

year_summ <- year_summ %>% left_join(regular_season_start_end_date, by = "yearSeason")

# NBA data ----
# 1. Win rate ----
# Pre-change : win_rate <= 0.49
# post_change: win_rate >= 0.51

plot(as.Date(CLE_dat$dateGame), CLE_dat$win_rate_month, pch=20,
     xlab = "Date", ylab = "win rate", main = "Monthly win rate with seasonaly average")
for (i in 1:nrow(year_summ)){
  year_date_range <- year_summ[i,c("start_date", "end_date")] %>% as.character() %>% as.Date()
  win_rate_year <- year_summ[i, "win_rate_year"] %>% as.numeric
  lines(x = year_date_range, y = rep(win_rate_year, 2), col = 2, lwd = 2)
}


alpha <- 1e-3 # Inverse of target ARL
p_pre <- 0.49

# Compute parameters
base_param <- compute_baseline(
  alpha = alpha,
  delta_lower = 0.02,
  delta_upper = 0.3,
  psi_fn_list = generate_sub_B_fn(p = p_pre),
  v_min = 1,
  k_max = 1000
)

# Compute e-detectors
log_base_fn_list <- sapply(base_param$lambda,
                           generate_log_base_fn,
                           psi_fn = base_param$psi_fn_list$psi)

# Compute mixture of SR-type e-detectors.
mix_SR <- update_log_mix_e_detectors(CLE_dat$isWin - p_pre,
                                     base_param$omega,
                                     log_base_fn_list)

# Stopping time of the mixture of SR procedure.
mix_SR_stop <- min(
  CLE_dat$dateGame[which(mix_SR$log_mix_e_detect_val > log(1/alpha))]
) %>% as.Date()

plot(as.Date(CLE_dat$dateGame), mix_SR$log_mix_e_detect_val, type = "l",
     xlab = "Date", ylab = "log e-detectors", main = paste0("CP detected at ",mix_SR_stop))
abline(h = log(1/alpha), col = 2)
abline(v = as.Date(regular_season_start_end_date$start_date) , col = 1, lty = 2)
abline(v = as.Date(regular_season_start_end_date$end_date) , col = 1, lty = 3)
abline(v = mix_SR_stop, col = 2, lty = 2)

plot(as.Date(CLE_dat$dateGame), CLE_dat$win_rate_month, pch=20,
     xlab = "Date", ylab = "win rate", main = "Monthly win rate with seasonaly average")
for (i in 1:nrow(year_summ)){
  year_date_range <- year_summ[i,c("start_date", "end_date")] %>% as.character() %>% as.Date()
  win_rate_year <- year_summ[i, "win_rate_year"] %>% as.numeric
  lines(x = year_date_range, y = rep(win_rate_year, 2), col = 2, lwd = 2)
}
abline(v = mix_SR_stop, col = 2, lty = 2, lwd = 2)


# 2. +/-  ----
# Pre-change : +/- =< -1
# Post-change: +/- > 1
# Assume +/- for each game is always between -100 and 100

plot(as.Date(CLE_dat$dateGame), CLE_dat$plusminusTeam, pch=20,
     xlab = "Game Date", ylab = "+/-", main = "Plus-Minus of the Cavaliers")
for (i in 1:nrow(year_summ)){
  year_date_range <- year_summ[i,c("start_date", "end_date")] %>% as.character() %>% as.Date()
  pm_year <- year_summ[i, "pm_year"] %>% as.numeric
  lines(x = year_date_range, y = rep(pm_year, 2), col = 2, lwd = 3)
}

lower <- -80
upper <- 80

CLE_dat <- CLE_dat %>% dplyr::mutate(normalized_pm = (plusminusTeam - lower) / (upper - lower))

alpha <- 1e-3 # Inverse of target ARL
m <- (-1 - lower) / (upper - lower) # Upper bound of mean of pre-change observations
d <- 2 / (upper - lower)  # Guess on the minimum gap between pre- and post-change means
E_fn_list <- generate_sub_E_fn()

# Compute parameters
base_param <- compute_baseline(
  alpha = alpha,
  delta_lower = m * d / (1/4 + (1-m)^2), # 0.012
  delta_upper = m * (1-m) / d^2,  # 1599.8
  psi_fn_list = generate_sub_E_fn(),
  v_min = 0,
  k_max = 1000
)

# Compute e-detectors
log_base_fn_list <- sapply(base_param$lambda,
                           generate_log_bounded_base_fn,
                           m = m)

# Compute mixture of SR-type e-detectors.
mix_SR <- update_log_mix_e_detectors(CLE_dat$normalized_pm,
                                     base_param$omega,
                                     log_base_fn_list)

# Stopping time of the mixture of SR procedure.
mix_SR_stop <- min(
  CLE_dat$dateGame[which(mix_SR$log_mix_e_detect_val > log(1/alpha))]
) %>% as.Date()

plot(as.Date(CLE_dat$dateGame), mix_SR$log_mix_e_detect_val, type = "l",
     xlab = "Game Date", ylab = "log e-detectors", main = paste0("CP detected at ",mix_SR_stop))
abline(h = log(1/alpha), col = 2)
# abline(v = as.Date(regular_season_start_end_date$start_date) , col = 1, lty = 2)
# abline(v = as.Date(regular_season_start_end_date$end_date) , col = 1, lty = 3)
abline(v = mix_SR_stop, col = 2, lty = 2)

plot(as.Date(CLE_dat$dateGame), CLE_dat$plusminusTeam, pch=20,
     xlab = "Date", ylab = "+/-", main = "+/- with detected CP")
for (i in 1:nrow(year_summ)){
  year_date_range <- year_summ[i,c("start_date", "end_date")] %>% as.character() %>% as.Date()
  pm_year <- year_summ[i, "pm_year"] %>% as.numeric
  lines(x = year_date_range, y = rep(pm_year, 2), col = 2, lwd = 2)
}

abline(v = mix_SR_stop, col = 2, lty = 2, lwd = 2)


# # 3. Score ----
# # Pre-change : PTS =< 99
# # Post-change: PTS >= 100
# # Assume points for each game is always between 50 and 200
#
# plot(as.Date(CLE_dat$dateGame), CLE_dat$ptsTeam, pch=20,
#      xlab = "Date", ylab = "X_n", main = "Game points with seasonal averages")
# for (i in 1:nrow(year_summ)){
#   year_date_range <- year_summ[i,c("start_date", "end_date")] %>% as.character() %>% as.Date()
#   pts_year <- year_summ[i, "pts_year"] %>% as.numeric
#   lines(x = year_date_range, y = rep(pts_year, 2), col = 2, lwd = 2)
# }
#
# CLE_dat <- CLE_dat %>% dplyr::mutate(normalized_PTS = (ptsTeam - 50) / (200 - 50))
#
# alpha <- 1e-3 # Inverse of target ARL
# m <- (99 - 50) / (200-50) # Upper bound of mean of pre-change observations
# d <- 1 / 150  # Guess on the minimum gap between pre- and post-change means
# E_fn_list <- generate_sub_E_fn()
#
# # Compute parameters
# base_param <- compute_baseline(
#   alpha = alpha,
#   delta_lower = m * d / (1/4 + (1-m)^2), # 0.003
#   delta_upper = m * (1-m) / d^2,  # 4949
#   psi_fn_list = generate_sub_E_fn(),
#   v_min = 0,
#   k_max = 1000
# )
#
# # Compute e-detectors
# log_base_fn_list <- sapply(base_param$lambda,
#                            generate_log_bounded_base_fn,
#                            m = m)
#
# # Compute mixture of SR-type e-detectors.
# mix_SR <- update_log_mix_e_detectors(CLE_dat$normalized_PTS,
#                                      base_param$omega,
#                                      log_base_fn_list)
#
# # Inferred change-point
# mix_SR_stop <- min(
#   CLE_dat$dateGame[which(mix_SR$log_mix_e_detect_val > log(1/alpha))]
# ) %>% as.Date()
#
# plot(as.Date(CLE_dat$dateGame), mix_SR$log_mix_e_detect_val, type = "l",
#      xlab = "Date", ylab = "log e-detectors", main = paste0("v_hat = ",mix_SR_stop))
# abline(h = log(1/alpha), col = 2)
# abline(v = as.Date(regular_season_start_end_date$start_date) , col = 1, lty = 2)
# abline(v = as.Date(regular_season_start_end_date$end_date) , col = 1, lty = 3)
# abline(v = mix_SR_stop, col = 2, lty = 2)
#
# plot(as.Date(CLE_dat$dateGame), CLE_dat$ptsTeam, pch=20,
#      xlab = "Date", ylab = "X_n", main = "Game points with the detected CP")
# for (i in 1:nrow(year_summ)){
#   year_date_range <- year_summ[i,c("start_date", "end_date")] %>% as.character() %>% as.Date()
#   pts_year <- year_summ[i, "pts_year"] %>% as.numeric
#   lines(x = year_date_range, y = rep(pts_year, 2), col = 2, lwd = 2)
# }
# abline(v = mix_SR_stop, col = 2, lty = 2, lwd = 2)
#
# # 3-2. Score change
# # Pre-change : PTS =< 99
# # Post-change: PTS >= 100
# # Assume points for each game is always between 50 and 200
# year_summ_after_2012 <- year_summ  %>%
#   mutate(
#     pts_pre_year = c(first(pts_year),
#                      pts_year[-length(pts_year)])
#   ) %>%
#   filter(yearSeason > 2011)
#
# CLE_dat_after_2012 <- CLE_dat %>%
#   left_join(year_summ_after_2012 %>%
#               select(yearSeason, pts_pre_year),
#             by = "yearSeason") %>%
#   mutate(
#     normalized_PTS_change = (ptsTeam - 50) / (pts_pre_year - 50),
#     pts_year_change = pts_year / pts_pre_year) %>%
#   filter(yearSeason > 2011)
#
#
#
# plot(as.Date(CLE_dat_after_2012$dateGame), CLE_dat_after_2012$normalized_PTS_change, pch=20,
#      xlab = "Date", ylab = "X_n", main = "Game points with seasonal averages")
# for (i in 1:nrow(year_summ_after_2012)){
#   year_date_range <- year_summ_after_2012[i,c("start_date", "end_date")] %>% as.character() %>% as.Date()
#   pts_year <- year_summ_after_2012[i, "pts_year"] %>% as.numeric()
#   pts_pre_year <- year_summ_after_2012[i, "pts_pre_year"] %>% as.numeric()
#   lines(x = year_date_range, y = rep(pts_year / pts_pre_year, 2), col = 2, lwd = 2)
# }
#
# alpha <- 1e-3 # Inverse of target ARL
# r_upper <- 1.02 # Upper bound of ratio of pre-change observations
# d <- 0.01  # Guess on the minimum gap between pre- and post-change means
# E_fn_list <- generate_sub_E_fn()
#
# # Compute parameters
# base_param <- compute_baseline(
#   alpha = alpha,
#   delta_lower = 0.01, # 0.003
#   delta_upper = 1,  # 4949
#   psi_fn_list = generate_sub_E_fn(),
#   v_min = 0,
#   k_max = 1000
# )
#
# generate_log_ratio_base_fn <- function(lambda, r = 1.03){
#   if (r <= 0) stop("Ratio parameter must be strictly positive.")
#   log_base_fn <- function(x){
#     log(1 + lambda * (x / r  - 1))
#   }
#   return(log_base_fn)
# }
#
# # Compute e-detectors
# log_base_fn_list <- sapply(base_param$lambda,
#                            generate_log_ratio_base_fn,
#                            r = r_upper)
#
# # Compute mixture of SR-type e-detectors.
# mix_SR <- update_log_mix_e_detectors(CLE_dat_after_2012$normalized_PTS_change,
#                                      base_param$omega,
#                                      log_base_fn_list)
#
# # Inferred change-point
# mix_SR_stop <- min(
#   CLE_dat_after_2012$dateGame[which(mix_SR$log_mix_e_detect_val > log(1/alpha))]
# ) %>% as.Date()
#
# plot(as.Date(CLE_dat_after_2012$dateGame), mix_SR$log_mix_e_detect_val, type = "l",
#      xlab = "Date", ylab = "log e-detectors", main = paste0("v_hat = ",mix_SR_stop))
# abline(h = log(1/alpha), col = 2)
# abline(v = as.Date(regular_season_start_end_date$start_date) , col = 1, lty = 2)
# abline(v = as.Date(regular_season_start_end_date$end_date) , col = 1, lty = 3)
# abline(v = mix_SR_stop, col = 2, lty = 2)
#
# plot(as.Date(CLE_dat_after_2012$dateGame), CLE_dat_after_2012$ptsTeam, pch=20,
#      xlab = "Date", ylab = "X_n", main = "Game points with the detected CP")
# for (i in 1:nrow(year_summ)){
#   year_date_range <- year_summ[i,c("start_date", "end_date")] %>% as.character() %>% as.Date()
#   pts_year <- year_summ[i, "pts_year"] %>% as.numeric
#   lines(x = year_date_range, y = rep(pts_year, 2), col = 2, lwd = 2)
# }
# abline(v = mix_SR_stop, col = 2, lty = 2, lwd = 2)
#
#
