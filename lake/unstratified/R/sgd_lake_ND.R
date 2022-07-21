# *************************
# Rn mass balance calculation for time series lake radon budget
# Author: Peter Fuleky
# Date: 2022-05-22
# notes:
# *************************

# *************************
#  setup ----
# *************************

# rm(list = ls(pattern = glob2rx("*__*")))

# load packages, functions, etc.
source(here::here("R", "setup.R"))

# study type (determines folder)
study_folder <- "lake/unstratified"

# input file name
csv_file_in <- "sgd_lake_unstratified_dataND1.csv"

# *************************
#  load data ----
# *************************

# load the time series
in_tbl <- read_csv(here(study_folder, "input", csv_file_in)) %>%
  # mutate(across(.cols = 1, .fns = ~ force_tz(., tzone = ""))) %>%
  mutate(across(.cols = 1, .fns = ~ lubridate::parse_date_time(., c("ymdHMS", "ymdHM", "mdyHM", "mdyHMS")))) %>%
  mutate(across(.cols = -1, .fns = as.numeric))
  # mutate(across(.cols = 2:4, .fns = ~ lubridate::parse_date_time(., c("ymdHMS", "ymdHM", "mdyHM", "mdyHMS")))) %>%
  # mutate(across(.cols = -(2:4), .fns = as.numeric))

# names(in_tbl)

# *************************
#  calculations ----
# *************************

# create data table
dat_tbl <- in_tbl %>%
  mutate(
    
    # calculates radon measurement time interval in hours based on provided measurement times, 
    # all other time series measurements provided by the user should be averaged to this interval
    time = `Rad Mid Time`,
    meas_t__h = (time %>% as.numeric() - lag(time %>% as.numeric())) / 60 / 60,
    tot_t__h = slider::slide_dbl(meas_t__h, ~sum(.x , na.rm = TRUE), .before = Inf),
    tot_t__day = tot_t__h / 24,
    wind__ms = `Wind speed (m/s)`,
    temp_wat__C = lag(`Water Temperature C`),
    temp_wat__K = temp_wat__C + 273.15,
    L = 10^(-(980 / (273 + temp_wat__C) + 1.59)),
    M = (999.842594 + (0.06793952) * temp_wat__C - (0.00909529) * temp_wat__C^2 + (0.0001001685) * temp_wat__C^3 - (0.000001120083) * temp_wat__C^4 + (0.000000006536332) * temp_wat__C^5) / 1000,
    N = 0.00002414 * 10^(247.8 / (temp_wat__K - 140)) * 1000 / 100,
    O = N / M,
    P = O / L,
    S = ((0.45 * (wind__ms^1.6) * ((P / 600)^-0.667)) / 100) / 60,
    Q = 0.105 + 0.405 * exp(-0.05027 * temp_wat__C),
    V = lag(`Rn in Air dpm/m3`),
    U = lag(`Water Level m`),
    user_Rn_atm__dpmm3 = 800,
    user_depth__m = 1.3,
    # decay constant of Rn in hours
    user_lambda__hr = log(2) / (3.84 * 24),
    user_Rn_ss__dpmm3 = 3000, # driver of optimization
    user_Rn_ini_inv__dpmm2 = user_depth__m * user_Rn_ss__dpmm3,
    user_Rn_ss_flux__dpmm2hr = user_lambda__hr * user_Rn_ini_inv__dpmm2, # calculate from components
    user_Rn_ss_inv__dpmm2 = user_Rn_ss_flux__dpmm2hr / user_lambda__hr,
    user_kin_visc__cm2s = 0.01004,
    user_mol_diff__cm2s = 0.000012,
    X = user_kin_visc__cm2s / O,
    Y = (34.6 * X * (user_mol_diff__cm2s^0.5)*(wind__ms^1.5)) / (24 * 60),
    W = c(`Rn in Water dpm/m3`[2], rep(NA_real_, length(time) - 1)), #`Rn in Water dpm/m3`,
  )

for(t_ind in seq_along(dat_tbl$time)[-1]) {
  dat_tbl$W[t_ind] = dat_tbl %$% {W[t_ind - 1] * exp(-user_lambda__hr[t_ind] * meas_t__h[t_ind]) + (user_Rn_ss_flux__dpmm2hr[t_ind] - (W[t_ind-1] - (Q[t_ind] * V[t_ind])) * S[t_ind] * 60) * (1 - exp(-user_lambda__hr[t_ind] * meas_t__h[t_ind])) / user_lambda__hr[t_ind] / U[t_ind]}
  }

dat_tbl <- dat_tbl %>%
  mutate(
    R = lag(W) - (Q*V),
    T = R * S * 60,
    Z = R * Y * 60,
    # there is an error in the xls: not using lambda in the next line:
    AB = c((lag(W) * exp(-user_lambda__hr * meas_t__h) * U)[1:2], rep(NA_real_, length(time) - 2)),
    AC = user_Rn_ss_flux__dpmm2hr * (1 - exp(-user_lambda__hr * meas_t__h)) / user_lambda__hr,
    AD = T * (1 - exp(-user_lambda__hr * meas_t__h)) / user_lambda__hr,
    AE = AB + (AC - AD),
  )

for(t_ind in seq_along(dat_tbl$time)[-(1:2)]) {
  dat_tbl$AB[t_ind] = dat_tbl %$% {AE[t_ind - 1] * exp(-user_lambda__hr[t_ind] * meas_t__h[t_ind])}
  dat_tbl$AE[t_ind] = dat_tbl %$% {AB[t_ind] + (AC[t_ind] - AD[t_ind])}
}

dat_tbl <- dat_tbl %>%
  mutate(
    AF = AE / user_depth__m,
    AG = Z * (1 - exp(-user_lambda__hr * meas_t__h)) / user_lambda__hr,
    AH = AB + (AC - AG),
    AI = AH / U,
    AJ = Q/Q[2],
    AK = c((AF * user_depth__m)[1:2], (AF * AJ)[-(1:2)]),
    AL = AI * AJ, 
    # AM = `Rn in Water dpm/m3`,
    AN = (AK - `Rn in Water dpm/m3`)^2,
    AO = (AL - `Rn in Water dpm/m3`)^2,
    
  ) %>% 
  
  # drop columns with no values (keep those with at least one value)
  select(where(~(!(.x %>% is.na())) %>% any()))

RMSE_vals <- dat_tbl %>% 
  summarize(
    RMSE_old = mean(AN, na.rm = TRUE) %>% sqrt(),
    RMSE_new = mean(AO, na.rm = TRUE) %>% sqrt(),
  )

# results saved in a csv file
write_csv(dat_tbl, here(study_folder, "output", "rn_budgetND.csv"), na = "")

# *************************
#  optimization ----
# *************************

# objective function to find steady state Rn flux
Rn_ss_calc <- function(optim_driver){
  #input: ini_Rn_flux__Bqm2hr - initial value for Rn flux via SGD
  
  # create data table
  dat_tbl <- in_tbl %>%
    mutate(
      
      # calculates radon measurement time interval in hours based on provided measurement times, 
      # all other time series measurements provided by the user should be averaged to this interval
      time = `Rad Mid Time`,
      meas_t__h = (time %>% as.numeric() - lag(time %>% as.numeric())) / 60 / 60,
      tot_t__h = slider::slide_dbl(meas_t__h, ~sum(.x , na.rm = TRUE), .before = Inf),
      tot_t__day = tot_t__h / 24,
      wind__ms = `Wind speed (m/s)`,
      temp_wat__C = lag(`Water Temperature C`),
      temp_wat__K = temp_wat__C + 273.15,
      L = 10^(-(980 / (273 + temp_wat__C) + 1.59)),
      M = (999.842594 + (0.06793952) * temp_wat__C - (0.00909529) * temp_wat__C^2 + (0.0001001685) * temp_wat__C^3 - (0.000001120083) * temp_wat__C^4 + (0.000000006536332) * temp_wat__C^5) / 1000,
      N = 0.00002414 * 10^(247.8 / (temp_wat__K - 140)) * 1000 / 100,
      O = N / M,
      P = O / L,
      S = ((0.45 * (wind__ms^1.6) * ((P / 600)^-0.667)) / 100) / 60,
      Q = 0.105 + 0.405 * exp(-0.05027 * temp_wat__C),
      V = lag(`Rn in Air dpm/m3`),
      U = lag(`Water Level m`),
      user_Rn_atm__dpmm3 = 800,
      user_depth__m = 1.3,
      # decay constant of Rn in hours
      user_lambda__hr = log(2) / (3.84 * 24),
      user_Rn_ss__dpmm3 = optim_driver, # driver of optimization
      user_Rn_ini_inv__dpmm2 = user_depth__m * user_Rn_ss__dpmm3,
      user_Rn_ss_flux__dpmm2hr = user_lambda__hr * user_Rn_ini_inv__dpmm2, # calculate from components
      user_Rn_ss_inv__dpmm2 = user_Rn_ss_flux__dpmm2hr / user_lambda__hr,
      user_kin_visc__cm2s = 0.01004,
      user_mol_diff__cm2s = 0.000012,
      X = user_kin_visc__cm2s / O,
      Y = (34.6 * X * (user_mol_diff__cm2s^0.5)*(wind__ms^1.5)) / (24 * 60),
      W = c(`Rn in Water dpm/m3`[2], rep(NA_real_, length(time) - 1)), #`Rn in Water dpm/m3`,
    )
  
  for(t_ind in seq_along(dat_tbl$time)[-1]) {
    dat_tbl$W[t_ind] = dat_tbl %$% {W[t_ind - 1] * exp(-user_lambda__hr[t_ind] * meas_t__h[t_ind]) + (user_Rn_ss_flux__dpmm2hr[t_ind] - (W[t_ind-1] - (Q[t_ind] * V[t_ind])) * S[t_ind] * 60) * (1 - exp(-user_lambda__hr[t_ind] * meas_t__h[t_ind])) / user_lambda__hr[t_ind] / U[t_ind]}
  }
  
  dat_tbl <- dat_tbl %>%
    mutate(
      R = lag(W) - (Q*V),
      T = R * S * 60,
      Z = R * Y * 60,
      # there is an error in the xls: not using lambda in the next line:
      AB = c((lag(W) * exp(-user_lambda__hr * meas_t__h) * U)[1:2], rep(NA_real_, length(time) - 2)),
      AC = user_Rn_ss_flux__dpmm2hr * (1 - exp(-user_lambda__hr * meas_t__h)) / user_lambda__hr,
      AD = T * (1 - exp(-user_lambda__hr * meas_t__h)) / user_lambda__hr,
      AE = AB + (AC - AD),
    )
  
  for(t_ind in seq_along(dat_tbl$time)[-(1:2)]) {
    dat_tbl$AB[t_ind] = dat_tbl %$% {AE[t_ind - 1] * exp(-user_lambda__hr[t_ind] * meas_t__h[t_ind])}
    dat_tbl$AE[t_ind] = dat_tbl %$% {AB[t_ind] + (AC[t_ind] - AD[t_ind])}
  }
  
  dat_tbl <- dat_tbl %>%
    mutate(
      AF = AE / user_depth__m,
      AG = Z * (1 - exp(-user_lambda__hr * meas_t__h)) / user_lambda__hr,
      AH = AB + (AC - AG),
      AI = AH / U,
      AJ = Q/Q[2],
      AK = c((AF * user_depth__m)[1:2], (AF * AJ)[-(1:2)]),
      AL = AI * AJ, 
      # AM = `Rn in Water dpm/m3`,
      AN = (AK - `Rn in Water dpm/m3`)^2,
      AO = (AL - `Rn in Water dpm/m3`)^2,
      
    ) %>% 
    
    # drop columns with no values (keep those with at least one value)
    select(where(~(!(.x %>% is.na())) %>% any()))
  
  RMSE_vals <- dat_tbl %>% 
    summarize(
      RMSE_old = mean(AN, na.rm = TRUE) %>% sqrt(),
      RMSE_new = mean(AO, na.rm = TRUE) %>% sqrt(),
    ) %>% 
    
    pull(RMSE_new)
  
}

# minimize objective function with respect to steady state Rn flux
# initial value is based on:
# 1) initial guess of steady state Rn concentration
# Rn_ss__Bqm3 = Rn_wat__Bqm3,
# 2) Rn inventory in steady state
# inv_Rn_ss__Bqm2 = Rn_ss__Bqm3 * depth__m,
# 3) initial guess of steady state Rn flux
# Rn_flux__Bqm2hr = lambda * inv_Rn_ss__Bqm2,
# combine formulas 1 -> 2 -> 3 and take average:
optim_start <- dat_tbl %>% 
  summarize(mean(user_Rn_ss__dpmm3, na.rm = TRUE)) %>% 
  pull()

# compare two optimization techniques
opt1 <- optim(optim_start, Rn_ss_calc, method = "BFGS")
opt2 <- optimize(Rn_ss_calc, c(0, 50000))
opt1$par - opt2$minimum


# *************************
#  calcs with optimal value ----
# *************************

optim_driver <- opt1$par

# create data table
dat_tbl <- in_tbl %>%
  mutate(
    
    # calculates radon measurement time interval in hours based on provided measurement times, 
    # all other time series measurements provided by the user should be averaged to this interval
    time = `Rad Mid Time`,
    meas_t__h = (time %>% as.numeric() - lag(time %>% as.numeric())) / 60 / 60,
    tot_t__h = slider::slide_dbl(meas_t__h, ~sum(.x , na.rm = TRUE), .before = Inf),
    tot_t__day = tot_t__h / 24,
    wind__ms = `Wind speed (m/s)`,
    temp_wat__C = lag(`Water Temperature C`),
    temp_wat__K = temp_wat__C + 273.15,
    L = 10^(-(980 / (273 + temp_wat__C) + 1.59)),
    M = (999.842594 + (0.06793952) * temp_wat__C - (0.00909529) * temp_wat__C^2 + (0.0001001685) * temp_wat__C^3 - (0.000001120083) * temp_wat__C^4 + (0.000000006536332) * temp_wat__C^5) / 1000,
    N = 0.00002414 * 10^(247.8 / (temp_wat__K - 140)) * 1000 / 100,
    O = N / M,
    P = O / L,
    S = ((0.45 * (wind__ms^1.6) * ((P / 600)^-0.667)) / 100) / 60,
    Q = 0.105 + 0.405 * exp(-0.05027 * temp_wat__C),
    V = lag(`Rn in Air dpm/m3`),
    U = lag(`Water Level m`),
    user_Rn_atm__dpmm3 = 800,
    user_depth__m = 1.3,
    # decay constant of Rn in hours
    user_lambda__hr = log(2) / (3.84 * 24),
    user_Rn_ss__dpmm3 = optim_driver, # driver of optimization
    user_Rn_ini_inv__dpmm2 = user_depth__m * user_Rn_ss__dpmm3,
    user_Rn_ss_flux__dpmm2hr = user_lambda__hr * user_Rn_ini_inv__dpmm2, # calculate from components
    user_Rn_ss_inv__dpmm2 = user_Rn_ss_flux__dpmm2hr / user_lambda__hr,
    user_kin_visc__cm2s = 0.01004,
    user_mol_diff__cm2s = 0.000012,
    X = user_kin_visc__cm2s / O,
    Y = (34.6 * X * (user_mol_diff__cm2s^0.5)*(wind__ms^1.5)) / (24 * 60),
    W = c(`Rn in Water dpm/m3`[2], rep(NA_real_, length(time) - 1)), #`Rn in Water dpm/m3`,
  )

for(t_ind in seq_along(dat_tbl$time)[-1]) {
  dat_tbl$W[t_ind] = dat_tbl %$% {W[t_ind - 1] * exp(-user_lambda__hr[t_ind] * meas_t__h[t_ind]) + (user_Rn_ss_flux__dpmm2hr[t_ind] - (W[t_ind-1] - (Q[t_ind] * V[t_ind])) * S[t_ind] * 60) * (1 - exp(-user_lambda__hr[t_ind] * meas_t__h[t_ind])) / user_lambda__hr[t_ind] / U[t_ind]}
}

dat_tbl <- dat_tbl %>%
  mutate(
    R = lag(W) - (Q*V),
    T = R * S * 60,
    Z = R * Y * 60,
    # there is an error in the xls: not using lambda in the next line:
    AB = c((lag(W) * exp(-user_lambda__hr * meas_t__h) * U)[1:2], rep(NA_real_, length(time) - 2)),
    AC = user_Rn_ss_flux__dpmm2hr * (1 - exp(-user_lambda__hr * meas_t__h)) / user_lambda__hr,
    AD = T * (1 - exp(-user_lambda__hr * meas_t__h)) / user_lambda__hr,
    AE = AB + (AC - AD),
  )

for(t_ind in seq_along(dat_tbl$time)[-(1:2)]) {
  dat_tbl$AB[t_ind] = dat_tbl %$% {AE[t_ind - 1] * exp(-user_lambda__hr[t_ind] * meas_t__h[t_ind])}
  dat_tbl$AE[t_ind] = dat_tbl %$% {AB[t_ind] + (AC[t_ind] - AD[t_ind])}
}

dat_tbl <- dat_tbl %>%
  mutate(
    AF = AE / user_depth__m,
    AG = Z * (1 - exp(-user_lambda__hr * meas_t__h)) / user_lambda__hr,
    AH = AB + (AC - AG),
    AI = AH / U,
    AJ = Q/Q[2],
    AK = c((AF * user_depth__m)[1:2], (AF * AJ)[-(1:2)]),
    AL = AI * AJ, 
    # AM = `Rn in Water dpm/m3`,
    AN = (AK - `Rn in Water dpm/m3`)^2,
    AO = (AL - `Rn in Water dpm/m3`)^2,
    
  ) %>% 
  
  # drop columns with no values (keep those with at least one value)
  select(where(~(!(.x %>% is.na())) %>% any()))

RMSE_vals <- dat_tbl %>% 
  summarize(
    RMSE_old = mean(AN, na.rm = TRUE) %>% sqrt(),
    RMSE_new = mean(AO, na.rm = TRUE) %>% sqrt(),
  )

# results saved in a csv file
write_csv(dat_tbl, here(study_folder, "output", "rn_budget.csv"), na = "")
  
# END OF RADON BUDGET CALCULATION

# plot
dat_tbl %>% 
  ggplot(mapping = aes(x = time)) + 
  geom_line(mapping = aes(y = `Rn in Water dpm/m3`), alpha = 0.6) +
  geom_point(mapping = aes(y = AL), color = "red") +
  labs(y = "Rn in water: actual (black) vs model (red)") +
  theme_light()
ggsave(here(study_folder, "output", "rn_budget.png"))
  
# *************************
#  end ----
# *************************
