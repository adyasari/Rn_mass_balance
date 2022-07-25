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

# names(in_tbl)

# *************************
#  calculations ----
# *************************

# data manipulations and calculations
Rn_ss_calc <- function(in_tbl, optim_driver = NULL){
  #input: in_tbl - observed (user supplied) data
  #input: optim_driver - initial value for Rn flux via SGD
  
  # create data table
  dat_tbl <- in_tbl %>%
    mutate(
      
      # calculates radon measurement time interval in hours based on provided measurement times, 
      # all other time series measurements provided by the user should be averaged to this interval
      time = `Rad Mid Time`,
      meas_t__h = (time %>% as.numeric() - lag(time %>% as.numeric())) / 60 / 60,
      # tot_t__h = slider::slide_dbl(meas_t__h, ~sum(.x , na.rm = TRUE), .before = Inf),
      # tot_t__day = tot_t__h / 24,
      wind__ms = `Wind speed (m/s)`,
      temp_wat__C = lag(`Water Temperature C`),
      temp_wat__K = temp_wat__C + 273.15,
      mole_diff_Rn__cm2s = 10^(-(980 / (273 + temp_wat__C) + 1.59)),
      dens_wat__gcm3 = (999.842594 + (0.06793952) * temp_wat__C - (0.00909529) * temp_wat__C^2 + (0.0001001685) * temp_wat__C^3 - (0.000001120083) * temp_wat__C^4 + (0.000000006536332) * temp_wat__C^5) / 1000,
      abs_visc_wat__gscm = 0.00002414 * 10^(247.8 / (temp_wat__K - 140)) * 1000 / 100,
      kin_visc_wat__cm2s = abs_visc_wat__gscm / dens_wat__gcm3,
      schm_num = kin_visc_wat__cm2s / mole_diff_Rn__cm2s,
      k_piston_Rn__mmin = ((0.45 * (wind__ms^1.6) * ((schm_num / 600)^-0.667)) / 100) / 60,
      kw_air = 0.105 + 0.405 * exp(-0.05027 * temp_wat__C),
      Rn_air__Bqm3 = lag(`Rn in Air dpm/m3`),
      depth__m = lag(`Water Level m`),
      # user_Rn_atm__dpmm3 = 800,
      # user_depth__m = 1.3,
      # decay constant of Rn in hours
      lambda__hr = log(2) / (3.84 * 24),
      user_Rn_ss__dpmm3 = if (is.null(optim_driver)) user_Rn_ss__dpmm3 else optim_driver, # driver of optimization
      user_Rn_ini_inv__dpmm2 = user_depth__m * user_Rn_ss__dpmm3,
      user_Rn_ss_flux__dpmm2hr = lambda__hr * user_Rn_ini_inv__dpmm2, # calculate from components
      user_Rn_ss_inv__dpmm2 = user_Rn_ss_flux__dpmm2hr / lambda__hr,
      # user_kin_visc__cm2s = 0.01004,
      # user_mol_diff__cm2s = 0.000012,
      # X = user_kin_visc__cm2s / kin_visc_wat__cm2s,
      # Y = (34.6 * X * (user_mol_diff__cm2s^0.5)*(wind__ms^1.5)) / (24 * 60),
      mc_Rn_water__dpmm3 = c(`Rn in Water dpm/m3`[2], rep(NA_real_, length(time) - 1)), #`Rn in Water dpm/m3`,
    )
  
  for(t_ind in seq_along(dat_tbl$time)[-1]) {
    dat_tbl$mc_Rn_water__dpmm3[t_ind] = dat_tbl %$% {mc_Rn_water__dpmm3[t_ind - 1] * exp(-lambda__hr[t_ind] * meas_t__h[t_ind]) + (user_Rn_ss_flux__dpmm2hr[t_ind] - (mc_Rn_water__dpmm3[t_ind-1] - (kw_air[t_ind] * Rn_air__Bqm3[t_ind])) * k_piston_Rn__mmin[t_ind] * 60) * (1 - exp(-lambda__hr[t_ind] * meas_t__h[t_ind])) / lambda__hr[t_ind] / depth__m[t_ind]}
  }
  
  dat_tbl <- dat_tbl %>%
    mutate(
      mc_Rn_air_exch__dpmm3 = lag(mc_Rn_water__dpmm3),
      grad_Rn_conc__dpmm3 = mc_Rn_air_exch__dpmm3 - (kw_air*Rn_air__Bqm3),
      f_Rn_atm__Bqm2hr = grad_Rn_conc__dpmm3 * k_piston_Rn__mmin * 60,
      # Z = grad_Rn_conc__dpmm3 * Y * 60,
      # there is an error in the xls: not using lambda in the next line:
      m_inv_Rn_decay_corr__dpmm2 = c((lag(mc_Rn_water__dpmm3) * exp(-lambda__hr * meas_t__h) * depth__m)[1:2], rep(NA_real_, length(time) - 2)),
      m_SGD_inv_Rn__dpmm2 = user_Rn_ss_flux__dpmm2hr * (1 - exp(-lambda__hr * meas_t__h)) / lambda__hr,
      m_air_loss_inv_Rn__dpmm2 = f_Rn_atm__Bqm2hr * (1 - exp(-lambda__hr * meas_t__h)) / lambda__hr,
      m_net_inv_Rn__dpmm2 = m_inv_Rn_decay_corr__dpmm2 + (m_SGD_inv_Rn__dpmm2 - m_air_loss_inv_Rn__dpmm2),
    )
  
  for(t_ind in seq_along(dat_tbl$time)[-(1:2)]) {
    dat_tbl$m_inv_Rn_decay_corr__dpmm2[t_ind] = dat_tbl %$% {m_net_inv_Rn__dpmm2[t_ind - 1] * exp(-lambda__hr[t_ind] * meas_t__h[t_ind])}
    dat_tbl$m_net_inv_Rn__dpmm2[t_ind] = dat_tbl %$% {m_inv_Rn_decay_corr__dpmm2[t_ind] + (m_SGD_inv_Rn__dpmm2[t_ind] - m_air_loss_inv_Rn__dpmm2[t_ind])}
  }
  
  dat_tbl <- dat_tbl %>%
    mutate(
      final_Rn_wat_no_temp_corr__dpmm3 = m_net_inv_Rn__dpmm2 / depth__m,
      # AG = Z * (1 - exp(-lambda__hr * meas_t__h)) / lambda__hr,
      # AH = m_inv_Rn_decay_corr__dpmm2 + (m_SGD_inv_Rn__dpmm2 - AG),
      # AI = AH / depth__m,
      R_tem_sol_corr = kw_air/kw_air[2],
      # final_mc_Rn__dpmm3 = c((final_Rn_wat_no_temp_corr__dpmm3 * user_depth__m)[1:2], (final_Rn_wat_no_temp_corr__dpmm3 * R_tem_sol_corr)[-(1:2)]),
      final_mc_Rn__dpmm3 = final_Rn_wat_no_temp_corr__dpmm3 * R_tem_sol_corr,
      # AL = AI * R_tem_sol_corr, 
      # AM = `Rn in Water dpm/m3`,
      sq_err = (final_mc_Rn__dpmm3 - `Rn in Water dpm/m3`)^2,
      # sq_err = (AL - `Rn in Water dpm/m3`)^2,
      
    ) %>% 
    
    # drop columns with no values (keep those with at least one value)
    select(where(~(!(.x %>% is.na())) %>% any()))
  
  RMSE_vals <- dat_tbl %>% 
    summarize(
      # RMSE_old = mean(sq_err, na.rm = TRUE) %>% sqrt(),
      RMSE_new = mean(sq_err, na.rm = TRUE) %>% sqrt(),
    ) %>% 
    pull(RMSE_new)
  
  return(list(dat = dat_tbl, rmse = RMSE_vals))
  
}

# data manipulations and calculations
dat_out <- Rn_ss_calc(in_tbl)
dat_tbl <- dat_out$dat

# results saved in a csv file
write_csv(dat_tbl, here(study_folder, "output", "rn_budgetND.csv"), na = "")

# *************************
#  optimization ----
# *************************

# objective function to find steady state Rn flux
Rn_ss_obj_fun <- function(in_tbl, optim_driver){
  Rn_ss_calc(in_tbl, optim_driver)$rmse
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
opt1 <- optim(optim_start, Rn_ss_obj_fun, in_tbl = in_tbl, method = "BFGS")
opt2 <- optimize(Rn_ss_obj_fun, in_tbl = in_tbl, c(0, 50000))
opt1$par - opt2$minimum

# *************************
#  calcs with optimal value ----
# *************************

optim_driver <- opt1$par

# data manipulations and calculations
dat_out <- Rn_ss_calc(in_tbl, optim_driver)
dat_tbl <- dat_out$dat

# results saved in a csv file
write_csv(dat_tbl, here(study_folder, "output", "rn_budget.csv"), na = "")
  
# END OF RADON BUDGET CALCULATION

# plot
dat_tbl %>% 
  ggplot(mapping = aes(x = time)) + 
  geom_line(mapping = aes(y = `Rn in Water dpm/m3`), alpha = 0.6) +
  geom_point(mapping = aes(y = final_mc_Rn__dpmm3), color = "red") +
  labs(y = "Rn in water: actual (black) vs model (red)") +
  theme_light()
ggsave(here(study_folder, "output", "rn_budget.png"))
  
# *************************
#  end ----
# *************************
