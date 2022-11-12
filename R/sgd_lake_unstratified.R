# *************************
# Rn mass balance calculation for time series lake radon budget based on Dimova and Burnett, 2011
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

# if folder hierarchy segmented by study type, specify it here
study_folder <- "."

# input file name
csv_file_in <- "sgd_lake_unstratified_data.csv"

# *************************
#  load data ----
# *************************

# load lake time series
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
      meas_t__h = (time %>% as.numeric() - lag(time %>% as.numeric())) / 60 / 60,
      
      # wind__ms = wind__ms, # not lagged
      temp_wat__C = lag(temp_wat__C), # lagged WHY??? parameters are taken from previous measurement step to account for radon measurement delay???
      temp_wat__K = temp_wat__C + 273.15,
      # mole_diff_Rn__cm2s = 10^(-(980 / (273 + temp_wat__C) + 1.59)),
      mole_diff_Rn__m2s = 10^(-(980 / (273 + temp_wat__C) + 1.59))/10000,
      # dens_wat__gcm3 = (999.842594 + (0.06793952) * temp_wat__C - (0.00909529) * temp_wat__C^2 + (0.0001001685) * temp_wat__C^3 - (0.000001120083) * temp_wat__C^4 + (0.000000006536332) * temp_wat__C^5) / 1000,
      dens_wat__kgm3 = (999.842594 + (0.06793952) * temp_wat__C - (0.00909529) * temp_wat__C^2 + (0.0001001685) * temp_wat__C^3 - (0.000001120083) * temp_wat__C^4 + (0.000000006536332) * temp_wat__C^5),
      # abs_visc_wat__gscm = 0.00002414 * 10^(247.8 / (temp_wat__K - 140)) * 1000 / 100,
      abs_visc_wat__kgsm = 0.00002414 * 10^(247.8 / (temp_wat__K - 140)),
      kin_visc_wat__m2s = abs_visc_wat__kgsm / dens_wat__kgm3,
      schm_num = kin_visc_wat__m2s / mole_diff_Rn__m2s,
      k_piston_Rn__ms = ((0.45 * (wind__ms^1.6) * ((schm_num / 600)^-0.667)) / 100) / 3600,
      kw_air = 0.105 + 0.405 * exp(-0.05027 * temp_wat__C),
      
      # if Rad-Aqua was used to collect radon data and radon in the exchanger (therefore in air) is provided 
      # it is converted to Rn in water in this step;
      # otherwise, the provided radon in water is used for further calculations
      # the condition checks if a Rn_wat__Bqm3 column exists and if it is non-empty
      Rn_wat__Bqm3 = if (("Rn_wat__Bqm3" %in% names(.)) && (!(Rn_wat__Bqm3 %>% is.na())) %>% any()) {
        Rn_wat__Bqm3
      } else {
        Rn_exch__Bqm3 * kw_air
      },
      
      Rn_air__Bqm3 = lag(Rn_air__Bqm3), # lagged WHY???
      depth__m = lag(depth__m), # lagged WHY???
      
      # decay constant of Rn in hours
      lambda__s = log(2) / (3.84 * 24 * 3600),
      
      # if an optimal value or the driver of the optimization process is not provided, 
      # then use the sample mean of Rn_wat__Bqm3
      Rn_ss__Bqm3 = if (is.null(optim_driver)) mean(Rn_wat__Bqm3, na.rm = TRUE) else optim_driver, # driver of optimization
      
      Rn_ini_inv__Bqm2 = depth__m * Rn_ss__Bqm3,
      Rn_ss_flux__Bqm2s = lambda__s * Rn_ini_inv__Bqm2,
      mc_Rn_water__Bqm3 = c(Rn_wat__Bqm3[2], rep(NA_real_, length(time) - 1)), # if first row deleted then start from 1 rather than 2
    )
  
  # fill in mc_Rn_water__Bqm3 using this recursion (mc_Rn_water__Bqm3 values depend on their own lags)
  for(t_ind in seq_along(dat_tbl$time)[-1]) {
    dat_tbl$mc_Rn_water__Bqm3[t_ind] = dat_tbl %$% {mc_Rn_water__Bqm3[t_ind - 1] * exp(-lambda__s[t_ind] * meas_t__h[t_ind]) + (Rn_ss_flux__Bqm2s[t_ind] - (mc_Rn_water__Bqm3[t_ind-1] - (kw_air[t_ind] * Rn_air__Bqm3[t_ind])) * k_piston_Rn__ms[t_ind]) * (1 - exp(-lambda__s[t_ind] * meas_t__h[t_ind])) / lambda__s[t_ind] / depth__m[t_ind]}
  }
  
  # fill in remaining values (they depend on mc_Rn_water__Bqm3)
  dat_tbl <- dat_tbl %>%
    mutate(
      mc_Rn_air_exch__Bqm3 = lag(mc_Rn_water__Bqm3),
      grad_Rn_conc__Bqm3 = mc_Rn_air_exch__Bqm3 - (kw_air*Rn_air__Bqm3),
      f_Rn_atm__Bqm2s = grad_Rn_conc__Bqm3 * k_piston_Rn__ms,
      m_inv_Rn_decay_corr__Bqm2 = c((lag(mc_Rn_water__Bqm3) * exp(-lambda__s * meas_t__h) * depth__m)[1:2], rep(NA_real_, length(time) - 2)),  # if first row deleted then start from 1 rather than 2
      m_SGD_inv_Rn__Bqm2 = Rn_ss_flux__Bqm2s * (1 - exp(-lambda__s * meas_t__h)) / lambda__s,
      m_air_loss_inv_Rn__Bqm2 = f_Rn_atm__Bqm2s * (1 - exp(-lambda__s * meas_t__h)) / lambda__s,
      m_net_inv_Rn__Bqm2 = m_inv_Rn_decay_corr__Bqm2 + (m_SGD_inv_Rn__Bqm2 - m_air_loss_inv_Rn__Bqm2),
    )
  
  # fill in m_inv_Rn_decay_corr__Bqm2 using this recursion (m_inv_Rn_decay_corr__Bqm2 values depend on the lag of m_net_inv_Rn__Bqm2)
  for(t_ind in seq_along(dat_tbl$time)[-(1:2)]) { # if first row deleted then start from 1 rather than 2
    dat_tbl$m_inv_Rn_decay_corr__Bqm2[t_ind] = dat_tbl %$% {m_net_inv_Rn__Bqm2[t_ind - 1] * exp(-lambda__s[t_ind] * meas_t__h[t_ind])}
    dat_tbl$m_net_inv_Rn__Bqm2[t_ind] = dat_tbl %$% {m_inv_Rn_decay_corr__Bqm2[t_ind] + (m_SGD_inv_Rn__Bqm2[t_ind] - m_air_loss_inv_Rn__Bqm2[t_ind])}
  }
  
  # fill in remaining series
  dat_tbl <- dat_tbl %>%
    mutate(
      final_Rn_wat_no_temp_corr__Bqm3 = m_net_inv_Rn__Bqm2 / depth__m,
      R_tem_sol_corr = kw_air/kw_air[2],  # if first row deleted then start from 1 rather than 2
      final_mc_Rn__Bqm3 = final_Rn_wat_no_temp_corr__Bqm3 * R_tem_sol_corr,
      sq_err = (final_mc_Rn__Bqm3 - Rn_wat__Bqm3)^2,

    ) %>% 
    
    # drop columns with no values (keep those with at least one value)
    select(where(~(!(.x %>% is.na())) %>% any()))
  
  # calculate the root mean squared error
  RMSE <- dat_tbl %>% 
    summarize(
      RMSE = mean(sq_err, na.rm = TRUE) %>% sqrt(),
    ) %>% 
    pull(RMSE)
  
  return(list(dat = dat_tbl, rmse = RMSE))
  
}

# data manipulations and calculations
dat_out <- Rn_ss_calc(in_tbl)
dat_tbl <- dat_out$dat

# intermediate results saved in a csv file
# write_csv(dat_tbl, here(study_folder, "output", "sgd_lake_unstratified_rn_budget.csv"), na = "")

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
  summarize(mean(Rn_ss__Bqm3, na.rm = TRUE)) %>% 
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
write_csv(dat_tbl, here(study_folder, "output", "sgd_lake_unstratified_rn_budget.csv"), na = "")
  
# END OF RADON BUDGET CALCULATION

# plot
dat_tbl %>% 
  ggplot(mapping = aes(x = time)) + 
  geom_line(mapping = aes(y = Rn_wat__Bqm3), alpha = 0.6) +
  geom_point(mapping = aes(y = final_mc_Rn__Bqm3), color = "red") +
  labs(y = "Rn in water: actual (black) vs model (red)") +
  theme_light()
ggsave(here(study_folder, "output", "sgd_lake_unstratified_rn_budget.png"))
  
# *************************
#  end ----
# *************************

