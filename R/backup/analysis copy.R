# *************************
# Rn mass balance calculation
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

# input file name
excel_file_in <- "sgd_ts_data.xlsx"

# input sheet names
single_sheet_in <- 1
timser_sheet_in <- 2

# *************************
#  load data ----
# *************************

# load the variables containing single observations
single_in <- readxl::read_excel(here("input", excel_file_in), sheet = single_sheet_in, range = NULL) %>% 
  mutate(across(.cols = everything(), .fns = as.numeric))

# load the time series
timser_in <- readxl::read_excel(here("input", excel_file_in), sheet = timser_sheet_in, range = NULL) %>% 
  mutate(across(.cols = 1, .fns = ~force_tz(., tzone = ""))) %>%
  mutate(across(.cols = -1, .fns = as.numeric)) %>% 
  ts_long() %>% 
  ts_xts()

# *************************
#  calculations ----
# *************************

# coastal radon measurement time interval in minutes, all other time series parameters provided by the user should be averaged to this interval
ts_dat <- timser_in %>% 
  ts_c(
    timser_in %>% 
      ts_tbl() %>% 
      ts_wide() %>% 
      mutate(meas_t__min = time) %>% 
      select(time, meas_t__min) %>% 
      ts_long() %>% 
      ts_diff() %>% 
      mutate(value = value %>% as.numeric / 60)
  )

# water/air partitioning coefficient Kw/air; calculations according to Schubert et al., 2012
ts_dat$temp_wat__K <- timser_in$temp_wat__C + 273.15

ts_dat$kw_air <- 
  exp(-76.14 + 120.36 * (100 / ts_dat$temp_wat__K) + 31.26 * log(ts_dat$temp_wat__K / 100) + timser_in$sal_wat) *
  (-0.2631 + 0.1673 * (ts_dat$temp_wat__K / 100) + (-0.0273 * (ts_dat$temp_wat__K / 100)^2)) * 
  ts_dat$temp_wat__K / 273.15

# if(timser_in$Rn_exch__Bqm3) isnumber <- Rn_wat__Bqm3 = timser_in$Rn_exch__Bqm3 * kw/air
# otherwise Rn_wat__Bqm3 = timser_in$Rn_wat__Bqm3
# 
# If timser_in$wind__ms > 3.6
# f_atm__Bqm2hr =
#   ((0.45*(timser_in$wind__ms^1.6)*((0.0086/(10^(-((980/temp_wat__K+1.59)))/600)^-(1/2)))/100)/60)*
#   (Rn_wat__Bqm3=- kw/air* timser_in$Rn_air__Bqm3) * 60
# 
# If timser_in$wind__ms >1.5 and <= 3.6
# f_atm__Bqm2hr =
#   ((0.45*(timser_in$wind__ms^1.6)*((0.0086/(10^(-((980/temp_wat__K+1.59)))/600)^-(2/3)))/100)/60)*
#   (Rn_wat__Bqm3=- kw/air* timser_in$Rn_air__Bqm3) * 60
# 
# If timser_in$wind__ms <= 1.5
# f_atm__Bqm2hr =
#   ((0.45*(1.5^1.6)*((0.0086/(10^(-((980/temp_wat__K+1.59)))/600)^-(2/3)))/100)/60)*
#   (Rn_wat__Bqm3=- kw/air* timser_in$Rn_air__Bqm3) * 60
# 
# 
# f_dif__Bqm2hr = (495 * single_in$Ra226_sed__Bqg * 60 + 18.2) / 24
# 
# exRn_wat__Bqm3 = Rn_wat__Bqm3 - single_in$Ra226_wat__Bqm3
# 
# exRn_wat_inv__Bqm2 = exRn_wat__Bqm3 * timser_in$depth__m
# 
# f_Rn_gross__Bqm2hr = (ts_lag(exRn_wat_inv__Bqm2, 1) - exRn_wat_inv__Bqm2) * 60/meas_t__min
# 
# if ts_lag(timser_in$depth__m, 1) - timser_in$depth__m < 0
# f_Rn_flood__Bqm2hr = ((ts_lag(timser_in$depth__m, 1) - timser_in$depth__m) * single_in$Rn_offshore__Bqm3) * 60/meas_t__min 
# 
# if ts_lag(timser_in$depth__m, 1) - timser_in$depth__m > 0
# f_Rn_ebb__Bqm2hr = ((ts_lag(timser_in$depth__m, 1) - timser_in$depth__m) * exRn_wat__Bqm3) * 60/meas_t__min
# 
# f_Rn_net__Bqm2hr = f_Rn_gross__Bqm2hr + f_atm__Bqm2hr + f_Rn_ebb__Bqm2hr + f_Rn_flood__Bqm2hr - f_dif__Bqm2hr + f_mix_exp__Bqm2hr
# 
# if f_Rn_net__Bqm2hr <0 then f_mix__Bqm2hr = - f_Rn_net__Bqm2hr otherwise f_mix__Bqm2hr = 0
# 
# f_Rn_total__Bqm2hr = f_Rn_net__Bqm2hr + f_mix__Bqm2hr
# 
# f_gw__m3m2d = (f_Rn_total__Bqm2hr / single_in$Rn_gw__Bqm3) * 24
# 
# 
# 
