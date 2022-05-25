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
csv_file_in <- "sgd_ts_data.csv"

# *************************
#  load data ----
# *************************

# load the time series
in_tbl <- read_csv(here("input", csv_file_in)) %>%
  mutate(across(.cols = 1, .fns = ~ force_tz(., tzone = ""))) %>%
  mutate(across(.cols = -1, .fns = as.numeric))

# *************************
#  calculations ----
# *************************

# explain
dat_tbl <- in_tbl %>%
  mutate(
    # add coastal radon measurement time interval in minutes,
    # all other time series parameters provided by the user should be averaged to this interval
    meas_t__min = (time %>% as.numeric() - lag(time %>% as.numeric())) / 60,

    # water/air partitioning coefficient Kw/air; calculations according to Schubert et al., 2012
    temp_wat__K = temp_wat__C + 273.15,

    # explain
    kw_air =
      exp(-76.14 + 120.36 * (100 / temp_wat__K) + 31.26 * log(temp_wat__K / 100) + sal_wat) *
        (-0.2631 + 0.1673 * (temp_wat__K / 100) + (-0.0273 * (temp_wat__K / 100)^2)) *
        temp_wat__K / 273.15,

    # explain
    Rn_wat__Bqm3 = if (!(Rn_exch__Bqm3 %>% is.null()) & !(Rn_exch__Bqm3 %>% is.na()) %>% any()) {
      Rn_exch__Bqm3 * kw_air
    } else {
      Rn_wat__Bqm3
    },

    # explain
    f_atm__Bqm2hr =
      case_when(
        wind__ms > 3.6 ~
          ((0.45 * (wind__ms^1.6) *
            ((0.0086 / (10^(-((980 / temp_wat__K + 1.59))) / 600)^(-(1 / 2)))) / 100) / 60) *
            (Rn_wat__Bqm3 - kw_air * Rn_air__Bqm3) * 60,
        wind__ms > 1.5 ~
          ((0.45 * (wind__ms^1.6) *
            ((0.0086 / (10^(-((980 / temp_wat__K + 1.59))) / 600)^(-(2 / 3)))) / 100) / 60) *
            (Rn_wat__Bqm3 - kw_air * Rn_air__Bqm3) * 60,
        TRUE ~
          ((0.45 * (1.5^1.6) *
            ((0.0086 / (10^(-((980 / temp_wat__K + 1.59))) / 600)^(-(2 / 3)))) / 100) / 60) *
            (Rn_wat__Bqm3 - kw_air * Rn_air__Bqm3) * 60
      ),

    # explain
    f_dif__Bqm2hr = (495 * Ra226_sed__Bqg * 60 + 18.2) / 24,

    # explain
    ex_Rn_wat__Bqm3 = Rn_wat__Bqm3 - Ra226_wat__Bqm3,

    # explain
    ex_Rn_wat_inv__Bqm2 = ex_Rn_wat__Bqm3 * depth__m,

    # explain
    f_Rn_gross__Bqm2hr = (ex_Rn_wat_inv__Bqm2 - lag(ex_Rn_wat_inv__Bqm2)) * 60 / meas_t__min,

    # explain
    f_Rn_flood__Bqm2hr = if_else(depth__m - lag(depth__m) > 0,
      ((depth__m - lag(depth__m)) * Rn_offshore__Bqm3) * 60 / meas_t__min,
      NA_real_,
      NA_real_
    ),

    # explain
    f_Rn_ebb__Bqm2hr = if_else(depth__m - lag(depth__m) < 0,
      ((depth__m - lag(depth__m)) * ex_Rn_wat__Bqm3) * 60 / meas_t__min,
      NA_real_,
      NA_real_
    ),

    # explain
    f_Rn_net__Bqm2hr = f_Rn_gross__Bqm2hr +
      f_atm__Bqm2hr +
      f_Rn_ebb__Bqm2hr +
      f_Rn_flood__Bqm2hr -
      f_dif__Bqm2hr +
      f_mix_exp__Bqm2hr,

    # explain
    f_mix__Bqm2hr <- if_else(f_Rn_net__Bqm2hr < 0,
      -f_Rn_net__Bqm2hr,
      0,
      NA_real_
    ),

    # explain
    f_Rn_total__Bqm2hr = f_Rn_net__Bqm2hr + f_mix__Bqm2hr,

    # explain
    f_gw__m3m2d = (f_Rn_total__Bqm2hr / Rn_gw__Bqm3) * 24
  )
