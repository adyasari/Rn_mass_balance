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
  mutate(across(.cols = 1, .fns = ~ force_tz(., tzone = ""))) %>%
  mutate(across(.cols = -1, .fns = as.numeric)) %>%
  ts_long() %>%
  ts_xts()

# *************************
#  calculations ----
# *************************

# add coastal radon measurement time interval in minutes, all other time series parameters provided by the user should be averaged to this interval
ts_dat <- timser_in %>%
  ts_c(
    timser_in %>%
      ts_tbl() %>%
      ts_wide() %>%
      mutate(meas_t__min = time) %>%
      select(time, meas_t__min) %>%
      ts_long() %>%
      ts_diff() %>%
      mutate(value = value %>% as.numeric() / 60)
  )

# alternative to combine all data into single xts
# temp <- ts_dat %>% ts_tbl() %>% ts_wide() %>% bind_cols(single_in) %>% ts_long() %>% ts_xts()

# water/air partitioning coefficient Kw/air; calculations according to Schubert et al., 2012
ts_dat$temp_wat__K <- ts_dat$temp_wat__C + 273.15

# explain
ts_dat$kw_air <-
  exp(-76.14 + 120.36 * (100 / ts_dat$temp_wat__K) + 31.26 * log(ts_dat$temp_wat__K / 100) + ts_dat$sal_wat) *
    (-0.2631 + 0.1673 * (ts_dat$temp_wat__K / 100) + (-0.0273 * (ts_dat$temp_wat__K / 100)^2)) *
    ts_dat$temp_wat__K / 273.15

# explain
if (!(timser_in$Rn_exch__Bqm3 %>% is.null()) & !(timser_in$Rn_exch__Bqm3 %>% is.na()) %>% sum()) {
  ts_dat$Rn_wat__Bqm3 <- ts_dat$Rn_exch__Bqm3 * ts_dat$kw_air
} # else {
#   ts_dat$Rn_wat__Bqm3 <- timser_in$Rn_wat__Bqm3
# }

# if timser_in$Rn_air__Bqm3 is missing altogether or an empty column, then populate with single_in$Rn_air__Bqm3
if (timser_in$Rn_air__Bqm3 %>% is.null() | !(!(timser_in$Rn_air__Bqm3 %>% is.na()) %>% sum())) {
  ts_dat$Rn_air__Bqm3 <- single_in$Rn_air__Bqm3
} # else {
#   ts_dat$Rn_air__Bqm3 <- timser_in$Rn_air__Bqm3
# }

# explain
ts_dat <- ts_dat %>%
  ts_tbl() %>%
  ts_wide() %>%
  mutate(
    f_atm__Bqm2hr =
      case_when(
        wind__ms > 3.6 ~
          ((0.45 * (wind__ms^1.6) * ((0.0086 / (10^(-((980 / temp_wat__K + 1.59))) / 600)^(-(1 / 2)))) / 100) / 60) *
            (Rn_wat__Bqm3 - kw_air * Rn_air__Bqm3) * 60,
        wind__ms > 1.5 ~
          ((0.45 * (wind__ms^1.6) * ((0.0086 / (10^(-((980 / temp_wat__K + 1.59))) / 600)^(-(2 / 3)))) / 100) / 60) *
            (Rn_wat__Bqm3 - kw_air * Rn_air__Bqm3) * 60,
        TRUE ~
          ((0.45 * (1.5^1.6) *      ((0.0086 / (10^(-((980 / temp_wat__K + 1.59))) / 600)^(-(2 / 3)))) / 100) / 60) *
            (Rn_wat__Bqm3 - kw_air * Rn_air__Bqm3) * 60
      )
  ) %>%
  ts_long() %>%
  ts_xts()

# explain
ts_dat$f_dif__Bqm2hr <- (495 * single_in$Ra226_sed__Bqg * 60 + 18.2) / 24

# explain
ts_dat$ex_Rn_wat__Bqm3 <- ts_dat$Rn_wat__Bqm3 - single_in$Ra226_wat__Bqm3

# explain
ts_dat$ex_Rn_wat_inv__Bqm2 <- ts_dat$ex_Rn_wat__Bqm3 * ts_dat$depth__m

# explain
ts_dat$f_Rn_gross__Bqm2hr <- ts_diff(ts_dat$ex_Rn_wat_inv__Bqm2) * 60 / ts_dat$meas_t__min

# PF edited through here (check ts_dat: missing values in Rn_wat__Bqm3 are propagated through the calcs.)

# explain
ts_dat$f_Rn_flood__Bqm2hr <- ifelse(ts_diff(ts_dat$depth__m) > 0, (ts_diff(ts_dat$depth__m) * single_in$Rn_offshore__Bqm3) * 60 / ts_dat$meas_t__min, NA)

# explain
ts_dat$f_Rn_ebb__Bqm2hr <- ifelse(ts_diff(ts_dat$depth__m) < 0, (ts_diff(ts_dat$depth__m) * ts_dat$ex_Rn_wat__Bqm3) * 60 / ts_dat$meas_t__min, NA)

# explain
ts_dat$f_Rn_net__Bqm2hr <- ts_dat$f_Rn_gross__Bqm2hr + ts_dat$f_atm__Bqm2hr + ts_dat$f_Rn_ebb__Bqm2hr + ts_dat$f_Rn_flood__Bqm2hr - ts_dat$f_dif__Bqm2hr + ts_dat$f_mix_exp__Bqm2hr

# explain
ts_dat$f_mix__Bqm2hr <- ifelse(ts_dat$f_Rn_net__Bqm2hr < 0, -ts_dat$f_Rn_net__Bqm2hr, 0)

# explain
ts_dat$f_Rn_total__Bqm2hr <- ts_dat$f_Rn_net__Bqm2hr + ts_dat$f_mix__Bqm2hr

# explain
ts_dat$f_gw__m3m2d <- (ts_dat$f_Rn_total__Bqm2hr / single_in$Rn_gw__Bqm3) * 24



