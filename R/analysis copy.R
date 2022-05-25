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
in_xts <- read_csv(here("input", csv_file_in)) %>%
  mutate(across(.cols = 1, .fns = ~ force_tz(., tzone = ""))) %>%
  mutate(across(.cols = -1, .fns = as.numeric)) %>%
  ts_long() %>%
  ts_xts()

# *************************
#  calculations ----
# *************************

# add coastal radon measurement time interval in minutes, 
# all other time series parameters provided by the user should be averaged to this interval
dat_xts <- in_xts %>%
  ts_c(
    in_xts %>%
      ts_tbl() %>%
      ts_wide() %>%
      mutate(meas_t__min = time) %>%
      select(time, meas_t__min) %>%
      ts_long() %>%
      ts_diff() %>%
      mutate(value = value %>% as.numeric() / 60)
  )

# water/air partitioning coefficient Kw/air; calculations according to Schubert et al., 2012
dat_xts$temp_wat__K <- dat_xts$temp_wat__C + 273.15

# explain
dat_xts$kw_air <-
  exp(-76.14 + 120.36 * (100 / dat_xts$temp_wat__K) + 31.26 * log(dat_xts$temp_wat__K / 100) + dat_xts$sal_wat) *
    (-0.2631 + 0.1673 * (dat_xts$temp_wat__K / 100) + (-0.0273 * (dat_xts$temp_wat__K / 100)^2)) *
    dat_xts$temp_wat__K / 273.15

# explain
if (!(in_xts$Rn_exch__Bqm3 %>% is.null()) & !(in_xts$Rn_exch__Bqm3 %>% is.na()) %>% sum()) {
  dat_xts$Rn_wat__Bqm3 <- dat_xts$Rn_exch__Bqm3 * dat_xts$kw_air
}

# explain
dat_xts <- dat_xts %>%
  ts_tbl() %>%
  ts_wide() %>%
  mutate(
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
      )
  ) %>%
  ts_long() %>%
  ts_xts()

# explain
dat_xts$f_dif__Bqm2hr <- (495 * dat_xts$Ra226_sed__Bqg * 60 + 18.2) / 24

# explain
dat_xts$ex_Rn_wat__Bqm3 <- dat_xts$Rn_wat__Bqm3 - dat_xts$Ra226_wat__Bqm3

# explain
dat_xts$ex_Rn_wat_inv__Bqm2 <- dat_xts$ex_Rn_wat__Bqm3 * dat_xts$depth__m

# explain
dat_xts$f_Rn_gross__Bqm2hr <- ts_diff(dat_xts$ex_Rn_wat_inv__Bqm2) * 60 / dat_xts$meas_t__min

# explain
dat_xts$f_Rn_flood__Bqm2hr <- ifelse(ts_diff(dat_xts$depth__m) > 0,
  (ts_diff(dat_xts$depth__m) * dat_xts$Rn_offshore__Bqm3) * 60 / dat_xts$meas_t__min,
  NA
)

# explain
dat_xts$f_Rn_ebb__Bqm2hr <- ifelse(ts_diff(dat_xts$depth__m) < 0,
  (ts_diff(dat_xts$depth__m) * dat_xts$ex_Rn_wat__Bqm3) * 60 / dat_xts$meas_t__min,
  NA
)

# explain
dat_xts$f_Rn_net__Bqm2hr <- dat_xts$f_Rn_gross__Bqm2hr +
  dat_xts$f_atm__Bqm2hr +
  dat_xts$f_Rn_ebb__Bqm2hr +
  dat_xts$f_Rn_flood__Bqm2hr -
  dat_xts$f_dif__Bqm2hr +
  dat_xts$f_mix_exp__Bqm2hr

# explain
dat_xts$f_mix__Bqm2hr <- ifelse(dat_xts$f_Rn_net__Bqm2hr < 0,
  -dat_xts$f_Rn_net__Bqm2hr,
  0
)

# explain
dat_xts$f_Rn_total__Bqm2hr <- dat_xts$f_Rn_net__Bqm2hr + dat_xts$f_mix__Bqm2hr

# explain
dat_xts$f_gw__m3m2d <- (dat_xts$f_Rn_total__Bqm2hr / dat_xts$Rn_gw__Bqm3) * 24
