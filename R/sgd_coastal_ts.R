# *************************
# Rn mass balance calculation for time series coastal radon budget
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
csv_file_in <- "sgd_coastal_ts_data.csv"

# *************************
#  load data ----
# *************************

# load the time series
in_tbl <- read_csv(here(study_folder, "input", csv_file_in)) %>%
  # mutate(across(.cols = 1, .fns = ~ force_tz(., tzone = ""))) %>%
  mutate(across(.cols = 1, .fns = ~ lubridate::parse_date_time(., c("ymdHMS", "ymdHM", "mdyHM", "mdyHMS")))) %>%
  mutate(across(.cols = -1, .fns = as.numeric))

# *************************
#  calculations ----
# *************************

# creates data table
dat_tbl <- in_tbl %>%
  mutate(
    
    # the loaded data should have a fixed periodicity
    # interpolate linearly when a single value is missing from the time series
    # (e.g. due to issues with the measuring device)
    Rn_wat__Bqm3 = if_else(is.na(Rn_wat__Bqm3), (lag(Rn_wat__Bqm3) + lead(Rn_wat__Bqm3))/2, Rn_wat__Bqm3),
    temp_wat__C = if_else(is.na(temp_wat__C), (lag(temp_wat__C) + lead(temp_wat__C))/2, temp_wat__C),
    sal_wat = if_else(is.na(sal_wat), (lag(sal_wat) + lead(sal_wat))/2, sal_wat),
    
    # calculates coastal radon measurement time interval in minutes based on provided measurement times, 
    # all other time series parameters provided by the user should be averaged to this interval
    meas_t__min = (time %>% as.numeric() - lag(time %>% as.numeric())) / 60,
    
    # change in water depth between two time stamps
    diff_owl__m = depth__m - lag(depth__m),

    # if radon mixing losses are determined via an independent method (current measurements, residence time estimates)
    # then f_mix_exp__Bqm2hr should be provided in the spreadsheet
    # otherwise, mixing losses will be estimated from the Rn mass balance, see definition of f_Rn_mix__Bqm2hr below
    # the condition checks if a f_mix_exp__Bqm2hr column exists and if it is non-empty
    f_mix_exp__Bqm2hr = if (("f_mix_exp__Bqm2hr" %in% names(.)) && (!(f_mix_exp__Bqm2hr %>% is.na())) %>% any()) {
      f_mix_exp__Bqm2hr %>% if_else(is.na(.), 0, .)
    } else {
      0
    },
    
    # water temperature converted from degrees Celsius to Kelvin
    temp_wat__K = temp_wat__C + 273.15,

    # water/air partitioning coefficient kw_air based on water temperature and salinity; 
    # calculations and coefficients from Schubert et al. 2012
    kw_air =
      exp(-76.14 + 120.36 * (100 / temp_wat__K) + 31.26 * log(temp_wat__K / 100) + sal_wat *
        (-0.2631 + 0.1673 * (temp_wat__K / 100) + (-0.0270 * (temp_wat__K / 100)^2))) *
        temp_wat__K / 273.15,

    # if Rad-Aqua was used to collect radon data and radon in the exchanger (therfore in air) is provided 
    # it is converted to Rn in water in this step;
    # otherwise, the provided radon in water is used for further calculations
    # the condition checks if a Rn_wat__Bqm3 column exists and if it is non-empty
    Rn_wat__Bqm3 = if (("Rn_wat__Bqm3" %in% names(.)) && (!(Rn_wat__Bqm3 %>% is.na())) %>% any()) {
      Rn_wat__Bqm3
    } else {
      Rn_exch__Bqm3 * kw_air
    },
    
    # Rn losses by evasion into the atmosphere are calculated according to MacIntyre et al. (1995)
    # for wind speeds above 3.6 m/s Sc^-1/2 and for wind speeds below 3.6 m/s Sc^-2/3 is applied (Turner et al., 1996);
    # for wind speeds below 1.5 m/s k is assumed to be constant and equivalent to wind speeds of 1.5 m/s (Ocampo-torres et al., 1994)
    # note that kinematic viscosity is considered constant, one can calculate more accurate values that account for salinity and temperature
    f_Rn_atm__Bqm2hr =
      case_when(
        wind__ms > 3.6 ~
          (0.45 * wind__ms^1.6 *
            ((0.0086 / 10^(-(980 / temp_wat__K + 1.59))) / 600)^(-1 / 2)) / 100 / 60 *
            (Rn_wat__Bqm3 - kw_air * Rn_air__Bqm3) * 60,
        wind__ms > 1.5 ~
          (0.45 * wind__ms^1.6 *
             ((0.0086 / 10^(-(980 / temp_wat__K + 1.59))) / 600)^(-2 / 3)) / 100 / 60 *
            (Rn_wat__Bqm3 - kw_air * Rn_air__Bqm3) * 60,
        TRUE ~
          (0.45 * 1.5^1.6 *
             ((0.0086 / 10^(-(980 / temp_wat__K + 1.59))) / 600)^(-2 / 3)) / 100 / 60 *
            (Rn_wat__Bqm3 - kw_air * Rn_air__Bqm3) * 60
      ),

    # excess Rn in water is calculated by subtracting dissolved 226Ra in water
    ex_Rn_wat__Bqm3 = Rn_wat__Bqm3 - Ra226_wat__Bqm3,

    # excess Rn inventory is excess Rn activity times water depth or depth of mixed layer/groundwater plume thickness
    ex_Rn_wat_inv__Bqm2 = ex_Rn_wat__Bqm3 * depth__m,

    # change in radon inventory is the difference between radon inventories in two consecutive time steps
    # the change in radon inventory is radon flux over the measurement time period and is converted to per hour flux
    f_Rn_gross__Bqm2hr = (ex_Rn_wat_inv__Bqm2 - lag(ex_Rn_wat_inv__Bqm2)) * 60 / meas_t__min,

    # Rn brought in by tides from offshore
    f_Rn_flood__Bqm2hr = if_else(depth__m - lag(depth__m) > 0,
      ((depth__m - lag(depth__m)) * Rn_offshore__Bqm3) * 60 / meas_t__min,
      0,
      NA_real_
    ),

    # Rn lost from coastal to offshore areas
    f_Rn_ebb__Bqm2hr = if_else(depth__m - lag(depth__m) < 0,
      (-(depth__m - lag(depth__m)) * ex_Rn_wat__Bqm3) * 60 / meas_t__min,
      0,
      NA_real_
    ),

    # all derived Rn fluxes summed together to get net Rn change per hour
    f_Rn_net__Bqm2hr = f_Rn_gross__Bqm2hr +
      f_Rn_atm__Bqm2hr +
      f_Rn_ebb__Bqm2hr -
      f_Rn_flood__Bqm2hr -
      f_dif__Bqm2hr +
      f_mix_exp__Bqm2hr,

    # if f_mix_exp__Bqm2hr has not been provided then losses by mixing are set to equal negative f_Rn_net__Bqm2hr
    # this is a conservative approach providing minimal estimate of mixing loss
    f_Rn_mix__Bqm2hr = if_else(f_Rn_net__Bqm2hr < 0,
      -f_Rn_net__Bqm2hr,
      0,
      NA_real_
    ),

    # total radon flux is the one corrected for mixing losses
    f_Rn_total__Bqm2hr = f_Rn_net__Bqm2hr + f_Rn_mix__Bqm2hr,

    # groundwater discharge f_gw__m3m2d equals the total Rn flux f_Rn_total__Bqm2hr 
    # divided by groundwater Rn activity
    f_gw__m3m2d = (f_Rn_total__Bqm2hr / Rn_gw__Bqm3) * 24
    
  ) %>% 

  # drop columns with no values (keep those with at least one value)
  select(where(~(!(.x %>% is.na())) %>% any()))

# results saved in a csv file
write_csv(dat_tbl, here(study_folder, "output", "sgd_coastal_ts_rn_budget.csv"), na = "")
  
# END OF RADON BUDGET CALCULATION

# *************************
#  end ----
# *************************
