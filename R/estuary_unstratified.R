# *************************
# Rn mass balance calculation for spatial survey estuarine radon budget
# Based off equations in Hussain et al., 1999 (10.1016/S0304-4203(99)00015-8) &
# Schwartz, 2003 (doi: 10.1016/S0272-7714(02)00118-X)
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
csv_file_in <- "estuary_unstratified_data.csv"

# *************************
#  load data ----
# *************************

# load the estuarine survey data
in_tbl <- read_csv(here(study_folder, "input", csv_file_in)) %>%
  mutate(across(.cols = 1, .fns = ~ lubridate::parse_date_time(., c("ymdHMS", "ymdHM", "mdyHM", "mdyHMS")))) %>%
  mutate(across(.cols = -1, .fns = as.numeric))

# *************************
#  calculations ----
# *************************

# creates data table
dat_tbl <- in_tbl %>%
  mutate(

    # decay constant of Rn in hours
    lambda__hr = log(2) / (3.84 * 24),

    # water temperature converted from degrees Celsius to Kelvin
    temp_wat__K = temp_wat__C + 273.15,

    # volume of the box
    v_box__m3 = a_box__m2 * depth__m, 
    
    # average of upstream and downstream values
    Rn_wat__Bqm3 = (Rn_wat_ups__Bqm3 + Rn_wat_dws__Bqm3) / 2, 
    sal_wat = (sal_wat_ups + sal_wat_dws) / 2, 

    # water/air partitioning coefficient kw_air based on water temperature and salinity;
    # calculations and coefficients from Schubert et al. 2012
    kw_air =
      exp(-76.14 + 120.36 * (100 / temp_wat__K) + 31.26 * log(temp_wat__K / 100) + sal_wat *
        (-0.2631 + 0.1673 * (temp_wat__K / 100) + (-0.0270 * (temp_wat__K / 100)^2))) *
        temp_wat__K / 273.15,
    
    # Rn losses by evasion into the atmosphere are calculated either according to Borges et al., 2004 or according to MacIntyre et al. (1995)
    # if wat_current__cms is provided in the spreadsheet, then currents and winds speed are used to estimate turbulence and hence k (Borges et al. 2004), otherwise only wind speed is used (MacIntyre et al. 1995).
    # for wind speeds above 3.6 m/s Sc^-1/2 and for wind speeds below 3.6 m/s Sc^-2/3 is applied (Turner et al., 1996);
    # for wind speeds below 1.5 m/s k is assumed to be constant and equivalent to wind speeds of 1.5 m/s (Ocampo-torres et al., 1994)
    # note that kinematic viscosity is considered constant, one can calculate more accurate values that account for salinity and temperature
    f_Rn_atm__Bqm2hr = if (("wat_current__cms" %in% names(.)) && (!(wat_current__cms %>% is.na())) %>% any()) {
      case_when(
        wind__ms > 3.6 ~
          ((1 + 1.719 * wat_current__cms^(1 / 2) * depth__m^(-1 / 2) + 2.58 * wind__ms) *
            ((0.0086 / 10^(-(980 / temp_wat__K + 1.59))) / 600)^(-1 / 2)) / 100 / 60 *
            (Rn_wat__Bqm3 - kw_air * Rn_air__Bqm3) * 60,
        wind__ms > 1.5 ~
          ((1 + 1.719 * wat_current__cms^(1 / 2) * depth__m^(-1 / 2) + 2.58 * wind__ms) *
            ((0.0086 / 10^(-(980 / temp_wat__K + 1.59))) / 600)^(-2 / 3)) / 100 / 60 *
            (Rn_wat__Bqm3 - kw_air * Rn_air__Bqm3) * 60,
        TRUE ~
          ((1 + 1.719 * wat_current__cms^(1 / 2) * depth__m^(-1 / 2) + 2.58 * 1.5) *
            ((0.0086 / 10^(-(980 / temp_wat__K + 1.59))) / 600)^(-2 / 3)) / 100 / 60 *
            (Rn_wat__Bqm3 - kw_air * Rn_air__Bqm3) * 60
      )
    } else {
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
      )
    },
    
    # volume of the box
    v_box__m3 = a_box__m2 * depth__m, 

    # groundwater discharge into the estuary
    q_gw__m3m2d = (q_dws__m3d * Rn_wat_dws__Bqm3 + # Advection - out (Bq/d)
             f_Rn_atm__Bqm2hr * 24 * a_box__m2 + # Rn flux via atmospheric evasion (Bq/m2/d and then converted to Bq/d)
               Rn_wat_dws__Bqm3 * lambda__hr / 24 * v_box__m3 - # Rn decay out of box (Bq/d)
             q_ups__m3d * Rn_wat_ups__Bqm3 - # Advection - in (Bq/d)
             Ra226_wat__Bqm3 * lambda__hr / 24 * v_box__m3 # Ra-226 production in the box (Bq/d) PF: Rn production????
           ) / a_box__m2 / Rn_gw__Bqm3,
    
  ) %>%
  
  # add a row with total gw for estuary 
  add_row(q_gw__m3m2d = sum(.$q_gw__m3m2d)) %>% 
  
  # drop columns with no values (keep those with at least one value)
  select(where(~ (!(.x %>% is.na())) %>% any()))

# results saved in a csv file
write_csv(dat_tbl, here(study_folder, "output", "estuary_unstratified_rn_budget.csv"), na = "")

# END OF RADON BUDGET AND GROUNDWATER FLUX CALCULATION

# *************************
#  end ----
# *************************
