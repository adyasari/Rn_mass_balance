# *************************
# Rn mass balance calculation for spatial survey estuarine radon budget
# Based off of eq. in Santos et al., 2010; doi:10.1016/j.marchem.2010.03.003
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
csv_file_in <- "sgd_lake_stratified_data.csv"

# *************************
#  load data ----
# *************************

# load the stratified lake data
in_tbl <- read_csv(here(study_folder, "input", csv_file_in)) %>%
  # mutate(across(.cols = 1, .fns = ~ force_tz(., tzone = ""))) %>%
  mutate(across(.cols = 1, .fns = ~ lubridate::parse_date_time(., c("ymdHMS", "ymdHM", "mdyHM", "mdyHMS")))) %>%
  mutate(across(.cols = -c(time, layerID), .fns = as.numeric))

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

    # water/air partitioning coefficient kw_air based on water temperature and salinity (add salinity 0.1 if unknown);
    # calculations and coefficients from Schubert et al. 2012
    kw_air =
      exp(-76.14 + 120.36 * (100 / temp_wat__K) + 31.26 * log(temp_wat__K / 100) + sal_wat *
        (-0.2631 + 0.1673 * (temp_wat__K / 100) + (-0.0270 * (temp_wat__K / 100)^2))) *
        temp_wat__K / 273.15,

    # if Rad-Aqua was used to collect radon data and radon in the exchanger (therefore in air) is provided 
    # it is converted to Rn in water in this step;
    # otherwise, the provided radon in water is used for further calculations
    # the condition checks if a Rn_exch__Bqm3 column exists and if it is non-empty
    Rn_wat__Bqm3 = if (("Rn_wat__Bqm3" %in% names(.)) && (!(Rn_wat__Bqm3 %>% is.na())) %>% any()) {
      Rn_wat__Bqm3
    } else {
      Rn_exch__Bqm3 * kw_air
    },
    
    # Rn losses by evasion into the atmosphere are calculated either according to MacIntyre et al. (1995)
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

    # groundwater discharge into the surface water layer above the pycnocline called epilmnion (check units)
    q_epi__m3d = if_else(str_detect(layerID, "epi"), (Rn_wat__Bqm3 * lambda__hr * d_box__m + f_Rn_atm__Bqm2hr) / Rn_gw__Bqm3 * a_box__m2, NA_real_),
    # groundwater discharge into the middle water layer called metalimnion isolated from the atmosphere by epilimnion and bottom sediments by hypolimnion (check units)
    q_meta__m3d = if_else(str_detect(layerID, "meta"), (Rn_wat__Bqm3 * lambda__hr * d_box__m) / Rn_gw__Bqm3 * a_box__m2, NA_real_),
    # groundwater discharge into the bottom water layer called hypolimnion below the pycnocline (check units)
    q_hypo__m3d = if_else(str_detect(layerID, "hypo"), (Rn_wat__Bqm3 * lambda__hr * d_box__m + f_dif__Bqm2hr) / Rn_gw__Bqm3 * a_box__m2, NA_real_),
    
    # total groundwater discharge into the lake accounting for inputs into all its layers (check units)
    q_gw__m3m2d = mean(q_epi__m3d, na.rm = TRUE) + mean(q_meta__m3d, na.rm = TRUE) + mean(q_hypo__m3d, na.rm = TRUE),
    
  ) %>%
  
  # drop columns with no values (keep those with at least one value)
  select(where(~ (!(.x %>% is.na())) %>% any()))

# results saved in a csv file
write_csv(dat_tbl, here(study_folder, "output", "sgd_lake_stratified_rn_budget.csv"), na = "")

# END OF RADON BUDGET AND GROUNDWATER FLUX CALCULATION

# *************************
#  end ----
# *************************
