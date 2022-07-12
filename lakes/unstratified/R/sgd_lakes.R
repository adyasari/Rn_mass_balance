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
study_folder <- "lakes/unstratified"

# input file name
csv_file_in <- "sgd_lake_data_Natasha.csv"

# decay constant of Rn in hours
lambda <- log(2)/(3.84 * 24)

# *************************
#  load data ----
# *************************

# load the time series
in_tbl <- read_csv(here(study_folder, "input", csv_file_in)) %>%
  # mutate(across(.cols = 1, .fns = ~ force_tz(., tzone = ""))) %>%
  mutate(across(.cols = 1, .fns = ~ parse_date_time(., c("ymdHMS", "mdyHM", "mdyHMS")))) %>%
  mutate(across(.cols = -1, .fns = as.numeric))

# *************************
#  calculations ----
# *************************

# create data table
dat_tbl <- in_tbl %>%
  mutate(
    
    # when a single value is missing interpolate linearly
    Rn_wat__Bqm3 = if_else(is.na(Rn_wat__Bqm3), (lag(Rn_wat__Bqm3) + lead(Rn_wat__Bqm3))/2, Rn_wat__Bqm3),
    temp_wat__C = if_else(is.na(temp_wat__C), (lag(temp_wat__C) + lead(temp_wat__C))/2, temp_wat__C),
    sal_wat = if_else(is.na(sal_wat), (lag(sal_wat) + lead(sal_wat))/2, sal_wat),
    
    # calculates coastal radon measurement time interval in minutes based on provided measurement times, 
    # all other time series parameters provided by the user should be averaged to this interval
    meas_t__min = (time %>% as.numeric() - lag(time %>% as.numeric())) / 60,
  
    # change in water depth between two time stamps
    diff_owl__m = depth__m - lag(depth__m),

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
    Rn_wat__Bqm3 = if (!(Rn_exch__Bqm3 %>% is.null()) & (!(Rn_exch__Bqm3 %>% is.na())) %>% any()) {
      Rn_exch__Bqm3 * kw_air
    } else {
      Rn_wat__Bqm3
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

    # The amount of radon diffusing from the bottom sediments can be estimated from an experimentally 
    # defined relationship between 226Ra content of sediments and the corresponding measured 222Rn flux 
    # by diffusion (Burnett et al., 2003).  That empirical relationship is based on experimental data 
    # from several different environments (both marine and fresh):Flux (dpm m-2 day-1) =  495 x 226Ra conc.(dpm g-1) + 18.2 
    # this can be written as f_dif__Bqm2hr = (495 x Ra226_sed__Bqg * 60 + 18.2) / 24
    f_dif__Bqm2hr = if (!(Ra226_sed__Bqg %>% is.null()) & (!(Ra226_sed__Bqg %>% is.na())) %>% any()) {
      if_else(is.na(Ra226_sed__Bqg), 0, (495 * Ra226_sed__Bqg * 60 + 18.2) / 24)
    } else {
      0
    },
    
    # inventory contributed by Rn diffusion from sediments
    inv_dif__Bqm2 = f_dif__Bqm2hr * meas_t__min / 60,
    
    # Rn inventory contributed by decay over measurement time
    inv_dec__Bqm2 = Rn_wat__Bqm3 * exp(-lambda * meas_t__min / 60) * depth__m, 
    
    # Rn inventory contributed by atmospheric evasion over measurement time
    inv_eva__Bqm2 = f_Rn_atm__Bqm2hr * (1 - exp(-lambda * meas_t__min / 60)) / lambda,
    
    # # excess Rn in water is calculated by sutracting dissolved 226Ra in water
    # ex_Rn_wat__Bqm3 = Rn_wat__Bqm3 - Ra226_wat__Bqm3,

  ) %>% 
  
  # drop columns with no values (keep those with at least one value)
  select(where(~(!(.x %>% is.na())) %>% any()))

# objective function to find steady state Rn flux
Rn_ss_calc <- function(ini_Rn_flux__Bqm2hr){
  #input: ini_Rn_flux__Bqm2hr - initial value for Rn flux via SGD
  
  # the calculation is based on transformations of the data
  obj_fun <- dat_tbl %>%
    mutate(
      
      # Rn inventory contributed by SGD over measurement time
      inv_sgd__Bqm2 = ini_Rn_flux__Bqm2hr * (1 - exp(-lambda * meas_t__min / 60)) / lambda,
      
      # Rn inventory contributed by all considered components
      inv_net__Bqm2 = inv_dec__Bqm2 + inv_sgd__Bqm2 - inv_eva__Bqm2 + inv_dif__Bqm2,
      
      # modeled Rn concentration
      mod_Rn_wat__Bqm3 = inv_net__Bqm2 * depth__m,
      
      # difference between actual and modeled Rn concentration
      diff_Rn_wat__Bqm3 = Rn_wat__Bqm3 - mod_Rn_wat__Bqm3,
      
      # root mean squared difference
      rmse_Rn_wat__Bqm3 = (diff_Rn_wat__Bqm3^2) %>% mean(na.rm = TRUE) %>% sqrt()
      
    ) %>% 
    
    # rmse_Rn_wat__Bqm3 is a column with identical values, get one of those
    slice(1) %>% 
    pull(rmse_Rn_wat__Bqm3)
  
}

# minimize objective function with respect to steady state Rn flux
# initial value is based  on:
# initial guess of steady state Rn concentration
# Rn_ss__Bqm3 = Rn_wat__Bqm3,
# Rn inventory in steady state
# inv_Rn_ss__Bqm2 = Rn_ss__Bqm3 * depth__m,
# initial guess of steady state Rn flux
# Rn_flux__Bqm2hr = lambda * inv_Rn_ss__Bqm2,
optim_start <- dat_tbl %>% 
  summarize(mean(lambda * Rn_wat__Bqm3 * depth__m, na.rm = TRUE)) %>% 
  pull()

# compare two optimization techniques
opt1 <- optim(optim_start, Rn_ss_calc, method = "BFGS")
opt2 <- optimize(Rn_ss_calc, c(0, 500))
opt1$par - opt2$minimum

# continue data table
dat_tbl_cont <- dat_tbl %>%
  mutate(

    # steady state Rn flux (from optimization above)
    Rn_flux__Bqm2hr = opt1$par,
    
    # Rn inventory contributed by SGD over measurement time
    inv_sgd__Bqm2 = Rn_flux__Bqm2hr * (1 - exp(-lambda * meas_t__min / 60)) / lambda,
    
    # Rn inventory contributed by all considered components
    inv_net__Bqm2 = inv_dec__Bqm2 + inv_sgd__Bqm2 - inv_eva__Bqm2 + inv_dif__Bqm2,
    
    # modeled Rn concentration
    mod_Rn_wat__Bqm3 = inv_net__Bqm2 * depth__m,
    
    # difference between actual and modeled Rn concentration
    diff_Rn_wat__Bqm3 = Rn_wat__Bqm3 - mod_Rn_wat__Bqm3,
    
    # ground water advection into lake in m/hr which is the same as m3/m2/hr
    drip_rate__mhr = Rn_flux__Bqm2hr / Rn_gw__Bqm3,
    
  ) %>% 
  
  # drop columns with no values (keep those with at least one value)
  select(where(~(!(.x %>% is.na())) %>% any()))

# results saved in a csv file
write_csv(dat_tbl_cont, here(study_folder, "output", "rn_budget.csv"), na = "")
  
# END OF RADON BUDGET CALCULATION

# plot
dat_tbl_cont %>% 
  ggplot(mapping = aes(x = time)) + 
  geom_line(mapping = aes(y = Rn_wat__Bqm3), alpha = 0.6) +
  geom_point(mapping = aes(y = mod_Rn_wat__Bqm3), color = "red") +
  labs(y = "Rn in water: actual (black) vs model (red)") +
  theme_light()
ggsave(here(study_folder, "output", "rn_budget.png"))
  
# *************************
#  end ----
# *************************
