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
  mutate(across(.cols = 1, .fns = ~as_datetime(., tz = Sys.timezone(location = TRUE)))) %>% 
  mutate(across(.cols = -1, .fns = as.numeric)) %>% 
  ts_long() %>% 
  ts_xts()

# *************************
#  calculations ----
# *************************

# sample expression 
temp1 <- ts_lag(timser_in$wind__ms, "1 day") * single_in$Rn_gw_err__Bqm3
temp2 <- ts_lag(timser_in$wind__ms, "2 min") * single_in$Rn_gw_err__Bqm3
temp2 <- ts_lag(timser_in$wind__ms, "0 min") * single_in$Rn_gw_err__Bqm3
temp3 <- timser_in$wind__ms * single_in$Rn_gw_err__Bqm3


