#
# Macroeconomic Forecasting for Hawaii
#
# Utility Functions for Time Series in the Model
#
# Author: Peter Fuleky
# Date: October 1, 2019

# generate help files for these functions
# d <- document::document(file_name = here::here("R", "util_funs.R"), output_directory = here::here("man"))

#' Download a single series from udaman using series name
#'
#' @param ser_id udaman series name
#' @param expand "true" or "raw" ("true" downloads formatted data, "raw" downloads raw units)
#'
#' @return time and data for a single series combined in a tibble
#' @export
#'
#' @examples
#' get_series_1("VISNS@HI.M")
get_series_1 <- function(ser_id, expand = "true") {
  # API call
  url <- str_c("https://api.uhero.hawaii.edu/v1.u/series?name=", ser_id, "&expand=", expand, "&u=uhero&nocache")
  req <- httr::GET(url, httr::add_headers(Authorization = str_c("Bearer ", Sys.getenv("udaman_token"))))
  json <- httr::content(req, as = "text")
  uhero_data <- jsonlite::fromJSON(json)
  # extract series info
  dates <- uhero_data$data$observations$transformationResults$dates[[1]]
  values <- uhero_data$data$observations$transformationResults$values[[1]]
  series <- bind_cols(time = ymd(dates), values = as.numeric(values))
  colnames(series) <- c("time", uhero_data$data$series$name)
  # colnames(series) <- c("time", uhero_data$data$series$name %>% str_replace_all(c("\\.[A-Z]" = "", "@" = "__", "OCUP%" = "OCUPP")))
  # return series
  return(series)
}


#' Download a set of series from udaman using series names
#'
#' @param ser_id_vec vector of series names
#' @param format "wide" (default) or "long" or "xts"
#' @param expand "true" (default) or "raw" ("true" downloads formatted data, "raw" downloads raw units)
#'
#' @return time and data for all series combined in an object specified by the format option
#' @export
#'
#' @examples
#' get_series(c("VISNS@HI.M", "VAPNS@HI.M"))
#' get_series(c("VISNS@HI.M", "VAPNS@HI.M"), format = "xts")
#' get_series(c("VISNS@HI.M"), format = "xts")
get_series <- function(ser_id_vec, format = "wide", expand = "true") {
  ser_tbl <- ser_id_vec %>%
    map(get_series_1, expand) %>%
    reduce(full_join, by = "time") %>%
    arrange(time)
  if (format == "wide") ser_out <- ser_tbl
  if (format == "long") ser_out <- ser_tbl %>% ts_long()
  if (format == "xts") ser_out <- ser_tbl %>% ts_long() %>% ts_xts()
  return(ser_out) 
}


#' Download series listed in an export table from udaman
#'
#' @param exp_id export id
#' @param format "wide" (default) or "long" or "xts"
#' @param save_loc location to save the csv of the retrieved data, set to NULL to avoid saving
#'
#' @return time and data for all series combined in a tibble
#' @export
#'
#' @examples
#' get_series_exp(74)
#' get_series_exp(74, format = "xts")
#' get_series_exp(74, save_loc = NULL)
get_series_exp <- function(exp_id, format = "wide", save_loc = "data/raw") {
  url <- "https://udaman.uhero.hawaii.edu/"
  dn_url <- str_c("https://udaman.uhero.hawaii.edu/exports/", exp_id, ".csv")
  session <- rvest::session(url)
  form <- rvest::html_form(session)[[1]]
  fl_fm <- rvest::html_form_set(form,
    `user[email]` = Sys.getenv("udaman_user"),
    `user[password]` = Sys.getenv("udaman_pwd")
  )
  main_page <- rvest::session_submit(session, fl_fm)
  download <- rvest::session_jump_to(main_page, dn_url)
  content <- readBin(download$response$content, what = character())[1]
  data_tbl <- read_csv(file = content) %>% 
    rename(time = date)
    # rename_with(~str_replace_all(., c("date" = "time", "\\.[A-Z]" = "", "@" = "__", "OCUP%" = "OCUPP")))
  if (!is.null(save_loc)) write_csv(data_tbl, here(save_loc, basename(dn_url)))
  if (format == "wide") data_out <- data_tbl
  if (format == "long") data_out <- data_tbl %>% ts_long()
  if (format == "xts") data_out <- data_tbl %>% ts_long() %>% ts_xts()
  return(data_out) 
}


#' Interpolate a single series from quarterly to monthly freq
#'
#' @param var_q vector containing a single variable at quarterly freq
#' @param ts_start starting period as c(year, quarter) e.g. c(2001, 1)
#' @param conv_type match the quarterly value via "first", "last", "sum", "average"
#'
#' @return vector containing a single variable at monthly freq
#' @export
#'
#' @examples
#' QtoM_1(test1, c(2010, 1), "average")
QtoM_1 <- function(var_q, ts_start, conv_type) {
  var_q_ts <- ts(var_q, frequency = 4, start = ts_start)
  tempdisagg::td(
    formula = var_q_ts ~ 1,
    conversion = conv_type,
    to = "monthly",
    method = "denton-cholette"
  ) %>%
    predict()
}


#' Interpolate a tibble of series from quaterly to monthly freq
#'
#' @param data_q tibble containing variables at quarterly freq
#           the first column of data_q named "time" contains dates
#' @param conv_type match the quarterly value via "first", "last", "sum", "average"
#'
#' @return tibble containing variables at monthly freq
#' @export
#'
#' @examples
#' `ncen@us.sola` <- ts(NA, start = 2016, end = 2021, freq = 1) %>% ts_xts()
#' `ncen@us.sola`["2016/2021"] <- c(323127513, 325511184, 327891911, 330268840, 332639102, 334998398)
#' test1 <- AtoQ(`ncen@us.sola`)
#' QtoM(ts_tbl(test1), "average")
#' ts_frequency(QtoM(ts_tbl(test1), "average") %>% ts_xts())
QtoM <- function(data_q, conv_type) {
  data_q_names <- colnames(data_q)
  data_q_dates <- ymd(data_q$time)
  data_q_first <- first(data_q_dates)
  data_q_last <- last(data_q_dates)
  data_q_start <- c(year(data_q_first), quarter(data_q_first))
  data_m <- data_q %>%
    select(-time) %>%
    map(QtoM_1, data_q_start, conv_type) %>%
    reduce(ts.union) %>%
    as_tibble()
  data_m_dates <- seq(data_q_first, data_q_last + months(2), by = "months") %>%
    enframe(name = NULL)
  data_m <- bind_cols(data_m_dates, data_m) %>%
    rename_all(~data_q_names)
  return(data_m)
}


#' Linear interpolation based on aremos command reference page 292
#'
#' @param ser_in the xts series to be interpolated (freq = a)
#' @param aggr interpolation method: aggregate via mean (default) or sum
#'
#' @return interpolated xts series (freq = q)
#' @export
#'
#' @examples
#' `ncen@us.sola` <- ts(NA, start = 2016, end = 2021, freq = 1) %>% ts_xts()
#' `ncen@us.sola`["2016/2021"] <- c(323127513, 325511184, 327891911, 330268840, 332639102, 334998398)
#' test1 <- AtoQ(`ncen@us.sola`)
AtoQ <- function(ser_in, aggr = "mean") {
  ser_out_name <- names(ser_in)
  ser_out_dates <- tibble(time = seq.Date(
    from = ts_summary(ser_in)$start,
    to = ts_summary(ser_in)$end %>% ceiling_date(unit = "year") %>% rollback(),
    by = "quarter"
  ))
  ser_out <- left_join(ser_out_dates, ser_in %>% ts_tbl()) %>%
    fill(value) %>%
    ts_xts() %>%
    magrittr::set_names(ser_out_name)
  dat_start <- ts_summary(ser_out)$start
  dat_end <- ts_summary(ser_out)$end
  increment <- (ts_lag(ser_out, -4) - ser_out) / 4
  increment <- increment %>% ts_bind(ts_lag(increment, 4)[str_c(dat_end - months(9), "/", dat_end)])
  ser_out[p(dat_start + months(3), dat_end)] <- (as.numeric(ser_out[dat_start]) + ts_lag(increment, 1)[p(dat_start + months(3), dat_end)] %>% cumsum()) %>% as.numeric()
  ser_out <- ser_out - 1.5 * increment
  if (aggr != "mean") ser_out <- ser_out / 4
  colnames(ser_out) <- colnames(ser_in) %>%
    gsub(".SOLA", ".SOLQ", .) %>%
    gsub(".A", ".Q", .)
  return(ser_out)
}

#' Conversion from quarterly to annual frequency
#'
#' @param ser_in the xts series to be converted (freq = q)
#' @param	aggr aggregate via mean (default) or sum
#'
#' @return converted xts series (freq = a)
#' @export
#'
#' @examples
#' `ncen@us.sola` <- ts(NA, start = 2016, end = 2021, freq = 1) %>% ts_xts()
#' `ncen@us.sola`["2016/2021"] <- c(323127513, 325511184, 327891911, 330268840, 332639102, 334998398)
#' test1 <- AtoQ(`ncen@us.sola`)
#' test2 <- QtoA(test1) # for stock type variables mean, for flow type variables sum
#' print(test1)
#' print(cbind(`ncen@us.sola`, test2))
QtoA <- function(ser_in, aggr = "mean") {
  ser_out <- ts_frequency(ser_in, to = "year", aggregate = aggr)
  colnames(ser_out) <- colnames(ser_in) %>%
    gsub(".SOLQ", ".SOLA", .) %>%
    gsub(".Q", ".A", .)
  return(ser_out)
}

#' Year to date sum or average
#'
#' @param long_tbl_in a long tibble of time series (produced by ts_long() for example)
#' @param avg if true, the year to date average, if false, the year to date sum
#'
#' @return a long tibble of time series containing year to date sum or average
#' @export
#'
#' @examples
#' get_series(c("VISNS@HI.M", "VAPNS@HI.M")) %>% ts_long() %>% ytd_cum()
ytd_cum <- function(long_tbl_in, avg = TRUE) {
  long_tbl_out <- long_tbl_in %>%
    mutate(yr = floor_date(time, "year")) %>% 
    group_by(id, yr) %>% 
    mutate(value = if(avg) cummean(value) else cumsum(value)) %>% 
    ungroup() %>% 
    select(-yr)
  return(long_tbl_out)
}

#' Year to date growth rate
#'
#' @param long_tbl_in a long tibble of time series (produced by ts_long() for example)
#' @param avg if true, the year to date average, if false, the year to date sum for calculation
#'
#' @return a long tibble of time series containing year to date growth rate
#' @export
#'
#' @examples
#' get_series(c("VISNS@HI.M", "VAPNS@HI.M")) %>% ts_long() %>% ytd_gr() %>% tail()
ytd_gr <- function(long_tbl_in, avg = TRUE) {
  long_tbl_out <- long_tbl_in %>%
    ytd_cum(avg) %>% 
    ts_pcy()
  return(long_tbl_out)
}

#' Month to date sum or average
#'
#' @param long_tbl_in a long tibble of time series (produced by ts_long() for example)
#' @param avg if true, the year to date average, if false, the year to date sum
#'
#' @return a long tibble of time series containing year to date sum or average
#' @export
#'
#' @examples
#' get_series(c("VISPNS@HI.D", "VAPNS@HI.D")) %>% ts_long() %>% mtd_cum()
#' test <- get_series("VAPNS@HI.D") %>% ts_long() %>% mtd_cum()
#' test %ts/% ts_lag(test, "3 years") %>% tail()
mtd_cum <- function(long_tbl_in, avg = TRUE) {
  long_tbl_out <- long_tbl_in %>%
    mutate(yrmo = floor_date(time, "month")) %>% 
    group_by(id, yrmo) %>% 
    mutate(value = if(avg) cummean(value) else cumsum(value)) %>% 
    ungroup() %>% 
    select(-yrmo)
  return(long_tbl_out)
}

#' Month to date growth rate
#'
#' @param long_tbl_in a long tibble of time series (produced by ts_long() for example)
#' @param avg if true, the year to date average, if false, the year to date sum for calculation
#'
#' @return a long tibble of time series containing year to date growth rate
#' @export
#'
#' @examples
#' get_series(c("VISPNS@HI.D", "VAPNS@HI.D")) %>% ts_long() %>% mtd_gr() %>% tail()
mtd_gr <- function(long_tbl_in, avg = TRUE) {
  long_tbl_out <- long_tbl_in %>%
    mtd_cum(avg) %>% 
    ts_pcy()
  return(long_tbl_out)
}

#' Period to date sum or average
#'
#' @param long_tbl_in a long tibble of time series (produced by ts_long() for example)
#' @param per unit of time supplied to floor_date() (for ytd per = "year" (default), for mtd per = "month")
#' @param avg if true (default), the year to date average, if false, the year to date sum
#'
#' @return a long tibble of time series containing year to date sum or average
#' @export
#'
#' @examples
#' get_series(c("VISNS@HI.M", "VAPNS@HI.M")) %>% ts_long() %>% ptd_cum()
ptd_cum <- function(long_tbl_in, per = "year", avg = TRUE) {
  long_tbl_out <- long_tbl_in %>%
    mutate(time_per = floor_date(time, per)) %>% 
    group_by(id, time_per) %>% 
    mutate(value = if(avg) cummean(value) else cumsum(value)) %>% 
    ungroup() %>% 
    select(-time_per)
  return(long_tbl_out)
}

#' Period to date growth rate
#'
#' @param long_tbl_in a long tibble of time series (produced by ts_long() for example)
#' @param per unit of time supplied to floor_date() (for ytd per = "year" (default), for mtd per = "month")
#' @param lag_length period over which growth is calculated (e.g. "1 year" (default), "3 years", etc. See ?ts_lag() for options)
#' @param avg if true, the year to date average, if false, the year to date sum for calculation
#'
#' @return a long tibble of time series containing year to date growth rate
#' @export
#'
#' @examples
#' get_series(c("VISNS@HI.M", "VAPNS@HI.M")) %>% ts_long() %>% ptd_gr() %>% tail()
#' get_series("VAPNS@HI.D") %>% ts_long() %>% ptd_gr(per = "month", lag_length = "3 years") %>% tail()
ptd_gr <- function(long_tbl_in, per = "year", lag_length = "1 year", avg = TRUE) {
  long_tbl_out <- long_tbl_in %>%
    ptd_cum(per, avg) %>% 
    {(. %ts/% ts_lag(., lag_length) %ts-% 1) %ts*% 100}
  return(long_tbl_out)
}

#' Splitting of xts matrix to individual xts vectors (don't use, pollutes global environment)
#'
#' @param xts_in the xts matrix to be split into individual xts vectors
#'
#' @return nothing (silently store split series in global environment)
#' @export
#'
#' @examples
#' get_series_exp(74, save_loc = NULL) %>%
#'   ts_long() %>%
#'   ts_xts() %>%
#'   explode_xts()
#' rm(list = ls(pattern = glob2rx("*@HI.Q")))
explode_xts <- function(xts_in) {
  temp <- xts_in %>%
    ts_tbl() %>%
    ts_wide()
  for (i in names(temp)[-1]) {
    assign(i, xts::xts(temp[i], temp %>% pull(time)), envir = .GlobalEnv)
  }
  return("OK. Done.")
}


#' Convert annualized growth to quarterly growth
#'
#' @param ser_in the series containing annualized growth (in percent)
#'
#' @return series containing quarterly growth (in percent)
#' @export
#'
#' @examples
#' `ncen@us.sola` <- ts(NA, start = 2016, end = 2021, freq = 1) %>% ts_xts()
#' `ncen@us.sola`["2016/2021"] <- c(323127513, 325511184, 327891911, 330268840, 332639102, 334998398)
#' test1 <- AtoQ(`ncen@us.sola`)
#' ts_c(test1 %>% ts_pca() %>% pca_to_pc(), test1 %>% ts_pc())
pca_to_pc <- function(ser_in) {
  ((1 + ser_in / 100)^0.25 - 1) * 100
}


#' Find the date of the first observation (NAs are dropped)
#'
#' @param ser_in an xts series
#'
#' @return date associated with first observation
#' @export
#'
#' @examples
#' `ncen@us.sola` <- ts(NA, start = 2016, end = 2021, freq = 1) %>% ts_xts()
#' `ncen@us.sola`["2017/2021"] <- c(325511184, 327891911, 330268840, 332639102, 334998398)
#' find_start(`ncen@us.sola`)
find_start <- function(ser_in) {
  # ser_in %>% na.omit() %>% start()
  ser_in %>%
    ts_summary() %>%
    pull(start)
}


#' Find the date of the last observation (NAs are dropped)
#'
#' @param ser_in an xts series
#'
#' @return date associated with last observation
#' @export
#'
#' @examples
#' `ncen@us.sola` <- ts(NA, start = 2016, end = 2060, freq = 1) %>% ts_xts()
#' `ncen@us.sola`["2016/2018"] <- c(323127513, 325511184, 327891911)
#' find_end(`ncen@us.sola`)
find_end <- function(ser_in) {
  # ser_in %>% na.omit() %>% end()
  ser_in %>%
    ts_summary() %>%
    pull(end)
}


#' Construct a series name from variable components and retrieve the series
#'
#' @param ser_in a variable name (string with substituted expressions)
#' @param env environment where the expression should be evaluated
#'
#' @return variable
#' @export
#'
#' @examples
#' ser_i <- "_NF"
#' cnty_i <- "HI"
#' get_series_exp(74, save_loc = NULL) %>%
#'   ts_long() %>%
#'   ts_xts() %$% get_var("E{ser_i}@{cnty_i}.Q")
get_var <- function(ser_in, env = parent.frame()) {
  return(ser_in %>% str_glue() %>% get(envir = env, inherits = TRUE))
}


#' Concatenate dates to obtain period
#'
#' @param dat1 date of period start (string: yyyy-mm-dd)
#' @param dat2 date of period end (string: yyyy-mm-dd)
#'
#' @return string containing date range
#' @export
#'
#' @examples
#' p("2010-01-01", "2020-01-01")
#' p(2010, 2020) # for annual period only
p <- function(dat1 = "", dat2 = "") {
  str_c(dat1, dat2, sep = "/")
}


#' Concatenate dates formatted as yyyyMm or yyyy.m to obtain period
#'
#' @param dat1 date of period start (string: yyyyMm or yyyy.m)
#' @param dat2 date of period end (string: yyyyMm or yyyy.m)
#'
#' @return string containing date range
#' @export
#'
#' @examples
#' pm("2010M1", "2020M4")
#' pm(2010.1, 2020.4)
#' pm(2010.1, )
#' pm(, 2010.1)
pm <- function(dat1 = "", dat2 = "") {
  str_c(if(dat1 != "") dat1 %>% ym(), "/", if(dat2 != "") dat2 %>% ym())
}


#' Concatenate dates formatted as yyyyQq or yyyy.q to obtain period
#'
#' @param dat1 date of period start (string: yyyyQq or yyyy.q)
#' @param dat2 date of period end (string: yyyyQq or yyyy.q)
#'
#' @return string containing date range
#' @export
#'
#' @examples
#' pq("2010Q1", "2020Q4")
#' pq(2010.1, 2020.4)
#' pq(2010.1, )
#' pq(, 2010.1)
pq <- function(dat1 = "", dat2 = "") {
  str_c(if(dat1 != "") dat1 %>% yq(), "/", if(dat2 != "") dat2 %>% yq())
}


#' Concatenate dates formatted as yyyy to obtain period
#'
#' @param dat1 year of period start (string or numeric: yyyy)
#' @param dat2 year of period end (string or numeric: yyyy)
#'
#' @return string containing date range
#' @export
#'
#' @examples
#' py("2010", "2020")
#' py(2010, 2020)
#' py(2010, )
#' py(, 2010)
py <- function(dat1 = "", dat2 = "") {
  str_c(if(dat1 != "") str_c(dat1, "-01-01"), "/", if(dat2 != "") str_c(dat2, "-01-01"))
}


#' Convert dates from yyyy-mm-dd to yyyyQqq format
#'
#' @param x dates (string: yyyy-mm-dd)
#'
#' @return formatted dates (string: yyyyQqq)
#' @export
#'
#' @examples
#' ymd_to_yQq(c("2010-01-01", "2020-10-01"))
#' ymd_to_yQq(c("2010-01-01", "2020-10-01")) %>% yq()
ymd_to_yQq <- function(x) {
  x %>% quarter(type = "year.quarter") %>% str_replace("\\.", "Q")
}


#' Create xts and fill with values
#'
#' @param start date of series start (string: "yyyy-mm-dd")
#' @param end date of series end (string: "yyyy-mm-dd")
#' @param per periodicity of series (string: "quarter", "year")
#' @param val values to fill in (scalar or vector)
#'
#' @return an xts series
#' @export
#'
#' @examples
#' make_xts()
#' make_xts(start = ymd("2010-01-01"), per = "quarter", val = 0)
make_xts <- function(start = bnk_start, end = bnk_end, per = "year", val = NA) {
  tibble(
    time = seq.Date(start, end, by = per),
    value = val
  ) %>%
    ts_xts()
}


#' Convert period in quarters to period months
#'
#' @param nr_quarters number of quarters in period (integer)
#'
#' @return number of months in period
#' @export
#'
#' @examples
#' qtrs(3)
#' ymd("2020-01-01") + qtrs(3)
qtrs <- function(nr_quarters) {
  nr_quarters * months(3)
}


#' Calculate multi-period average growth
#'
#' @param ser_in name of xts series for which growth is calculated
#' @param lag_in length of period over which growth is calculated
#'
#' @return series containing the average growth of ser_in (in percent)
#' @export
#'
#' @examples
#' `ncen@us.sola` <- ts(NA, start = 2016, end = 2021, freq = 1) %>% ts_xts()
#' `ncen@us.sola`["2016/2021"] <- c(323127513, 325511184, 327891911, 330268840, 332639102, 334998398)
#' test1 <- AtoQ(`ncen@us.sola`)
#' ts_c(pchmy(`ncen@us.sola`, lag_in = 3), ts_pc(`ncen@us.sola`))
#' ts_c(pchmy(test1, lag_in = 4), ts_pcy(test1), ts_pca(test1), ts_pc(test1))
pchmy <- function(ser_in, lag_in = 1) {
  ser_in <- ts_xts(ser_in)
  ser_out <- (((ser_in / ts_lag(ser_in, lag_in))^(1 / lag_in)) - 1) * 100
  return(ser_out)
}


#' Interactive plot of a single variable with level and growth rate
#'
#' @param ser time series to plot
#' @param rng_start start of zoom range ("YYYY-MM-DD")
#' @param height height of a single panel (px)
#' @param width width of a single panel (px)
#'
#' @return a dygraph plot
#' @export
#'
#' @examples
#' `ncen@us.sola` <- ts(NA, start = 2016, end = 2021, freq = 1) %>% ts_xts()
#' `ncen@us.sola`["2016/2021"] <- c(323127513, 325511184, 327891911, 330268840, 332639102, 334998398)
#' test1 <- AtoQ(`ncen@us.sola`)
#' plot_1(`ncen@us.sola`, rng_start = "2017-01-01")
#' plot_1(test1, rng_start = "2017-01-01")
plot_1 <- function(ser, rng_start = as.character(Sys.Date() - years(15)), height = 300, width = 900) {
  ser_name <- ser %>%
    ts_xts() %>%
    names()
  ser_plot <- ser %>%
    ts_xts() %>%
    ts_c(ts_pc(.)) %>%
    ts_tslist() %>%
    set_names(c(ser_name, str_glue("{ser_name}%"))) %>%
    ts_dygraphs(main = "", height = height, width = width) %>%
    dygraphs::dyAxis("y", label = "% change") %>%
    dygraphs::dyAxis("y2", label = "level") %>%
    dygraphs::dyLegend(width = width * 0.90) %>%
    dygraphs::dySeries(str_glue("{ser_name}"), axis = "y2") %>%
    dygraphs::dyBarSeries(str_glue("{ser_name}%"), axis = "y") %>%
    # dygraphs::dyOptions(colors = RColorBrewer::brewer.pal(length(ser_name), "Set1")) %>%
    dygraphs::dyRangeSelector(dateWindow = c(rng_start, as.character(Sys.Date())), height = 30, strokeColor = "red")
  ser_plot
}


#' Two-panel plot of levels, index, and growth rates
#'
#' @param sers a vector of series to plot
#' @param rng_start start of the zoom range ("YYYY-MM-DD")
#' @param rng_end end of the zoom range ("YYYY-MM-DD")
#' @param height height of a single panel (px)
#' @param width width of a single panel (px)
#'
#' @return a list with two dygraph plots (level, index, growth)
#' @export
#'
#' @examples
#' `ncen@us.sola` <- ts(NA, start = 2016, end = 2021, freq = 1) %>% ts_xts()
#' `ncen@us.sola`["2016/2021"] <- c(323127513, 325511184, 327891911, 330268840, 332639102, 334998398)
#' test1 <- AtoQ(`ncen@us.sola`)
#' plot_comp_2(ts_c(`ncen@us.sola`, test1), rng_start = "2017-01-01")
#' get_series_exp(74, save_loc = NULL) %>%
#'   ts_long() %>%
#'   ts_xts() %>%
#'   extract(, c("E_NF@HI.Q", "ECT@HI.Q", "EMN@HI.Q")) %>%
#'   plot_comp_2()
plot_comp_2 <- function(sers, rng_start = as.character(Sys.Date() - years(15)), rng_end = as.character(Sys.Date()), height = 300, width = 900) {
  ser_names <- sers %>%
    ts_xts() %>%
    names()
  plot_level <-
    sers %>%
    ts_xts() %>%
    ts_dygraphs(main = "Level", group = "comp", height = height, width = width) %>%
    dygraphs::dyLegend(width = width * 0.90) # %>%
  # dygraphs::dyOptions(colors = RColorBrewer::brewer.pal(length(ser_names), "Set2"))
  plot_growth <-
    sers %>%
    ts_xts() %>%
    ts_pcy() %>%
    ts_dygraphs(main = "Growth", group = "comp", height = height, width = width) %>%
    dygraphs::dyBarChart() %>%
    dygraphs::dyLegend(width = width * 0.90) %>%
    # dygraphs::dyOptions(colors = RColorBrewer::brewer.pal(length(ser_names), "Set2")) %>%
    dygraphs::dyRangeSelector(dateWindow = c(rng_start, rng_end), height = 30, strokeColor = "red")

  # render the dygraphs objects using htmltools
  list(plot_level, plot_growth) %>%
    htmltools::tagList() %>%
    htmltools::browsable()
}


#' Three-panel plot of levels, index, and growth rates
#'
#' @param sers a vector of series to plot
#' @param rng_start start of the zoom range ("YYYY-MM-DD")
#' @param rng_end end of the zoom range ("YYYY-MM-DD")
#' @param height height of a single panel (px)
#' @param width width of a single panel (px)
#'
#' @return a list with three dygraph plots (level, index, growth)
#' @export
#'
#' @examples
#' `ncen@us.sola` <- ts(NA, start = 2016, end = 2021, freq = 1) %>% ts_xts()
#' `ncen@us.sola`["2016/2021"] <- c(323127513, 325511184, 327891911, 330268840, 332639102, 334998398)
#' test1 <- AtoQ(`ncen@us.sola`)
#' plot_comp(ts_c(`ncen@us.sola`, test1), rng_start = "2017-01-01")
#' get_series_exp(74, save_loc = NULL) %>%
#'   ts_long() %>%
#'   ts_xts() %>%
#'   extract(, c("E_NF@HI.Q", "ECT@HI.Q", "EMN@HI.Q")) %>%
#'   plot_comp()
plot_comp_3 <- function(sers, rng_start = as.character(Sys.Date() - years(15)), rng_end = as.character(Sys.Date()), height = 300, width = 900) {
  ser_names <- sers %>%
    ts_xts() %>%
    names()
  plot_level <-
    sers %>%
    ts_xts() %>%
    ts_dygraphs(main = "Level", group = "comp", height = height, width = width) %>%
    dygraphs::dyLegend(width = width * 0.90) # %>%
  # dygraphs::dyOptions(colors = RColorBrewer::brewer.pal(length(ser_names), "Set2"))
  # plot_level[["elementId"]] <- ser_names %>% extract(1) %>% str_extract("^.*@")
  plot_index <-
    sers %>%
    ts_xts() %>%
    ts_dygraphs(main = "Index", group = "comp", height = height, width = width) %>%
    dygraphs::dyRebase(value = 100) %>%
    dygraphs::dyLegend(width = width * 0.90) # %>%
  # dygraphs::dyOptions(colors = RColorBrewer::brewer.pal(length(ser_names), "Set2"))
  plot_growth <-
    sers %>%
    ts_xts() %>%
    ts_pc() %>%
    ts_dygraphs(main = "Growth", group = "comp", height = height, width = width) %>%
    dygraphs::dyBarChart() %>%
    dygraphs::dyLegend(width = width * 0.90) %>%
    # dygraphs::dyOptions(colors = RColorBrewer::brewer.pal(length(ser_names), "Set2")) %>%
    dygraphs::dyRangeSelector(dateWindow = c(rng_start, rng_end), height = 30, strokeColor = "red")
  
  # render the dygraphs objects using htmltools
  list(plot_level, plot_index, plot_growth) %>%
    htmltools::tagList() %>%
    htmltools::browsable()
}


#' Parse lm() output and convert into bimets equation (GETS model development)
#'
#' @param model a model estimated by lm() (lm object)
#'
#' @return a character vector containing the estimated equation [1] and bimets components [2:4]
#' @export
#'
#' @examples
#' this function combines coefficient estimates and variable names into an equation in vector element 1 
#' and into bimets components in vector elements 2-4.
#' https://stats.stackexchange.com/questions/63600/how-to-translate-the-results-from-lm-to-an-equation
model_equation <- function(model, ...) { #   model =  est_lm   {model_equation(est_lm)[2:4]}
  format_args <- list(...)
  
  model_coeff <- model$coefficients
  format_args$x <- abs(model$coefficients)
  model_coeff_sign <- sign(model_coeff)
  model_coeff_prefix <- case_when(
    model_coeff_sign == -1 ~ " - ",
    model_coeff_sign == 1 ~ " + ",
    model_coeff_sign == 0 ~ " + "
  )
  
  # model_eqn <- paste(strsplit(as.character(model$call$formula), "~")[[2]], # 'y'
  # model_eqn <- paste(strsplit(as.character(model$full_formula), "~")[[2]], # 'y'
  model_eqn <- paste(
    strsplit(as.character(model$terms), "~")[[2]], # 'y'
    "=",
    paste(if_else(model_coeff[1] < 0, "- ", ""),
          do.call(format, format_args)[1],
          paste(model_coeff_prefix[-1],
                do.call(format, format_args)[-1],
                " * ",
                names(model_coeff[-1]),
                sep = "", collapse = ""
          ),
          sep = ""
    )
  )
  
  model_eqn_bim <- paste(
    strsplit(as.character(model$terms), "~")[[2]], # 'y'
    "=",
    paste("b0",
          paste(paste(" + b", 1:length(model_coeff_prefix[-1]), sep = ""),
                " * ",
                names(model_coeff[-1]),
                sep = "", collapse = ""
          ),
          sep = ""
    )
  )
  
  model_eqn_beh <- str_extract(model_eqn_bim, "\\w*") %>% # extract the target variable name 
    gsub("DL_([_.[:alnum:]]+)", "TSDELTALOG(\\1)", .) %>% # and parse DL_, L_, D_ to bimets code
    gsub( "L_([_.[:alnum:]]+)", "LOG(\\1)", .) %>% 
    gsub( "D_([_.[:alnum:]]+)", "TSDELTA(\\1)", .) 
  model_eqn_bim <- model_eqn_bim %>% 
    gsub("DL_([_.[:alnum:]]+)", "TSDELTALOG(\\1)", .) %>% # replace DL_, L_, D_ with TSDELTALOG(), LOG(), TSDELTA()
    gsub( "L_([_.[:alnum:]]+)", "LOG(\\1)", .) %>% 
    gsub( "D_([_.[:alnum:]]+)", "TSDELTA(\\1)", .) 
  model_eqn_bim <- gsub("([_[:alpha:]]+)(\\.)([[:digit:]]{1})", "TSLAG(\\1\\, \\3\\)", model_eqn_bim) # replace dot notation for lags with TSLAG()
  model_eqn_coe <- str_extract_all(model_eqn_bim, "(b[[:digit:]]+)", simplify = TRUE) %>% paste(collapse = " ") # # extract coefficients
  model_eqn_beh <- gsub("^", "BEHAVIORAL> ", model_eqn_beh) # add a line for BEHAVIORAL>
  model_eqn_bim <- gsub("^", "EQ> ", model_eqn_bim) # add EQ> to start of line
  model_eqn_coe <- gsub("^", "COEFF> ", model_eqn_coe) # add COEFF> to start of line
  
  # model_eqn_bim <- gsub("([_[:alpha:]]+)(\\.)([[:digit:]]{1})", "TSLAG(\\1\\, \\3\\)", model_eqn_bim) # replace dot notation for lags with TSLAG()
  # model_eqn_beh <- str_extract(model_eqn_bim, "\\w*") %>% gsub("([_[:alpha:]]+)(_QL)", "\\1_Q", .) # extract the target variable name
  # model_eqn_bim <- gsub("([_[:alpha:]]+)(_QL)", "LOG(\\1_Q)", model_eqn_bim) # replace _QL with LOG()
  # model_eqn_coe <- str_extract_all(model_eqn_bim, "(b[[:digit:]]+)", simplify = TRUE) %>% paste(collapse = " ") # # extract coefficients
  # model_eqn_beh <- gsub("^", "BEHAVIORAL> ", model_eqn_beh) # add a line for BEHAVIORAL>
  # model_eqn_bim <- gsub("^", "EQ> ", model_eqn_bim) # add EQ> to start of line
  # model_eqn_coe <- gsub("^", "COEFF> ", model_eqn_coe) ## add COEFF> to start of line
  
  return(c(model_eqn, model_eqn_beh, model_eqn_bim, model_eqn_coe))
}

#' Parse gets() output and extract underlying data (GETS model development)
#'
#' @param model_in a model estimated by arx(), isat(), or getsm()
#'
#' @return an xts containing the model variables
#' @export
#'
#' @examples
#' save the data associated with a gets model
extract_data <- function(model_in) {
  data_out <- eviews(model_in, print = FALSE, return = TRUE)$data %>%
    select(-c) %>%
    rename(!!sym(yvar_name) := "y") %>%
    rename_with(~ str_replace(., "ar", str_glue("{yvar_name}."))) %>%
    rename_with(~ str_replace_all(., c("iis" = "IIS_", "sis" = "SIS_"))) %>%
    rename_with(~ str_replace_all(., "-", "_")) %>%
    ts_long() %>%
    ts_xts()
  
  return(data_out)
}


