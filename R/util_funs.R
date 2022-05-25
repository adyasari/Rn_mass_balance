#
# Utility Functions
#
# Author: Peter Fuleky
# Date: May 24, 2022

# generate help files for these functions
# d <- document::document(file_name = here::here("R", "util_funs.R"), output_directory = here::here("man"))

#' Interactive plot of two variables (levels)
#'
#' @param ser time series to plot (xts format)
#' @param rng_start start of zoom range ("YYYY-MM-DD")
#' @param rng_end end of the zoom range ("YYYY-MM-DD")
#' @param height height of a single panel (px)
#' @param width width of a single panel (px)
#'
#' @return a dygraph plot
#' @export
#'
#' @examples
#' dat_tbl %>% 
#' ts_long() %>%
#'   ts_xts() %>%
#'   extract(, c("f_Rn_flood__Bqm2hr", "depth__m")) %T>% 
#'   {ser_summary <<- ts_summary(.)} %>% 
#'   plot_2(rng_start = ser_summary %>% pull(start) %>% extract(1), rng_end = ser_summary %>% pull(end) %>% extract(1), height = 300, width = 600)
plot_2 <- function(ser, rng_start = as.character(Sys.Date() - years(5)), rng_end = as.character(Sys.Date()), height = 300, width = 600) {
  ser_names <- ser %>%
    ts_xts() %>%
    names()
  ser_plot <- ser %>%
    ts_xts() %>%
    ts_tslist() %>%
    set_names(ser_names) %>%
    ts_dygraphs(main = "", height = height, width = width) %>%
    dygraphs::dyAxis("y", label = ser_names[1]) %>%
    dygraphs::dyAxis("y2", label = ser_names[2]) %>%
    dygraphs::dyLegend(width = width * 0.90) %>%
    dygraphs::dySeries(ser_names[1], axis = "y") %>%
    dygraphs::dySeries(ser_names[2], axis = "y2") %>%
    dygraphs::dyHighlight(highlightSeriesBackgroundAlpha = 0.2, hideOnMouseOut = FALSE) %>% 
    dygraphs::dyRangeSelector(dateWindow = c(rng_start, rng_end), height = 30, strokeColor = "red")
  ser_plot
}

