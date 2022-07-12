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

# input file name
# csv_file_in <- "sgd_ts_data_RADAquaMixDif.csv"
csv_file_in <- "Kiholo_seldata.csv"

# *************************
#  load data ----
# *************************

# load the time series
in_tbl <- read_csv(here("input", csv_file_in)) %>%
  # mutate(across(.cols = 1, .fns = ~ force_tz(., tzone = ""))) %>%
  mutate(across(.cols = 1, .fns = ~ parse_date_time(., c("ymdHMS", "mdyHM", "mdyHMS")))) %>%
  mutate(across(.cols = -1, .fns = as.numeric))

# *************************
#  calculations ----
# *************************

# creates data table
dat_tbl <- in_tbl %>%
  mutate(
    
    # when a single value is missing interpolate linearly
    Rn_wat__Bqm3 = if_else(is.na(Rn_wat__Bqm3), (lag(Rn_wat__Bqm3) + lead(Rn_wat__Bqm3))/2, Rn_wat__Bqm3),
    temp_wat__C = if_else(is.na(temp_wat__C), (lag(temp_wat__C) + lead(temp_wat__C))/2, temp_wat__C),
    sal_wat = if_else(is.na(sal_wat), (lag(sal_wat) + lead(sal_wat))/2, sal_wat),
    
    # calculates coastal radon measurement time interval in minutes based on provided measurement times, 
    # all other time series parameters provided by the user should be averaged to this interval
    meas_t__min = (time %>% as.numeric() - lag(time %>% as.numeric())) / 60,
    
    # EXPLAIN
    diff_owl__m = depth__m - lag(depth__m),

    # if radon mixing losses are determined via an independent method (current measurements, residence time estimates)
    # then f_mix_exp__Bqm2hr should be provided in the spreadsheet
    # otherwise, mixing losses will be estimated from the Rn mass balance
    f_mix_exp__Bqm2hr = if (!(f_mix_exp__Bqm2hr %>% is.null()) & (!(f_mix_exp__Bqm2hr %>% is.na())) %>% any()) {
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
 
    # excess Rn in water is calculated by sutracting dissolved 226Ra in water
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
  select(where(~(!(.x %>% is.na())) %>% any()))

# results saved in a csv file
write_csv(dat_tbl, here("output", "rn_budget.csv"), na = "")
  
# END OF RADON BUDGET CALCULATION

# *************************
#  diagnostics ----
# *************************

# interactive time series plot of two variables
dat_tbl %>% 
  ts_long() %>%
  ts_xts() %>%
  extract(, c("f_gw__m3m2d", "depth__m")) %T>% 
  {ser_summary <<- ts_summary(.)} %>% 
  plot_2(rng_start = ser_summary %>% pull(start) %>% min(), 
         rng_end = ser_summary %>% pull(end) %>% max(), 
         height = 300, width = 600)

# find start of second continuous chunk
start_2 <- dat_tbl %>% 
  ts_long() %>%
  filter(id == "f_gw__m3m2d") %>% 
  mutate(value = is.na(value) - is.na(value) %>% lag(default = 0)) %>% 
  filter(value %>% as.logical()) %>% 
  slice(n()) %>% 
  pull(time)

# # Compute features
# dat_features <- dat_tbl %>% 
#   drop_na() %>% 
#   ts_long() %>%
#   ts_tsibble() %>% 
#   fabletools::features(value, fabletools::feature_set(pkgs="feasts"))
# # fabletools::features(value, fabletools::feature_set(tags = "autocorrelation"))

# interactive time series plot of two variables
dat_tbl %>% 
  ts_long() %>%
  ts_xts() %>%
  extract(, c("f_gw__m3m2d", "depth__m")) %T>% 
  {ser_summary <<- ts_summary(.)} %>% 
  plot_2(rng_start = ser_summary %>% pull(start) %>% min(), 
         rng_end = ser_summary %>% pull(end) %>% max(), 
         height = 300, width = 600)
  
# proportional area plot of flux components
plot_out <- dat_tbl %>% 
  select(time, f_Rn_flood__Bqm2hr, f_Rn_atm__Bqm2hr, f_Rn_ebb__Bqm2hr, f_Rn_mix__Bqm2hr, f_Rn_gross__Bqm2hr) %>% 
  filter(time >= start_2) %>% 
  filter(time >= ymd("2015-04-05") & time < ymd("2016-04-08")) %>% 
  drop_na() %>% 
  ts_long() %>%
  mutate(id = factor(id, levels = c("f_Rn_flood__Bqm2hr", "f_Rn_atm__Bqm2hr", "f_Rn_ebb__Bqm2hr", "f_Rn_mix__Bqm2hr", "f_Rn_gross__Bqm2hr"))) %>%
  ts_regular() %>% 
  mutate(value = abs(value)) %>% 
  group_by(time) %>% 
  mutate(perc = value / sum(value)) %>% 
  ungroup() %>% 
  ggplot(aes(x = time, y = perc, fill = id)) + 
  geom_area(alpha = 1) + 
  scale_x_datetime("Time") + 
  scale_y_continuous("Rn Flux Components") + 
  theme_classic()
ggsave(here("output", "stack1.pdf"), plot_out)

# stacked area plot of flux components
plot_out <- dat_tbl %>% 
  select(time, f_Rn_flood__Bqm2hr, f_Rn_atm__Bqm2hr, f_Rn_ebb__Bqm2hr, f_Rn_mix__Bqm2hr, f_Rn_gross__Bqm2hr) %>% 
  filter(time >= start_2) %>% 
  filter(time >= ymd("2015-04-05") & time < ymd("2015-04-08")) %>% 
  drop_na() %>% 
  ts_long() %>% 
  mutate(id = factor(id, levels = c("f_Rn_flood__Bqm2hr", "f_Rn_atm__Bqm2hr", "f_Rn_ebb__Bqm2hr", "f_Rn_mix__Bqm2hr", "f_Rn_gross__Bqm2hr"))) %>%
  ts_regular() %>% 
  mutate(value = abs(value)) %>% 
  ggplot(aes(x = time, y = value, fill = id)) + 
  geom_area(alpha = 1) + 
  scale_x_datetime("Time") + 
  scale_y_continuous("Rn Flux Components", limits = c(0, 1000)) + 
  theme_classic()
ggsave(here("output", "stack2.pdf"), plot_out)

# *************************
#  spectral analysis ----
# *************************

# get the periodogram of individual series
spctrm <- dat_tbl %>% 
  ts_long() %>%
  filter(time >= start_2) %>% 
  ts_xts() %>%
  extract(, c("f_gw__m3m2d")) %>%
  # spec.pgram(log = "no") %>% 
  spectrum(log = "no") %>% 
  inset2("per", 1/.$freq / 60 / 60)

spec_df <- data.frame(freq = spctrm$freq, spec = spctrm$spec, per = spctrm$per)
# Create a vector of periods to label on the graph, units are in years
hrs_period <- rev(c(6, 24, 72, 672, 4380))
hrs_labels <- rev(c("6 hours", "1 day", "3 days", "28 days", "6 months"))
ggplot(data = spec_df) + 
  geom_line(aes(x = per, y = spec)) + 
  scale_x_log10("Period", breaks = hrs_period, labels = hrs_labels) +
  # scale_x_continuous("Period (hours)") + 
  scale_y_continuous("Spectrum") + 
  theme_minimal()

# most important periods
tibble(per = spctrm$per, spec = spctrm$spec) %>% 
  arrange(desc(spec)) %>% 
  print(n = 50)

# # decomposition of the time series
# # https://www.r-bloggers.com/2013/09/wheres-the-magic-emd-and-ssa-in-r/
# # library (Rssa)
# # library (EMD)
# # library (hht)
# 
# # singular spectrum analysis
# dat_ssa <- Rssa::ssa(dat_tbl %>% 
#                  ts_long() %>%
#                    filter(time >= start_2) %>% 
#                    ts_ts() %>%
#                  extract(, c("f_gw__m3m2d")))
# 
# plot(dat_ssa)
# summary(dat_ssa)
# plot(dat_ssa, type = "series", groups = 1:16)
# 
# dat_rec <- Rssa::reconstruct(dat_ssa, list(C1=1, C2=2:9, C3=10:11, C4=12:13))
# plot(dat_rec)
# plot(dat_tbl %>% 
#        ts_long() %>%
#        filter(time >= start_2) %>% 
#        ts_ts() %>%
#        extract(, c("f_gw__m3m2d")), col="gray")
# lines(dat_rec$C1)
# 
# # empirical mode decomposition
# ee <- EMD::emd(dat_tbl %>% filter(time >= start_2) %>% pull("f_gw__m3m2d"), 
#                 dat_tbl %>% filter(time >= start_2) %>% pull("time") %>% as.numeric())
# 
# plot(dat_tbl %>% filter(time >= start_2) %>% pull("time"), ee$imf[,1], type="l")
# lines(dat_tbl %>% filter(time >= start_2) %>% select("time", "f_gw__m3m2d"), col="red")
# plot(dat_tbl %>% filter(time >= start_2) %>% pull("time"), ee$imf[,2], type="l")
# lines(dat_tbl %>% filter(time >= start_2) %>% select("time", "f_gw__m3m2d"), col="red")
# plot(dat_tbl %>% filter(time >= start_2) %>% pull("time"), ee$imf[,3], type="l")
# lines(dat_tbl %>% filter(time >= start_2) %>% select("time", "f_gw__m3m2d"), col="red")
# 
# # ensemble empirical mode decomposition (takes a long time!)
# ee <- hht::EEMD(dat_tbl %>% pull("f_gw__m3m2d"), 
#            dat_tbl %>% pull("time") %>% as.numeric(), 100, 100, 6, "trials")
# eec <- hht::EEMDCompile ("trials", 100, 6)


# *************************
# EXPERIMENTAL CODE BELOW ----
# *************************


# *************************
#  regressions ----
# *************************

# QUESTION: should the variables be standardized?

# consider only variables that are non-NA
dat_est <- dat_tbl %>% 
  select(where(~ !(is.na(.x) %>% all()))) %>% 
  select(where(~ (sum(.x - lag(.x), na.rm = TRUE) != 0))) %>% 
  # filter(time >= start_2) %>% 
  # select(wind__ms, depth__m, winddir, waveht__m, f_gw__m3m2d) %>% 
  select(-1) %>%
  mutate(across(everything(), ~ scale(.x) %>% extract(,1)))

# # best subset regression
# # library(lmSubsets)
# training_top <- lmSubsets::lmSelect(f_gw__m3m2d ~ ., data = dat_est, nbest = 20, penalty = "AIC", include = NULL, exclude = NULL)
# 
# summary(training_top)
# image(training_top)
# plot(training_top)
# 
# training_top %>% formula(best = 1)
# training_top %>% variable.names(best = 1)
# 
# training_best <- training_top %>% 
#   formula(best = 8) %>% 
#   as.formula() %>% 
#   lm(data = dat_est)
# # refit()
# broom::tidy(training_best)
# broom::glance(training_best)
# # augment(training_best)
# training_best_coef <- broom::tidy(training_best) %>% 
#   mutate(conf.low = estimate - 2 * std.error, conf.high = estimate + 2 * std.error) #%>% 
# # arrange(term)
# training_best_coef_tab <- training_best_coef %>% 
#   select(-contains("conf"))
# 
# # library(car)
# pdf(here("output", "partial_plot_state.pdf"), width = 11, height = 8.5)
# # partial_plot_data <- car::avPlots(lm(f_gw__m3m2d ~ ., data = dat_est), 
# #                                   main=paste("Partial Regression Plot")
# # )
# partial_plot_data <- car::avPlots(lm(as.formula(training_best$terms), data = dat_est), 
#                                   main=paste("Partial Regression Plot")
# )
# dev.off()

# cor_mat <- dat_est %>% select(-time) %>% drop_na() %>% cor() %>% as_tibble() %>% round(2)
cor_mat <- dat_est %>% drop_na() %>% cor()
corrplot::corrplot.mixed(cor_mat, lower = "ellipse", upper = "number")
cor_mat <- cor_mat %>% as_tibble() %>% round(2)
rownames(cor_mat) <- colnames(cor_mat)
cor_mat <- cor_mat %>% rownames_to_column()
write_csv(cor_mat, here("output", "cor_mat.csv"), )

# scatter plot matrix
dat_est %>% names()
dat_scat_1 <- dat_est %>% 
  filter(!is.na(Rn_wat__Bqm3)) %>% 
  # select(atm = f_Rn_atm__Bqm2hr,kw_air, Rn_wat = Rn_wat__Bqm3, temp = temp_wat__C, sal = sal_wat, wind = wind__ms, owl = depth__m, wave = waveht__m, Rn_gross = f_Rn_gross__Bqm2hr)
  select(Rn_gross = f_Rn_gross__Bqm2hr, atm = f_Rn_atm__Bqm2hr, flood = f_Rn_flood__Bqm2hr, ebb = f_Rn_ebb__Bqm2hr, mix = f_Rn_mix__Bqm2hr, kw_air, Rn_wat = Rn_wat__Bqm3, temp = temp_wat__C, sal = sal_wat, wind = wind__ms, owl = depth__m, dowl = diff_owl__m, wave = waveht__m)
# pdf(here("output", "scatter_plot_matrix_1a.pdf"), width = 11, height = 8.5)
# car::scatterplotMatrix(dat_scat_1, smooth = FALSE)
# dev.off()
# pdf(here("output", "scatter_plot_matrix_1b.pdf"), width = 11, height = 8.5)
jpeg(here("output", "scatter_plot_matrix_1b.jpg"), width = 680, height = 680, quality = 95)
PerformanceAnalytics::chart.Correlation(dat_scat_1, histogram=TRUE, pch=19)
dev.off()
# pdf(here("output", "scatter_plot_matrix_1c.pdf"), width = 11, height = 8.5)
# corrplot::corrplot.mixed(dat_scat_1 %>% drop_na() %>% cor(), lower = "ellipse", upper = "number", tl.pos = "n")
# dev.off()

reg_out <- lm(Rn_gross ~ ., data = dat_scat_1 %>% select(-(atm:mix)))
reg_out %>% broom::tidy()

pdf(here("output", "partial_plot_gross.pdf"), width = 11, height = 8.5)
partial_plot_data <- car::avPlots(reg_out,
                                  main=paste("Partial Regression Plot")
)
dev.off()

reg_out <- lm(atm ~ ., data = dat_scat_1 %>% select(-c(Rn_gross, flood:mix)))
reg_out %>% broom::tidy()

pdf(here("output", "partial_plot_atm.pdf"), width = 11, height = 8.5)
partial_plot_data <- car::avPlots(reg_out,
                                  main=paste("Partial Regression Plot")
)
dev.off()

reg_out <- lm(mix ~ ., data = dat_scat_1 %>% select(-c(Rn_gross:ebb)))
reg_out %>% broom::tidy()

pdf(here("output", "partial_plot_mix.pdf"), width = 11, height = 8.5)
partial_plot_data <- car::avPlots(reg_out,
                                  main=paste("Partial Regression Plot")
)
dev.off()

# dat_scat_2 <- dat_est %>% 
#   filter(!is.na(Rn_wat__Bqm3)) %>% 
#   select(Rn_flood = f_Rn_flood__Bqm2hr, Rn_wat = Rn_wat__Bqm3, temp = temp_wat__C, sal = sal_wat, wind = wind__ms, owl = depth__m, wave = waveht__m, Rn_gross = f_Rn_gross__Bqm2hr)
# jpeg(here("output", "scatter_plot_matrix_2b.jpg"), width = 680, height = 680, quality = 95)
# PerformanceAnalytics::chart.Correlation(dat_scat_2, histogram=TRUE, pch=19)
# dev.off()
# 
# dat_scat_3 <- 
#   dat_est %>% 
#   filter(!is.na(Rn_wat__Bqm3)) %>% 
#   select(Rn_ebb = f_Rn_ebb__Bqm2hr, Rn_wat = Rn_wat__Bqm3, temp = temp_wat__C, sal = sal_wat, wind = wind__ms, owl = depth__m, wave = waveht__m, Rn_gross = f_Rn_gross__Bqm2hr)
# jpeg(here("output", "scatter_plot_matrix_3b.jpg"), width = 680, height = 680, quality = 95)
# PerformanceAnalytics::chart.Correlation(dat_scat_3, histogram=TRUE, pch=19)
# dev.off()
# 
# dat_scat_4 <- dat_est %>% filter(!is.na(Rn_wat__Bqm3)) %>% select(f_Rn_mix__Bqm2hr, Rn_wat__Bqm3, temp_wat__C, sal_wat, wind__ms, winddir, depth__m, waveht__m, f_Rn_gross__Bqm2hr)
# pdf(here("output", "scatter_plot_matrix_4a.pdf"), width = 11, height = 8.5)
# car::scatterplotMatrix(dat_scat_4, smooth = FALSE)
# dev.off()
# 
# dat_scat_5 <- dat_est %>% filter(!is.na(Rn_wat__Bqm3)) %>% select(f_Rn_total__Bqm2hr, f_Rn_gross__Bqm2hr, f_Rn_atm__Bqm2hr, f_Rn_flood__Bqm2hr, f_Rn_ebb__Bqm2hr, f_Rn_mix__Bqm2hr)
# pdf(here("output", "scatter_plot_matrix_5a.pdf"), width = 11, height = 8.5)
# car::scatterplotMatrix(dat_scat_5, smooth = FALSE)
# dev.off()

# *************************
#  shrinkage/penalized methods ----
# *************************

# # dat1 <- dat_est %>% drop_na() %>% select(-time) %>% scale() %>% as_tibble()
# dat1 <- dat_est %>% drop_na() %>% scale() %>% as_tibble()
# 
# # obtain automatic report about the PCA
# # library(FactoMineR)
# # library(FactoInvestigate)
# 
# res <- FactoMineR::PCA(dat1, graph=FALSE, scale.unit = TRUE, ncp = 5)
# FactoInvestigate::Investigate(res, file = here("output", "FactoInvestigateState.Rmd"), document = c("html_document"), keepRmd = TRUE)
# dimdesc(res, axes = 1:4, proba = 0.05)
# corr_table <-  dimdesc(res, axes = 1:1)$Dim.1$quanti %>% as.data.frame() %>% rownames_to_column() %>% as_tibble %>% rename_with(~c('Variable Name', 'Correlation', 'p-value'))
# res$var
# res$ind
# res$ind$coord
# res$ind$coord[,1] %>% scale()
# 
# dat1reg <- dat1 %>%
#   mutate(PC1 = res$ind$coord[, 1] %>% scale()) %>%
#   scale() %>%
#   as.data.frame()
# lm(f_gw__m3m2d ~ ., dat1reg) %>% summary()
# cor_w_HI <- cor(dat1reg)
# # scatter plot matrix
# pdf(here("output", "scatter_plot_matrix4.pdf"), width = 11, height = 8.5)
# car::scatterplotMatrix(dat1reg, smooth = FALSE, id = TRUE)
# dev.off()
# 
# cor_out1 <- cor(dat1reg) %>% round(2) %>% #select(-contains("Intercept"))
#   kableExtra::kable(format = "latex", booktabs = T, caption = "Correlation Matrix") %>% 
#   kableExtra::footnote(general = "Here is a general comments of the table.", threeparttable = T)
# 
# cor_w_HI[cor(dat1reg) %>% round(2) %>% lower.tri()] <- cor(dat1reg)[cor(dat1reg) %>% round(2) %>% lower.tri()]
# cor_out2 <- cor_w_HI %>% round(2) %>% #select(-contains("Intercept"))
#   kableExtra::kable(format = "latex", booktabs = T, caption = "Correlation Matrix") %>% 
#   # kableExtra::row_spec(3:5, background = "grey")
#   kableExtra::footnote(general = "Here is a general comments of the table.", threeparttable = T)
# 
# # # run PCA
# # res <- prcomp(dat1)
# # summary(res)
# # res$rotation
# # res$x
# # as.matrix(dat1) %*% res$rotation - res$x
# # res$x[,1]/res$sd[1]
# 
# # partial least squares (need to load package, namespace reference not enough)
# library(pls)
# 
# pls.model <- plsr(f_Rn_atm__Bqm2hr ~ ., data = dat_scat_1, validation = "CV")
# 
# # Find the number of dimensions with the lowest cross validation error
# cv <- RMSEP(pls.model)
# best.dims <- which.min(cv$val[estimate = "adjCV", , ]) - 1
# plot(RMSEP(pls.model), legendpos = "topright")
# plot(pls.model, plottype = "validation")
# 
# selectNcomp(pls.model, method = "onesigma", plot = TRUE)
# selectNcomp(pls.model, method = "randomization", plot = TRUE)
# temp <- crossval(pls.model, segments = 10)
# 
# # Rerun the model
# pls.model <- plsr(f_Rn_atm__Bqm2hr ~ ., data = dat_scat_1 %>% select(-c()), ncomp = best.dims)
# pls.model <- plsr(f_Rn_atm__Bqm2hr ~ ., data = dat_scat_1 %>% select(-c()), ncomp = best.dims, validation = "CV", jackknife = TRUE)
# temp <- coefplot(pls.model, se.whiskers = TRUE, labels = "names")
# pls.model <- plsr(f_Rn_atm__Bqm2hr ~ ., data = dat_scat_1 %>% select(-c()), ncomp = 4)
# 
# pls.model$validation
# pls.model %>% broom::tidy()
# 
# coefficients <- coef(pls.model)
# sum.coef <- sum(sapply(coefficients, abs))
# coefficients <- coefficients * 100 / sum.coef
# coefficients <- sort(coefficients[, 1 , 1])
# barplot(tail(coefficients, 6))
# barplot((coefficients))
# 
# scores(pls.model)
# temp <- pls::loadings(pls.model)
# str(temp)
# temp[[1]]
# 
# plot(pls.model, plottype = "coef", ncomp=1:7, legendpos = "bottomleft", labels = "names", xlab = "variables")
# plot(pls.model, ncomp = 7, asp = 1, line = TRUE)
# plot(pls.model, plottype = "scores", comps = 1:7)
# explvar(pls.model)
# plot(pls.model, plottype = "loadings", comps = 1:7, labels = "names", xlab = "variables")
# abline(h = 0)
# plot(pls.model, plottype = "correlation", comps = 1:4)
# plot(pls.model, plottype = "validation")
# loading.weights(pls.model)
# 
# 
# library(caret)
# 
# # Split the data into training and test set
# set.seed(123)
# training.samples <- dat1$f_gw__m3m2d %>%
#   createDataPartition(p = 0.8, list = FALSE)
# train.data  <- dat1[training.samples, ]
# test.data <- dat1[-training.samples, ]
# 
# 
# library(glmnet)
# 
# # Predictor variables
# x <- model.matrix(f_gw__m3m2d ~., train.data)[,-1]
# # Outcome variable
# y <- train.data$f_gw__m3m2d
# 
# # Ridge regression
# # Find the best lambda using cross-validation
# set.seed(123) 
# cv <- cv.glmnet(x, y, alpha = 0)
# # Display the best lambda value
# cv$lambda.min
# 
# # Fit the final model on the training data
# model_RR <- glmnet(x, y, alpha = 0, lambda = cv$lambda.min)
# # Display ridge regression coefficients
# coef(model_RR)
# 
# # LASSO
# glmnet(x, y, alpha = 1, lambda = NULL)
# 
# # Find the best lambda using cross-validation
# set.seed(123) 
# cv <- cv.glmnet(x, y, alpha = 1)
# # Display the best lambda value
# cv$lambda.min
# 
# # Fit the final model on the training data
# model_LASSO <- glmnet(x, y, alpha = 1, lambda = cv$lambda.min)
# # Dsiplay regression coefficients
# coef(model_LASSO)
# summary(model_LASSO)
# coef_table <- left_join(broom::tidy(model_RR) %>% select(term, estimate), broom::tidy(model_LASSO) %>% select(term, estimate), by = "term", suffix = c("RIDGE", "LASSO")) %>% 
#   left_join(bind_cols(coef(pls.model) %>% rownames(), coef(pls.model) %>% as_tibble()), by = c("term" = "...1")) %>% 
#   rename("estimatePLS" = "f_gw__m3m2d.11 comps") %>% 
#   left_join(broom::tidy(training_best)) %>% 
#   rename("estimateBEST" = "estimate")
# 
# # check weight of 100% filter out intercept.
# coef_sums <- coef_table %>% filter(term != "(Intercept)") %>% select(contains("estimate")) %>% 
#   abs() %>% colSums(na.rm = TRUE)
# coef_sums <- rep(1,4)
# colMeans(dat1)
# options(knitr.kable.NA = '')
# coef_table %>% mutate(estimateRIDGE = estimateRIDGE / coef_sums[1],
#                       estimateLASSO = estimateLASSO / coef_sums[2],
#                       estimatePLS = estimatePLS / coef_sums[3],
#                       estimateBEST = estimateBEST / coef_sums[4]) %>% 
#   filter(term != "(Intercept)") %>% 
#   select(-c(std.error, statistic)) %>% 
#   arrange(desc(estimateRIDGE)) %>%
#   mutate(across(where(is.numeric), function(x) round(x, 2))) %>%
#   kableExtra::kable(format = "latex", booktabs = T, caption = "Estimates", col.names = c("Predictor", "Ridge", "Lasso", "PLS", "Best", "p-value")) %>% 
#   # kableExtra::add_header_above(c(" ", "Estimate" = 3, "t-statistic" = 3, "p-value" = 3)) %>% 
#   kableExtra::footnote(general = "Here is a general comments of the table.", threeparttable = T)
# 
# 
#   
#   
#   
# # *************************
# # TEST CODE BELOW THIS LINE ----
# # *************************
# 
#   # spectal decomposition
#   # seasonal decomposition
#   # rising, falling tide
#   # correlogram lags, leads
#   # precip and gw in single time plot
#   # Rn and sal.ctd in single time plot
#   # Rn and water.owl in single time plot zoom into 1 week
#   
# 
# # don't use:
#   
#   
#   library(fpp3)
#   aus_cafe <- aus_retail %>%
#     filter(
#       Industry == "Cafes, restaurants and takeaway food services",
#       year(Month) %in% 2004:2018
#     ) %>%
#     summarise(Turnover = sum(Turnover))
#   
#   dat_harm <- dat_tbl %>% 
#     ts_long() %>%
#     ts_tsibble() %>% 
#     filter(
#       id == "f_gw__m3m2d",
#     )
#   
#   fit <- model(dat_harm,
#                `K = 1` = ARIMA(value ~ fourier(K=1) + PDQ(0,0,0)),
#                `K = 2` = ARIMA(value ~ fourier(K=2) + PDQ(0,0,0)),
#                `K = 3` = ARIMA(value ~ fourier(K=3) + PDQ(0,0,0)),
#                `K = 4` = ARIMA(value ~ fourier(K=4) + PDQ(0,0,0)),
#                `K = 5` = ARIMA(value ~ fourier(K=5) + PDQ(0,0,0)),
#                `K = 6` = ARIMA(value ~ fourier(K=6) + PDQ(0,0,0))
#   )
#   
#   fit %>%
#     forecast(h = "2 years") %>%
#     autoplot(aus_cafe, level = 95) +
#     facet_wrap(vars(.model), ncol = 2) +
#     guides(colour = "none", fill = "none", level = "none") +
#     geom_label(
#       aes(x = yearmonth("2007 Jan"), y = 4250,
#           label = paste0("AICc = ", format(AICc))),
#       data = glance(fit)
#     ) +
#     labs(title= "Total monthly eating-out expenditure",
#          y="$ billions")
#   
#   
