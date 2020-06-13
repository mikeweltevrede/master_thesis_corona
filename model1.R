#### Model 1 - Within-region spread ####
# We start with a simple model ignoring effects across regions:
# I_rt = alpha_within*I_rt-tau*S_rt-tau + X_rt*delta + nu_rt

#### Setup ####
# source("clean_full.R") # May error; if so, run it by hand

# Import standard variables
source("config.R")

library(tidyverse)
library(glue)
library(lmtest)
library(aTSA)

df_long = readr::read_csv(path_full_long, col_types = do.call(
  cols, list(Date = col_date(format = "%Y-%m-%d"))))

#### Data preprocessing ####
# Add weekend and weekday effect
df_long = df_long %>%
  mutate(weekNumber = lubridate::week(df_long$Date)) %>%
  mutate(weekend = lubridate::wday(df_long$Date, label = TRUE)
         %in% c("Sat", "Sun") %>% as.integer %>% as.factor)

X_regressors = c("weekend", "weekNumber", "medianAge")
base_vars = c("Date", "Code", "incidenceRate", "susceptibleRate", X_regressors)

#### Run models ####
lag = 5 # Incubation period
  
# Construct formula
fm = paste("incidenceRate ~ ",
           glue("lag(incidenceRate, {lag}):lag(susceptibleRate, {lag})+"),
           paste(X_regressors, collapse="+")) %>%
  paste("+factor(Code)") %>%
  as.formula

# Run model
lsdv = lm(fm, data=df_long)
residual = df_long$incidenceRate[-1:-lag] - lsdv$fitted.values

# TODO: Replace DW test by another unit root test
dw_results = dwtest(fm, data=df_long)
dw_stat = dw_results$statistic
dw_pval = dw_results$p.value

t_results = t.test(residual)
t_stat = t_results$statistic
t_pval = t_results$p.value

# Test
lag_col_counts = table(lag_col)
index = vector()
for (len in lag_col_counts){
  index = c(index, 1:len)
}

tbl = tibble(index, residuals = residual, lag=lag_col, dw_col, dw_p_col,
             t_col, t_p_col)
tbl$lag = factor(tbl$lag, levels=unique(tbl$lag),
                 labels=paste0("Lag: ", lags,
                               "\n Durbin-Watson statistic: ",
                               signif(unique(tbl$dw_col), 5),
                               " (p=", signif(unique(tbl$dw_p_col), 5), ")",
                               "\n t-statistic for mean 0: ",
                               signif(unique(tbl$t_col), 5),
                               " (p=", signif(unique(tbl$t_p_col), 5), ")"))

ggplot(tbl, aes(x=index, y=resids_col)) +
  facet_wrap(~lag) +
  theme(plot.title = element_text(face = "bold")) +
  geom_point(alpha=0.6, color='firebrick')
ggsave("model1_residuals_plot_lag_1_6.png", path=output_path)
# One Sample t-test for zero-mean
t.test(residual) # p=1: true mean is not equal to 0

# Tests for autocorrelation
lmtest::dwtest(fm, data=df_long) # Durbin-Watson test: p=0.76: ac > 0
lmtest::bgtest(fm, data=df_long) # Breusch-Godfrey test: p<2.2e-16: ac > 0(?)
Box.test(residual, type = "Ljung-Box") # Ljung-Box test: p<2.2e-16: ac > 0

# Tests for stationarity
aTSA::stationary.test(residual, method = "adf") # ADF: stationary
aTSA::stationary.test(residual, method = "pp") # Phillips-Perron: stationary
aTSA::stationary.test(residual, method = "kpss") # KPSS: nonstationary
