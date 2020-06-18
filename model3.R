#### Model 3 - Within and between-region spread ####
# I_rt = alpha_within*I_rt-tau*S_rt-tau +
#        alpha_within*S_rt-tau*\sum_{c!=r}I_ct-tau +
#        X_rt*delta + nu_rt

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

df_wide = readr::read_csv(path_full_wide, col_types = do.call(
  cols, list(Date = col_date(format = "%Y-%m-%d"))))

lag = 5 # Incubation period
codes = df_long$Code %>% unique
df_sumInc = tibble(Date = as.Date(NA),
                   Code = character(),
                   sumIncidenceRate = numeric())
for (region in codes){
  df_sumInc = df_sumInc %>%
    bind_rows(
      df_wide %>%
      select(map(codes[codes != region], starts_with, vars = colnames(.)) %>%
               unlist()) %>%
      select(ends_with("incidenceRate")) %>%
      mutate_all(dplyr::lag, n=lag) %>%
      transmute(
        Date = df_wide$Date,
        Code = region,
        sumIncidenceRate = rowSums(.)))
}

rm(df_wide)

df_long = df_long %>%
  left_join(df_sumInc, by = c("Code", "Date"))
rm(df_sumInc)

#### Data preprocessing ####
# Add weekend and weekday effect
df_long = df_long %>%
  mutate(weekNumber = lubridate::week(df_long$Date)) %>%
  mutate(weekend = lubridate::wday(df_long$Date, label = TRUE)
         %in% c("Sat", "Sun") %>% as.integer %>% as.factor)

X_regressors = c("weekend", "weekNumber", "medianAge")

#### Run model 3 ####
# Construct formula
fm = paste("incidenceRate ~ ",
           glue("lag(incidenceRate, {lag}):lag(susceptibleRate, {lag})+"),
           glue("lag(susceptibleRate, {lag}):sumIncidenceRate +"),
           paste(X_regressors, collapse="+")) %>%
  # paste("+ factor(Code)") %>%
  as.formula

# Run model
lsdv = lm(fm, data=df_long) # Multicolinearity regions and medianAge
summary(lsdv)
alias(lsdv) # Looks at the linearly dependent terms

png(glue("{output_path}/model3_lmplot_lag_{lag}.png"))
par(mfrow=c(2,2))
plot(lsdv)
par(mfrow=c(1,1))
dev.off()

residual = df_long$incidenceRate[-1:-lag] - lsdv$fitted.values

# Plot residuals
tibble(index = 1:length(residual), residuals = residual) %>%
  ggplot(aes(x=index, y=residual)) +
  theme(plot.title = element_text(face = "bold")) +
  geom_point(alpha=0.6, color='firebrick')
ggsave(glue("model3_residuals_plot_lag_{lag}.png"), path=output_path)

# One Sample t-test for zero-mean
t.test(residual) # p=0.7036: true mean is not equal to 0

# Tests for autocorrelation
lmtest::dwtest(fm, data=df_long) # Durbin-Watson test: p=0.7604: ac > 0
lmtest::bgtest(fm, data=df_long) # Breusch-Godfrey test: p=0.0007456: ac > 0(?)
Box.test(residual, type = "Ljung-Box") # Ljung-Box test: p=0.4148: ac not 0(?)

# Tests for stationarity
aTSA::stationary.test(residual, method = "adf") # ADF: stationary
aTSA::stationary.test(residual, method = "pp") # Phillips-Perron: stationary
aTSA::stationary.test(residual, method = "kpss") # KPSS: nonstationary
