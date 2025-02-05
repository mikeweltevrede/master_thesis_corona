#### Model 3 - Within and between-region spread ####
# I_rt = alpha_within*I_rt-tau*S_rt-tau +
#        alpha_within*S_rt-tau*\sum_{c!=r}I_ct-tau +
#        X_rt*delta + nu_rt

#### Setup ####
# Import standard variables
source("config.R")

library(tidyverse)
library(glue)
library(lmtest)
library(aTSA)

df_long = readr::read_csv(path_full_long, col_types = do.call(
  cols, list(date = col_date(format = "%Y-%m-%d"))))

df_wide = readr::read_csv(path_full_wide, col_types = do.call(
  cols, list(date = col_date(format = "%Y-%m-%d"))))

lag = 5 # Incubation period
codes = df_long$code %>% unique
df_sumInc = tibble(date = as.Date(NA),
                   code = character(),
                   sumInfectives = numeric())
for (region in codes){
  df_sumInc = df_sumInc %>%
    bind_rows(
      df_wide %>%
      select(map(codes[codes != region], starts_with, vars = colnames(.)) %>%
               unlist()) %>%
      select(ends_with("infectives")) %>%
      mutate_all(dplyr::lag, n=lag) %>%
      transmute(
        date = df_wide$date,
        code = region,
        sumInfectives = rowSums(.)))
}

rm(df_wide)

df_long = df_long %>%
  left_join(df_sumInc, by = c("code", "date"))
rm(df_sumInc)

#### Data preprocessing ####
# Add weekend and weekday effect
df_long = df_long %>%
  mutate(weekend = lubridate::wday(df_long$date, label = TRUE)
         %in% c("Sat", "Sun") %>% as.integer %>% as.factor)

X_regressors = c("weekend")

#### Run model 3 ####
# Construct formula
fm = paste("infectives ~ ",
           glue("lag(infectives, {lag}):lag(susceptibleRate, {lag})+"),
           glue("lag(susceptibleRate, {lag}):sumInfectives +"),
           paste(X_regressors, collapse="+")) %>%
  paste("+ factor(code)") %>%
  as.formula

# Run model
lsdv = lm(fm, data=df_long) # Multicollinearity regions and medianAge
summary(lsdv)

png(glue("{output_path}/model3_lag{lag}_lmplot.png"))
par(mfrow=c(2,2))
plot(lsdv)
par(mfrow=c(1,1))
dev.off()

# We remove the first `lag` dates for the residuals
df_long_sub = df_long[!df_long$date %in% (df_long$date %>% unique() %>% .[1:lag]), ]

# Compute the residuals
residual = df_long_sub$infectives - lsdv$fitted.values

# Plot residuals
tibble(index = 1:length(residual), residuals = residual) %>%
  ggplot(aes(x=index, y=residual)) +
  theme(plot.title = element_text(face = "bold")) +
  geom_point(alpha=0.6, color='firebrick')
ggsave(glue("model3_lag{lag}_residuals.png"), path=output_path)

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
