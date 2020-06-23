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
library(latex2exp)

df_long = readr::read_csv(path_full_long, col_types = do.call(
  cols, list(Date = col_date(format = "%Y-%m-%d"))))

df_wide = readr::read_csv(path_cleaned_wide, col_types = do.call(
  cols, list(Date=col_date(format="%Y-%m-%d"))))

df_wide_full = readr::read_csv(path_full_wide, col_types = do.call(
  cols, list(Date=col_date(format="%Y-%m-%d"))))

#### Data preprocessing ####
# Add nationwide sums
susceptible = df_wide_full %>%
  select(ends_with("susceptiblePopulation")) %>%
  rowSums
total = df_wide_full %>%
  select(ends_with("totalPopulation")) %>%
  rowSums
confirmed = df_wide_full %>% 
  select(ends_with("_Confirmed")) %>% 
  rowSums
df_wide_full = df_wide_full %>%
  mutate(susceptibleRate_total = susceptible/total,
         confirmed_total = confirmed)

# Add weekend and weekday effect
df_wide_full = df_wide_full %>%
  mutate(weekNumber = lubridate::week(df_wide_full$Date)) %>%
  mutate(weekend = lubridate::wday(df_wide_full$Date, label = TRUE)
         %in% c("Sat", "Sun") %>% as.integer %>% as.factor)

X_regressors = c("weekend", "weekNumber") # Multicolinearity regions and medianAge

#### Run models ####
lag = 5 # Incubation period

# Construct formula
fm = paste("confirmed_total ~ ",
           glue("lag(confirmed_total, {lag}):lag(susceptibleRate_total, {lag})+"), # TODO: Susceptible population nationwide
           paste(X_regressors, collapse="+")) %>%
  # paste("+factor(Code)") %>%
  as.formula

# Run model
lsdv = lm(fm, data=df_wide_full)
summary(lsdv)

png(glue("{output_path}/model1_lag{lag}_lmplot_national.png"))
par(mfrow=c(2,2))
plot(lsdv)
par(mfrow=c(1,1))
dev.off()

residual = df_wide_full$confirmed_total[-1:-lag] - lsdv$fitted.values

# Plot residuals
tibble(index = 1:length(residual), residuals = residual) %>%
  ggplot(aes(x=index, y=residual)) +
  theme(plot.title = element_text(face = "bold")) +
  geom_point(alpha=0.6, color='firebrick')
ggsave(glue("model1_lag{lag}_residuals_national.png"), path=output_path)

# One Sample t-test for zero-mean
t.test(residual) # p=1: true mean is equal to 0

# Tests for autocorrelation
lmtest::dwtest(fm, data=df_wide_full) # Durbin-Watson test: p-value < 2.2e-16: ac > 0
lmtest::bgtest(fm, data=df_wide_full) # Breusch-Godfrey test: p<2.2e-16: ac > 0 (?)
Box.test(residual, type = "Ljung-Box") # Ljung-Box test: p<2.2e-16: ac > 0

# Tests for stationarity
aTSA::stationary.test(residual, method = "adf") # ADF: (not sure about the results)
aTSA::stationary.test(residual, method = "pp") # Phillips-Perron: stationary
aTSA::stationary.test(residual, method = "kpss") # KPSS: nonstationary

#### Plot alpha over time ####
tbl = tibble(Date = as.Date(NA), Alpha=numeric(0))
start = 31 # starting index - we want at least this number of observations

alphas = vector("double")
dates = vector("character")
for (t in start:nrow(df_wide_full)){
  lsdv = lm(fm, data=df_wide_full[1:t, ])
  alpha = lsdv$coefficients[[glue("lag(confirmed_total, {lag}):lag(susceptibleRate_total, {lag})")]]
  alphas = c(alphas, alpha)
}

tbl = tibble(Date = df_wide_full$Date[start:nrow(df_wide_full)],
                   Alpha = alphas) %>%
  drop_na()

# Get the mean Alpha of the last `days_smooth` days per region
days_smooth = 14
mean_alphas = mean(tbl$Alpha %>% tail(days_smooth))

g= ggplot(tbl, aes(Date, Alpha)) + 
  geom_point(col="firebrick") +
  # geom_path(alpha=0.2) +
  geom_smooth(method="loess", span=0.3, se=FALSE, col="firebrick") +
  # coord_cartesian(ylim = c(0, 1)) +
  xlab("") +
  ylab(TeX("$\\alpha_{within}$"))
print(g)
ggsave(glue("model1_lag{lag}_alphawithin_national.png"), path=output_path)
