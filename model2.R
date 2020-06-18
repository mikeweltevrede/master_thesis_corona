#### Model 2 - Weighted Within-region spread ####
# We start with a simple model ignoring effects across regions:
# I_rt = I_rt-tau*S_rt-tau*sum_k alpha_within^k*W_rt-tau^k + X_rt*delta + nu_rt

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
# Transform into proportions
regressors_prop = c("airPassengersArrived", "touristArrivals",
                    "dischargeRateDiabetes", "dischargeRateRespiratory",
                    "dischargeRateHypertension", "dischargeRateCancer",
                    "dischargeRateChd", "dischargeRatePneumonia",
                    "dischargeRateTB", "availableBeds",
                    "maritimePassengersDisembarked", "railTravelers")

make_prop = function(x, na.rm = FALSE) { x / sum(x, na.rm = na.rm) }
df_long = df_long %>% 
  group_by(Date) %>% 
  mutate_at(regressors_prop, make_prop) %>% 
  ungroup

# Add weekend and weekday effect
df_long = df_long %>%
  mutate(weekNumber = lubridate::week(df_long$Date)) %>%
  mutate(weekend = lubridate::wday(df_long$Date, label = TRUE)
         %in% c("Sat", "Sun") %>% as.integer %>% as.factor)

cor(df_long$airPassengersArrived, df_long$touristArrivals) # = 0.5835016

regressors = c("touristArrivals", "broadbandAccess", "dischargeRateDiabetes",
               "dischargeRateHypertension", "dischargeRateCancer",
               "dischargeRateChd", "dischargeRateTB", "availableBeds",
               "riskOfPovertyOrSocialExclusion")
base_vars = c("Date", "Code", "incidenceRate", "susceptibleRate", X_regressors)
X_regressors = c("weekend", "weekNumber, medianAge")

#### Run models ####
lag = 5 # Incubation period

# Construct formula
fm = paste("incidenceRate ~",
           
           # The Eurostat regressors are multiplied by the lagged Inc and S
           # TODO: Investigate more closely if some variables should not be
           # multiplied in this way
           paste0(glue("dplyr::lag(incidenceRate, {lag})",
                       ":dplyr::lag(susceptibleRate, {lag}):")) %>%
             paste0(glue("dplyr::lag({regressors}, {lag})")) %>%
             paste(collapse="+"), "+",
           
           # These include the weekend and weekNumber effect
           paste(X_regressors, collapse="+")) %>%
  # paste("+factor(Code)") %>%
  as.formula

# Run model
lsdv = lm(fm, data=df_long)
summary(lsdv)
alias(lsdv) # Looks at the linearly dependent terms

png(glue("{output_path}/model1_lmplot_lag_{lag}.png"))
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
ggsave(glue("model2_residuals_plot_lag_{lag}.png"), path=output_path)

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

# Backwards stepwise selection with BIC
step(lsdv, direction = "backward", k=log(nrow(df_long)))
# In order, it kicks out:
# 1. touristArrivals
# 2. weekNumber
# 3. dischargeRateChd
# 4. weekend
# 5. dischargeRateHypertension

# Remaining:
# broadbandAccess, dischargeRateDiabetes, dischargeRateCancer, dischargeRateTB,
# availableBeds, riskOfPovertyOrSocialExclusion, and the regional effects
