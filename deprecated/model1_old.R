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

#### Data preprocessing ####
# Add weekend and weekday effect
df_long = df_long %>%
  mutate(weekNumber = lubridate::week(df_long$Date)) %>%
  mutate(weekend = lubridate::wday(df_long$Date, label = TRUE)
         %in% c("Sat", "Sun") %>% as.integer %>% as.factor)

X_regressors = c("weekend", "weekNumber", "medianAge") # Multicolinearity regions and medianAge

#### Run models ####
lag = 5 # Incubation period

# Construct formula
fm = paste("Confirmed ~ ",
           glue("lag(Confirmed, {lag}):lag(susceptibleRate, {lag})+"),
           paste(X_regressors, collapse="+")) %>%
  # paste("+factor(Code)") %>%
  as.formula

# Run model - Note: now pools all observations
lsdv = lm(fm, data=df_long)
summary(lsdv)

png(glue("{output_path}/model1_lag{lag}_lmplot.png"))
par(mfrow=c(2,2))
plot(lsdv)
par(mfrow=c(1,1))
dev.off()

# Note that df_long$incidenceRate[-1:-lag] may need to be replaced depending
# on whether the LSDV is run or a simple linear regression
residual = df_long$incidenceRate[-1:-lag] - lsdv$fitted.values

# Plot residuals
tibble(index = 1:length(residual), residuals = residual) %>%
  ggplot(aes(x=index, y=residual)) +
  theme(plot.title = element_text(face = "bold")) +
  geom_point(alpha=0.6, color='firebrick')
ggsave(glue("model1_lag{lag}_residuals.png"), path=output_path)

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

#### Plot alpha over time ####
df_meta = readxl::read_xlsx(path_wiki, sheet = "Metadata")
tbl = tibble(Date = as.Date(NA), Alpha=numeric(0), Code=character(0))
start = 31 # starting index - we want at least this number of observations

for (region in df_meta$Code){
  alphas = vector("double")
  dates = vector("character")
  data = df_long %>% filter(Code == !!region)
  for (t in start:nrow(data)){
    lsdv = lm(fm, data=data[1:t, ])
    alpha = lsdv$coefficients[["lag(Confirmed, 5):lag(susceptibleRate, 5)"]]
    alphas = c(alphas, alpha)
  }
  
  tbl = tbl %>%
    bind_rows(tibble(Date = data$Date[start:nrow(data)],
                     Alpha = alphas,
                     Code = region) %>%
                drop_na())
}

tbl = tbl %>%
  left_join(df_meta %>% select(c(Region, Code, Direction)), on=Code)

mean_alphas = vector("list")
days_smooth = 14
for (sub_tbl in split(tbl, tbl$Direction)){
  direc = sub_tbl$Direction[1]
  g = ggplot(sub_tbl, aes(Date, Alpha, colour = Region)) + 
    geom_point() +
    # geom_path(alpha=0.2) +
    geom_smooth(method="loess", span=0.3, se=FALSE) +
    # coord_cartesian(ylim = c(0, 1)) +
    xlab("") +
    ylab(TeX("$\\alpha_{within}$")) +
    scale_colour_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#0072B2",
                                 "#D55E00", "#CC79A7"))
  print(g)
  ggsave(glue("model1_lag{lag}_alphawithin_{direc}.png"), path=output_path)
  
  # Get the mean Alpha of the last `days_smooth` days per region
  mean_alphas[[direc]] = lapply(
    split(sub_tbl, sub_tbl$Code),
    function(x){mean(x$Alpha %>% tail(days_smooth))})
}

#### Burn-in period ####
data = df_long %>% filter(Code == "ABR")
ssr = vector()
for (burnin in 0:floor(nrow(data)/2)){
  lsdv = lm(fm, data=data[(1+burnin):nrow(data), ])
  ssr = c(ssr, sum(lsdv$residuals**2))
}

# TODO: Work on this; see a lot of difference in model specification per region
for (region in df_meta$Code){
  print(glue("#############{region}#############"))
  data = df_long %>% filter(Code == !!region)
  print(step(lm(fm, data=data), direction = "both",
             k=log(nrow(data)), trace=0))
}
