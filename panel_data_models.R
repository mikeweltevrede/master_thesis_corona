### Setup ####
rm(list=ls())

library(plm)
library(tidyverse)

# Import standard variables
source("config.R")

df_long = readr::read_csv(path_full_long, col_types = do.call(
  cols, list(date = col_date(format = "%Y-%m-%d")))) %>%
  mutate(proportionTested = totalTested / totalPopulation) %>%
  group_by(code) %>%
  mutate(
    infectivesTotal = cumsum(infectives),
    infectivesRateTotal = infectivesTotal / totalPopulation,
    infectivesLag = infectivesTotal - dplyr::lag(infectivesTotal, 1)) %>%
  ungroup()



fm = paste("infectivesLag ~ ",
           "dplyr::lag(infectivesTotal, 1):dplyr::lag(susceptibleRate, 1) + ",
           "dplyr::lag(infectivesTotal, 1)") %>%
  as.formula

re = plm(fm, data=df_long, index="code", model="random")
fe = plm(fm, data=df_long, index="code", model="within")
pooled = plm(fm, data=df_long, index="code", model="pooling")

summary(re)
summary(fe)
summary(pooled)

phtest(fe, re) # p-value < 2.2e-16 => Prefer FE

#### Make plot of R0 over time
dates = unique(df_long$date)
betas = numeric(0)
gammas = numeric(0)
Rs = numeric(0)

for (t in 30:length(dates)) {
  data = filter(df_long, date <= dates[t])
  
  fe = plm(fm, data=data, index="code",
              model="within")
  re = plm(fm, data=data, index="code",
           model="random")
  test = phtest(fe, re)
  
  if (test$p.value <= 0.05) {
    model = fe
  } else {
    model = re
  }
  
  model_summ = summary(model)
  
  beta = coef(model_summ)[
    "dplyr::lag(infectivesTotal, 1):dplyr::lag(susceptibleRate, 1)", "Estimate"]
  gamma = -coef(model_summ)["dplyr::lag(infectivesTotal, 1)", "Estimate"]
  
  betas = c(betas, beta)
  gammas = c(gammas, gamma)
  Rs = c(Rs, beta/gamma)
}

plot(betas, type="l", col="red")
lines(gammas, col="blue")

plot(Rs, type="l")

sum((betas-gammas)^2)