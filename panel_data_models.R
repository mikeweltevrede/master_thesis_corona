rm(list=ls())

#### Panel data models ####
# We estimate three panel data models: Pooled OLS, Fixed Effects, and Random
# Effects. This is done for the following discretized SIR model:
# I{rt} - I_{r,t-1} = beta*S_{r,t-tau}*I_{r,t-tau} - gamma*I_{r,t-1} + eta_{rt}

### Setup ####

# Import standard variables
# source("config.R")

# Import packages
library(plm)
library(tidyverse)

# Import data
path_full_long = "data/data_long_full.csv"
df_long = readr::read_csv(path_full_long, col_types = do.call(
  cols, list(date = col_date(format = "%Y-%m-%d"))))

#### Decide the parameters ####
# You need to adapt, if desired, the following parameters:
# lag: int, the latent period
# rolling: boolean, whether to apply a rolling window_size
# window_size: int, if using a rolling window_size, how large?
# form: str, the form of undocumented infections to model with (if any)

# Latent period; Incubation period has median value 5; latent period is
# estimated to be 2 days shorter: 5-3=2
tau = 3

# Do we want to use a rolling window_size, i.e. only use the most recent `window_size`
# observations?
rolling = TRUE
window_size = 100

if (rolling) {
  rolling_flag = "_rolling"
} else {
  rolling_flag = ""
}

# Determine if we want to model undocumented infectives and, if so, by which
# method. Note that infective_variable is the number of new cases, i.e. Delta X.
form = "Quadratic" %>%
  to_upper_camel_case

if (form %in% c("Linear", "Quadratic", "DownwardsVertex", "UpwardsVertex",
                "Cubic")){
  print(glue("#### Running models while modelling undocumented infections with",
             " the {form} functional form! ####"))
  infective_variable = glue("infectives{form}")
  undoc_flag = glue("_Undoc{form}")
  
} else if (form == ""){
  # Then do not use the undocumented infections modelling
  print("#### Running models WITHOUT modelling undocumented infections! ####")
  infective_variable = "infectives"
  undoc_flag = ""
  
} else {
  sprintf(paste("The variable `form` is %s but it should be one of %s.",
                "Choosing 'infectives' as the infectives variable, so *not*",
                "modelling undocumented infections."),
          form,
          paste(c("Linear", "Quadratic", "DownwardsVertex",
                  "UpwardsVertex", "Cubic", ""), collapse=", ")) %>%
    print
  
  infective_variable = "infectives"
  undoc_flag = ""
}

#### Data preparation ####
# infectivesRateTotal is the number of active infectives divided by the total population
# infectivesRate is the number of new cases divided by the total population
# susceptibleRate is the total number of susceptibles divided by the total population

df = df_long %>%
  group_by() %>%
  transmute(
    date = date,
    code = code,
    total := .data[[glue("infectivesRateTotal{form}")]],
    dependent := .data[[glue("infectivesRate{form}")]],
    SI := dplyr::lag(.data[[glue("susceptibleRate{form}")]], tau)*
      dplyr::lag(.data[[glue("infectivesRateTotal{form}")]], tau),
    Ilag := dplyr::lag(.data[[glue("infectivesRateTotal{form}")]], 1)) %>%
  ungroup() %>%
  drop_na()

if (rolling) {
  df = df %>% 
    group_by(code) %>% 
    slice(tail(row_number(), window_size)) %>%
    ungroup()
}

#### Formula ####
# I{r,t} - I_{r,t-1} = beta*S_{r,t-tau}*I_{r,t-tau} - gamma*I_{r,t-1} + eta_{rt}
fm = as.formula("dependent ~ -1 + SI + Ilag")

#### Pooled OLS ####
pols = plm(fm, data=df, model = "pooling", index = c("code", "date"))
summary(pols) # NB: gamma (Ilag) is positive but should be negative?

#### Fixed Effects ####
fe = plm(fm, data=df, model = "within", index = c("code", "date"))
summary(fe) # NB: beta (SI) is negative but should be positive?
fixef(fe) # Regional fixed effects

#### Random Effects ####
re = plm(fm, data=df, model = "random", index = c("code", "date"))
summary(re) # NB: gamma (Ilag) is positive but should be negative?

phtest(fe, re) # p-value < 2.2e-16 => Prefer FE

#### OLD ####
# df_long = readr::read_csv(path_full_long, col_types = do.call(
#   cols, list(date = col_date(format = "%Y-%m-%d")))) %>%
#   mutate(proportionTested = totalTested / totalPopulation) %>%
#   group_by(code) %>%
#   mutate(
#     infectivesTotal = cumsum(infectives),
#     infectivesRateTotal = infectivesTotal / totalPopulation,
#     infectivesLag = infectivesTotal - dplyr::lag(infectivesTotal, 1)) %>%
#   ungroup()
# 
# fm = paste("infectivesLag ~ -1+",
#            "dplyr::lag(infectivesRateTotal, 1):dplyr::lag(susceptibleRate, 1) + ",
#            "dplyr::lag(infectivesRateTotal, 2)") %>% # Delay for recovery
#   as.formula
# 
# re = plm(fm, data=df_long, index=c("code", "date"), model="random")
# fe = plm(fm, data=df_long, index=c("code", "date"), model="within")
# pooled = plm(fm, data=df_long, index=c("code", "date"), model="pooling")
# 
# summary(re)
# summary(fe)
# summary(pooled)
# 
# dates = unique(df_long$date)
# betas = numeric(0)
# gammas = numeric(0)
# 
# for (t in 100:length(dates)) {
#   fm = paste("infectivesLag ~ -1+",
#              "dplyr::lag(infectivesRateTotal, 1):dplyr::lag(susceptibleRate, 1)") %>% # Delay for recovery
#     as.formula
#   
#   pooled = plm(fm, data=filter(df_long, date <= dates[t]), index="code", model="pooling")
#   # gamma = coef(summary(pooled))["dplyr::lag(infectivesRateTotal, 2)",1]
#   beta = coef(summary(pooled))["dplyr::lag(infectivesRateTotal, 1):dplyr::lag(susceptibleRate, 1)",1]
#   
#   betas = c(betas, beta)
#   # gammas = c(gammas, gamma)
# }
# 
# plot(betas, type="l", col="blue")
# lines(-gammas, col = "red")
# 
# phtest(fe, re) # p-value < 2.2e-16 => Prefer FE
# 
# #### Make plot of R0 over time
# dates = unique(df_long$date)
# betas = numeric(0)
# gammas = numeric(0)
# Rs = numeric(0)
# 
# for (t in 30:length(dates)) {
#   data = filter(df_long, date <= dates[t])
#   
#   fe = plm(fm, data=data, index="code",
#               model="within")
#   re = plm(fm, data=data, index="code",
#            model="random")
#   test = phtest(fe, re)
#   
#   if (test$p.value <= 0.05) {
#     model = fe
#   } else {
#     model = re
#   }
#   
#   model_summ = summary(model)
#   
#   beta = coef(model_summ)[
#     "dplyr::lag(infectivesTotal, 1):dplyr::lag(susceptibleRate, 1)", "Estimate"]
#   gamma = -coef(model_summ)["dplyr::lag(infectivesTotal, 1)", "Estimate"]
#   
#   betas = c(betas, beta)
#   gammas = c(gammas, gamma)
#   Rs = c(Rs, beta/gamma)
# }
# 
# plot(betas, type="l", col="red")
# lines(gammas, col="blue")
# 
# plot(Rs, type="l")
# 
# sum((betas-gammas)^2)