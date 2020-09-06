rm(list=ls())

#### Panel data models ####
# We estimate two panel data models: Pooled OLS and Random Effects. This is done
# for the following discretized SIR model:
# I{rt} - I_{r,t-1} = beta*S_{r,t-tau}*I_{r,t-tau} - gamma*I_{r,t-1} + eta_{rt}
# Fixed Effects is not estimated because the SIR model does not include an
# intercept, also not a regional effect. Since, in essence, this is what fixed
# effects takes into account, it does not make sense in this case.

#### Setup ####
# Import standard variables
source("config.R")

# Import packages
library(snakecase)
library(plm)
library(tidyverse)

# Import data
df_long = readr::read_csv(path_full_long, col_types = do.call(
  cols, list(date = col_date(format = "%Y-%m-%d"))))

#### Decide the parameters ####
# You need to adapt, if desired, the following parameters:
# tau: int, the latent period
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
  
  df_long = df_long %>% 
    group_by(code) %>% 
    slice(tail(row_number(), window_size)) %>%
    ungroup()
  
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
tau = 1 # TODO: Test both 1 and 3; dependent on the location!

df_long = df_long %>%
  group_by(code) %>%
  transmute(
    date = date,
    code = code,
    
    # Density-dependent
    X = .data[[glue("susceptiblePopulation{form}")]],
    Xlag = dplyr::lag(X, 1),
    dX = X - Xlag,
    Y = .data[[glue("infectivesTotal{form}")]],
    Ylag = dplyr::lag(Y, 1),
    dY = Y - Ylag,
    Z = cumsum(recovered) + cumsum(deaths),
    Zlag = dplyr::lag(Z, 1),
    dZ = Z - Zlag,
    
    # Frequency dependent
    S = X / totalPopulation,
    Slag = dplyr::lag(S, 1),
    dS = S - Slag,
    I = Y / totalPopulation,
    Ilag = dplyr::lag(I, 1),
    dI = I - Ilag,
    R = Z / totalPopulation,
    Rlag = dplyr::lag(R, 1),
    dR = R - Rlag) %>%
  ungroup() %>%
  drop_na()

# Density dependent
print("#### Density dependent (X, Y, Z) ####")
for (method in c("pooling", "random")) {
  print(glue("\n\n## {toupper(method)} ##"))
  
  #### Estimate beta and gamma directly ####
  model_s = plm(dX ~ -1 + Xlag:Ylag, data=df_long, model = method,
               index = c("code", "date"))
  model_r = plm(dZ ~ -1 + Ylag, data=df_long, model = method,
               index = c("code", "date"))
  
  beta = -unname(coef(model_s))
  gamma = unname(coef(model_r))
  
  if (beta < 0) {
    print("Warning: Beta is negative!")
  }
  if (gamma < 0) {
    print("Warning: Gamma is negative!")
  }
  
  print(glue("Gamma: {gamma}"))
  
  r = beta / gamma
  print(glue("Normal   | Beta: {beta}, R: {r}"))
  
  # Estimate beta with gamma as estimated above
  df_long = df_long %>%
    mutate(dependent := dY + !!gamma*Ylag)
  
  model_i = plm(dependent ~ -1 + Xlag:Ylag, data=df_long, model = method,
               index = c("code", "date"))
  
  beta_twostep = unname(coef(model_i))
  
  if (beta_twostep < 0) {
    print("Warning: Beta is negative!")
  }
  
  r_twostep = beta_twostep / gamma
  print(glue("Two-step | Beta: {beta_twostep}, R: {r_twostep}"))
}

# Frequency dependent
print("\n\n#### Frequency dependent (S, I, R) ####")
for (method in c("pooling", "random")) {
  print(glue("\n\n## {toupper(method)} ##"))
  
  #### Estimate beta and gamma directly ####
  model_s = plm(dS ~ -1 + Slag:Ilag, data=df_long, model = method,
                index = c("code", "date"))
  model_r = plm(dR ~ -1 + Ilag, data=df_long, model = method,
                index = c("code", "date"))
  
  beta = -unname(coef(model_s))
  gamma = unname(coef(model_r))
  
  if (beta < 0) {
    print("Warning: Beta is negative!")
  }
  if (gamma < 0) {
    print("Warning: Gamma is negative!")
  }
  
  print(glue("Gamma: {gamma}"))
  
  r = beta / gamma
  print(glue("Normal   | Beta: {beta}, R: {r}"))
  
  # Estimate beta with gamma as estimated above
  df_long = df_long %>%
    mutate(dependent := dI + !!gamma*Ilag)
  
  model_i = plm(dependent ~ -1 + Slag:Ilag, data=df_long, model = method,
                index = c("code", "date"))
  
  beta_twostep = unname(coef(model_i))
  
  if (beta_twostep < 0) {
    print("Warning: Beta is negative!")
  }
  
  r_twostep = beta_twostep / gamma
  print(glue("Two-step | Beta: {beta_twostep}, R: {r_twostep}"))
}

#### NLS? ####
removed = function(infectives, gamma) {
  gamma*infectives
}

nls(dY ~ -1 + Xlag:Ylag + removed(Ylag, gamma),
    data=filter(df_long, code == "LOM"), start = list(gamma = 0.01), trace = TRUE)

# Google: nls r with function
# https://rpubs.com/RobinLovelace/nls-function
