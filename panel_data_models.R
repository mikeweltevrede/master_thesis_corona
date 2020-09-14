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
library(xtable)

# Import data
df_long = readr::read_csv(path_full_long, col_types = do.call(
  cols, list(date = col_date(format = "%Y-%m-%d"))))
df_wide = readr::read_csv(path_full_wide, col_types = do.call(
  cols, list(date = col_date(format = "%Y-%m-%d"))))

#### Decide the parameters ####
# You need to adapt, if desired, the following parameters:
# rolling: boolean, whether to apply a rolling window_size
# window_size: int, if using a rolling window_size, how large?
# form: str, the form of undocumented infections to model with (if any)

# Do we want to use a rolling window_size, i.e. only use the most recent `window_size`
# observations?
rolling = TRUE
window_size = 100

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
# Latent period; Incubation period has median value 5; latent period is
# estimated to be 2 days shorter: 5-3=2
tau = 3

deathsNational = df_wide %>% 
  select(ends_with("_deaths")) %>% 
  rowSums
recoveredNational = df_wide %>% 
  select(ends_with("_recovered")) %>% 
  rowSums

df_wide = df_wide %>%
  transmute(
    date = date,
    
    # Density-dependent
    X = .data[[glue("susceptiblePopulationNational{form}")]],
    XlagOne = dplyr::lag(X, 1),
    XlagTau = dplyr::lag(X, tau),
    dX = X - dplyr::lag(X, 1),
    Y = cumsum(.data[[glue("infectivesTotalNational{form}")]]),
    YlagOne = dplyr::lag(Y, 1),
    YlagTau = dplyr::lag(Y, tau),
    dY = Y - YlagOne,
    Z := cumsum(!!recoveredNational) + cumsum(!!deathsNational),
    Zlag = dplyr::lag(Z, 1),
    dZ = Z - Zlag,
    
    # Frequency dependent
    S = X / totalPopulationNational,
    SlagOne = dplyr::lag(S, 1),
    SlagTau = dplyr::lag(S, tau),
    dS = S - SlagOne,
    I = Y / totalPopulationNational,
    IlagOne = dplyr::lag(I, 1),
    IlagTau = dplyr::lag(I, tau),
    dI = I - IlagOne,
    R = Z / totalPopulationNational,
    Rlag = dplyr::lag(R, 1),
    dR = R - Rlag) %>%
  drop_na()

df_long = df_long %>%
  group_by(code) %>%
  transmute(
    date = date,
    code = code,
    
    # Density-dependent
    X = .data[[glue("susceptiblePopulation{form}")]],
    XlagOne = dplyr::lag(X, 1),
    XlagTau = dplyr::lag(X, tau),
    dX = X - dplyr::lag(X, 1),
    Y = .data[[glue("infectivesTotal{form}")]],
    YlagOne = dplyr::lag(Y, 1),
    YlagTau = dplyr::lag(Y, tau),
    dY = Y - YlagOne,
    Z = cumsum(recovered) + cumsum(deaths),
    Zlag = dplyr::lag(Z, 1),
    dZ = Z - Zlag,
    
    # Frequency dependent
    S = X / totalPopulation,
    SlagOne = dplyr::lag(S, 1),
    SlagTau = dplyr::lag(S, tau),
    dS = S - SlagOne,
    I = Y / totalPopulation,
    IlagOne = dplyr::lag(I, 1),
    IlagTau = dplyr::lag(I, tau),
    dI = I - IlagOne,
    R = Z / totalPopulation,
    Rlag = dplyr::lag(R, 1),
    dR = R - Rlag) %>%
  ungroup() %>%
  drop_na()


if (rolling) {
  rolling_flag = "_rolling"
  
  df_long = df_long %>% 
    group_by(code) %>% 
    slice(tail(row_number(), window_size)) %>%
    ungroup()
  
  df_wide = tail(df_wide, window_size)
  
} else {
  rolling_flag = ""
}

#### Models ####
output_for_table = function(model, method="ols", significance=6){
  
  get_stars = function(pval) {
    if (pval < 0.01) {
      stars = "***"
    } else if (pval < 0.05) {
      stars = "**"
    } else if (pval < 0.1) {
      stars = "*"
    } else {
      stars = ""
    }
  }
  
  if (method == "ols") {
    tvar = "t"
    tvardash = paste0(tvar, " ")
    
  } else {
    
    if (method == "random") {
      tvar = "z"
    } else {
      tvar = "t"
    }
    
    tvardash = paste0(tvar, "-")
  }
  
  stars = coef(summary(model))[, glue("Pr(>|{tvar}|)")] %>%
    sapply(get_stars)
  estimates = coefficients(model) %>%
    signif(significance) %>%
    paste0(stars)
  names(estimates) = names(coefficients(model))
  
  tvals = coef(summary(model))[, glue("{tvardash}value")] %>%
    signif(significance) %>%
    sapply(function(x){paste0("(", x, ")")})
  names(tvals) = names(coefficients(model))
  
  return(list("estimates"=estimates, "tvals"=tvals))
}

transmissions = c("frequency", "density")

#### Regional models ####
regions = unique(df_long$code)

taus = c("One", "Tau")
vars = c("Beta", "BetaTwostep", "Gamma")
suffix = c("", "_Z")
prefix_table = rep(c("Left", "Right"), each=length(vars))
vars_table = rep(vars, times=length(suffix))

table_headers = paste0(prefix_table, vars_table)
results_table = tibble(variables = table_headers)

## National model
column = character(0)
column_z = character(0)

for (transmission in transmissions) {
  tau_var = "One"
  
  X_var = ifelse(transmission == "frequency", "S", "X")
  Y_var = ifelse(transmission == "frequency", "I", "Y")
  Z_var = ifelse(transmission == "frequency", "R", "Z")
  
  fm_X = glue("d{X_var} ~ -1 + {X_var}lag{tau_var}:{Y_var}lag{tau_var}") %>%
    as.formula
  fm_Y_twostep = glue("d{Y_var}_twostep ~ -1 + ",
                      "{X_var}lag{tau_var}:{Y_var}lag{tau_var}") %>%
    as.formula
  fm_Z = glue("d{Z_var} ~ -1 + {Y_var}lagOne") %>%
    as.formula
  
  model_s = lm(fm_X, data=df_wide)
  model_r = lm(fm_Z, data=df_wide)
  
  gamma = unname(coef(model_r))
  
  df_wide[[glue("d{Y_var}_twostep")]] = 
    df_wide[[glue("d{Y_var}")]] + gamma*df_wide[[glue("{Y_var}lagOne")]]
  
  model_i = lm(fm_Y_twostep, data=df_wide)
  
  output_model_s = output_for_table(model_s)
  output_model_i = output_for_table(model_i)
  output_model_r = output_for_table(model_r)
  
  beta_table = unname(output_model_s$estimates)
  beta_z_table = unname(output_model_s$tvals)
  
  if (substring(beta_table, 1, 1) == "-") {
    beta_table = substring(beta_table, 2)
    beta_z_table = paste0("(", substring(beta_z_table, 3))
    
  } else {
    beta_table = paste0("-", beta_table)
    beta_z_table = paste0("(-", substring(beta_z_table, 3))
  }
  
  betas_twostep = unname(output_model_i$estimates)
  betas_twostep_z = unname(output_model_i$tvals)
  gammas = unname(output_model_r$estimates)
  gammas_z = unname(output_model_r$tvals)
  
  column = c(column, beta_table, betas_twostep, gammas)
  column_z = c(column_z, beta_z_table, betas_twostep_z, gammas_z)
}

# Insert parameter estimates in the results table
results_table = results_table %>%
  left_join(tibble("variables" = table_headers,
                   "National" = column,
                   "National_tvals" = column_z),
            by="variables")

for (region in regions) {
  column = character(0)
  column_z = character(0)
  
  data = df_long %>%
    filter(code == !!region)
  
  for (transmission in transmissions) {
    tau_var = "One"
    
    X_var = ifelse(transmission == "frequency", "S", "X")
    Y_var = ifelse(transmission == "frequency", "I", "Y")
    Z_var = ifelse(transmission == "frequency", "R", "Z")
    
    fm_X = glue("d{X_var} ~ -1 + {X_var}lag{tau_var}:{Y_var}lag{tau_var}") %>%
      as.formula
    fm_Y_twostep = glue("d{Y_var}_twostep ~ -1 + ",
                        "{X_var}lag{tau_var}:{Y_var}lag{tau_var}") %>%
      as.formula
    fm_Z = glue("d{Z_var} ~ -1 + {Y_var}lagOne") %>%
      as.formula
    
    model_s = lm(fm_X, data=data)
    model_r = lm(fm_Z, data=data)
    
    gamma = unname(coef(model_r))
    
    data[[glue("d{Y_var}_twostep")]] = 
      data[[glue("d{Y_var}")]] + gamma*data[[glue("{Y_var}lagOne")]]
    
    model_i = lm(fm_Y_twostep, data=data)
    
    output_model_s = output_for_table(model_s)
    output_model_i = output_for_table(model_i)
    output_model_r = output_for_table(model_r)
    
    beta_table = unname(output_model_s$estimates)
    beta_z_table = unname(output_model_s$tvals)
    
    if (substring(beta_table, 1, 1) == "-") {
      beta_table = substring(beta_table, 2)
      beta_z_table = paste0("(", substring(beta_z_table, 3))
      
    } else {
      beta_table = paste0("-", beta_table)
      beta_z_table = paste0("(-", substring(beta_z_table, 3))
    }
    
    betas_twostep = unname(output_model_i$estimates)
    betas_twostep_z = unname(output_model_i$tvals)
    gammas = unname(output_model_r$estimates)
    gammas_z = unname(output_model_r$tvals)
    
    column = c(column, beta_table, betas_twostep, gammas)
    column_z = c(column_z, beta_z_table, betas_twostep_z, gammas_z)
  }
  
  # Insert parameter estimates in the results table
  results_table = results_table %>%
    left_join(tibble("variables" = table_headers,
                     !!glue("{region}") := column,
                     !!glue("{region}_tvals") := column_z),
              by="variables")
  
}

results_table = results_table %>%
  gather(region, val, 2:ncol(results_table)) %>%
  spread(names(results_table)[1], val) %>%
  select(region, all_of(table_headers))

results_table = rbind(
  results_table %>% filter(str_detect(region, "National")),
  results_table %>% filter(!str_detect(region, "National"))) %>%
  column_to_rownames("region")

table_transmission = xtable(results_table)
print(table_transmission,
      file=glue("{output_path}/discrete_SIR_transmission{undoc_flag}",
                "{rolling_flag}.txt"))

### Select frequency-dependent transmission and get results for different taus
results_table = tibble(variables = table_headers)

transmission = "frequency"
X_var = ifelse(transmission == "frequency", "S", "X")
Y_var = ifelse(transmission == "frequency", "I", "Y")
Z_var = ifelse(transmission == "frequency", "R", "Z")

## National model
column = character(0)
column_z = character(0)

for (tau_var in taus) {
  fm_X = glue("d{X_var} ~ -1 + {X_var}lag{tau_var}:{Y_var}lag{tau_var}") %>%
    as.formula
  fm_Y_twostep = glue("d{Y_var}_twostep ~ -1 + ",
                      "{X_var}lag{tau_var}:{Y_var}lag{tau_var}") %>%
    as.formula
  fm_Z = glue("d{Z_var} ~ -1 + {Y_var}lagOne") %>%
    as.formula
  
  model_s = lm(fm_X, data=df_wide)
  model_r = lm(fm_Z, data=df_wide)
  
  gamma = unname(coef(model_r))
  
  df_wide[[glue("d{Y_var}_twostep")]] = 
    df_wide[[glue("d{Y_var}")]] + gamma*df_wide[[glue("{Y_var}lagOne")]]
  
  model_i = lm(fm_Y_twostep, data=df_wide)
  
  output_model_s = output_for_table(model_s)
  output_model_i = output_for_table(model_i)
  output_model_r = output_for_table(model_r)
  
  beta_table = unname(output_model_s$estimates)
  beta_z_table = unname(output_model_s$tvals)
  
  if (substring(beta_table, 1, 1) == "-") {
    beta_table = substring(beta_table, 2)
    beta_z_table = paste0("(", substring(beta_z_table, 3))
    
  } else {
    beta_table = paste0("-", beta_table)
    beta_z_table = paste0("(-", substring(beta_z_table, 3))
  }
  
  betas_twostep = unname(output_model_i$estimates)
  betas_twostep_z = unname(output_model_i$tvals)
  gammas = unname(output_model_r$estimates)
  gammas_z = unname(output_model_r$tvals)
  
  column = c(column, beta_table, betas_twostep, gammas)
  column_z = c(column_z, beta_z_table, betas_twostep_z, gammas_z)
}

# Insert parameter estimates in the results table
results_table = results_table %>%
  left_join(tibble("variables" = table_headers,
                   "National" = column,
                   "National_tvals" = column_z),
            by="variables")

for (region in regions) {
  column = character(0)
  column_z = character(0)
  
  data = df_long %>%
    filter(code == !!region)
  
  for (tau_var in taus) {
    
    fm_X = glue("d{X_var} ~ -1 + {X_var}lag{tau_var}:{Y_var}lag{tau_var}") %>%
      as.formula
    fm_Y_twostep = glue("d{Y_var}_twostep ~ -1 + ",
                        "{X_var}lag{tau_var}:{Y_var}lag{tau_var}") %>%
      as.formula
    fm_Z = glue("d{Z_var} ~ -1 + {Y_var}lagOne") %>%
      as.formula
    
    model_s = lm(fm_X, data=data)
    model_r = lm(fm_Z, data=data)
    
    gamma = unname(coef(model_r))
    
    data[[glue("d{Y_var}_twostep")]] = 
      data[[glue("d{Y_var}")]] + gamma*data[[glue("{Y_var}lagOne")]]
    
    model_i = lm(fm_Y_twostep, data=data)
    
    output_model_s = output_for_table(model_s)
    output_model_i = output_for_table(model_i)
    output_model_r = output_for_table(model_r)
    
    beta_table = unname(output_model_s$estimates)
    beta_z_table = unname(output_model_s$tvals)
    
    if (substring(beta_table, 1, 1) == "-") {
      beta_table = substring(beta_table, 2)
      beta_z_table = paste0("(", substring(beta_z_table, 3))
      
    } else {
      beta_table = paste0("-", beta_table)
      beta_z_table = paste0("(-", substring(beta_z_table, 3))
    }
    
    betas_twostep = unname(output_model_i$estimates)
    betas_twostep_z = unname(output_model_i$tvals)
    gammas = unname(output_model_r$estimates)
    gammas_z = unname(output_model_r$tvals)
    
    column = c(column, beta_table, betas_twostep, gammas)
    column_z = c(column_z, beta_z_table, betas_twostep_z, gammas_z)
  }
  
  # Insert parameter estimates in the results table
  results_table = results_table %>%
    left_join(tibble("variables" = table_headers,
                     !!glue("{region}") := column,
                     !!glue("{region}_tvals") := column_z),
              by="variables")
  
}

results_table = results_table %>%
  gather(region, val, 2:ncol(results_table)) %>%
  spread(names(results_table)[1], val) %>%
  select(region, all_of(table_headers))

results_table = rbind(
  results_table %>% filter(str_detect(region, "National")),
  results_table %>% filter(!str_detect(region, "National"))) %>%
  column_to_rownames("region")

table_tau = xtable(results_table)
print(table_tau,
      file=glue("{output_path}/discrete_SIR_taus{undoc_flag}",
                "{rolling_flag}.txt"))

#### Panel Data Models ####
for (transmission in transmissions) {
  results_table = tibble(model = c("POLS", "RE"))
  
  print(glue("\n\n### Transmission: {transmission} ###"))
  
  X_var = ifelse(transmission == "frequency", "S", "X")
  Y_var = ifelse(transmission == "frequency", "I", "Y")
  Z_var = ifelse(transmission == "frequency", "R", "Z")
  
  for (tau_var in taus) {
    print(glue("\n\n#### Lag: {tau_var} ####"))
  
    fm = as.formula("dY ~ -1 + YlagOne:SlagOne")
    
    df_long = 
    
    fm_X = glue("d{X_var} ~ -1 + {X_var}lag{tau_var}:{Y_var}lag{tau_var}") %>%
      as.formula
    fm_Y_twostep = glue("d{Y_var}_twostep ~ -1 + ",
                        "{X_var}lag{tau_var}:{Y_var}lag{tau_var}") %>%
      as.formula
    fm_Z = glue("d{Z_var} ~ -1 + {Y_var}lagOne") %>%
      as.formula
    
    betas = character(0)
    betas_z = character(0)
    betas_twostep = character(0)
    betas_twostep_z = character(0)
    gammas = character(0)
    gammas_z = character(0)
    
    for (method in c("pooling", "random")) {
      print(glue("\n\n## {toupper(method)} ##"))
      
      #### Estimate beta and gamma directly ####
      model_s = plm(fm_X, data=df_long, model = method,
                    index = c("code", "date"))
      model_r = plm(fm_Z, data=df_long, model = method,
                    index = c("code", "date"))
      
      gamma = unname(coef(model_r))
      # beta = -unname(coef(model_s))
      # r = beta / gamma
      
      # Estimate beta with gamma as estimated above
      df_long[[glue("d{Y_var}_twostep")]] = 
        df_long[[glue("d{Y_var}")]] + gamma*df_long[[glue("{Y_var}lagOne")]]
      
      model_i = plm(fm_Y_twostep, data=df_long, model = method,
                    index = c("code", "date"))
      
      # beta_twostep = unname(coef(model_i))
      # r_twostep = beta_twostep / gamma
      
      output_model_s = output_for_table(model_s, method)
      output_model_i = output_for_table(model_i, method)
      output_model_r = output_for_table(model_r, method)
      
      beta_table = unname(output_model_s$estimates)
      beta_z_table = unname(output_model_s$tvals)
      
      if (substring(beta_table, 1, 1) == "-") {
        beta_table = substring(beta_table, 2)
        beta_z_table = paste0("(", substring(beta_z_table, 3))
        
      } else {
        beta_table = paste0("-", beta_table)
        beta_z_table = paste0("(-", substring(beta_z_table, 3))
      }
      
      betas = c(betas, beta_table)
      betas_z = c(betas_z, beta_z_table)
      betas_twostep = c(betas_twostep, unname(output_model_i$estimates))
      betas_twostep_z = c(betas_twostep_z, unname(output_model_i$tvals))
      gammas = c(gammas, unname(output_model_r$estimates))
      gammas_z = c(gammas_z, unname(output_model_r$tvals))
      
    } # end method
    
    results_table = results_table %>%
      left_join(tibble("model" = c("POLS", "RE"),
                       !!glue("{tau_var}Beta") := betas,
                       !!glue("{tau_var}Beta_Z") := betas_z,
                       !!glue("{tau_var}BetaTwostep") := betas_twostep,
                       !!glue("{tau_var}BetaTwostep_Z") := betas_twostep_z,
                       !!glue("{tau_var}Gamma") := gammas,
                       !!glue("{tau_var}Gamma_Z") := gammas_z
                       ),
                by="model")
    
  } # end tau
  
  # Transpose tibble and set row names
  results_table = results_table %>%
    gather(col, val, 2:ncol(results_table)) %>%
    spread(names(results_table)[1], val) %>%
    column_to_rownames("col")
  
  print(xtable(results_table),
        file=glue("{output_path}/panel_table{undoc_flag}{rolling_flag}",
                  "_{transmission}.txt"))
  
} # end transmission

# TODO: Make plots over time

#### NLS - TODO ####
# removed = function(infectives, gamma) {
#   gamma*infectives
# }
# 
# nls(dR ~ -1 + removed(IlagOne, gamma),
#     data=filter(df_long, code == "LOM"), start = list(gamma = 0.001), trace = TRUE)
# 
# # Google: nls r with function
# # https://rpubs.com/RobinLovelace/nls-function
