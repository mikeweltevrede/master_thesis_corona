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
library(latex2exp)
library(plm)
library(snakecase)
library(tidyverse)
library(xtable)

# Import data
df_long = readr::read_csv(path_full_long, col_types = do.call(
  cols, list(date = col_date(format = "%Y-%m-%d"))))
df_wide = readr::read_csv(path_full_wide, col_types = do.call(
  cols, list(date = col_date(format = "%Y-%m-%d"))))

#### Decide the parameters ####
# Maximum incubation period
tau = 14

# Do we want to use a rolling window_size, i.e. only use the most recent
# `window_size` observations?
rolling = TRUE
window_size = 100 + tau

# Determine if we want to model undocumented infectives and, if so, by which
# method. Note that infective_variable is the number of new cases, i.e. Delta X.
form = "" %>%
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
  
  df_long_trimmed = df_long %>% 
    group_by(code) %>% 
    slice(tail(row_number(), window_size)) %>%
    ungroup()
  
  df_wide_trimmed = tail(df_wide, window_size)
  
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
  
  if (rolling) {
    data = df_long_trimmed %>%
      filter(code == !!region)
  } else {
    data = df_long %>%
      filter(code == !!region)
  }
  
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
  
  if (rolling) {
    data = df_long_trimmed %>%
      filter(code == !!region)
  } else {
    data = df_long %>%
      filter(code == !!region)
  }
  
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
      if (rolling) {
        data = df_long_trimmed
      } else {
        data = df_long
      }
      
      model_s = plm(fm_X, data=data, model = method,
                    index = c("code", "date"))
      model_r = plm(fm_Z, data=data, model = method,
                    index = c("code", "date"))
      
      gamma = unname(coef(model_r))
      
      # Estimate beta with gamma as estimated above
      data[[glue("d{Y_var}_twostep")]] = data[[glue("d{Y_var}")]] +
        gamma*data[[glue("{Y_var}lagOne")]]
      
      model_i = plm(fm_Y_twostep, data=data, model = method,
                    index = c("code", "date"))
      
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
      
      print(glue("Beta: {beta_table}"))
      print(glue("Beta two-step: {unname(output_model_i$estimates)}"))
      print(glue("Gamma: {unname(output_model_r$estimates)}"))
      
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

#### Make plots over time ####
tau_var = "Tau"

for (tau_var in taus) {
  fm_X = glue("dS ~ -1 + Slag{tau_var}:Ilag{tau_var}") %>%
    as.formula
  fm_Y_twostep = glue("dI_twostep ~ -1 + Slag{tau_var}:Ilag{tau_var}") %>%
    as.formula
  fm_Z = glue("dR ~ -1 + IlagOne") %>%
    as.formula
  
  betas = numeric(0)
  betas_twostep = numeric(0)
  gammas = numeric(0)
  
  unique_dates = unique(df_long$date)
  
  for (t in window_size:length(unique_dates)){
    start_date = ifelse(rolling, unique_dates[(t-window_size+1)], unique_dates[1])
    
    data = df_long %>%
      filter(date >= start_date & date <= unique_dates[t])
    
    model_s = plm(fm_X, data = data, model = "pooling", index = c("code", "date"))
    model_r = plm(fm_Z, data = data, model = "pooling", index = c("code", "date"))
    
    gamma = unname(coef(model_r))
    
    data[[glue("dI_twostep")]] = data[[glue("dI")]] + gamma*data[[glue("IlagOne")]]
    
    model_i = plm(fm_Y_twostep, data = data, model = "pooling",
                  index = c("code", "date"))
    
    beta = -unname(coef(model_s))
    beta_twostep = unname(coef(model_i))
    
    gammas = c(gammas, gamma)
    betas = c(betas, beta)
    betas_twostep = c(betas_twostep, beta_twostep)
  }
  
  rs = betas/gammas
  rs_twostep = betas_twostep/gammas
  
  tbl_plot = tibble(date = unique_dates[window_size:length(unique_dates)],
                    betas = betas, betas_twostep = betas_twostep,
                    gammas = gammas, rs = rs, rs_twostep = rs_twostep)
  
  firstColor = "#0072B2" # Dark blue
  secondColor = "#D55E00" # Orange-brown
  
  # Plot of Rs
  tbl_plot %>%
    ggplot(aes(x = date)) + 
    geom_point(aes(y = rs, color = firstColor)) +
    geom_smooth(aes(y = rs, color = firstColor), method="loess",
                span=0.3, se=FALSE) +
    geom_point(aes(y = rs_twostep, color = secondColor)) +
    geom_smooth(aes(y = rs_twostep, color = secondColor), method="loess",
                span=0.3, se=FALSE) +
    xlab("") +
    ylab(TeX("$R_{eff}$")) +
    scale_colour_manual(name = "", 
                        labels = c("Regular", "Two-Step"),
                        values = c(firstColor, secondColor)) +
    theme(
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    )
  
  ggsave(glue("panel_data_lag{tau}_Reff{undoc_flag}{rolling_flag}.pdf"),
         path=output_path, width = 8.08, height = 6.24, units = "in")
  
  # Plot of betas
  tbl_plot %>%
    ggplot(aes(x = date)) +
    geom_point(aes(y = betas, color = firstColor)) +
    geom_smooth(aes(y = betas, color = firstColor), method="loess",
                span=0.3, se=FALSE) +
    geom_point(aes(y = betas_twostep, color = secondColor)) +
    geom_smooth(aes(y = betas_twostep, color = secondColor), method="loess",
                span=0.3, se=FALSE) +
    xlab("") +
    ylab(TeX("$\\beta$")) +
    scale_colour_manual(name = "", 
                        labels = c("Regular", "Two-Step"),
                        values = c(firstColor, secondColor)) +
    theme(
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    )
  
  ggsave(glue("panel_data_lag{tau}_betas{undoc_flag}{rolling_flag}.pdf"),
         path=output_path, width = 8.08, height = 6.24, units = "in")
  
  # Plot of gamma
  tbl_plot %>%
    ggplot(aes(x = date)) +
    geom_point(aes(y = gammas, color = firstColor)) +
    geom_smooth(aes(y = gammas, color = firstColor), method="loess",
                span=0.3, se=FALSE) +
    xlab("") +
    ylab(TeX("$\\gamma$")) +
    scale_colour_manual(values = firstColor) +
    theme(legend.position = "none")
  
  ggsave(glue("panel_data_lag{tau}_gammas{undoc_flag}{rolling_flag}.pdf"),
         path=output_path, width = 8.08, height = 6.24, units = "in")
}
