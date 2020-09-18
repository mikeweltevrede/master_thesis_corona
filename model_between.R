rm(list=ls())

#### Model 3 - Within and between-region spread ####
# \Delta i_rt = beta_within*S_rt-tau*\Delta i_rt-tau +
#               beta_between*S_rt-tau*\sum_{c!=r}\Delta i_ct-tau +
#               delta*X_rt + nu_rt
# i_rt is the absolute number of new cases!

#### Set-up ####
# Import standard variables and activate Python environment
source("config.R")

# Import packages
library(glue)
library(gridExtra)
library(latex2exp)
library(snakecase)
library(tidyverse)
library(xtable)

# Import data
df_meta = readxl::read_xlsx(path_metadata, sheet = "Metadata")

df_long = readr::read_csv(path_full_long, col_types = do.call(
  cols, list(date = col_date(format = "%Y-%m-%d")))) %>%
  left_join(df_meta %>% select(code, regionGH, direction), by="code")
df_wide = readr::read_csv(path_full_wide, col_types = do.call(
  cols, list(date = col_date(format = "%Y-%m-%d"))))

# Get all region abbreviations
regions = unique(df_long$code)

# Select regressors
M_regressors = c("weekend")

#### Decide the parameters ####
# You need to adapt, if desired, the following parameters:
# tau: int, maximum incubation period
# rolling: boolean, whether to apply a rolling window
# window_size: int, if using a rolling window, how large?
# form: str, the form of undocumented infections to model with (if any)

# Maximum incubation period
tau = 14

# Do we want to use a rolling window_size, i.e. only use the most recent
# `window_size` observations?
rolling = TRUE
window_size = 100 + tau

if (rolling) {
  rolling_flag = "_rolling"
} else {
  rolling_flag = ""
}

# Determine if we want to model undocumented infectives and, if so, by which
# method. Note that infective_variable is the number of new cases, i.e. Delta i.
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

# Retrieve parameter estimates
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

#### Data preparation ####
df_sumInc = tibble(date = as.Date(NA),
                   code = character(),
                   sumInfectives = numeric())

for (region in regions){
  df_sumInc = df_sumInc %>%
    bind_rows(
      df_wide %>%
        select(
          map(regions[regions != region], starts_with, vars = colnames(.)) %>%
            unlist()) %>%
        select(ends_with(infective_variable)) %>%
        mutate_all(dplyr::lag, n=tau) %>%
        transmute(
          date = df_wide$date,
          code = region,
          sumInfectives = rowSums(.)))
}

df_long = df_long %>%
  left_join(df_sumInc, by = c("code", "date"))

#### Models without model selection ####
fm = glue("{infective_variable} ~ -1 +",
          "dplyr::lag({infective_variable}, {tau}):",
          "dplyr::lag(susceptibleRate, {tau})+",
          "dplyr::lag(susceptibleRate, {tau}):sumInfectives +",
          paste(M_regressors, collapse="+")) %>%
  as.formula

# Create results table
results_table = tibble(variables = c(M_regressors, "beta_w", "beta_b"))

for (region in regions){
  # Select only the data for the relevant region
  data = df_long %>%
    filter(code == !!region)
  
  # Estimate the model by OLS
  if (rolling) {
    model = lm(fm, data=tail(data, window_size))
  } else {
    model = lm(fm, data=data)
  }
  
  # Retrieve parameter estimates
  table_output = output_for_table(model)
  estimates = table_output$estimates
  tvals = table_output$tvals
  
  # Insert parameter estimates in the results table
  results_table = results_table %>%
    left_join(tibble("variables" = c(M_regressors, "beta_w", "beta_b"),
                     !!glue("{region}") := unname(
                       c(estimates[M_regressors],
                         estimates[
                           glue("dplyr::lag({infective_variable}, {tau}):",
                                "dplyr::lag(susceptibleRate, {tau})")],
                         estimates[glue("dplyr::lag(susceptibleRate, {tau}):",
                                        "sumInfectives")])),
                     !!glue("{region}_tvals") := unname(
                       c(tvals[M_regressors],
                         tvals[glue("dplyr::lag({infective_variable}, {tau}):",
                                    "dplyr::lag(susceptibleRate, {tau})")],
                         tvals[glue("dplyr::lag(susceptibleRate, {tau}):",
                                    "sumInfectives")]))),
              by="variables")
}

# Transpose and reorder columns
results_table = results_table %>%
  gather(region, val, 2:ncol(results_table)) %>%
  spread(names(results_table)[1], val) %>%
  select(region, beta_w, beta_b, everything()) %>%
  column_to_rownames("region")

# Return a LaTeX table
table_no_ms = xtable(results_table, math.style.exponents = TRUE)
print(table_no_ms,
      file=glue("{output_path}/model_between_table_no_ms{undoc_flag}",
                "{rolling_flag}.txt"))

#### Models with model selection (AIC) ####
results_table_aic = tibble(variables = c(M_regressors, "beta_w", "beta_b"))

# Construct formula
fm = glue("{infective_variable} ~ -1 +",
          "dplyr::lag({infective_variable}, {tau}):",
          "dplyr::lag(susceptibleRate, {tau})+",
          "dplyr::lag(susceptibleRate, {tau}):sumInfectives +",
          paste(M_regressors, collapse="+")) %>%
  as.formula

for (region in regions){
  # Select only the data for the relevant region
  data = df_long %>%
    filter(code == !!region)
  
  # Use AIC for model selection
  if (rolling) {
    model = step(lm(fm, data=tail(data, window_size)), k=2, trace=0,
                 scope=list(
                   "lower" = glue(
                     "{infective_variable} ~ -1 +",
                     "dplyr::lag({infective_variable}, {tau}):",
                     "dplyr::lag(susceptibleRate, {tau})+",
                     "dplyr::lag(susceptibleRate, {tau}):sumInfectives") %>%
                     as.formula,
                   "upper" = fm))
  } else {
    model = step(lm(fm, data=data), k=2, trace=0,
                 scope=list(
                   "lower" = glue(
                     "{infective_variable} ~ -1 +",
                     "dplyr::lag({infective_variable}, {tau}):",
                     "dplyr::lag(susceptibleRate, {tau})+",
                     "dplyr::lag(susceptibleRate, {tau}):sumInfectives") %>%
                     as.formula,
                   "upper" = fm))
  }
  
  # Retrieve parameter estimates
  table_output = output_for_table(model)
  estimates = table_output$estimates
  tvals = table_output$tvals
  
  # Insert parameter estimates in the results table
  results_table_aic = results_table_aic %>%
    left_join(tibble("variables" = c(M_regressors, "beta_w", "beta_b"),
                     !!glue("{region}") := unname(
                       c(estimates[M_regressors],
                         estimates[
                           glue("dplyr::lag({infective_variable}, {tau}):",
                                "dplyr::lag(susceptibleRate, {tau})")],
                         estimates[glue("dplyr::lag(susceptibleRate, {tau}):",
                                        "sumInfectives")])),
                     !!glue("{region}_tvals") := unname(
                       c(tvals[M_regressors],
                         tvals[glue("dplyr::lag({infective_variable}, {tau}):",
                                    "dplyr::lag(susceptibleRate, {tau})")],
                         tvals[glue("dplyr::lag(susceptibleRate, {tau}):",
                                    "sumInfectives")]))),
              by="variables")
}

# Transpose and reorder columns
results_table_aic = results_table_aic %>%
  gather(region, val, 2:ncol(results_table_aic)) %>%
  spread(names(results_table_aic)[1], val) %>%
  select(region, beta_w, beta_b, everything()) %>%
  column_to_rownames("region")

# Return a LaTeX table
table_aic = xtable(results_table_aic, math.style.exponents = TRUE)
print(table_aic,
      file=glue("{output_path}/model_between_table_aic{undoc_flag}",
                "{rolling_flag}.txt"))

#### Models with model selection (BIC) ####
results_table_bic = tibble(variables = c(M_regressors, "beta_w", "beta_b"))

# Construct formula
fm = glue("{infective_variable} ~ -1 +",
          "dplyr::lag({infective_variable}, {tau}):",
          "dplyr::lag(susceptibleRate, {tau})+",
          "dplyr::lag(susceptibleRate, {tau}):sumInfectives +",
          paste(M_regressors, collapse="+")) %>%
  as.formula

for (region in regions){
  # Select only the data for the relevant region
  data = df_long %>%
    filter(code == !!region)
  
  # Use BIC for model selection
  if (rolling) {
    model = step(lm(fm, data=tail(data, window_size)), k=log(window_size),
                 trace=0, scope=list(
                   "lower" = glue(
                     "{infective_variable} ~ -1 +",
                     "dplyr::lag({infective_variable}, {tau}):",
                     "dplyr::lag(susceptibleRate, {tau})+",
                     "dplyr::lag(susceptibleRate, {tau}):sumInfectives") %>%
                     as.formula,
                   "upper" = fm))
  } else {
    model = step(lm(fm, data=data), k=log(nrow(data)), trace=0,
                 scope=list(
                   "lower" = glue(
                     "{infective_variable} ~ -1 +",
                     "dplyr::lag({infective_variable}, {tau}):",
                     "dplyr::lag(susceptibleRate, {tau})+",
                     "dplyr::lag(susceptibleRate, {tau}):sumInfectives") %>%
                     as.formula,
                   "upper" = fm))
  }
  
  # Retrieve parameter estimates
  table_output = output_for_table(model)
  estimates = table_output$estimates
  tvals = table_output$tvals
  
  # Insert parameter estimates in the results table
  results_table_bic = results_table_bic %>%
    left_join(tibble("variables" = c(M_regressors, "beta_w", "beta_b"),
                     !!glue("{region}") := unname(
                       c(estimates[M_regressors],
                         estimates[
                           glue("dplyr::lag({infective_variable}, {tau}):",
                                "dplyr::lag(susceptibleRate, {tau})")],
                         estimates[glue("dplyr::lag(susceptibleRate, {tau}):",
                                        "sumInfectives")])),
                     !!glue("{region}_tvals") := unname(
                       c(tvals[M_regressors],
                         tvals[glue("dplyr::lag({infective_variable}, {tau}):",
                                    "dplyr::lag(susceptibleRate, {tau})")],
                         tvals[glue("dplyr::lag(susceptibleRate, {tau}):",
                                    "sumInfectives")]))),
              by="variables")
}

# Transpose and reorder columns
results_table_bic = results_table_bic %>%
  gather(region, val, 2:ncol(results_table_bic)) %>%
  spread(names(results_table_bic)[1], val) %>%
  select(region, beta_w, beta_b, everything()) %>%
  column_to_rownames("region")

# Return a LaTeX table
table_bic = xtable(results_table_bic, math.style.exponents = TRUE)
print(table_bic,
      file=glue("{output_path}/model_between_table_bic{undoc_flag}",
                "{rolling_flag}.txt"))

#### Plot beta over time ####
fm = glue("{infective_variable} ~ -1 +",
          "dplyr::lag({infective_variable}, {tau}):",
          "dplyr::lag(susceptibleRate, {tau})+",
          "dplyr::lag(susceptibleRate, {tau}):sumInfectives +",
          paste(M_regressors, collapse="+")) %>%
  as.formula

withinColor = "#0072B2" # Dark blue
betweenColor = "#D55E00" # Orange-brown

#### Without model selection ####
tbl_beta = tibble(date = as.Date(NA), Within=numeric(0), Between=numeric(0),
                  code=character(0))

# Find the estimates of beta per region over time
for (region in regions){
  betas_w = vector("double")
  betas_b = vector("double")
  dates = vector("character")
  
  # Select only the data for the relevant region
  data = df_long %>%
    filter(code == !!region)
  
  for (t in window_size:nrow(data)){
    # Estimate the model by OLS
    if (rolling) {
      model = lm(fm, data=data[(t-window_size+1):t, ])
    } else {
      model = lm(fm, data=head(data, t))
    }
    
    # Retrieve the beta estimate and append this to the list of betas
    beta_w = coef(model)[[glue("dplyr::lag({infective_variable}, {tau}):",
                               "dplyr::lag(susceptibleRate, {tau})")]]
    beta_b = coef(model)[[
      glue("dplyr::lag(susceptibleRate, {tau}):sumInfectives")]]
    betas_w = c(betas_w, beta_w)
    betas_b = c(betas_b, beta_b)
  }
  
  # Append the results to the table
  tbl_beta = tbl_beta %>%
    bind_rows(tibble(date = data$date[window_size:nrow(data)],
                     Within = betas_w,
                     Between = betas_b,
                     code = region) %>%
                drop_na())
}

# Add region and direction to the table
tbl_beta = tbl_beta %>%
  left_join(df_meta %>% select(c(regionGH, code, direction)), by="code")

# Make a plot per direction; facet_wrap by regions
for (sub_tbl in split(tbl_beta, tbl_beta$direction)){
  direc = sub_tbl$direction[1]
  
  plots = vector("list")
  
  sub_tbl = sub_tbl %>% 
    pivot_longer(all_of(c("Within", "Between")))
  
  for (region in unique(sub_tbl$code)) {
    region_full = filter(sub_tbl, code == !!region)$regionGH[1]
    
    plots[[region]] = 
      sub_tbl %>%
      filter(code == !!region) %>% 
      ggplot(aes(x = date)) + 
      geom_point(aes(y = value, color = name)) +
      geom_smooth(aes(y = value, color = name), method="loess", span=0.3,
                  se=FALSE)  +
      facet_wrap(vars(name), ncol=1, scales="free_y") +
      xlab("") +
      ylab("") +
      ggtitle(region_full) +
      scale_colour_manual(name = "Beta", 
                          labels = c("Within", "Between"),
                          values = c(withinColor, betweenColor)) +
      theme(legend.position = "none")
  }
  
  g = do.call("grid.arrange", c(plots, ncol=floor(sqrt(length(plots)))))
  
  ggsave(
    glue("model_between_lag{tau}_betas_{direc}{undoc_flag}{rolling_flag}.pdf"),
    plot = g, path = output_path)
}

#### With model selection (AIC) ####
tbl_beta = tibble(date = as.Date(NA), Within=numeric(0), Between=numeric(0),
                  code=character(0))

# Find the estimates of beta per region over time
for (region in regions){
  betas_w = vector("double")
  betas_b = vector("double")
  dates = vector("character")
  
  # Select only the data for the relevant region
  data = df_long %>%
    filter(code == !!region)
  
  for (t in window_size:nrow(data)){
    # Use AIC for model selection - scope says we want to always keep
    # beta_within in
    if (rolling) {
      model = step(lm(fm, data=data[(t-window_size+1):t, ]), k=2, trace=0,
                   scope=list(
                     "lower" = glue(
                       "{infective_variable} ~ -1 +",
                       "dplyr::lag({infective_variable}, {tau}):",
                       "dplyr::lag(susceptibleRate, {tau})+",
                       "dplyr::lag(susceptibleRate, {tau}):sumInfectives") %>%
                       as.formula,
                     "upper" = fm))
    } else {
      model = step(lm(fm, data=head(data, t)), k=2, trace=0,
                   scope=list(
                     "lower" = glue(
                       "{infective_variable} ~ -1 +",
                       "dplyr::lag({infective_variable}, {tau}):",
                       "dplyr::lag(susceptibleRate, {tau})+",
                       "dplyr::lag(susceptibleRate, {tau}):sumInfectives") %>%
                       as.formula,
                     "upper" = fm))
    }
    
    # Retrieve the beta estimate and append this to the list of betas
    beta_w = coef(model)[[glue("dplyr::lag({infective_variable}, {tau}):",
                               "dplyr::lag(susceptibleRate, {tau})")]]
    beta_b = coef(model)[[
      glue("dplyr::lag(susceptibleRate, {tau}):sumInfectives")]]
    betas_w = c(betas_w, beta_w)
    betas_b = c(betas_b, beta_b)
  }
  
  # Append the results to the table
  tbl_beta = tbl_beta %>%
    bind_rows(tibble(date = data$date[window_size:nrow(data)],
                     Within = betas_w,
                     Between = betas_b,
                     code = region) %>%
                drop_na())
}

# Add region and direction to the table
tbl_beta = tbl_beta %>%
  left_join(df_meta %>% select(c(regionGH, code, direction)), by="code")

# Make a plot per direction; facet_wrap by regions
for (sub_tbl in split(tbl_beta, tbl_beta$direction)){
  direc = sub_tbl$direction[1]
  coeff = median(sub_tbl$Within / sub_tbl$Between)
  
  plots = vector("list")
  
  sub_tbl = sub_tbl %>% 
    pivot_longer(all_of(c("Within", "Between")))
  
  for (region in unique(sub_tbl$code)) {
    region_full = filter(sub_tbl, code == !!region)$regionGH[1]
    
    plots[[region]] = 
      sub_tbl %>%
      filter(code == !!region) %>% 
      ggplot(aes(x = date)) + 
      geom_point(aes(y = value, color = name)) +
      geom_smooth(aes(y = value, color = name), method="loess", span=0.3,
                  se=FALSE)  +
      facet_wrap(vars(name), ncol=1, scales="free_y") +
      xlab("") +
      ylab("") +
      ggtitle(region_full) +
      scale_colour_manual(name = "Beta", 
                          labels = c("Within", "Between"),
                          values = c(withinColor, betweenColor)) +
      theme(legend.position = "none")
  }
  
  g = do.call("grid.arrange", c(plots, ncol=floor(sqrt(length(plots)))))
  
  ggsave(
    glue("model_between_lag{tau}_betas_{direc}_aic{undoc_flag}{rolling_flag}.pdf"),
    plot = g, path = output_path)
}

#### With model selection (BIC) ####
tbl_beta = tibble(date = as.Date(NA), Within=numeric(0), Between=numeric(0),
                  code=character(0))

# Find the estimates of beta per region over time
for (region in regions){
  betas_w = vector("double")
  betas_b = vector("double")
  dates = vector("character")
  
  # Select only the data for the relevant region
  data = df_long %>%
    filter(code == !!region)
  
  for (t in window_size:nrow(data)){
    # Use BIC for model selection - scope says we want to always keep
    # beta_within in
    if (rolling) {
      model = step(lm(fm, data=data[(t-window_size+1):t, ]), k=log(window_size),
                   trace=0, scope=list(
                     "lower" = glue(
                       "{infective_variable} ~ -1 +",
                       "dplyr::lag({infective_variable}, {tau}):",
                       "dplyr::lag(susceptibleRate, {tau})+",
                       "dplyr::lag(susceptibleRate, {tau}):sumInfectives") %>%
                       as.formula,
                     "upper" = fm))
    } else {
      model = step(lm(fm, data=head(data, t)), k=log(t), trace=0,
                   scope=list(
                     "lower" = glue(
                       "{infective_variable} ~ -1 +",
                       "dplyr::lag({infective_variable}, {tau}):",
                       "dplyr::lag(susceptibleRate, {tau})+",
                       "dplyr::lag(susceptibleRate, {tau}):sumInfectives") %>%
                       as.formula,
                     "upper" = fm))
    }
    
    # Retrieve the beta estimate and append this to the list of betas
    beta_w = coef(model)[[glue("dplyr::lag({infective_variable}, {tau}):",
                               "dplyr::lag(susceptibleRate, {tau})")]]
    beta_b = coef(model)[[
      glue("dplyr::lag(susceptibleRate, {tau}):sumInfectives")]]
    betas_w = c(betas_w, beta_w)
    betas_b = c(betas_b, beta_b)
  }
  
  # Append the results to the table
  tbl_beta = tbl_beta %>%
    bind_rows(tibble(date = data$date[window_size:nrow(data)],
                     Within = betas_w,
                     Between = betas_b,
                     code = region) %>%
                drop_na())
}

# Add region and direction to the table
tbl_beta = tbl_beta %>%
  left_join(df_meta %>% select(c(regionGH, code, direction)), by="code")

# Make a plot per direction; facet_wrap by regions
for (sub_tbl in split(tbl_beta, tbl_beta$direction)){
  direc = sub_tbl$direction[1]
  coeff = median(sub_tbl$Within / sub_tbl$Between)
  
  plots = vector("list")
  
  sub_tbl = sub_tbl %>% 
    pivot_longer(all_of(c("Within", "Between")))
  
  for (region in unique(sub_tbl$code)) {
    region_full = filter(sub_tbl, code == !!region)$regionGH[1]
    
    plots[[region]] = 
      sub_tbl %>%
      filter(code == !!region) %>% 
      ggplot(aes(x = date)) + 
      geom_point(aes(y = value, color = name)) +
      geom_smooth(aes(y = value, color = name), method="loess", span=0.3,
                  se=FALSE)  +
      facet_wrap(vars(name), ncol=1, scales="free_y") +
      xlab("") +
      ylab("") +
      ggtitle(region_full) +
      scale_colour_manual(name = "Beta", 
                          labels = c("Within", "Between"),
                          values = c(withinColor, betweenColor)) +
      theme(legend.position = "none")
  }
  
  g = do.call("grid.arrange", c(plots, ncol=floor(sqrt(length(plots)))))
  
  ggsave(
    glue("model_between_lag{tau}_betas_{direc}_bic{undoc_flag}{rolling_flag}.pdf"),
    plot = g, path = output_path)
}
