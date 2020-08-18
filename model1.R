rm(list=ls())

#### Model 1 - Within-region spread ####
# Delta Y_rt = beta_within*Delta Y_rt-tau*S_rt-tau + delta*M_rt + nu_rt
# Y_rt is the absolute number of cases!

#### Set-up ####
# Import standard variables and activate Python environment
source("config.R")

# Import packages
library(glue)
library(latex2exp)
library(snakecase)
library(tidyverse)
library(xtable)

# Import data
df_long = readr::read_csv(path_full_long, col_types = do.call(
  cols, list(date = col_date(format = "%Y-%m-%d"))))
df_wide = readr::read_csv(path_full_wide, col_types = do.call(
  cols, list(date = col_date(format = "%Y-%m-%d"))))

# Get all region abbreviations
regions = unique(df_long$code)

# Select regressors
M_regressors = c("weekend")
all_variables = c("(Intercept)", M_regressors)

#### Decide the parameters ####
# You need to adapt, if desired, the following parameters:
# tau: int, the latent period
# rolling: boolean, whether to apply a rolling window
# window_size: int, if using a rolling window, how large?
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

#### Models without model selection ####
# Create results table
results_table = tibble(variables = c(all_variables, "beta"))

#### National model ####
# Construct formula
fm = glue("infectivesNational{form} ~ ",
          "lag(infectivesNational{form}, {tau}):",
          "lag(susceptibleRateNational, {tau}) +",
  paste(M_regressors, collapse="+")) %>%
  as.formula

if (rolling) {
  model = lm(fm, data=tail(df_wide, window_size))
} else {
  model = lm(fm, data=df_wide)
}
summary(model)

#### Make parameter table ####
# Retrieve parameter estimates
output_for_table = function(model, significance=4){
  
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
  
  stars = coef(summary(model))[, "Pr(>|t|)"] %>%
    sapply(get_stars)
  estimates = coef(summary(model))[, "Estimate"] %>%
    signif(significance) %>%
    paste0(stars)
  names(estimates) = names(coef(summary(model))[, "Estimate"])
  
  tvals = coef(summary(model))[, "t value"] %>%
    signif(significance) %>%
    sapply(function(x){paste0("(", x, ")")})
  
  return(list("estimates"=estimates, "tvals"=tvals))
}

table_output = output_for_table(model)
estimates = table_output$estimates
tvals = table_output$tvals

# Insert parameter estimates in the results table
results_table = results_table %>%
  left_join(tibble("variables" = c(all_variables, "beta"),
                   "National" = unname(
                     c(estimates[all_variables],
                       estimates[glue("lag(infectivesNational{form}, {tau}):",
                                      "lag(susceptibleRateNational, {tau})")])),
                   "National_tvals" = unname(
                     c(tvals[all_variables],
                       tvals[glue("lag(infectivesNational{form}, {tau}):",
                                  "lag(susceptibleRateNational, {tau})")]))),
            by="variables")

#### Regional models ####
# Construct formula
fm = glue("{infective_variable} ~ ",
          "lag({infective_variable}, {tau}):lag(susceptibleRate, {tau})+",
          paste(M_regressors, collapse="+")) %>%
  as.formula

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
    left_join(tibble("variables" = c(all_variables, "beta"),
                     !!glue("{region}") := unname(
                       c(estimates[all_variables],
                         estimates[glue("lag({infective_variable}, {tau}):",
                                        "lag(susceptibleRate, {tau})")])),
                     !!glue("{region}_tvals") := unname(
                       c(tvals[all_variables],
                         tvals[glue("lag({infective_variable}, {tau}):",
                                    "lag(susceptibleRate, {tau})")]))),
              by="variables")
}

# Transpose and reorder columns
results_table = results_table %>%
  gather(region, val, 2:ncol(results_table)) %>%
  spread(names(results_table)[1], val) %>%
  select(region, beta, everything())

# Put the national results at the top
results_table = rbind(
  results_table %>% filter(str_detect(region, "National")),
  results_table %>% filter(!str_detect(region, "National"))) %>%
  column_to_rownames("region")

# Return a LaTeX table
table_no_ms = xtable(results_table, math.style.exponents = TRUE)
print(table_no_ms,
      file=glue("{output_path}/model1_table_no_ms{undoc_flag}{rolling_flag}.txt"))

#### Models with model selection (AIC) ####
results_table_aic = tibble(variables = c(all_variables, "beta"))

#### National model ####
# Construct formula
fm = glue("infectivesNational{form} ~ ",
          "lag(infectivesNational{form}, {tau}):",
          "lag(susceptibleRateNational, {tau})+",
          paste(M_regressors, collapse="+")) %>%
  as.formula

# Use AIC for model selection - scope says we want to always keep beta_within in
if (rolling) {
  model = step(lm(fm, data=tail(df_wide, window_size)), k=2, trace=0,
               scope=list(
                 "lower" = glue("infectivesNational{form} ~ ",
                                "lag(infectivesNational{form}, {tau}):",
                                "lag(susceptibleRateNational, {tau})") %>%
                   as.formula,
                 "upper" = fm))
} else {
  model = step(lm(fm, data=df_wide), k=2, trace=0,
               scope=list(
                 "lower" = glue("infectivesNational{form} ~ ",
                                "lag(infectivesNational{form}, {tau}):",
                                "lag(susceptibleRateNational, {tau})") %>%
                   as.formula,
                 "upper" = fm))
}

# Retrieve parameter estimates
table_output = output_for_table(model)
estimates = table_output$estimates
tvals = table_output$tvals

# Insert parameter estimates in the results table
results_table_aic = results_table_aic %>%
  left_join(tibble("variables" = c(all_variables, "beta"),
                   "National" = unname(
                     c(estimates[all_variables],
                       estimates[glue("lag(infectivesNational{form}, {tau}):",
                                      "lag(susceptibleRateNational, {tau})")])),
                   "National_tvals" = unname(
                     c(tvals[all_variables],
                       tvals[glue("lag(infectivesNational{form}, {tau}):",
                                  "lag(susceptibleRateNational, {tau})")]))),
            by="variables")

#### Regional models ####
# Construct formula
fm = glue("{infective_variable} ~ ",
          "lag({infective_variable}, {tau}):lag(susceptibleRate, {tau})+",
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
                   "lower" = glue("{infective_variable} ~ ",
                                  "lag({infective_variable}, {tau}):",
                                  "lag(susceptibleRate, {tau})") %>%
                     as.formula,
                   "upper" = fm))
  } else {
    model = step(lm(fm, data=data), k=2, trace=0,
                 scope=list(
                   "lower" = glue("{infective_variable} ~ ",
                                  "lag({infective_variable}, {tau}):",
                                  "lag(susceptibleRate, {tau})") %>%
                     as.formula,
                   "upper" = fm))
  }
  
  # Retrieve parameter estimates
  table_output = output_for_table(model)
  estimates = table_output$estimates
  tvals = table_output$tvals
  
  # Insert parameter estimates in the results table
  results_table_aic = results_table_aic %>%
    left_join(tibble("variables" = c(all_variables, "beta"),
                     !!glue("{region}") := unname(
                       c(estimates[all_variables],
                         estimates[glue("lag({infective_variable}, {tau}):",
                                        "lag(susceptibleRate, {tau})")])),
                     !!glue("{region}_tvals") := unname(
                       c(tvals[all_variables],
                         tvals[glue("lag({infective_variable}, {tau}):",
                                    "lag(susceptibleRate, {tau})")]))),
              by="variables")
}

# Transpose and reorder columns
results_table_aic = results_table_aic %>%
  gather(region, val, 2:ncol(results_table_aic)) %>%
  spread(names(results_table_aic)[1], val) %>%
  select(region, beta, everything())

# Put the national results at the top
results_table_aic = rbind(
  results_table_aic %>% filter(str_detect(region, "National")),
  results_table_aic %>% filter(!str_detect(region, "National"))) %>%
  column_to_rownames("region")

# Return a LaTeX table
table_ms_aic = xtable(results_table_aic, math.style.exponents = TRUE)
print(table_ms_aic,
      file=glue("{output_path}/model1_table_aic{undoc_flag}{rolling_flag}.txt"))

#### Models with model selection (BIC) ####
results_table_bic = tibble(variables = c(all_variables, "beta"))

#### National model ####
# Construct formula
fm = glue("infectivesNational{form} ~ ",
          "lag(infectivesNational{form}, {tau}):",
          "lag(susceptibleRateNational, {tau})+",
          paste(M_regressors, collapse="+")) %>%
  as.formula

# Use BIC for model selection - scope says we want to always keep beta_within in
if (rolling) {
  model = step(lm(fm, data=tail(df_wide, window_size)), k=log(window_size), trace=0,
               scope=list(
                 "lower" = glue("infectivesNational{form} ~ ",
                                "lag(infectivesNational{form}, {tau}):",
                                "lag(susceptibleRateNational, {tau})") %>%
                   as.formula,
                 "upper" = fm))
} else {
  model = step(lm(fm, data=df_wide), k=log(nrow(df_wide)), trace=0,
               scope=list(
                 "lower" = glue("infectivesNational{form} ~ ",
                                "lag(infectivesNational{form}, {tau}):",
                                "lag(susceptibleRateNational, {tau})") %>%
                   as.formula,
                 "upper" = fm))
}

# Retrieve parameter estimates
table_output = output_for_table(model)
estimates = table_output$estimates
tvals = table_output$tvals

# Insert parameter estimates in the results table
results_table_bic = results_table_bic %>%
  left_join(tibble("variables" = c(all_variables, "beta"),
                   "National" = unname(
                     c(estimates[all_variables],
                       estimates[glue("lag(infectivesNational{form}, {tau}):",
                                      "lag(susceptibleRateNational, {tau})")])),
                   "National_tvals" = unname(
                     c(tvals[all_variables],
                       tvals[glue("lag(infectivesNational{form}, {tau}):",
                                  "lag(susceptibleRateNational, {tau})")]))),
            by="variables")

#### Regional models ####
# Construct formula
fm = glue("{infective_variable} ~ ",
          "lag({infective_variable}, {tau}):lag(susceptibleRate, {tau})+",
          paste(M_regressors, collapse="+")) %>%
  as.formula

for (region in regions){
  # Select only the data for the relevant region
  data = df_long %>%
    filter(code == !!region)
  
  # Use BIC for model selection
  if (rolling) {
    model = step(lm(fm, data=tail(data, window_size)), k=log(window_size), trace=0,
                 scope=list(
                   "lower" = glue("{infective_variable} ~ ",
                                  "lag({infective_variable}, {tau}):",
                                  "lag(susceptibleRate, {tau})") %>%
                     as.formula,
                   "upper" = fm))
  } else {
    model = step(lm(fm, data=data), k=log(nrow(data)), trace=0,
                 scope=list(
                   "lower" = glue("{infective_variable} ~ ",
                                  "lag({infective_variable}, {tau}):",
                                  "lag(susceptibleRate, {tau})") %>%
                     as.formula,
                   "upper" = fm))
  }
  
  # Retrieve parameter estimates
  table_output = output_for_table(model)
  estimates = table_output$estimates
  tvals = table_output$tvals
  
  # Insert parameter estimates in the results table
  results_table_bic = results_table_bic %>%
    left_join(tibble("variables" = c(all_variables, "beta"),
                     !!glue("{region}") := unname(
                       c(estimates[all_variables],
                         estimates[glue("lag({infective_variable}, {tau}):",
                                        "lag(susceptibleRate, {tau})")])),
                     !!glue("{region}_tvals") := unname(
                       c(tvals[all_variables],
                         tvals[glue("lag({infective_variable}, {tau}):",
                                    "lag(susceptibleRate, {tau})")]))),
              by="variables")
}

# Transpose and reorder columns
results_table_bic = results_table_bic %>%
  gather(region, val, 2:ncol(results_table_bic)) %>%
  spread(names(results_table_bic)[1], val) %>%
  select(region, beta, everything())

# Put the national results at the top
results_table_bic = rbind(
  results_table_bic %>% filter(str_detect(region, "National")),
  results_table_bic %>% filter(!str_detect(region, "National"))) %>%
  column_to_rownames("region")

# Return a LaTeX table
table_ms_bic = xtable(results_table_bic, math.style.exponents = TRUE)
print(table_ms_bic,
      file=glue("{output_path}/model1_table_bic{undoc_flag}{rolling_flag}.txt"))

#### Plot beta over time ####
df_meta = readxl::read_xlsx(path_metadata, sheet = "Metadata")

fm = glue("{infective_variable} ~ ",
          "lag({infective_variable}, {tau}):lag(susceptibleRate, {tau})+",
          paste(M_regressors, collapse="+")) %>%
  as.formula

#### Without model selection ####
tbl_beta = tibble(date = as.Date(NA), betas=numeric(0), code=character(0))

# Find the estimates of beta per region over time
for (region in regions){
  betas = vector("double")
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
    beta = coef(model)[[glue("lag({infective_variable}, {tau}):",
                              "lag(susceptibleRate, {tau})")]]
    betas = c(betas, beta)
  }
  
  # Append the results to the table
  tbl_beta = tbl_beta %>%
    bind_rows(tibble(date = data$date[window_size:nrow(data)],
                     betas = betas,
                     code = region) %>%
                drop_na())
}

# Add region and direction to the table
tbl_beta = tbl_beta %>%
  left_join(df_meta %>% select(c(region, code, direction)), by="code")

# Make a plot per direction
for (sub_tbl in split(tbl_beta, tbl_beta$direction)){
  direc = sub_tbl$direction[1]
  g = ggplot(sub_tbl, aes(date, betas, colour = region)) + 
    geom_point() +
    geom_smooth(method="loess", span=0.3, se=FALSE) +
    xlab("") +
    ylab(TeX("$\\beta_{within}$")) +
    scale_colour_manual(values=c("#0072B2", # Dark blue
                                 "#D55E00", # Orange-brown
                                 "#CC79A7", # Pink
                                 "#009E73", # Green
                                 "#56B4E9", # Light blue
                                 "#E69F00")) # Yellow
  print(g)
  ggsave(
    glue("model1_lag{tau}_betawithin_{direc}{undoc_flag}{rolling_flag}.pdf"),
    path=output_path, width = 10.8, height = 6.62, units = "in")
}

#### With model selection (AIC) ####
tbl_beta = tibble(date = as.Date(NA), betas=numeric(0), code=character(0))

# Find the estimates of beta per region over time
for (region in regions){
  betas = vector("double")
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
                     "lower" = glue("{infective_variable} ~ ",
                                    "lag({infective_variable}, {tau}):",
                                    "lag(susceptibleRate, {tau})") %>%
                       as.formula,
                     "upper" = fm))
    } else {
      model = step(lm(fm, data=head(data, t)), k=2, trace=0,
                   scope=list(
                     "lower" = glue("{infective_variable} ~ ",
                                    "lag({infective_variable}, {tau}):",
                                    "lag(susceptibleRate, {tau})") %>%
                       as.formula,
                     "upper" = fm))
    }
    
    # Retrieve the beta estimate and append this to the list of betas
    beta = coef(model)[[glue("lag({infective_variable}, {tau}):",
                              "lag(susceptibleRate, {tau})")]]
    betas = c(betas, beta)
  }
  
  # Append the results to the table
  tbl_beta = tbl_beta %>%
    bind_rows(tibble(date = data$date[window_size:nrow(data)],
                     betas = betas,
                     code = region) %>%
                drop_na())
}

# Add region and direction to the table
tbl_beta = tbl_beta %>%
  left_join(df_meta %>% select(c(region, code, direction)), by="code")

# Make a plot per direction
for (sub_tbl in split(tbl_beta, tbl_beta$direction)){
  direc = sub_tbl$direction[1]
  g = ggplot(sub_tbl, aes(date, betas, colour = region)) + 
    geom_point() +
    geom_smooth(method="loess", span=0.3, se=FALSE) +
    xlab("") +
    ylab(TeX("$\\beta_{within}$")) +
    scale_colour_manual(values=c("#0072B2", # Dark blue
                                 "#D55E00", # Orange-brown
                                 "#CC79A7", # Pink
                                 "#009E73", # Green
                                 "#56B4E9", # Light blue
                                 "#E69F00")) # Yellow
  print(g)
  ggsave(
    glue("model1_lag{tau}_betawithin_{direc}_aic{undoc_flag}{rolling_flag}.pdf"),
    path=output_path, width = 10.8, height = 6.62, units = "in")
}

#### With model selection (BIC) ####
tbl_beta = tibble(date = as.Date(NA), betas=numeric(0), code=character(0))

# Find the estimates of beta per region over time
for (region in regions){
  betas = vector("double")
  dates = vector("character")
  
  # Select only the data for the relevant region
  data = df_long %>% filter(code == !!region)
  
  for (t in window_size:nrow(data)){
    # Use BIC for model selection - scope says we want to always keep
    # beta_within in
    if (rolling) {
      model = step(lm(fm, data=data[(t-window_size+1):t, ]), k=log(window_size), trace=0,
                   scope=list(
                     "lower" = glue("{infective_variable} ~ ",
                                    "lag({infective_variable}, {tau}):",
                                    "lag(susceptibleRate, {tau})") %>%
                       as.formula,
                     "upper" = fm))
    } else {
      model = step(lm(fm, data=head(data, t)), k=log(t), trace=0,
                   scope=list(
                     "lower" = glue("{infective_variable} ~ ",
                                    "lag({infective_variable}, {tau}):",
                                    "lag(susceptibleRate, {tau})") %>%
                       as.formula,
                     "upper" = fm))
      
    }
    # Retrieve the beta estimate and append this to the list of betas
    beta = coef(model)[[glue("lag({infective_variable}, {tau}):",
                              "lag(susceptibleRate, {tau})")]]
    betas = c(betas, beta)
  }
  
  # Append the results to the table
  tbl_beta = tbl_beta %>%
    bind_rows(tibble(date = data$date[window_size:nrow(data)],
                     betas = betas,
                     code = region) %>%
                drop_na())
}

# Add region and direction to the table
tbl_beta = tbl_beta %>%
  left_join(df_meta %>% select(c(region, code, direction)), by="code")

# Make a plot per direction
for (sub_tbl in split(tbl_beta, tbl_beta$direction)){
  direc = sub_tbl$direction[1]
  g = ggplot(sub_tbl, aes(date, betas, colour = region)) + 
    geom_point() +
    geom_smooth(method="loess", span=0.3, se=FALSE) +
    xlab("") +
    ylab(TeX("$\\beta_{within}$")) +
    scale_colour_manual(values=c("#0072B2", # Dark blue
                                 "#D55E00", # Orange-brown
                                 "#CC79A7", # Pink
                                 "#009E73", # Green
                                 "#56B4E9", # Light blue
                                 "#E69F00")) # Yellow
  print(g)
  ggsave(
    glue("model1_lag{tau}_betawithin_{direc}_bic{undoc_flag}{rolling_flag}.pdf"),
    path=output_path, width = 10.8, height = 6.62, units = "in")
}
