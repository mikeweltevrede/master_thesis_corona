rm(list=ls())

#### Model 1 - Within-region spread ####
# I_rt = alpha_within*I_rt-tau*S_rt-tau + X_rt*delta + nu_rt

#### Setup ####
# Import standard variables
source("config.R")

# Import packages
library(glue)
library(latex2exp)
library(restriktor)
library(snakecase)
library(tidyverse)
library(xtable)

# Import data
df_long = readr::read_csv(path_full_long, col_types = do.call(
  cols, list(date = col_date(format = "%Y-%m-%d"))))

pivot_to_df_wide = function(df_long) {
  # Turn long data into wide data
  df_wide = df_long %>%
    pivot_wider(names_from = code,
                values_from = all_of(colnames(df_long)[-c(1,2)]))
  
  # Column names are now of the form "variable_regionCode". We want to have them
  # in the form "regionCode_variable".
  colnames(df_wide) = colnames(df_wide) %>%
    sapply(function(s) {
      s %>%
        str_split("_", simplify = TRUE) %>%
        rev %>%
        paste(collapse="_")
    }, USE.NAMES=FALSE)
  
  # Order the columns alphabetically, keeping date at the start
  df_wide = df_wide[, c("date", sort(colnames(df_wide[-1])))]
  
  return(df_wide)
}

# Incubation period
lag = 5

# Select regressors
X_regressors = c("weekend")
all_variables = c("(Intercept)", paste0(X_regressors, 1))

# Should we apply the restriktor package to ensure that the alphas are positive?
restrict = FALSE

#### Data preprocessing ####
# Transform the variable to include undocumented infections, if applicable. Note
# that this does not depend on the region specifically but only the values at
# that moment of time. As such, we do not need to loop and can simply apply the
# function to each row.
form = "" %>%
  to_upper_camel_case

if (form %in% c("Linear", "Quadratic", "DownwardsVertex", "UpwardsVertex",
                "Cubic")){
  print(glue("####Running models while modelling undocumented infections with ",
             "the {form} functional form!####"))
  infective_variable = glue("infectivesRate{form}")
  undoc_flag = glue("_Undoc{form}")
  
} else if (form == ""){
  # Then do not use the undocumented infections modelling
  print("####Running models WITHOUT modelling undocumented infections!####")
  infective_variable = "infectivesRate"
  undoc_flag = ""
  
} else {
  sprintf(paste("The variable `form` is %s but it should be one of %s.",
                "Choosing 'infectives' as the infectives variable, so *not*",
                "modelling undocumented infections."),
          form,
          paste(c("Linear", "Quadratic", "DownwardsVertex",
                  "UpwardsVertex", "Cubic", ""), collapse=", ")) %>%
    print

  infective_variable = "infectivesRate"
  undoc_flag = ""
}

#### Adapt the lag to 1 and make the lagged dependent variable ####
df_wide[["infectivesRateNationalLag"]] = df_wide$infectivesRateNational -
  lag(df_wide$infectivesRateNational, 1)

df_long = df_long %>%
  group_by(code) %>%
  mutate(infectivesRateLag :=
           !!sym(infective_variable) - lag(!!sym(infective_variable), 1)) %>%
  ungroup()

# Get all region abbreviations
regions = df_long$code %>% unique

#### Models without model selection ####
results_table = tibble(variables = c(all_variables, "alpha"))

#### National model ####
# Construct formula
model_summ = glue("infectivesNational ~ ",
          "lag(infectivesNational, {lag}):lag(susceptibleRateNational, {lag})+",
          paste(X_regressors, collapse="+")) %>%
  as.formula %>%
  lm(data=df_wide) %>%
  summary
print(model_summ)

model_summ = glue("infectivesNational ~ ",
          "lag(infectivesNational, {lag}):lag(susceptibleRateNational, {lag})+",
          "lag(infectivesNational, {lag})+",
          paste(X_regressors, collapse="+")) %>%
  as.formula %>%
  lm(data=df_wide) %>%
  summary
print(model_summ)
print(glue("R0: {-coef(model_summ)[4,1]/coef(model_summ)[2,1]}"))

model_summ = glue("infectivesRateNational ~ ",
          "lag(infectivesRateNational, {lag}):lag(susceptibleRateNational, {lag})+",
          paste(X_regressors, collapse="+")) %>%
  as.formula %>%
  lm(data=df_wide) %>%
  summary
print(model_summ)

model_summ = glue(
  "infectivesRateNational ~ ",
  "lag(infectivesRateNational, {lag}):lag(susceptibleRateNational, {lag})+",
  "lag(infectivesRateNational, {lag})+",
  paste(X_regressors, collapse="+")) %>%
  as.formula %>%
  lm(data=df_wide) %>%
  summary
print(model_summ)
print(glue("R0: {-coef(model_summ)[4,1]/coef(model_summ)[2,1]}"))
 
# Estimate the model by OLS
model = lm(fm, data=df_wide)
summary(model)

# pdf(glue("{output_path}/model1_lag{lag}_lmplot_national{undoc_flag}.pdf"))
# par(mfrow=c(2,2))
# plot(model)
# par(mfrow=c(1,1))
# dev.off()

# TODO: Check the maths behind this
if (restrict) {
  model = restriktor(model,
                     constraints = c(0, 1, 0), # alpha_within > 0
                     rhs = 0)
}

summary(model)

# Retrieve parameter estimates
estimates = coef(summary(model))[, "Estimate"]
pvals = coef(summary(model))[, "Pr(>|t|)"] # TODO: SE with stars

# Insert parameter estimates in the results table
results_table = results_table %>%
  left_join(tibble("variables" = c(all_variables, "alpha"),
                   "National" = unname(
                     c(estimates[all_variables], estimates[
                       glue("lag(infectivesRateNational, {lag}):",
                            "lag(susceptibleRateNational, {lag})")])),
                   "National_pvals" = unname(
                     c(pvals[all_variables], pvals[
                       glue("lag(infectivesRateNational, {lag}):",
                            "lag(susceptibleRateNational, {lag})")]))),
            by="variables")

#### Regional models ####
# Construct formula
# TODO: Make left-hand side lagged
fm = glue("{infective_variable} ~ ",
          "lag({infective_variable}, {lag}):lag(susceptibleRate, {lag})+",
          paste(X_regressors, collapse="+")) %>%
  as.formula

for (region in regions){
  # Select only the data for the relevant region
  data = df_long %>% filter(code == !!region)
  
  # Estimate the model by OLS
  model = lm(fm, data=data)
  
  # pdf(glue("{output_path}/model1_lag{lag}_lmplot_{region}{undoc_flag}.pdf"))
  # par(mfrow=c(2,2))
  # plot(model)
  # par(mfrow=c(1,1))
  # dev.off()
  
  # TODO: Check the maths behind this
  if (restrict) {
    model = restriktor(model,
                       constraints = c(0, 1, 0), # alpha_within > 0
                       rhs = 0)
  }
  
  # Retrieve parameter estimates
  estimates = coef(summary(model))[, "Estimate"]
  pvals = coef(summary(model))[, "Pr(>|t|)"] # TODO: SE with stars
  
  # Insert parameter estimates in the results table
  results_table = results_table %>%
    left_join(tibble("variables" = c(all_variables, "alpha"),
                     !!glue("{region}") := unname(
                       c(estimates[all_variables], estimates[
                         glue("lag({infective_variable}, {lag}):",
                              "lag(susceptibleRate, {lag})")])),
                     !!glue("{region}_pvals") := unname(
                       c(pvals[all_variables], pvals[
                         glue("lag({infective_variable}, {lag}):",
                              "lag(susceptibleRate, {lag})")]))),
              by="variables")
}

# Transpose and reorder columns
results_table = results_table %>%
  gather(region, val, 2:ncol(results_table)) %>%
  spread(names(results_table)[1], val) %>%
  select(region, alpha, everything())

# Put the national results at the top
results_table = rbind(
  results_table %>% filter(str_detect(region, "National")),
  results_table %>% filter(!str_detect(region, "National"))) %>%
  column_to_rownames("region")

# Return a LaTeX table
table_no_ms = xtable(results_table, math.style.exponents = TRUE)
table_no_ms

#### Models with model selection (AIC) ####
results_table_aic = tibble(variables = c(all_variables, "alpha"))

#### National model ####
# Construct formula
fm = glue("infectivesNational ~ ",
          "lag(infectivesNational, {lag}):lag(susceptibleRateNational, {lag})+",
          paste(X_regressors, collapse="+")) %>%
  as.formula

# Use AIC for model selection - scope says we want to always keep
# alpha_within in
model = step(lm(fm, data=df_wide), k=2, trace=0,
             scope=list("lower" = glue("infectivesNational ~ ",
                                       "lag(infectivesNational, {lag}):",
                                       "lag(susceptibleRateNational, {lag})") %>%
                          as.formula,
                        "upper" = fm))

# pdf(glue("{output_path}/model1_lag{lag}_lmplot_national_aic{undoc_flag}.pdf"))
# par(mfrow=c(2,2))
# plot(model)
# par(mfrow=c(1,1))
# dev.off()

# TODO: Check the maths behind this
if (restrict) {
  model = restriktor(model,
                     constraints = c(0, 1, 0), # alpha_within > 0
                     rhs = 0)
}

# Retrieve parameter estimates
estimates = coef(summary(model))[, "Estimate"]
pvals = coef(summary(model))[, "Pr(>|t|)"] # TODO: SE with stars

# Insert parameter estimates in the results table
results_table_aic = results_table_aic %>%
  left_join(tibble("variables" = c(all_variables, "alpha"),
                   "National" = unname(
                     c(estimates[all_variables], estimates[
                       glue("lag(infectivesNational, {lag}):",
                            "lag(susceptibleRateNational, {lag})")])),
                   "National_pvals" = unname(
                     c(pvals[all_variables], pvals[
                       glue("lag(infectivesNational, {lag}):",
                            "lag(susceptibleRateNational, {lag})")]))),
            by="variables")

#### Regional models ####
# Construct formula
fm = glue("{infective_variable} ~ ",
          "lag({infective_variable}, {lag}):lag(susceptibleRate, {lag})+",
          paste(X_regressors, collapse="+")) %>%
  as.formula

for (region in regions){
  # Select only the data for the relevant region
  data = df_long %>% filter(code == !!region)
  
  # Use AIC for model selection
  model = step(lm(fm, data=data), k=2, trace=0,
               scope=list("lower" = glue("{infective_variable} ~ ",
                                         "lag({infective_variable}, {lag}):",
                                         "lag(susceptibleRate, {lag})") %>%
                            as.formula,
                          "upper" = fm))
  
  # pdf(glue("{output_path}/model1_lag{lag}_lmplot_{region}_aic{undoc_flag}.pdf"))
  # par(mfrow=c(2,2))
  # plot(model)
  # par(mfrow=c(1,1))
  # dev.off()
  
  # TODO: Check the maths behind this
  if (restrict) {
    model = restriktor(model,
                       constraints = c(0, 1, 0), # alpha_within > 0
                       rhs = 0)
  }
  
  # Retrieve parameter estimates
  estimates = coef(summary(model))[, "Estimate"]
  pvals = coef(summary(model))[, "Pr(>|t|)"] # TODO: SE with stars
  
  # Insert parameter estimates in the results table
  results_table_aic = results_table_aic %>%
    left_join(tibble("variables" = c(all_variables, "alpha"),
                     !!glue("{region}") := unname(
                       c(estimates[all_variables], estimates[
                         glue("lag({infective_variable}, {lag}):",
                              "lag(susceptibleRate, {lag})")])),
                     !!glue("{region}_pvals") := unname(
                       c(pvals[all_variables], pvals[
                         glue("lag({infective_variable}, {lag}):",
                              "lag(susceptibleRate, {lag})")]))),
              by="variables")
}

# Transpose and reorder columns
results_table_aic = results_table_aic %>%
  gather(region, val, 2:ncol(results_table_aic)) %>%
  spread(names(results_table_aic)[1], val) %>%
  select(region, alpha, everything())

# Put the national results at the top
results_table_aic = rbind(
  results_table_aic %>% filter(str_detect(region, "National")),
  results_table_aic %>% filter(!str_detect(region, "National"))) %>%
  column_to_rownames("region")

# Return a LaTeX table
table_ms_aic = xtable(results_table_aic, math.style.exponents = TRUE)
table_ms_aic

#### Models with model selection (BIC) ####
results_table_bic = tibble(variables = c(all_variables, "alpha"))

#### National model ####
# Construct formula
fm = glue("infectivesNational ~ ",
          "lag(infectivesNational, {lag}):lag(susceptibleRateNational, {lag})+",
          paste(X_regressors, collapse="+")) %>%
  as.formula

# Use BIC for model selection - scope says we want to always keep
# alpha_within in
model = step(lm(fm, data=df_wide), k=log(nrow(df_wide)), trace=0,
             scope=list("lower" = glue("infectivesNational ~ ",
                                       "lag(infectivesNational, {lag}):",
                                       "lag(susceptibleRateNational, {lag})") %>%
                          as.formula,
                        "upper" = fm))

# pdf(glue("{output_path}/model1_lag{lag}_lmplot_national_bic{undoc_flag}.pdf"))
# par(mfrow=c(2,2))
# plot(model)
# par(mfrow=c(1,1))
# dev.off()

# TODO: Check the maths behind this
if (restrict) {
  model = restriktor(model,
                     constraints = c(0, 1, 0), # alpha_within > 0
                     rhs = 0)
}

# Retrieve parameter estimates
estimates = coef(summary(model))[, "Estimate"]
pvals = coef(summary(model))[, "Pr(>|t|)"] # TODO: SE with stars

# Insert parameter estimates in the results table
results_table_bic = results_table_bic %>%
  left_join(tibble("variables" = c(all_variables, "alpha"),
                   "National" = unname(
                     c(estimates[all_variables], estimates[
                       glue("lag(infectivesNational, {lag}):",
                            "lag(susceptibleRateNational, {lag})")])),
                   "National_pvals" = unname(
                     c(pvals[all_variables], pvals[
                       glue("lag(infectivesNational, {lag}):",
                            "lag(susceptibleRateNational, {lag})")]))),
            by="variables")

#### Regional models ####
# Construct formula
fm = glue("{infective_variable} ~ ",
          "lag({infective_variable}, {lag}):lag(susceptibleRate, {lag})+",
          paste(X_regressors, collapse="+")) %>%
  as.formula

for (region in regions){
  # Select only the data for the relevant region
  data = df_long %>% filter(code == !!region)
  
  # Use BIC for model selection
  model = step(lm(fm, data=data), k=log(nrow(data)), trace=0,
               scope=list("lower" = glue("{infective_variable} ~ ",
                                         "lag({infective_variable}, {lag}):",
                                         "lag(susceptibleRate, {lag})") %>%
                            as.formula,
                          "upper" = fm))
  
  # pdf(glue("{output_path}/model1_lag{lag}_lmplot_{region}_bic{undoc_flag}.pdf"))
  # par(mfrow=c(2,2))
  # plot(model)
  # par(mfrow=c(1,1))
  # dev.off()
  
  # TODO: Check the maths behind this
  if (restrict) {
    model = restriktor(model,
                       constraints = c(0, 1, 0), # alpha_within > 0
                       rhs = 0)
  }
  
  # Retrieve parameter estimates
  estimates = coef(summary(model))[, "Estimate"]
  pvals = coef(summary(model))[, "Pr(>|t|)"] # TODO: SE with stars
  
  # Insert parameter estimates in the results table
  results_table_bic = results_table_bic %>%
    left_join(tibble("variables" = c(all_variables, "alpha"),
                     !!glue("{region}") := unname(
                       c(estimates[all_variables], estimates[
                         glue("lag({infective_variable}, {lag}):",
                              "lag(susceptibleRate, {lag})")])),
                     !!glue("{region}_pvals") := unname(
                       c(pvals[all_variables], pvals[
                         glue("lag({infective_variable}, {lag}):",
                              "lag(susceptibleRate, {lag})")]))),
              by="variables")
}

# Transpose and reorder columns
results_table_bic = results_table_bic %>%
  gather(region, val, 2:ncol(results_table_bic)) %>%
  spread(names(results_table_bic)[1], val) %>%
  select(region, alpha, everything())

# Put the national results at the top
results_table_bic = rbind(
  results_table_bic %>% filter(str_detect(region, "National")),
  results_table_bic %>% filter(!str_detect(region, "National"))) %>%
  column_to_rownames("region")

# Return a LaTeX table
table_ms_bic = xtable(results_table_bic, math.style.exponents = TRUE)
table_ms_bic

#### Plot alpha over time ####
df_meta = readxl::read_xlsx(path_metadata, sheet = "Metadata")

# Starting index - we want at least this number of observations
start = 50

fm = glue("{infective_variable} ~ ",
          "lag({infective_variable}, {lag}):lag(susceptibleRate, {lag})+",
          paste(X_regressors, collapse="+")) %>%
  as.formula

#### Without model selection ####
tbl_alpha = tibble(date = as.Date(NA), alphas=numeric(0), code=character(0))

# Find the estimates of alpha per region over time
for (region in regions){
  alphas = vector("double")
  dates = vector("character")
  
  # Select only the data for the relevant region
  data = df_long %>% filter(code == !!region)
  
  for (t in start:nrow(data)){
    # Estimate the model by OLS
    model = lm(fm, data=data[1:t, ])
    
    # TODO: Check the maths behind this
    if (restrict) {
      model = restriktor(model,
                         constraints = c(0, 1, 0), # alpha_within > 0
                         rhs = 0)
    }
    
    # Retrieve the alpha estimate and append this to the list of alphas
    alpha = coef(model)[[glue("lag({infective_variable}, {lag}):",
                              "lag(susceptibleRate, {lag})")]]
    alphas = c(alphas, alpha)
  }
  
  # Append the results to the table
  tbl_alpha = tbl_alpha %>%
    bind_rows(tibble(date = data$date[start:nrow(data)],
                     alphas = alphas,
                     code = region) %>%
                drop_na())
}

# Add region and direction to the table
tbl_alpha = tbl_alpha %>%
  left_join(df_meta %>% select(c(region, code, direction)), by="code")

# Make a plot per direction
for (sub_tbl in split(tbl_alpha, tbl_alpha$direction)){
  direc = sub_tbl$direction[1]
  g = ggplot(sub_tbl, aes(date, alphas, colour = region)) + 
    geom_point() +
    geom_smooth(method="loess", span=0.3, se=FALSE) +
    xlab("") +
    ylab(TeX("$\\alpha_{within}$")) +
    scale_colour_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#0072B2",
                                 "#D55E00", "#CC79A7"))
  print(g)
  ggsave(glue("model1_lag{lag}_alphawithin_{direc}{undoc_flag}.pdf"),
         path=output_path, width = 10.8, height = 6.62, units = "in")
}

#### With model selection (AIC) ####
tbl_alpha = tibble(date = as.Date(NA), alphas=numeric(0), code=character(0))

# Find the estimates of alpha per region over time
for (region in regions){
  alphas = vector("double")
  dates = vector("character")
  
  # Select only the data for the relevant region
  data = df_long %>% filter(code == !!region)
  
  for (t in start:nrow(data)){
    # Use AIC for model selection - scope says we want to always keep
    # alpha_within in
    model = step(lm(fm, data=data[1:t, ]), k=2, trace=0,
                 scope=list("lower" = glue("{infective_variable} ~ ",
                                           "lag({infective_variable}, {lag}):",
                                           "lag(susceptibleRate, {lag})") %>%
                              as.formula,
                            "upper" = fm))
    
    # TODO: Check the maths behind this
    if (restrict) {
      model = restriktor(model,
                         constraints = c(0, 1, 0), # alpha_within > 0
                         rhs = 0)
    }
    
    # Retrieve the alpha estimate and append this to the list of alphas
    alpha = coef(model)[[glue("lag({infective_variable}, {lag}):",
                              "lag(susceptibleRate, {lag})")]]
    alphas = c(alphas, alpha)
  }
  
  # Append the results to the table
  tbl_alpha = tbl_alpha %>%
    bind_rows(tibble(date = data$date[start:nrow(data)],
                     alphas = alphas,
                     code = region) %>%
                drop_na())
}

# Add region and direction to the table
tbl_alpha = tbl_alpha %>%
  left_join(df_meta %>% select(c(region, code, direction)), by="code")

# Make a plot per direction
for (sub_tbl in split(tbl_alpha, tbl_alpha$direction)){
  direc = sub_tbl$direction[1]
  g = ggplot(sub_tbl, aes(date, alphas, colour = region)) + 
    geom_point() +
    geom_smooth(method="loess", span=0.3, se=FALSE) +
    xlab("") +
    ylab(TeX("$\\alpha_{within}$")) +
    scale_colour_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#0072B2",
                                 "#D55E00", "#CC79A7"))
  print(g)
  ggsave(glue("model1_lag{lag}_alphawithin_{direc}_aic{undoc_flag}.pdf"),
         path=output_path, width = 10.8, height = 6.62, units = "in")
}

#### With model selection (BIC) ####
tbl_alpha = tibble(date = as.Date(NA), alphas=numeric(0), code=character(0))

# Find the estimates of alpha per region over time
for (region in regions){
  alphas = vector("double")
  dates = vector("character")
  
  # Select only the data for the relevant region
  data = df_long %>% filter(code == !!region)
  
  for (t in start:nrow(data)){
    # Use BIC for model selection - scope says we want to always keep
    # alpha_within in
    model = step(lm(fm, data=data[1:t, ]), k=log(t), trace=0,
                 scope=list("lower" = glue("{infective_variable} ~ ",
                                           "lag({infective_variable}, {lag}):",
                                           "lag(susceptibleRate, {lag})") %>%
                              as.formula,
                            "upper" = fm))
    # TODO: Check the maths behind this
    if (restrict) {
      model = restriktor(model,
                         constraints = c(0, 1, 0), # alpha_within > 0
                         rhs = 0)
    }
    
    # Retrieve the alpha estimate and append this to the list of alphas
    alpha = coef(model)[[glue("lag({infective_variable}, {lag}):",
                              "lag(susceptibleRate, {lag})")]]
    alphas = c(alphas, alpha)
  }
  
  # Append the results to the table
  tbl_alpha = tbl_alpha %>%
    bind_rows(tibble(date = data$date[start:nrow(data)],
                     alphas = alphas,
                     code = region) %>%
                drop_na())
}

# Add region and direction to the table
tbl_alpha = tbl_alpha %>%
  left_join(df_meta %>% select(c(region, code, direction)), by="code")

# Make a plot per direction
for (sub_tbl in split(tbl_alpha, tbl_alpha$direction)){
  direc = sub_tbl$direction[1]
  g = ggplot(sub_tbl, aes(date, alphas, colour = region)) + 
    geom_point() +
    geom_smooth(method="loess", span=0.3, se=FALSE) +
    xlab("") +
    ylab(TeX("$\\alpha_{within}$")) +
    scale_colour_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#0072B2",
                                 "#D55E00", "#CC79A7"))
  print(g)
  ggsave(glue("model1_lag{lag}_alphawithin_{direc}_bic{undoc_flag}.pdf"),
         path=output_path, width = 10.8, height = 6.62, units = "in")
}
