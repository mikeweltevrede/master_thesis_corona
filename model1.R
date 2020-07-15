#### Model 1 - Within-region spread ####
# We start with a simple model ignoring effects across regions:
# I_rt = alpha_within*I_rt-tau*S_rt-tau + X_rt*delta + nu_rt

#### Setup ####
# Import standard variables
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
  cols, list(date=col_date(format="%Y-%m-%d"))))

# Incubation period
lag = 5

X_regressors = c("weekend", "weekNumber")
all_variables = c("(Intercept)", X_regressors) %>%
  str_replace("weekend", "weekend1")

# Get all region abbreviations
regions = df_long$code %>% unique

#### Data preprocessing ####
# Transform the variable to include undocumented infections, if applicable. Note
# that this does not depend on the region specifically but only the values at
# that moment of time. As such, we do not need to loop and can simply apply the
# function to each row.
form = "quadratic" %>%
  to_upper_camel_case

if (form %in% c("Linear", "Quadratic", "DownwardsVertex", "UpwardsVertex",
                "Cubic")){
  infective_variable = glue("infectives{form}") %>%
    snakecase::to_lower_camel_case()
  
  undoc_flag = glue("_Undoc{form}")
  
} else if (form == ""){
  # Then do not use the undocumented infections modelling
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

# Add nationwide variables by summing the individual regions' variables
susceptibleTotal = df_wide %>%
  select(ends_with("susceptiblePopulation")) %>%
  rowSums
total = df_wide %>%
  select(ends_with("totalPopulation")) %>%
  rowSums
infectivesTotal = df_wide %>% 
  select(ends_with(glue("_{infective_variable}"))) %>% 
  rowSums
df_wide = df_wide %>%
  mutate(susceptibleRateTotal = susceptibleTotal/total,
         infectivesTotal = infectivesTotal)

# Add weekend and weekday effect
df_wide = df_wide %>%
  mutate(weekNumber = lubridate::week(df_wide$date)) %>%
  mutate(weekend = lubridate::wday(df_wide$date, label = TRUE)
         %in% c("Sat", "Sun") %>% as.integer %>% as.factor)

df_long = df_long %>%
  mutate(weekNumber = lubridate::week(df_long$date)) %>%
  mutate(weekend = lubridate::wday(df_long$date, label = TRUE)
         %in% c("Sat", "Sun") %>% as.integer %>% as.factor)

#### Least Squares Dummy Variables (LSDV) regression ####
fm = glue("{infective_variable} ~ lag({infective_variable}, {lag}):",
          "lag(susceptibleRate, {lag})+",
           paste(X_regressors, collapse="+")) %>%
  paste("+factor(code)") %>%
  as.formula

model = lm(fm, data=df_long)
summary(model)

pdf(glue("{output_path}/model1_lag{lag}_lmplot_lsdv{undoc_flag}.pdf"))
par(mfrow=c(2,2))
plot(model)
par(mfrow=c(1,1))
dev.off()

#### Models without model selection ####
results_table = tibble(variables = c(all_variables, "alpha"))

#### National model ####
# Construct formula
fm = glue("infectivesTotal ~ lag(infectivesTotal, {lag}):",
          "lag(susceptibleRateTotal, {lag})+",
          paste(X_regressors, collapse="+")) %>%
  as.formula

# Estimate the model by OLS
model = lm(fm, data=df_wide)
summary(model)

pdf(glue("{output_path}/model1_lag{lag}_lmplot_national{undoc_flag}.pdf"))
par(mfrow=c(2,2))
plot(model)
par(mfrow=c(1,1))
dev.off()

# Retrieve parameter estimates
estimates = coef(summary(model))[, "Estimate"]
pvals = coef(summary(model))[, "Pr(>|t|)"] # TODO: SE with stars

# Insert parameter estimates in the results table
# TODO: Idem dito
results_table = results_table %>%
  left_join(tibble("variables" = c(all_variables, "alpha"),
                   "National" = unname(
                     c(estimates[all_variables], estimates[
                       glue("lag(infectivesTotal, {lag}):",
                            "lag(susceptibleRateTotal, {lag})")])),
                   "National_pvals" = unname(
                     c(pvals[all_variables], pvals[
                       glue("lag(infectivesTotal, {lag}):",
                            "lag(susceptibleRateTotal, {lag})")]))),
            by="variables")

#### Regional models ####
# Construct formula
fm = glue("{infective_variable} ~ lag({infective_variable}, {lag}):",
          "lag(susceptibleRate, {lag})+",
          paste(X_regressors, collapse="+")) %>%
  as.formula

for (region in regions){
  # Select only the data for the relevant region
  data = df_long %>% filter(code == !!region)
  
  # Estimate the model by OLS
  model = lm(fm, data=data)
  
  pdf(glue("{output_path}/model1_lag{lag}_lmplot_{region}{undoc_flag}.pdf"))
  par(mfrow=c(2,2))
  plot(model)
  par(mfrow=c(1,1))
  dev.off()
  
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

#### Models with model selection (BIC) ####
results_table_ms = tibble(variables = c(all_variables, "alpha"))

#### National model ####
# Construct formula
# TODO: Idem dito
fm = glue("infectivesTotal ~ lag(infectivesTotal, {lag}):",
          "lag(susceptibleRateTotal, {lag})+",
           paste(X_regressors, collapse="+")) %>%
  as.formula

# Use BIC for model selection - scope says we want to always keep
# alpha_within in
             scope=list("lower" = glue("infectivesTotal ~ ",
                                       "lag(infectivesTotal, {lag}):",
                                       "lag(susceptibleRateTotal, {lag})") %>%
                          as.formula,
                        "upper" = fm))

pdf(glue("{output_path}/model1_lag{lag}_lmplot_national_aic{undoc_flag}.pdf"))
par(mfrow=c(2,2))
plot(model)
par(mfrow=c(1,1))
dev.off()

# Retrieve parameter estimates
estimates = coef(summary(model))[, "Estimate"]
pvals = coef(summary(model))[, "Pr(>|t|)"] # TODO: SE with stars

# Insert parameter estimates in the results table
# TODO: Idem dito
results_table_ms_aic = results_table_ms_aic %>%
  left_join(tibble("variables" = c(all_variables, "alpha"),
                   "National" = unname(
                     c(estimates[all_variables], estimates[
                       glue("lag(infectivesTotal, {lag}):",
                            "lag(susceptibleRateTotal, {lag})")])),
                   "National_pvals" = unname(
                     c(pvals[all_variables], pvals[
                       glue("lag(infectivesTotal, {lag}):",
                            "lag(susceptibleRateTotal, {lag})")]))),
            by="variables")

#### Regional models ####
# Construct formula
fm = glue("{infective_variable} ~ lag({infective_variable}, {lag}):",
          "lag(susceptibleRate, {lag})+",
          paste(X_regressors, collapse="+")) %>%
  as.formula

for (region in regions){
  # Select only the data for the relevant region
  data = df_long %>% filter(code == !!region)
  
  # Use BIC for model selection
  # Use AIC for model selection
               scope=list("lower" = glue("{infective_variable} ~ ",
                                         "lag({infective_variable}, {lag}):",
                                         "lag(susceptibleRate, {lag})") %>%
                            as.formula,
                          "upper" = fm))
  
  pdf(glue("{output_path}/model1_lag{lag}_lmplot_{region}_aic{undoc_flag}.pdf"))
  par(mfrow=c(2,2))
  plot(model)
  par(mfrow=c(1,1))
  dev.off()
  
  # Retrieve parameter estimates
  estimates = coef(summary(model))[, "Estimate"]
  pvals = coef(summary(model))[, "Pr(>|t|)"] # TODO: SE with stars
  
  # Insert parameter estimates in the results table
  results_table_ms_aic = results_table_ms_aic %>%
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
results_table_ms_aic = results_table_ms_aic %>%
  gather(region, val, 2:ncol(results_table_ms_aic)) %>%
  spread(names(results_table_ms_aic)[1], val) %>%
  select(region, alpha, everything())

# Put the national results at the top
results_table_ms_aic = rbind(
  results_table_ms_aic %>% filter(str_detect(region, "National")),
  results_table_ms_aic %>% filter(!str_detect(region, "National"))) %>%
  column_to_rownames("region")

# Return a LaTeX table
table_ms_aic = xtable(results_table_ms_aic, math.style.exponents = TRUE)
table_ms_aic

#### Models with model selection (BIC) ####
results_table_ms_bic = tibble(variables = c(all_variables, "alpha"))

#### National model ####
# Construct formula
# TODO: Idem dito
fm = glue("infectivesTotal ~ lag(infectivesTotal, {lag}):",
          "lag(susceptibleRateTotal, {lag})+",
          paste(X_regressors, collapse="+")) %>%
  as.formula

# Use BIC for model selection - scope says we want to always keep
# alpha_within in
model = step(lm(fm, data=df_wide), k=log(nrow(df_wide)), trace=0,
             scope=list("lower" = glue("infectivesTotal ~ ",
                                       "lag(infectivesTotal, {lag}):",
                                       "lag(susceptibleRateTotal, {lag})") %>%
                          as.formula,
                        "upper" = fm))

pdf(glue("{output_path}/model1_lag{lag}_lmplot_national_bic{undoc_flag}.pdf"))
par(mfrow=c(2,2))
plot(model)
par(mfrow=c(1,1))
dev.off()

# Retrieve parameter estimates
estimates = coef(summary(model))[, "Estimate"]
pvals = coef(summary(model))[, "Pr(>|t|)"] # TODO: SE with stars

# Insert parameter estimates in the results table
# TODO: Idem dito
results_table_ms_bic = results_table_ms_bic %>%
  left_join(tibble("variables" = c(all_variables, "alpha"),
                   "National" = unname(
                     c(estimates[all_variables], estimates[
                       glue("lag(infectivesTotal, {lag}):",
                            "lag(susceptibleRateTotal, {lag})")])),
                   "National_pvals" = unname(
                     c(pvals[all_variables], pvals[
                       glue("lag(infectivesTotal, {lag}):",
                            "lag(susceptibleRateTotal, {lag})")]))),
            by="variables")

#### Regional models ####
# Construct formula
fm = glue("{infective_variable} ~ lag({infective_variable}, {lag}):",
          "lag(susceptibleRate, {lag})+",
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
  
  pdf(glue("{output_path}/model1_lag{lag}_lmplot_{region}_bic{undoc_flag}.pdf"))
  par(mfrow=c(2,2))
  plot(model)
  par(mfrow=c(1,1))
  dev.off()
  
  # Retrieve parameter estimates
  estimates = coef(summary(model))[, "Estimate"]
  pvals = coef(summary(model))[, "Pr(>|t|)"] # TODO: SE with stars
  
  # Insert parameter estimates in the results table
  results_table_ms_bic = results_table_ms_bic %>%
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
results_table_ms_bic = results_table_ms_bic %>%
  gather(region, val, 2:ncol(results_table_ms_bic)) %>%
  spread(names(results_table_ms_bic)[1], val) %>%
  select(region, alpha, everything())

# Put the national results at the top
results_table_ms_bic = rbind(
  results_table_ms_bic %>% filter(str_detect(region, "National")),
  results_table_ms_bic %>% filter(!str_detect(region, "National"))) %>%
  column_to_rownames("region")

# Return a LaTeX table
table_ms_bic = xtable(results_table_ms_bic, math.style.exponents = TRUE)
table_ms_bic

#### Plot alpha over time ####
df_meta = readxl::read_xlsx(path_wiki, sheet = "Metadata")

# Starting index - we want at least this number of observations
start = 50

fm = glue("{infective_variable} ~ lag({infective_variable}, {lag}):",
          "lag(susceptibleRate, {lag})+",
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
    
    # Retrieve the alpha estimate and append this to the list of alphas
    # TODO: Unhardcode the lag
    alpha = model$coefficients[[glue("lag({infective_variable}, {lag}):",
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

#### With model selection BIC ####
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
    model
    
    # Retrieve the alpha estimate and append this to the list of alphas
    alpha = model$coefficients[[glue("lag({infective_variable}, {lag}):",
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
