# https://rpubs.com/choisy/sir
rm(list=ls())

#### SIR Model ####
#### Setup ####
# Import standard variables
source("config.R")

# Import packages
library(deSolve)
library(glue)
library(tidyverse)

# Import data
df_long = readr::read_csv(path_full_long, col_types = do.call(
  cols, list(date = col_date(format = "%Y-%m-%d"))))
df_wide = readr::read_csv(path_full_wide, col_types = do.call(
  cols, list(date=col_date(format="%Y-%m-%d"))))

# Add nationwide variables by summing the individual regions' variables
infective_variable = "infectivesQuadratic"

susceptiblePopulationNational = df_wide %>%
  select(ends_with("susceptiblePopulation")) %>%
  rowSums
totalPopulationNational = df_wide %>%
  select(ends_with("totalPopulation")) %>%
  rowSums
infectivesNational = df_wide %>% 
  select(ends_with(glue("_{infective_variable}"))) %>% 
  rowSums
recoveredNational = df_wide %>% 
  select(ends_with(glue("_recovered"))) %>% 
  rowSums
areaNational = df_wide %>%
  select(ends_with("area")) %>%
  rowSums
df_wide = df_wide %>%
  mutate(susceptiblePopulationNational = susceptiblePopulationNational,
         susceptibleRateNational =
           susceptiblePopulationNational/totalPopulationNational,
         infectivesNational = infectivesNational,
         infectivesRateNational = infectivesNational / totalPopulationNational,
         recoveredNational = recoveredNational,
         totalPopulationNational = totalPopulationNational,
         populationDensityNational = totalPopulationNational/areaNational)

df_long = df_long %>%
  mutate(recoveredRate = recovered/totalPopulation)

####
sir_equations <- function(time, variables, parameters) {
  with(as.list(c(variables, parameters)), {
    dS = -beta * I * S
    dI = beta * I * S - gamma * I
    dR = gamma * I
    return(list(c(dS, dI, dR)))
  })
}

df = df_long %>% filter(code == "LOM")
infective_variable = "infectives"
susceptible_variable = "susceptiblePopulation"
recovered_variable = "recovered"

# We need to start with a non-zero value of I, so we drop the rows till we
# encounter the first infectives
df = df[which(df[[infective_variable]] > 0)[1]:nrow(df), ]

initial_values <- c(
  S = df[[susceptible_variable]][1], # number of susceptibles at time = 0
  I = df[[infective_variable]][1], # number of infectious at time = 0
  R = df[[recovered_variable]][1] # number of recovered (and immune) at time = 0
)

time_values <- seq(0, nrow(df)) # days

sir_1 <- function(beta, gamma, S0, I0, R0, times) {
  require(deSolve) # for the "ode" function
  
  # the differential equations:
  sir_equations <- function(time, variables, parameters) {
    with(as.list(c(variables, parameters)), {
      dS <- -beta * I * S
      dI <-  beta * I * S - gamma * I
      dR <-  gamma * I
      return(list(c(dS, dI, dR)))
    })
  }
  
  # the parameters values:
  parameters_values <- c(beta  = beta, gamma = gamma)
  
  # the initial values of variables:
  initial_values <- c(S = S0, I = I0, R = R0)
  
  # solving
  out <- ode(initial_values, times, sir_equations, parameters_values)
  
  # returning the output:
  return(as.data.frame(out))
}

ss <- function(beta, gamma, data=df){
  N = data$totalPopulation[1]
  I0 = data[["infectives"]][1]
  S0 = N - I0
  R0 = data[["recovered"]][1]
  times = seq(0, nrow(data)-1)
  predictions = sir_1(beta = beta, gamma = gamma,    # parameters
                      S0 = S0, I0 = I0, R0 = R0,     # variables initial values
                      times = times)                 # time points
  return(sum((predictions$I[-1] - data[["infectives"]][-1])^2))
}

ss_print <- function(beta, gamma, data=df){
  N = data$totalPopulation[1]
  I0 = data[["infectives"]][1]
  S0 = N - I0
  R0 = data[["recovered"]][1]
  times = seq(0, nrow(data)-1)
  predictions = sir_1(beta = beta, gamma = gamma,    # parameters
                      S0 = S0, I0 = I0, R0 = R0,     # variables initial values
                      times = times)                 # time points
  print(predictions)
  return(sum((predictions$I[-1] - data[["infectives"]][-1])^2))
}

# Given gamma = 0.5, compute the beta for which the SSR is minimal
beta_val = seq(from = 2, to = 3, le = 500)
ss_val = sapply(beta_val, ss, gamma = 0.5)
beta_hat = beta_val[ss_val == min(ss_val)]
plot(beta_val, ss_val, type = "l", lwd = 2,
     xlab = expression(paste("infectious contact rate ", beta)),
     ylab = "sum of squares")
abline(h = min(ss_val), lty = 2, col = "grey")
abline(v = beta_hat, lty = 2, col = "grey")

ss_print(beta_hat, gamma = 0.5)

# Given beta = beta_hat, compute the gamma for which the SSR is minimal
gamma_val = seq(from = 0.01, to = 0.9, le = 100)
ss_val = sapply(gamma_val,
                 function(x){ss(beta_hat, x)})
plot(gamma_val, ss_val, type = "l", lwd = 2,
     xlab = expression(paste("recovery rate ", gamma)),
     ylab = "sum of squares")
abline(h = min(ss_val), lty = 2, col = "grey")
abline(v = gamma_val[ss_val == min(ss_val)], lty = 2, col = "grey")

# Try both at the same time
n <- 20 # number of parameter values to try
beta_val <- seq(from = 0.01, to = 3, le = n)
gamma_val <- seq(from = 0.01, to = 8, le = n)
param_val <- expand.grid(beta_val, gamma_val)
ss_val <- with(param_val, Map(ss, Var1, Var2))
ss_val <- matrix(unlist(ss_val), n)
persp(beta_val, gamma_val, ss_val, theta = 40, phi = 30,
      xlab = "beta", ylab = "gamma", zlab = "sum of squares")


model_fit <- function(beta, gamma, data, ...) {
  I0 <- data$cases[1] # initial number of infected (from data)
  times <- data$day   # time points (from data)
  # model's predictions:
  predictions <- sir_1(beta = beta, gamma = gamma,   # parameters
                       S0 = N - I0, I0 = I0, R0 = 0, # variables' initial values
                       times = times)                # time points
  # plotting the observed prevalences:
  with(data, plot(day, cases, ...))
  # adding the model-predicted prevalence:
  with(predictions, lines(time, I, col = "red"))
}
