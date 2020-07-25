rm(list=ls())

#### SIR Model ####
#### Setup ####
# Import standard variables
source("config.R")

# Import packages
library(deSolve)
library(glue)
library(latex2exp)
library(restriktor)
library(snakecase)
library(tidyverse)
library(xtable)

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
areaNational = df_wide %>%
  select(ends_with("area")) %>%
  rowSums
df_wide = df_wide %>%
  mutate(susceptibleRateNational =
           susceptiblePopulationNational/totalPopulationNational,
         infectivesNational = infectivesNational,
         populationDensityNational = totalPopulationNational/areaNational)

# Beta = alpha and gamma = beta
## National
infective_variable = "infectivesNational"
Infected = df_wide[[infective_variable]]
day = 0:(nrow(df_wide)-1)
N = df_wide %>%
  select(ends_with("totalPopulation")) %>%
  rowSums() %>%
  .[1]

plot(df_wide$date, Infected, type = "l")

## Regional
infective_variable = "infectives"
df = df_long %>% filter(code == "ABR")
df = df[which(df$infectives > 0)[1]:nrow(df), ]
Infected = df[[infective_variable]]
day = 0:(nrow(df)-1)
N = df$totalPopulation[1]

SIR <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dS <- -beta/N * I * S
    dI <- beta/N * I * S - gamma * I
    dR <- gamma * I
    list(c(dS, dI, dR))
  })
}

init <- c(S = N-Infected[1], I = Infected[1], R = 0)

RSS <- function(parameters) {
  names(parameters) <- c("beta", "gamma")
  out <- ode(y = init, times = day, func = SIR, parms = parameters)
  fit <- out[ , 3]
  sum((Infected - fit)^2)
}

betas = numeric(0)
gammas = numeric(0)
for (beta in seq(0,1,0.1)) {
  for (gamma in seq(0,1,0.1)) {
    Opt = optim(c(beta, gamma), RSS, method = "L-BFGS-B", lower = c(0, 0),
                upper = c(1, 1))
    
    print(glue("Beta: {beta}, Gamma: {gamma}, RSS: {RSS}"))
    betas = c(betas, Opt$par[1])
    gammas = c(gammas, Opt$par[2])
  }
}

betas_seq = seq(0, 1, 0.1)
gammas_seq = seq(0, 1, 0.1)
beta_gamma_opt = tibble(
  beta = rep(betas_seq, each=length(gammas_seq)),
  gamma = rep(gammas_seq, times=length(betas_seq)),
  beta_opt = betas,
  gamma_opt = gammas
)

# TODO: Make heatmap of combinations of beta and gamma or something similar

betas_seq = seq(0,1,0.01)
gamma = 0.5
plot(betas_seq,
     sapply(lapply(betas_seq, function(x){c(x, gamma)}), RSS))

print(glue("Mean beta: {mean(betas)}; Mean gamma: {mean(gammas)}"))

Opt <- optim(c(0.7, 0.5), RSS, method = "L-BFGS-B", lower = c(0, 0),
             upper = c(1, 1)) 

Opt_par <- setNames(Opt$par, c("beta", "gamma"))
Opt_par

t <- 1:nrow(df_wide) # time in days

fit <- data.frame(ode(y = init, times = t, func = SIR, parms = Opt_par))


###############################

init = c(S = N-Infected[1],
         I = Infected[1],
         R = df_wide %>%
           select(ends_with("recovered")) %>%
           rowSums() %>%
           .[1])
plot(day, Infected)

SIR <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  ####edit 2; use equally scaled variables 
  with(par, { dS <- -beta * (S/N) * I
  dI <- beta * (S/N) * I - gamma * I
  dR <- gamma * I
  list(c(dS, dI, dR))
  })
}

SIR2 <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  ####
  #### use as change of variables variable
  #### const = (beta-gamma)
  #### delta = gamma/beta
  #### R0 = beta/gamma > 1 
  #### 
  #### beta-gamma = beta*(1-delta)
  #### beta-gamma = beta*(1-1/R0)
  #### gamma = beta/R0
  with(par, { 
    beta  <- const/(1-1/R0)  
    gamma <- const/(R0-1)  
    dS <- -(beta * (S/N)      ) * I 
    dI <-  (beta * (S/N)-gamma) * I 
    dR <-  (             gamma) * I
    list(c(dS, dI, dR))
  })
}

RSS.SIR2 <- function(parameters) {
  names(parameters) <- c("const", "R0")
  out <- ode(y = init, times = day, func = SIR2, parms = parameters)
  fit <- out[ , 3]
  RSS <- sum((Infected - fit)^2)
  return(RSS)
}

### plotting different values R0
RSS.SIR <- function(parameters) {
  names(parameters) <- c("beta", "gamma")
  out <- ode(y = init, times = day, func = SIR, parms = parameters)
  fit <- out[ , 3]
  RSS <- sum((Infected - fit)^2)
  return(RSS)
}

lower = c(0, 0)
upper = c(1, 1)  ###adjust limit because different scale 1/N

### edit: get a good starting condition
mod <- nlsLM(Infected ~ a*exp(b*day), 
           start = list(a = Infected[1],
                        b = log(Infected[2]/Infected[1])))
optimsstart <- c(2,1)*coef(mod)[2]

# use the ordinary exponential model to determine const = beta - gamma
const <- coef(mod)[2]

set.seed(12)
Opt <- optim(optimsstart, RSS.SIR, method = "L-BFGS-B", lower = lower,
             upper = upper,
             hessian = TRUE)
Opt

sigest <- sqrt(Opt$value/(length(Infected)-1))
solve(1/(2*sigest^2)*Opt$hessian) 

optimsstart <- c(coef(mod)[2],5)
lower = c(0, 1)
upper = c(1, 10^3)  ### adjust limit because we use R0 now which should be >1

set.seed(12)
Opt2 <- optim(optimsstart, RSS.SIR2, method = "L-BFGS-B",lower=lower, upper=upper,
              hessian = TRUE, control = list(maxit = 1000, 
                                             parscale = c(10^-3,1)))
Opt2

sigest <- sqrt(Opt2$value/(length(Infected)-1))
solve(1/(2*sigest^2)*Opt2$hessian)

const <- coef(mod)[2]
R0 <- 1.1

# graph
plot(-100,-100, xlim=c(0,131), ylim = c(1,N),
     ylab = "infected", xlab = "days", yaxt = "n")
# axis(2, las=2, at=c(0:5)*1000,
#      labels=c(expression(0),
#               expression(1000),
#               expression(2000),
#               expression(3000),
#               expression(4000),
#               expression(5000)))
# axis(2, at=rep(c(2:9),9)*rep(10^c(0:8),each=8), labels=rep("",8*9),tck=-0.02)
title(bquote(paste("scenario's for different ", R[0])), cex.main = 1)

# time
t <- seq(0,131,0.1)

# plot model with different R0
R0s = c(0.8,1.2,4)
for (R0 in R0s) {
  fit <- data.frame(ode(y = init, times = t, func = SIR2, parms = c(const,R0)))$I
  
  if (R0s == R0s[1]) {
    plot(t,fit, type="l")
  } else {
    lines.default(t,fit)
  }
  
  text(t[601],fit[601],
       bquote(paste(R[0], " = ",.(R0))),
       cex=0.7,pos=4)
}

# plot observations
points(day,Infected)
