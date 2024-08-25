# Load Libraries ----
library(knitr)
library(data.table)
library(MASS)
library(rjags)
library(plyr)
library(stringr)
library(tictoc)

library(caret)
library(bayesplot)
library(BayesFactor)
library(brms)
library(rstanarm)

library(tidyverse)

# Read in data ----
## Clean data ----
Stormdata <- fread("E2_data.csv")
Stormdata$basin <- ifelse(Stormdata$basin == "atlantic", 1, 0)

## Create training and test data sets ----
training_data <- Stormdata |> filter(complete.cases(VMAX))
test_data <- Stormdata |> filter(!complete.cases(VMAX))


# 1. Model description ----
## Load data ----
Y <- training_data$VMAX
Z <- log(Y)
X <- training_data |>
  select(-c(
    StormID,
    Date,
    lead_time,
    VMAX,
    basin,
    LON,
    SST,
    CAPE1,
    SHTFL2,
    TCOND7002,
    CP1,
    VMAX_OP_T0
  ))
Xscale <- scale(X)
Xscale <- cbind("Intercept" = rep(1,length(Y)), Xscale)
HWRF <- training_data$HWRF

### Create StormID indicator for random effects model ----
StormIDs <- training_data |> select(StormID)
Storms <- unique(StormIDs$StormID)
IDs <- 1:length(Storms)
Storms <- data.frame(
  StormID = Storms,
  ID = IDs
)
StormIDs <- left_join(StormIDs, Storms)
StormID <- StormIDs$ID

n <- length(Y)
p <- ncol(Xscale)
N <- length(IDs)

burn <- 10000
iters <- 20000
chains <- 2
thin <- 5

# 2. Old Analysis ----
outputDF <- fread(file = "PollardTyler.csv")
colnames(outputDF) <- c("HWRF", "VMAXtlyer", "Ltyler", "Utyler")

ActualY <- read.csv("Actual Y.csv")

TotalY <- cbind(outputDF, ActualY)
TotalYpreds <- TotalY |> filter(!is.na(x))

MAD1 <- mean(abs(TotalYpreds$Y - TotalYpreds$x))
COV1 <- mean( (TotalYpreds$L <= TotalYpreds$x) & (TotalYpreds$x <= TotalYpreds$U))

MADtyler <- mean(abs(TotalYpreds$VMAXtlyer - TotalYpreds$x))
COVtyler <- mean( (TotalYpreds$Ltyler <= TotalYpreds$x) & (TotalYpreds$x <= TotalYpreds$Utyler))

MAD1
MADtyler

COV1
COVtyler


# 3. New Analysis ----
Y <- Stormdata$VMAX[complete.cases(Stormdata$VMAX)]
newDF <- cbind(Xscale, Y) |> as.data.frame()
fit1 <- generalTestBF(Y ~ .,
                      data = newDF)
fit2 <- stan_glm(Y ~ .,
                 data = newDF)
fit2
summary(fit2, digits = 3, probs = c(0.025, 0.975))
car::vif(fit2)

pp_check(fit2, nreps = 1000)
ppc_bars(Y, posterior_predict(fit2))

fit2 <- stan_glm(log(Y) ~ .,
                 data = newDF,
                 family = gaussian())
fit2
summary(fit2, digits = 3, probs = c(0.025, 0.975))
car::vif(fit2)

