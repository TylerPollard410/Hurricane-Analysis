# Load Libraries ----
library(knitr)
library(data.table)
library(MASS)
library(rjags)
library(plyr)
library(stringr)
library(lubridate)
library(tictoc)
library(cowplot)
library(caret)
library(splines)
library(DescTools)
library(bayesplot)
library(BayesFactor)
library(rstanarm)
library(tidybayes)
library(brms)
library(performance)
library(tidyverse)

# Read in data ----
## Clean data ----
Stormdata_raw <- fread("~/Desktop/Hurricane Analysis/_data/E2_data.csv")
Stormdata <- Stormdata_raw |>
  mutate(
    StormID = factor(StormID),
    basin = factor(basin),
    Date = as_datetime(Date, tz = "UTC")
  )

# str(Stormdata)

# Timediff <- Stormdata |>
#   group_by(StormID) |>
#   summarise(
#     difftime(Date, lag(Date, default = Date[1]), units = "hours")
#   )
# 
# nondense <- which(Timediff$`difftime(Date, lag(Date, default = Date[1]), units = "hours")` > 6)
# nondenseID <- Stormdata |> slice(nondense)

## Create training and test data sets ----
StormdataTrain <- Stormdata |> filter(complete.cases(VMAX))
StormdataTest <- Stormdata |> filter(!complete.cases(VMAX))

## Transform data ----
# Remove not varying 
StormdataTrain2 <- StormdataTrain |>
  select(-lead_time)

# Create Date vars 
dataYears <- year(StormdataTrain2$Date)
dataMonths <- month(StormdataTrain2$Date, label = TRUE)
dataDays <- day(StormdataTrain2$Date)

StormdataTrain3 <- StormdataTrain2 |>
  mutate(
    Year = factor(dataYears, ordered = TRUE),
    Month = dataMonths
  ) |>
  group_by(StormID) |>
  mutate(
    StormElapsedTime = as.numeric(difftime(Date, min(Date), units = "hours"))
  ) |>
  select(
    StormID,
    Date,
    Year,
    Month,
    StormElapsedTime,
    everything()
  ) |>
  ungroup()

summary(StormdataTrain3$VMAX)

## Plot VMAX ----
### Histogram ----
ggplot(data = StormdataTrain3) +
  geom_histogram(
    aes(x = VMAX, after_stat(density)),
    color = "#99c7c7", fill = "#bcdcdc",
    bins = 100) +
  geom_density(#data = final_data3,
    aes(x = VMAX),
    color = "#007C7C", 
    linewidth = 1) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  #scale_x_continuous(limits = c(-0.5,1)) +
  labs(title = "Density Plot",
       subtitle = "Data",
       x = "VMAX",
       y = "Density") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

### Time ----
ggplot(data = StormdataTrain3) +
  geom_point(aes(x = Date, y = VMAX)) +
  scale_x_datetime(date_breaks = "month", 
                   date_minor_breaks = "day", 
                   date_labels = "%b-%Y") + 
  theme(
    axis.text.x = element_text(angle = 90)
  )


### Map ----
# Do by time next
world_coordinates <- map_data("world") 
ggplot() + 
  # geom_map() function takes world coordinates  
  # as input to plot world map 
  geom_map( 
    data = world_coordinates, map = world_coordinates, 
    aes(x = long, y = lat, map_id = region) 
  ) + 
  geom_point(
    data = StormdataTrain3,
    aes(x = LON-360, y = LAT, 
        color = VMAX)
  ) +
  xlim(c(-180,0)) +
  ylim(c(0,60)) +
  scale_color_continuous(low = "green", high = "red") +
  theme_bw()

# Fit model ----
## Create Data ----
StormdataTrain4 <- StormdataTrain3 |>
  select(
    "StormID",
    #Date,
    Year,
    Month,
    StormElapsedTime,
    "basin",
    "LAT",
    "LON",
    "MINSLP",
    "SHR_MAG",
    "STM_SPD",
    "SST",
    "RHLO",
    "CAPE1",
    "CAPE3",
    "SHTFL2",
    "TCOND7002",
    "INST2",
    "CP1",
    "TCONDSYM2",
    "COUPLSYM3",
    "HWFI",
    "VMAX_OP_T0",
    "HWRF",
    "VMAX"
  ) |>
  mutate(
    StormID = droplevels(StormID)
  )
str(StormdataTrain4)

linFit1 <- brm(
  formula = VMAX ~ 
    Year +
    Month +
    StormElapsedTime + 
    #basin + 
    LAT +
    LON +
    MINSLP +
    SHR_MAG +
    STM_SPD +
    SST +
    RHLO +
    CAPE1 +
    CAPE3 +
    SHTFL2 +
    TCOND7002 +
    INST2 +
    CP1 +
    TCONDSYM2 +
    COUPLSYM3 +
    HWFI +
    VMAX_OP_T0 +
    HWRF +
    (1|StormID),
  data = StormdataTrain4, 
  family = gaussian(), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52, 
  warmup = 1000
)
prior_summary(linFit1)
posterior_summary(linFit1)
linFit1

print(linFit1, digits = 4)
plot(linFit1)
loo(linFit1)
performance::check_distribution(linFit1)
performance::check_outliers(linFit1)
performance::check_heteroskedasticity(linFit1)
performance_rmse(linFit1)
performance_mae(linFit1)
variance_decomposition(linFit1)
fixef(linFit1)
ranef(linFit1)



StormdataTrain5 <- StormdataTrain3 |>
  select(
    "StormID",
    #Date,
    Year,
    Month,
    StormElapsedTime,
    "basin",
    "LAT",
    "LON",
    "MINSLP",
    "SHR_MAG",
    "STM_SPD",
    "SST",
    "RHLO",
    "CAPE1",
    "CAPE3",
    "SHTFL2",
    "TCOND7002",
    "INST2",
    "CP1",
    "TCONDSYM2",
    "COUPLSYM3",
    "HWFI",
    "VMAX_OP_T0",
    "HWRF",
    "VMAX"
  ) |>
  mutate(
    StormID = droplevels(StormID),
    Year = factor(Year, ordered = FALSE),
    Month = factor(Month, ordered = FALSE)
  )
str(StormdataTrain5)

linFit2 <- brm(
  formula = VMAX ~ 
    Year +
    Month +
    StormElapsedTime + 
    #basin + 
    LAT +
    LON +
    MINSLP +
    SHR_MAG +
    STM_SPD +
    SST +
    RHLO +
    CAPE1 +
    CAPE3 +
    SHTFL2 +
    TCOND7002 +
    INST2 +
    CP1 +
    TCONDSYM2 +
    COUPLSYM3 +
    HWFI +
    VMAX_OP_T0 +
    HWRF +
    (1|StormID),
  data = StormdataTrain5, 
  family = gaussian(), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52, 
  warmup = 1000
)
prior_summary(linFit2)
posterior_summary(linFit2)
linFit2

print(linFit2, digits = 4)
plot(linFit2)
loo(linFit2)
performance::check_distribution(linFit2)
performance::check_outliers(linFit2)
performance::check_heteroskedasticity(linFit2)
performance_rmse(linFit2)
performance_mae(linFit2)
variance_decomposition(linFit2)
fixef(linFit2)
ranef(linFit2)

bayes_factor(linFit1, linFit2)


StormdataTrain6 <- StormdataTrain3 |>
  select(
    "StormID",
    #Date,
    Year,
    Month,
    StormElapsedTime,
    "basin",
    "LAT",
    "LON",
    "MINSLP",
    "SHR_MAG",
    "STM_SPD",
    "SST",
    "RHLO",
    "CAPE1",
    "CAPE3",
    "SHTFL2",
    "TCOND7002",
    "INST2",
    "CP1",
    "TCONDSYM2",
    "COUPLSYM3",
    "HWFI",
    "VMAX_OP_T0",
    "HWRF",
    "VMAX"
  ) |>
  mutate(
    StormID = droplevels(StormID),
    Year = factor(Year, ordered = FALSE),
    Month = factor(Month, ordered = FALSE)
  ) |>
  mutate(
    across(where(is.numeric) & !VMAX, function(x){scale(x)})
  )
str(StormdataTrain6)

linFit3 <- brm(
  formula = VMAX ~ 
    Year +
    Month +
    StormElapsedTime + 
    #basin + 
    LAT +
    LON +
    MINSLP +
    SHR_MAG +
    STM_SPD +
    SST +
    RHLO +
    CAPE1 +
    CAPE3 +
    SHTFL2 +
    TCOND7002 +
    INST2 +
    CP1 +
    TCONDSYM2 +
    COUPLSYM3 +
    HWFI +
    VMAX_OP_T0 +
    HWRF +
    (1|StormID),
  data = StormdataTrain6, 
  family = gaussian(), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52, 
  warmup = 1000
)
prior_summary(linFit3)
posterior_summary(linFit3)
linFit3

print(linFit3, digits = 4)
plot(linFit3)
loo(linFit3)
performance::check_distribution(linFit3)
performance::check_outliers(linFit3)
performance::check_heteroskedasticity(linFit3)
performance_rmse(linFit3)
performance_mae(linFit3)
variance_decomposition(linFit3)
fixef(linFit3)
ranef(linFit3)

bayes_factor(linFit2, linFit3)



linFit4 <- brm(
  formula = VMAX ~ 
    Year +
    Month +
    StormElapsedTime + 
    I(StormElapsedTime^2) +
    #basin + 
    LAT +
    LON +
    MINSLP +
    SHR_MAG +
    STM_SPD +
    SST +
    RHLO +
    CAPE1 +
    CAPE3 +
    SHTFL2 +
    TCOND7002 +
    INST2 +
    CP1 +
    TCONDSYM2 +
    COUPLSYM3 +
    HWFI +
    VMAX_OP_T0 +
    HWRF +
    (1|StormID),
  data = StormdataTrain6, 
  family = gaussian(), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52, 
  warmup = 1000
)
prior_summary(linFit4)
posterior_summary(linFit4)
linFit4

print(linFit4, digits = 4)
plot(linFit4)
loo(linFit4)
performance::check_distribution(linFit4)
performance::check_outliers(linFit4)
performance::check_heteroskedasticity(linFit4)
performance_rmse(linFit4)
performance_mae(linFit4)
variance_decomposition(linFit4)
fixef(linFit4)
ranef(linFit4)

bayes_factor(linFit4, linFit3)


linFit5 <- brm(
  formula = VMAX ~ 
    Year +
    Month +
    StormElapsedTime + 
    I(StormElapsedTime^2) +
    basin + 
    LAT +
    LON +
    MINSLP +
    SHR_MAG +
    STM_SPD +
    SST +
    RHLO +
    CAPE1 +
    CAPE3 +
    SHTFL2 +
    TCOND7002 +
    INST2 +
    CP1 +
    TCONDSYM2 +
    COUPLSYM3 +
    HWFI +
    VMAX_OP_T0 +
    HWRF +
    (1|StormID),
  data = StormdataTrain6, 
  family = gaussian(), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52, 
  warmup = 1000
)
prior_summary(linFit5)
posterior_summary(linFit5)
linFit5

print(linFit5, digits = 4)
plot(linFit5)
loo(linFit5)
performance::check_distribution(linFit5)
performance::check_outliers(linFit5)
performance::check_heteroskedasticity(linFit5)
performance_rmse(linFit5)
performance_mae(linFit5)
variance_decomposition(linFit5)
fixef(linFit5)
ranef(linFit5)

bayes_factor(linFit5, linFit4)


linFit6 <- brm(
  formula = VMAX ~ 
    #Year +
    Month +
    StormElapsedTime + 
    I(StormElapsedTime^2) +
    basin + 
    LAT +
    LON +
    MINSLP +
    SHR_MAG +
    STM_SPD +
    SST +
    RHLO +
    CAPE1 +
    CAPE3 +
    SHTFL2 +
    TCOND7002 +
    INST2 +
    CP1 +
    TCONDSYM2 +
    COUPLSYM3 +
    HWFI +
    VMAX_OP_T0 +
    HWRF +
    (1|StormID),
  data = StormdataTrain6, 
  family = gaussian(), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52, 
  warmup = 1000
)
prior_summary(linFit6)
posterior_summary(linFit6)
linFit6

print(linFit6, digits = 4)
plot(linFit6)
loo(linFit6)
performance::check_distribution(linFit6)
performance::check_outliers(linFit6)
performance::check_heteroskedasticity(linFit6)
performance_rmse(linFit6)
performance_mae(linFit6)
variance_decomposition(linFit6)
fixef(linFit6)
ranef(linFit6)

bayes_factor(linFit6, linFit4)
bayes_factor(linFit6, linFit5)


linFit7 <- brm(
  formula = VMAX ~ 
    #Year +
    Month +
    s(StormElapsedTime) + 
    #I(StormElapsedTime^2) +
    basin + 
    t2(LAT, LON) +
    #LON +
    MINSLP +
    SHR_MAG +
    STM_SPD +
    SST +
    RHLO +
    CAPE1 +
    CAPE3 +
    SHTFL2 +
    TCOND7002 +
    INST2 +
    CP1 +
    TCONDSYM2 +
    COUPLSYM3 +
    HWFI +
    VMAX_OP_T0 +
    HWRF +
    (1|StormID),
  data = StormdataTrain6, 
  family = gaussian(), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52, 
  warmup = 1000
)
prior_summary(linFit7)
posterior_summary(linFit7)
linFit7

print(linFit7, digits = 4)
plot(linFit7)
loo(linFit7)
performance::check_distribution(linFit7)
performance::check_outliers(linFit7)
performance::check_heteroskedasticity(linFit7)
performance_rmse(linFit7)
performance_mae(linFit7)
variance_decomposition(linFit7)
fixef(linFit7)
ranef(linFit7)

bayes_factor(linFit7, linFit4)
bayes_factor(linFit7, linFit5)
bayes_factor(linFit7, linFit6)


linFit8 <- brm(
  formula = log(VMAX) ~ 
    #Year +
    Month +
    s(StormElapsedTime) + 
    #I(StormElapsedTime^2) +
    basin + 
    t2(LAT, LON) +
    #LON +
    MINSLP +
    SHR_MAG +
    STM_SPD +
    SST +
    RHLO +
    CAPE1 +
    CAPE3 +
    SHTFL2 +
    TCOND7002 +
    INST2 +
    CP1 +
    TCONDSYM2 +
    COUPLSYM3 +
    HWFI +
    VMAX_OP_T0 +
    HWRF +
    (1|StormID),
  data = StormdataTrain6, 
  family = gaussian(), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52, 
  warmup = 1000
)
prior_summary(linFit8)
posterior_summary(linFit8)
linFit8

print(linFit8, digits = 4)
plot(linFit8)
loo(linFit8)
performance::check_distribution(linFit8)
performance::check_outliers(linFit8)
performance::check_heteroskedasticity(linFit8)
performance_rmse(linFit8)
performance_mae(linFit8)
variance_decomposition(linFit8)
fixef(linFit8)
ranef(linFit8)

bayes_factor(linFit8, linFit4)
bayes_factor(linFit8, linFit5)
bayes_factor(linFit8, linFit6)
bayes_factor(linFit8, linFit7)



# L <- 10 # number of knots
# B1   <- bs(modelData$LON, df=2*L, intercept=TRUE) # Longitude basis functions
# B2   <- bs(modelData$LAT, df=L, intercept=TRUE)   # Latitude basis functions
# X    <- NULL
# for(j in 1:ncol(B1)){
#   for(k in 1:ncol(B2)){
#     X <- cbind(X,B1[,j]*B2[,k])  # Products
#   }
# }
# X    <- X[,apply(X,2,max)>0.1]  # Remove basis function that are near zero for all sites
# X    <- ifelse(X>0.001,X,0)
# p    <- ncol(X)
# 
# plot(B1[,1])
# 
# index <- sample(1:p, 1)
# spline_data <- modelData |>
#   select(
#     LON, 
#     LAT
#   ) |>
#   add_column(spline = X[,index])
# #world_coordinates <- map_data("world") 
# ggplot() + 
#   # geom_map() function takes world coordinates  
#   # as input to plot world map 
#   geom_map( 
#     data = world_coordinates, map = world_coordinates, 
#     aes(x = long, y = lat, map_id = region) 
#   ) + 
#   geom_point(
#     data = spline_data,
#     aes(x = LON-360, y = LAT, 
#         color = spline)
#   ) +
#   xlim(c(-180,0)) +
#   ylim(c(0,60)) +
#   scale_color_continuous(low = "green", high = "red") +
#   theme_bw()
# 
# gamMod1 <- gam::gam(VMAX ~ gam::s(LON, df = 3), data = modelData)
# summary(gamMod1)
# 
# east <- north <- 1:10
# Grid <- expand.grid(east, north)
# K <- nrow(Grid)
# 
# # set up distance and neighbourhood matrices
# distance <- as.matrix(dist(Grid))
# W <- array(0, c(K, K))
# W[distance == 1] <- 1
# 
# # generate the covariates and response data
# x1 <- rnorm(K)
# x2 <- rnorm(K)
# theta <- rnorm(K, sd = 0.05)
# phi <- brms::rmulti_normal(
#   1, mu = rep(0, K), Sigma = 0.4 * exp(-0.1 * distance)
# )
# eta <- x1 + x2 + phi
# prob <- exp(eta) / (1 + exp(eta))
# size <- rep(50, K)
# y <- rbinom(n = K, size = size, prob = prob)
# dat <- data.frame(y, size, x1, x2)
# 
# # fit a CAR model
# fit <- brm(y | trials(size) ~ x1 + x2 + car(W),
#            data = dat, data2 = list(W = W),
#            family = binomial())
# summary(fit)
# 
# spaceTimefit <- brm(
#   formula = 
# )


# 5. Prediction ====
completeStormdata <- Stormdata_raw |>
  mutate(
    Observed = ifelse(is.na(VMAX), "Predict", "Observed")
  )

Actual_Y <- Stormdata_raw |> 
  select(-VMAX) |>
  cbind(Actual_Y)

for(i in 1:nrow(completeStormdata)){
  if(is.na(completeStormdata$VMAX[i])){
    completeStormdata$VMAX[i] <- Actual_Y$VMAX[i]
  }
}

completeStormdata |>
  filter(Observed == "Observed") |>
  summarise(
    mean(abs(HWRF - VMAX))
  )






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
    basin,
    lead_time,
    LAT,
    LON,
    MINSLP,
    SHR_MAG,
    STM_SPD,
    SST,
    RHLO,
    CAPE1,
    CAPE3,
    SHTFL2,
    TCOND7002,
    INST2,
    CP1,
    TCONDSYM2,
    HWFI,
    VMAX_OP_T0,
    VMAX
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

