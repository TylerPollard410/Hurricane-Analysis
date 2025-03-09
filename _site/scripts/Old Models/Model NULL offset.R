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
library(GGally)
library(patchwork)
library(fitdistrplus)
library(caret)
library(splines)
library(mgcv)
library(DescTools)
library(car)
library(bayesplot)
library(BayesFactor)
library(rstanarm)
library(tidybayes)
library(loo)
library(brms)
library(performance)
library(tidyverse)

# Read in data ----
Stormdata_raw <- fread("_data/E2_data.csv")
Stormdata <- Stormdata_raw |>
  mutate(
    Obs = 1:nrow(Stormdata_raw),
    StormID = factor(StormID),
    basin = factor(basin),
    Date = as_datetime(Date, tz = "UTC")
  )


## Clean data ----
### Training ----
StormdataTrain <- Stormdata |> filter(complete.cases(VMAX))

# Remove not varying 
StormdataTrain2 <- StormdataTrain |>
  select(-lead_time)

# Create Date vars 
dataTrainYears <- year(StormdataTrain2$Date)
dataTrainMonths <- month(StormdataTrain2$Date, label = TRUE)
dataTrainDays <- day(StormdataTrain2$Date)
dataTrainYearDay <- yday(StormdataTrain2$Date)

#### Train 3 ----
StormdataTrain3 <- StormdataTrain2 |>
  mutate(
    Year = factor(dataTrainYears, ordered = TRUE),
    Month = dataTrainMonths,
    Day = dataTrainYearDay
  ) |>
  group_by(StormID) |>
  mutate(
    StormElapsedTime = as.numeric(difftime(Date, min(Date), units = "hours")),
    StormElapsedTime2 = StormElapsedTime/6
  ) |>
  select(
    StormID,
    Date,
    Year,
    Month,
    Day,
    StormElapsedTime,
    StormElapsedTime2,
    everything()
  ) |>
  ungroup()

#### Train 7 ----
StormdataTrain7 <- StormdataTrain3 |>
  select(
    "StormID",
    #Date,
    Year,
    Month,
    Day,
    StormElapsedTime,
    StormElapsedTime2,
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
    across(where(is.numeric) & !c(VMAX, 
                                  Day,
                                  StormElapsedTime,StormElapsedTime2, 
                                  LAT, LON),
           function(x){scale(x)})
  )
str(StormdataTrain7)

#### Train 7scale ----
StormdataTrain7scale <- StormdataTrain3 |>
  select(
    "StormID",
    #Date,
    Year,
    Month,
    Day,
    StormElapsedTime,
    StormElapsedTime2,
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
    across(where(is.numeric) & !c(VMAX, StormElapsedTime2),
           function(x){scale(x)})
  )
str(StormdataTrain7scale)

#### Train 8 ----
StormdataTrain8 <- StormdataTrain3 |>
  select(
    "StormID",
    #Date,
    Year,
    Month,
    Day,
    StormElapsedTime,
    StormElapsedTime2,
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
    Month = factor(Month, ordered = FALSE),
  ) |>
  mutate(
    across(where(is.numeric) & !c(VMAX, HWRF, StormElapsedTime2),
           function(x){scale(x)})
  ) |>
  mutate(
    propVMAX = VMAX/HWRF
  )
str(StormdataTrain8)

### Test ----
StormdataTest <- Stormdata |> filter(!complete.cases(VMAX))

# Remove not varying 
StormdataTest2 <- StormdataTest |>
  select(-lead_time)

# Create Date vars 
dataTestYears <- year(StormdataTest2$Date)
dataTestMonths <- month(StormdataTest2$Date, label = TRUE)
dataTestDays <- day(StormdataTest2$Date)
dataTestYearDay <- yday(StormdataTest2$Date)

#### Test 3 ----
StormdataTest3 <- StormdataTest2 |>
  mutate(
    Year = factor(dataTestYears, ordered = TRUE),
    Month = dataTestMonths,
    Day = dataTestYearDay
  ) |>
  group_by(StormID) |>
  mutate(
    StormElapsedTime = as.numeric(difftime(Date, min(Date), units = "hours")),
    StormElapsedTime2 = StormElapsedTime/6
  ) |>
  select(
    StormID,
    Date,
    Year,
    Month,
    Day,
    StormElapsedTime,
    StormElapsedTime2,
    everything()
  ) |>
  ungroup()

#### Test Final ----
StormdataTest7 <- StormdataTest3 |>
  select(
    "StormID",
    #Date,
    Year,
    Month,
    Day,
    StormElapsedTime,
    StormElapsedTime2,
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
    across(where(is.numeric) & !c(VMAX, Day,
                                  StormElapsedTime,StormElapsedTime2, 
                                  LAT, LON),
           function(x){scale(x,
                             center = attr(StormdataTrain7 |> pull(x), "scaled:center"),
                             scale = attr(StormdataTrain7 |> pull(x), "scaled:scale"))
           })
  )
str(StormdataTest7)

#### Test Final Scale----
StormdataTest7scale <- StormdataTest3 |>
  select(
    "StormID",
    #Date,
    Year,
    Month,
    Day,
    StormElapsedTime,
    StormElapsedTime2,
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
    across(where(is.numeric) & !c(VMAX, StormElapsedTime2),
           function(x){scale(x,
                             center = attr(StormdataTrain7scale |> pull(x), "scaled:center"),
                             scale = attr(StormdataTrain7scale |> pull(x), "scaled:scale"))
           })
  )
str(StormdataTest7scale)

#### Test Final Scale2 ----
StormdataTest8 <- StormdataTest3 |>
  select(
    "StormID",
    #Date,
    Year,
    Month,
    Day,
    StormElapsedTime,
    StormElapsedTime2,
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
    across(where(is.numeric) & !c(VMAX, HWRF, StormElapsedTime2),
           function(x){scale(x,
                             center = attr(StormdataTrain8 |> pull(x), "scaled:center"),
                             scale = attr(StormdataTrain8 |> pull(x), "scaled:scale"))
           })
  ) |>
  mutate(
    propVMAX = VMAX/HWRF
  )
str(StormdataTest8)

### Actual ----
Actual_Y <- fread("_data/Actual Y.csv")
Actual_Yvec <- Actual_Y |> filter(complete.cases(x)) |> pull(x)

#### Total Data ----
StormdataTest3final <- bind_cols(
  StormdataTest3,
  Actual = Actual_Yvec
) |>
  mutate(VMAX = Actual) |>
  select(-Actual)

StormComplete <- bind_rows(
  StormdataTrain3,
  StormdataTest3final
) |> arrange(Obs)

HWRF_MAE <- mean(abs(StormComplete$VMAX - StormComplete$HWRF))

# Plot VMAX ----
ggpairs(StormdataTrain8 |>
          select(-StormID, -Date, -Month, -Year))

## Scatter ----
ggplot(data = StormdataTrain3, aes(x = HWRF, y = VMAX)) +
  geom_point() +
  geom_smooth()

mu1 <- mean(StormdataTrain3$VMAX)
mu2 <- mean(StormdataTrain3$HWRF)
sd1 <- sd(StormdataTrain3$VMAX)
sd2 <- sd(StormdataTrain3$HWRF)

mean(StormdataTrain3$VMAX - StormdataTrain3$HWRF)
sd(StormdataTrain3$VMAX - StormdataTrain3$HWRF)

mean(StormdataTrain3$VMAX/StormdataTrain3$HWRF)
sd(StormdataTrain3$VMAX/StormdataTrain3$HWRF)

mean(log(StormdataTrain3$VMAX/StormdataTrain3$HWRF))
sd(log(StormdataTrain3$VMAX/StormdataTrain3$HWRF))

## Histogram ----
ggplot(geom_histogram(
  aes(rlnorm(1705, meanlog = 2.786386, sdlog = 13.69776))
))
plot(hist(rlnorm(1705, meanlog = log(mu1), sdlog = log(sd1))))

ggplot(data = StormdataTrain3) +
  # geom_histogram(
  #   aes(x = VMAX, after_stat(density)),
  #   color = "#99c7c7", fill = "#bcdcdc",
  #   bins = 100) +
  geom_density(#data = final_data3,
    aes(x = VMAX),
    color = "#007C7C", 
    linewidth = 1) +
  geom_density(#data = final_data3,
    aes(x = HWRF),
    color = "blue", 
    linewidth = 1) +
  geom_density(#data = final_data3,
    aes(x = HWRF-1),
    color = "blue", 
    linewidth = 1) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  scale_x_continuous(limits = c(-25,200)) +
  labs(title = "Density Plot",
       subtitle = "Data",
       x = "VMAX",
       y = "Density") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

ggplot(data = StormdataTrain3) +
  geom_histogram(
    aes(x = log(VMAX), after_stat(density)),
    color = "#99c7c7", fill = "#bcdcdc",
    bins = 100) +
  geom_density(#data = final_data3,
    aes(x = log(VMAX)),
    color = "#007C7C", 
    linewidth = 1) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  #scale_x_continuous(limits = c(-0.5,1)) +
  labs(title = "Density Plot",
       subtitle = "Data",
       x = "VMAX/HWRF",
       y = "Density") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

ggplot(data = StormdataTrain8) +
  geom_histogram(
    aes(x = propVMAX, after_stat(density)),
    color = "#99c7c7", fill = "#bcdcdc",
    bins = 100) +
  geom_density(#data = final_data3,
    aes(x = propVMAX),
    color = "#007C7C", 
    linewidth = 1) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  #scale_x_continuous(limits = c(-0.5,1)) +
  labs(title = "Density Plot",
       subtitle = "Data",
       x = "VMAX/HWRF",
       y = "Density") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

ggplot(data = StormdataTrain8) +
  geom_histogram(
    aes(x = log(propVMAX), after_stat(density)),
    color = "#99c7c7", fill = "#bcdcdc",
    bins = 100) +
  geom_density(#data = final_data3,
    aes(x = log(propVMAX)),
    color = "#007C7C", 
    linewidth = 1) +
  # geom_density(#data = final_data3,
  #   aes(x = log(VMAX) - log(HWRF)),
  #   color = "red", 
  #   linewidth = 1) +
  geom_density(#data = final_data3,
    aes(x = rnorm(1705, 1.049843 - 1, 0.239613/sqrt(2))),
    color = "blue", 
    linewidth = 1) +
  # geom_density(#data = final_data3,
  #   aes(x = exp(rnorm(1705, 1.049843, 0.239613))),
  #   color = "green", 
  #   linewidth = 1) +
  geom_density(#data = final_data3,
    aes(x = log(rlnorm(1705, 1.049843, 0.239613/sqrt(2)))+1),
    color = "gold", 
    linewidth = 1) +
  geom_density(#data = final_data3,
    aes(x = log(rlnorm(1705, exp(1.049843-1), 0.239613))),
    color = "red", 
    linewidth = 1) +
  geom_density(#data = final_data3,
    aes(x = log(rlnorm(1705, 1.049843-0.239613^2/2, 0.239613))),
    color = "purple", 
    linewidth = 1) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  #scale_x_continuous(limits = c(-0.5,1)) +
  labs(title = "Density Plot",
       subtitle = "Data",
       x = "log(VMAX/HWRF)",
       y = "Density") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )



## Time ----
ggplot(data = StormdataTrain3) +
  geom_point(aes(x = Date, y = VMAX)) +
  scale_x_datetime(date_breaks = "month", 
                   date_minor_breaks = "day", 
                   date_labels = "%b-%Y") + 
  theme(
    axis.text.x = element_text(angle = 90)
  )

ggplot(data = StormdataTrain3) +
  geom_line(aes(x = StormElapsedTime, y = VMAX)) +
  scale_x_continuous(breaks = seq(0,402,24)) +
  facet_wrap(vars(StormID))

ggplot(data = StormdataTrain3) +
  geom_line(aes(x = HWRF, y = VMAX)) +
  #scale_x_continuous(breaks = seq(0,402,24)) +
  facet_wrap(vars(StormID))


## Map ----
# Do by time next
# [1] "StormID"           "Date"              "Year"              "Month"             "Day"              
# [6] "StormElapsedTime"  "StormElapsedTime2" "basin"             "LAT"               "LON"              
# [11] "MINSLP"            "SHR_MAG"           "STM_SPD"           "SST"               "RHLO"             
# [16] "CAPE1"             "CAPE3"             "SHTFL2"            "TCOND7002"         "INST2"            
# [21] "CP1"               "TCONDSYM2"         "COUPLSYM3"         "HWFI"              "VMAX_OP_T0"       
# [26] "HWRF"              "VMAX"

range(StormdataTrain3$StormElapsedTime)
range(StormdataTrain3$LAT)
range(StormdataTrain3$LON)
range(StormdataTest7$StormElapsedTime)
range(StormdataTest7$LAT)
which.min(StormdataTest7$LAT)
range(StormdataTest7$LON)
which.min(StormdataTest7$LON)

world_coordinates <- map_data("world") 
ggplot() + 
  # geom_map() function takes world coordinates  
  # as input to plot world map 
  geom_map( 
    data = world_coordinates, map = world_coordinates, 
    aes(map_id = region) 
  ) + 
  geom_point(
    data = StormdataTest7,
    aes(x = LON-360, y = LAT, 
        color = StormElapsedTime)
  ) +
  xlim(c(-190,0)) +
  ylim(c(0,60)) +
  scale_color_continuous(low = "green", high = "red") +
  theme_bw()

# Find Distributions ----

ddist <- descdist(StormdataTrain3$VMAX, discrete = FALSE)
summary(ddist)

ddistLog <- descdist(log(StormdataTrain3$VMAX), discrete = FALSE)
summary(ddist)

ddistLogProp <- descdist(log(StormdataTrain3$VMAX/StormdataTrain3$HWRF), discrete = FALSE)
summary(ddist)

fitnorm1 <- fitdist(StormdataTrain3$VMAX, distr = "norm", method = "mle")
summary(fitnorm1)
fitdist(StormdataTrain3$VMAX, distr = "norm", method = "mme")
fitdist(StormdataTrain3$VMAX, distr = "norm", method = "mle")
fitdist(StormdataTrain3$VMAX, distr = "norm", method = "mle")

# Fit model ----
load(file = "_data/gammaFit1.RData")
load(file = "_data/logNormalFit1.RData")
load(file = "_data/logNormalFit2A.RData")
load(file = "_data/propLinFit1.RData")
load(file = "_data/propStudentFit1.RData")
load(file = "_data/propStudentFit2.RData")
load(file = "_data/propStudentFit5.RData")
load(file = "_data/propGammaFit1.RData")
load(file = "_data/logpropLinFit1.RData")
load(file = "_data/logpropStudentFit1.RData")
load(file = "_data/logpropStudentFit2.RData")

## NULL MODELS ----
### GAUSSIAN ----
#### Model 1 ----
LinFitNULL <- brm(
  bf(
    # VMAX ~ offset(HWRF), family = gaussian()
    # VMAX ~ HWRF, family = gaussian(link = "log")
    VMAX ~ offset(HWRF), family = gaussian()
    # VMAX ~ log(HWRF), family = gaussian(link = "log)
  ),
  data = StormdataTrain8, 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)

LinFitNULL2 <- brm(
  bf(
    #VMAX ~ offset(HWRF), family = gaussian()
    # VMAX ~ offset(HWRF), family = gaussian(link = "log")
    # VMAX ~ log(HWRF), family = gaussian()
    VMAX ~ offset(log(HWRF)), family = gaussian(link = "log")
  ),
  data = StormdataTrain8, 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)

#save(LinFitNULL, file = "_data/LinFitNULL.RData")
prior_summary(LinFitNULL)
round(posterior_summary(LinFitNULL, probs = c(0.025, 0.975)))
LinFitNULL

print(LinFitNULL, digits = 4)
plot(LinFitNULL)
LinFitNULLppcFit <- pp_check(LinFitNULL, ndraws = 100) + 
  labs(title = "LinFitNULL Fit PPC") +
  theme_bw()
LinFitNULLppcFit
LinFitNULLloo <- loo(LinFitNULL)
waic(LinFitNULL)
performance::check_distribution(LinFitNULL)
performance::check_outliers(LinFitNULL)
performance::check_heteroskedasticity(LinFitNULL)
performance_rmse(LinFitNULL)
performance_mae(LinFitNULL)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(LinFitNULL)

print(LinFitNULL2, digits = 4)
plot(LinFitNULL2)
LinFitNULL2ppcFit <- pp_check(LinFitNULL2, ndraws = 100) + 
  labs(title = "LinFitNULL2 Fit PPC") +
  theme_bw()
LinFitNULL2ppcFit
LinFitNULL2loo <- loo(LinFitNULL2)
waic(LinFitNULL2)
performance::check_distribution(LinFitNULL2)
performance::check_outliers(LinFitNULL2)
performance::check_heteroskedasticity(LinFitNULL2)
performance_rmse(LinFitNULL2)
performance_mae(LinFitNULL2)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(LinFitNULL2)

LinFitNULLppcFit/LinFitNULL2ppcFit


variance_decomposition(LinFitNULL)
exp(fixef(LinFitNULL))
ranef(LinFitNULL)

bayes_R2(LinFitNULL)

bayes_factor(LinFitNULL2, LinFitNULL)
loo_compare(LinFitNULLloo,LinFitNULL2loo)

LinFitNULLsmooths <- conditional_smooths(LinFitNULL)
plot(LinFitNULLsmooths, stype = "raster", ask = FALSE)
LinFitNULLeffects <- conditional_effects(LinFitNULL, 
                                               method = "posterior_predict",
                                               robust = FALSE,
                                               re_formula = NULL)
LinFitNULLeffects <- conditional_effects(LinFitNULL)
plot(LinFitNULLeffects, points = TRUE, ask = FALSE)

##### Prediction ----
## Fitted
LinFitNULLfinalFit <- posterior_predict(LinFitNULL)
LinFitNULLfinalFitMean <- colMeans(LinFitNULLfinalFit)
LinFitNULLfinalFitMed <- apply(LinFitNULLfinalFit, 2, function(x){quantile(x, 0.5)})
LinFitNULLfinalFitLCB <- apply(LinFitNULLfinalFit, 2, function(x){quantile(x, 0.025)})
LinFitNULLfinalFitUCB <- apply(LinFitNULLfinalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
LinFitNULLfinalPreds <- posterior_predict(LinFitNULL, 
                                                newdata = StormdataTest8,
                                                allow_new_levels = TRUE)
LinFitNULLfinalPredsMean <- colMeans(LinFitNULLfinalPreds)
LinFitNULLfinalPredsMed <- apply(LinFitNULLfinalPreds, 2, function(x){quantile(x, 0.5)})
LinFitNULLfinalPredsLCB <- apply(LinFitNULLfinalPreds, 2, function(x){quantile(x, 0.025)})
LinFitNULLfinalPredsUCB <- apply(LinFitNULLfinalPreds, 2, function(x){quantile(x, 0.975)})

LinFitNULLpredMetrics <- tibble(
  Fit = "LinFitNULL",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(LinFitNULLfinalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(LinFitNULLfinalFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < LinFitNULLfinalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(LinFitNULLfinalPredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(LinFitNULLfinalPredsMed - Actual_Yvec)),
  COV_pred = mean(LinFitNULLfinalPredsLCB < Actual_Yvec & Actual_Yvec < LinFitNULLfinalPredsUCB)
)
LinFitNULLpredMetrics

## Fitted
LinFitNULL2finalFit <- posterior_predict(LinFitNULL2)
LinFitNULL2finalFitMean <- colMeans(LinFitNULL2finalFit)
LinFitNULL2finalFitMed <- apply(LinFitNULL2finalFit, 2, function(x){quantile(x, 0.5)})
LinFitNULL2finalFitLCB <- apply(LinFitNULL2finalFit, 2, function(x){quantile(x, 0.025)})
LinFitNULL2finalFitUCB <- apply(LinFitNULL2finalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
LinFitNULL2finalPreds <- posterior_predict(LinFitNULL2, 
                                          newdata = StormdataTest8,
                                          allow_new_levels = TRUE)
LinFitNULL2finalPredsMean <- colMeans(LinFitNULL2finalPreds)
LinFitNULL2finalPredsMed <- apply(LinFitNULL2finalPreds, 2, function(x){quantile(x, 0.5)})
LinFitNULL2finalPredsLCB <- apply(LinFitNULL2finalPreds, 2, function(x){quantile(x, 0.025)})
LinFitNULL2finalPredsUCB <- apply(LinFitNULL2finalPreds, 2, function(x){quantile(x, 0.975)})

LinFitNULL2predMetrics <- tibble(
  Fit = "LinFitNULL2",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(LinFitNULL2finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(LinFitNULL2finalFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < LinFitNULL2finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(LinFitNULL2finalPredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(LinFitNULL2finalPredsMed - Actual_Yvec)),
  COV_pred = mean(LinFitNULL2finalPredsLCB < Actual_Yvec & Actual_Yvec < LinFitNULL2finalPredsUCB)
)
LinFitNULLpredMetrics
LinFitNULL2predMetrics

##### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = LinFitNULLfinalPreds) +
  labs(title = "LinFitNULL Predict") +
  theme_bw()

LinFitNULLFitDF <- bind_cols(
  StormdataTrain3,
  LCB = LinFitNULLfinalFitLCB,
  Mean = LinFitNULLfinalFitMean,
  Med = LinFitNULLfinalFitMed,
  UCB = LinFitNULLfinalFitUCB
)

LinFitNULLstormsFitplot <- ggplot(data = LinFitNULLFitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))
LinFitNULLstormsFitplot

## Prediction
LinFitNULLPredDF <- bind_cols(
  StormdataTest3,
  LCB = LinFitNULLfinalPredsLCB,
  Mean = LinFitNULLfinalPredsMean,
  Med = LinFitNULLfinalPredsMed,
  UCB = LinFitNULLfinalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

LinFitNULLstormsPredplot <- ggplot(data = LinFitNULLPredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX, color = "Observed")) +
  geom_line(aes(y = Mean, color = "PPD Mean")) +
  facet_wrap(vars(StormID), ncol = 6)+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50)) +
  labs(title = "LinFitNULL PPD Mean vs Observed VMAX",
       subtitle = "95% Credible Interval about PPD Mean") +
  scale_color_manual(name = NULL, values = c("black","red")) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1))
  ) +
  theme_bw()
LinFitNULLstormsPredplot

##### PPC ----
###### Quantile 2.5 
LinFitNULLLCBsims <- apply(LinFitNULLfinalFit, 
                                 MARGIN = 1,
                                 function(x){
                                   quantile(x, 0.025)
                                 })
LinFitNULLLCBpvalueVec <- LinFitNULLLCBsims < quantile(StormdataTrain3$VMAX, 0.025)
LinFitNULLLCBpvalue <- sum(LinFitNULLLCBpvalueVec)
LinFitNULLLCBpvalue <- round(LinFitNULLLCBpvalue/4000, 3)
LinFitNULLLCBpvalue <- min(LinFitNULLLCBpvalue, 1 - LinFitNULLLCBpvalue)

LinFitNULL_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           LinFitNULLfinalFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", LinFitNULLLCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#LinFitNULL_ppcLCB

###### Quantile 97.5 
LinFitNULLUCBsims <- apply(LinFitNULLfinalFit, 
                                 MARGIN = 1,
                                 function(x){
                                   quantile(x, 0.975)
                                 })
LinFitNULLUCBpvalueVec <- LinFitNULLUCBsims < quantile(StormdataTrain3$VMAX, 0.975)
LinFitNULLUCBpvalue <- as.numeric(sum(LinFitNULLUCBpvalueVec))
LinFitNULLUCBpvalue <- round(LinFitNULLUCBpvalue/4000, 3)
LinFitNULLUCBpvalue <- min(LinFitNULLUCBpvalue, 1 - LinFitNULLUCBpvalue)

LinFitNULL_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           LinFitNULLfinalFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", LinFitNULLUCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#LinFitNULL_ppcUCB

###### Mean 
LinFitNULLMEANsims <- apply(LinFitNULLfinalFit, 
                                  MARGIN = 1,
                                  function(x){
                                    mean(x)
                                  })
LinFitNULLMEANpvalueVec <- LinFitNULLMEANsims < mean(StormdataTrain3$VMAX)
LinFitNULLMEANpvalue <- sum(LinFitNULLMEANpvalueVec)
LinFitNULLMEANpvalue <- round(LinFitNULLMEANpvalue/4000, 3)
LinFitNULLMEANpvalue <- min(LinFitNULLMEANpvalue, 1 - LinFitNULLMEANpvalue)

LinFitNULL_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           LinFitNULLfinalFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", LinFitNULLMEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#LinFitNULL_ppcMEAN

###### Med 
LinFitNULLMEDsims <- apply(LinFitNULLfinalFit, 
                                 MARGIN = 1,
                                 function(x){
                                   quantile(x, 0.5)
                                 })
LinFitNULLMEDpvalueVec <- LinFitNULLMEDsims < quantile(StormdataTrain3$VMAX, 0.5)
LinFitNULLMEDpvalue <- sum(LinFitNULLMEDpvalueVec)
LinFitNULLMEDpvalue <- round(LinFitNULLMEDpvalue/4000, 3)
LinFitNULLMEDpvalue <- min(LinFitNULLMEDpvalue, 1 - LinFitNULLMEDpvalue)

LinFitNULL_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           LinFitNULLfinalFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", LinFitNULLMEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#LinFitNULL_ppcMED

###### SD 
LinFitNULLSDsims <- apply(LinFitNULLfinalFit, 
                                MARGIN = 1,
                                function(x){
                                  sd(x)
                                })
LinFitNULLSDpvalueVec <- LinFitNULLSDsims < sd(StormdataTrain3$VMAX)
LinFitNULLSDpvalue <- sum(LinFitNULLSDpvalueVec)
LinFitNULLSDpvalue <- round(LinFitNULLSDpvalue/4000, 3)
LinFitNULLSDpvalue <- min(LinFitNULLSDpvalue, 1 - LinFitNULLSDpvalue)

LinFitNULL_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           LinFitNULLfinalFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", LinFitNULLSDpvalue, ")")) +
  theme_bw() +
  legend_none()
#LinFitNULL_ppcSD

###### Range 
LinFitNULLRANGEsims <- apply(LinFitNULLfinalFit, 
                                   MARGIN = 1,
                                   function(x){
                                     max(x)-min(x)
                                   })
LinFitNULLRANGEpvalueVec <- LinFitNULLRANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
LinFitNULLRANGEpvalue <- sum(LinFitNULLRANGEpvalueVec)
LinFitNULLRANGEpvalue <- round(LinFitNULLRANGEpvalue/4000, 3)
LinFitNULLRANGEpvalue <- min(LinFitNULLRANGEpvalue, 1 - LinFitNULLRANGEpvalue)

LinFitNULL_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           LinFitNULLfinalFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", LinFitNULLRANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#LinFitNULL_ppcRANGE

##### Combined Plot ----
LinFitNULL_ppcComb <- 
  #LinFitNULLppcFit /
  LinFitNULL_ppcLCB | LinFitNULL_ppcMED | LinFitNULL_ppcUCB |
  LinFitNULL_ppcRANGE | LinFitNULL_ppcMEAN | LinFitNULL_ppcSD
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
LinFitNULL_ppcComb

##### Bayes p-values ----
LinFitNULLpvalues <- tibble(
  Fit = "LinFitNULL",
  LCB = LinFitNULLLCBpvalue,
  Median = LinFitNULLMEDpvalue,
  UCB = LinFitNULLUCBpvalue,
  Range = LinFitNULLRANGEpvalue,
  Mean = LinFitNULLMEANpvalue,
  SD = LinFitNULLSDpvalue
)
LinFitNULLpvalues

##### PPC ----
###### Quantile 2.5 
LinFitNULL2LCBsims <- apply(LinFitNULL2finalFit, 
                           MARGIN = 1,
                           function(x){
                             quantile(x, 0.025)
                           })
LinFitNULL2LCBpvalueVec <- LinFitNULL2LCBsims < quantile(StormdataTrain3$VMAX, 0.025)
LinFitNULL2LCBpvalue <- sum(LinFitNULL2LCBpvalueVec)
LinFitNULL2LCBpvalue <- round(LinFitNULL2LCBpvalue/4000, 3)
LinFitNULL2LCBpvalue <- min(LinFitNULL2LCBpvalue, 1 - LinFitNULL2LCBpvalue)

LinFitNULL2_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           LinFitNULL2finalFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", LinFitNULL2LCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#LinFitNULL2_ppcLCB

###### Quantile 97.5 
LinFitNULL2UCBsims <- apply(LinFitNULL2finalFit, 
                           MARGIN = 1,
                           function(x){
                             quantile(x, 0.975)
                           })
LinFitNULL2UCBpvalueVec <- LinFitNULL2UCBsims < quantile(StormdataTrain3$VMAX, 0.975)
LinFitNULL2UCBpvalue <- as.numeric(sum(LinFitNULL2UCBpvalueVec))
LinFitNULL2UCBpvalue <- round(LinFitNULL2UCBpvalue/4000, 3)
LinFitNULL2UCBpvalue <- min(LinFitNULL2UCBpvalue, 1 - LinFitNULL2UCBpvalue)

LinFitNULL2_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           LinFitNULL2finalFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", LinFitNULL2UCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#LinFitNULL2_ppcUCB

###### Mean 
LinFitNULL2MEANsims <- apply(LinFitNULL2finalFit, 
                            MARGIN = 1,
                            function(x){
                              mean(x)
                            })
LinFitNULL2MEANpvalueVec <- LinFitNULL2MEANsims < mean(StormdataTrain3$VMAX)
LinFitNULL2MEANpvalue <- sum(LinFitNULL2MEANpvalueVec)
LinFitNULL2MEANpvalue <- round(LinFitNULL2MEANpvalue/4000, 3)
LinFitNULL2MEANpvalue <- min(LinFitNULL2MEANpvalue, 1 - LinFitNULL2MEANpvalue)

LinFitNULL2_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           LinFitNULL2finalFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", LinFitNULL2MEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#LinFitNULL2_ppcMEAN

###### Med 
LinFitNULL2MEDsims <- apply(LinFitNULL2finalFit, 
                           MARGIN = 1,
                           function(x){
                             quantile(x, 0.5)
                           })
LinFitNULL2MEDpvalueVec <- LinFitNULL2MEDsims < quantile(StormdataTrain3$VMAX, 0.5)
LinFitNULL2MEDpvalue <- sum(LinFitNULL2MEDpvalueVec)
LinFitNULL2MEDpvalue <- round(LinFitNULL2MEDpvalue/4000, 3)
LinFitNULL2MEDpvalue <- min(LinFitNULL2MEDpvalue, 1 - LinFitNULL2MEDpvalue)

LinFitNULL2_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           LinFitNULL2finalFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", LinFitNULL2MEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#LinFitNULL2_ppcMED

###### SD 
LinFitNULL2SDsims <- apply(LinFitNULL2finalFit, 
                          MARGIN = 1,
                          function(x){
                            sd(x)
                          })
LinFitNULL2SDpvalueVec <- LinFitNULL2SDsims < sd(StormdataTrain3$VMAX)
LinFitNULL2SDpvalue <- sum(LinFitNULL2SDpvalueVec)
LinFitNULL2SDpvalue <- round(LinFitNULL2SDpvalue/4000, 3)
LinFitNULL2SDpvalue <- min(LinFitNULL2SDpvalue, 1 - LinFitNULL2SDpvalue)

LinFitNULL2_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           LinFitNULL2finalFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", LinFitNULL2SDpvalue, ")")) +
  theme_bw() +
  legend_none()
#LinFitNULL2_ppcSD

###### Range 
LinFitNULL2RANGEsims <- apply(LinFitNULL2finalFit, 
                             MARGIN = 1,
                             function(x){
                               max(x)-min(x)
                             })
LinFitNULL2RANGEpvalueVec <- LinFitNULL2RANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
LinFitNULL2RANGEpvalue <- sum(LinFitNULL2RANGEpvalueVec)
LinFitNULL2RANGEpvalue <- round(LinFitNULL2RANGEpvalue/4000, 3)
LinFitNULL2RANGEpvalue <- min(LinFitNULL2RANGEpvalue, 1 - LinFitNULL2RANGEpvalue)

LinFitNULL2_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           LinFitNULL2finalFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", LinFitNULL2RANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#LinFitNULL2_ppcRANGE

##### Combined Plot ----
LinFitNULL2_ppcComb <- 
  #LinFitNULL2ppcFit /
  LinFitNULL2_ppcLCB | LinFitNULL2_ppcMED | LinFitNULL2_ppcUCB |
  LinFitNULL2_ppcRANGE | LinFitNULL2_ppcMEAN | LinFitNULL2_ppcSD
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
LinFitNULL2_ppcComb

##### Bayes p-values ----
LinFitNULL2pvalues <- tibble(
  Fit = "LinFitNULL2",
  LCB = LinFitNULL2LCBpvalue,
  Median = LinFitNULL2MEDpvalue,
  UCB = LinFitNULL2UCBpvalue,
  Range = LinFitNULL2RANGEpvalue,
  Mean = LinFitNULL2MEANpvalue,
  SD = LinFitNULL2SDpvalue
)
LinFitNULL2pvalues

LinFitNULLfiles <- ls()[str_detect(ls(), pattern = "LinFitNULL")]
LinFitNULLfilesRM <- LinFitNULLfiles[!(LinFitNULLfiles %in% c("LinFitNULL_ppcComb",
                                                                                "LinFitNULLloo",
                                                                                "LinFitNULLpredMetrics",
                                                                                "LinFitNULLpvalues",
                                                                                "LinFitNULLstormsPredplot"))]

##### CV ----
LinFitNULLkfold <- kfold(LinFitNULL,
                         save_fits = TRUE,
                         chains = 1, 
                         K = 5)
LinFitNULLkfoldPreds <- kfold_predict(LinFitNULLkfold)
LinFitNULLkfoldPredsDat <- LinFitNULLkfoldPreds$yrep
LinFitNULLkfoldPredsMean <- colMeans(LinFitNULLkfoldPredsDat)
LinFitNULLkfoldPredsMed <- apply(LinFitNULLkfoldPredsDat, 2, function(x){quantile(x, 0.5)})
LinFitNULLkfoldPredsLCB <- apply(LinFitNULLkfoldPredsDat, 2, function(x){quantile(x, 0.025)})
LinFitNULLkfoldPredsUCB <- apply(LinFitNULLkfoldPredsDat, 2, function(x){quantile(x, 0.975)})

LinFitNULLkfoldMetrics <- tibble(
  Fit = "LinFitNULL",
  MAE_kfold = mean(abs(LinFitNULLkfoldPredsMean - LinFitNULLkfoldPreds$y)),
  MAD_kfold = mean(abs(LinFitNULLkfoldPredsMed - LinFitNULLkfoldPreds$y)),
  COV_kfold = mean(LinFitNULLkfoldPredsLCB < LinFitNULLkfoldPreds$y & LinFitNULLkfoldPreds$y < LinFitNULLkfoldPredsUCB)
)
LinFitNULLkfoldMetrics

LinFitNULL2kfold <- kfold(LinFitNULL2,
                         save_fits = TRUE,
                         chains = 1, 
                         K = 5)

LinFitNULL2kfoldPreds <- kfold_predict(LinFitNULL2kfold)
LinFitNULL2kfoldPredsDat <- LinFitNULL2kfoldPreds$yrep
LinFitNULL2kfoldPredsMean <- colMeans(LinFitNULL2kfoldPredsDat)
LinFitNULL2kfoldPredsMed <- apply(LinFitNULL2kfoldPredsDat, 2, function(x){quantile(x, 0.5)})
LinFitNULL2kfoldPredsLCB <- apply(LinFitNULL2kfoldPredsDat, 2, function(x){quantile(x, 0.025)})
LinFitNULL2kfoldPredsUCB <- apply(LinFitNULL2kfoldPredsDat, 2, function(x){quantile(x, 0.975)})

LinFitNULL2kfoldMetrics <- tibble(
  Fit = "LinFitNULL2",
  MAE_kfold = mean(abs(LinFitNULL2kfoldPredsMean - LinFitNULL2kfoldPreds$y)),
  MAD_kfold = mean(abs(LinFitNULL2kfoldPredsMed - LinFitNULL2kfoldPreds$y)),
  COV_kfold = mean(LinFitNULL2kfoldPredsLCB < LinFitNULL2kfoldPreds$y & LinFitNULL2kfoldPreds$y < LinFitNULL2kfoldPredsUCB)
)
LinFitNULL2kfoldMetrics

rm(list = LinFitNULLfilesRM)
rm(LinFitNULLfilesRM)

### SKEW GAUSSIAN ----
#### Model 1 ----
skewLinFitNULL <- brm(
  bf(
    # VMAX ~ offset(HWRF), family = gaussian()
    # VMAX ~ HWRF, family = gaussian(link = "log")
    VMAX ~ offset(HWRF), family = skew_normal()
    # VMAX ~ log(HWRF), family = gaussian(link = "log)
  ),
  data = StormdataTrain8, 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)

skewLinFitNULL2 <- brm(
  bf(
    #VMAX ~ offset(HWRF), family = gaussian()
    # VMAX ~ offset(HWRF), family = gaussian(link = "log")
    # VMAX ~ log(HWRF), family = gaussian()
    VMAX ~ offset(log(HWRF)), family = gaussian(link = "log")
  ),
  data = StormdataTrain8, 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)

#save(skewLinFitNULL, file = "_data/skewLinFitNULL.RData")
prior_summary(skewLinFitNULL)
round(posterior_summary(skewLinFitNULL, probs = c(0.025, 0.975)))
skewLinFitNULL

print(skewLinFitNULL, digits = 4)
plot(skewLinFitNULL)
skewLinFitNULLppcFit <- pp_check(skewLinFitNULL, ndraws = 100) + 
  labs(title = "skewLinFitNULL Fit PPC") +
  theme_bw()
skewLinFitNULLppcFit
skewLinFitNULLloo <- loo(skewLinFitNULL)
waic(skewLinFitNULL)
performance::check_distribution(skewLinFitNULL)
performance::check_outliers(skewLinFitNULL)
performance::check_heteroskedasticity(skewLinFitNULL)
performance_rmse(skewLinFitNULL)
performance_mae(skewLinFitNULL)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(skewLinFitNULL)

print(skewLinFitNULL2, digits = 4)
plot(skewLinFitNULL2)
skewLinFitNULL2ppcFit <- pp_check(skewLinFitNULL2, ndraws = 100) + 
  labs(title = "skewLinFitNULL2 Fit PPC") +
  theme_bw()
skewLinFitNULL2ppcFit
skewLinFitNULL2loo <- loo(skewLinFitNULL2)
waic(skewLinFitNULL2)
performance::check_distribution(skewLinFitNULL2)
performance::check_outliers(skewLinFitNULL2)
performance::check_heteroskedasticity(skewLinFitNULL2)
performance_rmse(skewLinFitNULL2)
performance_mae(skewLinFitNULL2)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(skewLinFitNULL2)

skewLinFitNULLppcFit/skewLinFitNULL2ppcFit


variance_decomposition(skewLinFitNULL)
exp(fixef(skewLinFitNULL))
ranef(skewLinFitNULL)

bayes_R2(skewLinFitNULL)

bayes_factor(skewLinFitNULL, skewLinFitNULL2)
loo_compare(skewLinFitNULLloo,skewLinFitNULL2loo)

skewLinFitNULLsmooths <- conditional_smooths(skewLinFitNULL)
plot(skewLinFitNULLsmooths, stype = "raster", ask = FALSE)
skewLinFitNULLeffects <- conditional_effects(skewLinFitNULL, 
                                         method = "posterior_predict",
                                         robust = FALSE,
                                         re_formula = NULL)
skewLinFitNULLeffects <- conditional_effects(skewLinFitNULL)
plot(skewLinFitNULLeffects, points = TRUE, ask = FALSE)

##### Prediction ----
## Fitted
skewLinFitNULLfinalFit <- posterior_predict(skewLinFitNULL)
skewLinFitNULLfinalFitMean <- colMeans(skewLinFitNULLfinalFit)
skewLinFitNULLfinalFitMed <- apply(skewLinFitNULLfinalFit, 2, function(x){quantile(x, 0.5)})
skewLinFitNULLfinalFitLCB <- apply(skewLinFitNULLfinalFit, 2, function(x){quantile(x, 0.025)})
skewLinFitNULLfinalFitUCB <- apply(skewLinFitNULLfinalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
skewLinFitNULLfinalPreds <- posterior_predict(skewLinFitNULL, 
                                          newdata = StormdataTest8,
                                          allow_new_levels = TRUE)
skewLinFitNULLfinalPredsMean <- colMeans(skewLinFitNULLfinalPreds)
skewLinFitNULLfinalPredsMed <- apply(skewLinFitNULLfinalPreds, 2, function(x){quantile(x, 0.5)})
skewLinFitNULLfinalPredsLCB <- apply(skewLinFitNULLfinalPreds, 2, function(x){quantile(x, 0.025)})
skewLinFitNULLfinalPredsUCB <- apply(skewLinFitNULLfinalPreds, 2, function(x){quantile(x, 0.975)})

skewLinFitNULLpredMetrics <- tibble(
  Fit = "skewLinFitNULL",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(skewLinFitNULLfinalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(skewLinFitNULLfinalFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < skewLinFitNULLfinalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(skewLinFitNULLfinalPredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(skewLinFitNULLfinalPredsMed - Actual_Yvec)),
  COV_pred = mean(skewLinFitNULLfinalPredsLCB < Actual_Yvec & Actual_Yvec < skewLinFitNULLfinalPredsUCB)
)
skewLinFitNULLpredMetrics

## Fitted
skewLinFitNULL2finalFit <- posterior_predict(skewLinFitNULL2)
skewLinFitNULL2finalFitMean <- colMeans(skewLinFitNULL2finalFit)
skewLinFitNULL2finalFitMed <- apply(skewLinFitNULL2finalFit, 2, function(x){quantile(x, 0.5)})
skewLinFitNULL2finalFitLCB <- apply(skewLinFitNULL2finalFit, 2, function(x){quantile(x, 0.025)})
skewLinFitNULL2finalFitUCB <- apply(skewLinFitNULL2finalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
skewLinFitNULL2finalPreds <- posterior_predict(skewLinFitNULL2, 
                                           newdata = StormdataTest8,
                                           allow_new_levels = TRUE)
skewLinFitNULL2finalPredsMean <- colMeans(skewLinFitNULL2finalPreds)
skewLinFitNULL2finalPredsMed <- apply(skewLinFitNULL2finalPreds, 2, function(x){quantile(x, 0.5)})
skewLinFitNULL2finalPredsLCB <- apply(skewLinFitNULL2finalPreds, 2, function(x){quantile(x, 0.025)})
skewLinFitNULL2finalPredsUCB <- apply(skewLinFitNULL2finalPreds, 2, function(x){quantile(x, 0.975)})

skewLinFitNULL2predMetrics <- tibble(
  Fit = "skewLinFitNULL2",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(skewLinFitNULL2finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(skewLinFitNULL2finalFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < skewLinFitNULL2finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(skewLinFitNULL2finalPredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(skewLinFitNULL2finalPredsMed - Actual_Yvec)),
  COV_pred = mean(skewLinFitNULL2finalPredsLCB < Actual_Yvec & Actual_Yvec < skewLinFitNULL2finalPredsUCB)
)
skewLinFitNULLpredMetrics
skewLinFitNULL2predMetrics

##### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = skewLinFitNULLfinalPreds) +
  labs(title = "skewLinFitNULL Predict") +
  theme_bw()

skewLinFitNULLFitDF <- bind_cols(
  StormdataTrain3,
  LCB = skewLinFitNULLfinalFitLCB,
  Mean = skewLinFitNULLfinalFitMean,
  Med = skewLinFitNULLfinalFitMed,
  UCB = skewLinFitNULLfinalFitUCB
)

skewLinFitNULLstormsFitplot <- ggplot(data = skewLinFitNULLFitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))
skewLinFitNULLstormsFitplot

## Prediction
skewLinFitNULLPredDF <- bind_cols(
  StormdataTest3,
  LCB = skewLinFitNULLfinalPredsLCB,
  Mean = skewLinFitNULLfinalPredsMean,
  Med = skewLinFitNULLfinalPredsMed,
  UCB = skewLinFitNULLfinalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

skewLinFitNULLstormsPredplot <- ggplot(data = skewLinFitNULLPredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX, color = "Observed")) +
  geom_line(aes(y = Mean, color = "PPD Mean")) +
  facet_wrap(vars(StormID), ncol = 6)+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50)) +
  labs(title = "skewLinFitNULL PPD Mean vs Observed VMAX",
       subtitle = "95% Credible Interval about PPD Mean") +
  scale_color_manual(name = NULL, values = c("black","red")) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1))
  ) +
  theme_bw()
skewLinFitNULLstormsPredplot

##### PPC ----
###### Quantile 2.5 
skewLinFitNULLLCBsims <- apply(skewLinFitNULLfinalFit, 
                           MARGIN = 1,
                           function(x){
                             quantile(x, 0.025)
                           })
skewLinFitNULLLCBpvalueVec <- skewLinFitNULLLCBsims < quantile(StormdataTrain3$VMAX, 0.025)
skewLinFitNULLLCBpvalue <- sum(skewLinFitNULLLCBpvalueVec)
skewLinFitNULLLCBpvalue <- round(skewLinFitNULLLCBpvalue/4000, 3)
skewLinFitNULLLCBpvalue <- min(skewLinFitNULLLCBpvalue, 1 - skewLinFitNULLLCBpvalue)

skewLinFitNULL_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           skewLinFitNULLfinalFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", skewLinFitNULLLCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#skewLinFitNULL_ppcLCB

###### Quantile 97.5 
skewLinFitNULLUCBsims <- apply(skewLinFitNULLfinalFit, 
                           MARGIN = 1,
                           function(x){
                             quantile(x, 0.975)
                           })
skewLinFitNULLUCBpvalueVec <- skewLinFitNULLUCBsims < quantile(StormdataTrain3$VMAX, 0.975)
skewLinFitNULLUCBpvalue <- as.numeric(sum(skewLinFitNULLUCBpvalueVec))
skewLinFitNULLUCBpvalue <- round(skewLinFitNULLUCBpvalue/4000, 3)
skewLinFitNULLUCBpvalue <- min(skewLinFitNULLUCBpvalue, 1 - skewLinFitNULLUCBpvalue)

skewLinFitNULL_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           skewLinFitNULLfinalFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", skewLinFitNULLUCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#skewLinFitNULL_ppcUCB

###### Mean 
skewLinFitNULLMEANsims <- apply(skewLinFitNULLfinalFit, 
                            MARGIN = 1,
                            function(x){
                              mean(x)
                            })
skewLinFitNULLMEANpvalueVec <- skewLinFitNULLMEANsims < mean(StormdataTrain3$VMAX)
skewLinFitNULLMEANpvalue <- sum(skewLinFitNULLMEANpvalueVec)
skewLinFitNULLMEANpvalue <- round(skewLinFitNULLMEANpvalue/4000, 3)
skewLinFitNULLMEANpvalue <- min(skewLinFitNULLMEANpvalue, 1 - skewLinFitNULLMEANpvalue)

skewLinFitNULL_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           skewLinFitNULLfinalFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", skewLinFitNULLMEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#skewLinFitNULL_ppcMEAN

###### Med 
skewLinFitNULLMEDsims <- apply(skewLinFitNULLfinalFit, 
                           MARGIN = 1,
                           function(x){
                             quantile(x, 0.5)
                           })
skewLinFitNULLMEDpvalueVec <- skewLinFitNULLMEDsims < quantile(StormdataTrain3$VMAX, 0.5)
skewLinFitNULLMEDpvalue <- sum(skewLinFitNULLMEDpvalueVec)
skewLinFitNULLMEDpvalue <- round(skewLinFitNULLMEDpvalue/4000, 3)
skewLinFitNULLMEDpvalue <- min(skewLinFitNULLMEDpvalue, 1 - skewLinFitNULLMEDpvalue)

skewLinFitNULL_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           skewLinFitNULLfinalFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", skewLinFitNULLMEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#skewLinFitNULL_ppcMED

###### SD 
skewLinFitNULLSDsims <- apply(skewLinFitNULLfinalFit, 
                          MARGIN = 1,
                          function(x){
                            sd(x)
                          })
skewLinFitNULLSDpvalueVec <- skewLinFitNULLSDsims < sd(StormdataTrain3$VMAX)
skewLinFitNULLSDpvalue <- sum(skewLinFitNULLSDpvalueVec)
skewLinFitNULLSDpvalue <- round(skewLinFitNULLSDpvalue/4000, 3)
skewLinFitNULLSDpvalue <- min(skewLinFitNULLSDpvalue, 1 - skewLinFitNULLSDpvalue)

skewLinFitNULL_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           skewLinFitNULLfinalFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", skewLinFitNULLSDpvalue, ")")) +
  theme_bw() +
  legend_none()
#skewLinFitNULL_ppcSD

###### Range 
skewLinFitNULLRANGEsims <- apply(skewLinFitNULLfinalFit, 
                             MARGIN = 1,
                             function(x){
                               max(x)-min(x)
                             })
skewLinFitNULLRANGEpvalueVec <- skewLinFitNULLRANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
skewLinFitNULLRANGEpvalue <- sum(skewLinFitNULLRANGEpvalueVec)
skewLinFitNULLRANGEpvalue <- round(skewLinFitNULLRANGEpvalue/4000, 3)
skewLinFitNULLRANGEpvalue <- min(skewLinFitNULLRANGEpvalue, 1 - skewLinFitNULLRANGEpvalue)

skewLinFitNULL_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           skewLinFitNULLfinalFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", skewLinFitNULLRANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#skewLinFitNULL_ppcRANGE

##### Combined Plot ----
skewLinFitNULL_ppcComb <- 
  #skewLinFitNULLppcFit /
  skewLinFitNULL_ppcLCB | skewLinFitNULL_ppcMED | skewLinFitNULL_ppcUCB |
  skewLinFitNULL_ppcRANGE | skewLinFitNULL_ppcMEAN | skewLinFitNULL_ppcSD
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
skewLinFitNULL_ppcComb

##### Bayes p-values ----
skewLinFitNULLpvalues <- tibble(
  Fit = "skewLinFitNULL",
  LCB = skewLinFitNULLLCBpvalue,
  Median = skewLinFitNULLMEDpvalue,
  UCB = skewLinFitNULLUCBpvalue,
  Range = skewLinFitNULLRANGEpvalue,
  Mean = skewLinFitNULLMEANpvalue,
  SD = skewLinFitNULLSDpvalue
)
skewLinFitNULLpvalues

##### PPC ----
###### Quantile 2.5 
skewLinFitNULL2LCBsims <- apply(skewLinFitNULL2finalFit, 
                            MARGIN = 1,
                            function(x){
                              quantile(x, 0.025)
                            })
skewLinFitNULL2LCBpvalueVec <- skewLinFitNULL2LCBsims < quantile(StormdataTrain3$VMAX, 0.025)
skewLinFitNULL2LCBpvalue <- sum(skewLinFitNULL2LCBpvalueVec)
skewLinFitNULL2LCBpvalue <- round(skewLinFitNULL2LCBpvalue/4000, 3)
skewLinFitNULL2LCBpvalue <- min(skewLinFitNULL2LCBpvalue, 1 - skewLinFitNULL2LCBpvalue)

skewLinFitNULL2_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           skewLinFitNULL2finalFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", skewLinFitNULL2LCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#skewLinFitNULL2_ppcLCB

###### Quantile 97.5 
skewLinFitNULL2UCBsims <- apply(skewLinFitNULL2finalFit, 
                            MARGIN = 1,
                            function(x){
                              quantile(x, 0.975)
                            })
skewLinFitNULL2UCBpvalueVec <- skewLinFitNULL2UCBsims < quantile(StormdataTrain3$VMAX, 0.975)
skewLinFitNULL2UCBpvalue <- as.numeric(sum(skewLinFitNULL2UCBpvalueVec))
skewLinFitNULL2UCBpvalue <- round(skewLinFitNULL2UCBpvalue/4000, 3)
skewLinFitNULL2UCBpvalue <- min(skewLinFitNULL2UCBpvalue, 1 - skewLinFitNULL2UCBpvalue)

skewLinFitNULL2_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           skewLinFitNULL2finalFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", skewLinFitNULL2UCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#skewLinFitNULL2_ppcUCB

###### Mean 
skewLinFitNULL2MEANsims <- apply(skewLinFitNULL2finalFit, 
                             MARGIN = 1,
                             function(x){
                               mean(x)
                             })
skewLinFitNULL2MEANpvalueVec <- skewLinFitNULL2MEANsims < mean(StormdataTrain3$VMAX)
skewLinFitNULL2MEANpvalue <- sum(skewLinFitNULL2MEANpvalueVec)
skewLinFitNULL2MEANpvalue <- round(skewLinFitNULL2MEANpvalue/4000, 3)
skewLinFitNULL2MEANpvalue <- min(skewLinFitNULL2MEANpvalue, 1 - skewLinFitNULL2MEANpvalue)

skewLinFitNULL2_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           skewLinFitNULL2finalFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", skewLinFitNULL2MEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#skewLinFitNULL2_ppcMEAN

###### Med 
skewLinFitNULL2MEDsims <- apply(skewLinFitNULL2finalFit, 
                            MARGIN = 1,
                            function(x){
                              quantile(x, 0.5)
                            })
skewLinFitNULL2MEDpvalueVec <- skewLinFitNULL2MEDsims < quantile(StormdataTrain3$VMAX, 0.5)
skewLinFitNULL2MEDpvalue <- sum(skewLinFitNULL2MEDpvalueVec)
skewLinFitNULL2MEDpvalue <- round(skewLinFitNULL2MEDpvalue/4000, 3)
skewLinFitNULL2MEDpvalue <- min(skewLinFitNULL2MEDpvalue, 1 - skewLinFitNULL2MEDpvalue)

skewLinFitNULL2_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           skewLinFitNULL2finalFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", skewLinFitNULL2MEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#skewLinFitNULL2_ppcMED

###### SD 
skewLinFitNULL2SDsims <- apply(skewLinFitNULL2finalFit, 
                           MARGIN = 1,
                           function(x){
                             sd(x)
                           })
skewLinFitNULL2SDpvalueVec <- skewLinFitNULL2SDsims < sd(StormdataTrain3$VMAX)
skewLinFitNULL2SDpvalue <- sum(skewLinFitNULL2SDpvalueVec)
skewLinFitNULL2SDpvalue <- round(skewLinFitNULL2SDpvalue/4000, 3)
skewLinFitNULL2SDpvalue <- min(skewLinFitNULL2SDpvalue, 1 - skewLinFitNULL2SDpvalue)

skewLinFitNULL2_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           skewLinFitNULL2finalFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", skewLinFitNULL2SDpvalue, ")")) +
  theme_bw() +
  legend_none()
#skewLinFitNULL2_ppcSD

###### Range 
skewLinFitNULL2RANGEsims <- apply(skewLinFitNULL2finalFit, 
                              MARGIN = 1,
                              function(x){
                                max(x)-min(x)
                              })
skewLinFitNULL2RANGEpvalueVec <- skewLinFitNULL2RANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
skewLinFitNULL2RANGEpvalue <- sum(skewLinFitNULL2RANGEpvalueVec)
skewLinFitNULL2RANGEpvalue <- round(skewLinFitNULL2RANGEpvalue/4000, 3)
skewLinFitNULL2RANGEpvalue <- min(skewLinFitNULL2RANGEpvalue, 1 - skewLinFitNULL2RANGEpvalue)

skewLinFitNULL2_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           skewLinFitNULL2finalFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", skewLinFitNULL2RANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#skewLinFitNULL2_ppcRANGE

##### Combined Plot ----
skewLinFitNULL2_ppcComb <- 
  #skewLinFitNULL2ppcFit /
  skewLinFitNULL2_ppcLCB | skewLinFitNULL2_ppcMED | skewLinFitNULL2_ppcUCB |
  skewLinFitNULL2_ppcRANGE | skewLinFitNULL2_ppcMEAN | skewLinFitNULL2_ppcSD
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
skewLinFitNULL2_ppcComb

##### Bayes p-values ----
skewLinFitNULL2pvalues <- tibble(
  Fit = "skewLinFitNULL2",
  LCB = skewLinFitNULL2LCBpvalue,
  Median = skewLinFitNULL2MEDpvalue,
  UCB = skewLinFitNULL2UCBpvalue,
  Range = skewLinFitNULL2RANGEpvalue,
  Mean = skewLinFitNULL2MEANpvalue,
  SD = skewLinFitNULL2SDpvalue
)
skewLinFitNULL2pvalues


##### CV ----
skewLinFitNULLkfold <- kfold(skewLinFitNULL,
                         save_fits = TRUE,
                         chains = 1, 
                         K = 5)
skewLinFitNULLkfoldPreds <- kfold_predict(skewLinFitNULLkfold)
skewLinFitNULLkfoldPredsDat <- skewLinFitNULLkfoldPreds$yrep
skewLinFitNULLkfoldPredsMean <- colMeans(skewLinFitNULLkfoldPredsDat)
skewLinFitNULLkfoldPredsMed <- apply(skewLinFitNULLkfoldPredsDat, 2, function(x){quantile(x, 0.5)})
skewLinFitNULLkfoldPredsLCB <- apply(skewLinFitNULLkfoldPredsDat, 2, function(x){quantile(x, 0.025)})
skewLinFitNULLkfoldPredsUCB <- apply(skewLinFitNULLkfoldPredsDat, 2, function(x){quantile(x, 0.975)})

skewLinFitNULLkfoldMetrics <- tibble(
  Fit = "skewLinFitNULL",
  MAE_kfold = mean(abs(skewLinFitNULLkfoldPredsMean - skewLinFitNULLkfoldPreds$y)),
  MAD_kfold = mean(abs(skewLinFitNULLkfoldPredsMed - skewLinFitNULLkfoldPreds$y)),
  COV_kfold = mean(skewLinFitNULLkfoldPredsLCB < skewLinFitNULLkfoldPreds$y & skewLinFitNULLkfoldPreds$y < skewLinFitNULLkfoldPredsUCB)
)
skewLinFitNULLkfoldMetrics

skewLinFitNULL2kfold <- kfold(skewLinFitNULL2,
                          save_fits = TRUE,
                          chains = 1, 
                          K = 5)

skewLinFitNULL2kfoldPreds <- kfold_predict(skewLinFitNULL2kfold)
skewLinFitNULL2kfoldPredsDat <- skewLinFitNULL2kfoldPreds$yrep
skewLinFitNULL2kfoldPredsMean <- colMeans(skewLinFitNULL2kfoldPredsDat)
skewLinFitNULL2kfoldPredsMed <- apply(skewLinFitNULL2kfoldPredsDat, 2, function(x){quantile(x, 0.5)})
skewLinFitNULL2kfoldPredsLCB <- apply(skewLinFitNULL2kfoldPredsDat, 2, function(x){quantile(x, 0.025)})
skewLinFitNULL2kfoldPredsUCB <- apply(skewLinFitNULL2kfoldPredsDat, 2, function(x){quantile(x, 0.975)})

skewLinFitNULL2kfoldMetrics <- tibble(
  Fit = "skewLinFitNULL2",
  MAE_kfold = mean(abs(skewLinFitNULL2kfoldPredsMean - skewLinFitNULL2kfoldPreds$y)),
  MAD_kfold = mean(abs(skewLinFitNULL2kfoldPredsMed - skewLinFitNULL2kfoldPreds$y)),
  COV_kfold = mean(skewLinFitNULL2kfoldPredsLCB < skewLinFitNULL2kfoldPreds$y & skewLinFitNULL2kfoldPreds$y < skewLinFitNULL2kfoldPredsUCB)
)
skewLinFitNULL2kfoldMetrics


skewLinFitNULLfiles <- ls()[str_detect(ls(), pattern = "skewLinFitNULL")]
skewLinFitNULLfilesRM <- skewLinFitNULLfiles[!(skewLinFitNULLfiles %in% c("skewLinFitNULL_ppcComb",
                                                              "skewLinFitNULLloo",
                                                              "skewLinFitNULLpredMetrics",
                                                              "skewLinFitNULLpvalues",
                                                              "skewLinFitNULLstormsPredplot"))]

rm(list = skewLinFitNULLfilesRM)
rm(skewLinFitNULLfilesRM)

### STUDENT ----
StudentFitNULL <- brm(
  bf(
    # VMAX ~ offset(HWRF), family = gaussian()
    # VMAX ~ HWRF, family = gaussian(link = "log")
    VMAX ~ offset(HWRF), family = student()
    # VMAX ~ log(HWRF), family = gaussian(link = "log)
  ),
  data = StormdataTrain8, 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)

StudentFitNULL2 <- brm(
  bf(
    # VMAX ~ offset(HWRF), family = gaussian()
    # VMAX ~ HWRF, family = gaussian(link = "log")
    VMAX ~ offset(log(HWRF)), family = student(link = "log")
    # VMAX ~ log(HWRF), family = gaussian(link = "log)
  ),
  data = StormdataTrain8, 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)

StudentFitNULL3 <- brm(
  bf(
    # VMAX ~ offset(HWRF), family = gaussian()
    # VMAX ~ HWRF, family = gaussian(link = "log")
    VMAX ~ offset(log(HWRF)), family = student(link = "log")
    # VMAX ~ log(HWRF), family = gaussian(link = "log)
  ),
  data = StormdataTrain8, 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)

#save(StudentFitNULL, file = "_data/StudentFitNULL.RData")
prior_summary(StudentFitNULL)
round(posterior_summary(StudentFitNULL, probs = c(0.025, 0.975)))
StudentFitNULL

print(StudentFitNULL3, digits = 4)
plot(StudentFitNULL)
StudentFitNULLppcFit <- pp_check(StudentFitNULL3, ndraws = 100) + 
  labs(title = "StudentFitNULL Fit3 PPC") +
  theme_bw()
StudentFitNULLppcFit
StudentFitNULLloo <- loo(StudentFitNULL)
waic(StudentFitNULL)
performance::check_distribution(StudentFitNULL)
performance::check_outliers(StudentFitNULL)
performance::check_heteroskedasticity(StudentFitNULL)
performance_rmse(StudentFitNULL)
performance_mae(StudentFitNULL3)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(StudentFitNULL)

print(StudentFitNULL2, digits = 4)
plot(StudentFitNULL2)
StudentFitNULL2ppcFit <- pp_check(StudentFitNULL2, ndraws = 100) + 
  labs(title = "StudentFitNULL2 Fit PPC") +
  theme_bw()
StudentFitNULL2ppcFit
StudentFitNULL2loo <- loo(StudentFitNULL2)
waic(StudentFitNULL2)
performance::check_distribution(StudentFitNULL2)
performance::check_outliers(StudentFitNULL2)
performance::check_heteroskedasticity(StudentFitNULL2)
performance_rmse(StudentFitNULL2)
performance_mae(StudentFitNULL2)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(StudentFitNULL2)


variance_decomposition(StudentFitNULL)
exp(fixef(StudentFitNULL))
ranef(StudentFitNULL)

bayes_R2(StudentFitNULL)

bayes_factor(StudentFitNULL2, StudentFitNULL)
bayes_factor(StudentFitNULL, gammaFit2)
bayes_factor(StudentFitNULL, gammaFit3)
bayes_factor(StudentFitNULL, studentFit1)
bayes_factor(StudentFitNULL, studentFit2)
bayes_factor(StudentFitNULL, studentFit3)
bayes_factor(StudentFitNULL, linFit11)
bayes_factor(StudentFitNULL, propFit1)
bayes_factor(StudentFitNULL, logPropFit1)
loo(StudentFitNULL, gammaFit3)

StudentFitNULLsmooths <- conditional_smooths(StudentFitNULL)
plot(StudentFitNULLsmooths, stype = "raster", ask = FALSE)
StudentFitNULLeffects <- conditional_effects(StudentFitNULL, 
                                         method = "posterior_predict",
                                         robust = FALSE,
                                         re_formula = NULL)
StudentFitNULLeffects <- conditional_effects(StudentFitNULL)
plot(StudentFitNULLeffects, points = TRUE, ask = FALSE)

##### Prediction ----
## Fitted
StudentFitNULLfinalFit <- posterior_predict(StudentFitNULL3)
StudentFitNULLfinalFitMean <- colMeans(StudentFitNULLfinalFit)
StudentFitNULLfinalFitMed <- apply(StudentFitNULLfinalFit, 2, function(x){quantile(x, 0.5)})
StudentFitNULLfinalFitLCB <- apply(StudentFitNULLfinalFit, 2, function(x){quantile(x, 0.025)})
StudentFitNULLfinalFitUCB <- apply(StudentFitNULLfinalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
StudentFitNULLfinalPreds <- posterior_predict(StudentFitNULL3, 
                                          newdata = StormdataTest8,
                                          allow_new_levels = TRUE)
StudentFitNULLfinalPredsMean <- colMeans(StudentFitNULLfinalPreds)
StudentFitNULLfinalPredsMed <- apply(StudentFitNULLfinalPreds, 2, function(x){quantile(x, 0.5)})
StudentFitNULLfinalPredsLCB <- apply(StudentFitNULLfinalPreds, 2, function(x){quantile(x, 0.025)})
StudentFitNULLfinalPredsUCB <- apply(StudentFitNULLfinalPreds, 2, function(x){quantile(x, 0.975)})

StudentFitNULLpredMetrics <- tibble(
  Fit = "StudentFitNULL",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(StudentFitNULLfinalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(StudentFitNULLfinalFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < StudentFitNULLfinalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(StudentFitNULLfinalPredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(StudentFitNULLfinalPredsMed - Actual_Yvec)),
  COV_pred = mean(StudentFitNULLfinalPredsLCB < Actual_Yvec & Actual_Yvec < StudentFitNULLfinalPredsUCB)
)
StudentFitNULLpredMetrics

## Fitted
StudentFitNULL2finalFit <- posterior_predict(StudentFitNULL2)
StudentFitNULL2finalFitMean <- colMeans(StudentFitNULL2finalFit)
StudentFitNULL2finalFitMed <- apply(StudentFitNULL2finalFit, 2, function(x){quantile(x, 0.5)})
StudentFitNULL2finalFitLCB <- apply(StudentFitNULL2finalFit, 2, function(x){quantile(x, 0.025)})
StudentFitNULL2finalFitUCB <- apply(StudentFitNULL2finalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
StudentFitNULL2finalPreds <- posterior_predict(StudentFitNULL2, 
                                              newdata = StormdataTest8,
                                              allow_new_levels = TRUE)
StudentFitNULL2finalPredsMean <- colMeans(StudentFitNULL2finalPreds)
StudentFitNULL2finalPredsMed <- apply(StudentFitNULL2finalPreds, 2, function(x){quantile(x, 0.5)})
StudentFitNULL2finalPredsLCB <- apply(StudentFitNULL2finalPreds, 2, function(x){quantile(x, 0.025)})
StudentFitNULL2finalPredsUCB <- apply(StudentFitNULL2finalPreds, 2, function(x){quantile(x, 0.975)})

StudentFitNULL2predMetrics <- tibble(
  Fit = "StudentFitNULL2",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(StudentFitNULL2finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(StudentFitNULL2finalFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < StudentFitNULL2finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(StudentFitNULL2finalPredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(StudentFitNULL2finalPredsMed - Actual_Yvec)),
  COV_pred = mean(StudentFitNULL2finalPredsLCB < Actual_Yvec & Actual_Yvec < StudentFitNULL2finalPredsUCB)
)
StudentFitNULLpredMetrics
StudentFitNULL2predMetrics

##### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = StudentFitNULLfinalPreds) +
  labs(title = "StudentFitNULL Predict") +
  theme_bw()

StudentFitNULLFitDF <- bind_cols(
  StormdataTrain3,
  LCB = StudentFitNULLfinalFitLCB,
  Mean = StudentFitNULLfinalFitMean,
  Med = StudentFitNULLfinalFitMed,
  UCB = StudentFitNULLfinalFitUCB
)

StudentFitNULLstormsFitplot <- ggplot(data = StudentFitNULLFitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))
StudentFitNULLstormsFitplot

## Prediction
StudentFitNULLPredDF <- bind_cols(
  StormdataTest3,
  LCB = StudentFitNULLfinalPredsLCB,
  Mean = StudentFitNULLfinalPredsMean,
  Med = StudentFitNULLfinalPredsMed,
  UCB = StudentFitNULLfinalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

StudentFitNULLstormsPredplot <- ggplot(data = StudentFitNULLPredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX, color = "Observed")) +
  geom_line(aes(y = Mean, color = "PPD Mean")) +
  facet_wrap(vars(StormID), ncol = 6)+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50)) +
  labs(title = "StudentFitNULL PPD Mean vs Observed VMAX",
       subtitle = "95% Credible Interval about PPD Mean") +
  scale_color_manual(name = NULL, values = c("black","red")) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1))
  ) +
  theme_bw()
StudentFitNULLstormsPredplot

##### PPC ----
###### Quantile 2.5 
StudentFitNULLLCBsims <- apply(StudentFitNULLfinalFit, 
                           MARGIN = 1,
                           function(x){
                             quantile(x, 0.025)
                           })
StudentFitNULLLCBpvalueVec <- StudentFitNULLLCBsims < quantile(StormdataTrain3$VMAX, 0.025)
StudentFitNULLLCBpvalue <- sum(StudentFitNULLLCBpvalueVec)
StudentFitNULLLCBpvalue <- round(StudentFitNULLLCBpvalue/4000, 3)
StudentFitNULLLCBpvalue <- min(StudentFitNULLLCBpvalue, 1 - StudentFitNULLLCBpvalue)

StudentFitNULL_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           StudentFitNULLfinalFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", StudentFitNULLLCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#StudentFitNULL_ppcLCB

###### Quantile 97.5 
StudentFitNULLUCBsims <- apply(StudentFitNULLfinalFit, 
                           MARGIN = 1,
                           function(x){
                             quantile(x, 0.975)
                           })
StudentFitNULLUCBpvalueVec <- StudentFitNULLUCBsims < quantile(StormdataTrain3$VMAX, 0.975)
StudentFitNULLUCBpvalue <- as.numeric(sum(StudentFitNULLUCBpvalueVec))
StudentFitNULLUCBpvalue <- round(StudentFitNULLUCBpvalue/4000, 3)
StudentFitNULLUCBpvalue <- min(StudentFitNULLUCBpvalue, 1 - StudentFitNULLUCBpvalue)

StudentFitNULL_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           StudentFitNULLfinalFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", StudentFitNULLUCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#StudentFitNULL_ppcUCB

###### Mean 
StudentFitNULLMEANsims <- apply(StudentFitNULLfinalFit, 
                            MARGIN = 1,
                            function(x){
                              mean(x)
                            })
StudentFitNULLMEANpvalueVec <- StudentFitNULLMEANsims < mean(StormdataTrain3$VMAX)
StudentFitNULLMEANpvalue <- sum(StudentFitNULLMEANpvalueVec)
StudentFitNULLMEANpvalue <- round(StudentFitNULLMEANpvalue/4000, 3)
StudentFitNULLMEANpvalue <- min(StudentFitNULLMEANpvalue, 1 - StudentFitNULLMEANpvalue)

StudentFitNULL_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           StudentFitNULLfinalFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", StudentFitNULLMEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#StudentFitNULL_ppcMEAN

###### Med 
StudentFitNULLMEDsims <- apply(StudentFitNULLfinalFit, 
                           MARGIN = 1,
                           function(x){
                             quantile(x, 0.5)
                           })
StudentFitNULLMEDpvalueVec <- StudentFitNULLMEDsims < quantile(StormdataTrain3$VMAX, 0.5)
StudentFitNULLMEDpvalue <- sum(StudentFitNULLMEDpvalueVec)
StudentFitNULLMEDpvalue <- round(StudentFitNULLMEDpvalue/4000, 3)
StudentFitNULLMEDpvalue <- min(StudentFitNULLMEDpvalue, 1 - StudentFitNULLMEDpvalue)

StudentFitNULL_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           StudentFitNULLfinalFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", StudentFitNULLMEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#StudentFitNULL_ppcMED

###### SD 
StudentFitNULLSDsims <- apply(StudentFitNULLfinalFit, 
                          MARGIN = 1,
                          function(x){
                            sd(x)
                          })
StudentFitNULLSDpvalueVec <- StudentFitNULLSDsims < sd(StormdataTrain3$VMAX)
StudentFitNULLSDpvalue <- sum(StudentFitNULLSDpvalueVec)
StudentFitNULLSDpvalue <- round(StudentFitNULLSDpvalue/4000, 3)
StudentFitNULLSDpvalue <- min(StudentFitNULLSDpvalue, 1 - StudentFitNULLSDpvalue)

StudentFitNULL_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           StudentFitNULLfinalFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", StudentFitNULLSDpvalue, ")")) +
  theme_bw() +
  legend_none()
#StudentFitNULL_ppcSD

###### Range 
StudentFitNULLRANGEsims <- apply(StudentFitNULLfinalFit, 
                             MARGIN = 1,
                             function(x){
                               max(x)-min(x)
                             })
StudentFitNULLRANGEpvalueVec <- StudentFitNULLRANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
StudentFitNULLRANGEpvalue <- sum(StudentFitNULLRANGEpvalueVec)
StudentFitNULLRANGEpvalue <- round(StudentFitNULLRANGEpvalue/4000, 3)
StudentFitNULLRANGEpvalue <- min(StudentFitNULLRANGEpvalue, 1 - StudentFitNULLRANGEpvalue)

StudentFitNULL_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           StudentFitNULLfinalFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", StudentFitNULLRANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#StudentFitNULL_ppcRANGE

##### Combined Plot ----
StudentFitNULL_ppcComb <- 
  #StudentFitNULLppcFit /
  (StudentFitNULL_ppcLCB | StudentFitNULL_ppcMED | StudentFitNULL_ppcUCB |
  StudentFitNULL_ppcRANGE | StudentFitNULL_ppcMEAN | StudentFitNULL_ppcSD)
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
StudentFitNULL_ppcComb

##### Bayes p-values ----
StudentFitNULLpvalues <- tibble(
  Fit = "StudentFitNULL",
  LCB = StudentFitNULLLCBpvalue,
  Median = StudentFitNULLMEDpvalue,
  UCB = StudentFitNULLUCBpvalue,
  Range = StudentFitNULLRANGEpvalue,
  Mean = StudentFitNULLMEANpvalue,
  SD = StudentFitNULLSDpvalue
)
StudentFitNULLpvalues

##### PPC ----
###### Quantile 2.5 
StudentFitNULL2LCBsims <- apply(StudentFitNULL2finalFit, 
                               MARGIN = 1,
                               function(x){
                                 quantile(x, 0.025)
                               })
StudentFitNULL2LCBpvalueVec <- StudentFitNULL2LCBsims < quantile(StormdataTrain3$VMAX, 0.025)
StudentFitNULL2LCBpvalue <- sum(StudentFitNULL2LCBpvalueVec)
StudentFitNULL2LCBpvalue <- round(StudentFitNULL2LCBpvalue/4000, 3)
StudentFitNULL2LCBpvalue <- min(StudentFitNULL2LCBpvalue, 1 - StudentFitNULL2LCBpvalue)

StudentFitNULL2_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           StudentFitNULL2finalFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", StudentFitNULL2LCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#StudentFitNULL2_ppcLCB

###### Quantile 97.5 
StudentFitNULL2UCBsims <- apply(StudentFitNULL2finalFit, 
                               MARGIN = 1,
                               function(x){
                                 quantile(x, 0.975)
                               })
StudentFitNULL2UCBpvalueVec <- StudentFitNULL2UCBsims < quantile(StormdataTrain3$VMAX, 0.975)
StudentFitNULL2UCBpvalue <- as.numeric(sum(StudentFitNULL2UCBpvalueVec))
StudentFitNULL2UCBpvalue <- round(StudentFitNULL2UCBpvalue/4000, 3)
StudentFitNULL2UCBpvalue <- min(StudentFitNULL2UCBpvalue, 1 - StudentFitNULL2UCBpvalue)

StudentFitNULL2_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           StudentFitNULL2finalFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", StudentFitNULL2UCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#StudentFitNULL2_ppcUCB

###### Mean 
StudentFitNULL2MEANsims <- apply(StudentFitNULL2finalFit, 
                                MARGIN = 1,
                                function(x){
                                  mean(x)
                                })
StudentFitNULL2MEANpvalueVec <- StudentFitNULL2MEANsims < mean(StormdataTrain3$VMAX)
StudentFitNULL2MEANpvalue <- sum(StudentFitNULL2MEANpvalueVec)
StudentFitNULL2MEANpvalue <- round(StudentFitNULL2MEANpvalue/4000, 3)
StudentFitNULL2MEANpvalue <- min(StudentFitNULL2MEANpvalue, 1 - StudentFitNULL2MEANpvalue)

StudentFitNULL2_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           StudentFitNULL2finalFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", StudentFitNULL2MEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#StudentFitNULL2_ppcMEAN

###### Med 
StudentFitNULL2MEDsims <- apply(StudentFitNULL2finalFit, 
                               MARGIN = 1,
                               function(x){
                                 quantile(x, 0.5)
                               })
StudentFitNULL2MEDpvalueVec <- StudentFitNULL2MEDsims < quantile(StormdataTrain3$VMAX, 0.5)
StudentFitNULL2MEDpvalue <- sum(StudentFitNULL2MEDpvalueVec)
StudentFitNULL2MEDpvalue <- round(StudentFitNULL2MEDpvalue/4000, 3)
StudentFitNULL2MEDpvalue <- min(StudentFitNULL2MEDpvalue, 1 - StudentFitNULL2MEDpvalue)

StudentFitNULL2_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           StudentFitNULL2finalFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", StudentFitNULL2MEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#StudentFitNULL2_ppcMED

###### SD 
StudentFitNULL2SDsims <- apply(StudentFitNULL2finalFit, 
                              MARGIN = 1,
                              function(x){
                                sd(x)
                              })
StudentFitNULL2SDpvalueVec <- StudentFitNULL2SDsims < sd(StormdataTrain3$VMAX)
StudentFitNULL2SDpvalue <- sum(StudentFitNULL2SDpvalueVec)
StudentFitNULL2SDpvalue <- round(StudentFitNULL2SDpvalue/4000, 3)
StudentFitNULL2SDpvalue <- min(StudentFitNULL2SDpvalue, 1 - StudentFitNULL2SDpvalue)

StudentFitNULL2_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           StudentFitNULL2finalFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", StudentFitNULL2SDpvalue, ")")) +
  theme_bw() +
  legend_none()
#StudentFitNULL2_ppcSD

###### Range 
StudentFitNULL2RANGEsims <- apply(StudentFitNULL2finalFit, 
                                 MARGIN = 1,
                                 function(x){
                                   max(x)-min(x)
                                 })
StudentFitNULL2RANGEpvalueVec <- StudentFitNULL2RANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
StudentFitNULL2RANGEpvalue <- sum(StudentFitNULL2RANGEpvalueVec)
StudentFitNULL2RANGEpvalue <- round(StudentFitNULL2RANGEpvalue/4000, 3)
StudentFitNULL2RANGEpvalue <- min(StudentFitNULL2RANGEpvalue, 1 - StudentFitNULL2RANGEpvalue)

StudentFitNULL2_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           StudentFitNULL2finalFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", StudentFitNULL2RANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#StudentFitNULL2_ppcRANGE

##### Combined Plot ----
StudentFitNULL2_ppcComb <- 
  #StudentFitNULL2ppcFit /
  (StudentFitNULL2_ppcLCB | StudentFitNULL2_ppcMED | StudentFitNULL2_ppcUCB |
     StudentFitNULL2_ppcRANGE | StudentFitNULL2_ppcMEAN | StudentFitNULL2_ppcSD)
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
StudentFitNULL2_ppcComb

##### Bayes p-values ----
StudentFitNULL2pvalues <- tibble(
  Fit = "StudentFitNULL2",
  LCB = StudentFitNULL2LCBpvalue,
  Median = StudentFitNULL2MEDpvalue,
  UCB = StudentFitNULL2UCBpvalue,
  Range = StudentFitNULL2RANGEpvalue,
  Mean = StudentFitNULL2MEANpvalue,
  SD = StudentFitNULL2SDpvalue
)
StudentFitNULL2pvalues

##### CV ----
StudentFitNULLkfold <- kfold(StudentFitNULL,
                         save_fits = TRUE,
                         chains = 1, 
                         K = 5)
StudentFitNULLkfoldPreds <- kfold_predict(StudentFitNULLkfold)
StudentFitNULLkfoldPredsDat <- StudentFitNULLkfoldPreds$yrep
StudentFitNULLkfoldPredsMean <- colMeans(StudentFitNULLkfoldPredsDat)
StudentFitNULLkfoldPredsMed <- apply(StudentFitNULLkfoldPredsDat, 2, function(x){quantile(x, 0.5)})
StudentFitNULLkfoldPredsLCB <- apply(StudentFitNULLkfoldPredsDat, 2, function(x){quantile(x, 0.025)})
StudentFitNULLkfoldPredsUCB <- apply(StudentFitNULLkfoldPredsDat, 2, function(x){quantile(x, 0.975)})

StudentFitNULLkfoldMetrics <- tibble(
  Fit = "StudentFitNULL",
  MAE_kfold = mean(abs(StudentFitNULLkfoldPredsMean - StudentFitNULLkfoldPreds$y)),
  MAD_kfold = mean(abs(StudentFitNULLkfoldPredsMed - StudentFitNULLkfoldPreds$y)),
  COV_kfold = mean(StudentFitNULLkfoldPredsLCB < StudentFitNULLkfoldPreds$y & StudentFitNULLkfoldPreds$y < StudentFitNULLkfoldPredsUCB)
)
StudentFitNULLkfoldMetrics

StudentFitNULL2kfold <- kfold(StudentFitNULL2,
                          save_fits = TRUE,
                          chains = 1, 
                          K = 5)

StudentFitNULL2kfoldPreds <- kfold_predict(StudentFitNULL2kfold)
StudentFitNULL2kfoldPredsDat <- StudentFitNULL2kfoldPreds$yrep
StudentFitNULL2kfoldPredsMean <- colMeans(StudentFitNULL2kfoldPredsDat)
StudentFitNULL2kfoldPredsMed <- apply(StudentFitNULL2kfoldPredsDat, 2, function(x){quantile(x, 0.5)})
StudentFitNULL2kfoldPredsLCB <- apply(StudentFitNULL2kfoldPredsDat, 2, function(x){quantile(x, 0.025)})
StudentFitNULL2kfoldPredsUCB <- apply(StudentFitNULL2kfoldPredsDat, 2, function(x){quantile(x, 0.975)})

StudentFitNULL2kfoldMetrics <- tibble(
  Fit = "StudentFitNULL2",
  MAE_kfold = mean(abs(StudentFitNULL2kfoldPredsMean - StudentFitNULL2kfoldPreds$y)),
  MAD_kfold = mean(abs(StudentFitNULL2kfoldPredsMed - StudentFitNULL2kfoldPreds$y)),
  COV_kfold = mean(StudentFitNULL2kfoldPredsLCB < StudentFitNULL2kfoldPreds$y & StudentFitNULL2kfoldPreds$y < StudentFitNULL2kfoldPredsUCB)
)
StudentFitNULL2kfoldMetrics

StudentFitNULLfiles <- ls()[str_detect(ls(), pattern = "StudentFitNULL")]
StudentFitNULLfilesRM <- StudentFitNULLfiles[!(StudentFitNULLfiles %in% c("StudentFitNULL_ppcComb",
                                                              "StudentFitNULLloo",
                                                              "StudentFitNULLpredMetrics",
                                                              "StudentFitNULLpvalues",
                                                              "StudentFitNULLstormsPredplot"))]

rm(list = StudentFitNULLfilesRM)
rm(StudentFitNULLfilesRM)

### LOGNORMAL ----
logNormalFitNULL <- brm(
  bf(
    VMAX ~ offset(HWRF), family = lognormal()
    # VMAX ~ log(HWRF), family = lognormal(),
    # VMAX ~ log(HWRF), family = lognormal(),
    # VMAX ~ log(HWRF), family = lognormal(),
  ),
  data = StormdataTrain8, 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)

logNormalFitNULL2 <- brm(
  bf(
    VMAX ~ offset(log(HWRF)), family = lognormal()
    # VMAX ~ log(HWRF), family = lognormal(),
    # VMAX ~ log(HWRF), family = lognormal(),
    # VMAX ~ log(HWRF), family = lognormal(),
  ),
  data = StormdataTrain8, 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)

#save(logNormalFitNULL, file = "_data/logNormalFitNULL.RData")
prior_summary(logNormalFitNULL)
round(posterior_summary(logNormalFitNULL, probs = c(0.025, 0.975)))
logNormalFitNULL

print(logNormalFitNULL2, digits = 4)
plot(logNormalFitNULL)
logNormalFitNULLppcFit <- pp_check(logNormalFitNULL, ndraws = 100) + 
  labs(title = "logNormalFitNULL Fit PPC") +
  theme_bw()
logNormalFitNULLppcFit
logNormalFitNULLloo <- loo(logNormalFitNULL)
waic(logNormalFitNULL)
performance::check_distribution(logNormalFitNULL)
performance::check_outliers(logNormalFitNULL)
performance::check_heteroskedasticity(logNormalFitNULL)
performance_rmse(logNormalFitNULL)
performance_mae(logNormalFitNULL)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(logNormalFitNULL)

print(logNormalFitNULL2, digits = 4)
plot(logNormalFitNULL2)
logNormalFitNULL2ppcFit <- pp_check(logNormalFitNULL2, ndraws = 100) + 
  labs(title = "logNormalFitNULL2 Fit PPC") +
  theme_bw()
logNormalFitNULL2ppcFit
logNormalFitNULL2loo <- loo(logNormalFitNULL2)
waic(logNormalFitNULL2)
performance::check_distribution(logNormalFitNULL2)
performance::check_outliers(logNormalFitNULL2)
performance::check_heteroskedasticity(logNormalFitNULL2)
performance_rmse(logNormalFitNULL2)
performance_mae(logNormalFitNULL2)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(logNormalFitNULL2)


variance_decomposition(logNormalFitNULL)
exp(fixef(logNormalFitNULL))
ranef(logNormalFitNULL)

bayes_R2(logNormalFitNULL)

bayes_factor(logNormalFitNULL2, logNormalFitNULL)

logNormalFitNULLsmooths <- conditional_smooths(logNormalFitNULL)
plot(logNormalFitNULLsmooths, stype = "raster", ask = FALSE)
logNormalFitNULLeffects <- conditional_effects(logNormalFitNULL, 
                                             method = "posterior_predict",
                                             robust = FALSE,
                                             re_formula = NULL)
logNormalFitNULLeffects <- conditional_effects(logNormalFitNULL)
plot(logNormalFitNULLeffects, points = TRUE, ask = FALSE)

##### Prediction ----
## Fitted
logNormalFitNULLfinalFit <- posterior_predict(logNormalFitNULL)
logNormalFitNULLfinalFitMean <- colMeans(logNormalFitNULLfinalFit)
logNormalFitNULLfinalFitMed <- apply(logNormalFitNULLfinalFit, 2, function(x){quantile(x, 0.5)})
logNormalFitNULLfinalFitLCB <- apply(logNormalFitNULLfinalFit, 2, function(x){quantile(x, 0.025)})
logNormalFitNULLfinalFitUCB <- apply(logNormalFitNULLfinalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
logNormalFitNULLfinalPreds <- posterior_predict(logNormalFitNULL, 
                                              newdata = StormdataTest8,
                                              allow_new_levels = TRUE)
logNormalFitNULLfinalPredsMean <- colMeans(logNormalFitNULLfinalPreds)
logNormalFitNULLfinalPredsMed <- apply(logNormalFitNULLfinalPreds, 2, function(x){quantile(x, 0.5)})
logNormalFitNULLfinalPredsLCB <- apply(logNormalFitNULLfinalPreds, 2, function(x){quantile(x, 0.025)})
logNormalFitNULLfinalPredsUCB <- apply(logNormalFitNULLfinalPreds, 2, function(x){quantile(x, 0.975)})

logNormalFitNULLpredMetrics <- tibble(
  Fit = "logNormalFitNULL",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(logNormalFitNULLfinalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(logNormalFitNULLfinalFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < logNormalFitNULLfinalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(logNormalFitNULLfinalPredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(logNormalFitNULLfinalPredsMed - Actual_Yvec)),
  COV_pred = mean(logNormalFitNULLfinalPredsLCB < Actual_Yvec & Actual_Yvec < logNormalFitNULLfinalPredsUCB)
)
logNormalFitNULLpredMetrics

## Fitted
logNormalFitNULL2finalFit <- posterior_predict(logNormalFitNULL2)
logNormalFitNULL2finalFitMean <- colMeans(logNormalFitNULL2finalFit)
logNormalFitNULL2finalFitMed <- apply(logNormalFitNULL2finalFit, 2, function(x){quantile(x, 0.5)})
logNormalFitNULL2finalFitLCB <- apply(logNormalFitNULL2finalFit, 2, function(x){quantile(x, 0.025)})
logNormalFitNULL2finalFitUCB <- apply(logNormalFitNULL2finalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
logNormalFitNULL2finalPreds <- posterior_predict(logNormalFitNULL2, 
                                                newdata = StormdataTest8,
                                                allow_new_levels = TRUE)
logNormalFitNULL2finalPredsMean <- colMeans(logNormalFitNULL2finalPreds)
logNormalFitNULL2finalPredsMed <- apply(logNormalFitNULL2finalPreds, 2, function(x){quantile(x, 0.5)})
logNormalFitNULL2finalPredsLCB <- apply(logNormalFitNULL2finalPreds, 2, function(x){quantile(x, 0.025)})
logNormalFitNULL2finalPredsUCB <- apply(logNormalFitNULL2finalPreds, 2, function(x){quantile(x, 0.975)})

logNormalFitNULL2predMetrics <- tibble(
  Fit = "logNormalFitNULL2",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(logNormalFitNULL2finalFitMean - StormdataTrain3$VMAX)),
  MAD_fit = mean(abs(logNormalFitNULL2finalFitMed - StormdataTrain3$VMAX)),
  COV_fit = mean(logNormalFitNULL2finalFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < logNormalFitNULL2finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(logNormalFitNULL2finalPredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(logNormalFitNULL2finalPredsMed - Actual_Yvec)),
  COV_pred = mean(logNormalFitNULL2finalPredsLCB < Actual_Yvec & Actual_Yvec < logNormalFitNULL2finalPredsUCB)
)
logNormalFitNULL2predMetrics

##### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = logNormalFitNULLfinalPreds) +
  labs(title = "logNormalFitNULL Predict") +
  theme_bw()

logNormalFitNULLFitDF <- bind_cols(
  StormdataTrain3,
  LCB = logNormalFitNULLfinalFitLCB,
  Mean = logNormalFitNULLfinalFitMean,
  Med = logNormalFitNULLfinalFitMed,
  UCB = logNormalFitNULLfinalFitUCB
)

logNormalFitNULLstormsFitplot <- ggplot(data = logNormalFitNULLFitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))
logNormalFitNULLstormsFitplot

## Prediction
logNormalFitNULLPredDF <- bind_cols(
  StormdataTest3,
  LCB = logNormalFitNULLfinalPredsLCB,
  Mean = logNormalFitNULLfinalPredsMean,
  Med = logNormalFitNULLfinalPredsMed,
  UCB = logNormalFitNULLfinalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

logNormalFitNULLstormsPredplot <- ggplot(data = logNormalFitNULLPredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX, color = "Observed")) +
  geom_line(aes(y = Mean, color = "PPD Mean")) +
  facet_wrap(vars(StormID), ncol = 6)+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50)) +
  labs(title = "logNormalFitNULL PPD Mean vs Observed VMAX",
       subtitle = "95% Credible Interval about PPD Mean") +
  scale_color_manual(name = NULL, values = c("black","red")) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1))
  ) +
  theme_bw()
logNormalFitNULLstormsPredplot

##### PPC ----
###### Quantile 2.5 
logNormalFitNULL2LCBsims <- apply(logNormalFitNULL2finalFit, 
                               MARGIN = 1,
                               function(x){
                                 quantile(x, 0.025)
                               })
logNormalFitNULL2LCBpvalueVec <- logNormalFitNULL2LCBsims < quantile(StormdataTrain3$VMAX, 0.025)
logNormalFitNULL2LCBpvalue <- sum(logNormalFitNULL2LCBpvalueVec)
logNormalFitNULL2LCBpvalue <- round(logNormalFitNULL2LCBpvalue/4000, 3)
logNormalFitNULL2LCBpvalue <- min(logNormalFitNULL2LCBpvalue, 1 - logNormalFitNULL2LCBpvalue)

logNormalFitNULL2_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULL2finalFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", logNormalFitNULL2LCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULL2_ppcLCB

###### Quantile 97.5 
logNormalFitNULL2UCBsims <- apply(logNormalFitNULL2finalFit, 
                               MARGIN = 1,
                               function(x){
                                 quantile(x, 0.975)
                               })
logNormalFitNULL2UCBpvalueVec <- logNormalFitNULL2UCBsims < quantile(StormdataTrain3$VMAX, 0.975)
logNormalFitNULL2UCBpvalue <- as.numeric(sum(logNormalFitNULL2UCBpvalueVec))
logNormalFitNULL2UCBpvalue <- round(logNormalFitNULL2UCBpvalue/4000, 3)
logNormalFitNULL2UCBpvalue <- min(logNormalFitNULL2UCBpvalue, 1 - logNormalFitNULL2UCBpvalue)

logNormalFitNULL2_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULL2finalFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", logNormalFitNULL2UCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULL2_ppcUCB

###### Mean 
logNormalFitNULL2MEANsims <- apply(logNormalFitNULL2finalFit, 
                                MARGIN = 1,
                                function(x){
                                  mean(x)
                                })
logNormalFitNULL2MEANpvalueVec <- logNormalFitNULL2MEANsims < mean(StormdataTrain3$VMAX)
logNormalFitNULL2MEANpvalue <- sum(logNormalFitNULL2MEANpvalueVec)
logNormalFitNULL2MEANpvalue <- round(logNormalFitNULL2MEANpvalue/4000, 3)
logNormalFitNULL2MEANpvalue <- min(logNormalFitNULL2MEANpvalue, 1 - logNormalFitNULL2MEANpvalue)

logNormalFitNULL2_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULL2finalFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", logNormalFitNULL2MEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULL2_ppcMEAN

###### Med 
logNormalFitNULL2MEDsims <- apply(logNormalFitNULL2finalFit, 
                               MARGIN = 1,
                               function(x){
                                 quantile(x, 0.5)
                               })
logNormalFitNULL2MEDpvalueVec <- logNormalFitNULL2MEDsims < quantile(StormdataTrain3$VMAX, 0.5)
logNormalFitNULL2MEDpvalue <- sum(logNormalFitNULL2MEDpvalueVec)
logNormalFitNULL2MEDpvalue <- round(logNormalFitNULL2MEDpvalue/4000, 3)
logNormalFitNULL2MEDpvalue <- min(logNormalFitNULL2MEDpvalue, 1 - logNormalFitNULL2MEDpvalue)

logNormalFitNULL2_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULL2finalFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", logNormalFitNULL2MEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULL2_ppcMED

###### SD 
logNormalFitNULL2SDsims <- apply(logNormalFitNULL2finalFit, 
                              MARGIN = 1,
                              function(x){
                                sd(x)
                              })
logNormalFitNULL2SDpvalueVec <- logNormalFitNULL2SDsims < sd(StormdataTrain3$VMAX)
logNormalFitNULL2SDpvalue <- sum(logNormalFitNULL2SDpvalueVec)
logNormalFitNULL2SDpvalue <- round(logNormalFitNULL2SDpvalue/4000, 3)
logNormalFitNULL2SDpvalue <- min(logNormalFitNULL2SDpvalue, 1 - logNormalFitNULL2SDpvalue)

logNormalFitNULL2_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULL2finalFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", logNormalFitNULL2SDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULL2_ppcSD

###### Range 
logNormalFitNULL2RANGEsims <- apply(logNormalFitNULL2finalFit, 
                                 MARGIN = 1,
                                 function(x){
                                   max(x)-min(x)
                                 })
logNormalFitNULL2RANGEpvalueVec <- logNormalFitNULL2RANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
logNormalFitNULL2RANGEpvalue <- sum(logNormalFitNULL2RANGEpvalueVec)
logNormalFitNULL2RANGEpvalue <- round(logNormalFitNULL2RANGEpvalue/4000, 3)
logNormalFitNULL2RANGEpvalue <- min(logNormalFitNULL2RANGEpvalue, 1 - logNormalFitNULL2RANGEpvalue)

logNormalFitNULL2_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULL2finalFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", logNormalFitNULL2RANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULL2_ppcRANGE

##### Combined Plot ----
logNormalFitNULL2_ppcComb <- 
  #logNormalFitNULL2ppcFit /
  (logNormalFitNULL2_ppcLCB | logNormalFitNULL2_ppcMED | logNormalFitNULL2_ppcUCB |
  logNormalFitNULL2_ppcRANGE | logNormalFitNULL2_ppcMEAN | logNormalFitNULL2_ppcSD)
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
logNormalFitNULL2_ppcComb

##### Bayes p-values ----
logNormalFitNULL2pvalues <- tibble(
  Fit = "logNormalFitNULL2",
  LCB = logNormalFitNULL2LCBpvalue,
  Median = logNormalFitNULL2MEDpvalue,
  UCB = logNormalFitNULL2UCBpvalue,
  Range = logNormalFitNULL2RANGEpvalue,
  Mean = logNormalFitNULL2MEANpvalue,
  SD = logNormalFitNULL2SDpvalue
)
logNormalFitNULL2pvalues

##### CV ----
logNormalFitNULL2kfold <- kfold(logNormalFitNULL2,
                         save_fits = TRUE,
                         chains = 1, 
                         K = 5)
logNormalFitNULL2kfoldPreds <- kfold_predict(logNormalFitNULL2kfold)
logNormalFitNULL2kfoldPredsDat <- logNormalFitNULL2kfoldPreds$yrep
logNormalFitNULL2kfoldPredsMean <- colMeans(logNormalFitNULL2kfoldPredsDat)
logNormalFitNULL2kfoldPredsMed <- apply(logNormalFitNULL2kfoldPredsDat, 2, function(x){quantile(x, 0.5)})
logNormalFitNULL2kfoldPredsLCB <- apply(logNormalFitNULL2kfoldPredsDat, 2, function(x){quantile(x, 0.025)})
logNormalFitNULL2kfoldPredsUCB <- apply(logNormalFitNULL2kfoldPredsDat, 2, function(x){quantile(x, 0.975)})

logNormalFitNULL2kfoldMetrics <- tibble(
  Fit = "logNormalFitNULL2",
  MAE_kfold = mean(abs(logNormalFitNULL2kfoldPredsMean - logNormalFitNULL2kfoldPreds$y)),
  MAD_kfold = mean(abs(logNormalFitNULL2kfoldPredsMed - logNormalFitNULL2kfoldPreds$y)),
  COV_kfold = mean(logNormalFitNULL2kfoldPredsLCB < logNormalFitNULL2kfoldPreds$y & logNormalFitNULL2kfoldPreds$y < logNormalFitNULL2kfoldPredsUCB)
)
logNormalFitNULL2kfoldMetrics


logNormalFitNULLfiles <- ls()[str_detect(ls(), pattern = "logNormalFitNULL")]
logNormalFitNULLfilesRM <- logNormalFitNULLfiles[!(logNormalFitNULLfiles %in% c("logNormalFitNULL_ppcComb",
                                                                          "logNormalFitNULLloo",
                                                                          "logNormalFitNULLpredMetrics",
                                                                          "logNormalFitNULLpvalues",
                                                                          "logNormalFitNULLstormsPredplot"))]

rm(list = logNormalFitNULLfilesRM)
rm(logNormalFitNULLfilesRM)

### GAMMA ----
gammaFitNULL <- brm(
  bf(
    VMAX ~ offset(HWRF), family = Gamma(link = "inverse")
    # VMAX ~ log(HWRF), family = lognormal(),
    # VMAX ~ log(HWRF), family = lognormal(),
    # VMAX ~ log(HWRF), family = lognormal(),
  ),
  data = StormdataTrain8,
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)

gammaFitNULL2 <- brm(
  bf(
    VMAX ~ offset(log(HWRF)) + (1|StormID), family = Gamma(link = "log")
    # VMAX ~ log(HWRF), family = lognormal(),
    # VMAX ~ log(HWRF), family = lognormal(),
    # VMAX ~ log(HWRF), family = lognormal(),
  ),
  data = StormdataTrain8, 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)

#save(gammaFitNULL, file = "_data/gammaFitNULL.RData")
prior_summary(gammaFitNULL)
round(posterior_summary(gammaFitNULL, probs = c(0.025, 0.975)))
gammaFitNULL

print(gammaFitNULL, digits = 4)
plot(gammaFitNULL)
gammaFitNULLppcFit <- pp_check(gammaFitNULL, ndraws = 100) + 
  labs(title = "gammaFitNULL Fit PPC") +
  theme_bw()
gammaFitNULLppcFit
gammaFitNULLloo <- loo(gammaFitNULL)
waic(gammaFitNULL)
performance::check_distribution(gammaFitNULL)
performance::check_outliers(gammaFitNULL)
performance::check_heteroskedasticity(gammaFitNULL)
performance_rmse(gammaFitNULL)
performance_mae(gammaFitNULL)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(gammaFitNULL)


print(gammaFitNULL2, digits = 4)
plot(gammaFitNULL2)
gammaFitNULL2ppcFit <- pp_check(gammaFitNULL2, ndraws = 100) + 
  labs(title = "gammaFitNULL2 Fit PPC") +
  theme_bw()
gammaFitNULL2ppcFit
gammaFitNULL2loo <- loo(gammaFitNULL2)
waic(gammaFitNULL2)
performance::check_distribution(gammaFitNULL2)
performance::check_outliers(gammaFitNULL2)
performance::check_heteroskedasticity(gammaFitNULL2)
performance_rmse(gammaFitNULL2)
performance_mae(gammaFitNULL2)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(gammaFitNULL2)

variance_decomposition(gammaFitNULL)
exp(fixef(gammaFitNULL))
ranef(gammaFitNULL)

bayes_R2(gammaFitNULL)

bayes_factor(gammaFitNULL2, gammaFitNULL)
bayes_factor(gammaFitNULL2, logNormalFitNULL2)
bayes_factor(logNormalFitNULL2, gammaFitNULL2)
bayes_factor(gammaFitNULL, studentFit1)
bayes_factor(gammaFitNULL, studentFit2)
bayes_factor(gammaFitNULL, studentFit3)
bayes_factor(gammaFitNULL, linFit11)
bayes_factor(gammaFitNULL, propFit1)
bayes_factor(gammaFitNULL, logPropFit1)
loo(gammaFitNULL, gammaFit3)

gammaFitNULLsmooths <- conditional_smooths(gammaFitNULL)
plot(gammaFitNULLsmooths, stype = "raster", ask = FALSE)
gammaFitNULLeffects <- conditional_effects(gammaFitNULL, 
                                               method = "posterior_predict",
                                               robust = FALSE,
                                               re_formula = NULL)
gammaFitNULLeffects <- conditional_effects(gammaFitNULL)
plot(gammaFitNULLeffects, points = TRUE, ask = FALSE)

##### Prediction ----
## Fitted
gammaFitNULLfinalFit <- posterior_predict(gammaFitNULL)
gammaFitNULLfinalFitMean <- colMeans(gammaFitNULLfinalFit)
gammaFitNULLfinalFitMed <- apply(gammaFitNULLfinalFit, 2, function(x){quantile(x, 0.5)})
gammaFitNULLfinalFitLCB <- apply(gammaFitNULLfinalFit, 2, function(x){quantile(x, 0.025)})
gammaFitNULLfinalFitUCB <- apply(gammaFitNULLfinalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
gammaFitNULLfinalPreds <- posterior_predict(gammaFitNULL, 
                                                newdata = StormdataTest8,
                                                allow_new_levels = TRUE)
gammaFitNULLfinalPredsMean <- colMeans(gammaFitNULLfinalPreds)
gammaFitNULLfinalPredsMed <- apply(gammaFitNULLfinalPreds, 2, function(x){quantile(x, 0.5)})
gammaFitNULLfinalPredsLCB <- apply(gammaFitNULLfinalPreds, 2, function(x){quantile(x, 0.025)})
gammaFitNULLfinalPredsUCB <- apply(gammaFitNULLfinalPreds, 2, function(x){quantile(x, 0.975)})

gammaFitNULLpredMetrics <- tibble(
  Fit = "gammaFitNULL",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(gammaFitNULLfinalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(gammaFitNULLfinalFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < gammaFitNULLfinalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(gammaFitNULLfinalPredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(gammaFitNULLfinalPredsMed - Actual_Yvec)),
  COV_pred = mean(gammaFitNULLfinalPredsLCB < Actual_Yvec & Actual_Yvec < gammaFitNULLfinalPredsUCB)
)
gammaFitNULLpredMetrics

## Fitted
gammaFitNULL2finalFit <- posterior_predict(gammaFitNULL2)
gammaFitNULL2finalFitMean <- colMeans(gammaFitNULL2finalFit)
gammaFitNULL2finalFitMed <- apply(gammaFitNULL2finalFit, 2, function(x){quantile(x, 0.5)})
gammaFitNULL2finalFitLCB <- apply(gammaFitNULL2finalFit, 2, function(x){quantile(x, 0.025)})
gammaFitNULL2finalFitUCB <- apply(gammaFitNULL2finalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
gammaFitNULL2finalPreds <- posterior_predict(gammaFitNULL2, 
                                            newdata = StormdataTest8,
                                            allow_new_levels = TRUE)
gammaFitNULL2finalPredsMean <- colMeans(gammaFitNULL2finalPreds)
gammaFitNULL2finalPredsMed <- apply(gammaFitNULL2finalPreds, 2, function(x){quantile(x, 0.5)})
gammaFitNULL2finalPredsLCB <- apply(gammaFitNULL2finalPreds, 2, function(x){quantile(x, 0.025)})
gammaFitNULL2finalPredsUCB <- apply(gammaFitNULL2finalPreds, 2, function(x){quantile(x, 0.975)})

gammaFitNULL2predMetrics <- tibble(
  Fit = "gammaFitNULL2",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(gammaFitNULL2finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(gammaFitNULL2finalFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < gammaFitNULL2finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(gammaFitNULL2finalPredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(gammaFitNULL2finalPredsMed - Actual_Yvec)),
  COV_pred = mean(gammaFitNULL2finalPredsLCB < Actual_Yvec & Actual_Yvec < gammaFitNULL2finalPredsUCB)
)
gammaFitNULL2predMetrics

##### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = gammaFitNULLfinalPreds) +
  labs(title = "gammaFitNULL Predict") +
  theme_bw()

gammaFitNULLFitDF <- bind_cols(
  StormdataTrain3,
  LCB = gammaFitNULLfinalFitLCB,
  Mean = gammaFitNULLfinalFitMean,
  Med = gammaFitNULLfinalFitMed,
  UCB = gammaFitNULLfinalFitUCB
)

gammaFitNULLstormsFitplot <- ggplot(data = gammaFitNULLFitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))
gammaFitNULLstormsFitplot

## Prediction
gammaFitNULLPredDF <- bind_cols(
  StormdataTest3,
  LCB = gammaFitNULLfinalPredsLCB,
  Mean = gammaFitNULLfinalPredsMean,
  Med = gammaFitNULLfinalPredsMed,
  UCB = gammaFitNULLfinalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

gammaFitNULLstormsPredplot <- ggplot(data = gammaFitNULLPredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX, color = "Observed")) +
  geom_line(aes(y = Mean, color = "PPD Mean")) +
  facet_wrap(vars(StormID), ncol = 6)+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50)) +
  labs(title = "gammaFitNULL PPD Mean vs Observed VMAX",
       subtitle = "95% Credible Interval about PPD Mean") +
  scale_color_manual(name = NULL, values = c("black","red")) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1))
  ) +
  theme_bw()
gammaFitNULLstormsPredplot

##### PPC ----
###### Quantile 2.5 
gammaFitNULL2LCBsims <- apply(gammaFitNULL2finalFit, 
                                 MARGIN = 1,
                                 function(x){
                                   quantile(x, 0.025)
                                 })
gammaFitNULL2LCBpvalueVec <- gammaFitNULL2LCBsims < quantile(StormdataTrain3$VMAX, 0.025)
gammaFitNULL2LCBpvalue <- sum(gammaFitNULL2LCBpvalueVec)
gammaFitNULL2LCBpvalue <- round(gammaFitNULL2LCBpvalue/4000, 3)
gammaFitNULL2LCBpvalue <- min(gammaFitNULL2LCBpvalue, 1 - gammaFitNULL2LCBpvalue)

gammaFitNULL2_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULL2finalFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", gammaFitNULL2LCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULL2_ppcLCB

###### Quantile 97.5 
gammaFitNULL2UCBsims <- apply(gammaFitNULL2finalFit, 
                                 MARGIN = 1,
                                 function(x){
                                   quantile(x, 0.975)
                                 })
gammaFitNULL2UCBpvalueVec <- gammaFitNULL2UCBsims < quantile(StormdataTrain3$VMAX, 0.975)
gammaFitNULL2UCBpvalue <- as.numeric(sum(gammaFitNULL2UCBpvalueVec))
gammaFitNULL2UCBpvalue <- round(gammaFitNULL2UCBpvalue/4000, 3)
gammaFitNULL2UCBpvalue <- min(gammaFitNULL2UCBpvalue, 1 - gammaFitNULL2UCBpvalue)

gammaFitNULL2_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULL2finalFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", gammaFitNULL2UCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULL2_ppcUCB

###### Mean 
gammaFitNULL2MEANsims <- apply(gammaFitNULL2finalFit, 
                                  MARGIN = 1,
                                  function(x){
                                    mean(x)
                                  })
gammaFitNULL2MEANpvalueVec <- gammaFitNULL2MEANsims < mean(StormdataTrain3$VMAX)
gammaFitNULL2MEANpvalue <- sum(gammaFitNULL2MEANpvalueVec)
gammaFitNULL2MEANpvalue <- round(gammaFitNULL2MEANpvalue/4000, 3)
gammaFitNULL2MEANpvalue <- min(gammaFitNULL2MEANpvalue, 1 - gammaFitNULL2MEANpvalue)

gammaFitNULL2_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULL2finalFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", gammaFitNULL2MEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULL2_ppcMEAN

###### Med 
gammaFitNULL2MEDsims <- apply(gammaFitNULL2finalFit, 
                                 MARGIN = 1,
                                 function(x){
                                   quantile(x, 0.5)
                                 })
gammaFitNULL2MEDpvalueVec <- gammaFitNULL2MEDsims < quantile(StormdataTrain3$VMAX, 0.5)
gammaFitNULL2MEDpvalue <- sum(gammaFitNULL2MEDpvalueVec)
gammaFitNULL2MEDpvalue <- round(gammaFitNULL2MEDpvalue/4000, 3)
gammaFitNULL2MEDpvalue <- min(gammaFitNULL2MEDpvalue, 1 - gammaFitNULL2MEDpvalue)

gammaFitNULL2_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULL2finalFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", gammaFitNULL2MEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULL2_ppcMED

###### SD 
gammaFitNULL2SDsims <- apply(gammaFitNULL2finalFit, 
                                MARGIN = 1,
                                function(x){
                                  sd(x)
                                })
gammaFitNULL2SDpvalueVec <- gammaFitNULL2SDsims < sd(StormdataTrain3$VMAX)
gammaFitNULL2SDpvalue <- sum(gammaFitNULL2SDpvalueVec)
gammaFitNULL2SDpvalue <- round(gammaFitNULL2SDpvalue/4000, 3)
gammaFitNULL2SDpvalue <- min(gammaFitNULL2SDpvalue, 1 - gammaFitNULL2SDpvalue)

gammaFitNULL2_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULL2finalFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", gammaFitNULL2SDpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULL2_ppcSD

###### Range 
gammaFitNULL2RANGEsims <- apply(gammaFitNULL2finalFit, 
                                   MARGIN = 1,
                                   function(x){
                                     max(x)-min(x)
                                   })
gammaFitNULL2RANGEpvalueVec <- gammaFitNULL2RANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
gammaFitNULL2RANGEpvalue <- sum(gammaFitNULL2RANGEpvalueVec)
gammaFitNULL2RANGEpvalue <- round(gammaFitNULL2RANGEpvalue/4000, 3)
gammaFitNULL2RANGEpvalue <- min(gammaFitNULL2RANGEpvalue, 1 - gammaFitNULL2RANGEpvalue)

gammaFitNULL2_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULL2finalFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", gammaFitNULL2RANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULL2_ppcRANGE

##### Combined Plot ----
gammaFitNULL2_ppcComb <- 
  #gammaFitNULL2ppcFit /
  (gammaFitNULL2_ppcLCB | gammaFitNULL2_ppcMED | gammaFitNULL2_ppcUCB |
  gammaFitNULL2_ppcRANGE | gammaFitNULL2_ppcMEAN | gammaFitNULL2_ppcSD)
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
gammaFitNULL2_ppcComb

##### Bayes p-values ----
gammaFitNULL2pvalues <- tibble(
  Fit = "gammaFitNULL2",
  LCB = gammaFitNULL2LCBpvalue,
  Median = gammaFitNULL2MEDpvalue,
  UCB = gammaFitNULL2UCBpvalue,
  Range = gammaFitNULL2RANGEpvalue,
  Mean = gammaFitNULL2MEANpvalue,
  SD = gammaFitNULL2SDpvalue
)
gammaFitNULL2pvalues

##### CV ----
gammaFitNULL2kfold <- kfold(gammaFitNULL2,
                                save_fits = TRUE,
                                chains = 1, 
                                K = 5)
gammaFitNULL2kfoldPreds <- kfold_predict(gammaFitNULL2kfold)
gammaFitNULL2kfoldPredsDat <- gammaFitNULL2kfoldPreds$yrep
gammaFitNULL2kfoldPredsMean <- colMeans(gammaFitNULL2kfoldPredsDat)
gammaFitNULL2kfoldPredsMed <- apply(gammaFitNULL2kfoldPredsDat, 2, function(x){quantile(x, 0.5)})
gammaFitNULL2kfoldPredsLCB <- apply(gammaFitNULL2kfoldPredsDat, 2, function(x){quantile(x, 0.025)})
gammaFitNULL2kfoldPredsUCB <- apply(gammaFitNULL2kfoldPredsDat, 2, function(x){quantile(x, 0.975)})

gammaFitNULL2kfoldMetrics <- tibble(
  Fit = "gammaFitNULL2",
  MAE_kfold = mean(abs(gammaFitNULL2kfoldPredsMean - gammaFitNULL2kfoldPreds$y)),
  MAD_kfold = mean(abs(gammaFitNULL2kfoldPredsMed - gammaFitNULL2kfoldPreds$y)),
  COV_kfold = mean(gammaFitNULL2kfoldPredsLCB < gammaFitNULL2kfoldPreds$y & gammaFitNULL2kfoldPreds$y < gammaFitNULL2kfoldPredsUCB)
)
gammaFitNULL2kfoldMetrics

gammaFitNULLfiles <- ls()[str_detect(ls(), pattern = "gammaFitNULL")]
gammaFitNULLfilesRM <- gammaFitNULLfiles[!(gammaFitNULLfiles %in% c("gammaFitNULL_ppcComb",
                                                                                "gammaFitNULLloo",
                                                                                "gammaFitNULLpredMetrics",
                                                                                "gammaFitNULLpvalues",
                                                                                "gammaFitNULLstormsPredplot"))]

rm(list = gammaFitNULLfilesRM)
rm(gammaFitNULLfilesRM)


## Compare NULL Models ----
### PPC ----
ppcCompPlot <- (LinFitNULLppcFit + scale_x_continuous(limits = c(-300, 400))) /
(LinFitNULL2ppcFit  + scale_x_continuous(limits = c(-300, 400))) /
(skewLinFitNULLppcFit + scale_x_continuous(limits = c(-300, 400))) /
(skewLinFitNULL2ppcFit  + scale_x_continuous(limits = c(-300, 400))) /
(StudentFitNULLppcFit + scale_x_continuous(limits = c(-300, 400))) /
(StudentFitNULL2ppcFit + scale_x_continuous(limits = c(-300, 400))) /
(logNormalFitNULL2ppcFit + scale_x_continuous(limits = c(-300, 400))) /
(gammaFitNULL2ppcFit + scale_x_continuous(limits = c(-300, 400))) +
  plot_layout(guides = "collect")

### Bayes pvalue ----
LinFitNULL_ppcComb /
LinFitNULL2_ppcComb /
skewLinFitNULL_ppcComb /
skewLinFitNULL2_ppcComb /
StudentFitNULL_ppcComb /
StudentFitNULL2_ppcComb /
logNormalFitNULL2_ppcComb /
gammaFitNULL2_ppcComb 

pvalueCompNULL <- bind_rows(
  LinFitNULLpvalues, 
  LinFitNULL2pvalues,
  skewLinFitNULLpvalues,
  skewLinFitNULL2pvalues,
  StudentFitNULLpvalues, 
  StudentFitNULL2pvalues,
  logNormalFitNULL2pvalues,
  gammaFitNULL2pvalues
)
pvalueCompNULL
save(pvalueCompNULL, file = "~/Desktop/Temp Hurricane Model Data/pvalueCompNULL.RData")

### LOO ----
looCompNULL <- loo_compare(LinFitNULLloo, LinFitNULL2loo,
                           skewLinFitNULLloo, skewLinFitNULL2loo,
                           StudentFitNULLloo, StudentFitNULL2loo,
                           logNormalFitNULL2loo,
                           gammaFitNULL2loo)
looCompNULL
save(looCompNULL, file = "~/Desktop/Temp Hurricane Model Data/looCompNULL.RData")

### CV ----
cvCompNULL <- bind_rows(
  LinFitNULLkfoldMetrics, LinFitNULL2kfoldMetrics,
  skewLinFitNULLkfoldMetrics, skewLinFitNULL2kfoldMetrics,
  StudentFitNULLkfoldMetrics, StudentFitNULL2kfoldMetrics,
  logNormalFitNULL2kfoldMetrics,
  gammaFitNULL2kfoldMetrics
)
cvCompNULL2 <- cvCompNULL |> arrange(MAE_kfold)
save(cvCompNULL2, file = "~/Desktop/Temp Hurricane Model Data/cvCompNULL2.RData")


### Preds ----
predCompNULL <- bind_rows(
  LinFitNULLpredMetrics, LinFitNULL2predMetrics,
  skewLinFitNULLpredMetrics, skewLinFitNULL2predMetrics,
  StudentFitNULLpredMetrics, StudentFitNULL2predMetrics,
  logNormalFitNULL2predMetrics,
  gammaFitNULL2predMetrics
)
predCompNULL2 <- predCompNULL |> arrange(MAE_pred)
save(predCompNULL2, file = "~/Desktop/Temp Hurricane Model Data/predCompNULL2.RData")


### Bayes R2
bayes_R2(LinFitNULL)
bayes_R2(LinFitNULL2)
bayes_R2(skewLinFitNULL)
bayes_R2(skewLinFitNULL2)
bayes_R2(StudentFitNULL)
bayes_R2(StudentFitNULL2)
bayes_R2(logNormalFitNULL2)
bayes_R2(gammaFitNULL2)

bayes_factor(LinFitNULL, LinFitNULL2)
bayes_factor(skewLinFitNULL, skewLinFitNULL2)
bayes_factor(StudentFitNULL, StudentFitNULL2)
bayes_factor(LinFitNULL, skewLinFitNULL)
bayes_factor(skewLinFitNULL, LinFitNULL)
bayes_factor(skewLinFitNULL, StudentFitNULL)
bayes_factor(StudentFitNULL, skewLinFitNULL)
bayes_factor(StudentFitNULL, logNormalFitNULL2)
bayes_factor(logNormalFitNULL2, StudentFitNULL)
bayes_factor(logNormalFitNULL2, gammaFitNULL2)

save(LinFitNULL, file = "~/Desktop/Temp Hurricane Model Data/LinFitNULL.RData")
save(LinFitNULL2, file = "~/Desktop/Temp Hurricane Model Data/LinFitNULL2.RData")
save(skewLinFitNULL, file = "~/Desktop/Temp Hurricane Model Data/skewLinFitNULL.RData")
save(skewLinFitNULL2, file = "~/Desktop/Temp Hurricane Model Data/skewLinFitNULL2.RData")
save(StudentFitNULL, file = "~/Desktop/Temp Hurricane Model Data/StudentFitNULL.RData")
save(StudentFitNULL2, file = "~/Desktop/Temp Hurricane Model Data/StudentFitNULL2.RData")
save(logNormalFitNULL2, file = "~/Desktop/Temp Hurricane Model Data/logNormalFitNULL2.RData")
save(gammaFitNULL2, file = "~/Desktop/Temp Hurricane Model Data/gammaFitNULL2.RData")

