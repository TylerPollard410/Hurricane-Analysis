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
Stormdata_raw <- fread("~/Desktop/Hurricane Analysis/_data/E2_data.csv")
Stormdata <- Stormdata_raw |>
  mutate(
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
    across(where(is.numeric) & !c(VMAX, Day,
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
plot(x = StormdataTrain7scale$Day, y = StormdataTrain7scale$VMAX)
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
StormdataTestFinal <- StormdataTest3 |>
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
str(StormdataTestFinal)

#### Test Final Scale----
StormdataTestFinalscale <- StormdataTest3 |>
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
str(StormdataTestFinalscale)

#### Test Final Scale2 ----
StormdataTestFinalscale2 <- StormdataTest3 |>
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
str(StormdataTestFinalscale2)

### Actual ----
Actual_Y <- fread("_data/Actual Y.csv")
Actual_Yvec <- Actual_Y |> filter(complete.cases(x)) |> pull(x)


# Plot VMAX ----
ggpairs(StormdataTrain8 |>
          select(-StormID, -Date, -Month, -Year))

## Scatter ----
ggplot(data = StormdataTrain3, aes(x = HWFI, y = VMAX)) +
  geom_point() +
  geom_smooth()


## Histogram ----
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
range(StormdataTestFinal$StormElapsedTime)
range(StormdataTestFinal$LAT)
which.min(StormdataTestFinal$LAT)
range(StormdataTestFinal$LON)
which.min(StormdataTestFinal$LON)

world_coordinates <- map_data("world") 
ggplot() + 
  # geom_map() function takes world coordinates  
  # as input to plot world map 
  geom_map( 
    data = world_coordinates, map = world_coordinates, 
    aes(map_id = region) 
  ) + 
  geom_point(
    data = StormdataTestFinal,
    aes(x = LON-360, y = LAT, 
        color = StormElapsedTime)
  ) +
  xlim(c(-190,0)) +
  ylim(c(0,60)) +
  scale_color_continuous(low = "green", high = "red") +
  theme_bw()

# Fit model ----
load(file = "_data/gammaFit0.RData")
load(file = "_data/gammaFit1.RData")
load(file = "_data/gammaFit2.RData")
load(file = "_data/gammaFit3.RData")
load(file = "_data/gammaFit4.RData")
load(file = "_data/gammaFit5.RData")

## VMAX ----
### Model 1 ----
gammaFit5 <- brm(
  formula = VMAX ~
    #Year +
    #Month +
    basin + 
    s(Day) +
    s(StormElapsedTime) + 
    t2(LON, LAT) +
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
  data = StormdataTrain7scale, 
  family = brmsfamily("Gamma", link = "log"), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)
save(gammaFit5, file = "_data/gammaFit5.RData")
prior_summary(gammaFit5)
round(posterior_summary(gammaFit5, probs = c(0.025, 0.975)))
gammaFit5

print(gammaFit5, digits = 4)
plot(gammaFit5)
pp_check(gammaFit5, ndraws = 100)
loo(gammaFit5)
waic(gammaFit5)
performance::check_distribution(gammaFit5)
performance::check_outliers(gammaFit5)
performance::check_heteroskedasticity(gammaFit5)
performance_rmse(gammaFit5)
performance_mae(gammaFit5)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(gammaFit5)


variance_decomposition(gammaFit5)
exp(fixef(gammaFit5))
ranef(gammaFit5)

bayes_R2(gammaFit5)

bayes_factor(gammaFit5, gammaFit1)
bayes_factor(gammaFit5, gammaFit2)
bayes_factor(gammaFit5, gammaFit3)
bayes_factor(gammaFit5, studentFit1)
bayes_factor(gammaFit5, studentFit2)
bayes_factor(gammaFit5, studentFit3)
bayes_factor(gammaFit5, linFit11)
bayes_factor(gammaFit5, propFit1)
bayes_factor(gammaFit5, logPropFit1)
loo(gammaFit5, gammaFit3)

gammaFit5smooths <- conditional_smooths(gammaFit5)
plot(gammaFit5smooths, stype = "raster", ask = FALSE)
gammaFit5effects <- conditional_effects(gammaFit5, 
                                        method = "posterior_predict",
                                        robust = FALSE,
                                        re_formula = NULL)
gammaFit5effects <- conditional_effects(gammaFit5)
plot(gammaFit5effects, points = TRUE, ask = FALSE)

#### Prediction ----
## Fitted
gammaFit5finalFit <- posterior_predict(gammaFit5)
gammaFit5finalFitMean <- colMeans(gammaFit5finalFit)
gammaFit5finalFitMed <- apply(gammaFit5finalFit, 2, function(x){quantile(x, 0.5)})
gammaFit5finalFitLCB <- apply(gammaFit5finalFit, 2, function(x){quantile(x, 0.025)})
gammaFit5finalFitUCB <- apply(gammaFit5finalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
gammaFit5finalPreds <- posterior_predict(gammaFit5, 
                                         newdata = StormdataTestFinalscale,
                                         allow_new_levels = TRUE)
gammaFit5finalPreds2 <- colMeans(gammaFit5finalPreds)
gammaFit5finalPredsMed <- apply(gammaFit5finalPreds, 2, function(x){quantile(x, 0.5)})
gammaFit5finalPredsLCB <- apply(gammaFit5finalPreds, 2, function(x){quantile(x, 0.025)})
gammaFit5finalPredsUCB <- apply(gammaFit5finalPreds, 2, function(x){quantile(x, 0.975)})

gammaFit5predMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(gammaFit5finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(gammaFit5finalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < gammaFit5finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(gammaFit5finalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(gammaFit5finalPredsMed - Actual_Yvec)),
  COV_pred = mean(gammaFit5finalPredsLCB < Actual_Yvec & Actual_Yvec < gammaFit5finalPredsUCB)
)
gammaFit5predMetrics

#### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = gammaFit5finalPreds) +
  labs(title = "GammaFit5 Predict")

gammaFit5FitDF <- bind_cols(
  StormdataTrain3,
  LCB = gammaFit5finalFitLCB,
  Mean = gammaFit5finalFitMean,
  Med = gammaFit5finalFitMed,
  UCB = gammaFit5finalFitUCB
) |>
  filter(StormID %in% c(5712015))

ggplot(data = gammaFit5FitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))

## Prediction
gammaFit5PredDF <- bind_cols(
  StormdataTest3,
  LCB = gammaFit5finalPredsLCB,
  Mean = gammaFit5finalPreds2,
  Med = gammaFit5finalPredsMed,
  UCB = gammaFit5finalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

ggplot(data = gammaFit5PredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))

rm(gammaFit5finalFit)
rm(gammaFit5finalPreds)

## VMAX/HWRF ----
### LINEAR ----
#### Model 1 ----
propLinFit1 <- brm(
  bf(
    VMAX/HWRF ~
      #Year +
      Month +
      basin + 
      LON +
      LAT +
      #s(Day) +
      #s(StormElapsedTime) + 
      #t2(LON, LAT) +
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
      #HWRF +
      (1|StormID)
  ),
  data = StormdataTrain8, 
  family = gaussian(link = "log"), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)
save(propLinFit1, file = "_data/propLinFit1.RData")
prior_summary(propLinFit1)
round(posterior_summary(propLinFit1, probs = c(0.025, 0.975)))
propLinFit1

print(propLinFit1, digits = 4)
plot(propLinFit1)
pp_check(propLinFit1, ndraws = 100) + labs(title = "propLinFit1")
loo(propLinFit1)
waic(propLinFit1)
performance::check_distribution(propLinFit1)
performance::check_outliers(propLinFit1)
performance::check_heteroskedasticity(propLinFit1)
performance_rmse(propLinFit1)
performance_mae(propLinFit1)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(propLinFit1)


variance_decomposition(propLinFit1)
exp(fixef(propLinFit1))
ranef(propLinFit1)

bayes_R2(propLinFit1)

bayes_factor(propLinFit1, propLinFit1)
bayes_factor(propLinFit1, gammaFit1)
bayes_factor(propLinFit1, gammaFit1)
bayes_factor(propLinFit1, studentFit1)
bayes_factor(propLinFit1, studentFit1)
bayes_factor(propLinFit1, studentFit1)
bayes_factor(propLinFit1, linFit11)
bayes_factor(propLinFit1, propFit1)
bayes_factor(propLinFit1, logPropFit1)
loo(propLinFit1, gammaFit1)

propLinFit1smooths <- conditional_smooths(propLinFit1)
plot(propLinFit1smooths, stype = "raster", ask = FALSE)
propLinFit1effects <- conditional_effects(propLinFit1, 
                                          method = "posterior_predict",
                                          robust = FALSE,
                                          re_formula = NULL)
propLinFit1effects <- conditional_effects(propLinFit1)
plot(propLinFit1effects, points = TRUE, ask = FALSE)

##### Prediction ----
## Fitted
propLinFit1finalFit <- posterior_predict(propLinFit1)
propLinFit1finalFit2 <- t(t(propLinFit1finalFit)*StormdataTrain3$HWRF)
#propLinFit1finalFitMean <- colMeans(propLinFit1finalFit)*StormdataTrain3$HWRF
propLinFit1finalFitMean <- colMeans(propLinFit1finalFit2)
propLinFit1finalFitMed <- apply(propLinFit1finalFit2, 2, function(x){quantile(x, 0.5)})
propLinFit1finalFitLCB <- apply(propLinFit1finalFit2, 2, function(x){quantile(x, 0.025)})
propLinFit1finalFitUCB <- apply(propLinFit1finalFit2, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
propLinFit1finalPreds <- posterior_predict(propLinFit1, 
                                           newdata = StormdataTestFinalscale2,
                                           allow_new_levels = TRUE)
propLinFit1finalPreds2 <- t(t(propLinFit1finalPreds)*StormdataTest3$HWRF)
propLinFit1finalPreds3 <- colMeans(propLinFit1finalPreds2)
propLinFit1finalPredsMed <- apply(propLinFit1finalPreds2, 2, function(x){quantile(x, 0.5)})
propLinFit1finalPredsLCB <- apply(propLinFit1finalPreds2, 2, function(x){quantile(x, 0.025)})
propLinFit1finalPredsUCB <- apply(propLinFit1finalPreds2, 2, function(x){quantile(x, 0.975)})

propLinFit1predMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(propLinFit1finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(propLinFit1finalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < propLinFit1finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(propLinFit1finalPreds3 - Actual_Yvec)),
  MAD_pred = mean(abs(propLinFit1finalPredsMed - Actual_Yvec)),
  COV_pred = mean(propLinFit1finalPredsLCB < Actual_Yvec & Actual_Yvec < propLinFit1finalPredsUCB)
)
propLinFit1predMetrics

##### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = propLinFit1finalPreds) +
  labs(title = "propFit1 Predict")

propLinFit1FitDF <- bind_cols(
  StormdataTrain3,
  LCB = propLinFit1finalFitLCB,
  Mean = propLinFit1finalFitMean,
  Med = propLinFit1finalFitMed,
  UCB = propLinFit1finalFitUCB
) 

propLinFit1ppcDraws <- as_draws_matrix(propLinFit1)

propLinFit1basePPCplot <- ggplot(data = propLinFit1FitDF) +
  # geom_histogram(
  #   aes(x = VMAX/HWRF, after_stat(density)),
  #   color = "#99c7c7", fill = "#bcdcdc",
  #   bins = 100) +
  geom_density(#data = final_data3,
    aes(x = VMAX/HWRF),
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
propLinFit1basePPCplot

set.seed(52)
propLinFit1PPCindex <- sample(1:1705, 100)
for(i in propLinFit1PPCindex){
  propLinFit1basePPCplot <- propLinFit1basePPCplot +
    geom_density(
      aes(x = propLinFit1finalFit[i,])
    )
}
propLinFit1basePPCplot

ggplot(data = propLinFit1FitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))

## Prediction
propLinFit1PredDF <- bind_cols(
  StormdataTest3,
  LCB = propLinFit1finalPredsLCB,
  Mean = propLinFit1finalPreds2,
  Med = propLinFit1finalPredsMed,
  UCB = propLinFit1finalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

ggplot(data = propLinFit1PredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))

rm(propLinFit1finalFit)
rm(propLinFit1finalFit2)
rm(propLinFit1finalPreds)
rm(propLinFit1finalPreds2)

#### Model 2 ----
propLinFit2 <- brm(
  bf(
    VMAX/HWRF ~
      #Year +
      #Month +
      basin + 
      LON +
      LAT +
      Day +
      StormElapsedTime + 
      t2(LON, LAT) +
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
      #HWRF +
      (1|StormID)
  ),
  data = StormdataTrain8, 
  family = gaussian(link = "log"), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)
save(propLinFit3B, file = "_data/propLinFit3B.RData")
prior_summary(propLinFit3B)
round(posterior_summary(propLinFit3B, probs = c(0.025, 0.975)))
propLinFit3B

print(propLinFit3B, digits = 4)
plot(propLinFit3B)
pp_check(propLinFit3B, ndraws = 100)
loo(propLinFit3B)
waic(propLinFit3B)
performance::check_distribution(propLinFit3B)
performance::check_outliers(propLinFit3B)
performance::check_heteroskedasticity(propLinFit3B)
performance_rmse(propLinFit3B)
performance_mae(propLinFit3B)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(propLinFit3B)


variance_decomposition(propLinFit3B)
exp(fixef(propLinFit3B))
ranef(propLinFit3B)

bayes_R2(propLinFit3B)

bayes_factor(propLinFit3B, propLinFit1)
bayes_factor(propLinFit3B, gammaFit3B)
bayes_factor(propLinFit3B, gammaFit3B)
bayes_factor(propLinFit3B, studentFit3B)
bayes_factor(propLinFit3B, studentFit3B)
bayes_factor(propLinFit3B, studentFit3B)
bayes_factor(propLinFit3B, linFit3B1)
bayes_factor(propLinFit3B, propFit3B)
bayes_factor(propLinFit3B, logPropFit3B)
loo(propLinFit3B, gammaFit3B)

propLinFit3Bsmooths <- conditional_smooths(propLinFit3B)
plot(propLinFit3Bsmooths, stype = "raster", ask = FALSE)
propLinFit3Beffects <- conditional_effects(propLinFit3B, 
                                           method = "posterior_predict",
                                           robust = FALSE,
                                           re_formula = NULL)
propLinFit3Beffects <- conditional_effects(propLinFit3B)
plot(propLinFit3Beffects, points = TRUE, ask = FALSE)

##### Prediction ----
## Fitted
propLinFit3BfinalFit <- posterior_predict(propLinFit3B)
propLinFit3BfinalFit <- t(t(propLinFit3BfinalFit)*StormdataTrain3$HWRF)
#propLinFit3BfinalFitMean <- colMeans(propLinFit3BfinalFit)*StormdataTrain3$HWRF
propLinFit3BfinalFitMean <- colMeans(propLinFit3BfinalFit)
propLinFit3BfinalFitMed <- apply(propLinFit3BfinalFit, 2, function(x){quantile(x, 0.5)})
propLinFit3BfinalFitLCB <- apply(propLinFit3BfinalFit, 2, function(x){quantile(x, 0.025)})
propLinFit3BfinalFitUCB <- apply(propLinFit3BfinalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
propLinFit3BfinalPreds <- posterior_predict(propLinFit3B, 
                                            newdata = StormdataTestFinalscale2,
                                            allow_new_levels = TRUE)
propLinFit3BfinalPreds <- t(t(propLinFit3BfinalPreds)*StormdataTest3$HWRF)
propLinFit3BfinalPreds2 <- colMeans(propLinFit3BfinalPreds)
propLinFit3BfinalPredsMed <- apply(propLinFit3BfinalPreds, 2, function(x){quantile(x, 0.5)})
propLinFit3BfinalPredsLCB <- apply(propLinFit3BfinalPreds, 2, function(x){quantile(x, 0.025)})
propLinFit3BfinalPredsUCB <- apply(propLinFit3BfinalPreds, 2, function(x){quantile(x, 0.975)})

propLinFit3BpredMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(propLinFit3BfinalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(propLinFit3BfinalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < propLinFit3BfinalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(propLinFit3BfinalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(propLinFit3BfinalPredsMed - Actual_Yvec)),
  COV_pred = mean(propLinFit3BfinalPredsLCB < Actual_Yvec & Actual_Yvec < propLinFit3BfinalPredsUCB)
)
propLinFit3BpredMetrics

##### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = propLinFit3BfinalPreds) +
  labs(title = "GammaFit5 Predict")

propLinFit3BFitDF <- bind_cols(
  StormdataTrain3,
  LCB = propLinFit3BfinalFitLCB,
  Mean = propLinFit3BfinalFitMean,
  Med = propLinFit3BfinalFitMed,
  UCB = propLinFit3BfinalFitUCB
) 

ggplot(data = propLinFit3BFitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))

## Prediction
propLinFit3BPredDF <- bind_cols(
  StormdataTest3,
  LCB = propLinFit3BfinalPredsLCB,
  Mean = propLinFit3BfinalPreds2,
  Med = propLinFit3BfinalPredsMed,
  UCB = propLinFit3BfinalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

ggplot(data = propLinFit3BPredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))

rm(propLinFit3BfinalFit)
rm(propLinFit3BfinalPreds)

### STUDENT----
#### Model 1 ----
propStudentFit1 <- brm(
  bf(
    VMAX/HWRF ~
      #Year +
      Month +
      basin + 
      LON +
      LAT +
      #s(Day) +
      #s(StormElapsedTime) + 
      #t2(LON, LAT) +
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
      #HWRF +
      (1|StormID)
  ),
  data = StormdataTrain8, 
  family = student(link = "log"), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)
save(propStudentFit1, file = "_data/propStudentFit1.RData")
prior_summary(propStudentFit1)
round(posterior_summary(propStudentFit1, probs = c(0.025, 0.975)))
propStudentFit1

print(propStudentFit1, digits = 4)
plot(propStudentFit1)
pp_check(propStudentFit1, ndraws = 100) + labs(title = "propStudentFit1 PPC")
loo(propStudentFit1)
waic(propStudentFit1)
performance::check_distribution(propStudentFit1)
performance::check_outliers(propStudentFit1)
performance::check_heteroskedasticity(propStudentFit1)
performance_rmse(propStudentFit1)
performance_mae(propStudentFit1)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(propStudentFit1)

variance_decomposition(propStudentFit1)
exp(fixef(propStudentFit1))
ranef(propStudentFit1)

bayes_R2(propStudentFit1)
loo_R2(propStudentFit1)

bayes_factor(propStudentFit1, propLinFit1)
bayes_factor(propStudentFit1, gammaFit3C)
bayes_factor(propStudentFit1, gammaFit3C)
bayes_factor(propStudentFit1, studentFit3C)
bayes_factor(propStudentFit1, studentFit3C)
bayes_factor(propStudentFit1, studentFit3C)
bayes_factor(propStudentFit1, linFit3C1)
bayes_factor(propStudentFit1, propFit3C)
bayes_factor(propStudentFit1, logPropFit3C)
loo(propStudentFit1, gammaFit3C)

propStudentFit1smooths <- conditional_smooths(propStudentFit1)
plot(propStudentFit1smooths, stype = "raster", ask = FALSE)
propStudentFit1effects <- conditional_effects(propStudentFit1, 
                                           method = "posterior_predict",
                                           robust = FALSE,
                                           re_formula = NULL)
propStudentFit1effects <- conditional_effects(propStudentFit1)
plot(propStudentFit1effects, points = TRUE, ask = FALSE)

##### Prediction ----
## Fitted
propStudentFit1finalFit <- posterior_predict(propStudentFit1)
propStudentFit1finalFit2 <- t(t(propStudentFit1finalFit)*StormdataTrain3$HWRF)
#propStudentFit1finalFitMean <- colMeans(propStudentFit1finalFit)*StormdataTrain3$HWRF
propStudentFit1finalFitMean <- colMeans(propStudentFit1finalFit2)
propStudentFit1finalFitMed <- apply(propStudentFit1finalFit2, 2, function(x){quantile(x, 0.5)})
propStudentFit1finalFitLCB <- apply(propStudentFit1finalFit2, 2, function(x){quantile(x, 0.025)})
propStudentFit1finalFitUCB <- apply(propStudentFit1finalFit2, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
propStudentFit1finalPreds <- posterior_predict(propStudentFit1, 
                                            newdata = StormdataTestFinalscale2,
                                            allow_new_levels = TRUE)
propStudentFit1finalPreds <- t(t(propStudentFit1finalPreds)*StormdataTest3$HWRF)
propStudentFit1finalPreds2 <- colMeans(propStudentFit1finalPreds)
propStudentFit1finalPredsMed <- apply(propStudentFit1finalPreds, 2, function(x){quantile(x, 0.5)})
propStudentFit1finalPredsLCB <- apply(propStudentFit1finalPreds, 2, function(x){quantile(x, 0.025)})
propStudentFit1finalPredsUCB <- apply(propStudentFit1finalPreds, 2, function(x){quantile(x, 0.975)})

propStudentFit1predMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(propStudentFit1finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(propStudentFit1finalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < propStudentFit1finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(propStudentFit1finalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(propStudentFit1finalPredsMed - Actual_Yvec)),
  COV_pred = mean(propStudentFit1finalPredsLCB < Actual_Yvec & Actual_Yvec < propStudentFit1finalPredsUCB)
)
propStudentFit1predMetrics

##### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = propStudentFit1finalPreds) +
  labs(title = "GammaFit5 Predict")

propStudentFit1FitDF <- bind_cols(
  StormdataTrain3,
  LCB = propStudentFit1finalFitLCB,
  Mean = propStudentFit1finalFitMean,
  Med = propStudentFit1finalFitMed,
  UCB = propStudentFit1finalFitUCB
) 

ggplot(data = propStudentFit1FitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))

## Prediction
propStudentFit1PredDF <- bind_cols(
  StormdataTest3,
  LCB = propStudentFit1finalPredsLCB,
  Mean = propStudentFit1finalPreds2,
  Med = propStudentFit1finalPredsMed,
  UCB = propStudentFit1finalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

ggplot(data = propStudentFit1PredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))

rm(propStudentFit1finalFit)
rm(propStudentFit1finalFit2)
rm(propStudentFit1finalPreds)

#### Model 2 ----
propStudentFit2 <- brm(
  bf(
    VMAX/HWRF ~
      #Year +
      #Month +
      basin + 
      LON +
      LAT +
      s(Day) +
      #s(StormElapsedTime) + 
      #t2(LON, LAT) +
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
      #HWRF +
      (1|StormID)
  ),
  data = StormdataTrain8, 
  family = student(link = "log"), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)
save(propStudentFit2, file = "_data/propStudentFit2.RData")
prior_summary(propStudentFit2)
round(posterior_summary(propStudentFit2, probs = c(0.025, 0.975)))
propStudentFit2

print(propStudentFit2, digits = 4)
plot(propStudentFit2)
pp_check(propStudentFit2, ndraws = 100) + labs(title = "propStudentFit2 PPC")
loo(propStudentFit2)
waic(propStudentFit2)
performance::check_distribution(propStudentFit2)
performance::check_outliers(propStudentFit2)
performance::check_heteroskedasticity(propStudentFit2)
performance_rmse(propStudentFit2)
performance_mae(propStudentFit2)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(propStudentFit2)

variance_decomposition(propStudentFit2)
exp(fixef(propStudentFit2))
ranef(propStudentFit2)

bayes_R2(propStudentFit2)
loo_R2(propStudentFit2)

bayes_factor(propStudentFit2, propLinFit2)
bayes_factor(propStudentFit2, gammaFit3C)
bayes_factor(propStudentFit2, gammaFit3C)
bayes_factor(propStudentFit2, studentFit3C)
bayes_factor(propStudentFit2, studentFit3C)
bayes_factor(propStudentFit2, studentFit3C)
bayes_factor(propStudentFit2, linFit3C1)
bayes_factor(propStudentFit2, propFit3C)
bayes_factor(propStudentFit2, logPropFit3C)
loo(propStudentFit2, gammaFit3C)

propStudentFit2smooths <- conditional_smooths(propStudentFit2)
plot(propStudentFit2smooths, stype = "raster", ask = FALSE)
propStudentFit2effects <- conditional_effects(propStudentFit2, 
                                              method = "posterior_predict",
                                              robust = FALSE,
                                              re_formula = NULL)
propStudentFit2effects <- conditional_effects(propStudentFit2)
plot(propStudentFit2effects, points = TRUE, ask = FALSE)

##### Prediction ----
## Fitted
propStudentFit2finalFit <- posterior_predict(propStudentFit2)
propStudentFit2finalFit2 <- t(t(propStudentFit2finalFit)*StormdataTrain3$HWRF)
#propStudentFit2finalFitMean <- colMeans(propStudentFit2finalFit)*StormdataTrain3$HWRF
propStudentFit2finalFitMean <- colMeans(propStudentFit2finalFit2)
propStudentFit2finalFitMed <- apply(propStudentFit2finalFit2, 2, function(x){quantile(x, 0.5)})
propStudentFit2finalFitLCB <- apply(propStudentFit2finalFit2, 2, function(x){quantile(x, 0.025)})
propStudentFit2finalFitUCB <- apply(propStudentFit2finalFit2, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
propStudentFit2finalPreds <- posterior_predict(propStudentFit2, 
                                               newdata = StormdataTestFinalscale2,
                                               allow_new_levels = TRUE)
propStudentFit2finalPreds <- t(t(propStudentFit2finalPreds)*StormdataTest3$HWRF)
propStudentFit2finalPreds2 <- colMeans(propStudentFit2finalPreds)
propStudentFit2finalPredsMed <- apply(propStudentFit2finalPreds, 2, function(x){quantile(x, 0.5)})
propStudentFit2finalPredsLCB <- apply(propStudentFit2finalPreds, 2, function(x){quantile(x, 0.025)})
propStudentFit2finalPredsUCB <- apply(propStudentFit2finalPreds, 2, function(x){quantile(x, 0.975)})

propStudentFit2predMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(propStudentFit2finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(propStudentFit2finalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < propStudentFit2finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(propStudentFit2finalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(propStudentFit2finalPredsMed - Actual_Yvec)),
  COV_pred = mean(propStudentFit2finalPredsLCB < Actual_Yvec & Actual_Yvec < propStudentFit2finalPredsUCB)
)
propStudentFit2predMetrics

##### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = propStudentFit2finalPreds) +
  labs(title = "GammaFit5 Predict")

propStudentFit2FitDF <- bind_cols(
  StormdataTrain3,
  LCB = propStudentFit2finalFitLCB,
  Mean = propStudentFit2finalFitMean,
  Med = propStudentFit2finalFitMed,
  UCB = propStudentFit2finalFitUCB
) 

ggplot(data = propStudentFit2FitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))

## Prediction
propStudentFit2PredDF <- bind_cols(
  StormdataTest3,
  LCB = propStudentFit2finalPredsLCB,
  Mean = propStudentFit2finalPreds2,
  Med = propStudentFit2finalPredsMed,
  UCB = propStudentFit2finalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

ggplot(data = propStudentFit2PredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))

rm(propStudentFit2finalFit)
rm(propStudentFit2finalFit2)
rm(propStudentFit2finalPreds)

#### Model 3 ----
propStudentFit3 <- brm(
  bf(
    VMAX/HWRF ~
      #Year +
      #Month +
      basin + 
      LON +
      LAT +
      s(Day) +
      s(StormElapsedTime) + 
      #t2(LON, LAT) +
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
      #HWRF +
      (1|StormID)
  ),
  data = StormdataTrain8, 
  family = student(link = "log"), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)
save(propStudentFit3, file = "_data/propStudentFit3.RData")
prior_summary(propStudentFit3)
round(posterior_summary(propStudentFit3, probs = c(0.025, 0.975)))
propStudentFit3

print(propStudentFit3, digits = 4)
plot(propStudentFit3)
pp_check(propStudentFit3, ndraws = 100) + labs(title = "propStudentFit3 PPC")
pp_check(propStudentFit3, ndraws = 100,
         fun = "ppc_bars", stat = "max")
loo(propStudentFit3)
waic(propStudentFit3)
performance::check_distribution(propStudentFit3)
performance::check_outliers(propStudentFit3)
performance::check_heteroskedasticity(propStudentFit3)
performance_rmse(propStudentFit3)
performance_mae(propStudentFit3)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(propStudentFit3)

variance_decomposition(propStudentFit3)
exp(fixef(propStudentFit3))
ranef(propStudentFit3)

bayes_R2(propStudentFit3)
loo_R2(propStudentFit3)

bayes_factor(propStudentFit3, propLinFit3)
bayes_factor(propStudentFit3, gammaFit3C)
bayes_factor(propStudentFit3, gammaFit3C)
bayes_factor(propStudentFit3, studentFit3C)
bayes_factor(propStudentFit3, studentFit3C)
bayes_factor(propStudentFit3, studentFit3C)
bayes_factor(propStudentFit3, linFit3C1)
bayes_factor(propStudentFit3, propFit3C)
bayes_factor(propStudentFit3, logPropFit3C)
loo(propStudentFit3, gammaFit3C)

propStudentFit3smooths <- conditional_smooths(propStudentFit3)
plot(propStudentFit3smooths, stype = "raster", ask = FALSE)
propStudentFit3effects <- conditional_effects(propStudentFit3, 
                                              method = "posterior_predict",
                                              robust = FALSE,
                                              re_formula = NULL)
propStudentFit3effects <- conditional_effects(propStudentFit3)
plot(propStudentFit3effects, points = TRUE, ask = FALSE)

##### Prediction ----
## Fitted
propStudentFit3finalFit <- posterior_predict(propStudentFit3)
propStudentFit3finalFit3 <- t(t(propStudentFit3finalFit)*StormdataTrain3$HWRF)
#propStudentFit3finalFitMean <- colMeans(propStudentFit3finalFit)*StormdataTrain3$HWRF
propStudentFit3finalFitMean <- colMeans(propStudentFit3finalFit3)
propStudentFit3finalFitMed <- apply(propStudentFit3finalFit3, 2, function(x){quantile(x, 0.5)})
propStudentFit3finalFitLCB <- apply(propStudentFit3finalFit3, 2, function(x){quantile(x, 0.025)})
propStudentFit3finalFitUCB <- apply(propStudentFit3finalFit3, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
propStudentFit3finalPreds <- posterior_predict(propStudentFit3, 
                                               newdata = StormdataTestFinalscale2,
                                               allow_new_levels = TRUE)
propStudentFit3finalPreds <- t(t(propStudentFit3finalPreds)*StormdataTest3$HWRF)
propStudentFit3finalPreds2 <- colMeans(propStudentFit3finalPreds)
propStudentFit3finalPredsMed <- apply(propStudentFit3finalPreds, 2, function(x){quantile(x, 0.5)})
propStudentFit3finalPredsLCB <- apply(propStudentFit3finalPreds, 2, function(x){quantile(x, 0.025)})
propStudentFit3finalPredsUCB <- apply(propStudentFit3finalPreds, 2, function(x){quantile(x, 0.975)})

propStudentFit3predMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(propStudentFit3finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(propStudentFit3finalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < propStudentFit3finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(propStudentFit3finalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(propStudentFit3finalPredsMed - Actual_Yvec)),
  COV_pred = mean(propStudentFit3finalPredsLCB < Actual_Yvec & Actual_Yvec < propStudentFit3finalPredsUCB)
)
propStudentFit3predMetrics

##### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = propStudentFit3finalPreds) +
  labs(title = "GammaFit5 Predict")

propStudentFit3FitDF <- bind_cols(
  StormdataTrain3,
  LCB = propStudentFit3finalFitLCB,
  Mean = propStudentFit3finalFitMean,
  Med = propStudentFit3finalFitMed,
  UCB = propStudentFit3finalFitUCB
) 

ggplot(data = propStudentFit3FitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))

## Prediction
propStudentFit3PredDF <- bind_cols(
  StormdataTest3,
  LCB = propStudentFit3finalPredsLCB,
  Mean = propStudentFit3finalPreds2,
  Med = propStudentFit3finalPredsMed,
  UCB = propStudentFit3finalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

ggplot(data = propStudentFit3PredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))

ppc_q2.5_plot2 <- 
  ppc_stat(Y2train, YppcSamps2Comb, stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = "2.5% Quantile") +
  theme_bw() +
  legend_none()

ppcL2dens <- density(ppcL2)
ppcL2dens <- cbind(ppcL2dens$x, ppcL2dens$y)
ppcL2densB <- ppcL2dens[between(ppcL2dens[,1], quantile(ppcL2, 0.025), quantile(ppcL2, 0.975)), ] 
ppc_q2.5_plot2B <- ggplot() +
  # geom_histogram(aes(x = ppcL2,  after_stat(density)),
  #                fill = "#bcdcdc", color = "#99c7c7") +
  geom_ribbon(aes(x = ppcL2densB[,1], ymin = 0, ymax = ppcL2densB[,2]),
              fill = "#bcdcdc") +
  geom_density(aes(x = ppcL2), color = "#99c7c7", linewidth = .75) +
  #geom_vline(aes(xintercept = quantile(ppcL2, 0.975)), color = "#007C7C", linewidth = 2) +
  geom_vline(aes(xintercept = DppcL2), color = "#007C7C", linewidth = 1) +
  geom_text(aes(label = paste0("p-value = ", round(mean(ppcL2 > DppcL2), 3))),
            x = 0.84*max(ppcL2), y = max(ppcL2dens[,2]), size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  labs(title = "2.5% Quantile",
       x = NULL,
       y = "Posterior Density") +
  theme_bw() +
  theme(
    plot.title = element_text(size = rel(1)))
ppc_q2.5_plot2B
rm(propStudentFit3finalFit)
rm(propStudentFit3finalFit3)
rm(propStudentFit3finalPreds)

#### Model 4 ----
propStudentFit4 <- brm(
  bf(
    VMAX/HWRF ~
      #Year +
      #Month +
      basin + 
      #LON +
      #LAT +
      s(Day) +
      s(StormElapsedTime) + 
      t2(LON, LAT) +
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
      #HWRF +
      (1|StormID)
  ),
  data = StormdataTrain8, 
  family = student(link = "log"), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)
save(propStudentFit4, file = "_data/propStudentFit4.RData")
prior_summary(propStudentFit4)
round(posterior_summary(propStudentFit4, probs = c(0.025, 0.975)))
propStudentFit4

print(propStudentFit4, digits = 4)
plot(propStudentFit4)
pp_check(propStudentFit4, ndraws = 100) + labs(title = "propStudentFit4 PPC")
pp_check(propStudentFit4, ndraws = 100,
         fun = "ppc_bars", stat = "max")
loo(propStudentFit4)
waic(propStudentFit4)
performance::check_distribution(propStudentFit4)
performance::check_outliers(propStudentFit4)
performance::check_heteroskedasticity(propStudentFit4)
performance_rmse(propStudentFit4)
performance_mae(propStudentFit4)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(propStudentFit4)

variance_decomposition(propStudentFit4)
exp(fixef(propStudentFit4))
ranef(propStudentFit4)

bayes_R2(propStudentFit4)
loo_R2(propStudentFit4)

bayes_factor(propStudentFit4, propStudentFit1)
bayes_factor(propStudentFit4, propStudentFit2)
bayes_factor(propStudentFit4, propStudentFit3)
bayes_factor(propStudentFit4, propStudentFit5)
bayes_factor(propStudentFit4, studentFit4C)
bayes_factor(propStudentFit4, studentFit4C)
bayes_factor(propStudentFit4, linFit4C1)
bayes_factor(propStudentFit4, propFit4C)
bayes_factor(propStudentFit4, logPropFit4C)
loo(propStudentFit4, gammaFit4C)

propStudentFit4smooths <- conditional_smooths(propStudentFit4)
plot(propStudentFit4smooths, stype = "raster", ask = FALSE)
propStudentFit4effects <- conditional_effects(propStudentFit4, 
                                              method = "posterior_predict",
                                              robust = FALSE,
                                              re_formula = NULL)
propStudentFit4effects <- conditional_effects(propStudentFit4)
plot(propStudentFit4effects, points = TRUE, ask = FALSE)

##### Prediction ----
## Fitted
propStudentFit4finalFit <- posterior_predict(propStudentFit4)
propStudentFit4finalFit4 <- t(t(propStudentFit4finalFit)*StormdataTrain3$HWRF)
#propStudentFit4finalFitMean <- colMeans(propStudentFit4finalFit)*StormdataTrain3$HWRF
propStudentFit4finalFitMean <- colMeans(propStudentFit4finalFit4)
propStudentFit4finalFitMed <- apply(propStudentFit4finalFit4, 2, function(x){quantile(x, 0.5)})
propStudentFit4finalFitLCB <- apply(propStudentFit4finalFit4, 2, function(x){quantile(x, 0.025)})
propStudentFit4finalFitUCB <- apply(propStudentFit4finalFit4, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
propStudentFit4finalPreds <- posterior_predict(propStudentFit4, 
                                               newdata = StormdataTestFinalscale2,
                                               allow_new_levels = TRUE)
propStudentFit4finalPreds <- t(t(propStudentFit4finalPreds)*StormdataTest3$HWRF)
propStudentFit4finalPreds2 <- colMeans(propStudentFit4finalPreds)
propStudentFit4finalPredsMed <- apply(propStudentFit4finalPreds, 2, function(x){quantile(x, 0.5)})
propStudentFit4finalPredsLCB <- apply(propStudentFit4finalPreds, 2, function(x){quantile(x, 0.025)})
propStudentFit4finalPredsUCB <- apply(propStudentFit4finalPreds, 2, function(x){quantile(x, 0.975)})

propStudentFit4predMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(propStudentFit4finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(propStudentFit4finalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < propStudentFit4finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(propStudentFit4finalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(propStudentFit4finalPredsMed - Actual_Yvec)),
  COV_pred = mean(propStudentFit4finalPredsLCB < Actual_Yvec & Actual_Yvec < propStudentFit4finalPredsUCB)
)
propStudentFit4predMetrics

##### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = propStudentFit4finalPreds) +
  labs(title = "GammaFit5 Predict")

propStudentFit4FitDF <- bind_cols(
  StormdataTrain3,
  LCB = propStudentFit4finalFitLCB,
  Mean = propStudentFit4finalFitMean,
  Med = propStudentFit4finalFitMed,
  UCB = propStudentFit4finalFitUCB
) 

ggplot(data = propStudentFit4FitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))

## Prediction
propStudentFit4PredDF <- bind_cols(
  StormdataTest3,
  LCB = propStudentFit4finalPredsLCB,
  Mean = propStudentFit4finalPreds2,
  Med = propStudentFit4finalPredsMed,
  UCB = propStudentFit4finalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

ggplot(data = propStudentFit4PredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))

ppc_q2.5_plot2 <- 
  ppc_stat(Y2train, YppcSamps2Comb, stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = "2.5% Quantile") +
  theme_bw() +
  legend_none()

ppcL2dens <- density(ppcL2)
ppcL2dens <- cbind(ppcL2dens$x, ppcL2dens$y)
ppcL2densB <- ppcL2dens[between(ppcL2dens[,1], quantile(ppcL2, 0.025), quantile(ppcL2, 0.975)), ] 
ppc_q2.5_plot2B <- ggplot() +
  # geom_histogram(aes(x = ppcL2,  after_stat(density)),
  #                fill = "#bcdcdc", color = "#99c7c7") +
  geom_ribbon(aes(x = ppcL2densB[,1], ymin = 0, ymax = ppcL2densB[,2]),
              fill = "#bcdcdc") +
  geom_density(aes(x = ppcL2), color = "#99c7c7", linewidth = .75) +
  #geom_vline(aes(xintercept = quantile(ppcL2, 0.975)), color = "#007C7C", linewidth = 2) +
  geom_vline(aes(xintercept = DppcL2), color = "#007C7C", linewidth = 1) +
  geom_text(aes(label = paste0("p-value = ", round(mean(ppcL2 > DppcL2), 3))),
            x = 0.84*max(ppcL2), y = max(ppcL2dens[,2]), size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  labs(title = "2.5% Quantile",
       x = NULL,
       y = "Posterior Density") +
  theme_bw() +
  theme(
    plot.title = element_text(size = rel(1)))
ppc_q2.5_plot2B
rm(propStudentFit4finalFit)
rm(propStudentFit4finalFit4)
rm(propStudentFit4finalPreds)

#### Model 5 ----
propStudentFit5 <- brm(
  bf(
    VMAX/HWRF ~
      #Year +
      #Month +
      basin + 
      LON +
      LAT +
      s(Day) +
      StormElapsedTime + 
      #t2(LON, LAT) +
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
      #HWRF +
      (1+StormElapsedTime|StormID)
  ),
  data = StormdataTrain8, 
  family = student(link = "log"), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)
save(propStudentFit5, file = "_data/propStudentFit5.RData")
prior_summary(propStudentFit5)
round(posterior_summary(propStudentFit5, probs = c(0.025, 0.975)))
propStudentFit5

print(propStudentFit5, digits = 4)
plot(propStudentFit5)
pp_check(propStudentFit5, ndraws = 100) + labs(title = "propStudentFit5 PPC")
pp_check(propStudentFit5, ndraws = 100,
         fun = "ppc_bars", stat = "max")
loo(propStudentFit5)
waic(propStudentFit5)
performance::check_distribution(propStudentFit5)
performance::check_outliers(propStudentFit5)
performance::check_heteroskedasticity(propStudentFit5)
performance_rmse(propStudentFit5)
performance_mae(propStudentFit5)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(propStudentFit5)

variance_decomposition(propStudentFit5)
exp(fixef(propStudentFit5))
ranef(propStudentFit5)

bayes_R2(propStudentFit5)
loo_R2(propStudentFit5)

bayes_factor(propStudentFit5, propLinFit5)
bayes_factor(propStudentFit5, gammaFit5C)
bayes_factor(propStudentFit5, gammaFit5C)
bayes_factor(propStudentFit5, studentFit5C)
bayes_factor(propStudentFit5, studentFit5C)
bayes_factor(propStudentFit5, studentFit5C)
bayes_factor(propStudentFit5, linFit5C1)
bayes_factor(propStudentFit5, propFit5C)
bayes_factor(propStudentFit5, logPropFit5C)
loo(propStudentFit5, gammaFit5C)

propStudentFit5smooths <- conditional_smooths(propStudentFit5)
plot(propStudentFit5smooths, stype = "raster", ask = FALSE)
propStudentFit5effects <- conditional_effects(propStudentFit5, 
                                              method = "posterior_predict",
                                              robust = FALSE,
                                              re_formula = NULL)
propStudentFit5effects <- conditional_effects(propStudentFit5)
plot(propStudentFit5effects, points = TRUE, ask = FALSE)

##### Prediction ----
## Fitted
propStudentFit5finalFit <- posterior_predict(propStudentFit5)
propStudentFit5finalFit5 <- t(t(propStudentFit5finalFit)*StormdataTrain3$HWRF)
#propStudentFit5finalFitMean <- colMeans(propStudentFit5finalFit)*StormdataTrain3$HWRF
propStudentFit5finalFitMean <- colMeans(propStudentFit5finalFit5)
propStudentFit5finalFitMed <- apply(propStudentFit5finalFit5, 2, function(x){quantile(x, 0.5)})
propStudentFit5finalFitLCB <- apply(propStudentFit5finalFit5, 2, function(x){quantile(x, 0.025)})
propStudentFit5finalFitUCB <- apply(propStudentFit5finalFit5, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
propStudentFit5finalPreds <- posterior_predict(propStudentFit5, 
                                               newdata = StormdataTestFinalscale2,
                                               allow_new_levels = TRUE)
propStudentFit5finalPreds <- t(t(propStudentFit5finalPreds)*StormdataTest3$HWRF)
propStudentFit5finalPreds2 <- colMeans(propStudentFit5finalPreds)
propStudentFit5finalPredsMed <- apply(propStudentFit5finalPreds, 2, function(x){quantile(x, 0.5)})
propStudentFit5finalPredsLCB <- apply(propStudentFit5finalPreds, 2, function(x){quantile(x, 0.025)})
propStudentFit5finalPredsUCB <- apply(propStudentFit5finalPreds, 2, function(x){quantile(x, 0.975)})

propStudentFit5predMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(propStudentFit5finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(propStudentFit5finalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < propStudentFit5finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(propStudentFit5finalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(propStudentFit5finalPredsMed - Actual_Yvec)),
  COV_pred = mean(propStudentFit5finalPredsLCB < Actual_Yvec & Actual_Yvec < propStudentFit5finalPredsUCB)
)
propStudentFit5predMetrics

##### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = propStudentFit5finalPreds) +
  labs(title = "GammaFit5 Predict")

propStudentFit5FitDF <- bind_cols(
  StormdataTrain3,
  LCB = propStudentFit5finalFitLCB,
  Mean = propStudentFit5finalFitMean,
  Med = propStudentFit5finalFitMed,
  UCB = propStudentFit5finalFitUCB
) 

ggplot(data = propStudentFit5FitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))

## Prediction
propStudentFit5PredDF <- bind_cols(
  StormdataTest3,
  LCB = propStudentFit5finalPredsLCB,
  Mean = propStudentFit5finalPreds2,
  Med = propStudentFit5finalPredsMed,
  UCB = propStudentFit5finalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

ggplot(data = propStudentFit5PredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))

ppc_q2.5_plot2 <- 
  ppc_stat(Y2train, YppcSamps2Comb, stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = "2.5% Quantile") +
  theme_bw() +
  legend_none()

ppcL2dens <- density(ppcL2)
ppcL2dens <- cbind(ppcL2dens$x, ppcL2dens$y)
ppcL2densB <- ppcL2dens[between(ppcL2dens[,1], quantile(ppcL2, 0.025), quantile(ppcL2, 0.975)), ] 
ppc_q2.5_plot2B <- ggplot() +
  # geom_histogram(aes(x = ppcL2,  after_stat(density)),
  #                fill = "#bcdcdc", color = "#99c7c7") +
  geom_ribbon(aes(x = ppcL2densB[,1], ymin = 0, ymax = ppcL2densB[,2]),
              fill = "#bcdcdc") +
  geom_density(aes(x = ppcL2), color = "#99c7c7", linewidth = .75) +
  #geom_vline(aes(xintercept = quantile(ppcL2, 0.975)), color = "#007C7C", linewidth = 2) +
  geom_vline(aes(xintercept = DppcL2), color = "#007C7C", linewidth = 1) +
  geom_text(aes(label = paste0("p-value = ", round(mean(ppcL2 > DppcL2), 3))),
            x = 0.84*max(ppcL2), y = max(ppcL2dens[,2]), size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  labs(title = "2.5% Quantile",
       x = NULL,
       y = "Posterior Density") +
  theme_bw() +
  theme(
    plot.title = element_text(size = rel(1)))
ppc_q2.5_plot2B
rm(propStudentFit5finalFit)
rm(propStudentFit5finalFit5)
rm(propStudentFit5finalPreds)

### GAMMA ----
#### Model 1 ----
propGammaFit1 <- brm(
  bf(
    VMAX/HWRF ~
      #Year +
      Month +
      basin + 
      LON +
      LAT +
      #s(Day) +
      #s(StormElapsedTime) + 
      #t2(LON, LAT) +
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
      #HWRF +
      (1|StormID)
  ),
  data = StormdataTrain8, 
  family = brmsfamily("Gamma", link = "log"), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)
save(propGammaFit1, file = "_data/propGammaFit1.RData")
prior_summary(propGammaFit1)
round(posterior_summary(propGammaFit1, probs = c(0.025, 0.975)))
propGammaFit1

print(propGammaFit1, digits = 4)
pp_check(propGammaFit1, ndraws = 100)
plot(propGammaFit1)
loo(propGammaFit1)
waic(propGammaFit1)
performance::check_distribution(propGammaFit1)
performance::check_outliers(propGammaFit1)
performance::check_heteroskedasticity(propGammaFit1)
performance_rmse(propGammaFit1)
performance_mae(propGammaFit1)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(propGammaFit1)


variance_decomposition(propGammaFit1)
exp(fixef(propGammaFit1))
ranef(propGammaFit1)

bayes_R2(propGammaFit1)
loo_R2(propGammaFit1)

bayes_factor(propGammaFit1, propLinFit1)
bayes_factor(propGammaFit1, gammaFit7)
bayes_factor(propGammaFit1, gammaFit7)
bayes_factor(propGammaFit1, studentFit7)
bayes_factor(propGammaFit1, studentFit7)
bayes_factor(propGammaFit1, studentFit7)
bayes_factor(propGammaFit1, linFit71)
bayes_factor(propGammaFit1, propFit7)
bayes_factor(propGammaFit1, logPropFit7)
loo(propGammaFit1, gammaFit7)

propGammaFit1smooths <- conditional_smooths(propGammaFit1)
plot(propGammaFit1smooths, stype = "raster", ask = FALSE)
propGammaFit1effects <- conditional_effects(propGammaFit1, 
                                          method = "posterior_predict",
                                          robust = FALSE,
                                          re_formula = NULL)
propGammaFit1effects <- conditional_effects(propGammaFit1)
plot(propGammaFit1effects, points = TRUE, ask = FALSE)

##### Prediction ----
## Fitted
propGammaFit1finalFit <- posterior_predict(propGammaFit1)
propGammaFit1finalFit2 <- t(t(propGammaFit1finalFit)*StormdataTrain3$HWRF)
#propGammaFit1finalFitMean <- colMeans(propGammaFit1finalFit)*StormdataTrain3$HWRF
propGammaFit1finalFitMean <- colMeans(propGammaFit1finalFit2)
propGammaFit1finalFitMed <- apply(propGammaFit1finalFit2, 2, function(x){quantile(x, 0.5)})
propGammaFit1finalFitLCB <- apply(propGammaFit1finalFit2, 2, function(x){quantile(x, 0.025)})
propGammaFit1finalFitUCB <- apply(propGammaFit1finalFit2, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
propGammaFit1finalPreds <- posterior_predict(propGammaFit1, 
                                           newdata = StormdataTestFinalscale2,
                                           allow_new_levels = TRUE)
propGammaFit1finalPreds <- t(t(propGammaFit1finalPreds)*StormdataTest3$HWRF)
propGammaFit1finalPreds2 <- colMeans(propGammaFit1finalPreds)
propGammaFit1finalPredsMed <- apply(propGammaFit1finalPreds, 2, function(x){quantile(x, 0.5)})
propGammaFit1finalPredsLCB <- apply(propGammaFit1finalPreds, 2, function(x){quantile(x, 0.025)})
propGammaFit1finalPredsUCB <- apply(propGammaFit1finalPreds, 2, function(x){quantile(x, 0.975)})

propGammaFit1predMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(propGammaFit1finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(propGammaFit1finalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < propGammaFit1finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(propGammaFit1finalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(propGammaFit1finalPredsMed - Actual_Yvec)),
  COV_pred = mean(propGammaFit1finalPredsLCB < Actual_Yvec & Actual_Yvec < propGammaFit1finalPredsUCB)
)
propGammaFit1predMetrics

##### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = propGammaFit1finalPreds) +
  labs(title = "propGammaFit1 Predict")

propGammaFit1FitDF <- bind_cols(
  StormdataTrain3,
  LCB = propGammaFit1finalFitLCB,
  Mean = propGammaFit1finalFitMean,
  Med = propGammaFit1finalFitMed,
  UCB = propGammaFit1finalFitUCB
) 

ggplot(data = propGammaFit1FitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))

## Prediction
propGammaFit1PredDF <- bind_cols(
  StormdataTest3,
  LCB = propGammaFit1finalPredsLCB,
  Mean = propGammaFit1finalPreds2,
  Med = propGammaFit1finalPredsMed,
  UCB = propGammaFit1finalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

ggplot(data = propGammaFit1PredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))

rm(propGammaFit1finalFit)
rm(propGammaFit1finalFit2)
rm(propGammaFit1finalPreds)

## log(VMAX/HWRF) ----
### LINEAR ----
#### Model 1 ----
logpropLinFit1 <- brm(
  bf(
    log(VMAX/HWRF) ~
      #Year +
      Month +
      basin + 
      LON +
      LAT +
      #s(Day) +
      #s(StormElapsedTime) + 
      #t2(LON, LAT) +
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
      #HWRF +
      (1|StormID)
  ),
  data = StormdataTrain8, 
  family = gaussian(), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)
save(logpropLinFit1, file = "_data/logpropLinFit1.RData")
prior_summary(logpropLinFit1)
round(posterior_summary(logpropLinFit1, probs = c(0.025, 0.975)))
logpropLinFit1

print(logpropLinFit1, digits = 4)
plot(logpropLinFit1)
pp_check(logpropLinFit1, ndraws = 100) + labs(title = "logpropLinFit1 PPC")
loo(logpropLinFit1)
waic(logpropLinFit1)
performance::check_distribution(logpropLinFit1)
performance::check_outliers(logpropLinFit1)
performance::check_heteroskedasticity(logpropLinFit1)
performance_rmse(logpropLinFit1)
performance_mae(logpropLinFit1)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(logpropLinFit1)


variance_decomposition(logpropLinFit1)
exp(fixef(logpropLinFit1))
ranef(logpropLinFit1)

bayes_R2(logpropLinFit1)

bayes_factor(logpropLinFit1, propLinFit1)
bayes_factor(logpropLinFit1, gammaFit5)
bayes_factor(logpropLinFit1, gammaFit5)
bayes_factor(logpropLinFit1, studentFit5)
bayes_factor(logpropLinFit1, studentFit5)
bayes_factor(logpropLinFit1, studentFit5)
bayes_factor(logpropLinFit1, linFit51)
bayes_factor(logpropLinFit1, propFit5)
bayes_factor(logpropLinFit1, logPropFit5)
loo(logpropLinFit1, gammaFit5)

logpropLinFit1smooths <- conditional_smooths(logpropLinFit1)
plot(logpropLinFit1smooths, stype = "raster", ask = FALSE)
logpropLinFit1effects <- conditional_effects(logpropLinFit1, 
                                          method = "posterior_predict",
                                          robust = FALSE,
                                          re_formula = NULL)
logpropLinFit1effects <- conditional_effects(logpropLinFit1)
plot(logpropLinFit1effects, points = TRUE, ask = FALSE)

##### Prediction ----
## Fitted
logpropLinFit1finalFit <- posterior_predict(logpropLinFit1)
logpropLinFit1finalFit <- t(t(exp(logpropLinFit1finalFit))*StormdataTrain3$HWRF)
#logpropLinFit1finalFitMean <- colMeans(logpropLinFit1finalFit)*StormdataTrain3$HWRF
logpropLinFit1finalFitMean <- colMeans(logpropLinFit1finalFit)
logpropLinFit1finalFitMed <- apply(logpropLinFit1finalFit, 2, function(x){quantile(x, 0.5)})
logpropLinFit1finalFitLCB <- apply(logpropLinFit1finalFit, 2, function(x){quantile(x, 0.025)})
logpropLinFit1finalFitUCB <- apply(logpropLinFit1finalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
logpropLinFit1finalPreds <- posterior_predict(logpropLinFit1, 
                                           newdata = StormdataTestFinalscale2,
                                           allow_new_levels = TRUE)
logpropLinFit1finalPreds <- t(t(exp(logpropLinFit1finalPreds))*StormdataTest3$HWRF)
logpropLinFit1finalPreds2 <- colMeans(logpropLinFit1finalPreds)
logpropLinFit1finalPredsMed <- apply(logpropLinFit1finalPreds, 2, function(x){quantile(x, 0.5)})
logpropLinFit1finalPredsLCB <- apply(logpropLinFit1finalPreds, 2, function(x){quantile(x, 0.025)})
logpropLinFit1finalPredsUCB <- apply(logpropLinFit1finalPreds, 2, function(x){quantile(x, 0.975)})

logpropLinFit1predMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(logpropLinFit1finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(logpropLinFit1finalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < logpropLinFit1finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(logpropLinFit1finalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(logpropLinFit1finalPredsMed - Actual_Yvec)),
  COV_pred = mean(logpropLinFit1finalPredsLCB < Actual_Yvec & Actual_Yvec < logpropLinFit1finalPredsUCB)
)
logpropLinFit1predMetrics

##### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = logpropLinFit1finalPreds) +
  labs(title = "GammaFit5 Predict")

logpropLinFit1FitDF <- bind_cols(
  StormdataTrain3,
  LCB = logpropLinFit1finalFitLCB,
  Mean = logpropLinFit1finalFitMean,
  Med = logpropLinFit1finalFitMed,
  UCB = logpropLinFit1finalFitUCB
) 

ggplot(data = logpropLinFit1FitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))

## Prediction
logpropLinFit1PredDF <- bind_cols(
  StormdataTest3,
  LCB = logpropLinFit1finalPredsLCB,
  Mean = logpropLinFit1finalPreds2,
  Med = logpropLinFit1finalPredsMed,
  UCB = logpropLinFit1finalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

ggplot(data = logpropLinFit1PredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))

rm(logpropLinFit1finalFit)
rm(logpropLinFit1finalPreds)

### STUDENT ----
#### Model 1 ----
logpropStudentFit1 <- brm(
  bf(
    log(VMAX/HWRF) ~
      #Year +
      Month +
      basin + 
      LON +
      LAT +
      #s(Day) +
      #s(StormElapsedTime) + 
      #t2(LON, LAT) +
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
      #HWRF +
      (1|StormID)
  ),
  data = StormdataTrain8, 
  family = student(), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)
save(logpropStudentFit1, file = "_data/logpropStudentFit1.RData")
prior_summary(logpropStudentFit1)
round(posterior_summary(logpropStudentFit1, probs = c(0.025, 0.975)))
logpropStudentFit1

print(logpropStudentFit1, digits = 4)
pp_check(logpropStudentFit1, ndraws = 100) + labs(title = "logpropStudentFit1 PPC")
plot(logpropStudentFit1)
loo(logpropStudentFit1)
waic(logpropStudentFit1)
performance::check_distribution(logpropStudentFit1)
performance::check_outliers(logpropStudentFit1)
performance::check_heteroskedasticity(logpropStudentFit1)
performance_rmse(logpropStudentFit1)
performance_mae(logpropStudentFit1)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(logpropStudentFit1)


variance_decomposition(logpropStudentFit1)
exp(fixef(logpropStudentFit1))
ranef(logpropStudentFit1)

bayes_R2(logpropStudentFit1)
loo_R2(logpropStudentFit1)

bayes_factor(logpropStudentFit1, propLinFit1)
bayes_factor(logpropStudentFit1, gammaFit6)
bayes_factor(logpropStudentFit1, gammaFit6)
bayes_factor(logpropStudentFit1, studentFit6)
bayes_factor(logpropStudentFit1, studentFit6)
bayes_factor(logpropStudentFit1, studentFit6)
bayes_factor(logpropStudentFit1, linFit61)
bayes_factor(logpropStudentFit1, propFit6)
bayes_factor(logpropStudentFit1, logPropFit6)
loo(logpropStudentFit1, gammaFit6)

logpropStudentFit1smooths <- conditional_smooths(logpropStudentFit1)
plot(logpropStudentFit1smooths, stype = "raster", ask = FALSE)
logpropStudentFit1effects <- conditional_effects(logpropStudentFit1, 
                                          method = "posterior_predict",
                                          robust = FALSE,
                                          re_formula = NULL)
logpropStudentFit1effects <- conditional_effects(logpropStudentFit1)
plot(logpropStudentFit1effects, points = TRUE, ask = FALSE)

##### Prediction ----
## Fitted
logpropStudentFit1finalFit <- posterior_predict(logpropStudentFit1)
logpropStudentFit1finalFit2 <- t(t(exp(logpropStudentFit1finalFit))*StormdataTrain3$HWRF)
#logpropStudentFit1finalFitMean <- colMeans(logpropStudentFit1finalFit)*StormdataTrain3$HWRF
logpropStudentFit1finalFitMean <- colMeans(logpropStudentFit1finalFit2)
logpropStudentFit1finalFitMed <- apply(logpropStudentFit1finalFit2, 2, function(x){quantile(x, 0.5)})
logpropStudentFit1finalFitLCB <- apply(logpropStudentFit1finalFit2, 2, function(x){quantile(x, 0.025)})
logpropStudentFit1finalFitUCB <- apply(logpropStudentFit1finalFit2, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
logpropStudentFit1finalPreds <- posterior_predict(logpropStudentFit1, 
                                           newdata = StormdataTestFinalscale2,
                                           allow_new_levels = TRUE)
logpropStudentFit1finalPreds <- t(t(exp(logpropStudentFit1finalPreds))*StormdataTest3$HWRF)
logpropStudentFit1finalPreds2 <- colMeans(logpropStudentFit1finalPreds)
logpropStudentFit1finalPredsMed <- apply(logpropStudentFit1finalPreds, 2, function(x){quantile(x, 0.5)})
logpropStudentFit1finalPredsLCB <- apply(logpropStudentFit1finalPreds, 2, function(x){quantile(x, 0.025)})
logpropStudentFit1finalPredsUCB <- apply(logpropStudentFit1finalPreds, 2, function(x){quantile(x, 0.975)})

logpropStudentFit1predMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(logpropStudentFit1finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(logpropStudentFit1finalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < logpropStudentFit1finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(logpropStudentFit1finalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(logpropStudentFit1finalPredsMed - Actual_Yvec)),
  COV_pred = mean(logpropStudentFit1finalPredsLCB < Actual_Yvec & Actual_Yvec < logpropStudentFit1finalPredsUCB)
)
logpropStudentFit1predMetrics

##### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = logpropStudentFit1finalPreds) +
  labs(title = "logpropStudentFit1 Predict")

logpropStudentFit1FitDF <- bind_cols(
  StormdataTrain3,
  LCB = logpropStudentFit1finalFitLCB,
  Mean = logpropStudentFit1finalFitMean,
  Med = logpropStudentFit1finalFitMed,
  UCB = logpropStudentFit1finalFitUCB
) 

ggplot(data = logpropStudentFit1FitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))

## Prediction
logpropStudentFit1PredDF <- bind_cols(
  StormdataTest3,
  LCB = logpropStudentFit1finalPredsLCB,
  Mean = logpropStudentFit1finalPreds2,
  Med = logpropStudentFit1finalPredsMed,
  UCB = logpropStudentFit1finalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

ggplot(data = logpropStudentFit1PredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))

rm(logpropStudentFit1finalFit)
rm(logpropStudentFit1finalFit2)
rm(logpropStudentFit1finalPreds)

#### Model 2 ----
logpropStudentFit2 <- brm(
  bf(
    log(VMAX/HWRF) ~
      #Month +
      basin + 
      LON +
      LAT +
      Day +
      StormElapsedTime + 
      #t2(LON, LAT) +
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
      #HWRF +
      (1|StormID)
  ),
  data = StormdataTrain8, 
  family = student(), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)
save(logpropStudentFit2, file = "_data/logpropStudentFit2.RData")
prior_summary(logpropStudentFit2)
round(posterior_summary(logpropStudentFit2, probs = c(0.025, 0.975)))
logpropStudentFit2

print(logpropStudentFit2, digits = 4)
pp_check(logpropStudentFit2, ndraws = 100) + labs(title = "logpropStudentFit2 PPC")
plot(logpropStudentFit2)
loo(logpropStudentFit2)
waic(logpropStudentFit2)
performance::check_distribution(logpropStudentFit2)
performance::check_outliers(logpropStudentFit2)
performance::check_heteroskedasticity(logpropStudentFit2)
performance_rmse(logpropStudentFit2)
performance_mae(logpropStudentFit2)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(logpropStudentFit2)


variance_decomposition(logpropStudentFit2)
exp(fixef(logpropStudentFit2))
ranef(logpropStudentFit2)

bayes_R2(logpropStudentFit2)
loo_R2(logpropStudentFit2)

bayes_factor(logpropStudentFit2, propLinFit2)
bayes_factor(logpropStudentFit2, gammaFit6)
bayes_factor(logpropStudentFit2, gammaFit6)
bayes_factor(logpropStudentFit2, studentFit6)
bayes_factor(logpropStudentFit2, studentFit6)
bayes_factor(logpropStudentFit2, studentFit6)
bayes_factor(logpropStudentFit2, linFit61)
bayes_factor(logpropStudentFit2, propFit6)
bayes_factor(logpropStudentFit2, logPropFit6)
loo(logpropStudentFit2, gammaFit6)

logpropStudentFit2smooths <- conditional_smooths(logpropStudentFit2)
plot(logpropStudentFit2smooths, stype = "raster", ask = FALSE)
logpropStudentFit2effects <- conditional_effects(logpropStudentFit2, 
                                                 method = "posterior_predict",
                                                 robust = FALSE,
                                                 re_formula = NULL)
logpropStudentFit2effects <- conditional_effects(logpropStudentFit2)
plot(logpropStudentFit2effects, points = TRUE, ask = FALSE)

##### Prediction ----
## Fitted
logpropStudentFit2finalFit <- posterior_predict(logpropStudentFit2)
logpropStudentFit2finalFit2 <- t(t(exp(logpropStudentFit2finalFit))*StormdataTrain3$HWRF)
#logpropStudentFit2finalFitMean <- colMeans(logpropStudentFit2finalFit)*StormdataTrain3$HWRF
logpropStudentFit2finalFitMean <- colMeans(logpropStudentFit2finalFit2)
logpropStudentFit2finalFitMed <- apply(logpropStudentFit2finalFit2, 2, function(x){quantile(x, 0.5)})
logpropStudentFit2finalFitLCB <- apply(logpropStudentFit2finalFit2, 2, function(x){quantile(x, 0.025)})
logpropStudentFit2finalFitUCB <- apply(logpropStudentFit2finalFit2, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
logpropStudentFit2finalPreds <- posterior_predict(logpropStudentFit2, 
                                                  newdata = StormdataTestFinalscale2,
                                                  allow_new_levels = TRUE)
logpropStudentFit2finalPreds <- t(t(exp(logpropStudentFit2finalPreds))*StormdataTest3$HWRF)
logpropStudentFit2finalPreds2 <- colMeans(logpropStudentFit2finalPreds)
logpropStudentFit2finalPredsMed <- apply(logpropStudentFit2finalPreds, 2, function(x){quantile(x, 0.5)})
logpropStudentFit2finalPredsLCB <- apply(logpropStudentFit2finalPreds, 2, function(x){quantile(x, 0.025)})
logpropStudentFit2finalPredsUCB <- apply(logpropStudentFit2finalPreds, 2, function(x){quantile(x, 0.975)})

logpropStudentFit2predMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(logpropStudentFit2finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(logpropStudentFit2finalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < logpropStudentFit2finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(logpropStudentFit2finalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(logpropStudentFit2finalPredsMed - Actual_Yvec)),
  COV_pred = mean(logpropStudentFit2finalPredsLCB < Actual_Yvec & Actual_Yvec < logpropStudentFit2finalPredsUCB)
)
logpropStudentFit2predMetrics

##### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = logpropStudentFit2finalPreds) +
  labs(title = "logpropStudentFit2 Predict")

logpropStudentFit2FitDF <- bind_cols(
  StormdataTrain3,
  LCB = logpropStudentFit2finalFitLCB,
  Mean = logpropStudentFit2finalFitMean,
  Med = logpropStudentFit2finalFitMed,
  UCB = logpropStudentFit2finalFitUCB
) 

ggplot(data = logpropStudentFit2FitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))

## Prediction
logpropStudentFit2PredDF <- bind_cols(
  StormdataTest3,
  LCB = logpropStudentFit2finalPredsLCB,
  Mean = logpropStudentFit2finalPreds2,
  Med = logpropStudentFit2finalPredsMed,
  UCB = logpropStudentFit2finalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

ggplot(data = logpropStudentFit2PredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))

rm(logpropStudentFit2finalFit)
rm(logpropStudentFit2finalFit2)
rm(logpropStudentFit2finalPreds)

#### Model 3 ----
logpropStudentFit3 <- brm(
  bf(
    log(VMAX/HWRF) ~
      #Month +
      basin + 
      LON +
      LAT +
      Day +
      s(StormElapsedTime) + 
      #s(StormElapsedTime, by = StormID, m = 1) +
      #t2(LON, LAT) +
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
      #HWRF +
      (1|StormID)
  ),
  data = StormdataTrain8, 
  family = student(), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)
save(logpropStudentFit3, file = "_data/logpropStudentFit3.RData")
prior_summary(logpropStudentFit3)
round(posterior_summary(logpropStudentFit3, probs = c(0.025, 0.975)))
logpropStudentFit3

print(logpropStudentFit3, digits = 4)
pp_check(logpropStudentFit3, ndraws = 100)
plot(logpropStudentFit3)
loo(logpropStudentFit3)
waic(logpropStudentFit3)
performance::check_distribution(logpropStudentFit3)
performance::check_outliers(logpropStudentFit3)
performance::check_heteroskedasticity(logpropStudentFit3)
performance_rmse(logpropStudentFit3)
performance_mae(logpropStudentFit3)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(logpropStudentFit3)


variance_decomposition(logpropStudentFit3)
exp(fixef(logpropStudentFit3))
ranef(logpropStudentFit3)

bayes_R2(logpropStudentFit3)
loo_R2(logpropStudentFit3)

bayes_factor(logpropStudentFit3, propLinFit2)
bayes_factor(logpropStudentFit3, gammaFit6)
bayes_factor(logpropStudentFit3, gammaFit6)
bayes_factor(logpropStudentFit3, studentFit6)
bayes_factor(logpropStudentFit3, studentFit6)
bayes_factor(logpropStudentFit3, studentFit6)
bayes_factor(logpropStudentFit3, linFit61)
bayes_factor(logpropStudentFit3, propFit6)
bayes_factor(logpropStudentFit3, logPropFit6)
loo(logpropStudentFit3, gammaFit6)

logpropStudentFit3smooths <- conditional_smooths(logpropStudentFit3)
plot(logpropStudentFit3smooths, stype = "raster", ask = FALSE)
logpropStudentFit3effects <- conditional_effects(logpropStudentFit3, 
                                                 method = "posterior_predict",
                                                 robust = FALSE,
                                                 re_formula = NULL)
logpropStudentFit3effects <- conditional_effects(logpropStudentFit3)
plot(logpropStudentFit3effects, points = TRUE, ask = FALSE)

##### Prediction ----
## Fitted
logpropStudentFit3finalFit <- posterior_predict(logpropStudentFit3)
logpropStudentFit3finalFit2 <- t(t(exp(logpropStudentFit3finalFit))*StormdataTrain3$HWRF)
#logpropStudentFit3finalFitMean <- colMeans(logpropStudentFit3finalFit)*StormdataTrain3$HWRF
logpropStudentFit3finalFitMean <- colMeans(logpropStudentFit3finalFit2)
logpropStudentFit3finalFitMed <- apply(logpropStudentFit3finalFit2, 2, function(x){quantile(x, 0.5)})
logpropStudentFit3finalFitLCB <- apply(logpropStudentFit3finalFit2, 2, function(x){quantile(x, 0.025)})
logpropStudentFit3finalFitUCB <- apply(logpropStudentFit3finalFit2, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
logpropStudentFit3finalPreds <- posterior_predict(logpropStudentFit3, 
                                                  newdata = StormdataTestFinalscale2,
                                                  allow_new_levels = TRUE)
logpropStudentFit3finalPreds <- t(t(exp(logpropStudentFit3finalPreds))*StormdataTest3$HWRF)
logpropStudentFit3finalPreds2 <- colMeans(logpropStudentFit3finalPreds)
logpropStudentFit3finalPredsMed <- apply(logpropStudentFit3finalPreds, 2, function(x){quantile(x, 0.5)})
logpropStudentFit3finalPredsLCB <- apply(logpropStudentFit3finalPreds, 2, function(x){quantile(x, 0.025)})
logpropStudentFit3finalPredsUCB <- apply(logpropStudentFit3finalPreds, 2, function(x){quantile(x, 0.975)})

logpropStudentFit3predMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(logpropStudentFit3finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(logpropStudentFit3finalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < logpropStudentFit3finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(logpropStudentFit3finalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(logpropStudentFit3finalPredsMed - Actual_Yvec)),
  COV_pred = mean(logpropStudentFit3finalPredsLCB < Actual_Yvec & Actual_Yvec < logpropStudentFit3finalPredsUCB)
)
logpropStudentFit3predMetrics

##### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = logpropStudentFit3finalPreds) +
  labs(title = "logpropStudentFit3 Predict")

logpropStudentFit3FitDF <- bind_cols(
  StormdataTrain3,
  LCB = logpropStudentFit3finalFitLCB,
  Mean = logpropStudentFit3finalFitMean,
  Med = logpropStudentFit3finalFitMed,
  UCB = logpropStudentFit3finalFitUCB
) 

ggplot(data = logpropStudentFit3FitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))

## Prediction
logpropStudentFit3PredDF <- bind_cols(
  StormdataTest3,
  LCB = logpropStudentFit3finalPredsLCB,
  Mean = logpropStudentFit3finalPreds2,
  Med = logpropStudentFit3finalPredsMed,
  UCB = logpropStudentFit3finalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

ggplot(data = logpropStudentFit3PredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))

rm(logpropStudentFit3finalFit)
rm(logpropStudentFit3finalFit2)
rm(logpropStudentFit3finalPreds)

#### Model 4 ----
logpropStudentFit4 <- brm(
  bf(
    log(VMAX/HWRF) ~
      #Month +
      basin + 
      #LON +
      #LAT +
      s(Day) +
      s(StormElapsedTime) + 
      #s(StormElapsedTime, by = StormID, m = 1) +
      t2(LON, LAT) +
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
      #HWRF +
      (1|StormID)
  ),
  data = StormdataTrain8, 
  family = student(), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)
save(logpropStudentFit4, file = "_data/logpropStudentFit4.RData")
prior_summary(logpropStudentFit4)
round(posterior_summary(logpropStudentFit4, probs = c(0.025, 0.975)))
logpropStudentFit4

print(logpropStudentFit4, digits = 4)
pp_check(logpropStudentFit4, ndraws = 100)
plot(logpropStudentFit4)
loo(logpropStudentFit4)
waic(logpropStudentFit4)
performance::check_distribution(logpropStudentFit4)
performance::check_outliers(logpropStudentFit4)
performance::check_heteroskedasticity(logpropStudentFit4)
performance_rmse(logpropStudentFit4)
performance_mae(logpropStudentFit4)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(logpropStudentFit4)


variance_decomposition(logpropStudentFit4)
exp(fixef(logpropStudentFit4))
ranef(logpropStudentFit4)

bayes_R2(logpropStudentFit4)
loo_R2(logpropStudentFit4)

bayes_factor(logpropStudentFit4, logpropStudentFit1)
bayes_factor(logpropStudentFit4, logpropStudentFit2)
bayes_factor(logpropStudentFit4, logpropStudentFit3)
bayes_factor(logpropStudentFit4, logpropStudentFit5)
bayes_factor(logpropStudentFit4, studentFit6)
bayes_factor(logpropStudentFit4, studentFit6)
bayes_factor(logpropStudentFit4, studentFit6)
bayes_factor(logpropStudentFit4, linFit61)
bayes_factor(logpropStudentFit4, propFit6)
bayes_factor(logpropStudentFit4, logPropFit6)
loo(logpropStudentFit4, gammaFit6)

logpropStudentFit4smooths <- conditional_smooths(logpropStudentFit4)
plot(logpropStudentFit4smooths, stype = "raster", ask = FALSE)
logpropStudentFit4effects <- conditional_effects(logpropStudentFit4, 
                                                 method = "posterior_predict",
                                                 robust = FALSE,
                                                 re_formula = NULL)
logpropStudentFit4effects <- conditional_effects(logpropStudentFit4)
plot(logpropStudentFit4effects, points = TRUE, ask = FALSE)

##### Prediction ----
## Fitted
logpropStudentFit4finalFit <- posterior_predict(logpropStudentFit4)
logpropStudentFit4finalFit2 <- t(t(exp(logpropStudentFit4finalFit))*StormdataTrain3$HWRF)
#logpropStudentFit4finalFitMean <- colMeans(logpropStudentFit4finalFit)*StormdataTrain3$HWRF
logpropStudentFit4finalFitMean <- colMeans(logpropStudentFit4finalFit2)
logpropStudentFit4finalFitMed <- apply(logpropStudentFit4finalFit2, 2, function(x){quantile(x, 0.5)})
logpropStudentFit4finalFitLCB <- apply(logpropStudentFit4finalFit2, 2, function(x){quantile(x, 0.025)})
logpropStudentFit4finalFitUCB <- apply(logpropStudentFit4finalFit2, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
logpropStudentFit4finalPreds <- posterior_predict(logpropStudentFit4, 
                                                  newdata = StormdataTestFinalscale2,
                                                  allow_new_levels = TRUE)
logpropStudentFit4finalPreds <- t(t(exp(logpropStudentFit4finalPreds))*StormdataTest3$HWRF)
logpropStudentFit4finalPreds2 <- colMeans(logpropStudentFit4finalPreds)
logpropStudentFit4finalPredsMed <- apply(logpropStudentFit4finalPreds, 2, function(x){quantile(x, 0.5)})
logpropStudentFit4finalPredsLCB <- apply(logpropStudentFit4finalPreds, 2, function(x){quantile(x, 0.025)})
logpropStudentFit4finalPredsUCB <- apply(logpropStudentFit4finalPreds, 2, function(x){quantile(x, 0.975)})

logpropStudentFit4predMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(logpropStudentFit4finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(logpropStudentFit4finalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < logpropStudentFit4finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(logpropStudentFit4finalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(logpropStudentFit4finalPredsMed - Actual_Yvec)),
  COV_pred = mean(logpropStudentFit4finalPredsLCB < Actual_Yvec & Actual_Yvec < logpropStudentFit4finalPredsUCB)
)
logpropStudentFit4predMetrics

##### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = logpropStudentFit4finalPreds) +
  labs(title = "logpropStudentFit4 Predict")

logpropStudentFit4FitDF <- bind_cols(
  StormdataTrain3,
  LCB = logpropStudentFit4finalFitLCB,
  Mean = logpropStudentFit4finalFitMean,
  Med = logpropStudentFit4finalFitMed,
  UCB = logpropStudentFit4finalFitUCB
) 

ggplot(data = logpropStudentFit4FitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))

## Prediction
logpropStudentFit4PredDF <- bind_cols(
  StormdataTest3,
  LCB = logpropStudentFit4finalPredsLCB,
  Mean = logpropStudentFit4finalPreds2,
  Med = logpropStudentFit4finalPredsMed,
  UCB = logpropStudentFit4finalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

ggplot(data = logpropStudentFit4PredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))

rm(logpropStudentFit4finalFit)
rm(logpropStudentFit4finalFit2)
rm(logpropStudentFit4finalPreds)

#### Model 5 ----
logpropStudentFit5 <- brm(
  bf(
    log(VMAX/HWRF) ~
      #Year +
      #Month +
      basin + 
      LON +
      LAT +
      s(Day) +
      StormElapsedTime + 
      #t2(LON, LAT) +
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
      #HWRF +
      (1+StormElapsedTime|StormID)
  ),
  data = StormdataTrain8, 
  family = student(), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)

save(logpropStudentFit5, file = "_data/logpropStudentFit5.RData")
prior_summary(logpropStudentFit5)
round(posterior_summary(logpropStudentFit5, probs = c(0.025, 0.975)))
logpropStudentFit5

print(logpropStudentFit5, digits = 4)
pp_check(logpropStudentFit5, ndraws = 100)
plot(logpropStudentFit5)
loo(logpropStudentFit5)
waic(logpropStudentFit5)
performance::check_distribution(logpropStudentFit5)
performance::check_outliers(logpropStudentFit5)
performance::check_heteroskedasticity(logpropStudentFit5)
performance_rmse(logpropStudentFit5)
performance_mae(logpropStudentFit5)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(logpropStudentFit5)


variance_decomposition(logpropStudentFit5)
exp(fixef(logpropStudentFit5))
ranef(logpropStudentFit5)

bayes_R2(logpropStudentFit5)
loo_R2(logpropStudentFit5)

bayes_factor(logpropStudentFit5, logpropStudentFit1)
bayes_factor(logpropStudentFit5, logpropStudentFit2)
bayes_factor(logpropStudentFit5, logpropStudentFit3)
bayes_factor(logpropStudentFit5, 
             logpropStudentFit4)
bayes_factor(logpropStudentFit5, studentFit6)
bayes_factor(logpropStudentFit5, studentFit6)
bayes_factor(logpropStudentFit5, linFit61)
bayes_factor(logpropStudentFit5, propFit6)
bayes_factor(logpropStudentFit5, logPropFit6)
loo(logpropStudentFit5, gammaFit6)

logpropStudentFit5smooths <- conditional_smooths(logpropStudentFit5)
plot(logpropStudentFit5smooths, stype = "raster", ask = FALSE)
logpropStudentFit5effects <- conditional_effects(logpropStudentFit5, 
                                                 method = "posterior_predict",
                                                 robust = FALSE,
                                                 re_formula = NULL)
logpropStudentFit5effects <- conditional_effects(logpropStudentFit5)
plot(logpropStudentFit5effects, points = TRUE, ask = FALSE)

##### Prediction ----
## Fitted
logpropStudentFit5finalFit <- posterior_predict(logpropStudentFit5)
logpropStudentFit5finalFit2 <- t(t(exp(logpropStudentFit5finalFit))*StormdataTrain3$HWRF)
#logpropStudentFit5finalFitMean <- colMeans(logpropStudentFit5finalFit)*StormdataTrain3$HWRF
logpropStudentFit5finalFitMean <- colMeans(logpropStudentFit5finalFit2)
logpropStudentFit5finalFitMed <- apply(logpropStudentFit5finalFit2, 2, function(x){quantile(x, 0.5)})
logpropStudentFit5finalFitLCB <- apply(logpropStudentFit5finalFit2, 2, function(x){quantile(x, 0.025)})
logpropStudentFit5finalFitUCB <- apply(logpropStudentFit5finalFit2, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
logpropStudentFit5finalPreds <- posterior_predict(logpropStudentFit5, 
                                                  newdata = StormdataTestFinalscale2,
                                                  allow_new_levels = TRUE)
logpropStudentFit5finalPreds <- t(t(exp(logpropStudentFit5finalPreds))*StormdataTest3$HWRF)
logpropStudentFit5finalPreds2 <- colMeans(logpropStudentFit5finalPreds)
logpropStudentFit5finalPredsMed <- apply(logpropStudentFit5finalPreds, 2, function(x){quantile(x, 0.5)})
logpropStudentFit5finalPredsLCB <- apply(logpropStudentFit5finalPreds, 2, function(x){quantile(x, 0.025)})
logpropStudentFit5finalPredsUCB <- apply(logpropStudentFit5finalPreds, 2, function(x){quantile(x, 0.975)})

logpropStudentFit5predMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(logpropStudentFit5finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(logpropStudentFit5finalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < logpropStudentFit5finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(logpropStudentFit5finalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(logpropStudentFit5finalPredsMed - Actual_Yvec)),
  COV_pred = mean(logpropStudentFit5finalPredsLCB < Actual_Yvec & Actual_Yvec < logpropStudentFit5finalPredsUCB)
)
logpropStudentFit5predMetrics

##### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = logpropStudentFit5finalPreds) +
  labs(title = "logpropStudentFit5 Predict")

logpropStudentFit5FitDF <- bind_cols(
  StormdataTrain3,
  LCB = logpropStudentFit5finalFitLCB,
  Mean = logpropStudentFit5finalFitMean,
  Med = logpropStudentFit5finalFitMed,
  UCB = logpropStudentFit5finalFitUCB
) 

ggplot(data = logpropStudentFit5FitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))

## Prediction
logpropStudentFit5PredDF <- bind_cols(
  StormdataTest3,
  LCB = logpropStudentFit5finalPredsLCB,
  Mean = logpropStudentFit5finalPreds2,
  Med = logpropStudentFit5finalPredsMed,
  UCB = logpropStudentFit5finalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

ggplot(data = logpropStudentFit5PredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))

rm(logpropStudentFit5finalFit)
rm(logpropStudentFit5finalFit2)
rm(logpropStudentFit5finalPreds)

# Compare Predictions ----
propLinFit1loo <- loo(propLinFit1)
propStudentFit1loo <- loo(propStudentFit1)
propStudentFit2loo <- loo(propStudentFit2)
propStudentFit3loo <- loo(propStudentFit3)
propStudentFit4loo <- loo(propStudentFit4)
propStudentFit5loo <- loo(propStudentFit5)
propGammaFit1loo <- loo(propGammaFit1)
logpropLinFit1loo <- loo(logpropLinFit1)
logpropStudentFit1loo <- loo(logpropStudentFit1)
logpropStudentFit2loo <- loo(logpropStudentFit2)
logpropStudentFit3loo <- loo(logpropStudentFit3)
logpropStudentFit4loo <- loo(logpropStudentFit4)
logpropStudentFit5loo <- loo(logpropStudentFit5)

looComp <- loo_compare(propLinFit1loo,
                       propStudentFit1loo,
                       propStudentFit2loo,
                       propStudentFit3loo,
                       propStudentFit4loo,
                       propStudentFit5loo,
                       propGammaFit1loo,
                       logpropLinFit1loo,
                       logpropStudentFit1loo,
                       logpropStudentFit2loo,
                       logpropStudentFit3loo,
                       logpropStudentFit4loo,
                       logpropStudentFit5loo)
looComp
save(looComp, file = "_data/looCompFinal.RData")

predCompMetrics <- bind_rows(
  propLinFit1predMetrics |> bind_cols(Fit = "propLinFit1"),
  propStudentFit1predMetrics |> bind_cols(Fit = "propStudentFit1"),
  propStudentFit2predMetrics |> bind_cols(Fit = "propStudentFit2"),
  propStudentFit3predMetrics |> bind_cols(Fit = "propStudentFit3"),
  propStudentFit4predMetrics |> bind_cols(Fit = "propStudentFit4"),
  propStudentFit5predMetrics |> bind_cols(Fit = "propStudentFit5"),
  propGammaFit1predMetrics |> bind_cols(Fit = "propGammaFit1"),
  logpropLinFit1predMetrics |> bind_cols(Fit = "logpropLinFit1"),
  logpropStudentFit1predMetrics |> bind_cols(Fit = "logpropStudentFit1"),
  logpropStudentFit2predMetrics |> bind_cols(Fit = "logpropStudentFit2"),
  logpropStudentFit3predMetrics |> bind_cols(Fit = "logpropStudentFit3"),
  logpropStudentFit4predMetrics |> bind_cols(Fit = "logpropStudentFit4"),
  logpropStudentFit5predMetrics |> bind_cols(Fit = "logpropStudentFit5"),
)
predCompMetrics <- predCompMetrics |> arrange(MAE_pred)

rankComps1 <- predCompMetrics |> 
  select(Fit) |> 
  mutate(RankPred = 1:nrow(predCompMetrics))

rankComps2 <- looComp |> 
  as_tibble() |>
  mutate(
    RankLoo = 1:nrow(looComp),
    Fit = rownames(looComp)
  ) |>
  select(Fit, RankLoo)

rankComps <- left_join(rankComps1, rankComps2) |>
  mutate(
    OverallRank = RankPred + RankLoo
  ) |>
  arrange(OverallRank)

save(rankComps, file = "_data/rankCompsFinal.RData")


