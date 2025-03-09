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
  facet_wrap(vars(StormID)) +
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

# Fit model ----
## Load Models ----
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

## VMAX ----
### GAMMA ----
#### Model 1 ----
gammaFit1 <- brm(
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
#save(gammaFit1, file = "_data/gammaFit1.RData")
prior_summary(gammaFit1)
round(posterior_summary(gammaFit1, probs = c(0.025, 0.975)))
gammaFit1

print(gammaFit1, digits = 4)
plot(gammaFit1)
gammaFit1ppcFit <- pp_check(gammaFit1, ndraws = 100) + 
  labs(title = "gammaFit1 Fit PPC") +
  theme_bw()
gammaFit1ppcFit
gammaFit1loo <- loo(gammaFit1)
waic(gammaFit1)
performance::check_distribution(gammaFit1)
performance::check_outliers(gammaFit1)
performance::check_heteroskedasticity(gammaFit1)
performance_rmse(gammaFit1)
performance_mae(gammaFit1)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(gammaFit1)


variance_decomposition(gammaFit1)
exp(fixef(gammaFit1))
ranef(gammaFit1)

bayes_R2(gammaFit1)

bayes_factor(gammaFit1, gammaFit1)
bayes_factor(gammaFit1, gammaFit2)
bayes_factor(gammaFit1, gammaFit3)
bayes_factor(gammaFit1, studentFit1)
bayes_factor(gammaFit1, studentFit2)
bayes_factor(gammaFit1, studentFit3)
bayes_factor(gammaFit1, linFit11)
bayes_factor(gammaFit1, propFit1)
bayes_factor(gammaFit1, logPropFit1)
loo(gammaFit1, gammaFit3)

gammaFit1smooths <- conditional_smooths(gammaFit1)
plot(gammaFit1smooths, stype = "raster", ask = FALSE)
gammaFit1effects <- conditional_effects(gammaFit1, 
                                        method = "posterior_predict",
                                        robust = FALSE,
                                        re_formula = NULL)
gammaFit1effects <- conditional_effects(gammaFit1)
plot(gammaFit1effects, points = TRUE, ask = FALSE)

##### Prediction ----
## Fitted
gammaFit1finalFit <- posterior_predict(gammaFit1)
gammaFit1finalFitMean <- colMeans(gammaFit1finalFit)
gammaFit1finalFitMed <- apply(gammaFit1finalFit, 2, function(x){quantile(x, 0.5)})
gammaFit1finalFitLCB <- apply(gammaFit1finalFit, 2, function(x){quantile(x, 0.025)})
gammaFit1finalFitUCB <- apply(gammaFit1finalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
gammaFit1finalPreds <- posterior_predict(gammaFit1, 
                                         newdata = StormdataTest7scale,
                                         allow_new_levels = TRUE)
gammaFit1finalPreds2 <- colMeans(gammaFit1finalPreds)
gammaFit1finalPredsMed <- apply(gammaFit1finalPreds, 2, function(x){quantile(x, 0.5)})
gammaFit1finalPredsLCB <- apply(gammaFit1finalPreds, 2, function(x){quantile(x, 0.025)})
gammaFit1finalPredsUCB <- apply(gammaFit1finalPreds, 2, function(x){quantile(x, 0.975)})

gammaFit1predMetrics <- tibble(
  Fit = "gammaFit1",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(gammaFit1finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(gammaFit1finalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < gammaFit1finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(gammaFit1finalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(gammaFit1finalPredsMed - Actual_Yvec)),
  COV_pred = mean(gammaFit1finalPredsLCB < Actual_Yvec & Actual_Yvec < gammaFit1finalPredsUCB)
)
gammaFit1predMetrics

##### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = gammaFit1finalPreds) +
  labs(title = "gammaFit1 Predict") +
  theme_bw()

gammaFit1FitDF <- bind_cols(
  StormdataTrain3,
  LCB = gammaFit1finalFitLCB,
  Mean = gammaFit1finalFitMean,
  Med = gammaFit1finalFitMed,
  UCB = gammaFit1finalFitUCB
)

gammaFit1stormsFitplot <- ggplot(data = gammaFit1FitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))
gammaFit1stormsFitplot

## Prediction
gammaFit1PredDF <- bind_cols(
  StormdataTest3,
  LCB = gammaFit1finalPredsLCB,
  Mean = gammaFit1finalPreds2,
  Med = gammaFit1finalPredsMed,
  UCB = gammaFit1finalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

gammaFit1stormsPredplot <- ggplot(data = gammaFit1PredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX, color = "Observed")) +
  geom_line(aes(y = Mean, color = "PPD Mean")) +
  facet_wrap(vars(StormID), ncol = 6)+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50)) +
  labs(title = "gammaFit1 PPD Mean vs Observed VMAX",
       subtitle = "95% Credible Interval about PPD Mean") +
  scale_color_manual(name = NULL, values = c("black","red")) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1))
  ) +
  theme_bw()
gammaFit1stormsPredplot

##### PPC ----
###### Quantile 2.5 
gammaFit1LCBsims <- apply(gammaFit1finalFit, 
                          MARGIN = 1,
                          function(x){
                            quantile(x, 0.025)
                          })
gammaFit1LCBpvalueVec <- gammaFit1LCBsims < quantile(StormdataTrain3$VMAX, 0.025)
gammaFit1LCBpvalue <- sum(gammaFit1LCBpvalueVec)
gammaFit1LCBpvalue <- round(gammaFit1LCBpvalue/4000, 3)
gammaFit1LCBpvalue <- min(gammaFit1LCBpvalue, 1 - gammaFit1LCBpvalue)

gammaFit1_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFit1finalFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", gammaFit1LCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFit1_ppcLCB

###### Quantile 97.5
gammaFit1UCBsims <- apply(gammaFit1finalFit, 
                          MARGIN = 1,
                          function(x){
                            quantile(x, 0.975)
                          })
gammaFit1UCBpvalueVec <- gammaFit1UCBsims < quantile(StormdataTrain3$VMAX, 0.975)
gammaFit1UCBpvalue <- as.numeric(sum(gammaFit1UCBpvalueVec))
gammaFit1UCBpvalue <- round(gammaFit1UCBpvalue/4000, 3)
gammaFit1UCBpvalue <- min(gammaFit1UCBpvalue, 1 - gammaFit1UCBpvalue)

gammaFit1_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFit1finalFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", gammaFit1UCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFit1_ppcUCB

###### Mean 
gammaFit1MEANsims <- apply(gammaFit1finalFit, 
                           MARGIN = 1,
                           function(x){
                             mean(x)
                           })
gammaFit1MEANpvalueVec <- gammaFit1MEANsims < mean(StormdataTrain3$VMAX)
gammaFit1MEANpvalue <- sum(gammaFit1MEANpvalueVec)
gammaFit1MEANpvalue <- round(gammaFit1MEANpvalue/4000, 3)
gammaFit1MEANpvalue <- min(gammaFit1MEANpvalue, 1 - gammaFit1MEANpvalue)

gammaFit1_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFit1finalFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", gammaFit1MEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFit1_ppcMEAN

###### Med 
gammaFit1MEDsims <- apply(gammaFit1finalFit, 
                          MARGIN = 1,
                          function(x){
                            quantile(x, 0.5)
                          })
gammaFit1MEDpvalueVec <- gammaFit1MEDsims < quantile(StormdataTrain3$VMAX, 0.5)
gammaFit1MEDpvalue <- sum(gammaFit1MEDpvalueVec)
gammaFit1MEDpvalue <- round(gammaFit1MEDpvalue/4000, 3)
gammaFit1MEDpvalue <- min(gammaFit1MEDpvalue, 1 - gammaFit1MEDpvalue)

gammaFit1_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFit1finalFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", gammaFit1MEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFit1_ppcMED

###### SD 
gammaFit1SDsims <- apply(gammaFit1finalFit, 
                         MARGIN = 1,
                         function(x){
                           sd(x)
                         })
gammaFit1SDpvalueVec <- gammaFit1SDsims < sd(StormdataTrain3$VMAX)
gammaFit1SDpvalue <- sum(gammaFit1SDpvalueVec)
gammaFit1SDpvalue <- round(gammaFit1SDpvalue/4000, 3)
gammaFit1SDpvalue <- min(gammaFit1SDpvalue, 1 - gammaFit1SDpvalue)

gammaFit1_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFit1finalFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", gammaFit1SDpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFit1_ppcSD

###### Range 
gammaFit1RANGEsims <- apply(gammaFit1finalFit, 
                            MARGIN = 1,
                            function(x){
                              max(x)-min(x)
                            })
gammaFit1RANGEpvalueVec <- gammaFit1RANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
gammaFit1RANGEpvalue <- sum(gammaFit1RANGEpvalueVec)
gammaFit1RANGEpvalue <- round(gammaFit1RANGEpvalue/4000, 3)
gammaFit1RANGEpvalue <- min(gammaFit1RANGEpvalue, 1 - gammaFit1RANGEpvalue)

gammaFit1_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFit1finalFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", gammaFit1RANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFit1_ppcRANGE

##### Combined Plot ----
gammaFit1_ppcComb <- 
  gammaFit1ppcFit /
  (gammaFit1_ppcLCB | gammaFit1_ppcMED | gammaFit1_ppcUCB) /
  (gammaFit1_ppcRANGE | gammaFit1_ppcMEAN | gammaFit1_ppcSD)
  # plot_layout(ncol = 3, 
  #             widths = c(3,1,1,1,1,1,1), 
  #             heights = c(1,1,1,1,1,1,1),
  #             byrow = FALSE)
gammaFit1_ppcComb

##### Bayes p-values ----
gammaFit1pvalues <- tibble(
  Fit = "gammaFit1",
  LCB = gammaFit1LCBpvalue,
  Median = gammaFit1MEDpvalue,
  UCB = gammaFit1UCBpvalue,
  Range = gammaFit1RANGEpvalue,
  Mean = gammaFit1MEANpvalue,
  SD = gammaFit1SDpvalue
)
gammaFit1pvalues

gammaFit1files <- ls()[str_detect(ls(), pattern = "gammaFit1")]
gammaFit1filesRM <- gammaFit1files[!(gammaFit1files %in% c("gammaFit1_ppcComb",
                                     "gammaFit1loo",
                                     "gammaFit1predMetrics",
                                     "gammaFit1pvalues",
                                     "gammaFit1stormsPredplot"))]

rm(list = gammaFit1filesRM)
rm(gammaFit1filesRM)

### LOGNORMAL ----
#### Model 1 ----
logNormalFit1 <- brm(
  bf(
    VMAX ~ log(HWRF) + eta,
    eta ~
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
      (1+StormElapsedTime|StormID),
    nl = TRUE
  ),
  data = StormdataTrain8, 
  family = lognormal(), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)

#save(logNormalFit1, file = "_data/logNormalFit1.RData")
prior_summary(logNormalFit1)
round(posterior_summary(logNormalFit1, probs = c(0.025, 0.975)))
logNormalFit1

print(logNormalFit1, digits = 4)
plot(logNormalFit1)
logNormalFit1ppcFit <- pp_check(logNormalFit1, ndraws = 100) + 
  labs(title = "logNormalFit1 Fit PPC") +
  theme_bw()
logNormalFit1ppcFit
logNormalFit1loo <- loo(logNormalFit1)
waic(logNormalFit1)
performance::check_distribution(logNormalFit1)
performance::check_outliers(logNormalFit1)
performance::check_heteroskedasticity(logNormalFit1)
performance_rmse(logNormalFit1)
performance_mae(logNormalFit1)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(logNormalFit1)


variance_decomposition(logNormalFit1)
exp(fixef(logNormalFit1))
ranef(logNormalFit1)

bayes_R2(logNormalFit1)

bayes_factor(logNormalFit1, logNormalFit1)
bayes_factor(logNormalFit1, gammaFit2)
bayes_factor(logNormalFit1, gammaFit3)
bayes_factor(logNormalFit1, studentFit1)
bayes_factor(logNormalFit1, studentFit2)
bayes_factor(logNormalFit1, studentFit3)
bayes_factor(logNormalFit1, linFit11)
bayes_factor(logNormalFit1, propFit1)
bayes_factor(logNormalFit1, logPropFit1)
loo(logNormalFit1, gammaFit3)

logNormalFit1smooths <- conditional_smooths(logNormalFit1)
plot(logNormalFit1smooths, stype = "raster", ask = FALSE)
logNormalFit1effects <- conditional_effects(logNormalFit1, 
                                        method = "posterior_predict",
                                        robust = FALSE,
                                        re_formula = NULL)
logNormalFit1effects <- conditional_effects(logNormalFit1)
plot(logNormalFit1effects, points = TRUE, ask = FALSE)

##### Prediction ----
## Fitted
logNormalFit1finalFit <- posterior_predict(logNormalFit1)
logNormalFit1finalFitMean <- colMeans(logNormalFit1finalFit)
logNormalFit1finalFitMed <- apply(logNormalFit1finalFit, 2, function(x){quantile(x, 0.5)})
logNormalFit1finalFitLCB <- apply(logNormalFit1finalFit, 2, function(x){quantile(x, 0.025)})
logNormalFit1finalFitUCB <- apply(logNormalFit1finalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
logNormalFit1finalPreds <- posterior_predict(logNormalFit1, 
                                         newdata = StormdataTest8,
                                         allow_new_levels = TRUE)
logNormalFit1finalPredsMean <- colMeans(logNormalFit1finalPreds)
logNormalFit1finalPredsMed <- apply(logNormalFit1finalPreds, 2, function(x){quantile(x, 0.5)})
logNormalFit1finalPredsLCB <- apply(logNormalFit1finalPreds, 2, function(x){quantile(x, 0.025)})
logNormalFit1finalPredsUCB <- apply(logNormalFit1finalPreds, 2, function(x){quantile(x, 0.975)})

logNormalFit1predMetrics <- tibble(
  Fit = "logNormalFit1",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(logNormalFit1finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(logNormalFit1finalFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < logNormalFit1finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(logNormalFit1finalPredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(logNormalFit1finalPredsMed - Actual_Yvec)),
  COV_pred = mean(logNormalFit1finalPredsLCB < Actual_Yvec & Actual_Yvec < logNormalFit1finalPredsUCB)
)
logNormalFit1predMetrics

##### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = logNormalFit1finalPreds) +
  labs(title = "logNormalFit1 Predict") +
  theme_bw()

logNormalFit1FitDF <- bind_cols(
  StormdataTrain3,
  LCB = logNormalFit1finalFitLCB,
  Mean = logNormalFit1finalFitMean,
  Med = logNormalFit1finalFitMed,
  UCB = logNormalFit1finalFitUCB
)

logNormalFit1stormsFitplot <- ggplot(data = logNormalFit1FitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))
logNormalFit1stormsFitplot

## Prediction
logNormalFit1PredDF <- bind_cols(
  StormdataTest3,
  LCB = logNormalFit1finalPredsLCB,
  Mean = logNormalFit1finalPredsMean,
  Med = logNormalFit1finalPredsMed,
  UCB = logNormalFit1finalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

logNormalFit1stormsPredplot <- ggplot(data = logNormalFit1PredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX, color = "Observed")) +
  geom_line(aes(y = Mean, color = "PPD Mean")) +
  facet_wrap(vars(StormID), ncol = 6)+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50)) +
  labs(title = "logNormalFit1 PPD Mean vs Observed VMAX",
       subtitle = "95% Credible Interval about PPD Mean") +
  scale_color_manual(name = NULL, values = c("black","red")) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1))
  ) +
  theme_bw()
logNormalFit1stormsPredplot

##### PPC ----
###### Quantile 2.5 
logNormalFit1LCBsims <- apply(logNormalFit1finalFit, 
                          MARGIN = 1,
                          function(x){
                            quantile(x, 0.025)
                          })
logNormalFit1LCBpvalueVec <- logNormalFit1LCBsims < quantile(StormdataTrain3$VMAX, 0.025)
logNormalFit1LCBpvalue <- sum(logNormalFit1LCBpvalueVec)
logNormalFit1LCBpvalue <- round(logNormalFit1LCBpvalue/4000, 3)
logNormalFit1LCBpvalue <- min(logNormalFit1LCBpvalue, 1 - logNormalFit1LCBpvalue)

logNormalFit1_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFit1finalFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", logNormalFit1LCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFit1_ppcLCB

###### Quantile 97.5 
logNormalFit1UCBsims <- apply(logNormalFit1finalFit, 
                          MARGIN = 1,
                          function(x){
                            quantile(x, 0.975)
                          })
logNormalFit1UCBpvalueVec <- logNormalFit1UCBsims < quantile(StormdataTrain3$VMAX, 0.975)
logNormalFit1UCBpvalue <- as.numeric(sum(logNormalFit1UCBpvalueVec))
logNormalFit1UCBpvalue <- round(logNormalFit1UCBpvalue/4000, 3)
logNormalFit1UCBpvalue <- min(logNormalFit1UCBpvalue, 1 - logNormalFit1UCBpvalue)

logNormalFit1_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFit1finalFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", logNormalFit1UCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFit1_ppcUCB

###### Mean 
logNormalFit1MEANsims <- apply(logNormalFit1finalFit, 
                           MARGIN = 1,
                           function(x){
                             mean(x)
                           })
logNormalFit1MEANpvalueVec <- logNormalFit1MEANsims < mean(StormdataTrain3$VMAX)
logNormalFit1MEANpvalue <- sum(logNormalFit1MEANpvalueVec)
logNormalFit1MEANpvalue <- round(logNormalFit1MEANpvalue/4000, 3)
logNormalFit1MEANpvalue <- min(logNormalFit1MEANpvalue, 1 - logNormalFit1MEANpvalue)

logNormalFit1_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFit1finalFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", logNormalFit1MEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFit1_ppcMEAN

###### Med 
logNormalFit1MEDsims <- apply(logNormalFit1finalFit, 
                          MARGIN = 1,
                          function(x){
                            quantile(x, 0.5)
                          })
logNormalFit1MEDpvalueVec <- logNormalFit1MEDsims < quantile(StormdataTrain3$VMAX, 0.5)
logNormalFit1MEDpvalue <- sum(logNormalFit1MEDpvalueVec)
logNormalFit1MEDpvalue <- round(logNormalFit1MEDpvalue/4000, 3)
logNormalFit1MEDpvalue <- min(logNormalFit1MEDpvalue, 1 - logNormalFit1MEDpvalue)

logNormalFit1_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFit1finalFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", logNormalFit1MEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFit1_ppcMED

###### SD 
logNormalFit1SDsims <- apply(logNormalFit1finalFit, 
                         MARGIN = 1,
                         function(x){
                           sd(x)
                         })
logNormalFit1SDpvalueVec <- logNormalFit1SDsims < sd(StormdataTrain3$VMAX)
logNormalFit1SDpvalue <- sum(logNormalFit1SDpvalueVec)
logNormalFit1SDpvalue <- round(logNormalFit1SDpvalue/4000, 3)
logNormalFit1SDpvalue <- min(logNormalFit1SDpvalue, 1 - logNormalFit1SDpvalue)

logNormalFit1_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFit1finalFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", logNormalFit1SDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFit1_ppcSD

###### Range 
logNormalFit1RANGEsims <- apply(logNormalFit1finalFit, 
                            MARGIN = 1,
                            function(x){
                              max(x)-min(x)
                            })
logNormalFit1RANGEpvalueVec <- logNormalFit1RANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
logNormalFit1RANGEpvalue <- sum(logNormalFit1RANGEpvalueVec)
logNormalFit1RANGEpvalue <- round(logNormalFit1RANGEpvalue/4000, 3)
logNormalFit1RANGEpvalue <- min(logNormalFit1RANGEpvalue, 1 - logNormalFit1RANGEpvalue)

logNormalFit1_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFit1finalFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", logNormalFit1RANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFit1_ppcRANGE

##### Combined Plot ----
logNormalFit1_ppcComb <- 
  logNormalFit1ppcFit /
  (logNormalFit1_ppcLCB | logNormalFit1_ppcMED | logNormalFit1_ppcUCB) /
  (logNormalFit1_ppcRANGE | logNormalFit1_ppcMEAN | logNormalFit1_ppcSD)
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
logNormalFit1_ppcComb

##### Bayes p-values ----
logNormalFit1pvalues <- tibble(
  Fit = "logNormalFit1",
  LCB = logNormalFit1LCBpvalue,
  Median = logNormalFit1MEDpvalue,
  UCB = logNormalFit1UCBpvalue,
  Range = logNormalFit1RANGEpvalue,
  Mean = logNormalFit1MEANpvalue,
  SD = logNormalFit1SDpvalue
)
logNormalFit1pvalues

logNormalFit1files <- ls()[str_detect(ls(), pattern = "logNormalFit1")]
logNormalFit1filesRM <- logNormalFit1files[!(logNormalFit1files %in% c("logNormalFit1_ppcComb",
                                                           "logNormalFit1loo",
                                                           "logNormalFit1predMetrics",
                                                           "logNormalFit1pvalues",
                                                           "logNormalFit1stormsPredplot"))]

rm(list = logNormalFit1filesRM)
rm(logNormalFit1filesRM)
rm(logNormalFit1files)

#### Model 2 ----
logNormalFit2 <- brm(
  bf(
    VMAX ~ 
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
      HWRF +
      (1|StormID)
  ),
  data = StormdataTrain7scale, 
  family = lognormal(), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)
#save(logNormalFit2, file = "_data/logNormalFit2.RData")
prior_summary(logNormalFit2)
round(posterior_summary(logNormalFit2, probs = c(0.025, 0.975)))
logNormalFit2

print(logNormalFit2, digits = 4)
plot(logNormalFit2)
logNormalFit2ppcFit <- pp_check(logNormalFit2, ndraws = 100) + 
  labs(title = "logNormalFit2 Fit PPC") +
  theme_bw()
logNormalFit2ppcFit
logNormalFit2loo <- loo(logNormalFit2)
waic(logNormalFit2)
performance::check_distribution(logNormalFit2)
performance::check_outliers(logNormalFit2)
performance::check_heteroskedasticity(logNormalFit2)
performance_rmse(logNormalFit2)
performance_mae(logNormalFit2)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(logNormalFit2)


variance_decomposition(logNormalFit2)
exp(fixef(logNormalFit2))
ranef(logNormalFit2)

bayes_R2(logNormalFit2)

bayes_factor(logNormalFit2, logNormalFit2)
bayes_factor(logNormalFit2, gammaFit2)
bayes_factor(logNormalFit2, gammaFit3)
bayes_factor(logNormalFit2, studentFit1)
bayes_factor(logNormalFit2, studentFit2)
bayes_factor(logNormalFit2, studentFit3)
bayes_factor(logNormalFit2, linFit11)
bayes_factor(logNormalFit2, propFit1)
bayes_factor(logNormalFit2, logPropFit1)
loo(logNormalFit2, gammaFit3)

logNormalFit2smooths <- conditional_smooths(logNormalFit2)
plot(logNormalFit2smooths, stype = "raster", ask = FALSE)
logNormalFit2effects <- conditional_effects(logNormalFit2, 
                                            method = "posterior_predict",
                                            robust = FALSE,
                                            re_formula = NULL)
logNormalFit2effects <- conditional_effects(logNormalFit2)
plot(logNormalFit2effects, points = TRUE, ask = FALSE)

##### Prediction ----
## Fitted
logNormalFit2finalFit <- posterior_predict(logNormalFit2)
logNormalFit2finalFitMean <- colMeans(logNormalFit2finalFit)
logNormalFit2finalFitMed <- apply(logNormalFit2finalFit, 2, function(x){quantile(x, 0.5)})
logNormalFit2finalFitLCB <- apply(logNormalFit2finalFit, 2, function(x){quantile(x, 0.025)})
logNormalFit2finalFitUCB <- apply(logNormalFit2finalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
logNormalFit2finalPreds <- posterior_predict(logNormalFit2, 
                                             newdata = StormdataTest8,
                                             allow_new_levels = TRUE)
logNormalFit2finalPredsMean <- colMeans(logNormalFit2finalPreds)
logNormalFit2finalPredsMed <- apply(logNormalFit2finalPreds, 2, function(x){quantile(x, 0.5)})
logNormalFit2finalPredsLCB <- apply(logNormalFit2finalPreds, 2, function(x){quantile(x, 0.025)})
logNormalFit2finalPredsUCB <- apply(logNormalFit2finalPreds, 2, function(x){quantile(x, 0.975)})

logNormalFit2predMetrics <- tibble(
  Fit = "logNormalFit2",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(logNormalFit2finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(logNormalFit2finalFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < logNormalFit2finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(logNormalFit2finalPredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(logNormalFit2finalPredsMed - Actual_Yvec)),
  COV_pred = mean(logNormalFit2finalPredsLCB < Actual_Yvec & Actual_Yvec < logNormalFit2finalPredsUCB)
)
logNormalFit2predMetrics

##### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = logNormalFit2finalPreds) +
  labs(title = "logNormalFit2 Predict") +
  theme_bw()

logNormalFit2FitDF <- bind_cols(
  StormdataTrain3,
  LCB = logNormalFit2finalFitLCB,
  Mean = logNormalFit2finalFitMean,
  Med = logNormalFit2finalFitMed,
  UCB = logNormalFit2finalFitUCB
)

logNormalFit2stormsFitplot <- ggplot(data = logNormalFit2FitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))
logNormalFit2stormsFitplot

## Prediction
logNormalFit2PredDF <- bind_cols(
  StormdataTest3,
  LCB = logNormalFit2finalPredsLCB,
  Mean = logNormalFit2finalPredsMean,
  Med = logNormalFit2finalPredsMed,
  UCB = logNormalFit2finalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

logNormalFit2stormsPredplot <- ggplot(data = logNormalFit2PredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX, color = "Observed")) +
  geom_line(aes(y = Mean, color = "PPD Mean")) +
  facet_wrap(vars(StormID), ncol = 6)+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50)) +
  labs(title = "logNormalFit2 PPD Mean vs Observed VMAX",
       subtitle = "95% Credible Interval about PPD Mean") +
  scale_color_manual(name = NULL, values = c("black","red")) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1))
  ) +
  theme_bw()
logNormalFit2stormsPredplot

##### PPC ----
###### Quantile 2.5 
logNormalFit2LCBsims <- apply(logNormalFit2finalFit, 
                              MARGIN = 1,
                              function(x){
                                quantile(x, 0.025)
                              })
logNormalFit2LCBpvalueVec <- logNormalFit2LCBsims < quantile(StormdataTrain3$VMAX, 0.025)
logNormalFit2LCBpvalue <- sum(logNormalFit2LCBpvalueVec)
logNormalFit2LCBpvalue <- round(logNormalFit2LCBpvalue/4000, 3)
logNormalFit2LCBpvalue <- min(logNormalFit2LCBpvalue, 1 - logNormalFit2LCBpvalue)

logNormalFit2_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFit2finalFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", logNormalFit2LCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFit2_ppcLCB

###### Quantile 97.5 
logNormalFit2UCBsims <- apply(logNormalFit2finalFit, 
                              MARGIN = 1,
                              function(x){
                                quantile(x, 0.975)
                              })
logNormalFit2UCBpvalueVec <- logNormalFit2UCBsims < quantile(StormdataTrain3$VMAX, 0.975)
logNormalFit2UCBpvalue <- as.numeric(sum(logNormalFit2UCBpvalueVec))
logNormalFit2UCBpvalue <- round(logNormalFit2UCBpvalue/4000, 3)
logNormalFit2UCBpvalue <- min(logNormalFit2UCBpvalue, 1 - logNormalFit2UCBpvalue)

logNormalFit2_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFit2finalFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", logNormalFit2UCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFit2_ppcUCB

###### Mean 
logNormalFit2MEANsims <- apply(logNormalFit2finalFit, 
                               MARGIN = 1,
                               function(x){
                                 mean(x)
                               })
logNormalFit2MEANpvalueVec <- logNormalFit2MEANsims < mean(StormdataTrain3$VMAX)
logNormalFit2MEANpvalue <- sum(logNormalFit2MEANpvalueVec)
logNormalFit2MEANpvalue <- round(logNormalFit2MEANpvalue/4000, 3)
logNormalFit2MEANpvalue <- min(logNormalFit2MEANpvalue, 1 - logNormalFit2MEANpvalue)

logNormalFit2_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFit2finalFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", logNormalFit2MEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFit2_ppcMEAN

###### Med 
logNormalFit2MEDsims <- apply(logNormalFit2finalFit, 
                              MARGIN = 1,
                              function(x){
                                quantile(x, 0.5)
                              })
logNormalFit2MEDpvalueVec <- logNormalFit2MEDsims < quantile(StormdataTrain3$VMAX, 0.5)
logNormalFit2MEDpvalue <- sum(logNormalFit2MEDpvalueVec)
logNormalFit2MEDpvalue <- round(logNormalFit2MEDpvalue/4000, 3)
logNormalFit2MEDpvalue <- min(logNormalFit2MEDpvalue, 1 - logNormalFit2MEDpvalue)

logNormalFit2_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFit2finalFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", logNormalFit2MEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFit2_ppcMED

###### SD 
logNormalFit2SDsims <- apply(logNormalFit2finalFit, 
                             MARGIN = 1,
                             function(x){
                               sd(x)
                             })
logNormalFit2SDpvalueVec <- logNormalFit2SDsims < sd(StormdataTrain3$VMAX)
logNormalFit2SDpvalue <- sum(logNormalFit2SDpvalueVec)
logNormalFit2SDpvalue <- round(logNormalFit2SDpvalue/4000, 3)
logNormalFit2SDpvalue <- min(logNormalFit2SDpvalue, 1 - logNormalFit2SDpvalue)

logNormalFit2_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFit2finalFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", logNormalFit2SDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFit2_ppcSD

###### Range 
logNormalFit2RANGEsims <- apply(logNormalFit2finalFit, 
                                MARGIN = 1,
                                function(x){
                                  max(x)-min(x)
                                })
logNormalFit2RANGEpvalueVec <- logNormalFit2RANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
logNormalFit2RANGEpvalue <- sum(logNormalFit2RANGEpvalueVec)
logNormalFit2RANGEpvalue <- round(logNormalFit2RANGEpvalue/4000, 3)
logNormalFit2RANGEpvalue <- min(logNormalFit2RANGEpvalue, 1 - logNormalFit2RANGEpvalue)

logNormalFit2_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFit2finalFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", logNormalFit2RANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFit2_ppcRANGE

##### Combined Plot ----
logNormalFit2_ppcComb <- 
  logNormalFit2ppcFit /
  (logNormalFit2_ppcLCB | logNormalFit2_ppcMED | logNormalFit2_ppcUCB) /
  (logNormalFit2_ppcRANGE | logNormalFit2_ppcMEAN | logNormalFit2_ppcSD)
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
logNormalFit2_ppcComb

##### Bayes p-values ----
logNormalFit2pvalues <- tibble(
  Fit = "logNormalFit2",
  LCB = logNormalFit2LCBpvalue,
  Median = logNormalFit2MEDpvalue,
  UCB = logNormalFit2UCBpvalue,
  Range = logNormalFit2RANGEpvalue,
  Mean = logNormalFit2MEANpvalue,
  SD = logNormalFit2SDpvalue
)
logNormalFit2pvalues

logNormalFit2files <- ls()[str_detect(ls(), pattern = "logNormalFit2")]
logNormalFit2filesRM <- logNormalFit2files[!(logNormalFit2files %in% c("logNormalFit2_ppcComb",
                                                                       "logNormalFit2loo",
                                                                       "logNormalFit2predMetrics",
                                                                       "logNormalFit2pvalues",
                                                                       "logNormalFit2stormsPredplot"))]

rm(list = logNormalFit2filesRM)
rm(logNormalFit2filesRM)

#### Model 2A ----
logNormalFit2A <- brm(
  bf(
    VMAX ~ 
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
      log(HWRF) +
      (1+StormElapsedTime|StormID)
  ),
  data = StormdataTrain8, 
  family = lognormal(), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)
#save(logNormalFit2A, file = "_data/logNormalFit2A.RData")
prior_summary(logNormalFit2A)
round(posterior_summary(logNormalFit2A), 4)
logNormalFit2A

print(logNormalFit2A, digits = 4)
plot(logNormalFit2A)
logNormalFit2AppcFit <- pp_check(logNormalFit2A, ndraws = 100) + 
  labs(title = "logNormalFit2A Fit PPC") +
  theme_bw()
logNormalFit2AppcFit
logNormalFit2Aloo <- loo(logNormalFit2A)
waic(logNormalFit2A)
performance::check_distribution(logNormalFit2A)
performance::check_outliers(logNormalFit2A)
performance::check_heteroskedasticity(logNormalFit2A)
performance_rmse(logNormalFit2A)
performance_mae(logNormalFit2A)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(logNormalFit2A)
logNormalFit2AlogLik <- log_lik(logNormalFit2A)

variance_decomposition(logNormalFit2A)
exp(fixef(logNormalFit2A))
ranef(logNormalFit2A)

bayes_R2(logNormalFit2A)

bayes_factor(logNormalFit2A, logNormalFit2A)
bayes_factor(logNormalFit2A, gammaFit2)
bayes_factor(logNormalFit2A, gammaFit3)
bayes_factor(logNormalFit2A, studentFit1)
bayes_factor(logNormalFit2A, studentFit2)
bayes_factor(logNormalFit2A, studentFit3)
bayes_factor(logNormalFit2A, linFit11)
bayes_factor(logNormalFit2A, propFit1)
bayes_factor(logNormalFit2A, logPropFit1)
loo(logNormalFit2A, gammaFit3)

logNormalFit2Asmooths <- conditional_smooths(logNormalFit2A)
plot(logNormalFit2Asmooths, stype = "raster", ask = FALSE)
logNormalFit2Aeffects <- conditional_effects(logNormalFit2A, 
                                            method = "posterior_predict",
                                            robust = FALSE,
                                            re_formula = NULL)
logNormalFit2Aeffects <- conditional_effects(logNormalFit2A)
plot(logNormalFit2Aeffects, points = TRUE, ask = FALSE)

##### Prediction ----
## Fitted
logNormalFit2AfinalFit <- posterior_predict(logNormalFit2A)
logNormalFit2AfinalFitMean <- colMeans(logNormalFit2AfinalFit)
logNormalFit2AfinalFitMed <- apply(logNormalFit2AfinalFit, 2, function(x){quantile(x, 0.5)})
logNormalFit2AfinalFitLCB <- apply(logNormalFit2AfinalFit, 2, function(x){quantile(x, 0.025)})
logNormalFit2AfinalFitUCB <- apply(logNormalFit2AfinalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
logNormalFit2AfinalPreds <- posterior_predict(logNormalFit2A, 
                                             newdata = StormdataTest8,
                                             allow_new_levels = TRUE)
logNormalFit2AfinalPredsMean <- colMeans(logNormalFit2AfinalPreds)
logNormalFit2AfinalPredsMed <- apply(logNormalFit2AfinalPreds, 2, function(x){quantile(x, 0.5)})
logNormalFit2AfinalPredsLCB <- apply(logNormalFit2AfinalPreds, 2, function(x){quantile(x, 0.025)})
logNormalFit2AfinalPredsUCB <- apply(logNormalFit2AfinalPreds, 2, function(x){quantile(x, 0.975)})

logNormalFit2ApredMetrics <- tibble(
  Fit = "logNormalFit2A",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(logNormalFit2AfinalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(logNormalFit2AfinalFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < logNormalFit2AfinalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(logNormalFit2AfinalPredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(logNormalFit2AfinalPredsMed - Actual_Yvec)),
  COV_pred = mean(logNormalFit2AfinalPredsLCB < Actual_Yvec & Actual_Yvec < logNormalFit2AfinalPredsUCB)
)
logNormalFit2ApredMetrics

##### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = logNormalFit2AfinalPreds) +
  labs(title = "logNormalFit2A Predict") +
  theme_bw()

logNormalFit2AFitDF <- bind_cols(
  StormdataTrain3,
  LCB = logNormalFit2AfinalFitLCB,
  Mean = logNormalFit2AfinalFitMean,
  Med = logNormalFit2AfinalFitMed,
  UCB = logNormalFit2AfinalFitUCB
)

logNormalFit2AstormsFitplot <- ggplot(data = logNormalFit2AFitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))
logNormalFit2AstormsFitplot

## Prediction
logNormalFit2APredDF <- bind_cols(
  StormdataTest3,
  LCB = logNormalFit2AfinalPredsLCB,
  Mean = logNormalFit2AfinalPredsMean,
  Med = logNormalFit2AfinalPredsMed,
  UCB = logNormalFit2AfinalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

logNormalFit2AstormsPredplot <- ggplot(data = logNormalFit2APredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX, color = "Observed")) +
  geom_line(aes(y = Mean, color = "PPD Mean")) +
  facet_wrap(vars(StormID), ncol = 6)+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50)) +
  labs(title = "logNormalFit2A PPD Mean vs Observed VMAX",
       subtitle = "95% Credible Interval about PPD Mean") +
  scale_color_manual(name = NULL, values = c("black","red")) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1))
  ) +
  theme_bw()
logNormalFit2AstormsPredplot

##### PPC ----
###### Quantile 2.5 
logNormalFit2ALCBsims <- apply(logNormalFit2AfinalFit, 
                              MARGIN = 1,
                              function(x){
                                quantile(x, 0.025)
                              })
logNormalFit2ALCBpvalueVec <- logNormalFit2ALCBsims < quantile(StormdataTrain3$VMAX, 0.025)
logNormalFit2ALCBpvalue <- sum(logNormalFit2ALCBpvalueVec)
logNormalFit2ALCBpvalue <- round(logNormalFit2ALCBpvalue/4000, 3)
logNormalFit2ALCBpvalue <- min(logNormalFit2ALCBpvalue, 1 - logNormalFit2ALCBpvalue)

logNormalFit2A_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFit2AfinalFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", logNormalFit2ALCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFit2A_ppcLCB

###### Quantile 97.5 
logNormalFit2AUCBsims <- apply(logNormalFit2AfinalFit, 
                              MARGIN = 1,
                              function(x){
                                quantile(x, 0.975)
                              })
logNormalFit2AUCBpvalueVec <- logNormalFit2AUCBsims < quantile(StormdataTrain3$VMAX, 0.975)
logNormalFit2AUCBpvalue <- as.numeric(sum(logNormalFit2AUCBpvalueVec))
logNormalFit2AUCBpvalue <- round(logNormalFit2AUCBpvalue/4000, 3)
logNormalFit2AUCBpvalue <- min(logNormalFit2AUCBpvalue, 1 - logNormalFit2AUCBpvalue)

logNormalFit2A_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFit2AfinalFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", logNormalFit2AUCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFit2A_ppcUCB

###### Mean 
logNormalFit2AMEANsims <- apply(logNormalFit2AfinalFit, 
                               MARGIN = 1,
                               function(x){
                                 mean(x)
                               })
logNormalFit2AMEANpvalueVec <- logNormalFit2AMEANsims < mean(StormdataTrain3$VMAX)
logNormalFit2AMEANpvalue <- sum(logNormalFit2AMEANpvalueVec)
logNormalFit2AMEANpvalue <- round(logNormalFit2AMEANpvalue/4000, 3)
logNormalFit2AMEANpvalue <- min(logNormalFit2AMEANpvalue, 1 - logNormalFit2AMEANpvalue)

logNormalFit2A_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFit2AfinalFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", logNormalFit2AMEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFit2A_ppcMEAN

###### Med 
logNormalFit2AMEDsims <- apply(logNormalFit2AfinalFit, 
                              MARGIN = 1,
                              function(x){
                                quantile(x, 0.5)
                              })
logNormalFit2AMEDpvalueVec <- logNormalFit2AMEDsims < quantile(StormdataTrain3$VMAX, 0.5)
logNormalFit2AMEDpvalue <- sum(logNormalFit2AMEDpvalueVec)
logNormalFit2AMEDpvalue <- round(logNormalFit2AMEDpvalue/4000, 3)
logNormalFit2AMEDpvalue <- min(logNormalFit2AMEDpvalue, 1 - logNormalFit2AMEDpvalue)

logNormalFit2A_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFit2AfinalFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", logNormalFit2AMEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFit2A_ppcMED

###### SD 
logNormalFit2ASDsims <- apply(logNormalFit2AfinalFit, 
                             MARGIN = 1,
                             function(x){
                               sd(x)
                             })
logNormalFit2ASDpvalueVec <- logNormalFit2ASDsims < sd(StormdataTrain3$VMAX)
logNormalFit2ASDpvalue <- sum(logNormalFit2ASDpvalueVec)
logNormalFit2ASDpvalue <- round(logNormalFit2ASDpvalue/4000, 3)
logNormalFit2ASDpvalue <- min(logNormalFit2ASDpvalue, 1 - logNormalFit2ASDpvalue)

logNormalFit2A_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFit2AfinalFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", logNormalFit2ASDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFit2A_ppcSD

###### Range 
logNormalFit2ARANGEsims <- apply(logNormalFit2AfinalFit, 
                                MARGIN = 1,
                                function(x){
                                  max(x)-min(x)
                                })
logNormalFit2ARANGEpvalueVec <- logNormalFit2ARANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
logNormalFit2ARANGEpvalue <- sum(logNormalFit2ARANGEpvalueVec)
logNormalFit2ARANGEpvalue <- round(logNormalFit2ARANGEpvalue/4000, 3)
logNormalFit2ARANGEpvalue <- min(logNormalFit2ARANGEpvalue, 1 - logNormalFit2ARANGEpvalue)

logNormalFit2A_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFit2AfinalFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", logNormalFit2ARANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFit2A_ppcRANGE

##### Combined Plot ----
logNormalFit2A_ppcComb <- 
  logNormalFit2AppcFit /
  (logNormalFit2A_ppcLCB | logNormalFit2A_ppcMED | logNormalFit2A_ppcUCB) /
  (logNormalFit2A_ppcRANGE | logNormalFit2A_ppcMEAN | logNormalFit2A_ppcSD)
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
logNormalFit2A_ppcComb

##### Bayes p-values ----
logNormalFit2Apvalues <- tibble(
  Fit = "logNormalFit2A",
  LCB = logNormalFit2ALCBpvalue,
  Median = logNormalFit2AMEDpvalue,
  UCB = logNormalFit2AUCBpvalue,
  Range = logNormalFit2ARANGEpvalue,
  Mean = logNormalFit2AMEANpvalue,
  SD = logNormalFit2ASDpvalue
)
logNormalFit2Apvalues

logNormalFit2Afiles <- ls()[str_detect(ls(), pattern = "logNormalFit2A")]
logNormalFit2AfilesRM <- logNormalFit2Afiles[!(logNormalFit2Afiles %in% c("logNormalFit2A_ppcComb",
                                                                       "logNormalFit2Aloo",
                                                                       "logNormalFit2ApredMetrics",
                                                                       "logNormalFit2Apvalues",
                                                                       "logNormalFit2AstormsPredplot"))]

rm(list = logNormalFit2AfilesRM)
rm(logNormalFit2AfilesRM)
rm(logNormalFit2Afiles)

#### Model 2B ----
logNormalFit2B <- brm(
  bf(
    VMAX ~ 
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
      log(HWRF) +
      (1|StormID)
  ),
  data = StormdataTrain8, 
  family = lognormal(), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)
#save(logNormalFit2B, file = "_data/logNormalFit2B.RData")
prior_summary(logNormalFit2B)
round(posterior_summary(logNormalFit2B, probs = c(0.025, 0.975)))
logNormalFit2B

print(logNormalFit2B, digits = 4)
plot(logNormalFit2B)
logNormalFit2BppcFit <- pp_check(logNormalFit2B, ndraws = 100) + 
  labs(title = "logNormalFit2B Fit PPC") +
  theme_bw()
logNormalFit2BppcFit
logNormalFit2Bloo <- loo(logNormalFit2B)
waic(logNormalFit2B)
performance::check_distribution(logNormalFit2B)
performance::check_outliers(logNormalFit2B)
performance::check_heteroskedasticity(logNormalFit2B)
performance_rmse(logNormalFit2B)
performance_mae(logNormalFit2B)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(logNormalFit2B)


variance_decomposition(logNormalFit2B)
exp(fixef(logNormalFit2B))
ranef(logNormalFit2B)

bayes_R2(logNormalFit2B)

bayes_factor(logNormalFit2B, logNormalFit2B)
bayes_factor(logNormalFit2B, gammaFit2)
bayes_factor(logNormalFit2B, gammaFit3)
bayes_factor(logNormalFit2B, studentFit1)
bayes_factor(logNormalFit2B, studentFit2)
bayes_factor(logNormalFit2B, studentFit3)
bayes_factor(logNormalFit2B, linFit11)
bayes_factor(logNormalFit2B, propFit1)
bayes_factor(logNormalFit2B, logPropFit1)
loo(logNormalFit2B, gammaFit3)

logNormalFit2Bsmooths <- conditional_smooths(logNormalFit2B)
plot(logNormalFit2Bsmooths, stype = "raster", ask = FALSE)
logNormalFit2Beffects <- conditional_effects(logNormalFit2B, 
                                             method = "posterior_predict",
                                             robust = FALSE,
                                             re_formula = NULL)
logNormalFit2Beffects <- conditional_effects(logNormalFit2B)
plot(logNormalFit2Beffects, points = TRUE, ask = FALSE)

##### Prediction ----
## Fitted
logNormalFit2BfinalFit <- posterior_predict(logNormalFit2B)
logNormalFit2BfinalFitMean <- colMeans(logNormalFit2BfinalFit)
logNormalFit2BfinalFitMed <- apply(logNormalFit2BfinalFit, 2, function(x){quantile(x, 0.5)})
logNormalFit2BfinalFitLCB <- apply(logNormalFit2BfinalFit, 2, function(x){quantile(x, 0.025)})
logNormalFit2BfinalFitUCB <- apply(logNormalFit2BfinalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
logNormalFit2BfinalPreds <- posterior_predict(logNormalFit2B, 
                                              newdata = StormdataTest8,
                                              allow_new_levels = TRUE)
logNormalFit2BfinalPredsMean <- colMeans(logNormalFit2BfinalPreds)
logNormalFit2BfinalPredsMed <- apply(logNormalFit2BfinalPreds, 2, function(x){quantile(x, 0.5)})
logNormalFit2BfinalPredsLCB <- apply(logNormalFit2BfinalPreds, 2, function(x){quantile(x, 0.025)})
logNormalFit2BfinalPredsUCB <- apply(logNormalFit2BfinalPreds, 2, function(x){quantile(x, 0.975)})

logNormalFit2BpredMetrics <- tibble(
  Fit = "logNormalFit2B",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(logNormalFit2BfinalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(logNormalFit2BfinalFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < logNormalFit2BfinalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(logNormalFit2BfinalPredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(logNormalFit2BfinalPredsMed - Actual_Yvec)),
  COV_pred = mean(logNormalFit2BfinalPredsLCB < Actual_Yvec & Actual_Yvec < logNormalFit2BfinalPredsUCB)
)
logNormalFit2BpredMetrics

##### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = logNormalFit2BfinalPreds) +
  labs(title = "logNormalFit2B Predict") +
  theme_bw()

logNormalFit2BFitDF <- bind_cols(
  StormdataTrain3,
  LCB = logNormalFit2BfinalFitLCB,
  Mean = logNormalFit2BfinalFitMean,
  Med = logNormalFit2BfinalFitMed,
  UCB = logNormalFit2BfinalFitUCB
)

logNormalFit2BstormsFitplot <- ggplot(data = logNormalFit2BFitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))
logNormalFit2BstormsFitplot

## Prediction
logNormalFit2BPredDF <- bind_cols(
  StormdataTest3,
  LCB = logNormalFit2BfinalPredsLCB,
  Mean = logNormalFit2BfinalPredsMean,
  Med = logNormalFit2BfinalPredsMed,
  UCB = logNormalFit2BfinalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

logNormalFit2BstormsPredplot <- ggplot(data = logNormalFit2BPredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX, color = "Observed")) +
  geom_line(aes(y = Mean, color = "PPD Mean")) +
  facet_wrap(vars(StormID), ncol = 6)+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50)) +
  labs(title = "logNormalFit2B PPD Mean vs Observed VMAX",
       subtitle = "95% Credible Interval about PPD Mean") +
  scale_color_manual(name = NULL, values = c("black","red")) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1))
  ) +
  theme_bw()
logNormalFit2BstormsPredplot

##### PPC ----
###### Quantile 2.5 
logNormalFit2BLCBsims <- apply(logNormalFit2BfinalFit, 
                               MARGIN = 1,
                               function(x){
                                 quantile(x, 0.025)
                               })
logNormalFit2BLCBpvalueVec <- logNormalFit2BLCBsims < quantile(StormdataTrain3$VMAX, 0.025)
logNormalFit2BLCBpvalue <- sum(logNormalFit2BLCBpvalueVec)
logNormalFit2BLCBpvalue <- round(logNormalFit2BLCBpvalue/4000, 3)
logNormalFit2BLCBpvalue <- min(logNormalFit2BLCBpvalue, 1 - logNormalFit2BLCBpvalue)

logNormalFit2B_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFit2BfinalFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", logNormalFit2BLCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFit2B_ppcLCB

###### Quantile 97.5 
logNormalFit2BUCBsims <- apply(logNormalFit2BfinalFit, 
                               MARGIN = 1,
                               function(x){
                                 quantile(x, 0.975)
                               })
logNormalFit2BUCBpvalueVec <- logNormalFit2BUCBsims < quantile(StormdataTrain3$VMAX, 0.975)
logNormalFit2BUCBpvalue <- as.numeric(sum(logNormalFit2BUCBpvalueVec))
logNormalFit2BUCBpvalue <- round(logNormalFit2BUCBpvalue/4000, 3)
logNormalFit2BUCBpvalue <- min(logNormalFit2BUCBpvalue, 1 - logNormalFit2BUCBpvalue)

logNormalFit2B_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFit2BfinalFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", logNormalFit2BUCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFit2B_ppcUCB

###### Mean 
logNormalFit2BMEANsims <- apply(logNormalFit2BfinalFit, 
                                MARGIN = 1,
                                function(x){
                                  mean(x)
                                })
logNormalFit2BMEANpvalueVec <- logNormalFit2BMEANsims < mean(StormdataTrain3$VMAX)
logNormalFit2BMEANpvalue <- sum(logNormalFit2BMEANpvalueVec)
logNormalFit2BMEANpvalue <- round(logNormalFit2BMEANpvalue/4000, 3)
logNormalFit2BMEANpvalue <- min(logNormalFit2BMEANpvalue, 1 - logNormalFit2BMEANpvalue)

logNormalFit2B_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFit2BfinalFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", logNormalFit2BMEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFit2B_ppcMEAN

###### Med 
logNormalFit2BMEDsims <- apply(logNormalFit2BfinalFit, 
                               MARGIN = 1,
                               function(x){
                                 quantile(x, 0.5)
                               })
logNormalFit2BMEDpvalueVec <- logNormalFit2BMEDsims < quantile(StormdataTrain3$VMAX, 0.5)
logNormalFit2BMEDpvalue <- sum(logNormalFit2BMEDpvalueVec)
logNormalFit2BMEDpvalue <- round(logNormalFit2BMEDpvalue/4000, 3)
logNormalFit2BMEDpvalue <- min(logNormalFit2BMEDpvalue, 1 - logNormalFit2BMEDpvalue)

logNormalFit2B_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFit2BfinalFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", logNormalFit2BMEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFit2B_ppcMED

###### SD 
logNormalFit2BSDsims <- apply(logNormalFit2BfinalFit, 
                              MARGIN = 1,
                              function(x){
                                sd(x)
                              })
logNormalFit2BSDpvalueVec <- logNormalFit2BSDsims < sd(StormdataTrain3$VMAX)
logNormalFit2BSDpvalue <- sum(logNormalFit2BSDpvalueVec)
logNormalFit2BSDpvalue <- round(logNormalFit2BSDpvalue/4000, 3)
logNormalFit2BSDpvalue <- min(logNormalFit2BSDpvalue, 1 - logNormalFit2BSDpvalue)

logNormalFit2B_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFit2BfinalFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", logNormalFit2BSDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFit2B_ppcSD

###### Range 
logNormalFit2BRANGEsims <- apply(logNormalFit2BfinalFit, 
                                 MARGIN = 1,
                                 function(x){
                                   max(x)-min(x)
                                 })
logNormalFit2BRANGEpvalueVec <- logNormalFit2BRANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
logNormalFit2BRANGEpvalue <- sum(logNormalFit2BRANGEpvalueVec)
logNormalFit2BRANGEpvalue <- round(logNormalFit2BRANGEpvalue/4000, 3)
logNormalFit2BRANGEpvalue <- min(logNormalFit2BRANGEpvalue, 1 - logNormalFit2BRANGEpvalue)

logNormalFit2B_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFit2BfinalFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", logNormalFit2BRANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFit2B_ppcRANGE

##### Combined Plot ----
logNormalFit2B_ppcComb <- 
  logNormalFit2BppcFit /
  (logNormalFit2B_ppcLCB | logNormalFit2B_ppcMED | logNormalFit2B_ppcUCB) /
  (logNormalFit2B_ppcRANGE | logNormalFit2B_ppcMEAN | logNormalFit2B_ppcSD)
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
logNormalFit2B_ppcComb

##### Bayes p-values ----
logNormalFit2Bpvalues <- tibble(
  Fit = "logNormalFit2B",
  LCB = logNormalFit2BLCBpvalue,
  Median = logNormalFit2BMEDpvalue,
  UCB = logNormalFit2BUCBpvalue,
  Range = logNormalFit2BRANGEpvalue,
  Mean = logNormalFit2BMEANpvalue,
  SD = logNormalFit2BSDpvalue
)
logNormalFit2Bpvalues

logNormalFit2Bfiles <- ls()[str_detect(ls(), pattern = "logNormalFit2B")]
logNormalFit2BfilesRM <- logNormalFit2Bfiles[!(logNormalFit2Bfiles %in% c("logNormalFit2B_ppcComb",
                                                                          "logNormalFit2Bloo",
                                                                          "logNormalFit2BpredMetrics",
                                                                          "logNormalFit2Bpvalues",
                                                                          "logNormalFit2BstormsPredplot"))]

rm(list = logNormalFit2BfilesRM)
rm(logNormalFit2BfilesRM)
rm(logNormalFit2Bfiles)


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
#save(propLinFit1, file = "_data/propLinFit1.RData")
prior_summary(propLinFit1)
round(posterior_summary(propLinFit1, probs = c(0.025, 0.975)))
propLinFit1

print(propLinFit1, digits = 4)
plot(propLinFit1)
propLinFit1ppcFit <- pp_check(propLinFit1, ndraws = 100) + 
  labs(title = "propLinFit1 Fit PPC") +
  theme_bw()
propLinFit1ppcFit
propLinFit1loo <- loo(propLinFit1)
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
bayes_factor(propLinFit1, gammaFit2)
bayes_factor(propLinFit1, gammaFit3)
bayes_factor(propLinFit1, studentFit1)
bayes_factor(propLinFit1, studentFit2)
bayes_factor(propLinFit1, studentFit3)
bayes_factor(propLinFit1, linFit11)
bayes_factor(propLinFit1, propFit1)
bayes_factor(propLinFit1, logPropFit1)
loo(propLinFit1, gammaFit3)

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
propLinFit1finalFit <- t(t(propLinFit1finalFit)*StormdataTrain3$HWRF)
propLinFit1finalFitMean <- colMeans(propLinFit1finalFit)
propLinFit1finalFitMed <- apply(propLinFit1finalFit, 2, function(x){quantile(x, 0.5)})
propLinFit1finalFitLCB <- apply(propLinFit1finalFit, 2, function(x){quantile(x, 0.025)})
propLinFit1finalFitUCB <- apply(propLinFit1finalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
propLinFit1finalPreds <- posterior_predict(propLinFit1, 
                                              newdata = StormdataTest8,
                                              allow_new_levels = TRUE)
propLinFit1finalPreds <- t(t(propLinFit1finalPreds)*StormdataTest3$HWRF)
propLinFit1finalPredsMean <- colMeans(propLinFit1finalPreds)
propLinFit1finalPredsMed <- apply(propLinFit1finalPreds, 2, function(x){quantile(x, 0.5)})
propLinFit1finalPredsLCB <- apply(propLinFit1finalPreds, 2, function(x){quantile(x, 0.025)})
propLinFit1finalPredsUCB <- apply(propLinFit1finalPreds, 2, function(x){quantile(x, 0.975)})

propLinFit1predMetrics <- tibble(
  Fit = "propLinFit1",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(propLinFit1finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(propLinFit1finalFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < propLinFit1finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(propLinFit1finalPredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(propLinFit1finalPredsMed - Actual_Yvec)),
  COV_pred = mean(propLinFit1finalPredsLCB < Actual_Yvec & Actual_Yvec < propLinFit1finalPredsUCB)
)
propLinFit1predMetrics

##### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = propLinFit1finalPreds) +
  labs(title = "propLinFit1 Predict") +
  theme_bw()

propLinFit1FitDF <- bind_cols(
  StormdataTrain3,
  LCB = propLinFit1finalFitLCB,
  Mean = propLinFit1finalFitMean,
  Med = propLinFit1finalFitMed,
  UCB = propLinFit1finalFitUCB
)

propLinFit1stormsFitplot <- ggplot(data = propLinFit1FitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))
propLinFit1stormsFitplot

## Prediction
propLinFit1PredDF <- bind_cols(
  StormdataTest3,
  LCB = propLinFit1finalPredsLCB,
  Mean = propLinFit1finalPredsMean,
  Med = propLinFit1finalPredsMed,
  UCB = propLinFit1finalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

propLinFit1stormsPredplot <- ggplot(data = propLinFit1PredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX, color = "Observed")) +
  geom_line(aes(y = Mean, color = "PPD Mean")) +
  facet_wrap(vars(StormID), ncol = 6)+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50)) +
  labs(title = "propLinFit1 PPD Mean vs Observed VMAX",
       subtitle = "95% Credible Interval about PPD Mean") +
  scale_color_manual(name = NULL, values = c("black","red")) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1))
  ) +
  theme_bw()
propLinFit1stormsPredplot

##### PPC ----
###### Quantile 2.5 
propLinFit1LCBsims <- apply(propLinFit1finalFit, 
                               MARGIN = 1,
                               function(x){
                                 quantile(x, 0.025)
                               })
propLinFit1LCBpvalueVec <- propLinFit1LCBsims < quantile(StormdataTrain3$VMAX, 0.025)
propLinFit1LCBpvalue <- sum(propLinFit1LCBpvalueVec)
propLinFit1LCBpvalue <- round(propLinFit1LCBpvalue/4000, 3)
propLinFit1LCBpvalue <- min(propLinFit1LCBpvalue, 1 - propLinFit1LCBpvalue)

propLinFit1_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           propLinFit1finalFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", propLinFit1LCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#propLinFit1_ppcLCB

###### Quantile 97.5 
propLinFit1UCBsims <- apply(propLinFit1finalFit, 
                               MARGIN = 1,
                               function(x){
                                 quantile(x, 0.975)
                               })
propLinFit1UCBpvalueVec <- propLinFit1UCBsims < quantile(StormdataTrain3$VMAX, 0.975)
propLinFit1UCBpvalue <- as.numeric(sum(propLinFit1UCBpvalueVec))
propLinFit1UCBpvalue <- round(propLinFit1UCBpvalue/4000, 3)
propLinFit1UCBpvalue <- min(propLinFit1UCBpvalue, 1 - propLinFit1UCBpvalue)

propLinFit1_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           propLinFit1finalFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", propLinFit1UCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#propLinFit1_ppcUCB

###### Mean 
propLinFit1MEANsims <- apply(propLinFit1finalFit, 
                                MARGIN = 1,
                                function(x){
                                  mean(x)
                                })
propLinFit1MEANpvalueVec <- propLinFit1MEANsims < mean(StormdataTrain3$VMAX)
propLinFit1MEANpvalue <- sum(propLinFit1MEANpvalueVec)
propLinFit1MEANpvalue <- round(propLinFit1MEANpvalue/4000, 3)
propLinFit1MEANpvalue <- min(propLinFit1MEANpvalue, 1 - propLinFit1MEANpvalue)

propLinFit1_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           propLinFit1finalFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", propLinFit1MEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#propLinFit1_ppcMEAN

###### Med 
propLinFit1MEDsims <- apply(propLinFit1finalFit, 
                               MARGIN = 1,
                               function(x){
                                 quantile(x, 0.5)
                               })
propLinFit1MEDpvalueVec <- propLinFit1MEDsims < quantile(StormdataTrain3$VMAX, 0.5)
propLinFit1MEDpvalue <- sum(propLinFit1MEDpvalueVec)
propLinFit1MEDpvalue <- round(propLinFit1MEDpvalue/4000, 3)
propLinFit1MEDpvalue <- min(propLinFit1MEDpvalue, 1 - propLinFit1MEDpvalue)

propLinFit1_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           propLinFit1finalFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", propLinFit1MEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#propLinFit1_ppcMED

###### SD 
propLinFit1SDsims <- apply(propLinFit1finalFit, 
                              MARGIN = 1,
                              function(x){
                                sd(x)
                              })
propLinFit1SDpvalueVec <- propLinFit1SDsims < sd(StormdataTrain3$VMAX)
propLinFit1SDpvalue <- sum(propLinFit1SDpvalueVec)
propLinFit1SDpvalue <- round(propLinFit1SDpvalue/4000, 3)
propLinFit1SDpvalue <- min(propLinFit1SDpvalue, 1 - propLinFit1SDpvalue)

propLinFit1_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           propLinFit1finalFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", propLinFit1SDpvalue, ")")) +
  theme_bw() +
  legend_none()
#propLinFit1_ppcSD

###### Range 
propLinFit1RANGEsims <- apply(propLinFit1finalFit, 
                                 MARGIN = 1,
                                 function(x){
                                   max(x)-min(x)
                                 })
propLinFit1RANGEpvalueVec <- propLinFit1RANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
propLinFit1RANGEpvalue <- sum(propLinFit1RANGEpvalueVec)
propLinFit1RANGEpvalue <- round(propLinFit1RANGEpvalue/4000, 3)
propLinFit1RANGEpvalue <- min(propLinFit1RANGEpvalue, 1 - propLinFit1RANGEpvalue)

propLinFit1_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           propLinFit1finalFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", propLinFit1RANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#propLinFit1_ppcRANGE

##### Combined Plot ----
propLinFit1_ppcComb <- 
  propLinFit1ppcFit /
  (propLinFit1_ppcLCB | propLinFit1_ppcMED | propLinFit1_ppcUCB) /
  (propLinFit1_ppcRANGE | propLinFit1_ppcMEAN | propLinFit1_ppcSD)
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
propLinFit1_ppcComb

##### Bayes p-values ----
propLinFit1pvalues <- tibble(
  Fit = "propLinFit1",
  LCB = propLinFit1LCBpvalue,
  Median = propLinFit1MEDpvalue,
  UCB = propLinFit1UCBpvalue,
  Range = propLinFit1RANGEpvalue,
  Mean = propLinFit1MEANpvalue,
  SD = propLinFit1SDpvalue
)
propLinFit1pvalues

propLinFit1files <- ls()[str_detect(ls(), pattern = "propLinFit1")]
propLinFit1filesRM <- propLinFit1files[!(propLinFit1files %in% c("propLinFit1_ppcComb",
                                                                          "propLinFit1loo",
                                                                          "propLinFit1predMetrics",
                                                                          "propLinFit1pvalues",
                                                                          "propLinFit1stormsPredplot"))]

rm(list = propLinFit1filesRM)
rm(propLinFit1filesRM)
rm(propLinFit1files)

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
#save(propStudentFit1, file = "_data/propStudentFit1.RData")
prior_summary(propStudentFit1)
round(posterior_summary(propStudentFit1, probs = c(0.025, 0.975)))
propStudentFit1

print(propStudentFit1, digits = 4)
plot(propStudentFit1)
propStudentFit1ppcFit <- pp_check(propStudentFit1, ndraws = 100) + 
  labs(title = "propStudentFit1 Fit PPC") +
  theme_bw()
propStudentFit1ppcFit
propStudentFit1loo <- loo(propStudentFit1)
waic(propStudentFit1)
performance::check_distribution(propStudentFit1)
performance::check_outliers(propStudentFit1)
performance::check_heteroskedasticity(propStudentFit1)
performance_rmse(propStudentFit1)
performance_mae(propStudentFit1)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(propStudentFit1)
propStudentFit1logLik <- log_lik(propStudentFit1)


variance_decomposition(propStudentFit1)
exp(fixef(propStudentFit1))
ranef(propStudentFit1)

bayes_R2(propStudentFit1)

bayes_factor(propStudentFit1, propStudentFit1)
bayes_factor(propStudentFit1, gammaFit2)
bayes_factor(propStudentFit1, gammaFit3)
bayes_factor(propStudentFit1, studentFit1)
bayes_factor(propStudentFit1, studentFit2)
bayes_factor(propStudentFit1, studentFit3)
bayes_factor(propStudentFit1, linFit11)
bayes_factor(propStudentFit1, propFit1)
bayes_factor(propStudentFit1, logPropFit1)
loo(propStudentFit1, gammaFit3)

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
propStudentFit1finalFit <- t(t(propStudentFit1finalFit)*StormdataTrain3$HWRF)
propStudentFit1finalFitMean <- colMeans(propStudentFit1finalFit)
propStudentFit1finalFitMed <- apply(propStudentFit1finalFit, 2, function(x){quantile(x, 0.5)})
propStudentFit1finalFitLCB <- apply(propStudentFit1finalFit, 2, function(x){quantile(x, 0.025)})
propStudentFit1finalFitUCB <- apply(propStudentFit1finalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
propStudentFit1finalPreds <- posterior_predict(propStudentFit1, 
                                           newdata = StormdataTest8,
                                           allow_new_levels = TRUE)
propStudentFit1finalPreds <- t(t(propStudentFit1finalPreds)*StormdataTest3$HWRF)
propStudentFit1finalPredsMean <- colMeans(propStudentFit1finalPreds)
propStudentFit1finalPredsMed <- apply(propStudentFit1finalPreds, 2, function(x){quantile(x, 0.5)})
propStudentFit1finalPredsLCB <- apply(propStudentFit1finalPreds, 2, function(x){quantile(x, 0.025)})
propStudentFit1finalPredsUCB <- apply(propStudentFit1finalPreds, 2, function(x){quantile(x, 0.975)})

propStudentFit1predMetrics <- tibble(
  Fit = "propStudentFit1",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(propStudentFit1finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(propStudentFit1finalFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < propStudentFit1finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(propStudentFit1finalPredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(propStudentFit1finalPredsMed - Actual_Yvec)),
  COV_pred = mean(propStudentFit1finalPredsLCB < Actual_Yvec & Actual_Yvec < propStudentFit1finalPredsUCB)
)
propStudentFit1predMetrics

##### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = propStudentFit1finalPreds) +
  labs(title = "propStudentFit1 Predict") +
  theme_bw()

propStudentFit1FitDF <- bind_cols(
  StormdataTrain3,
  LCB = propStudentFit1finalFitLCB,
  Mean = propStudentFit1finalFitMean,
  Med = propStudentFit1finalFitMed,
  UCB = propStudentFit1finalFitUCB
)

propStudentFit1stormsFitplot <- ggplot(data = propStudentFit1FitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))
propStudentFit1stormsFitplot

## Prediction
propStudentFit1PredDF <- bind_cols(
  StormdataTest3,
  LCB = propStudentFit1finalPredsLCB,
  Mean = propStudentFit1finalPredsMean,
  Med = propStudentFit1finalPredsMed,
  UCB = propStudentFit1finalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

propStudentFit1stormsPredplot <- ggplot(data = propStudentFit1PredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX, color = "Observed")) +
  geom_line(aes(y = Mean, color = "PPD Mean")) +
  facet_wrap(vars(StormID), ncol = 6)+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50)) +
  labs(title = "propStudentFit1 PPD Mean vs Observed VMAX",
       subtitle = "95% Credible Interval about PPD Mean") +
  scale_color_manual(name = NULL, values = c("black","red")) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1))
  ) +
  theme_bw()
propStudentFit1stormsPredplot

##### PPC ----
###### Quantile 2.5 
propStudentFit1LCBsims <- apply(propStudentFit1finalFit, 
                            MARGIN = 1,
                            function(x){
                              quantile(x, 0.025)
                            })
propStudentFit1LCBpvalueVec <- propStudentFit1LCBsims < quantile(StormdataTrain3$VMAX, 0.025)
propStudentFit1LCBpvalue <- sum(propStudentFit1LCBpvalueVec)
propStudentFit1LCBpvalue <- round(propStudentFit1LCBpvalue/4000, 3)
propStudentFit1LCBpvalue <- min(propStudentFit1LCBpvalue, 1 - propStudentFit1LCBpvalue)

propStudentFit1_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           propStudentFit1finalFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", propStudentFit1LCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#propStudentFit1_ppcLCB

###### Quantile 97.5 
propStudentFit1UCBsims <- apply(propStudentFit1finalFit, 
                            MARGIN = 1,
                            function(x){
                              quantile(x, 0.975)
                            })
propStudentFit1UCBpvalueVec <- propStudentFit1UCBsims < quantile(StormdataTrain3$VMAX, 0.975)
propStudentFit1UCBpvalue <- as.numeric(sum(propStudentFit1UCBpvalueVec))
propStudentFit1UCBpvalue <- round(propStudentFit1UCBpvalue/4000, 3)
propStudentFit1UCBpvalue <- min(propStudentFit1UCBpvalue, 1 - propStudentFit1UCBpvalue)

propStudentFit1_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           propStudentFit1finalFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", propStudentFit1UCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#propStudentFit1_ppcUCB

###### Mean 
propStudentFit1MEANsims <- apply(propStudentFit1finalFit, 
                             MARGIN = 1,
                             function(x){
                               mean(x)
                             })
propStudentFit1MEANpvalueVec <- propStudentFit1MEANsims < mean(StormdataTrain3$VMAX)
propStudentFit1MEANpvalue <- sum(propStudentFit1MEANpvalueVec)
propStudentFit1MEANpvalue <- round(propStudentFit1MEANpvalue/4000, 3)
propStudentFit1MEANpvalue <- min(propStudentFit1MEANpvalue, 1 - propStudentFit1MEANpvalue)

propStudentFit1_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           propStudentFit1finalFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", propStudentFit1MEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#propStudentFit1_ppcMEAN

###### Med 
propStudentFit1MEDsims <- apply(propStudentFit1finalFit, 
                            MARGIN = 1,
                            function(x){
                              quantile(x, 0.5)
                            })
propStudentFit1MEDpvalueVec <- propStudentFit1MEDsims < quantile(StormdataTrain3$VMAX, 0.5)
propStudentFit1MEDpvalue <- sum(propStudentFit1MEDpvalueVec)
propStudentFit1MEDpvalue <- round(propStudentFit1MEDpvalue/4000, 3)
propStudentFit1MEDpvalue <- min(propStudentFit1MEDpvalue, 1 - propStudentFit1MEDpvalue)

propStudentFit1_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           propStudentFit1finalFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", propStudentFit1MEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#propStudentFit1_ppcMED

###### SD 
propStudentFit1SDsims <- apply(propStudentFit1finalFit, 
                           MARGIN = 1,
                           function(x){
                             sd(x)
                           })
propStudentFit1SDpvalueVec <- propStudentFit1SDsims < sd(StormdataTrain3$VMAX)
propStudentFit1SDpvalue <- sum(propStudentFit1SDpvalueVec)
propStudentFit1SDpvalue <- round(propStudentFit1SDpvalue/4000, 3)
propStudentFit1SDpvalue <- min(propStudentFit1SDpvalue, 1 - propStudentFit1SDpvalue)

propStudentFit1_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           propStudentFit1finalFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", propStudentFit1SDpvalue, ")")) +
  theme_bw() +
  legend_none()
#propStudentFit1_ppcSD

###### Range 
propStudentFit1RANGEsims <- apply(propStudentFit1finalFit, 
                              MARGIN = 1,
                              function(x){
                                max(x)-min(x)
                              })
propStudentFit1RANGEpvalueVec <- propStudentFit1RANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
propStudentFit1RANGEpvalue <- sum(propStudentFit1RANGEpvalueVec)
propStudentFit1RANGEpvalue <- round(propStudentFit1RANGEpvalue/4000, 3)
propStudentFit1RANGEpvalue <- min(propStudentFit1RANGEpvalue, 1 - propStudentFit1RANGEpvalue)

propStudentFit1_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           propStudentFit1finalFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", propStudentFit1RANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#propStudentFit1_ppcRANGE

##### Combined Plot ----
propStudentFit1_ppcComb <- 
  propStudentFit1ppcFit /
  (propStudentFit1_ppcLCB | propStudentFit1_ppcMED | propStudentFit1_ppcUCB) /
  (propStudentFit1_ppcRANGE | propStudentFit1_ppcMEAN | propStudentFit1_ppcSD)
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
propStudentFit1_ppcComb

##### Bayes p-values ----
propStudentFit1pvalues <- tibble(
  Fit = "propStudentFit1",
  LCB = propStudentFit1LCBpvalue,
  Median = propStudentFit1MEDpvalue,
  UCB = propStudentFit1UCBpvalue,
  Range = propStudentFit1RANGEpvalue,
  Mean = propStudentFit1MEANpvalue,
  SD = propStudentFit1SDpvalue
)
propStudentFit1pvalues

propStudentFit1files <- ls()[str_detect(ls(), pattern = "propStudentFit1")]
propStudentFit1filesRM <- propStudentFit1files[!(propStudentFit1files %in% c("propStudentFit1_ppcComb",
                                                                 "propStudentFit1loo",
                                                                 "propStudentFit1predMetrics",
                                                                 "propStudentFit1pvalues",
                                                                 "propStudentFit1stormsPredplot"))]

rm(list = propStudentFit1filesRM)
rm(propStudentFit1filesRM)
rm(propStudentFit1files)

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
#save(propStudentFit2, file = "_data/propStudentFit2.RData")
prior_summary(propStudentFit2)
round(posterior_summary(propStudentFit2, probs = c(0.025, 0.975)))
propStudentFit2

print(propStudentFit2, digits = 4)
plot(propStudentFit2)
propStudentFit2ppcFit <- pp_check(propStudentFit2, ndraws = 100) + 
  labs(title = "propStudentFit2 Fit PPC") +
  theme_bw()
propStudentFit2ppcFit
propStudentFit2loo <- loo(propStudentFit2)
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

bayes_factor(propStudentFit2, propStudentFit2)
bayes_factor(propStudentFit2, gammaFit2)
bayes_factor(propStudentFit2, gammaFit3)
bayes_factor(propStudentFit2, studentFit1)
bayes_factor(propStudentFit2, studentFit2)
bayes_factor(propStudentFit2, studentFit3)
bayes_factor(propStudentFit2, linFit11)
bayes_factor(propStudentFit2, propFit1)
bayes_factor(propStudentFit2, logPropFit1)
loo(propStudentFit2, gammaFit3)

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
propStudentFit2finalFit <- t(t(propStudentFit2finalFit)*StormdataTrain3$HWRF)
propStudentFit2finalFitMean <- colMeans(propStudentFit2finalFit)
propStudentFit2finalFitMed <- apply(propStudentFit2finalFit, 2, function(x){quantile(x, 0.5)})
propStudentFit2finalFitLCB <- apply(propStudentFit2finalFit, 2, function(x){quantile(x, 0.025)})
propStudentFit2finalFitUCB <- apply(propStudentFit2finalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
propStudentFit2finalPreds <- posterior_predict(propStudentFit2, 
                                           newdata = StormdataTest8,
                                           allow_new_levels = TRUE)
propStudentFit2finalPreds <- t(t(propStudentFit2finalPreds)*StormdataTest3$HWRF)
propStudentFit2finalPredsMean <- colMeans(propStudentFit2finalPreds)
propStudentFit2finalPredsMed <- apply(propStudentFit2finalPreds, 2, function(x){quantile(x, 0.5)})
propStudentFit2finalPredsLCB <- apply(propStudentFit2finalPreds, 2, function(x){quantile(x, 0.025)})
propStudentFit2finalPredsUCB <- apply(propStudentFit2finalPreds, 2, function(x){quantile(x, 0.975)})

propStudentFit2predMetrics <- tibble(
  Fit = "propStudentFit2",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(propStudentFit2finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(propStudentFit2finalFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < propStudentFit2finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(propStudentFit2finalPredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(propStudentFit2finalPredsMed - Actual_Yvec)),
  COV_pred = mean(propStudentFit2finalPredsLCB < Actual_Yvec & Actual_Yvec < propStudentFit2finalPredsUCB)
)
propStudentFit2predMetrics

##### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = propStudentFit2finalPreds) +
  labs(title = "propStudentFit2 Predict") +
  theme_bw()

propStudentFit2FitDF <- bind_cols(
  StormdataTrain3,
  LCB = propStudentFit2finalFitLCB,
  Mean = propStudentFit2finalFitMean,
  Med = propStudentFit2finalFitMed,
  UCB = propStudentFit2finalFitUCB
)

propStudentFit2stormsFitplot <- ggplot(data = propStudentFit2FitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))
propStudentFit2stormsFitplot

## Prediction
propStudentFit2PredDF <- bind_cols(
  StormdataTest3,
  LCB = propStudentFit2finalPredsLCB,
  Mean = propStudentFit2finalPredsMean,
  Med = propStudentFit2finalPredsMed,
  UCB = propStudentFit2finalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

propStudentFit2stormsPredplot <- ggplot(data = propStudentFit2PredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX, color = "Observed")) +
  geom_line(aes(y = Mean, color = "PPD Mean")) +
  facet_wrap(vars(StormID), ncol = 6)+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50)) +
  labs(title = "propStudentFit2 PPD Mean vs Observed VMAX",
       subtitle = "95% Credible Interval about PPD Mean") +
  scale_color_manual(name = NULL, values = c("black","red")) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1))
  ) +
  theme_bw()
propStudentFit2stormsPredplot

##### PPC ----
###### Quantile 2.5 
propStudentFit2LCBsims <- apply(propStudentFit2finalFit, 
                            MARGIN = 1,
                            function(x){
                              quantile(x, 0.025)
                            })
propStudentFit2LCBpvalueVec <- propStudentFit2LCBsims < quantile(StormdataTrain3$VMAX, 0.025)
propStudentFit2LCBpvalue <- sum(propStudentFit2LCBpvalueVec)
propStudentFit2LCBpvalue <- round(propStudentFit2LCBpvalue/4000, 3)
propStudentFit2LCBpvalue <- min(propStudentFit2LCBpvalue, 1 - propStudentFit2LCBpvalue)

propStudentFit2_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           propStudentFit2finalFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", propStudentFit2LCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#propStudentFit2_ppcLCB

###### Quantile 97.5 
propStudentFit2UCBsims <- apply(propStudentFit2finalFit, 
                            MARGIN = 1,
                            function(x){
                              quantile(x, 0.975)
                            })
propStudentFit2UCBpvalueVec <- propStudentFit2UCBsims < quantile(StormdataTrain3$VMAX, 0.975)
propStudentFit2UCBpvalue <- as.numeric(sum(propStudentFit2UCBpvalueVec))
propStudentFit2UCBpvalue <- round(propStudentFit2UCBpvalue/4000, 3)
propStudentFit2UCBpvalue <- min(propStudentFit2UCBpvalue, 1 - propStudentFit2UCBpvalue)

propStudentFit2_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           propStudentFit2finalFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", propStudentFit2UCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#propStudentFit2_ppcUCB

###### Mean 
propStudentFit2MEANsims <- apply(propStudentFit2finalFit, 
                             MARGIN = 1,
                             function(x){
                               mean(x)
                             })
propStudentFit2MEANpvalueVec <- propStudentFit2MEANsims < mean(StormdataTrain3$VMAX)
propStudentFit2MEANpvalue <- sum(propStudentFit2MEANpvalueVec)
propStudentFit2MEANpvalue <- round(propStudentFit2MEANpvalue/4000, 3)
propStudentFit2MEANpvalue <- min(propStudentFit2MEANpvalue, 1 - propStudentFit2MEANpvalue)

propStudentFit2_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           propStudentFit2finalFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", propStudentFit2MEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#propStudentFit2_ppcMEAN

###### Med 
propStudentFit2MEDsims <- apply(propStudentFit2finalFit, 
                            MARGIN = 1,
                            function(x){
                              quantile(x, 0.5)
                            })
propStudentFit2MEDpvalueVec <- propStudentFit2MEDsims < quantile(StormdataTrain3$VMAX, 0.5)
propStudentFit2MEDpvalue <- sum(propStudentFit2MEDpvalueVec)
propStudentFit2MEDpvalue <- round(propStudentFit2MEDpvalue/4000, 3)
propStudentFit2MEDpvalue <- min(propStudentFit2MEDpvalue, 1 - propStudentFit2MEDpvalue)

propStudentFit2_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           propStudentFit2finalFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", propStudentFit2MEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#propStudentFit2_ppcMED

###### SD 
propStudentFit2SDsims <- apply(propStudentFit2finalFit, 
                           MARGIN = 1,
                           function(x){
                             sd(x)
                           })
propStudentFit2SDpvalueVec <- propStudentFit2SDsims < sd(StormdataTrain3$VMAX)
propStudentFit2SDpvalue <- sum(propStudentFit2SDpvalueVec)
propStudentFit2SDpvalue <- round(propStudentFit2SDpvalue/4000, 3)
propStudentFit2SDpvalue <- min(propStudentFit2SDpvalue, 1 - propStudentFit2SDpvalue)

propStudentFit2_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           propStudentFit2finalFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", propStudentFit2SDpvalue, ")")) +
  theme_bw() +
  legend_none()
#propStudentFit2_ppcSD

###### Range 
propStudentFit2RANGEsims <- apply(propStudentFit2finalFit, 
                              MARGIN = 1,
                              function(x){
                                max(x)-min(x)
                              })
propStudentFit2RANGEpvalueVec <- propStudentFit2RANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
propStudentFit2RANGEpvalue <- sum(propStudentFit2RANGEpvalueVec)
propStudentFit2RANGEpvalue <- round(propStudentFit2RANGEpvalue/4000, 3)
propStudentFit2RANGEpvalue <- min(propStudentFit2RANGEpvalue, 1 - propStudentFit2RANGEpvalue)

propStudentFit2_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           propStudentFit2finalFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", propStudentFit2RANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#propStudentFit2_ppcRANGE

##### Combined Plot ----
propStudentFit2_ppcComb <- 
  propStudentFit2ppcFit /
  (propStudentFit2_ppcLCB | propStudentFit2_ppcMED | propStudentFit2_ppcUCB) /
  (propStudentFit2_ppcRANGE | propStudentFit2_ppcMEAN | propStudentFit2_ppcSD)
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
propStudentFit2_ppcComb

##### Bayes p-values ----
propStudentFit2pvalues <- tibble(
  Fit = "propStudentFit2",
  LCB = propStudentFit2LCBpvalue,
  Median = propStudentFit2MEDpvalue,
  UCB = propStudentFit2UCBpvalue,
  Range = propStudentFit2RANGEpvalue,
  Mean = propStudentFit2MEANpvalue,
  SD = propStudentFit2SDpvalue
)
propStudentFit2pvalues

propStudentFit2files <- ls()[str_detect(ls(), pattern = "propStudentFit2")]
propStudentFit2filesRM <- propStudentFit2files[!(propStudentFit2files %in% c("propStudentFit2_ppcComb",
                                                                 "propStudentFit2loo",
                                                                 "propStudentFit2predMetrics",
                                                                 "propStudentFit2pvalues",
                                                                 "propStudentFit2stormsPredplot"))]

rm(list = propStudentFit2filesRM)
rm(propStudentFit2filesRM)
rm(propStudentFit2files)

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
                                               newdata = StormdataTest8,
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
                                               newdata = StormdataTest8,
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
#save(propStudentFit5, file = "_data/propStudentFit5.RData")
prior_summary(propStudentFit5)
round(posterior_summary(propStudentFit5, probs = c(0.025, 0.975)))
propStudentFit5

print(propStudentFit5, digits = 4)
plot(propStudentFit5)
propStudentFit5ppcFit <- pp_check(propStudentFit5, ndraws = 100) + 
  labs(title = "propStudentFit5 Fit PPC") +
  theme_bw()
propStudentFit5ppcFit
propStudentFit5loo <- loo(propStudentFit5)
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

bayes_factor(propStudentFit5, propStudentFit5)
bayes_factor(propStudentFit5, gammaFit2)
bayes_factor(propStudentFit5, gammaFit3)
bayes_factor(propStudentFit5, studentFit1)
bayes_factor(propStudentFit5, studentFit2)
bayes_factor(propStudentFit5, studentFit3)
bayes_factor(propStudentFit5, linFit11)
bayes_factor(propStudentFit5, propFit1)
bayes_factor(propStudentFit5, logPropFit1)
loo(propStudentFit5, gammaFit3)

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
propStudentFit5finalFit <- t(t(propStudentFit5finalFit)*StormdataTrain3$HWRF)
propStudentFit5finalFitMean <- colMeans(propStudentFit5finalFit)
propStudentFit5finalFitMed <- apply(propStudentFit5finalFit, 2, function(x){quantile(x, 0.5)})
propStudentFit5finalFitLCB <- apply(propStudentFit5finalFit, 2, function(x){quantile(x, 0.025)})
propStudentFit5finalFitUCB <- apply(propStudentFit5finalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
propStudentFit5finalPreds <- posterior_predict(propStudentFit5, 
                                           newdata = StormdataTest8,
                                           allow_new_levels = TRUE)
propStudentFit5finalPreds <- t(t(propStudentFit5finalPreds)*StormdataTest3$HWRF)
propStudentFit5finalPredsMean <- colMeans(propStudentFit5finalPreds)
propStudentFit5finalPredsMed <- apply(propStudentFit5finalPreds, 2, function(x){quantile(x, 0.5)})
propStudentFit5finalPredsLCB <- apply(propStudentFit5finalPreds, 2, function(x){quantile(x, 0.025)})
propStudentFit5finalPredsUCB <- apply(propStudentFit5finalPreds, 2, function(x){quantile(x, 0.975)})

propStudentFit5predMetrics <- tibble(
  Fit = "propStudentFit5",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(propStudentFit5finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(propStudentFit5finalFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < propStudentFit5finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(propStudentFit5finalPredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(propStudentFit5finalPredsMed - Actual_Yvec)),
  COV_pred = mean(propStudentFit5finalPredsLCB < Actual_Yvec & Actual_Yvec < propStudentFit5finalPredsUCB)
)
propStudentFit5predMetrics

##### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = propStudentFit5finalPreds) +
  labs(title = "propStudentFit5 Predict") +
  theme_bw()

propStudentFit5FitDF <- bind_cols(
  StormdataTrain3,
  LCB = propStudentFit5finalFitLCB,
  Mean = propStudentFit5finalFitMean,
  Med = propStudentFit5finalFitMed,
  UCB = propStudentFit5finalFitUCB
)

propStudentFit5stormsFitplot <- ggplot(data = propStudentFit5FitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))
propStudentFit5stormsFitplot

## Prediction
propStudentFit5PredDF <- bind_cols(
  StormdataTest3,
  LCB = propStudentFit5finalPredsLCB,
  Mean = propStudentFit5finalPredsMean,
  Med = propStudentFit5finalPredsMed,
  UCB = propStudentFit5finalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

propStudentFit5stormsPredplot <- ggplot(data = propStudentFit5PredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX, color = "Observed")) +
  geom_line(aes(y = Mean, color = "PPD Mean")) +
  facet_wrap(vars(StormID), ncol = 6)+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50)) +
  labs(title = "propStudentFit5 PPD Mean vs Observed VMAX",
       subtitle = "95% Credible Interval about PPD Mean") +
  scale_color_manual(name = NULL, values = c("black","red")) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1))
  ) +
  theme_bw()
propStudentFit5stormsPredplot

##### PPC ----
###### Quantile 2.5 
propStudentFit5LCBsims <- apply(propStudentFit5finalFit, 
                            MARGIN = 1,
                            function(x){
                              quantile(x, 0.025)
                            })
propStudentFit5LCBpvalueVec <- propStudentFit5LCBsims < quantile(StormdataTrain3$VMAX, 0.025)
propStudentFit5LCBpvalue <- sum(propStudentFit5LCBpvalueVec)
propStudentFit5LCBpvalue <- round(propStudentFit5LCBpvalue/4000, 3)
propStudentFit5LCBpvalue <- min(propStudentFit5LCBpvalue, 1 - propStudentFit5LCBpvalue)

propStudentFit5_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           propStudentFit5finalFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", propStudentFit5LCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#propStudentFit5_ppcLCB

###### Quantile 97.5 
propStudentFit5UCBsims <- apply(propStudentFit5finalFit, 
                            MARGIN = 1,
                            function(x){
                              quantile(x, 0.975)
                            })
propStudentFit5UCBpvalueVec <- propStudentFit5UCBsims < quantile(StormdataTrain3$VMAX, 0.975)
propStudentFit5UCBpvalue <- as.numeric(sum(propStudentFit5UCBpvalueVec))
propStudentFit5UCBpvalue <- round(propStudentFit5UCBpvalue/4000, 3)
propStudentFit5UCBpvalue <- min(propStudentFit5UCBpvalue, 1 - propStudentFit5UCBpvalue)

propStudentFit5_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           propStudentFit5finalFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", propStudentFit5UCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#propStudentFit5_ppcUCB

###### Mean 
propStudentFit5MEANsims <- apply(propStudentFit5finalFit, 
                             MARGIN = 1,
                             function(x){
                               mean(x)
                             })
propStudentFit5MEANpvalueVec <- propStudentFit5MEANsims < mean(StormdataTrain3$VMAX)
propStudentFit5MEANpvalue <- sum(propStudentFit5MEANpvalueVec)
propStudentFit5MEANpvalue <- round(propStudentFit5MEANpvalue/4000, 3)
propStudentFit5MEANpvalue <- min(propStudentFit5MEANpvalue, 1 - propStudentFit5MEANpvalue)

propStudentFit5_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           propStudentFit5finalFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", propStudentFit5MEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#propStudentFit5_ppcMEAN

###### Med 
propStudentFit5MEDsims <- apply(propStudentFit5finalFit, 
                            MARGIN = 1,
                            function(x){
                              quantile(x, 0.5)
                            })
propStudentFit5MEDpvalueVec <- propStudentFit5MEDsims < quantile(StormdataTrain3$VMAX, 0.5)
propStudentFit5MEDpvalue <- sum(propStudentFit5MEDpvalueVec)
propStudentFit5MEDpvalue <- round(propStudentFit5MEDpvalue/4000, 3)
propStudentFit5MEDpvalue <- min(propStudentFit5MEDpvalue, 1 - propStudentFit5MEDpvalue)

propStudentFit5_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           propStudentFit5finalFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", propStudentFit5MEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#propStudentFit5_ppcMED

###### SD 
propStudentFit5SDsims <- apply(propStudentFit5finalFit, 
                           MARGIN = 1,
                           function(x){
                             sd(x)
                           })
propStudentFit5SDpvalueVec <- propStudentFit5SDsims < sd(StormdataTrain3$VMAX)
propStudentFit5SDpvalue <- sum(propStudentFit5SDpvalueVec)
propStudentFit5SDpvalue <- round(propStudentFit5SDpvalue/4000, 3)
propStudentFit5SDpvalue <- min(propStudentFit5SDpvalue, 1 - propStudentFit5SDpvalue)

propStudentFit5_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           propStudentFit5finalFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", propStudentFit5SDpvalue, ")")) +
  theme_bw() +
  legend_none()
#propStudentFit5_ppcSD

###### Range 
propStudentFit5RANGEsims <- apply(propStudentFit5finalFit, 
                              MARGIN = 1,
                              function(x){
                                max(x)-min(x)
                              })
propStudentFit5RANGEpvalueVec <- propStudentFit5RANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
propStudentFit5RANGEpvalue <- sum(propStudentFit5RANGEpvalueVec)
propStudentFit5RANGEpvalue <- round(propStudentFit5RANGEpvalue/4000, 3)
propStudentFit5RANGEpvalue <- min(propStudentFit5RANGEpvalue, 1 - propStudentFit5RANGEpvalue)

propStudentFit5_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           propStudentFit5finalFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", propStudentFit5RANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#propStudentFit5_ppcRANGE

##### Combined Plot ----
propStudentFit5_ppcComb <- 
  propStudentFit5ppcFit /
  (propStudentFit5_ppcLCB | propStudentFit5_ppcMED | propStudentFit5_ppcUCB) /
  (propStudentFit5_ppcRANGE | propStudentFit5_ppcMEAN | propStudentFit5_ppcSD)
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
propStudentFit5_ppcComb

##### Bayes p-values ----
propStudentFit5pvalues <- tibble(
  Fit = "propStudentFit5",
  LCB = propStudentFit5LCBpvalue,
  Median = propStudentFit5MEDpvalue,
  UCB = propStudentFit5UCBpvalue,
  Range = propStudentFit5RANGEpvalue,
  Mean = propStudentFit5MEANpvalue,
  SD = propStudentFit5SDpvalue
)
propStudentFit5pvalues

propStudentFit5files <- ls()[str_detect(ls(), pattern = "propStudentFit5")]
propStudentFit5filesRM <- propStudentFit5files[!(propStudentFit5files %in% c("propStudentFit5_ppcComb",
                                                                 "propStudentFit5loo",
                                                                 "propStudentFit5predMetrics",
                                                                 "propStudentFit5pvalues",
                                                                 "propStudentFit5stormsPredplot"))]

rm(list = propStudentFit5filesRM)
rm(propStudentFit5filesRM)
rm(propStudentFit5files)

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
#save(propGammaFit1, file = "_data/propGammaFit1.RData")
prior_summary(propGammaFit1)
round(posterior_summary(propGammaFit1, probs = c(0.025, 0.975)))
propGammaFit1

print(propGammaFit1, digits = 4)
plot(propGammaFit1)
propGammaFit1ppcFit <- pp_check(propGammaFit1, ndraws = 100) + 
  labs(title = "propGammaFit1 Fit PPC") +
  theme_bw()
propGammaFit1ppcFit
propGammaFit1loo <- loo(propGammaFit1)
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

bayes_factor(propGammaFit1, propGammaFit1)
bayes_factor(propGammaFit1, gammaFit2)
bayes_factor(propGammaFit1, gammaFit3)
bayes_factor(propGammaFit1, studentFit1)
bayes_factor(propGammaFit1, studentFit2)
bayes_factor(propGammaFit1, studentFit3)
bayes_factor(propGammaFit1, linFit11)
bayes_factor(propGammaFit1, propFit1)
bayes_factor(propGammaFit1, logPropFit1)
loo(propGammaFit1, gammaFit3)

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
propGammaFit1finalFit <- t(t(propGammaFit1finalFit)*StormdataTrain3$HWRF)
propGammaFit1finalFitMean <- colMeans(propGammaFit1finalFit)
propGammaFit1finalFitMed <- apply(propGammaFit1finalFit, 2, function(x){quantile(x, 0.5)})
propGammaFit1finalFitLCB <- apply(propGammaFit1finalFit, 2, function(x){quantile(x, 0.025)})
propGammaFit1finalFitUCB <- apply(propGammaFit1finalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
propGammaFit1finalPreds <- posterior_predict(propGammaFit1, 
                                               newdata = StormdataTest8,
                                               allow_new_levels = TRUE)
propGammaFit1finalPreds <- t(t(propGammaFit1finalPreds)*StormdataTest3$HWRF)
propGammaFit1finalPredsMean <- colMeans(propGammaFit1finalPreds)
propGammaFit1finalPredsMed <- apply(propGammaFit1finalPreds, 2, function(x){quantile(x, 0.5)})
propGammaFit1finalPredsLCB <- apply(propGammaFit1finalPreds, 2, function(x){quantile(x, 0.025)})
propGammaFit1finalPredsUCB <- apply(propGammaFit1finalPreds, 2, function(x){quantile(x, 0.975)})

propGammaFit1predMetrics <- tibble(
  Fit = "propGammaFit1",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(propGammaFit1finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(propGammaFit1finalFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < propGammaFit1finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(propGammaFit1finalPredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(propGammaFit1finalPredsMed - Actual_Yvec)),
  COV_pred = mean(propGammaFit1finalPredsLCB < Actual_Yvec & Actual_Yvec < propGammaFit1finalPredsUCB)
)
propGammaFit1predMetrics

##### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = propGammaFit1finalPreds) +
  labs(title = "propGammaFit1 Predict") +
  theme_bw()

propGammaFit1FitDF <- bind_cols(
  StormdataTrain3,
  LCB = propGammaFit1finalFitLCB,
  Mean = propGammaFit1finalFitMean,
  Med = propGammaFit1finalFitMed,
  UCB = propGammaFit1finalFitUCB
)

propGammaFit1stormsFitplot <- ggplot(data = propGammaFit1FitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))
propGammaFit1stormsFitplot

## Prediction
propGammaFit1PredDF <- bind_cols(
  StormdataTest3,
  LCB = propGammaFit1finalPredsLCB,
  Mean = propGammaFit1finalPredsMean,
  Med = propGammaFit1finalPredsMed,
  UCB = propGammaFit1finalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

propGammaFit1stormsPredplot <- ggplot(data = propGammaFit1PredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX, color = "Observed")) +
  geom_line(aes(y = Mean, color = "PPD Mean")) +
  facet_wrap(vars(StormID), ncol = 6)+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50)) +
  labs(title = "propGammaFit1 PPD Mean vs Observed VMAX",
       subtitle = "95% Credible Interval about PPD Mean") +
  scale_color_manual(name = NULL, values = c("black","red")) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1))
  ) +
  theme_bw()
propGammaFit1stormsPredplot

##### PPC ----
###### Quantile 2.5 
propGammaFit1LCBsims <- apply(propGammaFit1finalFit, 
                                MARGIN = 1,
                                function(x){
                                  quantile(x, 0.025)
                                })
propGammaFit1LCBpvalueVec <- propGammaFit1LCBsims < quantile(StormdataTrain3$VMAX, 0.025)
propGammaFit1LCBpvalue <- sum(propGammaFit1LCBpvalueVec)
propGammaFit1LCBpvalue <- round(propGammaFit1LCBpvalue/4000, 3)
propGammaFit1LCBpvalue <- min(propGammaFit1LCBpvalue, 1 - propGammaFit1LCBpvalue)

propGammaFit1_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           propGammaFit1finalFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", propGammaFit1LCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#propGammaFit1_ppcLCB

###### Quantile 97.5 
propGammaFit1UCBsims <- apply(propGammaFit1finalFit, 
                                MARGIN = 1,
                                function(x){
                                  quantile(x, 0.975)
                                })
propGammaFit1UCBpvalueVec <- propGammaFit1UCBsims < quantile(StormdataTrain3$VMAX, 0.975)
propGammaFit1UCBpvalue <- as.numeric(sum(propGammaFit1UCBpvalueVec))
propGammaFit1UCBpvalue <- round(propGammaFit1UCBpvalue/4000, 3)
propGammaFit1UCBpvalue <- min(propGammaFit1UCBpvalue, 1 - propGammaFit1UCBpvalue)

propGammaFit1_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           propGammaFit1finalFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", propGammaFit1UCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#propGammaFit1_ppcUCB

###### Mean 
propGammaFit1MEANsims <- apply(propGammaFit1finalFit, 
                                 MARGIN = 1,
                                 function(x){
                                   mean(x)
                                 })
propGammaFit1MEANpvalueVec <- propGammaFit1MEANsims < mean(StormdataTrain3$VMAX)
propGammaFit1MEANpvalue <- sum(propGammaFit1MEANpvalueVec)
propGammaFit1MEANpvalue <- round(propGammaFit1MEANpvalue/4000, 3)
propGammaFit1MEANpvalue <- min(propGammaFit1MEANpvalue, 1 - propGammaFit1MEANpvalue)

propGammaFit1_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           propGammaFit1finalFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", propGammaFit1MEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#propGammaFit1_ppcMEAN

###### Med 
propGammaFit1MEDsims <- apply(propGammaFit1finalFit, 
                                MARGIN = 1,
                                function(x){
                                  quantile(x, 0.5)
                                })
propGammaFit1MEDpvalueVec <- propGammaFit1MEDsims < quantile(StormdataTrain3$VMAX, 0.5)
propGammaFit1MEDpvalue <- sum(propGammaFit1MEDpvalueVec)
propGammaFit1MEDpvalue <- round(propGammaFit1MEDpvalue/4000, 3)
propGammaFit1MEDpvalue <- min(propGammaFit1MEDpvalue, 1 - propGammaFit1MEDpvalue)

propGammaFit1_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           propGammaFit1finalFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", propGammaFit1MEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#propGammaFit1_ppcMED

###### SD 
propGammaFit1SDsims <- apply(propGammaFit1finalFit, 
                               MARGIN = 1,
                               function(x){
                                 sd(x)
                               })
propGammaFit1SDpvalueVec <- propGammaFit1SDsims < sd(StormdataTrain3$VMAX)
propGammaFit1SDpvalue <- sum(propGammaFit1SDpvalueVec)
propGammaFit1SDpvalue <- round(propGammaFit1SDpvalue/4000, 3)
propGammaFit1SDpvalue <- min(propGammaFit1SDpvalue, 1 - propGammaFit1SDpvalue)

propGammaFit1_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           propGammaFit1finalFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", propGammaFit1SDpvalue, ")")) +
  theme_bw() +
  legend_none()
#propGammaFit1_ppcSD

###### Range 
propGammaFit1RANGEsims <- apply(propGammaFit1finalFit, 
                                  MARGIN = 1,
                                  function(x){
                                    max(x)-min(x)
                                  })
propGammaFit1RANGEpvalueVec <- propGammaFit1RANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
propGammaFit1RANGEpvalue <- sum(propGammaFit1RANGEpvalueVec)
propGammaFit1RANGEpvalue <- round(propGammaFit1RANGEpvalue/4000, 3)
propGammaFit1RANGEpvalue <- min(propGammaFit1RANGEpvalue, 1 - propGammaFit1RANGEpvalue)

propGammaFit1_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           propGammaFit1finalFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", propGammaFit1RANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#propGammaFit1_ppcRANGE

##### Combined Plot ----
propGammaFit1_ppcComb <- 
  propGammaFit1ppcFit /
  (propGammaFit1_ppcLCB | propGammaFit1_ppcMED | propGammaFit1_ppcUCB) /
  (propGammaFit1_ppcRANGE | propGammaFit1_ppcMEAN | propGammaFit1_ppcSD)
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
propGammaFit1_ppcComb

##### Bayes p-values ----
propGammaFit1pvalues <- tibble(
  Fit = "propGammaFit1",
  LCB = propGammaFit1LCBpvalue,
  Median = propGammaFit1MEDpvalue,
  UCB = propGammaFit1UCBpvalue,
  Range = propGammaFit1RANGEpvalue,
  Mean = propGammaFit1MEANpvalue,
  SD = propGammaFit1SDpvalue
)
propGammaFit1pvalues

propGammaFit1files <- ls()[str_detect(ls(), pattern = "propGammaFit1")]
propGammaFit1filesRM <- propGammaFit1files[!(propGammaFit1files %in% c("propGammaFit1_ppcComb",
                                                                             "propGammaFit1loo",
                                                                             "propGammaFit1predMetrics",
                                                                             "propGammaFit1pvalues",
                                                                             "propGammaFit1stormsPredplot"))]

rm(list = propGammaFit1filesRM)
rm(propGammaFit1filesRM)
rm(propGammaFit1files)

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
#save(logpropLinFit1, file = "_data/logpropLinFit1.RData")
prior_summary(logpropLinFit1)
round(posterior_summary(logpropLinFit1, probs = c(0.025, 0.975)))
logpropLinFit1

print(logpropLinFit1, digits = 4)
plot(logpropLinFit1)
logpropLinFit1ppcFit <- pp_check(logpropLinFit1, ndraws = 100) + 
  labs(title = "logpropLinFit1 Fit PPC") +
  theme_bw()
logpropLinFit1ppcFit
logpropLinFit1loo <- loo(logpropLinFit1)
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

bayes_factor(logpropLinFit1, logpropLinFit1)
bayes_factor(logpropLinFit1, gammaFit2)
bayes_factor(logpropLinFit1, gammaFit3)
bayes_factor(logpropLinFit1, studentFit1)
bayes_factor(logpropLinFit1, studentFit2)
bayes_factor(logpropLinFit1, studentFit3)
bayes_factor(logpropLinFit1, linFit11)
bayes_factor(logpropLinFit1, propFit1)
bayes_factor(logpropLinFit1, logPropFit1)
loo(logpropLinFit1, gammaFit3)

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
logpropLinFit1finalFitMean <- colMeans(logpropLinFit1finalFit)
logpropLinFit1finalFitMed <- apply(logpropLinFit1finalFit, 2, function(x){quantile(x, 0.5)})
logpropLinFit1finalFitLCB <- apply(logpropLinFit1finalFit, 2, function(x){quantile(x, 0.025)})
logpropLinFit1finalFitUCB <- apply(logpropLinFit1finalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
logpropLinFit1finalPreds <- posterior_predict(logpropLinFit1, 
                                               newdata = StormdataTest8,
                                               allow_new_levels = TRUE)
logpropLinFit1finalPreds <- t(t(exp(logpropLinFit1finalPreds))*StormdataTest3$HWRF)
logpropLinFit1finalPredsMean <- colMeans(logpropLinFit1finalPreds)
logpropLinFit1finalPredsMed <- apply(logpropLinFit1finalPreds, 2, function(x){quantile(x, 0.5)})
logpropLinFit1finalPredsLCB <- apply(logpropLinFit1finalPreds, 2, function(x){quantile(x, 0.025)})
logpropLinFit1finalPredsUCB <- apply(logpropLinFit1finalPreds, 2, function(x){quantile(x, 0.975)})

logpropLinFit1predMetrics <- tibble(
  Fit = "logpropLinFit1",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(logpropLinFit1finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(logpropLinFit1finalFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < logpropLinFit1finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(logpropLinFit1finalPredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(logpropLinFit1finalPredsMed - Actual_Yvec)),
  COV_pred = mean(logpropLinFit1finalPredsLCB < Actual_Yvec & Actual_Yvec < logpropLinFit1finalPredsUCB)
)
logpropLinFit1predMetrics

##### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = logpropLinFit1finalPreds) +
  labs(title = "logpropLinFit1 Predict") +
  theme_bw()

logpropLinFit1FitDF <- bind_cols(
  StormdataTrain3,
  LCB = logpropLinFit1finalFitLCB,
  Mean = logpropLinFit1finalFitMean,
  Med = logpropLinFit1finalFitMed,
  UCB = logpropLinFit1finalFitUCB
)

logpropLinFit1stormsFitplot <- ggplot(data = logpropLinFit1FitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))
logpropLinFit1stormsFitplot

## Prediction
logpropLinFit1PredDF <- bind_cols(
  StormdataTest3,
  LCB = logpropLinFit1finalPredsLCB,
  Mean = logpropLinFit1finalPredsMean,
  Med = logpropLinFit1finalPredsMed,
  UCB = logpropLinFit1finalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

logpropLinFit1stormsPredplot <- ggplot(data = logpropLinFit1PredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX, color = "Observed")) +
  geom_line(aes(y = Mean, color = "PPD Mean")) +
  facet_wrap(vars(StormID), ncol = 6)+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50)) +
  labs(title = "logpropLinFit1 PPD Mean vs Observed VMAX",
       subtitle = "95% Credible Interval about PPD Mean") +
  scale_color_manual(name = NULL, values = c("black","red")) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1))
  ) +
  theme_bw()
logpropLinFit1stormsPredplot

##### PPC ----
###### Quantile 2.5 
logpropLinFit1LCBsims <- apply(logpropLinFit1finalFit, 
                                MARGIN = 1,
                                function(x){
                                  quantile(x, 0.025)
                                })
logpropLinFit1LCBpvalueVec <- logpropLinFit1LCBsims < quantile(StormdataTrain3$VMAX, 0.025)
logpropLinFit1LCBpvalue <- sum(logpropLinFit1LCBpvalueVec)
logpropLinFit1LCBpvalue <- round(logpropLinFit1LCBpvalue/4000, 3)
logpropLinFit1LCBpvalue <- min(logpropLinFit1LCBpvalue, 1 - logpropLinFit1LCBpvalue)

logpropLinFit1_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logpropLinFit1finalFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", logpropLinFit1LCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logpropLinFit1_ppcLCB

###### Quantile 97.5 
logpropLinFit1UCBsims <- apply(logpropLinFit1finalFit, 
                                MARGIN = 1,
                                function(x){
                                  quantile(x, 0.975)
                                })
logpropLinFit1UCBpvalueVec <- logpropLinFit1UCBsims < quantile(StormdataTrain3$VMAX, 0.975)
logpropLinFit1UCBpvalue <- as.numeric(sum(logpropLinFit1UCBpvalueVec))
logpropLinFit1UCBpvalue <- round(logpropLinFit1UCBpvalue/4000, 3)
logpropLinFit1UCBpvalue <- min(logpropLinFit1UCBpvalue, 1 - logpropLinFit1UCBpvalue)

logpropLinFit1_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logpropLinFit1finalFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", logpropLinFit1UCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logpropLinFit1_ppcUCB

###### Mean 
logpropLinFit1MEANsims <- apply(logpropLinFit1finalFit, 
                                 MARGIN = 1,
                                 function(x){
                                   mean(x)
                                 })
logpropLinFit1MEANpvalueVec <- logpropLinFit1MEANsims < mean(StormdataTrain3$VMAX)
logpropLinFit1MEANpvalue <- sum(logpropLinFit1MEANpvalueVec)
logpropLinFit1MEANpvalue <- round(logpropLinFit1MEANpvalue/4000, 3)
logpropLinFit1MEANpvalue <- min(logpropLinFit1MEANpvalue, 1 - logpropLinFit1MEANpvalue)

logpropLinFit1_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           logpropLinFit1finalFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", logpropLinFit1MEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#logpropLinFit1_ppcMEAN

###### Med 
logpropLinFit1MEDsims <- apply(logpropLinFit1finalFit, 
                                MARGIN = 1,
                                function(x){
                                  quantile(x, 0.5)
                                })
logpropLinFit1MEDpvalueVec <- logpropLinFit1MEDsims < quantile(StormdataTrain3$VMAX, 0.5)
logpropLinFit1MEDpvalue <- sum(logpropLinFit1MEDpvalueVec)
logpropLinFit1MEDpvalue <- round(logpropLinFit1MEDpvalue/4000, 3)
logpropLinFit1MEDpvalue <- min(logpropLinFit1MEDpvalue, 1 - logpropLinFit1MEDpvalue)

logpropLinFit1_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           logpropLinFit1finalFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", logpropLinFit1MEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logpropLinFit1_ppcMED

###### SD 
logpropLinFit1SDsims <- apply(logpropLinFit1finalFit, 
                               MARGIN = 1,
                               function(x){
                                 sd(x)
                               })
logpropLinFit1SDpvalueVec <- logpropLinFit1SDsims < sd(StormdataTrain3$VMAX)
logpropLinFit1SDpvalue <- sum(logpropLinFit1SDpvalueVec)
logpropLinFit1SDpvalue <- round(logpropLinFit1SDpvalue/4000, 3)
logpropLinFit1SDpvalue <- min(logpropLinFit1SDpvalue, 1 - logpropLinFit1SDpvalue)

logpropLinFit1_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           logpropLinFit1finalFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", logpropLinFit1SDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logpropLinFit1_ppcSD

###### Range 
logpropLinFit1RANGEsims <- apply(logpropLinFit1finalFit, 
                                  MARGIN = 1,
                                  function(x){
                                    max(x)-min(x)
                                  })
logpropLinFit1RANGEpvalueVec <- logpropLinFit1RANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
logpropLinFit1RANGEpvalue <- sum(logpropLinFit1RANGEpvalueVec)
logpropLinFit1RANGEpvalue <- round(logpropLinFit1RANGEpvalue/4000, 3)
logpropLinFit1RANGEpvalue <- min(logpropLinFit1RANGEpvalue, 1 - logpropLinFit1RANGEpvalue)

logpropLinFit1_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           logpropLinFit1finalFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", logpropLinFit1RANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#logpropLinFit1_ppcRANGE

##### Combined Plot ----
logpropLinFit1_ppcComb <- 
  logpropLinFit1ppcFit /
  (logpropLinFit1_ppcLCB | logpropLinFit1_ppcMED | logpropLinFit1_ppcUCB) /
  (logpropLinFit1_ppcRANGE | logpropLinFit1_ppcMEAN | logpropLinFit1_ppcSD)
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
logpropLinFit1_ppcComb

##### Bayes p-values ----
logpropLinFit1pvalues <- tibble(
  Fit = "logpropLinFit1",
  LCB = logpropLinFit1LCBpvalue,
  Median = logpropLinFit1MEDpvalue,
  UCB = logpropLinFit1UCBpvalue,
  Range = logpropLinFit1RANGEpvalue,
  Mean = logpropLinFit1MEANpvalue,
  SD = logpropLinFit1SDpvalue
)
logpropLinFit1pvalues

logpropLinFit1files <- ls()[str_detect(ls(), pattern = "logpropLinFit1")]
logpropLinFit1filesRM <- logpropLinFit1files[!(logpropLinFit1files %in% c("logpropLinFit1_ppcComb",
                                                                             "logpropLinFit1loo",
                                                                             "logpropLinFit1predMetrics",
                                                                             "logpropLinFit1pvalues",
                                                                             "logpropLinFit1stormsPredplot"))]

rm(list = logpropLinFit1filesRM)
rm(logpropLinFit1filesRM)
rm(logpropLinFit1files)

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
#save(logpropStudentFit1, file = "_data/logpropStudentFit1.RData")
prior_summary(logpropStudentFit1)
round(posterior_summary(logpropStudentFit1, probs = c(0.025, 0.975)))
logpropStudentFit1

print(logpropStudentFit1, digits = 4)
plot(logpropStudentFit1)
logpropStudentFit1ppcFit <- pp_check(logpropStudentFit1, ndraws = 100) + 
  labs(title = "logpropStudentFit1 Fit PPC") +
  theme_bw()
logpropStudentFit1ppcFit
logpropStudentFit1loo <- loo(logpropStudentFit1)
waic(logpropStudentFit1)
performance::check_distribution(logpropStudentFit1)
performance::check_outliers(logpropStudentFit1)
performance::check_heteroskedasticity(logpropStudentFit1)
performance_rmse(logpropStudentFit1)
performance_mae(logpropStudentFit1)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(logpropStudentFit1)
logpropStudentFit1logLik <- log_lik(logpropStudentFit1)

variance_decomposition(logpropStudentFit1)
exp(fixef(logpropStudentFit1))
ranef(logpropStudentFit1)

bayes_R2(logpropStudentFit1)

bayes_factor(logpropStudentFit1, logpropStudentFit1)
bayes_factor(logpropStudentFit1, gammaFit2)
bayes_factor(logpropStudentFit1, gammaFit3)
bayes_factor(logpropStudentFit1, studentFit1)
bayes_factor(logpropStudentFit1, studentFit2)
bayes_factor(logpropStudentFit1, studentFit3)
bayes_factor(logpropStudentFit1, linFit11)
bayes_factor(logpropStudentFit1, propFit1)
bayes_factor(logpropStudentFit1, logPropFit1)
loo(logpropStudentFit1, gammaFit3)

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
logpropStudentFit1finalFit <- t(t(exp(logpropStudentFit1finalFit))*StormdataTrain3$HWRF)
logpropStudentFit1finalFitMean <- colMeans(logpropStudentFit1finalFit)
logpropStudentFit1finalFitMed <- apply(logpropStudentFit1finalFit, 2, function(x){quantile(x, 0.5)})
logpropStudentFit1finalFitLCB <- apply(logpropStudentFit1finalFit, 2, function(x){quantile(x, 0.025)})
logpropStudentFit1finalFitUCB <- apply(logpropStudentFit1finalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
logpropStudentFit1finalPreds <- posterior_predict(logpropStudentFit1, 
                                              newdata = StormdataTest8,
                                              allow_new_levels = TRUE)
logpropStudentFit1finalPreds <- t(t(exp(logpropStudentFit1finalPreds))*StormdataTest3$HWRF)
logpropStudentFit1finalPredsMean <- colMeans(logpropStudentFit1finalPreds)
logpropStudentFit1finalPredsMed <- apply(logpropStudentFit1finalPreds, 2, function(x){quantile(x, 0.5)})
logpropStudentFit1finalPredsLCB <- apply(logpropStudentFit1finalPreds, 2, function(x){quantile(x, 0.025)})
logpropStudentFit1finalPredsUCB <- apply(logpropStudentFit1finalPreds, 2, function(x){quantile(x, 0.975)})

logpropStudentFit1predMetrics <- tibble(
  Fit = "logpropStudentFit1",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(logpropStudentFit1finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(logpropStudentFit1finalFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < logpropStudentFit1finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(logpropStudentFit1finalPredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(logpropStudentFit1finalPredsMed - Actual_Yvec)),
  COV_pred = mean(logpropStudentFit1finalPredsLCB < Actual_Yvec & Actual_Yvec < logpropStudentFit1finalPredsUCB)
)
logpropStudentFit1predMetrics

##### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = logpropStudentFit1finalPreds) +
  labs(title = "logpropStudentFit1 Predict") +
  theme_bw()

logpropStudentFit1FitDF <- bind_cols(
  StormdataTrain3,
  LCB = logpropStudentFit1finalFitLCB,
  Mean = logpropStudentFit1finalFitMean,
  Med = logpropStudentFit1finalFitMed,
  UCB = logpropStudentFit1finalFitUCB
)

logpropStudentFit1stormsFitplot <- ggplot(data = logpropStudentFit1FitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))
logpropStudentFit1stormsFitplot

## Prediction
logpropStudentFit1PredDF <- bind_cols(
  StormdataTest3,
  LCB = logpropStudentFit1finalPredsLCB,
  Mean = logpropStudentFit1finalPredsMean,
  Med = logpropStudentFit1finalPredsMed,
  UCB = logpropStudentFit1finalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

logpropStudentFit1stormsPredplot <- ggplot(data = logpropStudentFit1PredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX, color = "Observed")) +
  geom_line(aes(y = Mean, color = "PPD Mean")) +
  facet_wrap(vars(StormID), ncol = 6)+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50)) +
  labs(title = "logpropStudentFit1 PPD Mean vs Observed VMAX",
       subtitle = "95% Credible Interval about PPD Mean") +
  scale_color_manual(name = NULL, values = c("black","red")) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1))
  ) +
  theme_bw()
logpropStudentFit1stormsPredplot

##### PPC ----
###### Quantile 2.5 
logpropStudentFit1LCBsims <- apply(logpropStudentFit1finalFit, 
                               MARGIN = 1,
                               function(x){
                                 quantile(x, 0.025)
                               })
logpropStudentFit1LCBpvalueVec <- logpropStudentFit1LCBsims < quantile(StormdataTrain3$VMAX, 0.025)
logpropStudentFit1LCBpvalue <- sum(logpropStudentFit1LCBpvalueVec)
logpropStudentFit1LCBpvalue <- round(logpropStudentFit1LCBpvalue/4000, 3)
logpropStudentFit1LCBpvalue <- min(logpropStudentFit1LCBpvalue, 1 - logpropStudentFit1LCBpvalue)

logpropStudentFit1_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logpropStudentFit1finalFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", logpropStudentFit1LCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logpropStudentFit1_ppcLCB

###### Quantile 97.5 
logpropStudentFit1UCBsims <- apply(logpropStudentFit1finalFit, 
                               MARGIN = 1,
                               function(x){
                                 quantile(x, 0.975)
                               })
logpropStudentFit1UCBpvalueVec <- logpropStudentFit1UCBsims < quantile(StormdataTrain3$VMAX, 0.975)
logpropStudentFit1UCBpvalue <- as.numeric(sum(logpropStudentFit1UCBpvalueVec))
logpropStudentFit1UCBpvalue <- round(logpropStudentFit1UCBpvalue/4000, 3)
logpropStudentFit1UCBpvalue <- min(logpropStudentFit1UCBpvalue, 1 - logpropStudentFit1UCBpvalue)

logpropStudentFit1_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logpropStudentFit1finalFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", logpropStudentFit1UCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logpropStudentFit1_ppcUCB

###### Mean 
logpropStudentFit1MEANsims <- apply(logpropStudentFit1finalFit, 
                                MARGIN = 1,
                                function(x){
                                  mean(x)
                                })
logpropStudentFit1MEANpvalueVec <- logpropStudentFit1MEANsims < mean(StormdataTrain3$VMAX)
logpropStudentFit1MEANpvalue <- sum(logpropStudentFit1MEANpvalueVec)
logpropStudentFit1MEANpvalue <- round(logpropStudentFit1MEANpvalue/4000, 3)
logpropStudentFit1MEANpvalue <- min(logpropStudentFit1MEANpvalue, 1 - logpropStudentFit1MEANpvalue)

logpropStudentFit1_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           logpropStudentFit1finalFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", logpropStudentFit1MEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#logpropStudentFit1_ppcMEAN

###### Med 
logpropStudentFit1MEDsims <- apply(logpropStudentFit1finalFit, 
                               MARGIN = 1,
                               function(x){
                                 quantile(x, 0.5)
                               })
logpropStudentFit1MEDpvalueVec <- logpropStudentFit1MEDsims < quantile(StormdataTrain3$VMAX, 0.5)
logpropStudentFit1MEDpvalue <- sum(logpropStudentFit1MEDpvalueVec)
logpropStudentFit1MEDpvalue <- round(logpropStudentFit1MEDpvalue/4000, 3)
logpropStudentFit1MEDpvalue <- min(logpropStudentFit1MEDpvalue, 1 - logpropStudentFit1MEDpvalue)

logpropStudentFit1_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           logpropStudentFit1finalFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", logpropStudentFit1MEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logpropStudentFit1_ppcMED

###### SD 
logpropStudentFit1SDsims <- apply(logpropStudentFit1finalFit, 
                              MARGIN = 1,
                              function(x){
                                sd(x)
                              })
logpropStudentFit1SDpvalueVec <- logpropStudentFit1SDsims < sd(StormdataTrain3$VMAX)
logpropStudentFit1SDpvalue <- sum(logpropStudentFit1SDpvalueVec)
logpropStudentFit1SDpvalue <- round(logpropStudentFit1SDpvalue/4000, 3)
logpropStudentFit1SDpvalue <- min(logpropStudentFit1SDpvalue, 1 - logpropStudentFit1SDpvalue)

logpropStudentFit1_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           logpropStudentFit1finalFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", logpropStudentFit1SDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logpropStudentFit1_ppcSD

###### Range 
logpropStudentFit1RANGEsims <- apply(logpropStudentFit1finalFit, 
                                 MARGIN = 1,
                                 function(x){
                                   max(x)-min(x)
                                 })
logpropStudentFit1RANGEpvalueVec <- logpropStudentFit1RANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
logpropStudentFit1RANGEpvalue <- sum(logpropStudentFit1RANGEpvalueVec)
logpropStudentFit1RANGEpvalue <- round(logpropStudentFit1RANGEpvalue/4000, 3)
logpropStudentFit1RANGEpvalue <- min(logpropStudentFit1RANGEpvalue, 1 - logpropStudentFit1RANGEpvalue)

logpropStudentFit1_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           logpropStudentFit1finalFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", logpropStudentFit1RANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#logpropStudentFit1_ppcRANGE

##### Combined Plot ----
logpropStudentFit1_ppcComb <- 
  logpropStudentFit1ppcFit /
  (logpropStudentFit1_ppcLCB | logpropStudentFit1_ppcMED | logpropStudentFit1_ppcUCB) /
  (logpropStudentFit1_ppcRANGE | logpropStudentFit1_ppcMEAN | logpropStudentFit1_ppcSD)
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
logpropStudentFit1_ppcComb

##### Bayes p-values ----
logpropStudentFit1pvalues <- tibble(
  Fit = "logpropStudentFit1",
  LCB = logpropStudentFit1LCBpvalue,
  Median = logpropStudentFit1MEDpvalue,
  UCB = logpropStudentFit1UCBpvalue,
  Range = logpropStudentFit1RANGEpvalue,
  Mean = logpropStudentFit1MEANpvalue,
  SD = logpropStudentFit1SDpvalue
)
logpropStudentFit1pvalues

logpropStudentFit1files <- ls()[str_detect(ls(), pattern = "logpropStudentFit1")]
logpropStudentFit1filesRM <- logpropStudentFit1files[!(logpropStudentFit1files %in% c("logpropStudentFit1_ppcComb",
                                                                          "logpropStudentFit1loo",
                                                                          "logpropStudentFit1predMetrics",
                                                                          "logpropStudentFit1pvalues",
                                                                          "logpropStudentFit1stormsPredplot"))]

rm(list = logpropStudentFit1filesRM)
rm(logpropStudentFit1filesRM)
rm(logpropStudentFit1files)

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
#save(logpropStudentFit2, file = "_data/logpropStudentFit2.RData")
prior_summary(logpropStudentFit2)
round(posterior_summary(logpropStudentFit2, probs = c(0.025, 0.975)))
logpropStudentFit2

print(logpropStudentFit2, digits = 4)
plot(logpropStudentFit2)
logpropStudentFit2ppcFit <- pp_check(logpropStudentFit2, ndraws = 100) + 
  labs(title = "logpropStudentFit2 Fit PPC") +
  theme_bw()
logpropStudentFit2ppcFit
logpropStudentFit2loo <- loo(logpropStudentFit2)
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

bayes_factor(logpropStudentFit2, logpropStudentFit2)
bayes_factor(logpropStudentFit2, gammaFit2)
bayes_factor(logpropStudentFit2, gammaFit3)
bayes_factor(logpropStudentFit2, studentFit1)
bayes_factor(logpropStudentFit2, studentFit2)
bayes_factor(logpropStudentFit2, studentFit3)
bayes_factor(logpropStudentFit2, linFit11)
bayes_factor(logpropStudentFit2, propFit1)
bayes_factor(logpropStudentFit2, logPropFit1)
loo(logpropStudentFit2, gammaFit3)

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
logpropStudentFit2finalFit <- t(t(exp(logpropStudentFit2finalFit))*StormdataTrain3$HWRF)
logpropStudentFit2finalFitMean <- colMeans(logpropStudentFit2finalFit)
logpropStudentFit2finalFitMed <- apply(logpropStudentFit2finalFit, 2, function(x){quantile(x, 0.5)})
logpropStudentFit2finalFitLCB <- apply(logpropStudentFit2finalFit, 2, function(x){quantile(x, 0.025)})
logpropStudentFit2finalFitUCB <- apply(logpropStudentFit2finalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
logpropStudentFit2finalPreds <- posterior_predict(logpropStudentFit2, 
                                              newdata = StormdataTest8,
                                              allow_new_levels = TRUE)
logpropStudentFit2finalPreds <- t(t(exp(logpropStudentFit2finalPreds))*StormdataTest3$HWRF)
logpropStudentFit2finalPredsMean <- colMeans(logpropStudentFit2finalPreds)
logpropStudentFit2finalPredsMed <- apply(logpropStudentFit2finalPreds, 2, function(x){quantile(x, 0.5)})
logpropStudentFit2finalPredsLCB <- apply(logpropStudentFit2finalPreds, 2, function(x){quantile(x, 0.025)})
logpropStudentFit2finalPredsUCB <- apply(logpropStudentFit2finalPreds, 2, function(x){quantile(x, 0.975)})

logpropStudentFit2predMetrics <- tibble(
  Fit = "logpropStudentFit2",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(logpropStudentFit2finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(logpropStudentFit2finalFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < logpropStudentFit2finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(logpropStudentFit2finalPredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(logpropStudentFit2finalPredsMed - Actual_Yvec)),
  COV_pred = mean(logpropStudentFit2finalPredsLCB < Actual_Yvec & Actual_Yvec < logpropStudentFit2finalPredsUCB)
)
logpropStudentFit2predMetrics

##### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = logpropStudentFit2finalPreds) +
  labs(title = "logpropStudentFit2 Predict") +
  theme_bw()

logpropStudentFit2FitDF <- bind_cols(
  StormdataTrain3,
  LCB = logpropStudentFit2finalFitLCB,
  Mean = logpropStudentFit2finalFitMean,
  Med = logpropStudentFit2finalFitMed,
  UCB = logpropStudentFit2finalFitUCB
)

logpropStudentFit2stormsFitplot <- ggplot(data = logpropStudentFit2FitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))
logpropStudentFit2stormsFitplot

## Prediction
logpropStudentFit2PredDF <- bind_cols(
  StormdataTest3,
  LCB = logpropStudentFit2finalPredsLCB,
  Mean = logpropStudentFit2finalPredsMean,
  Med = logpropStudentFit2finalPredsMed,
  UCB = logpropStudentFit2finalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

logpropStudentFit2stormsPredplot <- ggplot(data = logpropStudentFit2PredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX, color = "Observed")) +
  geom_line(aes(y = Mean, color = "PPD Mean")) +
  facet_wrap(vars(StormID), ncol = 6)+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50)) +
  labs(title = "logpropStudentFit2 PPD Mean vs Observed VMAX",
       subtitle = "95% Credible Interval about PPD Mean") +
  scale_color_manual(name = NULL, values = c("black","red")) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1))
  ) +
  theme_bw()
logpropStudentFit2stormsPredplot

##### PPC ----
###### Quantile 2.5 
logpropStudentFit2LCBsims <- apply(logpropStudentFit2finalFit, 
                               MARGIN = 1,
                               function(x){
                                 quantile(x, 0.025)
                               })
logpropStudentFit2LCBpvalueVec <- logpropStudentFit2LCBsims < quantile(StormdataTrain3$VMAX, 0.025)
logpropStudentFit2LCBpvalue <- sum(logpropStudentFit2LCBpvalueVec)
logpropStudentFit2LCBpvalue <- round(logpropStudentFit2LCBpvalue/4000, 3)
logpropStudentFit2LCBpvalue <- min(logpropStudentFit2LCBpvalue, 1 - logpropStudentFit2LCBpvalue)

logpropStudentFit2_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logpropStudentFit2finalFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", logpropStudentFit2LCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logpropStudentFit2_ppcLCB

###### Quantile 97.5 
logpropStudentFit2UCBsims <- apply(logpropStudentFit2finalFit, 
                               MARGIN = 1,
                               function(x){
                                 quantile(x, 0.975)
                               })
logpropStudentFit2UCBpvalueVec <- logpropStudentFit2UCBsims < quantile(StormdataTrain3$VMAX, 0.975)
logpropStudentFit2UCBpvalue <- as.numeric(sum(logpropStudentFit2UCBpvalueVec))
logpropStudentFit2UCBpvalue <- round(logpropStudentFit2UCBpvalue/4000, 3)
logpropStudentFit2UCBpvalue <- min(logpropStudentFit2UCBpvalue, 1 - logpropStudentFit2UCBpvalue)

logpropStudentFit2_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logpropStudentFit2finalFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", logpropStudentFit2UCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logpropStudentFit2_ppcUCB

###### Mean 
logpropStudentFit2MEANsims <- apply(logpropStudentFit2finalFit, 
                                MARGIN = 1,
                                function(x){
                                  mean(x)
                                })
logpropStudentFit2MEANpvalueVec <- logpropStudentFit2MEANsims < mean(StormdataTrain3$VMAX)
logpropStudentFit2MEANpvalue <- sum(logpropStudentFit2MEANpvalueVec)
logpropStudentFit2MEANpvalue <- round(logpropStudentFit2MEANpvalue/4000, 3)
logpropStudentFit2MEANpvalue <- min(logpropStudentFit2MEANpvalue, 1 - logpropStudentFit2MEANpvalue)

logpropStudentFit2_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           logpropStudentFit2finalFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", logpropStudentFit2MEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#logpropStudentFit2_ppcMEAN

###### Med 
logpropStudentFit2MEDsims <- apply(logpropStudentFit2finalFit, 
                               MARGIN = 1,
                               function(x){
                                 quantile(x, 0.5)
                               })
logpropStudentFit2MEDpvalueVec <- logpropStudentFit2MEDsims < quantile(StormdataTrain3$VMAX, 0.5)
logpropStudentFit2MEDpvalue <- sum(logpropStudentFit2MEDpvalueVec)
logpropStudentFit2MEDpvalue <- round(logpropStudentFit2MEDpvalue/4000, 3)
logpropStudentFit2MEDpvalue <- min(logpropStudentFit2MEDpvalue, 1 - logpropStudentFit2MEDpvalue)

logpropStudentFit2_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           logpropStudentFit2finalFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", logpropStudentFit2MEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logpropStudentFit2_ppcMED

###### SD 
logpropStudentFit2SDsims <- apply(logpropStudentFit2finalFit, 
                              MARGIN = 1,
                              function(x){
                                sd(x)
                              })
logpropStudentFit2SDpvalueVec <- logpropStudentFit2SDsims < sd(StormdataTrain3$VMAX)
logpropStudentFit2SDpvalue <- sum(logpropStudentFit2SDpvalueVec)
logpropStudentFit2SDpvalue <- round(logpropStudentFit2SDpvalue/4000, 3)
logpropStudentFit2SDpvalue <- min(logpropStudentFit2SDpvalue, 1 - logpropStudentFit2SDpvalue)

logpropStudentFit2_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           logpropStudentFit2finalFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", logpropStudentFit2SDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logpropStudentFit2_ppcSD

###### Range 
logpropStudentFit2RANGEsims <- apply(logpropStudentFit2finalFit, 
                                 MARGIN = 1,
                                 function(x){
                                   max(x)-min(x)
                                 })
logpropStudentFit2RANGEpvalueVec <- logpropStudentFit2RANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
logpropStudentFit2RANGEpvalue <- sum(logpropStudentFit2RANGEpvalueVec)
logpropStudentFit2RANGEpvalue <- round(logpropStudentFit2RANGEpvalue/4000, 3)
logpropStudentFit2RANGEpvalue <- min(logpropStudentFit2RANGEpvalue, 1 - logpropStudentFit2RANGEpvalue)

logpropStudentFit2_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           logpropStudentFit2finalFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", logpropStudentFit2RANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#logpropStudentFit2_ppcRANGE

##### Combined Plot ----
logpropStudentFit2_ppcComb <- 
  logpropStudentFit2ppcFit /
  (logpropStudentFit2_ppcLCB | logpropStudentFit2_ppcMED | logpropStudentFit2_ppcUCB) /
  (logpropStudentFit2_ppcRANGE | logpropStudentFit2_ppcMEAN | logpropStudentFit2_ppcSD)
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
logpropStudentFit2_ppcComb

##### Bayes p-values ----
logpropStudentFit2pvalues <- tibble(
  Fit = "logpropStudentFit2",
  LCB = logpropStudentFit2LCBpvalue,
  Median = logpropStudentFit2MEDpvalue,
  UCB = logpropStudentFit2UCBpvalue,
  Range = logpropStudentFit2RANGEpvalue,
  Mean = logpropStudentFit2MEANpvalue,
  SD = logpropStudentFit2SDpvalue
)
logpropStudentFit2pvalues

logpropStudentFit2files <- ls()[str_detect(ls(), pattern = "logpropStudentFit2")]
logpropStudentFit2filesRM <- logpropStudentFit2files[!(logpropStudentFit2files %in% c("logpropStudentFit2_ppcComb",
                                                                          "logpropStudentFit2loo",
                                                                          "logpropStudentFit2predMetrics",
                                                                          "logpropStudentFit2pvalues",
                                                                          "logpropStudentFit2stormsPredplot"))]

rm(list = logpropStudentFit2filesRM)
rm(logpropStudentFit2filesRM)
rm(logpropStudentFit2files)

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
                                                  newdata = StormdataTest8,
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
                                                  newdata = StormdataTest8,
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
                                                  newdata = StormdataTest8,
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

###### Quantile 2.5 ----
logpropStudentFit5_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logpropStudentFit5finalFit2,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = "2.5% Quantile") +
  theme_bw() +
  legend_none()
logpropStudentFit5_ppcLCB

logpropStudentFit5_ppcLCBdens <- density(logpropStudentFit5finalFitLCB)
logpropStudentFit5_ppcLCBdens <- cbind(logpropStudentFit5_ppcLCBdens$x, 
                                       logpropStudentFit5_ppcLCBdens$y)
logpropStudentFit5_ppcLCBdensB <- logpropStudentFit5_ppcLCBdens[between(logpropStudentFit5_ppcLCBdens[,1], 
                                                                        quantile(logpropStudentFit5finalFitLCB, 0.025), 
                                                                        quantile(logpropStudentFit5finalFitLCB, 0.975)), ] 
logpropStudentFit5_ppcLCB2 <- ggplot() +
  # geom_histogram(aes(x = ppcL2,  after_stat(density)),
  #                fill = "#bcdcdc", color = "#99c7c7") +
  geom_ribbon(aes(x = logpropStudentFit5_ppcLCBdensB[,1],
                  ymin = 0, ymax = logpropStudentFit5_ppcLCBdensB[,2]),
              fill = "#bcdcdc") +
  geom_density(aes(x = logpropStudentFit5finalFitLCB), 
               color = "#99c7c7", linewidth = .75) +
  #geom_vline(aes(xintercept = quantile(ppcL2, 0.975)), color = "#007C7C", linewidth = 2) +
  geom_vline(aes(xintercept = DppcL2), color = "#007C7C", linewidth = 1) +
  geom_text(aes(label = paste0("p-value = ", 
                               round(mean(logpropStudentFit5finalFitLCB > DppcL2), 3))),
            x = 0.84*max(ppcL2), y = max(logpropStudentFit5_ppcLCBdens[,2]), size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  labs(title = "2.5% Quantile",
       x = NULL,
       y = "Posterior Density") +
  theme_bw() +
  theme(
    plot.title = element_text(size = rel(1)))
ppc_q2.5_plot2B

###### Quantile 97.5 ----
logpropStudentFit5_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logpropStudentFit5finalFit2,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = "97.5% Quantile") +
  theme_bw() +
  legend_none()
logpropStudentFit5_ppcUCB

###### Mean ----
logpropStudentFit5_ppcMean <- 
  ppc_stat(StormdataTrain3$VMAX,
           logpropStudentFit5finalFit2,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = "Mean") +
  theme_bw() +
  legend_none()
logpropStudentFit5_ppcMean

###### Med ----
logpropStudentFit5_ppcMed <- 
  ppc_stat(StormdataTrain3$VMAX,
           logpropStudentFit5finalFit2,
           stat = function(y) median(y), freq = FALSE) +
  labs(title = "Median") +
  theme_bw() +
  legend_none()
logpropStudentFit5_ppcMed

###### SD ----
logpropStudentFit5_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           logpropStudentFit5finalFit2,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = "SD") +
  theme_bw() +
  legend_none()
logpropStudentFit5_ppcSD

rm(logpropStudentFit5finalFit)
rm(logpropStudentFit5finalFit2)
rm(logpropStudentFit5finalPreds)

# Compare Predictions ----
## LOO ----
looComp <- loo_compare(gammaFit1loo,
                       logNormalFit1loo,
                       logNormalFit2Aloo,
                       #logNormalFit2Bloo,
                       propLinFit1loo,
                       propStudentFit1loo,
                       propStudentFit2loo,
                       #propStudentFit3loo,
                       #propStudentFit4loo,
                       propStudentFit5loo,
                       propGammaFit1loo,
                       logpropLinFit1loo,
                       logpropStudentFit1loo,
                       logpropStudentFit2loo
                       #logpropStudentFit3loo,
                       #logpropStudentFit4loo,
                       #logpropStudentFit5loo
)
looComp
save(looComp, file = "_data/looCompFinal.RData")

## P-Values ----
pvaluesDF <- bind_rows(
  gammaFit1pvalues,
  logNormalFit1pvalues,
  logNormalFit2Apvalues,
  propLinFit1pvalues,
  propStudentFit1pvalues,
  propStudentFit2pvalues,
  propStudentFit5pvalues,
  propGammaFit1pvalues,
  logpropLinFit1pvalues,
  logpropStudentFit1pvalues,
  logpropStudentFit2pvalues
)

## Pred Metrics ----
predCompMetrics <- bind_rows(
  gammaFit1predMetrics, # |> bind_cols(Fit = "gammaFit1"),
  logNormalFit1predMetrics, # |> bind_cols(Fit = "logNormalFit1"),
  logNormalFit2ApredMetrics, # |> bind_cols(Fit = "logNormalFit2A"),
  propLinFit1predMetrics, # |> bind_cols(Fit = "propLinFit1"),
  propStudentFit1predMetrics, # |> bind_cols(Fit = "propStudentFit1"),
  propStudentFit2predMetrics, # |> bind_cols(Fit = "propStudentFit2"),
  #propStudentFit3predMetrics |> bind_cols(Fit = "propStudentFit3"),
  #propStudentFit4predMetrics |> bind_cols(Fit = "propStudentFit4"),
  propStudentFit5predMetrics, # |> bind_cols(Fit = "propStudentFit5"),
  propGammaFit1predMetrics, # |> bind_cols(Fit = "propGammaFit1"),
  logpropLinFit1predMetrics, # |> bind_cols(Fit = "logpropLinFit1"),
  logpropStudentFit1predMetrics, # |> bind_cols(Fit = "logpropStudentFit1"),
  logpropStudentFit2predMetrics, # |> bind_cols(Fit = "logpropStudentFit2"),
  #logpropStudentFit3predMetrics |> bind_cols(Fit = "logpropStudentFit3"),
  #logpropStudentFit4predMetrics |> bind_cols(Fit = "logpropStudentFit4"),
  #logpropStudentFit5predMetrics |> bind_cols(Fit = "logpropStudentFit5"),
)
predCompMetrics <- predCompMetrics |> arrange(MAE_pred)

## Ranked Fits ----
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
  arrange(OverallRank) |>
  left_join(pvaluesDF)

save(rankComps, file = "_data/rankCompsFinal.RData")



