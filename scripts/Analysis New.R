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

# str(Stormdata)

# Timediff <- Stormdata |>
#   group_by(StormID) |>
#   summarise(
#     difftime(Date, lag(Date, default = Date[1]), units = "hours")
#   )
# nondense <- which(Timediff$`difftime(Date, lag(Date, default = Date[1]), units = "hours")` > 6)
# nondenseID <- Stormdata |> slice(nondense)


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

summary(StormdataTrain3$VMAX)

#### Train 4 ----
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


#### Train 5 ----
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

#### Train 6 ----
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
plot(x = StormdataTrain7$Day, y = StormdataTrain7$VMAX)
str(StormdataTrain7)

#### Train 7B ----
StormdataTrain7B <- StormdataTrain3 |>
  filter(
    StormID %in% c(602017, 3712016, 5412015, 3012014, 3402016, 702017)
  ) |>
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
str(StormdataTrain7B)

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
ggpairs(StormdataTrain3 |> 
          select()
)

## Scatter ----
# [1] "StormID"           "Date"              "Year"              "Month"             "Day"              
# [6] "StormElapsedTime"  "StormElapsedTime2" "basin"             "LAT"               "LON"              
# [11] "MINSLP"            "SHR_MAG"           "STM_SPD"           "SST"               "RHLO"             
# [16] "CAPE1"             "CAPE3"             "SHTFL2"            "TCOND7002"         "INST2"            
# [21] "CP1"               "TCONDSYM2"         "COUPLSYM3"         "HWFI"              "VMAX_OP_T0"       
# [26] "HWRF"              "VMAX"
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
load(file = "_data/studentFit1.RData")
load(file = "_data/studentFit2.RData")
load(file = "_data/studentFit3.RData")
load(file = "_data/linFit10.RData")
load(file = "_data/linFit11.RData")
load(file = "_data/linFit12.RData")

## GAUSSIAN ----
### Model 1 ----

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


### Model 2 ----

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

### Model 3 ----

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


### Model 4 ----
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

### Model 5 ----
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

### Model 6 ----
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

### Model 7 ----
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
pp_check(linFit7)
loo(linFit7)
performance::check_distribution(linFit7)
performance::check_outliers(linFit7)
performance::check_heteroskedasticity(linFit7)
performance_rmse(linFit7)
performance_mae(linFit7)

variance_decomposition(linFit7)
fixef(linFit7)
ranef(linFit7)

bayes_R2(linFit7)

bayes_factor(linFit7, linFit4)
bayes_factor(linFit7, linFit5)
bayes_factor(linFit7, linFit6)

conditional_smooths(linFit7)
conditional_effects(linFit7)



### Model 8 ----
linFit8 <- brm(
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
  data = StormdataTrain5, 
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
pp_check(linFit8)
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

conditional_smooths(linFit8)
conditional_effects(linFit8)

### Model 9 ----

linFit9 <- brm(
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
  data = StormdataTrain7, 
  family = gaussian(), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52, 
  warmup = 1000
)

prior_summary(linFit9)
posterior_summary(linFit9)
linFit9

print(linFit9, digits = 4)
plot(linFit9)
pp_check(linFit9)
loo(linFit9)
performance::check_distribution(linFit9)
performance::check_outliers(linFit9)
performance::check_heteroskedasticity(linFit9)
performance_rmse(linFit9)
performance_mae(linFit9)

variance_decomposition(linFit9)
fixef(linFit9)
ranef(linFit9)

bayes_R2(linFit9)

bayes_factor(linFit9, linFit4)
bayes_factor(linFit9, linFit5)
bayes_factor(linFit9, linFit6)

conditional_smooths(linFit9)
conditional_effects(linFit9)

### Model 10 ----
linFit10 <- brm(
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
  data = StormdataTrain7, 
  family = skew_normal(), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52, 
  warmup = 1000
)

prior_summary(linFit10)
posterior_summary(linFit10)
linFit10

print(linFit10, digits = 4)
plot(linFit10)
pp_check(linFit10, ndraws = 100)
loo(linFit10)
performance::check_distribution(linFit10)
performance::check_outliers(linFit10)
performance::check_heteroskedasticity(linFit10)
performance_rmse(linFit10)
performance_mae(linFit10)

variance_decomposition(linFit10)
fixef(linFit10)
ranef(linFit10)

bayes_R2(linFit10)

bayes_factor(linFit10, linFit4)
bayes_factor(linFit10, linFit5)
bayes_factor(linFit10, linFit9)

conditional_smooths(linFit10)
conditional_effects(linFit10)

#### Prediction ----
## Fitted
linFit10finalFit <- posterior_predict(linFit10)
linFit10finalFitMean <- colMeans(linFit10finalFit)
linFit10finalFitMed <- apply(linFit10finalFit, 2, function(x){quantile(x, 0.5)})
linFit10finalFitLCB <- apply(linFit10finalFit, 2, function(x){quantile(x, 0.025)})
linFit10finalFitUCB <- apply(linFit10finalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
linFit10finalPreds <- posterior_predict(linFit10, 
                                        newdata = StormdataTestFinal,
                                        allow_new_levels = TRUE)
linFit10finalPreds2 <- colMeans(linFit10finalPreds)
linFit10finalPredsMed <- apply(linFit10finalPreds, 2, function(x){quantile(x, 0.5)})
linFit10finalPredsLCB <- apply(linFit10finalPreds, 2, function(x){quantile(x, 0.025)})
linFit10finalPredsUCB <- apply(linFit10finalPreds, 2, function(x){quantile(x, 0.975)})

linFit10predMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(linFit10finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(linFit10finalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < linFit10finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(linFit10finalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(linFit10finalPredsMed - Actual_Yvec)),
  COV_pred = mean(linFit10finalPredsLCB < Actual_Yvec & Actual_Yvec < linFit10finalPredsUCB)
)
linFit10predMetrics

ppc_dens_overlay(y = Actual_Yvec, yrep = linFit10finalPreds) +
  labs(title = "GammaFit5 Predict")
rm(linFit10finalFit)
rm(linFit10finalPreds)

### Model 11 ----
linFit11 <- brm(
  formula = VMAX ~ 
    #Year +
    Month +
    #s(StormElapsedTime) + 
    #I(StormElapsedTime^2) +
    basin + 
    t2(LAT, LON, StormElapsedTime) +
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
  data = StormdataTrain7, 
  family = gaussian(), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52, 
  warmup = 1000
)

prior_summary(linFit11)
posterior_summary(linFit11)
linFit11

print(linFit11, digits = 4)
plot(linFit11)
pp_check(linFit11, ndraws = 100)
loo(linFit11, linFit10)
waic(linFit11)
performance::check_distribution(linFit11)
performance::check_outliers(linFit11)
performance::check_heteroskedasticity(linFit11)
performance_rmse(linFit11)
performance_mae(linFit11)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(linFit11)


variance_decomposition(linFit11)
exp(fixef(linFit11))
ranef(linFit11)

bayes_R2(linFit11)

bayes_factor(linFit11, gammaFit1)
bayes_factor(linFit11, studentFit1)
bayes_factor(linFit11, linFit10)
bayes_factor(linFit11, propFit1)
bayes_factor(linFit11, logPropFit1)

conditional_smooths(linFit11)
conditional_effects(linFit11)

#### Prediction ----
## Fitted
linFit11finalFit <- posterior_predict(linFit11)
linFit11finalFitMean <- colMeans(linFit11finalFit)
linFit11finalFitMed <- apply(linFit11finalFit, 2, function(x){quantile(x, 0.5)})
linFit11finalFitLCB <- apply(linFit11finalFit, 2, function(x){quantile(x, 0.025)})
linFit11finalFitUCB <- apply(linFit11finalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
linFit11finalPreds <- posterior_predict(linFit11, 
                                        newdata = StormdataTestFinal,
                                        allow_new_levels = TRUE)
linFit11finalPreds2 <- colMeans(linFit11finalPreds)
linFit11finalPredsMed <- apply(linFit11finalPreds, 2, function(x){quantile(x, 0.5)})
linFit11finalPredsLCB <- apply(linFit11finalPreds, 2, function(x){quantile(x, 0.025)})
linFit11finalPredsUCB <- apply(linFit11finalPreds, 2, function(x){quantile(x, 0.975)})

linFit11predMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(linFit11finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(linFit11finalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < linFit11finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(linFit11finalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(linFit11finalPredsMed - Actual_Yvec)),
  COV_pred = mean(linFit11finalPredsLCB < Actual_Yvec & Actual_Yvec < linFit11finalPredsUCB)
)
linFit11predMetrics

ppc_dens_overlay(y = Actual_Yvec, yrep = linFit11finalPreds) +
  labs(title = "GammaFit5 Predict")
rm(linFit11finalFit)
rm(linFit11finalPreds)

### Model 12 ----
linFit12 <- brm(
  formula = VMAX ~ 
    #Year +
    #Month +
    #basin + 
    t2(LAT, LON, StormElapsedTime, d = c(2,1), bs = c("tp", "cr")) +
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
  data = StormdataTrain7, 
  family = gaussian(), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52, 
  warmup = 1000
)
save(linFit12, file = "_data/linFit12.RData")
prior_summary(linFit12)
posterior_summary(linFit12)
linFit12

print(linFit12, digits = 4)
plot(linFit12)
pp_check(linFit12, ndraws = 100)
loo(linFit12)
waic(linFit12)
performance::check_distribution(linFit12)
performance::check_outliers(linFit12)
performance::check_heteroskedasticity(linFit12)
performance_rmse(linFit12)
performance_mae(linFit12)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(linFit12)


variance_decomposition(linFit12)
exp(fixef(linFit12))
ranef(linFit12)

bayes_R2(linFit12)

bayes_factor(linFit12, gammaFit1)
bayes_factor(linFit12, studentFit1)
bayes_factor(linFit12, linFit10)
bayes_factor(linFit12, propFit1)
bayes_factor(linFit12, logPropFit1)

conditional_smooths(linFit12)
conditional_effects(linFit12)

linFit12effects <- conditional_effects(linFit12, 
                                       method = "posterior_predict",
                                       robust = FALSE,
                                       re_formula = NULL)

plot(linFit12effects, points = TRUE)

#### Prediction ----
## Fitted
linFit12finalFit <- posterior_predict(linFit12)
linFit12finalFitMean <- colMeans(linFit12finalFit)
linFit12finalFitMed <- apply(linFit12finalFit, 2, function(x){quantile(x, 0.5)})
linFit12finalFitLCB <- apply(linFit12finalFit, 2, function(x){quantile(x, 0.025)})
linFit12finalFitUCB <- apply(linFit12finalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
linFit12finalPreds <- posterior_predict(linFit12, 
                                        newdata = StormdataTestFinal,
                                        allow_new_levels = TRUE)
linFit12finalPreds2 <- colMeans(linFit12finalPreds)
linFit12finalPredsMed <- apply(linFit12finalPreds, 2, function(x){quantile(x, 0.5)})
linFit12finalPredsLCB <- apply(linFit12finalPreds, 2, function(x){quantile(x, 0.025)})
linFit12finalPredsUCB <- apply(linFit12finalPreds, 2, function(x){quantile(x, 0.975)})

linFit12predMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(linFit12finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(linFit12finalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < linFit12finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(linFit12finalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(linFit12finalPredsMed - Actual_Yvec)),
  COV_pred = mean(linFit12finalPredsLCB < Actual_Yvec & Actual_Yvec < linFit12finalPredsUCB)
)
linFit12predMetrics

ppc_dens_overlay(y = Actual_Yvec, yrep = linFit12finalPreds) +
  labs(title = "GammaFit5 Predict")
rm(linFit12finalFit)
rm(linFit12finalPreds)

### Model 13 ----
linFit13 <- brm(
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
  family = skew_normal(), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000,
  knots = list(
    Day = c(scale(0,
                  center = attr(StormdataTrain7scale |> pull(Day), "scaled:center"),
                  scale = attr(StormdataTrain7scale |> pull(Day), "scaled:scale")),
            scale(365,
                  center = attr(StormdataTrain7scale |> pull(Day), "scaled:center"),
                  scale = attr(StormdataTrain7scale |> pull(Day), "scaled:scale"))
    ))
)
save(linFit13, file = "_data/linFit13.RData")
prior_summary(linFit13)
round(posterior_summary(linFit13, probs = c(0.025, 0.975)))
linFit13

print(linFit13, digits = 4)
plot(linFit13)
pp_check(linFit13, ndraws = 100)
loo(linFit13)
waic(linFit13)
performance::check_distribution(linFit13)
performance::check_outliers(linFit13)
performance::check_heteroskedasticity(linFit13)
performance_rmse(linFit13)
performance_mae(linFit13)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(linFit13)


variance_decomposition(linFit13)
exp(fixef(linFit13))
ranef(linFit13)

bayes_R2(linFit13)

bayes_factor(linFit13, gammaFit1)
bayes_factor(linFit13, gammaFit2)
bayes_factor(linFit13, gammaFit3)
bayes_factor(linFit13, studentFit1)
bayes_factor(linFit13, studentFit2)
bayes_factor(linFit13, studentFit3)
bayes_factor(linFit13, linFit11)
bayes_factor(linFit13, propFit1)
bayes_factor(linFit13, logPropFit1)
loo(linFit13, gammaFit3)

linFit13smooths <- conditional_smooths(linFit13)
plot(linFit13smooths, stype = "raster", ask = FALSE)
linFit13effects <- conditional_effects(linFit13, 
                                       method = "posterior_predict",
                                       robust = FALSE,
                                       re_formula = NULL)
linFit13effects <- conditional_effects(linFit13)
plot(linFit13effects, points = TRUE, ask = FALSE)

#### Prediction ----
## Fitted
linFit13finalFit <- posterior_predict(linFit13)
linFit13finalFitMean <- colMeans(linFit13finalFit)
linFit13finalFitMed <- apply(linFit13finalFit, 2, function(x){quantile(x, 0.5)})
linFit13finalFitLCB <- apply(linFit13finalFit, 2, function(x){quantile(x, 0.025)})
linFit13finalFitUCB <- apply(linFit13finalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
linFit13finalPreds <- posterior_predict(linFit13, 
                                        newdata = StormdataTestFinalscale,
                                        allow_new_levels = TRUE)
linFit13finalPreds2 <- colMeans(linFit13finalPreds)
linFit13finalPredsMed <- apply(linFit13finalPreds, 2, function(x){quantile(x, 0.5)})
linFit13finalPredsLCB <- apply(linFit13finalPreds, 2, function(x){quantile(x, 0.025)})
linFit13finalPredsUCB <- apply(linFit13finalPreds, 2, function(x){quantile(x, 0.975)})

linFit13predMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(linFit13finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(linFit13finalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < linFit13finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(linFit13finalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(linFit13finalPredsMed - Actual_Yvec)),
  COV_pred = mean(linFit13finalPredsLCB < Actual_Yvec & Actual_Yvec < linFit13finalPredsUCB)
)
linFit13predMetrics

#### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = linFit13finalPreds) +
  labs(title = "GammaFit5 Predict")

linFit13FitDF <- bind_cols(
  StormdataTrain3,
  LCB = linFit13finalFitLCB,
  Mean = linFit13finalFitMean,
  Med = linFit13finalFitMed,
  UCB = linFit13finalFitUCB
) 

ggplot(data = linFit13FitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))

## Prediction
linFit13PredDF <- bind_cols(
  StormdataTest3,
  LCB = linFit13finalPredsLCB,
  Mean = linFit13finalPreds2,
  Med = linFit13finalPredsMed,
  UCB = linFit13finalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

ggplot(data = linFit13PredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))

rm(linFit13finalFit)
rm(linFit13finalPreds)


## STUDENT T ----
### Model 1 ----
studentFit1 <- brm(
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
  data = StormdataTrain7, 
  family = student(), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52, 
  warmup = 1000
)

prior_summary(studentFit1)
posterior_summary(studentFit1)
studentFit1

print(studentFit1, digits = 4)
plot(studentFit1)
pp_check(studentFit1, ndraws = 100)
loo(studentFit1)
performance::check_distribution(studentFit1)
performance::check_outliers(studentFit1)
performance::check_heteroskedasticity(studentFit1)
performance_rmse(studentFit1)
performance_mae(studentFit1)

variance_decomposition(studentFit1)
fixef(studentFit1)
ranef(studentFit1)

bayes_R2(studentFit1)

bayes_factor(studentFit1, linFit9)
bayes_factor(studentFit1, linFit10)
bayes_factor(studentFit1, linFit11)

conditional_smooths(studentFit1)
conditional_effects(studentFit1)

#### Prediction ----
## Fitted
studentFit1finalFit <- posterior_predict(studentFit1)
studentFit1finalFitMean <- colMeans(studentFit1finalFit)
studentFit1finalFitMed <- apply(studentFit1finalFit, 2, function(x){quantile(x, 0.5)})
studentFit1finalFitLCB <- apply(studentFit1finalFit, 2, function(x){quantile(x, 0.025)})
studentFit1finalFitUCB <- apply(studentFit1finalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
studentFit1finalPreds <- posterior_predict(studentFit1, 
                                           newdata = StormdataTestFinal,
                                           allow_new_levels = TRUE)
studentFit1finalPreds2 <- colMeans(studentFit1finalPreds)
studentFit1finalPredsMed <- apply(studentFit1finalPreds, 2, function(x){quantile(x, 0.5)})
studentFit1finalPredsLCB <- apply(studentFit1finalPreds, 2, function(x){quantile(x, 0.025)})
studentFit1finalPredsUCB <- apply(studentFit1finalPreds, 2, function(x){quantile(x, 0.975)})

studentFit1predMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(studentFit1finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(studentFit1finalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < studentFit1finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(studentFit1finalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(studentFit1finalPredsMed - Actual_Yvec)),
  COV_pred = mean(studentFit1finalPredsLCB < Actual_Yvec & Actual_Yvec < studentFit1finalPredsUCB)
)
studentFit1predMetrics

#### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = studentFit1finalPreds) +
  labs(title = "studentFit1 Predict")

studentFit1FitDF <- bind_cols(
  StormdataTrain3,
  LCB = studentFit1finalFitLCB,
  Mean = studentFit1finalFitMean,
  Med = studentFit1finalFitMed,
  UCB = studentFit1finalFitUCB
)

ggplot(data = studentFit1FitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))

## Prediction
studentFit1PredDF <- bind_cols(
  StormdataTest3,
  LCB = studentFit1finalPredsLCB,
  Mean = studentFit1finalPreds2,
  Med = studentFit1finalPredsMed,
  UCB = studentFit1finalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  )

ggplot(data = studentFit1PredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))

rm(studentFit1finalFit)
rm(studentFit1finalPreds)

### Model 2 ----
studentFit2 <- brm(
  formula = VMAX ~ 
    #Year +
    Month +
    #s(StormElapsedTime) + 
    #I(StormElapsedTime^2) +
    basin + 
    t2(LAT, LON, StormElapsedTime) +
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
  data = StormdataTrain7, 
  family = student(link = "identity"), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52, 
  warmup = 1000
)

prior_summary(studentFit2)
posterior_summary(studentFit2)
studentFit2

print(studentFit2, digits = 4)
plot(studentFit2)
pp_check(studentFit2, ndraws = 100)
loo(studentFit2, linFit10)
waic(studentFit2)
performance::check_distribution(studentFit2)
performance::check_outliers(studentFit2)
performance::check_heteroskedasticity(studentFit2)
performance_rmse(studentFit2)
performance_mae(studentFit2)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(studentFit2)


variance_decomposition(studentFit2)
exp(fixef(studentFit2))
ranef(studentFit2)

bayes_R2(studentFit2)

bayes_factor(studentFit2, gammaFit1)
bayes_factor(studentFit2, studentFit1)
bayes_factor(studentFit2, linFit10)
bayes_factor(studentFit2, linFit11)
bayes_factor(studentFit2, propFit1)
bayes_factor(studentFit2, logPropFit1)

#### Prediction ----
## Fitted
studentFit2finalFit <- posterior_predict(studentFit2)
studentFit2finalFitMean <- colMeans(studentFit2finalFit)
studentFit2finalFitMed <- apply(studentFit2finalFit, 2, function(x){quantile(x, 0.5)})
studentFit2finalFitLCB <- apply(studentFit2finalFit, 2, function(x){quantile(x, 0.025)})
studentFit2finalFitUCB <- apply(studentFit2finalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
studentFit2finalPreds <- posterior_predict(studentFit2, 
                                           newdata = StormdataTestFinal,
                                           allow_new_levels = TRUE)
studentFit2finalPreds2 <- colMeans(studentFit2finalPreds)
studentFit2finalPredsMed <- apply(studentFit2finalPreds, 2, function(x){quantile(x, 0.5)})
studentFit2finalPredsLCB <- apply(studentFit2finalPreds, 2, function(x){quantile(x, 0.025)})
studentFit2finalPredsUCB <- apply(studentFit2finalPreds, 2, function(x){quantile(x, 0.975)})

studentFit2predMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(studentFit2finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(studentFit2finalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < studentFit2finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(studentFit2finalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(studentFit2finalPredsMed - Actual_Yvec)),
  COV_pred = mean(studentFit2finalPredsLCB < Actual_Yvec & Actual_Yvec < studentFit2finalPredsUCB)
)
studentFit2predMetrics

ppc_dens_overlay(y = Actual_Yvec, yrep = studentFit2finalPreds) +
  labs(title = "GammaFit5 Predict")
rm(studentFit2finalFit)
rm(studentFit2finalPreds)

### Model 3 ----
studentFit3 <- brm(
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
  family = student(link = "log"), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52, 
  warmup = 1000
)

prior_summary(studentFit3)
posterior_summary(studentFit3)
studentFit3

print(studentFit3, digits = 4)
plot(studentFit3)
pp_check(studentFit3, ndraws = 100) + xlim(c(-20, 300))
loo(studentFit3, studentFit2)
waic(studentFit3)
performance::check_distribution(studentFit3)
performance::check_outliers(studentFit3)
performance::check_heteroskedasticity(studentFit3)
performance_rmse(studentFit3)
performance_mae(studentFit3)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(studentFit3)


variance_decomposition(studentFit3)
exp(fixef(studentFit3))
ranef(studentFit3)

bayes_R2(studentFit3)

bayes_factor(studentFit3, studentFit1)
bayes_factor(studentFit3, studentFit2)

conditional_smooths(studentFit3)
conditional_effects(studentFit3)

#### Prediction ----
## Fitted
studentFit3finalFit <- posterior_predict(studentFit3)
studentFit3finalFitMean <- colMeans(studentFit3finalFit)
studentFit3finalFitMed <- apply(studentFit3finalFit, 2, function(x){quantile(x, 0.5)})
studentFit3finalFitLCB <- apply(studentFit3finalFit, 2, function(x){quantile(x, 0.025)})
studentFit3finalFitUCB <- apply(studentFit3finalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
studentFit3finalPreds <- posterior_predict(studentFit3, 
                                           newdata = StormdataTestFinal,
                                           allow_new_levels = TRUE)
studentFit3finalPreds2 <- colMeans(studentFit3finalPreds)
studentFit3finalPredsMed <- apply(studentFit3finalPreds, 2, function(x){quantile(x, 0.5)})
studentFit3finalPredsLCB <- apply(studentFit3finalPreds, 2, function(x){quantile(x, 0.025)})
studentFit3finalPredsUCB <- apply(studentFit3finalPreds, 2, function(x){quantile(x, 0.975)})

studentFit3predMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(studentFit3finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(studentFit3finalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < studentFit3finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(studentFit3finalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(studentFit3finalPredsMed - Actual_Yvec)),
  COV_pred = mean(studentFit3finalPredsLCB < Actual_Yvec & Actual_Yvec < studentFit3finalPredsUCB)
)
studentFit3predMetrics

ppc_dens_overlay(y = Actual_Yvec, yrep = studentFit3finalPreds) +
  labs(title = "GammaFit5 Predict")
rm(studentFit3finalFit)
rm(studentFit3finalPreds)

## GAMMA ----
### Model 0 ----
gammaFit0 <- brm(
  formula = VMAX ~ 
    #Year +
    Month +
    basin + 
    LAT +
    LON +
    StormElapsedTime +
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
  data = StormdataTrain7, 
  family = brmsfamily("Gamma", link = "log"), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)
save(gammaFit0, file = "_data/gammaFit0.RData")
prior_summary(gammaFit0)
posterior_summary(gammaFit0)
summary(gammaFit0)
gammaFit0

print(gammaFit0, digits = 4)
plot(gammaFit0)
pp_check(gammaFit0, ndraws = 100)
loo(gammaFit0)
waic(gammaFit0)
performance::check_distribution(gammaFit0)
performance::check_outliers(gammaFit0)
performance::check_heteroskedasticity(gammaFit0)
check_overdispersion(gammaFit0)
check_autocorrelation(gammaFit0)
performance_rmse(gammaFit0)
performance_mae(gammaFit0)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(gammaFit0)

gammaFit0residuals <- residuals(gammaFit0, method = "posterior_predict")


variance_decomposition(gammaFit0)
fixef(gammaFit0)
ranef(gammaFit0)

bayes_R2(gammaFit0)

bayes_factor(gammaFit0, gammaFit1)
bayes_factor(gammaFit0, gammaFit2)
bayes_factor(gammaFit0, studentFit1)
bayes_factor(gammaFit0, studentFit2)
bayes_factor(gammaFit0, studentFit0)
bayes_factor(gammaFit0, linFit11)
bayes_factor(gammaFit0, propFit1)
bayes_factor(gammaFit0, logPropFit1)
loo(gammaFit1, gammaFit2, gammaFit0)

gammaFit0smooths <- conditional_smooths(gammaFit0)
plot(gammaFit0smooths, stype = "raster")
gammaFit0effects <- conditional_effects(gammaFit0, 
                                        effects = c(
                                          "LON",
                                          "LAT",
                                          "StormElapsedTime"
                                        ),
                                        surface = TRUE)
gammaFit0effects
conditional_effects(gammaFit0, surface = TRUE, robust = FALSE)

#### Prediction ----
## Fitted
gammaFit0finalFit <- posterior_predict(gammaFit0)
gammaFit0finalFitMean <- colMeans(gammaFit0finalFit)
gammaFit0finalFitMed <- apply(gammaFit0finalFit, 2, function(x){quantile(x, 0.5)})
gammaFit0finalFitLCB <- apply(gammaFit0finalFit, 2, function(x){quantile(x, 0.025)})
gammaFit0finalFitUCB <- apply(gammaFit0finalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
gammaFit0finalPreds <- posterior_predict(gammaFit0, 
                                         newdata = StormdataTestFinal,
                                         allow_new_levels = TRUE)
gammaFit0finalPreds2 <- colMeans(gammaFit0finalPreds)
gammaFit0finalPredsMed <- apply(gammaFit0finalPreds, 2, function(x){quantile(x, 0.5)})
gammaFit0finalPredsLCB <- apply(gammaFit0finalPreds, 2, function(x){quantile(x, 0.025)})
gammaFit0finalPredsUCB <- apply(gammaFit0finalPreds, 2, function(x){quantile(x, 0.975)})

gammaFit0predMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(gammaFit0finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(gammaFit0finalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < gammaFit0finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(gammaFit0finalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(gammaFit0finalPredsMed - Actual_Yvec)),
  COV_pred = mean(gammaFit0finalPredsLCB < Actual_Yvec & Actual_Yvec < gammaFit0finalPredsUCB)
)
gammaFit0predMetrics

ppc_dens_overlay(y = Actual_Yvec, yrep = gammaFit0finalPreds) +
  labs(title = "GammaFit0 Predict")
rm(gammaFit0finalFit)
rm(gammaFit0finalPreds)

### Model 1 ----
gammaFit1 <- brm(
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
  data = StormdataTrain7, 
  family = brmsfamily("Gamma", link = "log"), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52, 
  warmup = 1000
)

prior_summary(gammaFit1)
posterior_summary(gammaFit1)
gammaFit1

print(gammaFit1, digits = 4)
plot(gammaFit1)
pp_check(gammaFit1, ndraws = 100)
loo(gammaFit1, studentFit1)
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

bayes_factor(gammaFit1, gammaFit4)
bayes_factor(gammaFit1, gammaFit5)
bayes_factor(gammaFit1, studentFit1)

conditional_smooths(gammaFit1)
conditional_effects(gammaFit1)

gammapostEPreds1 <- posterior_epred(gammaFit1)
gammapostEPreds2 <- colMeans(gammapostEPreds1)
mean(abs(gammapostEPreds2 - StormdataTrain3$VMAX))

#### Cross Validation ----
set.seed(52)
folds <- kfold_split_random(K = 5, N = nrow(StormdataTrain7))
foldFits <- list()
foldDataList <- list()
foldPredList <- list()
tic()
for(i in 1:5){
  foldIndex <- folds == i
  foldData <- StormdataTrain7 |> filter(!foldIndex)
  foldDataList[[i]] <- foldData
  foldPredList[[i]] <- StormdataTrain7 |> filter(foldIndex)
  foldFits[[i]] <- brm(
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
    data = foldData, 
    family = brmsfamily("Gamma", link = "log"), 
    save_pars = save_pars(all = TRUE), 
    chains = 4,
    iter = 2000,
    seed = 52, 
    warmup = 1000
  )
}
toc()

totalPreds <- matrix( , nrow = 4000, ncol = 0)
predMeans <- c()
predLower <- c()
predUpper <- c()
set.seed(52)
for(j in 1:5){
  predDat <- foldPredList[[j]] |>
    select(-VMAX)
  tempPred <- posterior_predict(foldFits[[j]], 
                                newdata = predDat,
                                allow_new_levels = TRUE)
  Obs_Y <- foldPredList[[j]] |> pull(VMAX)
  predMeansT <- abs(apply(tempPred, 2, function(x){mean(x)}) - Obs_Y)
  predLowerT <- apply(tempPred, 2, function(x){quantile(x, 0.025)}) < Obs_Y
  predUpperT <- apply(tempPred, 2, function(x){quantile(x, 0.975)}) > Obs_Y
  predMeans <- c(predMeans, predMeansT)
  predLower <- c(predLower, predLowerT)
  predUpper <- c(predUpper, predUpperT)
}

#### MAE and COV ----
cvMAE <- mean(predMeans)
cvCov <- mean(predLower & predUpper)

cvMetrics <- tibble(
  MAE_HWRF = mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF)),
  MAE_fit = cvMAE,
  COV = cvCov
)
cvMetrics


##### Using loo ----
logLikGamma <- log_lik(gammaFit1)

# foldsIndex1 <- which(folds == 1)
# foldsIndex2 <- which(folds == 2)
# foldsIndex3 <- which(folds == 3)
# foldsIndex4 <- which(folds == 4)
# foldsIndex5 <- which(folds == 5)
# 
# totalfoldIndex <- c(foldsIndex5,
#                     foldsIndex4,
#                     foldsIndex3,
#                     foldsIndex2,
#                     foldsIndex1)
# 
# totalPredsDF <- data.frame(totalPreds) |> select(totalfoldIndex)
# predMeans <- apply(totalPredsDF, 2, function(x){mean(x)})
# predLower <- apply(totalPredsDF, 2, function(x){quantile(x, 0.025)})
# predUpper <- apply(totalPredsDF, 2, function(x){quantile(x, 0.975)})
# 
# predDF <- StormdataTrain3 |>
#   select(HWRF, VMAX) |>
#   bind_cols(predMeans = predMeans,
#             predLower = predLower, 
#             predUpper = predUpper)
# 
# predSumDF <- predDF |>
#   summarise(
#     MAE_HWRF = mean(abs(VMAX - HWRF)),
#     MAE_fit = mean(abs(VMAX - predMeans)),
#     COV = mean(between(VMAX, predLower, predUpper))
#   )
# predSumDF

#### Variable Selection ----
varSelGamma <- varsel(gammaFit1)

#### Prediction ----
## Fitted
gammaFit1finalFit <- posterior_predict(gammaFit1)
gammaFit1finalFitMean <- colMeans(gammaFit1finalFit)
gammaFit1finalFitMed <- apply(gammaFit1finalFit, 2, function(x){quantile(x, 0.5)})
gammaFit1finalFitLCB <- apply(gammaFit1finalFit, 2, function(x){quantile(x, 0.025)})
gammaFit1finalFitUCB <- apply(gammaFit1finalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
gammaFit1finalPreds <- posterior_predict(gammaFit1, 
                                         newdata = StormdataTestFinal,
                                         allow_new_levels = TRUE)
gammaFit1finalPreds2 <- colMeans(gammaFit1finalPreds)
gammaFit1finalPredsMed <- apply(gammaFit1finalPreds, 2, function(x){quantile(x, 0.5)})
gammaFit1finalPredsLCB <- apply(gammaFit1finalPreds, 2, function(x){quantile(x, 0.025)})
gammaFit1finalPredsUCB <- apply(gammaFit1finalPreds, 2, function(x){quantile(x, 0.975)})

gammaFit1predMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(gammaFit1finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(gammaFit1finalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < gammaFit1finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(gammaFit1finalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(gammaFit1finalPredsMed - Actual_Yvec)),
  COV_pred = mean(gammaFit1finalPredsLCB < Actual_Yvec & Actual_Yvec < gammaFit1finalPredsUCB)
)
gammaFit1predMetrics

ppc_dens_overlay(y = Actual_Yvec, yrep = gammaFit1finalPreds) +
  labs(title = "GammaFit1 Predict")
rm(gammaFit1finalFit)
rm(gammaFit1finalPreds)



### Model 2 ----
gammaFit2 <- brm(
  formula = VMAX ~ 
    #Year +
    Month +
    #s(StormElapsedTime) + 
    #I(StormElapsedTime^2) +
    basin + 
    t2(LAT, LON, StormElapsedTime) +
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
  data = StormdataTrain7, 
  family = brmsfamily("Gamma", link = "log"), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52, 
  warmup = 1000
)

prior_summary(gammaFit2)
posterior_summary(gammaFit2)
gammaFit2

print(gammaFit2, digits = 4)
plot(gammaFit2)
pp_check(gammaFit2, ndraws = 100)
loo(gammaFit2, gammaFit1)
waic(gammaFit2)
performance::check_distribution(gammaFit2)
performance::check_outliers(gammaFit2)
performance::check_heteroskedasticity(gammaFit2)
performance_rmse(gammaFit2)
performance_mae(gammaFit2)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(gammaFit2)


variance_decomposition(gammaFit2)
exp(fixef(gammaFit2))
ranef(gammaFit2)

bayes_R2(gammaFit2)

bayes_factor(gammaFit2, gammaFit1)
bayes_factor(gammaFit2, studentFit1)
bayes_factor(gammaFit2, linFit11)
bayes_factor(gammaFit2, propFit1)
bayes_factor(gammaFit2, logPropFit1)

conditional_smooths(gammaFit2)
conditional_effects(gammaFit2)

#### Prediction ----
## Fitted
gammaFit2finalFit <- posterior_predict(gammaFit2)
gammaFit2finalFitMean <- colMeans(gammaFit2finalFit)
gammaFit2finalFitMed <- apply(gammaFit2finalFit, 2, function(x){quantile(x, 0.5)})
gammaFit2finalFitLCB <- apply(gammaFit2finalFit, 2, function(x){quantile(x, 0.025)})
gammaFit2finalFitUCB <- apply(gammaFit2finalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
gammaFit2finalPreds <- posterior_predict(gammaFit2, 
                                         newdata = StormdataTestFinal,
                                         allow_new_levels = TRUE)
gammaFit2finalPreds2 <- colMeans(gammaFit2finalPreds)
gammaFit2finalPredsMed <- apply(gammaFit2finalPreds, 2, function(x){quantile(x, 0.5)})
gammaFit2finalPredsLCB <- apply(gammaFit2finalPreds, 2, function(x){quantile(x, 0.025)})
gammaFit2finalPredsUCB <- apply(gammaFit2finalPreds, 2, function(x){quantile(x, 0.975)})

gammaFit2predMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(gammaFit2finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(gammaFit2finalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < gammaFit2finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(gammaFit2finalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(gammaFit2finalPredsMed - Actual_Yvec)),
  COV_pred = mean(gammaFit2finalPredsLCB < Actual_Yvec & Actual_Yvec < gammaFit2finalPredsUCB)
)
gammaFit2predMetrics

ppc_dens_overlay(y = Actual_Yvec, yrep = gammaFit2finalPreds) +
  labs(title = "GammaFit2 Predict")
rm(gammaFit2finalFit)
rm(gammaFit2finalPreds)

### Model 3 ----
gammaFit3 <- brm(
  formula = VMAX ~ 
    #Year +
    Month +
    basin + 
    t2(LON, LAT, StormElapsedTime, d = c(2,1)) +
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
  data = StormdataTrain7, 
  family = brmsfamily("Gamma", link = "log"), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)

prior_summary(gammaFit3)
posterior_summary(gammaFit3)
summary(gammaFit3)

print(gammaFit3, digits = 4)
plot(gammaFit3)
pp_check(gammaFit3, ndraws = 40)
loo(gammaFit3)
waic(gammaFit3)
performance::check_distribution(gammaFit3)
performance::check_outliers(gammaFit3)
performance::check_heteroskedasticity(gammaFit3)
performance_rmse(gammaFit3)
performance_mae(gammaFit3)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(gammaFit3)

gammaFit3residuals <- residuals(gammaFit3, method = "posterior_predict")


variance_decomposition(gammaFit3)
fixef(gammaFit3)
ranef(gammaFit3)

bayes_R2(gammaFit3)

bayes_factor(gammaFit3, gammaFit1)
bayes_factor(gammaFit3, gammaFit2)
bayes_factor(gammaFit3, studentFit1)
bayes_factor(gammaFit3, studentFit2)
bayes_factor(gammaFit3, studentFit3)
bayes_factor(gammaFit3, linFit11)
bayes_factor(gammaFit3, propFit1)
bayes_factor(gammaFit3, logPropFit1)
loo(gammaFit0, gammaFit1, gammaFit2, gammaFit3)

gammaFit3smooths <- conditional_smooths(gammaFit3)
plot(gammaFit3smooths, stype = "raster")
gammaFit3effects <- conditional_effects(gammaFit3, 
                                        effects = c(
                                          "LON",
                                          "LAT",
                                          "StormElapsedTime"
                                        ),
                                        surface = TRUE)
gammaFit3effects
conditional_effects(gammaFit3, surface = TRUE, robust = FALSE,)

#### Prediction ----
## Fitted
gammaFit3finalFit <- posterior_predict(gammaFit3)
gammaFit3finalFitMean <- colMeans(gammaFit3finalFit)
gammaFit3finalFitMed <- apply(gammaFit3finalFit, 2, function(x){quantile(x, 0.5)})
gammaFit3finalFitLCB <- apply(gammaFit3finalFit, 2, function(x){quantile(x, 0.025)})
gammaFit3finalFitUCB <- apply(gammaFit3finalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
gammaFit3finalPreds <- posterior_predict(gammaFit3, 
                                         newdata = StormdataTestFinal,
                                         allow_new_levels = TRUE)
gammaFit3finalPreds2 <- colMeans(gammaFit3finalPreds)
gammaFit3finalPredsMed <- apply(gammaFit3finalPreds, 2, function(x){quantile(x, 0.5)})
gammaFit3finalPredsLCB <- apply(gammaFit3finalPreds, 2, function(x){quantile(x, 0.025)})
gammaFit3finalPredsUCB <- apply(gammaFit3finalPreds, 2, function(x){quantile(x, 0.975)})

gammaFit3predMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(gammaFit3finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(gammaFit3finalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < gammaFit3finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(gammaFit3finalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(gammaFit3finalPredsMed - Actual_Yvec)),
  COV_pred = mean(gammaFit3finalPredsLCB < Actual_Yvec & Actual_Yvec < gammaFit3finalPredsUCB)
)
gammaFit3predMetrics

#### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = gammaFit3finalPreds) +
  labs(title = "GammaFit3 Predict")

gammaFit3FitDF <- bind_cols(
  StormdataTrain3,
  LCB = gammaFit3finalFitLCB,
  Mean = gammaFit3finalFitMean,
  Med = gammaFit3finalFitMed,
  UCB = gammaFit3finalFitUCB
) 

ggplot(data = gammaFit3FitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))

## Prediction
gammaFit3PredDF <- bind_cols(
  StormdataTest3,
  LCB = gammaFit3finalPredsLCB,
  Mean = gammaFit3finalPreds2,
  Med = gammaFit3finalPredsMed,
  UCB = gammaFit3finalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

ggplot(data = gammaFit3PredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))

rm(gammaFit3finalFit)
rm(gammaFit3finalPreds)

### Model 4 ----
gammaFit4 <- brm(
  formula = VMAX ~
    #Year +
    #Month +
    #basin + 
    t2(LON, LAT, StormElapsedTime, d = c(2,1), bs = c("tp", "cr")) +
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
  data = StormdataTrain7, 
  family = brmsfamily("Gamma", link = "log"), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)
save(gammaFit4, file = "_data/gammaFit4.RData")
prior_summary(gammaFit4)
posterior_summary(gammaFit4)
gammaFit4

print(gammaFit4, digits = 4)
plot(gammaFit4)
pp_check(gammaFit4, ndraws = 100)
loo(gammaFit4)
waic(gammaFit4)
performance::check_distribution(gammaFit4)
performance::check_outliers(gammaFit4)
performance::check_heteroskedasticity(gammaFit4)
performance_rmse(gammaFit4)
performance_mae(gammaFit4)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(gammaFit4)


variance_decomposition(gammaFit4)
exp(fixef(gammaFit4))
ranef(gammaFit4)

bayes_R2(gammaFit4)

bayes_factor(gammaFit4, gammaFit1)
bayes_factor(gammaFit4, gammaFit2)
bayes_factor(gammaFit4, gammaFit3)
bayes_factor(gammaFit4, studentFit1)
bayes_factor(gammaFit4, studentFit2)
bayes_factor(gammaFit4, studentFit3)
bayes_factor(gammaFit4, linFit11)
bayes_factor(gammaFit4, propFit1)
bayes_factor(gammaFit4, logPropFit1)
loo(gammaFit4, gammaFit3)

gammaFit4smooths <- conditional_smooths(gammaFit4, spaghetti = TRUE, ndraws = 100,)
plot(gammaFit4smooths, stype = "raster")
gammaFit4effects <- conditional_effects(gammaFit4, 
                                        method = "posterior_predict",
                                        robust = FALSE,
                                        re_formula = NULL)
#spaghetti = TRUE,
#ndraws = 100)
plot(gammaFit4effects, points = TRUE)

#### Prediction ----
## Fitted
gammaFit4finalFit <- posterior_predict(gammaFit4)
gammaFit4finalFitMean <- colMeans(gammaFit4finalFit)
gammaFit4finalFitMed <- apply(gammaFit4finalFit, 2, function(x){quantile(x, 0.5)})
gammaFit4finalFitLCB <- apply(gammaFit4finalFit, 2, function(x){quantile(x, 0.025)})
gammaFit4finalFitUCB <- apply(gammaFit4finalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
gammaFit4finalPreds <- posterior_predict(gammaFit4, 
                                         newdata = StormdataTestFinal,
                                         allow_new_levels = TRUE)
gammaFit4finalPreds2 <- colMeans(gammaFit4finalPreds)
gammaFit4finalPredsMed <- apply(gammaFit4finalPreds, 2, function(x){quantile(x, 0.5)})
gammaFit4finalPredsLCB <- apply(gammaFit4finalPreds, 2, function(x){quantile(x, 0.025)})
gammaFit4finalPredsUCB <- apply(gammaFit4finalPreds, 2, function(x){quantile(x, 0.975)})

gammaFit4predMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(gammaFit4finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(gammaFit4finalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < gammaFit4finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(gammaFit4finalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(gammaFit4finalPredsMed - Actual_Yvec)),
  COV_pred = mean(gammaFit4finalPredsLCB < Actual_Yvec & Actual_Yvec < gammaFit4finalPredsUCB)
)
gammaFit4predMetrics

ppc_dens_overlay(y = Actual_Yvec, yrep = gammaFit4finalPreds) +
  labs(title = "GammaFit4 Predict")
rm(gammaFit4finalFit)
rm(gammaFit4finalPreds)


### Model 5 ----
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

### Model 5B ----
gammaFit5B <- brm(
  formula = VMAX ~
    #Year +
    #Month +
    basin + 
    s(Day, bs = "cc") +
    s(StormElapsedTime) + 
    s(LON, LAT) +
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
  warmup = 1000,
  knots = list(
    Day = c(scale(0,
                  center = attr(StormdataTrain7scale |> pull(Day), "scaled:center"),
                  scale = attr(StormdataTrain7scale |> pull(Day), "scaled:scale")),
            scale(365,
                  center = attr(StormdataTrain7scale |> pull(Day), "scaled:center"),
                  scale = attr(StormdataTrain7scale |> pull(Day), "scaled:scale"))
    ))
)
save(gammaFit5B, file = "_data/gammaFit5B.RData")
prior_summary(gammaFit5B)
round(posterior_summary(gammaFit5B, probs = c(0.025, 0.975)))
gammaFit5B

print(gammaFit5B, digits = 4)
plot(gammaFit5B)
pp_check(gammaFit5B, ndraws = 100)
loo(gammaFit5B)
waic(gammaFit5B)
performance::check_distribution(gammaFit5B)
performance::check_outliers(gammaFit5B)
performance::check_heteroskedasticity(gammaFit5B)
performance_rmse(gammaFit5B)
performance_mae(gammaFit5B)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(gammaFit5B)


variance_decomposition(gammaFit5B)
exp(fixef(gammaFit5B))
ranef(gammaFit5B)

bayes_R2(gammaFit5B)

bayes_factor(gammaFit5B, gammaFit1)
bayes_factor(gammaFit5B, gammaFit2)
bayes_factor(gammaFit5B, gammaFit3)
bayes_factor(gammaFit5B, studentFit1)
bayes_factor(gammaFit5B, studentFit2)
bayes_factor(gammaFit5B, studentFit3)
bayes_factor(gammaFit5B, linFit11)
bayes_factor(gammaFit5B, propFit1)
bayes_factor(gammaFit5B, logPropFit1)
loo(gammaFit5B, gammaFit3)

gammaFit5Bsmooths <- conditional_smooths(gammaFit5B)
plot(gammaFit5Bsmooths, stype = "raster", ask = FALSE)
gammaFit5Beffects <- conditional_effects(gammaFit5B, 
                                         method = "posterior_predict",
                                         robust = FALSE,
                                         re_formula = NULL)
gammaFit5Beffects <- conditional_effects(gammaFit5B)
plot(gammaFit5Beffects, points = TRUE, ask = FALSE)

#### Prediction ----
## Fitted
gammaFit5BfinalFit <- posterior_predict(gammaFit5B)
gammaFit5BfinalFitMean <- colMeans(gammaFit5BfinalFit)
gammaFit5BfinalFitMed <- apply(gammaFit5BfinalFit, 2, function(x){quantile(x, 0.5)})
gammaFit5BfinalFitLCB <- apply(gammaFit5BfinalFit, 2, function(x){quantile(x, 0.025)})
gammaFit5BfinalFitUCB <- apply(gammaFit5BfinalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
gammaFit5BfinalPreds <- posterior_predict(gammaFit5B, 
                                          newdata = StormdataTestFinalscale,
                                          allow_new_levels = TRUE)
gammaFit5BfinalPreds2 <- colMeans(gammaFit5BfinalPreds)
gammaFit5BfinalPredsMed <- apply(gammaFit5BfinalPreds, 2, function(x){quantile(x, 0.5)})
gammaFit5BfinalPredsLCB <- apply(gammaFit5BfinalPreds, 2, function(x){quantile(x, 0.025)})
gammaFit5BfinalPredsUCB <- apply(gammaFit5BfinalPreds, 2, function(x){quantile(x, 0.975)})

gammaFit5BpredMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(gammaFit5BfinalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(gammaFit5BfinalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < gammaFit5BfinalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(gammaFit5BfinalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(gammaFit5BfinalPredsMed - Actual_Yvec)),
  COV_pred = mean(gammaFit5BfinalPredsLCB < Actual_Yvec & Actual_Yvec < gammaFit5BfinalPredsUCB)
)
gammaFit5BpredMetrics

#### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = gammaFit5BfinalPreds) +
  labs(title = "GammaFit5B Predict")

gammaFit5BFitDF <- bind_cols(
  StormdataTrain3,
  LCB = gammaFit5BfinalFitLCB,
  Mean = gammaFit5BfinalFitMean,
  Med = gammaFit5BfinalFitMed,
  UCB = gammaFit5BfinalFitUCB
) 

ggplot(data = gammaFit5BFitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))

## Prediction
gammaFit5BPredDF <- bind_cols(
  StormdataTest3,
  LCB = gammaFit5BfinalPredsLCB,
  Mean = gammaFit5BfinalPreds2,
  Med = gammaFit5BfinalPredsMed,
  UCB = gammaFit5BfinalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

ggplot(data = gammaFit5BPredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))

rm(gammaFit5BfinalFit)
rm(gammaFit5BfinalPreds)

### Model 5C ----
gammaFit5C <- brm(
  formula = VMAX ~
    #Year +
    #Month +
    basin + 
    s(Day, bs = "cc") +
    s(StormElapsedTime) + 
    s(LON, LAT) +
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
  warmup = 1000,
  knots = list(
    Day = c(scale(0,
                  center = attr(StormdataTrain7scale |> pull(Day), "scaled:center"),
                  scale = attr(StormdataTrain7scale |> pull(Day), "scaled:scale")),
            scale(365,
                  center = attr(StormdataTrain7scale |> pull(Day), "scaled:center"),
                  scale = attr(StormdataTrain7scale |> pull(Day), "scaled:scale"))
    ))
)
save(gammaFit5C, file = "_data/gammaFit5C.RData")
prior_summary(gammaFit5C)
round(posterior_summary(gammaFit5C, probs = c(0.025, 0.975)))
gammaFit5C

print(gammaFit5C, digits = 4)
plot(gammaFit5C)
pp_check(gammaFit5C, ndraws = 100)
loo(gammaFit5C)
waic(gammaFit5C)
performance::check_distribution(gammaFit5C)
performance::check_outliers(gammaFit5C)
performance::check_heteroskedasticity(gammaFit5C)
performance_rmse(gammaFit5C)
performance_mae(gammaFit5C)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(gammaFit5C)


variance_decomposition(gammaFit5C)
exp(fixef(gammaFit5C))
ranef(gammaFit5C)

bayes_R2(gammaFit5C)

bayes_factor(gammaFit5C, gammaFit1)
bayes_factor(gammaFit5C, gammaFit2)
bayes_factor(gammaFit5C, gammaFit3)
bayes_factor(gammaFit5C, studentFit1)
bayes_factor(gammaFit5C, studentFit2)
bayes_factor(gammaFit5C, studentFit3)
bayes_factor(gammaFit5C, linFit11)
bayes_factor(gammaFit5C, propFit1)
bayes_factor(gammaFit5C, logPropFit1)
loo(gammaFit5C, gammaFit3)

gammaFit5Csmooths <- conditional_smooths(gammaFit5C)
plot(gammaFit5Csmooths, stype = "raster", ask = FALSE)
gammaFit5Ceffects <- conditional_effects(gammaFit5C, 
                                         method = "posterior_predict",
                                         robust = FALSE,
                                         re_formula = NULL)
gammaFit5Ceffects <- conditional_effects(gammaFit5C)
plot(gammaFit5Ceffects, points = TRUE, ask = FALSE)

#### Prediction ----
## Fitted
gammaFit5CfinalFit <- posterior_predict(gammaFit5C)
gammaFit5CfinalFitMean <- colMeans(gammaFit5CfinalFit)
gammaFit5CfinalFitMed <- apply(gammaFit5CfinalFit, 2, function(x){quantile(x, 0.5)})
gammaFit5CfinalFitLCB <- apply(gammaFit5CfinalFit, 2, function(x){quantile(x, 0.025)})
gammaFit5CfinalFitUCB <- apply(gammaFit5CfinalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
gammaFit5CfinalPreds <- posterior_predict(gammaFit5C, 
                                          newdata = StormdataTestFinalscale,
                                          allow_new_levels = TRUE)
gammaFit5CfinalPreds2 <- colMeans(gammaFit5CfinalPreds)
gammaFit5CfinalPredsMed <- apply(gammaFit5CfinalPreds, 2, function(x){quantile(x, 0.5)})
gammaFit5CfinalPredsLCB <- apply(gammaFit5CfinalPreds, 2, function(x){quantile(x, 0.025)})
gammaFit5CfinalPredsUCB <- apply(gammaFit5CfinalPreds, 2, function(x){quantile(x, 0.975)})

gammaFit5CpredMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(gammaFit5CfinalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(gammaFit5CfinalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < gammaFit5CfinalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(gammaFit5CfinalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(gammaFit5CfinalPredsMed - Actual_Yvec)),
  COV_pred = mean(gammaFit5CfinalPredsLCB < Actual_Yvec & Actual_Yvec < gammaFit5CfinalPredsUCB)
)
gammaFit5CpredMetrics

#### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = gammaFit5CfinalPreds) +
  labs(title = "GammaFit5C Predict")

gammaFit5CFitDF <- bind_cols(
  StormdataTrain3,
  LCB = gammaFit5CfinalFitLCB,
  Mean = gammaFit5CfinalFitMean,
  Med = gammaFit5CfinalFitMed,
  UCB = gammaFit5CfinalFitUCB
) 

ggplot(data = gammaFit5CFitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))

## Prediction
gammaFit5CPredDF <- bind_cols(
  StormdataTest3,
  LCB = gammaFit5CfinalPredsLCB,
  Mean = gammaFit5CfinalPreds2,
  Med = gammaFit5CfinalPredsMed,
  UCB = gammaFit5CfinalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

ggplot(data = gammaFit5CPredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))

rm(gammaFit5CfinalFit)
rm(gammaFit5CfinalPreds)

### Model 5D ----
gammaFit5D <- brm(
  formula = VMAX ~
    #Year +
    #Month +
    basin + 
    s(Day) +
    s(StormElapsedTime) + 
    s(StormElapsedTime, StormID, bs = "fs") + 
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
  warmup = 1000,
  knots = list(
    Day = c(scale(0,
                  center = attr(StormdataTrain7scale |> pull(Day), "scaled:center"),
                  scale = attr(StormdataTrain7scale |> pull(Day), "scaled:scale")),
            scale(365,
                  center = attr(StormdataTrain7scale |> pull(Day), "scaled:center"),
                  scale = attr(StormdataTrain7scale |> pull(Day), "scaled:scale"))
    ))
)
save(gammaFit5D, file = "_data/gammaFit5D.RData")
prior_summary(gammaFit5D)
round(posterior_summary(gammaFit5D, probs = c(0.025, 0.975)))
gammaFit5D

print(gammaFit5D, digits = 4)
plot(gammaFit5D)
pp_check(gammaFit5D, ndraws = 100)
loo(gammaFit5D)
waic(gammaFit5D)
performance::check_distribution(gammaFit5D)
performance::check_outliers(gammaFit5D)
performance::check_heteroskedasticity(gammaFit5D)
performance_rmse(gammaFit5D)
performance_mae(gammaFit5D)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(gammaFit5D)


variance_decomposition(gammaFit5D)
exp(fixef(gammaFit5D))
ranef(gammaFit5D)

bayes_R2(gammaFit5D)

bayes_factor(gammaFit5D, gammaFit1)
bayes_factor(gammaFit5D, gammaFit2)
bayes_factor(gammaFit5D, gammaFit3)
bayes_factor(gammaFit5D, studentFit1)
bayes_factor(gammaFit5D, studentFit2)
bayes_factor(gammaFit5D, studentFit3)
bayes_factor(gammaFit5D, linFit11)
bayes_factor(gammaFit5D, propFit1)
bayes_factor(gammaFit5D, logPropFit1)
loo(gammaFit5D, gammaFit3)

gammaFit5Dsmooths <- conditional_smooths(gammaFit5D)
plot(gammaFit5Dsmooths, stype = "raster", ask = FALSE)
gammaFit5Deffects <- conditional_effects(gammaFit5D, 
                                         method = "posterior_predict",
                                         robust = FALSE,
                                         re_formula = NULL)
gammaFit5Deffects <- conditional_effects(gammaFit5D)
plot(gammaFit5Deffects, points = TRUE, ask = FALSE)

#### Prediction ----
## Fitted
gammaFit5DfinalFit <- posterior_predict(gammaFit5D)
gammaFit5DfinalFitMean <- colMeans(gammaFit5DfinalFit)
gammaFit5DfinalFitMed <- apply(gammaFit5DfinalFit, 2, function(x){quantile(x, 0.5)})
gammaFit5DfinalFitLCB <- apply(gammaFit5DfinalFit, 2, function(x){quantile(x, 0.025)})
gammaFit5DfinalFitUCB <- apply(gammaFit5DfinalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
gammaFit5DfinalPreds <- posterior_predict(gammaFit5D, 
                                          newdata = StormdataTestFinalscale,
                                          allow_new_levels = TRUE)
gammaFit5DfinalPreds2 <- colMeans(gammaFit5DfinalPreds)
gammaFit5DfinalPredsMed <- apply(gammaFit5DfinalPreds, 2, function(x){quantile(x, 0.5)})
gammaFit5DfinalPredsLCB <- apply(gammaFit5DfinalPreds, 2, function(x){quantile(x, 0.025)})
gammaFit5DfinalPredsUCB <- apply(gammaFit5DfinalPreds, 2, function(x){quantile(x, 0.975)})

gammaFit5DpredMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(gammaFit5DfinalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(gammaFit5DfinalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < gammaFit5DfinalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(gammaFit5DfinalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(gammaFit5DfinalPredsMed - Actual_Yvec)),
  COV_pred = mean(gammaFit5DfinalPredsLCB < Actual_Yvec & Actual_Yvec < gammaFit5DfinalPredsUCB)
)
gammaFit5DpredMetrics

#### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = gammaFit5DfinalPreds) +
  labs(title = "GammaFit5D Predict")

gammaFit5DFitDF <- bind_cols(
  StormdataTrain3,
  LCB = gammaFit5DfinalFitLCB,
  Mean = gammaFit5DfinalFitMean,
  Med = gammaFit5DfinalFitMed,
  UCB = gammaFit5DfinalFitUCB
) |>
  filter(StormID %in% c(5712015))

ggplot(data = gammaFit5DFitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))

## Prediction
gammaFit5DPredDF <- bind_cols(
  StormdataTest3,
  LCB = gammaFit5DfinalPredsLCB,
  Mean = gammaFit5DfinalPreds2,
  Med = gammaFit5DfinalPredsMed,
  UCB = gammaFit5DfinalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

ggplot(data = gammaFit5DPredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))

rm(gammaFit5DfinalFit)
rm(gammaFit5DfinalPreds)

### Model 6 ----
gammaFit6 <- brm(
  formula = VMAX ~
    #Year +
    #Month +
    #basin + 
    s(Day, bs = "cc") +
    #s(StormElapsedTime) + 
    #t2(LON, LAT) +
    t2(LON, LAT, StormElapsedTime, d = c(2,1)) +
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
  warmup = 1000,
  knots = list(
    Day = c(scale(0,
                  center = attr(StormdataTrain7scale |> pull(Day), "scaled:center"),
                  scale = attr(StormdataTrain7scale |> pull(Day), "scaled:scale")),
            scale(365,
                  center = attr(StormdataTrain7scale |> pull(Day), "scaled:center"),
                  scale = attr(StormdataTrain7scale |> pull(Day), "scaled:scale"))
    ))
  # LAT = c(scale(0,
  #               center = attr(StormdataTrain7scale |> pull(Day), "scaled:center"),
  #               scale = attr(StormdataTrain7scale |> pull(Day), "scaled:scale")),
  #         scale(365,
  #               center = attr(StormdataTrain7scale |> pull(Day), "scaled:center"),
  #               scale = attr(StormdataTrain7scale |> pull(Day), "scaled:scale"))
  # ),
  # LON = c(scale(0,
  #               center = attr(StormdataTrain7scale |> pull(Day), "scaled:center"),
  #               scale = attr(StormdataTrain7scale |> pull(Day), "scaled:scale")),
  #         scale(365,
  #               center = attr(StormdataTrain7scale |> pull(Day), "scaled:center"),
  #               scale = attr(StormdataTrain7scale |> pull(Day), "scaled:scale"))
  # ))
)
save(gammaFit6, file = "_data/gammaFit6.RData")
prior_summary(gammaFit6)
round(posterior_summary(gammaFit6, probs = c(0.025, 0.975)))
gammaFit6

print(gammaFit6, digits = 4)
plot(gammaFit6)
pp_check(gammaFit6, ndraws = 100)
loo(gammaFit6)
waic(gammaFit6)
performance::check_distribution(gammaFit6)
performance::check_outliers(gammaFit6)
performance::check_heteroskedasticity(gammaFit6)
performance_rmse(gammaFit6)
performance_mae(gammaFit6)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(gammaFit6)


variance_decomposition(gammaFit6)
exp(fixef(gammaFit6))
ranef(gammaFit6)

bayes_R2(gammaFit6)

bayes_factor(gammaFit6, gammaFit1)
bayes_factor(gammaFit6, gammaFit2)
bayes_factor(gammaFit6, gammaFit3)
bayes_factor(gammaFit6, studentFit1)
bayes_factor(gammaFit6, studentFit2)
bayes_factor(gammaFit6, studentFit3)
bayes_factor(gammaFit6, linFit11)
bayes_factor(gammaFit6, propFit1)
bayes_factor(gammaFit6, logPropFit1)
loo(gammaFit6, gammaFit3)

gammaFit6smooths <- conditional_smooths(gammaFit6)
plot(gammaFit6smooths, stype = "raster", ask = FALSE)
gammaFit6effects <- conditional_effects(gammaFit6, 
                                        method = "posterior_predict",
                                        robust = FALSE,
                                        re_formula = NULL)
gammaFit6effects <- conditional_effects(gammaFit6)
plot(gammaFit6effects, points = TRUE, ask = FALSE)

#### Prediction ----
## Fitted
gammaFit6finalFit <- posterior_predict(gammaFit6)
gammaFit6finalFitMean <- colMeans(gammaFit6finalFit)
gammaFit6finalFitMed <- apply(gammaFit6finalFit, 2, function(x){quantile(x, 0.5)})
gammaFit6finalFitLCB <- apply(gammaFit6finalFit, 2, function(x){quantile(x, 0.025)})
gammaFit6finalFitUCB <- apply(gammaFit6finalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
gammaFit6finalPreds <- posterior_predict(gammaFit6, 
                                         newdata = StormdataTestFinalscale,
                                         allow_new_levels = TRUE)
gammaFit6finalPreds2 <- colMeans(gammaFit6finalPreds)
gammaFit6finalPredsMed <- apply(gammaFit6finalPreds, 2, function(x){quantile(x, 0.5)})
gammaFit6finalPredsLCB <- apply(gammaFit6finalPreds, 2, function(x){quantile(x, 0.025)})
gammaFit6finalPredsUCB <- apply(gammaFit6finalPreds, 2, function(x){quantile(x, 0.975)})

gammaFit6predMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(gammaFit6finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(gammaFit6finalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < gammaFit6finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(gammaFit6finalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(gammaFit6finalPredsMed - Actual_Yvec)),
  COV_pred = mean(gammaFit6finalPredsLCB < Actual_Yvec & Actual_Yvec < gammaFit6finalPredsUCB)
)
gammaFit6predMetrics

ppc_dens_overlay(y = Actual_Yvec, yrep = gammaFit6finalPreds) +
  labs(title = "GammaFit6 Predict")
rm(gammaFit6finalFit)
rm(gammaFit6finalPreds)

### Model 7 ----
gammaFit7 <- brm(
  formula = VMAX ~
    #Year +
    #Month +
    #basin + 
    s(Day, bs = "cc") +
    #s(StormElapsedTime) + 
    #t2(LON, LAT) +
    t2(LON, LAT, StormElapsedTime, d = c(2,1)) +
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
  warmup = 1000,
  knots = list(
    Day = c(scale(0,
                  center = attr(StormdataTrain7scale |> pull(Day), "scaled:center"),
                  scale = attr(StormdataTrain7scale |> pull(Day), "scaled:scale")),
            scale(365,
                  center = attr(StormdataTrain7scale |> pull(Day), "scaled:center"),
                  scale = attr(StormdataTrain7scale |> pull(Day), "scaled:scale"))
    ))
  # LAT = c(scale(0,
  #               center = attr(StormdataTrain7scale |> pull(Day), "scaled:center"),
  #               scale = attr(StormdataTrain7scale |> pull(Day), "scaled:scale")),
  #         scale(365,
  #               center = attr(StormdataTrain7scale |> pull(Day), "scaled:center"),
  #               scale = attr(StormdataTrain7scale |> pull(Day), "scaled:scale"))
  # ),
  # LON = c(scale(0,
  #               center = attr(StormdataTrain7scale |> pull(Day), "scaled:center"),
  #               scale = attr(StormdataTrain7scale |> pull(Day), "scaled:scale")),
  #         scale(365,
  #               center = attr(StormdataTrain7scale |> pull(Day), "scaled:center"),
  #               scale = attr(StormdataTrain7scale |> pull(Day), "scaled:scale"))
  # ))
)
save(gammaFit7, file = "_data/gammaFit7.RData")
prior_summary(gammaFit7)
round(posterior_summary(gammaFit7, probs = c(0.025, 0.975)))
gammaFit7

print(gammaFit7, digits = 4)
plot(gammaFit7)
pp_check(gammaFit7, ndraws = 100)
loo(gammaFit7)
waic(gammaFit7)
performance::check_distribution(gammaFit7)
performance::check_outliers(gammaFit7)
performance::check_heteroskedasticity(gammaFit7)
performance_rmse(gammaFit7)
performance_mae(gammaFit7)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(gammaFit7)


variance_decomposition(gammaFit7)
exp(fixef(gammaFit7))
ranef(gammaFit7)

bayes_R2(gammaFit7)

bayes_factor(gammaFit7, gammaFit1)
bayes_factor(gammaFit7, gammaFit2)
bayes_factor(gammaFit7, gammaFit3)
bayes_factor(gammaFit7, studentFit1)
bayes_factor(gammaFit7, studentFit2)
bayes_factor(gammaFit7, studentFit3)
bayes_factor(gammaFit7, linFit11)
bayes_factor(gammaFit7, propFit1)
bayes_factor(gammaFit7, logPropFit1)
loo(gammaFit7, gammaFit3)

gammaFit7smooths <- conditional_smooths(gammaFit7)
plot(gammaFit7smooths, stype = "raster", ask = FALSE)
gammaFit7effects <- conditional_effects(gammaFit7, 
                                        method = "posterior_predict",
                                        robust = FALSE,
                                        re_formula = NULL)
gammaFit7effects <- conditional_effects(gammaFit7)
plot(gammaFit7effects, points = TRUE)

#### Prediction ----
## Fitted
gammaFit7finalFit <- posterior_predict(gammaFit7)
gammaFit7finalFitMean <- colMeans(gammaFit7finalFit)
gammaFit7finalFitMed <- apply(gammaFit7finalFit, 2, function(x){quantile(x, 0.5)})
gammaFit7finalFitLCB <- apply(gammaFit7finalFit, 2, function(x){quantile(x, 0.025)})
gammaFit7finalFitUCB <- apply(gammaFit7finalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
gammaFit7finalPreds <- posterior_predict(gammaFit7, 
                                         newdata = StormdataTestFinalscale,
                                         allow_new_levels = TRUE)
gammaFit7finalPreds2 <- colMeans(gammaFit7finalPreds)
gammaFit7finalPredsMed <- apply(gammaFit7finalPreds, 2, function(x){quantile(x, 0.5)})
gammaFit7finalPredsLCB <- apply(gammaFit7finalPreds, 2, function(x){quantile(x, 0.025)})
gammaFit7finalPredsUCB <- apply(gammaFit7finalPreds, 2, function(x){quantile(x, 0.975)})

gammaFit7predMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(gammaFit7finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(gammaFit7finalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < gammaFit7finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(gammaFit7finalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(gammaFit7finalPredsMed - Actual_Yvec)),
  COV_pred = mean(gammaFit7finalPredsLCB < Actual_Yvec & Actual_Yvec < gammaFit7finalPredsUCB)
)
gammaFit7predMetrics

ppc_dens_overlay(y = Actual_Yvec, yrep = gammaFit7finalPreds) +
  labs(title = "GammaFit7 Predict")
rm(gammaFit7finalFit)
rm(gammaFit7finalPreds)

### Model 5 ----
gammaFit8 <- brm(
  formula = VMAX ~
    #Year +
    #Month +
    basin + 
    s(Day) +
    t2(LON, LAT, StormElapsedTime, d = c(2,1), bs = c("tp", "cr")) +
    t2(LON, LAT, StormElapsedTime, StormID, d = c(2,1,1), bs = c("tp","cr", "re")) +
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

## MGCV ----
### Model 1 ----
set.seed(52)
mgcvFit1 <- gam(
  formula = VMAX ~
    #Year +
    Month +
    basin + 
    #t2(LON, LAT, StormElapsedTime, bs = c("tp","tp"), k = c(20,10)) +
    #s(StormElapsedTime, bs = "cr") +
    t2(LON, LAT, StormElapsedTime, bs = c("cr", "cr", "cr")) +
    #t2(LON, LAT, StormElapsedTime, StormID, bs = c("cr", "cr", "cr","re")) +
    #t2(LON, LAT, StormElapsedTime, StormID) +
    #bs = c("cr", "cr", "cr"), k = c(20, 10, 5)) +
    # t2(LON, LAT, StormElapsedTime, StormID,
    #    d = c(2,1,1), bs = c("tp", "cr", "re")) +
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
    HWRF,# +
  #s(StormID, bs = "re"),
  data = StormdataTrain7, 
  family = Gamma(link = "log"),
  method = "REML"
)

plot(gammaFit5, pages = 1, residuals = TRUE)
summary(gammaFit5)
anova.gam(gammaFit5)
anova(gammaFit8C, gammaFit8C1, test = "Chisq")
gam.check(gammaFit5)
k.check(gammaFit5)
performance_mae(gammaFit5)
gammaFit5Preds <- predict(gammaFit5, newdata = StormdataTestFinal, type = "response")
mean(abs(gammaFit5Preds - Actual_Yvec))


## GAUSSIAN PROP ----
### Model 1 ----
propLinFit1 <- brm(
  formula = propVMAX ~
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
    (1|StormID),
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
pp_check(propLinFit1, ndraws = 100)
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

bayes_factor(propLinFit1, gammaFit1)
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

#### Prediction ----
## Fitted
propLinFit1finalFit <- posterior_predict(propLinFit1)
propLinFit1finalFit <- t(t(propLinFit1finalFit)*StormdataTrain3$HWRF)
#propLinFit1finalFitMean <- colMeans(propLinFit1finalFit)*StormdataTrain3$HWRF
propLinFit1finalFitMean <- colMeans(propLinFit1finalFit)
propLinFit1finalFitMed <- apply(propLinFit1finalFit, 2, function(x){quantile(x, 0.5)})
propLinFit1finalFitLCB <- apply(propLinFit1finalFit, 2, function(x){quantile(x, 0.025)})
propLinFit1finalFitUCB <- apply(propLinFit1finalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
propLinFit1finalPreds <- posterior_predict(propLinFit1, 
                                           newdata = StormdataTestFinalscale,
                                           allow_new_levels = TRUE)
propLinFit1finalPreds <- t(t(propLinFit1finalPreds)*StormdataTest3$HWRF)
propLinFit1finalPreds2 <- colMeans(propLinFit1finalPreds)
propLinFit1finalPredsMed <- apply(propLinFit1finalPreds, 2, function(x){quantile(x, 0.5)})
propLinFit1finalPredsLCB <- apply(propLinFit1finalPreds, 2, function(x){quantile(x, 0.025)})
propLinFit1finalPredsUCB <- apply(propLinFit1finalPreds, 2, function(x){quantile(x, 0.975)})

propLinFit1predMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(propLinFit1finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(propLinFit1finalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < propLinFit1finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(propLinFit1finalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(propLinFit1finalPredsMed - Actual_Yvec)),
  COV_pred = mean(propLinFit1finalPredsLCB < Actual_Yvec & Actual_Yvec < propLinFit1finalPredsUCB)
)
propLinFit1predMetrics

#### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = propLinFit1finalPreds) +
  labs(title = "GammaFit5 Predict")

propLinFit1FitDF <- bind_cols(
  StormdataTrain3,
  LCB = propLinFit1finalFitLCB,
  Mean = propLinFit1finalFitMean,
  Med = propLinFit1finalFitMed,
  UCB = propLinFit1finalFitUCB
) 

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
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))

rm(propLinFit1finalFit)
rm(propLinFit1finalPreds)

### Model 2 ----
propLinFit2 <- brm(
  bf(
    VMAX ~ HWRF*exp(eta),
    eta ~
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
      (1|StormID),
    nl = TRUE
  ),
  data = StormdataTrain8, 
  family = gaussian(), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)
save(propLinFit2, file = "_data/propLinFit2.RData")
prior_summary(propLinFit2)
round(posterior_summary(propLinFit2, probs = c(0.025, 0.975)))
propLinFit2

print(propLinFit2, digits = 4)
plot(propLinFit2)
pp_check(propLinFit2, ndraws = 100)
loo(propLinFit2)
waic(propLinFit2)
performance::check_distribution(propLinFit2)
performance::check_outliers(propLinFit2)
performance::check_heteroskedasticity(propLinFit2)
performance_rmse(propLinFit2)
performance_mae(propLinFit2)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(propLinFit2)


variance_decomposition(propLinFit2)
exp(fixef(propLinFit2))
ranef(propLinFit2)

bayes_R2(propLinFit2)

bayes_factor(propLinFit2, propLinFit1)
bayes_factor(propLinFit2, gammaFit2)
bayes_factor(propLinFit2, gammaFit3)
bayes_factor(propLinFit2, studentFit2)
bayes_factor(propLinFit2, studentFit2)
bayes_factor(propLinFit2, studentFit3)
bayes_factor(propLinFit2, linFit21)
bayes_factor(propLinFit2, propFit2)
bayes_factor(propLinFit2, logPropFit2)
loo(propLinFit2, gammaFit3)

propLinFit2smooths <- conditional_smooths(propLinFit2)
plot(propLinFit2smooths, stype = "raster", ask = FALSE)
propLinFit2effects <- conditional_effects(propLinFit2, 
                                          method = "posterior_predict",
                                          robust = FALSE,
                                          re_formula = NULL)
propLinFit2effects <- conditional_effects(propLinFit2)
plot(propLinFit2effects, points = TRUE, ask = FALSE)

#### Prediction ----
## Fitted
propLinFit2finalFit <- posterior_predict(propLinFit2)
propLinFit2finalFit <- t(t(propLinFit2finalFit)*StormdataTrain3$HWRF)
#propLinFit2finalFitMean <- colMeans(propLinFit2finalFit)*StormdataTrain3$HWRF
propLinFit2finalFitMean <- colMeans(propLinFit2finalFit)
propLinFit2finalFitMed <- apply(propLinFit2finalFit, 2, function(x){quantile(x, 0.5)})
propLinFit2finalFitLCB <- apply(propLinFit2finalFit, 2, function(x){quantile(x, 0.025)})
propLinFit2finalFitUCB <- apply(propLinFit2finalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
propLinFit2finalPreds <- posterior_predict(propLinFit2, 
                                           newdata = StormdataTestFinalscale2,
                                           allow_new_levels = TRUE)
propLinFit2finalPreds <- t(t(propLinFit2finalPreds)*StormdataTest3$HWRF)
propLinFit2finalPreds2 <- colMeans(propLinFit2finalPreds)
propLinFit2finalPredsMed <- apply(propLinFit2finalPreds, 2, function(x){quantile(x, 0.5)})
propLinFit2finalPredsLCB <- apply(propLinFit2finalPreds, 2, function(x){quantile(x, 0.025)})
propLinFit2finalPredsUCB <- apply(propLinFit2finalPreds, 2, function(x){quantile(x, 0.975)})

propLinFit2predMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(propLinFit2finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(propLinFit2finalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < propLinFit2finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(propLinFit2finalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(propLinFit2finalPredsMed - Actual_Yvec)),
  COV_pred = mean(propLinFit2finalPredsLCB < Actual_Yvec & Actual_Yvec < propLinFit2finalPredsUCB)
)
propLinFit2predMetrics

#### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = propLinFit2finalPreds) +
  labs(title = "GammaFit5 Predict")

propLinFit2FitDF <- bind_cols(
  StormdataTrain3,
  LCB = propLinFit2finalFitLCB,
  Mean = propLinFit2finalFitMean,
  Med = propLinFit2finalFitMed,
  UCB = propLinFit2finalFitUCB
) 

ggplot(data = propLinFit2FitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))

## Prediction
propLinFit2PredDF <- bind_cols(
  StormdataTest3,
  LCB = propLinFit2finalPredsLCB,
  Mean = propLinFit2finalPreds2,
  Med = propLinFit2finalPredsMed,
  UCB = propLinFit2finalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

ggplot(data = propLinFit2PredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))

rm(propLinFit2finalFit)
rm(propLinFit2finalPreds)

### Model 2B ----
propLinFit2B <- brm(
  bf(
    VMAX ~ HWRF*exp(eta),
    eta ~
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
      (1|StormID),
    nl = TRUE
  ),
  data = StormdataTrain8, 
  family = lognormal(link_sigma = "log"), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)
save(propLinFit2B, file = "_data/propLinFit2B.RData")
prior_summary(propLinFit2B)
round(posterior_summary(propLinFit2B, probs = c(0.025, 0.975)))
propLinFit2B

print(propLinFit2B, digits = 4)
plot(propLinFit2B)
pp_check(propLinFit2B, ndraws = 100)
loo(propLinFit2B)
waic(propLinFit2B)
performance::check_distribution(propLinFit2B)
performance::check_outliers(propLinFit2B)
performance::check_heteroskedasticity(propLinFit2B)
performance_rmse(propLinFit2B)
performance_mae(propLinFit2B)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(propLinFit2B)


variance_decomposition(propLinFit2B)
exp(fixef(propLinFit2B))
ranef(propLinFit2B)

bayes_R2(propLinFit2B)

bayes_factor(propLinFit2B, propLinFit1)
bayes_factor(propLinFit2B, gammaFit2B)
bayes_factor(propLinFit2B, gammaFit3)
bayes_factor(propLinFit2B, studentFit2B)
bayes_factor(propLinFit2B, studentFit2B)
bayes_factor(propLinFit2B, studentFit3)
bayes_factor(propLinFit2B, linFit2B1)
bayes_factor(propLinFit2B, propFit2B)
bayes_factor(propLinFit2B, logPropFit2B)
loo(propLinFit2B, gammaFit3)

propLinFit2Bsmooths <- conditional_smooths(propLinFit2B)
plot(propLinFit2Bsmooths, stype = "raster", ask = FALSE)
propLinFit2Beffects <- conditional_effects(propLinFit2B, 
                                           method = "posterior_predict",
                                           robust = FALSE,
                                           re_formula = NULL)
propLinFit2Beffects <- conditional_effects(propLinFit2B)
plot(propLinFit2Beffects, points = TRUE, ask = FALSE)

#### Prediction ----
## Fitted
propLinFit2BfinalFit <- posterior_predict(propLinFit2B)
propLinFit2BfinalFit <- t(t(propLinFit2BfinalFit)*StormdataTrain3$HWRF)
#propLinFit2BfinalFitMean <- colMeans(propLinFit2BfinalFit)*StormdataTrain3$HWRF
propLinFit2BfinalFitMean <- colMeans(propLinFit2BfinalFit)
propLinFit2BfinalFitMed <- apply(propLinFit2BfinalFit, 2, function(x){quantile(x, 0.5)})
propLinFit2BfinalFitLCB <- apply(propLinFit2BfinalFit, 2, function(x){quantile(x, 0.025)})
propLinFit2BfinalFitUCB <- apply(propLinFit2BfinalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
propLinFit2BfinalPreds <- posterior_predict(propLinFit2B, 
                                            newdata = StormdataTestFinalscale2,
                                            allow_new_levels = TRUE)
propLinFit2BfinalPreds <- t(t(propLinFit2BfinalPreds)*StormdataTest3$HWRF)
propLinFit2BfinalPreds2 <- colMeans(propLinFit2BfinalPreds)
propLinFit2BfinalPredsMed <- apply(propLinFit2BfinalPreds, 2, function(x){quantile(x, 0.5)})
propLinFit2BfinalPredsLCB <- apply(propLinFit2BfinalPreds, 2, function(x){quantile(x, 0.025)})
propLinFit2BfinalPredsUCB <- apply(propLinFit2BfinalPreds, 2, function(x){quantile(x, 0.975)})

propLinFit2BpredMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(propLinFit2BfinalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(propLinFit2BfinalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < propLinFit2BfinalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(propLinFit2BfinalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(propLinFit2BfinalPredsMed - Actual_Yvec)),
  COV_pred = mean(propLinFit2BfinalPredsLCB < Actual_Yvec & Actual_Yvec < propLinFit2BfinalPredsUCB)
)
propLinFit2BpredMetrics

#### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = propLinFit2BfinalPreds) +
  labs(title = "GammaFit5 Predict")

propLinFit2BFitDF <- bind_cols(
  StormdataTrain3,
  LCB = propLinFit2BfinalFitLCB,
  Mean = propLinFit2BfinalFitMean,
  Med = propLinFit2BfinalFitMed,
  UCB = propLinFit2BfinalFitUCB
) 

ggplot(data = propLinFit2BFitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))

## Prediction
propLinFit2BPredDF <- bind_cols(
  StormdataTest3,
  LCB = propLinFit2BfinalPredsLCB,
  Mean = propLinFit2BfinalPreds2,
  Med = propLinFit2BfinalPredsMed,
  UCB = propLinFit2BfinalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

ggplot(data = propLinFit2BPredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))

rm(propLinFit2BfinalFit)
rm(propLinFit2BfinalPreds)

### Model 3 ----
propLinFit3 <- brm(
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
save(propLinFit3, file = "_data/propLinFit3.RData")
prior_summary(propLinFit3)
round(posterior_summary(propLinFit3, probs = c(0.025, 0.975)))
propLinFit3

print(propLinFit3, digits = 4)
plot(propLinFit3)
pp_check(propLinFit3, ndraws = 100) 
loo(propLinFit3)
waic(propLinFit3)
performance::check_distribution(propLinFit3)
performance::check_outliers(propLinFit3)
performance::check_heteroskedasticity(propLinFit3)
performance_rmse(propLinFit3)
performance_mae(propLinFit3)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(propLinFit3)


variance_decomposition(propLinFit3)
exp(fixef(propLinFit3))
ranef(propLinFit3)

bayes_R2(propLinFit3)

bayes_factor(propLinFit3, propLinFit1)
bayes_factor(propLinFit3, gammaFit3)
bayes_factor(propLinFit3, gammaFit3)
bayes_factor(propLinFit3, studentFit3)
bayes_factor(propLinFit3, studentFit3)
bayes_factor(propLinFit3, studentFit3)
bayes_factor(propLinFit3, linFit31)
bayes_factor(propLinFit3, propFit3)
bayes_factor(propLinFit3, logPropFit3)
loo(propLinFit3, gammaFit3)

propLinFit3smooths <- conditional_smooths(propLinFit3)
plot(propLinFit3smooths, stype = "raster", ask = FALSE)
propLinFit3effects <- conditional_effects(propLinFit3, 
                                          method = "posterior_predict",
                                          robust = FALSE,
                                          re_formula = NULL)
propLinFit3effects <- conditional_effects(propLinFit3)
plot(propLinFit3effects, points = TRUE, ask = FALSE)

#### Prediction ----
## Fitted
propLinFit3finalFit <- posterior_predict(propLinFit3)
propLinFit3finalFit2 <- t(t(propLinFit3finalFit)*StormdataTrain3$HWRF)
#propLinFit3finalFitMean <- colMeans(propLinFit3finalFit)*StormdataTrain3$HWRF
propLinFit3finalFitMean <- colMeans(propLinFit3finalFit2)
propLinFit3finalFitMed <- apply(propLinFit3finalFit2, 2, function(x){quantile(x, 0.5)})
propLinFit3finalFitLCB <- apply(propLinFit3finalFit2, 2, function(x){quantile(x, 0.025)})
propLinFit3finalFitUCB <- apply(propLinFit3finalFit2, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
propLinFit3finalPreds <- posterior_predict(propLinFit3, 
                                           newdata = StormdataTestFinalscale2,
                                           allow_new_levels = TRUE)
propLinFit3finalPreds2 <- t(t(propLinFit3finalPreds)*StormdataTest3$HWRF)
propLinFit3finalPreds3 <- colMeans(propLinFit3finalPreds2)
propLinFit3finalPredsMed <- apply(propLinFit3finalPreds2, 2, function(x){quantile(x, 0.5)})
propLinFit3finalPredsLCB <- apply(propLinFit3finalPreds2, 2, function(x){quantile(x, 0.025)})
propLinFit3finalPredsUCB <- apply(propLinFit3finalPreds2, 2, function(x){quantile(x, 0.975)})

propLinFit3predMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(propLinFit3finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(propLinFit3finalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < propLinFit3finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(propLinFit3finalPreds3 - Actual_Yvec)),
  MAD_pred = mean(abs(propLinFit3finalPredsMed - Actual_Yvec)),
  COV_pred = mean(propLinFit3finalPredsLCB < Actual_Yvec & Actual_Yvec < propLinFit3finalPredsUCB)
)
propLinFit3predMetrics

#### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = propLinFit3finalPreds) +
  labs(title = "propFit3 Predict")

propLinFit3FitDF <- bind_cols(
  StormdataTrain3,
  LCB = propLinFit3finalFitLCB,
  Mean = propLinFit3finalFitMean,
  Med = propLinFit3finalFitMed,
  UCB = propLinFit3finalFitUCB
) 

propLinFit3ppcDraws <- as_draws_matrix(propLinFit3)

propLinFit3basePPCplot <- ggplot(data = propLinFit3FitDF) +
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
propLinFit3basePPCplot

set.seed(52)
propLinFit3PPCindex <- sample(1:1705, 100)
for(i in propLinFit3PPCindex){
  propLinFit3basePPCplot <- propLinFit3basePPCplot +
    geom_density(
      aes(x = propLinFit3finalFit[i,])
    )
}
propLinFit3basePPCplot

ggplot(data = propLinFit3FitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))

## Prediction
propLinFit3PredDF <- bind_cols(
  StormdataTest3,
  LCB = propLinFit3finalPredsLCB,
  Mean = propLinFit3finalPreds2,
  Med = propLinFit3finalPredsMed,
  UCB = propLinFit3finalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

ggplot(data = propLinFit3PredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))

rm(propLinFit3finalFit)
rm(propLinFit3finalFit2)
rm(propLinFit3finalPreds)
rm(propLinFit3finalPreds2)

### Model 3B ----
propLinFit3B <- brm(
  bf(
    VMAX/HWRF ~
      #Year +
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

#### Prediction ----
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

#### Plotting ----
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

### Model 3C ----
propLinFit3C <- brm(
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
save(propLinFit3C, file = "_data/propLinFit3C.RData")
prior_summary(propLinFit3C)
round(posterior_summary(propLinFit3C, probs = c(0.025, 0.975)))
propLinFit3C

print(propLinFit3C, digits = 4)
plot(propLinFit3C)
pp_check(propLinFit3C, ndraws = 100)
loo(propLinFit3C)
waic(propLinFit3C)
performance::check_distribution(propLinFit3C)
performance::check_outliers(propLinFit3C)
performance::check_heteroskedasticity(propLinFit3C)
performance_rmse(propLinFit3C)
performance_mae(propLinFit3C)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(propLinFit3C)


variance_decomposition(propLinFit3C)
exp(fixef(propLinFit3C))
ranef(propLinFit3C)

bayes_R2(propLinFit3C)
loo_R2(propLinFit3C)

bayes_factor(propLinFit3C, propLinFit1)
bayes_factor(propLinFit3C, gammaFit3C)
bayes_factor(propLinFit3C, gammaFit3C)
bayes_factor(propLinFit3C, studentFit3C)
bayes_factor(propLinFit3C, studentFit3C)
bayes_factor(propLinFit3C, studentFit3C)
bayes_factor(propLinFit3C, linFit3C1)
bayes_factor(propLinFit3C, propFit3C)
bayes_factor(propLinFit3C, logPropFit3C)
loo(propLinFit3C, gammaFit3C)

propLinFit3Csmooths <- conditional_smooths(propLinFit3C)
plot(propLinFit3Csmooths, stype = "raster", ask = FALSE)
propLinFit3Ceffects <- conditional_effects(propLinFit3C, 
                                           method = "posterior_predict",
                                           robust = FALSE,
                                           re_formula = NULL)
propLinFit3Ceffects <- conditional_effects(propLinFit3C)
plot(propLinFit3Ceffects, points = TRUE, ask = FALSE)

#### Prediction ----
## Fitted
propLinFit3CfinalFit <- posterior_predict(propLinFit3C)
propLinFit3CfinalFit2 <- t(t(propLinFit3CfinalFit)*StormdataTrain3$HWRF)
#propLinFit3CfinalFitMean <- colMeans(propLinFit3CfinalFit)*StormdataTrain3$HWRF
propLinFit3CfinalFitMean <- colMeans(propLinFit3CfinalFit2)
propLinFit3CfinalFitMed <- apply(propLinFit3CfinalFit2, 2, function(x){quantile(x, 0.5)})
propLinFit3CfinalFitLCB <- apply(propLinFit3CfinalFit2, 2, function(x){quantile(x, 0.025)})
propLinFit3CfinalFitUCB <- apply(propLinFit3CfinalFit2, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
propLinFit3CfinalPreds <- posterior_predict(propLinFit3C, 
                                            newdata = StormdataTestFinalscale2,
                                            allow_new_levels = TRUE)
propLinFit3CfinalPreds <- t(t(propLinFit3CfinalPreds)*StormdataTest3$HWRF)
propLinFit3CfinalPreds2 <- colMeans(propLinFit3CfinalPreds)
propLinFit3CfinalPredsMed <- apply(propLinFit3CfinalPreds, 2, function(x){quantile(x, 0.5)})
propLinFit3CfinalPredsLCB <- apply(propLinFit3CfinalPreds, 2, function(x){quantile(x, 0.025)})
propLinFit3CfinalPredsUCB <- apply(propLinFit3CfinalPreds, 2, function(x){quantile(x, 0.975)})

propLinFit3CpredMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(propLinFit3CfinalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(propLinFit3CfinalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < propLinFit3CfinalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(propLinFit3CfinalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(propLinFit3CfinalPredsMed - Actual_Yvec)),
  COV_pred = mean(propLinFit3CfinalPredsLCB < Actual_Yvec & Actual_Yvec < propLinFit3CfinalPredsUCB)
)
propLinFit3CpredMetrics

#### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = propLinFit3CfinalPreds) +
  labs(title = "GammaFit5 Predict")

propLinFit3CFitDF <- bind_cols(
  StormdataTrain3,
  LCB = propLinFit3CfinalFitLCB,
  Mean = propLinFit3CfinalFitMean,
  Med = propLinFit3CfinalFitMed,
  UCB = propLinFit3CfinalFitUCB
) 

ggplot(data = propLinFit3CFitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))

## Prediction
propLinFit3CPredDF <- bind_cols(
  StormdataTest3,
  LCB = propLinFit3CfinalPredsLCB,
  Mean = propLinFit3CfinalPreds2,
  Med = propLinFit3CfinalPredsMed,
  UCB = propLinFit3CfinalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

ggplot(data = propLinFit3CPredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))

rm(propLinFit3CfinalFit)
rm(propLinFit3CfinalFit2)
rm(propLinFit3CfinalPreds)

### Model 4 ----
propLinFit4prior <- prior(constant(1), class = "b", coef = "HWRF")
propLinFit4 <- brm(
  bf(
    log(VMAX) ~ log(HWRF) + eta,
    eta ~
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
      (1|StormID),
    nl = TRUE
  ),
  data = StormdataTrain8, 
  family = gaussian(), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
  #prior = prior(constant(1), nlpar = "HWRF")
)
save(propLinFit4, file = "_data/propLinFit4.RData")
prior_summary(propLinFit4)
round(posterior_summary(propLinFit4, probs = c(0.025, 0.975)))
propLinFit4

print(propLinFit4, digits = 4)
plot(propLinFit4)
pp_check(propLinFit4, ndraws = 100)
loo(propLinFit4)
waic(propLinFit4)
performance::check_distribution(propLinFit4)
performance::check_outliers(propLinFit4)
performance::check_heteroskedasticity(propLinFit4)
performance_rmse(propLinFit4)
performance_mae(propLinFit4)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(propLinFit4)


variance_decomposition(propLinFit4)
exp(fixef(propLinFit4))
ranef(propLinFit4)

bayes_R2(propLinFit4)

bayes_factor(propLinFit4, propLinFit1)
bayes_factor(propLinFit4, gammaFit4)
bayes_factor(propLinFit4, gammaFit4)
bayes_factor(propLinFit4, studentFit4)
bayes_factor(propLinFit4, studentFit4)
bayes_factor(propLinFit4, studentFit4)
bayes_factor(propLinFit4, linFit41)
bayes_factor(propLinFit4, propFit4)
bayes_factor(propLinFit4, logPropFit4)
loo(propLinFit4, gammaFit4)

propLinFit4smooths <- conditional_smooths(propLinFit4)
plot(propLinFit4smooths, stype = "raster", ask = FALSE)
propLinFit4effects <- conditional_effects(propLinFit4, 
                                          method = "posterior_predict",
                                          robust = FALSE,
                                          re_formula = NULL)
propLinFit4effects <- conditional_effects(propLinFit4)
plot(propLinFit4effects, points = TRUE, ask = FALSE)

#### Prediction ----
## Fitted
propLinFit4finalFit <- posterior_predict(propLinFit4)
propLinFit4finalFit <- t(t(propLinFit4finalFit)*StormdataTrain3$HWRF)
#propLinFit4finalFitMean <- colMeans(propLinFit4finalFit)*StormdataTrain3$HWRF
propLinFit4finalFitMean <- colMeans(propLinFit4finalFit)
propLinFit4finalFitMed <- apply(propLinFit4finalFit, 2, function(x){quantile(x, 0.5)})
propLinFit4finalFitLCB <- apply(propLinFit4finalFit, 2, function(x){quantile(x, 0.025)})
propLinFit4finalFitUCB <- apply(propLinFit4finalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
propLinFit4finalPreds <- posterior_predict(propLinFit4, 
                                           newdata = StormdataTestFinalscale2,
                                           allow_new_levels = TRUE)
propLinFit4finalPreds <- t(t(propLinFit4finalPreds)*StormdataTest3$HWRF)
propLinFit4finalPreds2 <- colMeans(propLinFit4finalPreds)
propLinFit4finalPredsMed <- apply(propLinFit4finalPreds, 2, function(x){quantile(x, 0.5)})
propLinFit4finalPredsLCB <- apply(propLinFit4finalPreds, 2, function(x){quantile(x, 0.025)})
propLinFit4finalPredsUCB <- apply(propLinFit4finalPreds, 2, function(x){quantile(x, 0.975)})

propLinFit4predMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(propLinFit4finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(propLinFit4finalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < propLinFit4finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(propLinFit4finalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(propLinFit4finalPredsMed - Actual_Yvec)),
  COV_pred = mean(propLinFit4finalPredsLCB < Actual_Yvec & Actual_Yvec < propLinFit4finalPredsUCB)
)
propLinFit4predMetrics

#### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = propLinFit4finalPreds) +
  labs(title = "GammaFit5 Predict")

propLinFit4FitDF <- bind_cols(
  StormdataTrain3,
  LCB = propLinFit4finalFitLCB,
  Mean = propLinFit4finalFitMean,
  Med = propLinFit4finalFitMed,
  UCB = propLinFit4finalFitUCB
) 

ggplot(data = propLinFit4FitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))

## Prediction
propLinFit4PredDF <- bind_cols(
  StormdataTest3,
  LCB = propLinFit4finalPredsLCB,
  Mean = propLinFit4finalPreds2,
  Med = propLinFit4finalPredsMed,
  UCB = propLinFit4finalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

ggplot(data = propLinFit4PredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))

rm(propLinFit4finalFit)
rm(propLinFit4finalPreds)

### Model 5 ----
propLinFit5 <- brm(
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
save(propLinFit5, file = "_data/propLinFit5.RData")
prior_summary(propLinFit5)
round(posterior_summary(propLinFit5, probs = c(0.025, 0.975)))
propLinFit5

print(propLinFit5, digits = 4)
plot(propLinFit5)
pp_check(propLinFit5, ndraws = 100)
loo(propLinFit5)
waic(propLinFit5)
performance::check_distribution(propLinFit5)
performance::check_outliers(propLinFit5)
performance::check_heteroskedasticity(propLinFit5)
performance_rmse(propLinFit5)
performance_mae(propLinFit5)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(propLinFit5)


variance_decomposition(propLinFit5)
exp(fixef(propLinFit5))
ranef(propLinFit5)

bayes_R2(propLinFit5)

bayes_factor(propLinFit5, propLinFit1)
bayes_factor(propLinFit5, gammaFit5)
bayes_factor(propLinFit5, gammaFit5)
bayes_factor(propLinFit5, studentFit5)
bayes_factor(propLinFit5, studentFit5)
bayes_factor(propLinFit5, studentFit5)
bayes_factor(propLinFit5, linFit51)
bayes_factor(propLinFit5, propFit5)
bayes_factor(propLinFit5, logPropFit5)
loo(propLinFit5, gammaFit5)

propLinFit5smooths <- conditional_smooths(propLinFit5)
plot(propLinFit5smooths, stype = "raster", ask = FALSE)
propLinFit5effects <- conditional_effects(propLinFit5, 
                                          method = "posterior_predict",
                                          robust = FALSE,
                                          re_formula = NULL)
propLinFit5effects <- conditional_effects(propLinFit5)
plot(propLinFit5effects, points = TRUE, ask = FALSE)

#### Prediction ----
## Fitted
propLinFit5finalFit <- posterior_predict(propLinFit5)
propLinFit5finalFit <- t(t(exp(propLinFit5finalFit))*StormdataTrain3$HWRF)
#propLinFit5finalFitMean <- colMeans(propLinFit5finalFit)*StormdataTrain3$HWRF
propLinFit5finalFitMean <- colMeans(propLinFit5finalFit)
propLinFit5finalFitMed <- apply(propLinFit5finalFit, 2, function(x){quantile(x, 0.5)})
propLinFit5finalFitLCB <- apply(propLinFit5finalFit, 2, function(x){quantile(x, 0.025)})
propLinFit5finalFitUCB <- apply(propLinFit5finalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
propLinFit5finalPreds <- posterior_predict(propLinFit5, 
                                           newdata = StormdataTestFinalscale2,
                                           allow_new_levels = TRUE)
propLinFit5finalPreds <- t(t(exp(propLinFit5finalPreds))*StormdataTest3$HWRF)
propLinFit5finalPreds2 <- colMeans(propLinFit5finalPreds)
propLinFit5finalPredsMed <- apply(propLinFit5finalPreds, 2, function(x){quantile(x, 0.5)})
propLinFit5finalPredsLCB <- apply(propLinFit5finalPreds, 2, function(x){quantile(x, 0.025)})
propLinFit5finalPredsUCB <- apply(propLinFit5finalPreds, 2, function(x){quantile(x, 0.975)})

propLinFit5predMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(propLinFit5finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(propLinFit5finalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < propLinFit5finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(propLinFit5finalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(propLinFit5finalPredsMed - Actual_Yvec)),
  COV_pred = mean(propLinFit5finalPredsLCB < Actual_Yvec & Actual_Yvec < propLinFit5finalPredsUCB)
)
propLinFit5predMetrics

#### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = propLinFit5finalPreds) +
  labs(title = "GammaFit5 Predict")

propLinFit5FitDF <- bind_cols(
  StormdataTrain3,
  LCB = propLinFit5finalFitLCB,
  Mean = propLinFit5finalFitMean,
  Med = propLinFit5finalFitMed,
  UCB = propLinFit5finalFitUCB
) 

ggplot(data = propLinFit5FitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))

## Prediction
propLinFit5PredDF <- bind_cols(
  StormdataTest3,
  LCB = propLinFit5finalPredsLCB,
  Mean = propLinFit5finalPreds2,
  Med = propLinFit5finalPredsMed,
  UCB = propLinFit5finalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

ggplot(data = propLinFit5PredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))

rm(propLinFit5finalFit)
rm(propLinFit5finalPreds)

### Model 6 ----
propLinFit6 <- brm(
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
save(propLinFit6, file = "_data/propLinFit6.RData")
prior_summary(propLinFit6)
round(posterior_summary(propLinFit6, probs = c(0.025, 0.975)))
propLinFit6

print(propLinFit6, digits = 4)
pp_check(propLinFit6, ndraws = 1705)
plot(propLinFit6)
loo(propLinFit6)
waic(propLinFit6)
performance::check_distribution(propLinFit6)
performance::check_outliers(propLinFit6)
performance::check_heteroskedasticity(propLinFit6)
performance_rmse(propLinFit6)
performance_mae(propLinFit6)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(propLinFit6)


variance_decomposition(propLinFit6)
exp(fixef(propLinFit6))
ranef(propLinFit6)

bayes_R2(propLinFit6)
loo_R2(propLinFit6)

bayes_factor(propLinFit6, propLinFit1)
bayes_factor(propLinFit6, gammaFit6)
bayes_factor(propLinFit6, gammaFit6)
bayes_factor(propLinFit6, studentFit6)
bayes_factor(propLinFit6, studentFit6)
bayes_factor(propLinFit6, studentFit6)
bayes_factor(propLinFit6, linFit61)
bayes_factor(propLinFit6, propFit6)
bayes_factor(propLinFit6, logPropFit6)
loo(propLinFit6, gammaFit6)

propLinFit6smooths <- conditional_smooths(propLinFit6)
plot(propLinFit6smooths, stype = "raster", ask = FALSE)
propLinFit6effects <- conditional_effects(propLinFit6, 
                                           method = "posterior_predict",
                                           robust = FALSE,
                                           re_formula = NULL)
propLinFit6effects <- conditional_effects(propLinFit6)
plot(propLinFit6effects, points = TRUE, ask = FALSE)

#### Prediction ----
## Fitted
propLinFit6finalFit <- posterior_predict(propLinFit6)
propLinFit6finalFit2 <- t(t(exp(propLinFit6finalFit))*StormdataTrain3$HWRF)
#propLinFit6finalFitMean <- colMeans(propLinFit6finalFit)*StormdataTrain3$HWRF
propLinFit6finalFitMean <- colMeans(propLinFit6finalFit2)
propLinFit6finalFitMed <- apply(propLinFit6finalFit2, 2, function(x){quantile(x, 0.5)})
propLinFit6finalFitLCB <- apply(propLinFit6finalFit2, 2, function(x){quantile(x, 0.025)})
propLinFit6finalFitUCB <- apply(propLinFit6finalFit2, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
propLinFit6finalPreds <- posterior_predict(propLinFit6, 
                                            newdata = StormdataTestFinalscale2,
                                            allow_new_levels = TRUE)
propLinFit6finalPreds <- t(t(exp(propLinFit6finalPreds))*StormdataTest3$HWRF)
propLinFit6finalPreds2 <- colMeans(propLinFit6finalPreds)
propLinFit6finalPredsMed <- apply(propLinFit6finalPreds, 2, function(x){quantile(x, 0.5)})
propLinFit6finalPredsLCB <- apply(propLinFit6finalPreds, 2, function(x){quantile(x, 0.025)})
propLinFit6finalPredsUCB <- apply(propLinFit6finalPreds, 2, function(x){quantile(x, 0.975)})

propLinFit6predMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(propLinFit6finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(propLinFit6finalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < propLinFit6finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(propLinFit6finalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(propLinFit6finalPredsMed - Actual_Yvec)),
  COV_pred = mean(propLinFit6finalPredsLCB < Actual_Yvec & Actual_Yvec < propLinFit6finalPredsUCB)
)
propLinFit6predMetrics

#### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = propLinFit6finalPreds) +
  labs(title = "propLinFit6 Predict")

propLinFit6FitDF <- bind_cols(
  StormdataTrain3,
  LCB = propLinFit6finalFitLCB,
  Mean = propLinFit6finalFitMean,
  Med = propLinFit6finalFitMed,
  UCB = propLinFit6finalFitUCB
) 

ggplot(data = propLinFit6FitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))

## Prediction
propLinFit6PredDF <- bind_cols(
  StormdataTest3,
  LCB = propLinFit6finalPredsLCB,
  Mean = propLinFit6finalPreds2,
  Med = propLinFit6finalPredsMed,
  UCB = propLinFit6finalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

ggplot(data = propLinFit6PredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))

rm(propLinFit6finalFit)
rm(propLinFit6finalFit2)
rm(propLinFit6finalPreds)

### Model 7 ----
propLinFit7 <- brm(
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
save(propLinFit7, file = "_data/propLinFit7.RData")
prior_summary(propLinFit7)
round(posterior_summary(propLinFit7, probs = c(0.025, 0.975)))
propLinFit7

print(propLinFit7, digits = 4)
pp_check(propLinFit7, ndraws = 100)
plot(propLinFit7)
loo(propLinFit7)
waic(propLinFit7)
performance::check_distribution(propLinFit7)
performance::check_outliers(propLinFit7)
performance::check_heteroskedasticity(propLinFit7)
performance_rmse(propLinFit7)
performance_mae(propLinFit7)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(propLinFit7)


variance_decomposition(propLinFit7)
exp(fixef(propLinFit7))
ranef(propLinFit7)

bayes_R2(propLinFit7)
loo_R2(propLinFit7)

bayes_factor(propLinFit7, propLinFit1)
bayes_factor(propLinFit7, gammaFit7)
bayes_factor(propLinFit7, gammaFit7)
bayes_factor(propLinFit7, studentFit7)
bayes_factor(propLinFit7, studentFit7)
bayes_factor(propLinFit7, studentFit7)
bayes_factor(propLinFit7, linFit71)
bayes_factor(propLinFit7, propFit7)
bayes_factor(propLinFit7, logPropFit7)
loo(propLinFit7, gammaFit7)

propLinFit7smooths <- conditional_smooths(propLinFit7)
plot(propLinFit7smooths, stype = "raster", ask = FALSE)
propLinFit7effects <- conditional_effects(propLinFit7, 
                                          method = "posterior_predict",
                                          robust = FALSE,
                                          re_formula = NULL)
propLinFit7effects <- conditional_effects(propLinFit7)
plot(propLinFit7effects, points = TRUE, ask = FALSE)

#### Prediction ----
## Fitted
propLinFit7finalFit <- posterior_predict(propLinFit7)
propLinFit7finalFit2 <- t(t(propLinFit7finalFit)*StormdataTrain3$HWRF)
#propLinFit7finalFitMean <- colMeans(propLinFit7finalFit)*StormdataTrain3$HWRF
propLinFit7finalFitMean <- colMeans(propLinFit7finalFit2)
propLinFit7finalFitMed <- apply(propLinFit7finalFit2, 2, function(x){quantile(x, 0.5)})
propLinFit7finalFitLCB <- apply(propLinFit7finalFit2, 2, function(x){quantile(x, 0.025)})
propLinFit7finalFitUCB <- apply(propLinFit7finalFit2, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
propLinFit7finalPreds <- posterior_predict(propLinFit7, 
                                           newdata = StormdataTestFinalscale2,
                                           allow_new_levels = TRUE)
propLinFit7finalPreds <- t(t(propLinFit7finalPreds)*StormdataTest3$HWRF)
propLinFit7finalPreds2 <- colMeans(propLinFit7finalPreds)
propLinFit7finalPredsMed <- apply(propLinFit7finalPreds, 2, function(x){quantile(x, 0.5)})
propLinFit7finalPredsLCB <- apply(propLinFit7finalPreds, 2, function(x){quantile(x, 0.025)})
propLinFit7finalPredsUCB <- apply(propLinFit7finalPreds, 2, function(x){quantile(x, 0.975)})

propLinFit7predMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(propLinFit7finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(propLinFit7finalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < propLinFit7finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(propLinFit7finalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(propLinFit7finalPredsMed - Actual_Yvec)),
  COV_pred = mean(propLinFit7finalPredsLCB < Actual_Yvec & Actual_Yvec < propLinFit7finalPredsUCB)
)
propLinFit7predMetrics

#### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = propLinFit7finalPreds) +
  labs(title = "propLinFit7 Predict")

propLinFit7FitDF <- bind_cols(
  StormdataTrain3,
  LCB = propLinFit7finalFitLCB,
  Mean = propLinFit7finalFitMean,
  Med = propLinFit7finalFitMed,
  UCB = propLinFit7finalFitUCB
) 

ggplot(data = propLinFit7FitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))

## Prediction
propLinFit7PredDF <- bind_cols(
  StormdataTest3,
  LCB = propLinFit7finalPredsLCB,
  Mean = propLinFit7finalPreds2,
  Med = propLinFit7finalPredsMed,
  UCB = propLinFit7finalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

ggplot(data = propLinFit7PredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))

rm(propLinFit7finalFit)
rm(propLinFit7finalFit2)
rm(propLinFit7finalPreds)

## GAUSSIAN PROCESS ----
### Model 1 ----
gpFit1 <- brm(
  formula = VMAX ~
    #Year +
    #Month +
    basin + 
    s(Day) +
    gp(LON, LAT, StormElapsedTime, scale = TRUE) + 
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
  data = StormdataTrain7, 
  family = gaussian(), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)
save(gpFit1, file = "_data/gpFit1.RData")
prior_summary(gpFit1)
round(posterior_summary(gpFit1, probs = c(0.025, 0.975)))
gpFit1

print(gpFit1, digits = 4)
plot(gpFit1)
pp_check(gpFit1, ndraws = 100)
loo(gpFit1)
waic(gpFit1)
performance::check_distribution(gpFit1)
performance::check_outliers(gpFit1)
performance::check_heteroskedasticity(gpFit1)
performance_rmse(gpFit1)
performance_mae(gpFit1)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(gpFit1)


variance_decomposition(gpFit1)
exp(fixef(gpFit1))
ranef(gpFit1)

bayes_R2(gpFit1)

bayes_factor(gpFit1, gammaFit1)
bayes_factor(gpFit1, gammaFit2)
bayes_factor(gpFit1, gammaFit3)
bayes_factor(gpFit1, studentFit1)
bayes_factor(gpFit1, studentFit2)
bayes_factor(gpFit1, studentFit3)
bayes_factor(gpFit1, linFit11)
bayes_factor(gpFit1, propFit1)
bayes_factor(gpFit1, logPropFit1)
loo(gpFit1, gammaFit3)

gpFit1smooths <- conditional_smooths(gpFit1)
plot(gpFit1smooths, stype = "raster", ask = FALSE)
gpFit1effects <- conditional_effects(gpFit1, 
                                     method = "posterior_predict",
                                     robust = FALSE,
                                     re_formula = NULL)
gpFit1effects <- conditional_effects(gpFit1)
plot(gpFit1effects, points = TRUE, ask = FALSE)

#### Prediction ----
## Fitted
gpFit1finalFit <- posterior_predict(gpFit1)
gpFit1finalFitMean <- colMeans(gpFit1finalFit)
gpFit1finalFitMed <- apply(gpFit1finalFit, 2, function(x){quantile(x, 0.5)})
gpFit1finalFitLCB <- apply(gpFit1finalFit, 2, function(x){quantile(x, 0.025)})
gpFit1finalFitUCB <- apply(gpFit1finalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
gpFit1finalPreds <- posterior_predict(gpFit1, 
                                      newdata = StormdataTestFinalscale,
                                      allow_new_levels = TRUE)
gpFit1finalPreds2 <- colMeans(gpFit1finalPreds)
gpFit1finalPredsMed <- apply(gpFit1finalPreds, 2, function(x){quantile(x, 0.5)})
gpFit1finalPredsLCB <- apply(gpFit1finalPreds, 2, function(x){quantile(x, 0.025)})
gpFit1finalPredsUCB <- apply(gpFit1finalPreds, 2, function(x){quantile(x, 0.975)})

gpFit1predMetrics <- tibble(
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(gpFit1finalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(gpFit1finalFitLCB < StormdataTrain7$VMAX & StormdataTrain7$VMAX < gpFit1finalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(gpFit1finalPreds2 - Actual_Yvec)),
  MAD_pred = mean(abs(gpFit1finalPredsMed - Actual_Yvec)),
  COV_pred = mean(gpFit1finalPredsLCB < Actual_Yvec & Actual_Yvec < gpFit1finalPredsUCB)
)
gpFit1predMetrics

#### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = gpFit1finalPreds) +
  labs(title = "GammaFit5 Predict")

gpFit1FitDF <- bind_cols(
  StormdataTrain3,
  LCB = gpFit1finalFitLCB,
  Mean = gpFit1finalFitMean,
  Med = gpFit1finalFitMed,
  UCB = gpFit1finalFitUCB
) |>
  filter(StormID %in% c(5712015))

ggplot(data = gpFit1FitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))

## Prediction
gpFit1PredDF <- bind_cols(
  StormdataTest3,
  LCB = gpFit1finalPredsLCB,
  Mean = gpFit1finalPreds2,
  Med = gpFit1finalPredsMed,
  UCB = gpFit1finalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

ggplot(data = gpFit1PredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))

rm(gpFit1finalFit)
rm(gpFit1finalPreds)

# Compare Predictions ----
linFit10loo <- loo(linFit10)
linFit11loo <- loo(linFit11)
linFit12loo <- loo(linFit12)
linFit13loo <- loo(linFit13)
studentFit1loo <- loo(studentFit1)
studentFit2loo <- loo(studentFit2)
studentFit3loo <- loo(studentFit3)
gammaFit0loo <- loo(gammaFit0)
gammaFit1loo <- loo(gammaFit1)
gammaFit2loo <- loo(gammaFit2)
gammaFit3loo <- loo(gammaFit3)
gammaFit4loo <- loo(gammaFit4)
gammaFit5loo <- loo(gammaFit5)
# gammaFit5Bloo <- loo(gammaFit5B)
# gammaFit5Cloo <- loo(gammaFit5C)
# gammaFit6loo <- loo(gammaFit6)
# gammaFit7loo <- loo(gammaFit7)
propLinFit1loo <- loo(propLinFit1)
propLinFit2loo <- loo(propLinFit2)
propLinFit3loo <- loo(propLinFit3)
propLinFit3Bloo <- loo(propLinFit3B)
propLinFit3Cloo <- loo(propLinFit3C)
propLinFit4loo <- loo(propLinFit4)
propLinFit5loo <- loo(propLinFit5)
propLinFit6loo <- loo(propLinFit6)
propLinFit7loo <- loo(propLinFit7)

looComp <- loo_compare(linFit10loo,
                       linFit11loo,
                       linFit12loo,
                       linFit13loo,
                       studentFit1loo,
                       studentFit2loo,
                       studentFit3loo,
                       gammaFit0loo,
                       gammaFit1loo,
                       gammaFit2loo,
                       gammaFit3loo,
                       gammaFit4loo,
                       gammaFit5loo,
                       #gammaFit5Bloo,
                       #gammaFit5Cloo,
                       #gammaFit6loo,
                       #gammaFit7loo,
                       #propLinFit1loo,
                       propLinFit2loo,
                       propLinFit3loo,
                       propLinFit3Bloo,
                       propLinFit3Cloo,
                       #propLinFit4loo,
                       propLinFit5loo,
                       propLinFit6loo,
                       propLinFit7loo)
looComp
save(looComp, file = "_data/looComp.RData")

predCompMetrics <- bind_rows(
  linFit10predMetrics |> bind_cols(Fit = "linFit10"),
  linFit11predMetrics |> bind_cols(Fit = "linFit11"),
  linFit12predMetrics |> bind_cols(Fit = "linFit12"),
  linFit13predMetrics |> bind_cols(Fit = "linFit13"),
  studentFit1predMetrics |> bind_cols(Fit = "studentFit1"),
  studentFit2predMetrics |> bind_cols(Fit = "studentFit2"),
  studentFit3predMetrics |> bind_cols(Fit = "studentFit3"),
  gammaFit0predMetrics |> bind_cols(Fit = "gammaFit0"),
  gammaFit1predMetrics |> bind_cols(Fit = "gammaFit1"),
  gammaFit2predMetrics |> bind_cols(Fit = "gammaFit2"),
  gammaFit3predMetrics |> bind_cols(Fit = "gammaFit3"),
  gammaFit4predMetrics |> bind_cols(Fit = "gammaFit4"),
  gammaFit5predMetrics |> bind_cols(Fit = "gammaFit5"),
  # gammaFit5BpredMetrics |> bind_cols(Fit = "gammaFit5B"),
  # gammaFit5CpredMetrics |> bind_cols(Fit = "gammaFit5C"),
  # gammaFit6predMetrics |> bind_cols(Fit = "gammaFit6"),
  # gammaFit7predMetrics |> bind_cols(Fit = "gammaFit7"),
  # propLinFit1predMetrics |> bind_cols(Fit = "propLinFit1"),
  propLinFit2predMetrics |> bind_cols(Fit = "propLinFit2"),
  propLinFit3predMetrics |> bind_cols(Fit = "propLinFit3"),
  propLinFit3BpredMetrics |> bind_cols(Fit = "propLinFit3B"),
  propLinFit3CpredMetrics |> bind_cols(Fit = "propLinFit3C"),
  #propLinFit4predMetrics |> bind_cols(Fit = "propLinFit4"),
  propLinFit5predMetrics |> bind_cols(Fit = "propLinFit5"),
  propLinFit6predMetrics |> bind_cols(Fit = "propLinFit6"),
  propLinFit7predMetrics |> bind_cols(Fit = "propLinFit7")
)
predCompMetrics <- predCompMetrics |> arrange(MAE_pred)

save(linFit13, file = "~/Desktop/linFit13.RData")
save(propLinFit2, file = "~/Desktop/propLinFit2.RData")
save(propLinFit3, file = "~/Desktop/propLinFit3.RData")
save(propLinFit3B, file = "~/Desktop/propLinFit3B.RData")
save(propLinFit3C, file = "~/Desktop/propLinFit3C.RData")
save(propLinFit5, file = "~/Desktop/propLinFit5.RData")
save(predCompMetrics, file = "_data/predCompMetrics.RData")

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



