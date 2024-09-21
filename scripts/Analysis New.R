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
plot(hist(dataYearDay))

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
    across(where(is.numeric) & !c(VMAX, HWRF, StormElapsedTime, LAT, LON),
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

### Actual ----
Actual_Y <- fread("_data/Actual Y.csv")
Actual_Yvec <- Actual_Y |> filter(complete.cases(x)) |> pull(x)


# Plot VMAX ----
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
world_coordinates <- map_data("world") 
ggplot() + 
  # geom_map() function takes world coordinates  
  # as input to plot world map 
  geom_map( 
    data = world_coordinates, map = world_coordinates, 
    aes(map_id = region) 
  ) + 
  geom_point(
    data = StormdataTrain3,
    aes(x = LON-360, y = LAT, 
        color = StormElapsedTime)
  ) +
  xlim(c(-180,0)) +
  ylim(c(0,60)) +
  scale_color_continuous(low = "green", high = "red") +
  theme_bw()

# Fit model ----
load(file = "_data/gammaFit0.RData")
load(file = "_data/gammaFit1.RData")
load(file = "_data/gammaFit2.RData")
load(file = "_data/gammaFit3.RData")
load(file = "_data/gammaFit4.RData")
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

ppc_dens_overlay(y = Actual_Yvec, yrep = studentFit1finalPreds) +
  labs(title = "studentFit1 Predict")
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
    Month +
    #s(StormElapsedTime) + 
    #I(StormElapsedTime^2) +
    #s(LAT, LON) +
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

ppc_dens_overlay(y = Actual_Yvec, yrep = gammaFit3finalPreds) +
  labs(title = "GammaFit3 Predict")
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
  )

ggplot(data = gammaFit5PredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))

rm(gammaFit5finalFit)
rm(gammaFit5finalPreds)

### Model 6 ----
gammaFit6 <- brm(
  formula = VMAX ~
    #Year +
    #Month +
    #basin + 
    s(Day, bs = "cc") +
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
  warmup = 1000,
  knots = list(Day = c(scale(0,
                             center = attr(StormdataTrain7scale |> pull(Day), "scaled:center"),
                             scale = attr(StormdataTrain7scale |> pull(Day), "scaled:scale")),
                       scale(365,
                             center = attr(StormdataTrain7scale |> pull(Day), "scaled:center"),
                             scale = attr(StormdataTrain7scale |> pull(Day), "scaled:scale"))
                       ))
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

gammaFit6smooths <- conditional_smooths(gammaFit6, spaghetti = TRUE, ndraws = 100,)
plot(gammaFit6smooths, stype = "raster")
gammaFit6effects <- conditional_effects(gammaFit6, 
                                        method = "posterior_predict",
                                        robust = FALSE,
                                        re_formula = NULL)
#spaghetti = TRUE,
#ndraws = 100)
plot(gammaFit6effects, points = TRUE)

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

### Compare Predictions ----
linFit10loo <- loo(linFit10)
linFit11loo <- loo(linFit11)
linFit12loo <- loo(linFit12)
studentFit1loo <- loo(studentFit1)
studentFit2loo <- loo(studentFit2)
studentFit3loo <- loo(studentFit3)
gammaFit0loo <- loo(gammaFit0)
gammaFit1loo <- loo(gammaFit1)
gammaFit2loo <- loo(gammaFit2)
gammaFit3loo <- loo(gammaFit3)
gammaFit4loo <- loo(gammaFit4)
gammaFit5loo <- loo(gammaFit5)
gammaFit6loo <- loo(gammaFit6)

loo_compare(linFit10loo,
            linFit11loo,
            linFit12loo,
            studentFit1loo,
            studentFit2loo,
            studentFit3loo,
            gammaFit0loo,
            gammaFit1loo,
            gammaFit2loo,
            gammaFit3loo,
            gammaFit4loo,
            gammaFit5loo,
            gammaFit6loo)

predCompMetrics <- bind_rows(
  linFit10predMetrics |> bind_cols(Fit = "linFit10"),
  linFit11predMetrics |> bind_cols(Fit = "linFit11"),
  linFit12predMetrics |> bind_cols(Fit = "linFit12"),
  studentFit1predMetrics |> bind_cols(Fit = "studentFit1"),
  studentFit2predMetrics |> bind_cols(Fit = "studentFit2"),
  studentFit3predMetrics |> bind_cols(Fit = "studentFit3"),
  gammaFit0predMetrics |> bind_cols(Fit = "gammaFit0"),
  gammaFit1predMetrics |> bind_cols(Fit = "gammaFit1"),
  gammaFit2predMetrics |> bind_cols(Fit = "gammaFit2"),
  gammaFit3predMetrics |> bind_cols(Fit = "gammaFit3"),
  gammaFit4predMetrics |> bind_cols(Fit = "gammaFit4"),
  gammaFit5predMetrics |> bind_cols(Fit = "gammaFit5"),
  gammaFit6predMetrics |> bind_cols(Fit = "gammaFit6")
)
predCompMetrics |> arrange(MAE_pred)

### MGCV ----
### Model 5 ----
set.seed(52)
gammaFit5 <- gam(
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

### Model 6 ----
gammaFit6 <- gam(
  formula = VMAX ~
    #Year +
    Month +
    basin + 
    #t2(LON, LAT, StormElapsedTime, bs = c("tp","tp"), k = c(20,10)) +
    #s(StormElapsedTime, bs = "cr") +
    t2(LON, LAT, StormElapsedTime, bs = c("cr", "cr", "cr")) +
    s(StormID, bs = "re") +
    #t2(LON, LAT, StormElapsedTime, StormID, bs = c("cr", "cr", "cr","re")) +
    #t2(LON, LAT, StormElapsedTime, StormID) +
    #bs = c("cr", "cr", "cr"), k = c(20, 10, 6)) +
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

plot(gammaFit6, pages = 1, residuals = TRUE)
summary(gammaFit6)
anova.gam(gammaFit6)
anova(gammaFit6, gammaFit5, test = "Chisq")
gam.check(gammaFit6)
k.check(gammaFit6)
performance_mae(gammaFit6)
gammaFit6Preds <- predict(gammaFit6, newdata = StormdataTestFinal, type = "response")
mean(abs(gammaFit6Preds - Actual_Yvec))

### Model 7 ----
gammaFit7 <- gam(
  formula = VMAX ~
    #Year +
    Month +
    basin + 
    #t2(LON, LAT, StormElapsedTime, bs = c("tp","tp"), k = c(20,10)) +
    #s(StormElapsedTime, bs = "cr") +
    t2(LON, LAT, StormElapsedTime, bs = c("tp", "tp", "cr")) +
    s(StormID, bs = "re") +
    #t2(LON, LAT, StormElapsedTime, StormID, bs = c("cr", "cr", "cr","re")) +
    #t2(LON, LAT, StormElapsedTime, StormID) +
    #bs = c("cr", "cr", "cr"), k = c(20, 10, 7)) +
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

plot(gammaFit7)
summary(gammaFit7)
anova.gam(gammaFit7)
anova(gammaFit7, gammaFit5, test = "Chisq")
anova(gammaFit7, gammaFit6, test = "Chisq")
gam.check(gammaFit7)
k.check(gammaFit7)
performance_mae(gammaFit7)
gammaFit7Preds <- predict(gammaFit7, newdata = StormdataTestFinal, type = "response")
mean(abs(gammaFit7Preds - Actual_Yvec))

### Model 8 ----
set.seed(52)
gammaFit8 <- gam(
  formula = VMAX ~
    #Year +
    Month +
    basin + 
    #t2(LON, LAT, StormElapsedTime, bs = c("tp","tp"), k = c(20,10)) +
    #s(StormElapsedTime, bs = "cr") +
    t2(LON, LAT, StormElapsedTime, d = c(2,1), bs = c("tp", "tp")) +
    s(StormID, bs = "re") +
    #t2(LON, LAT, StormElapsedTime, StormID, bs = c("cr", "cr", "cr","re")) +
    #t2(LON, LAT, StormElapsedTime, StormID) +
    #bs = c("cr", "cr", "cr"), k = c(20, 10, 8)) +
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

plot(gammaFit8,residuals = TRUE)
summary(gammaFit8)
anova.gam(gammaFit8)
anova(gammaFit8, gammaFit5, test = "Chisq")
anova(gammaFit8, gammaFit6, test = "Chisq")
anova(gammaFit8, gammaFit7, test = "Chisq")
gam.check(gammaFit8)
k.check(gammaFit8)
performance_mae(gammaFit8)
gammaFit8Preds <- predict(gammaFit8, newdata = StormdataTestFinal, type = "response")
mean(abs(gammaFit8Preds - Actual_Yvec))


### Model 9 ----
set.seed(52)
gammaFit9 <- gam(
  formula = VMAX ~
    #Year +
    Month +
    basin + 
    t2(LON, LAT, StormElapsedTime, d = c(2,1), bs = c("tp", "cr")) +
    s(StormID, bs = "re") +
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

plot(gammaFit9)
summary(gammaFit9)
anova.gam(gammaFit9)
anova(gammaFit9, gammaFit5, test = "Chisq")
anova(gammaFit9, gammaFit6, test = "Chisq")
anova(gammaFit9, gammaFit7, test = "Chisq")
anova(gammaFit9, gammaFit8, test = "Chisq")
gam.check(gammaFit9)
k.check(gammaFit9)
performance_mae(gammaFit9)
gammaFit9Preds <- predict(gammaFit9, newdata = StormdataTestFinal, type = "response")
mean(abs(gammaFit9Preds - Actual_Yvec))

### Model 10 ----
set.seed(52)
gammaFit10 <- gam(
  formula = VMAX ~
    #Year +
    #Month +
    #basin + 
    t2(LON, LAT, StormElapsedTime, d = c(2,1), bs = c("tp", "cr")) +
    s(StormID, bs = "re") +
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

plot(gammaFit10)
summary(gammaFit10)
anova.gam(gammaFit10)
anova(gammaFit10, gammaFit7, test = "Chisq")
anova(gammaFit10, gammaFit8, test = "Chisq")
anova(gammaFit10, gammaFit9, test = "Chisq")
gam.check(gammaFit10)
k.check(gammaFit10)
performance_mae(gammaFit10)
gammaFit10Preds <- predict(gammaFit10, newdata = StormdataTestFinal, type = "response")
mean(abs(gammaFit10Preds - Actual_Yvec))

### Model 11 ----
set.seed(52)
gammaFit11 <- gam(
  formula = VMAX ~
    #Year +
    #Month +
    #basin + 
    t2(LON, LAT, StormElapsedTime, d = c(2,1), bs = c("tp", "cr"), k = c(30,5)) +
    s(StormID, bs = "re") +
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

plot(gammaFit11)
summary(gammaFit11)
anova.gam(gammaFit11)
anova(gammaFit11, gammaFit7, test = "Chisq")
anova(gammaFit11, gammaFit9, test = "Chisq")
anova(gammaFit11, gammaFit10, test = "Chisq")
gam.check(gammaFit11)
k.check(gammaFit11)
performance_mae(gammaFit11)
gammaFit11Preds <- predict(gammaFit11, newdata = StormdataTestFinal, type = "response")
mean(abs(gammaFit11Preds - Actual_Yvec))

### Model 12 ----
set.seed(52)
gammaFit12 <- gam(
  formula = VMAX ~
    #Year +
    #Month +
    #basin + 
    #t2(LON, LAT, StormElapsedTime, d = c(2,1), bs = c("tp", "cr")) +
    #s(StormElapsedTime, by = StormID, bs = "cr") +
    s(StormElapsedTime, bs = c("cr")) +
    t2(StormElapsedTime, StormID, bs = c("cr", "re")) +
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
  data = StormdataTrain7B, 
  family = Gamma(link = "log"),
  method = "REML"
)

plot(gammaFit12)
summary(gammaFit12)
anova.gam(gammaFit12)
anova(gammaFit12, gammaFit7, test = "Chisq")
anova(gammaFit12, gammaFit8, test = "Chisq")
anova(gammaFit12, gammaFit11, test = "Chisq")
gam.check(gammaFit12)
k.check(gammaFit12)
performance_mae(gammaFit12)
gammaFit12Preds <- predict(gammaFit12, newdata = StormdataTestFinal, type = "response")
mean(abs(gammaFit12Preds - Actual_Yvec))

### Model 13 ----
set.seed(52)
gammaFit13 <- gam(
  formula = VMAX ~
    #Year +
    #Month +
    #basin + 
    #t2(LON, LAT, StormElapsedTime, d = c(2,1), bs = c("tp", "cr")) +
    #s(StormElapsedTime, by = StormID, bs = "cr") +
    #s(StormElapsedTime, bs = c("cr")) +
    #s(StormElapsedTime, bs = "cr") +
    s(StormElapsedTime, by = StormID, bs = c("cr")) +
    s(StormID, bs = c("re")) +
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
  data = StormdataTrain7B, 
  family = Gamma(link = "log"),
  method = "REML"
)

plot(gammaFit13)
summary(gammaFit13)
anova.gam(gammaFit13)
anova(gammaFit13, gammaFit7, test = "Chisq")
anova(gammaFit13, gammaFit11, test = "Chisq")
anova(gammaFit13, gammaFit12, test = "Chisq")
gam.check(gammaFit13)
k.check(gammaFit13)
performance_mae(gammaFit13)
gammaFit13Preds <- predict(gammaFit13, newdata = StormdataTestFinal, type = "response")
mean(abs(gammaFit13Preds - Actual_Yvec))

### Model 14 ----
set.seed(52)
gammaFit14 <- gam(
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
    HWRF,# +
  #s(StormID, bs = "re"),
  data = StormdataTrain7B, 
  family = Gamma(link = "log"),
  method = "REML"
)

plot(gammaFit14)
summary(gammaFit14)
anova.gam(gammaFit14)
anova(gammaFit14, gammaFit11, test = "Chisq")
anova(gammaFit14, gammaFit12, test = "Chisq")
anova(gammaFit14, gammaFit13, test = "Chisq")
gam.check(gammaFit14)
k.check(gammaFit14)
performance_mae(gammaFit14)
gammaFit14Preds <- predict(gammaFit14, newdata = StormdataTestFinal, type = "response")
mean(abs(gammaFit14Preds - Actual_Yvec))

### Model 15 ----
set.seed(52)
gammaFit15 <- gam(
  formula = VMAX ~
    #Year +
    #Month +
    #basin + 
    t2(LON, LAT, StormElapsedTime, d = c(2,1), bs = c("tp", "cs")) +
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
  data = StormdataTrain7B, 
  family = Gamma(link = "log"),
  method = "REML"
)

plot(gammaFit15)
summary(gammaFit15)
anova.gam(gammaFit15)
anova(gammaFit15, gammaFit11, test = "Chisq")
anova(gammaFit15, gammaFit12, test = "Chisq")
anova(gammaFit15, gammaFit14, test = "Chisq")
gam.check(gammaFit15)
k.check(gammaFit15)
performance_mae(gammaFit15)
gammaFit15Preds <- predict(gammaFit15, newdata = StormdataTestFinal, type = "response")
mean(abs(gammaFit15Preds - Actual_Yvec))



## GAUSSIAN PROP ----

### Model 1 ----
propFit1 <- brm(
  formula = propVMAX ~ 
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
    #HWRF +
    (1|StormID),
  data = StormdataTrain8, 
  family = gaussian(), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52, 
  warmup = 1000
)

prior_summary(propFit1)
posterior_summary(propFit1)
propFit1

print(propFit1, digits = 4)
plot(propFit1)
pp_check(propFit1, ndraws = 100)
loo(propFit1, studentFit1)
waic(propFit1)
performance::check_distribution(propFit1)
performance::check_outliers(propFit1)
performance::check_heteroskedasticity(propFit1)
performance_rmse(propFit1)
performance_mae(propFit1)
mean(abs(StormdataTrain8$VMAX - StormdataTrain8$HWRF))
model_performance(propFit1)


variance_decomposition(propFit1)
exp(fixef(propFit1))
ranef(propFit1)

bayes_R2(propFit1)

bayes_factor(propFit1, propFit4)
bayes_factor(propFit1, propFit5)
bayes_factor(propFit1, studentFit1)

conditional_smooths(propFit1)
conditional_effects(propFit1)

proppostEPreds1 <- posterior_epred(propFit1)
proppostEPreds2 <- colMeans(proppostEPreds1)
mean(abs(proppostEPreds2 - StormdataTrain3$VMAX))

propPredsMean <- apply(proppostEPreds1, 2, function(x){mean(x)})
propPredsMed <- apply(proppostEPreds1, 2, function(x){quantile(x, 0.5)})
propPredsLCB <- apply(proppostEPreds1, 2, function(x){quantile(x, 0.025)})
propPredsUCB <- apply(proppostEPreds1, 2, function(x){quantile(x, 0.975)})

predMetrics <- tibble(
  MAE_HWRF = mean(abs(StormdataTest3$HWRF - Actual_Yvec)),
  MAE_fit = mean(abs(propPreds2 - Actual_Yvec)),
  MAD_fit = mean(abs(propPredsMed - Actual_Yvec)),
  COV = mean(propPredsLCB < Actual_Yvec & Actual_Yvec < propPredsUCB)
)
predMetrics

ppc_dens_overlay(y = Actual_Yvec, yrep = finalPreds)

### Model 2 (log) ----
logPropFit1 <- brm(
  formula = log(propVMAX) ~ 
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
    #HWRF +
    (1|StormID),
  data = StormdataTrain8, 
  family = gaussian(), 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52, 
  warmup = 1000
)
prior_summary(logPropFit1)
posterior_summary(logPropFit1)
logPropFit1

print(logPropFit1, digits = 4)
plot(logPropFit1)
pp_check(logPropFit1, ndraws = 100)
loo(logpropFit1)
loo(logPropFit1, gammaFit1)
waic(logPropFit1)
performance::check_distribution(logPropFit1)
performance::check_outliers(logPropFit1)
performance::check_heteroskedasticity(logPropFit1)
performance_rmse(logPropFit1)
performance_mae(logPropFit1)
mean(abs(StormdataTrain8$VMAX - StormdataTrain8$HWRF))
model_performance(logPropFit1)
logPropvs <- varsel(logPropFit1)
varImp(logPropFit1)

variance_decomposition(logPropFit1)
exp(fixef(logPropFit1))
ranef(logPropFit1)

bayes_R2(logPropFit1)

bayes_factor(logPropFit1, logPropFit4)
bayes_factor(logPropFit1, logPropFit5)
bayes_factor(logPropFit1, studentFit1)

conditional_smooths(logPropFit1)
conditional_effects(logPropFit1)

logPropEPreds1 <- posterior_epred(logPropFit1)
logPropEPreds2 <- colMeans(logPropEPreds1)
logPropEPreds3 <- exp(logPropEPreds2)*StormdataTrain8$HWRF
mean(abs(logPropEPreds3 - StormdataTrain3$VMAX))

logPropPredsMean <- apply(logPropEPreds1, 2, function(x){mean(x)})
logPropPredsMed <- apply(logPropEPreds1, 2, function(x){quantile(x, 0.5)})
logPropPredsLCB <- apply(logPropEPreds1, 2, function(x){quantile(x, 0.025)})
logPropPredsUCB <- apply(logPropEPreds1, 2, function(x){quantile(x, 0.975)})

predMetrics <- tibble(
  MAE_HWRF = mean(abs(StormdataTest3$HWRF - Actual_Yvec)),
  MAE_fit = mean(abs(logPropPreds2 - Actual_Yvec)),
  MAD_fit = mean(abs(logPropPredsMed - Actual_Yvec)),
  COV = mean(logPropPredsLCB < Actual_Yvec & Actual_Yvec < logPropPredsUCB)
)
predMetrics

ppc_dens_overlay(y = Actual_Yvec, yrep = finalPreds)



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



