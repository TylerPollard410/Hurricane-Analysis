# Load Libraries ----
library(knitr)
library(data.table)
library(MASS)
library(rjags)
library(plyr)
library(stringr)
library(lubridate)
library(sf)
library(spData)
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

## Clean data ----
Stormdata1 <- Stormdata_raw |>
  mutate(
    Obs = 1:nrow(Stormdata_raw),
    StormID = factor(StormID),
    basin = factor(basin),
    Date = as_datetime(Date, tz = "UTC"),
    Year = year(Date),
    Month = month(Date, label = TRUE),
    Day = yday(Date),
    LON2 = LON - 360,
    DataType = ifelse(is.na(VMAX), "Test", "Train")
  ) |>
  group_by(StormID) |>
  mutate(
    StormElapsedTime = as.numeric(difftime(Date, min(Date), units = "hours")),
    StormElapsedTime2 = StormElapsedTime/6
  ) |>
  ungroup() |>
  select(
    DataType,
    StormID,
    Date,
    Year,
    Month,
    Day,
    basin,
    StormElapsedTime,
    StormElapsedTime2,
    LAT,
    LON,
    LON2,
    everything(),
    -lead_time
  )

# Create Date vars
# dataYears <- year(Stormdata$Date)
# dataMonths <- month(Stormdata$Date, label = TRUE)
# dataDays <- day(Stormdata$Date)
# dataYearDay <- yday(Stormdata$Date)

### Land Indicator ----
pts <- st_as_sf(Stormdata1, # |> select(LON, LAT), 
                coords = c("LON2", "LAT"),
                crs = 4326
)
land_pts <- !is.na(as.numeric(st_intersects(pts, world)))

Stormdata <- Stormdata1 |>
  mutate(
    Land = factor(land_pts, 
                  labels = c("Water", "Land"))
  ) |>
  select(
    DataType,
    StormID,
    Date,
    Year,
    Month,
    Day,
    basin,
    StormElapsedTime,
    StormElapsedTime2,
    LAT,
    LON,
    LON2,
    Land,
    everything()
  )


## Training ----
StormdataTrain1 <- Stormdata |> 
  filter(DataType == "Train") |>
  select(-DataType)

#### Train 3 ----
StormdataTrain3 <- StormdataTrain1 |>
  mutate(
    StormID = droplevels(StormID),
    Year = factor(Year, ordered = TRUE)
    #Month = dataTrainMonths,
    #Day = dataTrainYearDay
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
  )

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
    "Land",
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
    #StormID = droplevels(StormID),
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
    "Land",
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
    #StormID = droplevels(StormID),
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
    "Land",
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
    #StormID = droplevels(StormID),
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

#### Train 9 ----
StormdataTrain9 <- StormdataTrain3 |>
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
    "Land",
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
    #StormID = droplevels(StormID),
    Year = factor(Year, ordered = FALSE),
    Month = factor(Month, ordered = FALSE),
    logHWRF = log(HWRF)
  ) |>
  mutate(
    across(where(is.numeric) & !c(VMAX, 
                                  HWRF, logHWRF,
                                  StormElapsedTime,StormElapsedTime2,
                                  LAT, LON),
           function(x){scale(x)}),
    logHWRFscale = scale(logHWRF, center = TRUE, scale = FALSE)
  )
str(StormdataTrain9)

## Test ----
StormdataTest1 <- Stormdata |> 
  filter(DataType == "Test") |>
  select(-DataType)

# # Remove not varying 
# StormdataTest2 <- StormdataTest1 |>
#   select(-lead_time)
# 
# # Create Date vars 
# dataTestYears <- year(StormdataTest2$Date)
# dataTestMonths <- month(StormdataTest2$Date, label = TRUE)
# dataTestDays <- day(StormdataTest2$Date)
# dataTestYearDay <- yday(StormdataTest2$Date)

#### Test 3 ----
StormdataTest3 <- StormdataTest1 |>
  mutate(
    StormID = droplevels(StormID),
    Year = factor(Year, ordered = TRUE)
    #Month = dataTestMonths,
    #Day = dataTestYearDay
  ) |>
  # group_by(StormID) |>
  # mutate(
  #   StormElapsedTime = as.numeric(difftime(Date, min(Date), units = "hours")),
  #   StormElapsedTime2 = StormElapsedTime/6
  # ) |>
  select(
    StormID,
    Date,
    Year,
    Month,
    Day,
    StormElapsedTime,
    StormElapsedTime2,
    everything()
  )

#### Test 7 ----
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
    "Land",
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
    # StormID = droplevels(StormID),
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

#### Test 7 Scale----
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
    "Land",
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
    #StormID = droplevels(StormID),
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

#### Test 8 ----
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
    "Land",
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
    #StormID = droplevels(StormID),
    Year = factor(Year, ordered = FALSE),
    Month = factor(Month, ordered = FALSE)
  ) |>
  mutate(
    across(where(is.numeric) & !c(VMAX, HWRF, StormElapsedTime2),
           function(x){scale(x,
                             center = attr(StormdataTrain8 |> pull(x), "scaled:center"),
                             scale = attr(StormdataTrain8 |> pull(x), "scaled:scale"))
           })
  )
str(StormdataTest8)

#### Test 9 ----
StormdataTest9 <- StormdataTest3 |>
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
    "Land",
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
    #StormID = droplevels(StormID),
    Year = factor(Year, ordered = FALSE),
    Month = factor(Month, ordered = FALSE),
    logHWRF = log(HWRF)
  ) |>
  mutate(
    across(where(is.numeric) & !c(VMAX, 
                                  HWRF, logHWRF,
                                  StormElapsedTime,StormElapsedTime2,
                                  LAT, LON),
           function(x){scale(x)}),
    logHWRFscale = scale(logHWRF, center = TRUE, scale = FALSE)
  )
str(StormdataTest9)

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
ggpairs(StormdataTrain3 |>
          select(
            #"StormID",
            #Date,
            #Year,
            #Month,
            Day,
            StormElapsedTime,
            #StormElapsedTime2,
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
            "VMAX"))

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
mean(log(StormdataTrain3$VMAX)/log(StormdataTrain3$HWRF))
sd(log(StormdataTrain3$VMAX/StormdataTrain3$HWRF))

## Histogram ----
ggplot(data = StormdataTrain3) +
  geom_histogram(
    aes(x = VMAX, after_stat(density)),
    color = "#99c7c7", fill = "#bcdcdc") +
  geom_vline(aes(xintercept = mean(VMAX)),
             color = "orange", linewidth = 2) +
  geom_vline(aes(xintercept = median(VMAX)),
             color = "orange3", linewidth = 2) +
  geom_density(#data = final_data3,
    aes(x = VMAX),
    color = "#007C7C", 
    linewidth = 1) +
  theme_bw()

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
    aes(x = log(VMAX) - log(HWRF), after_stat(density)),
    color = "#99c7c7", fill = "#bcdcdc",
    bins = 100) +
  geom_vline(aes(xintercept = mean(log(VMAX/HWRF))),
             color = "orange", linewidth = 2) +
  geom_vline(aes(xintercept = median(log(VMAX/HWRF))),
             color = "orange3", linewidth = 2) +
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
    aes(x = log(rlnorm(1705, 0.0242, 0.2203))),
    color = "gold", 
    linewidth = 1) +
  geom_density(#data = final_data3,
    aes(x = (rlnorm(1705, 0.0242, 0.2203))),
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

ggplot(data = StormdataTrain8) +
  geom_histogram(
    aes(x = VMAX, after_stat(density)),
    color = "#99c7c7", fill = "#bcdcdc", bins = 30) +
  geom_vline(aes(xintercept = mean(VMAX)),
             color = "orange", linewidth = 2) +
  geom_vline(aes(xintercept = median(VMAX)),
             color = "orange3", linewidth = 2) +
  geom_density(#data = final_data3,
    aes(x = VMAX),
    color = "#007C7C", 
    linewidth = 1) +
  geom_density(#data = final_data3,
    aes(x = exp(rnorm(1705, 1.049843 - 1, 0.239613/sqrt(2)) + log(HWRF))),
    color = "blue", 
    linewidth = 1) +
  geom_density(#data = final_data3,
    aes(x = exp(log(rlnorm(1705, 0.0242, 0.2203))+ log(HWRF))),
    color = "gold", 
    linewidth = 1) +
  geom_density(#data = final_data3,
    aes(x = exp((rlnorm(1705, 0.0242, 0.2203))+ log(HWRF))),
    color = "red", 
    linewidth = 1) +
  geom_density(#data = final_data3,
    aes(x = exp(log(rlnorm(1705, 1.049843-0.239613^2/2, 0.239613))+ log(HWRF))),
    color = "purple", 
    linewidth = 1) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  scale_x_continuous(limits = c(0,180)) +
  labs(title = "Density Plot",
       subtitle = "Data",
       x =  "VMAX",
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
  geom_point(aes(x = StormElapsedTime, y = VMAX)) +
  scale_x_continuous(breaks = seq(0,402,24))

ggplot(data = StormdataTrain3) +
  geom_line(aes(x = StormElapsedTime, y = log(VMAX))) +
  scale_x_continuous(breaks = seq(0,402,24)) +
  facet_wrap(vars(StormID))

ggplot(data = StormdataTrain3) +
  geom_line(aes(x = log(StormElapsedTime), y = log(VMAX) - log(HWRF))) +
  #scale_x_continuous(breaks = seq(0,402,24)) +
  facet_wrap(vars(StormID))

ggplot(data = StormdataTrain3) +
  geom_line(aes(x = StormElapsedTime, y = log(VMAX/HWRF))) +
  #scale_x_continuous(breaks = seq(0,402,24)) +
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
world_coordinates2 <- st_as_sf(world_coordinates,
                               coords = c("long", "lat"))
ggplot() + 
  # geom_map() function takes world coordinates  
  # as input to plot world map 
  geom_map( 
    data = world_coordinates, map = world_coordinates, 
    aes(map_id = region) 
  ) + 
  geom_point(
    data = Stormdata |> filter(complete.cases(VMAX)),
    aes(x = LON2, y = LAT, 
        color = VMAX,
        shape = Land
    )
  ) +
  xlim(c(-190,0)) +
  ylim(c(0,60)) +
  scale_color_continuous(low = "green", high = "red") +
  #facet_wrap(vars(DataType), nrow = 2) +
  theme_bw()

# ggplot() + 
#   # geom_map() function takes world coordinates  
#   # # as input to plot world map 
#   # geom_map( 
#   #   data = world_coordinates, map = world_coordinates, 
#   #   aes(map_id = region) 
#   # ) + 
#   geom_sf(
#     data = world_coordinates2
#     #aes(geometry = geometry)
#   ) +
#   geom_sf(
#     data = pts, # |> filter(complete.cases(VMAX)),
#     aes(geometry = geometry, 
#         color = VMAX
#     )
#   ) +
#   coord_sf(xlim = c(-190,0), ylim = c(0,60)) +
#   #xlim(c(-190,0)) +
#   #ylim(c(0,60)) +
#   #scale_color_continuous(low = "green", high = "red") +
#   #facet_wrap(vars(DataType), nrow = 2) +
#   theme_bw()

# Find Distributions ----

ddist <- descdist(StormdataTrain3$VMAX, discrete = FALSE)
summary(ddist)

ddistIdentity <- descdist(log(StormdataTrain3$VMAX), discrete = FALSE)
summary(ddist)

ddistIdentityProp <- descdist(log(StormdataTrain3$VMAX/StormdataTrain3$HWRF), discrete = FALSE)
summary(ddist)

fitnorm1 <- fitdist(StormdataTrain3$VMAX, distr = "norm", method = "mle")
summary(fitnorm1)
fitdist(StormdataTrain3$VMAX, distr = "norm", method = "mme")
fitdist(StormdataTrain3$VMAX, distr = "norm", method = "mle")
fitdist(StormdataTrain3$VMAX, distr = "norm", method = "mle")

# LOGNORMAL ----
## NULL Models ----
### No Scaling ----
#### Identity Link ----
logNormalFitNULL <- brm(
  bf(
    VMAX ~ HWRF
  ),
  data = StormdataTrain8, 
  family = lognormal(link = "identity"),
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)

##### Plot ----
print(logNormalFitNULL, digits = 4)
logNormalFitNULLppcFit <- pp_check(logNormalFitNULL, ndraws = 100) + 
  labs(title = "logNormalFitNULL Fit PPC") +
  theme_bw()
logNormalFitNULLppcFit

##### LOO ----
logNormalFitNULLloo <- loo(logNormalFitNULL)

##### Prediction ----
## Fitted
logNormalFitNULLFit <- posterior_predict(logNormalFitNULL)
logNormalFitNULLFitMean <- colMeans(logNormalFitNULLFit)
logNormalFitNULLFitMed <- apply(logNormalFitNULLFit, 2, function(x){quantile(x, 0.5)})
logNormalFitNULLFitLCB <- apply(logNormalFitNULLFit, 2, function(x){quantile(x, 0.025)})
logNormalFitNULLFitUCB <- apply(logNormalFitNULLFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
logNormalFitNULLPreds <- posterior_predict(logNormalFitNULL, 
                                           newdata = StormdataTest8,
                                           allow_new_levels = TRUE)
logNormalFitNULLPredsMean <- colMeans(logNormalFitNULLPreds)
logNormalFitNULLPredsMed <- apply(logNormalFitNULLPreds, 2, function(x){quantile(x, 0.5)})
logNormalFitNULLPredsLCB <- apply(logNormalFitNULLPreds, 2, function(x){quantile(x, 0.025)})
logNormalFitNULLPredsUCB <- apply(logNormalFitNULLPreds, 2, function(x){quantile(x, 0.975)})

logNormalFitNULLpredMetrics <- tibble(
  Fit = "logNormalFitNULL",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(logNormalFitNULLFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(logNormalFitNULLFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < logNormalFitNULLFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest3$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(logNormalFitNULLPredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(logNormalFitNULLPredsMed - Actual_Yvec)),
  COV_pred = mean(logNormalFitNULLPredsLCB < Actual_Yvec & Actual_Yvec < logNormalFitNULLPredsUCB)
)
logNormalFitNULLpredMetrics

##### PPC ----
###### Quantile 2.5 
logNormalFitNULLLCBsims <- apply(logNormalFitNULLFit, 
                                 MARGIN = 1,
                                 function(x){
                                   quantile(x, 0.025)
                                 })
logNormalFitNULLLCBpvalueVec <- logNormalFitNULLLCBsims < quantile(StormdataTrain3$VMAX, 0.025)
logNormalFitNULLLCBpvalue <- sum(logNormalFitNULLLCBpvalueVec)
logNormalFitNULLLCBpvalue <- round(logNormalFitNULLLCBpvalue/4000, 3)
logNormalFitNULLLCBpvalue <- min(logNormalFitNULLLCBpvalue, 1 - logNormalFitNULLLCBpvalue)

logNormalFitNULL_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", logNormalFitNULLLCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULL_ppcLCB

###### Quantile 97.5 
logNormalFitNULLUCBsims <- apply(logNormalFitNULLFit, 
                                 MARGIN = 1,
                                 function(x){
                                   quantile(x, 0.975)
                                 })
logNormalFitNULLUCBpvalueVec <- logNormalFitNULLUCBsims < quantile(StormdataTrain3$VMAX, 0.975)
logNormalFitNULLUCBpvalue <- as.numeric(sum(logNormalFitNULLUCBpvalueVec))
logNormalFitNULLUCBpvalue <- round(logNormalFitNULLUCBpvalue/4000, 3)
logNormalFitNULLUCBpvalue <- min(logNormalFitNULLUCBpvalue, 1 - logNormalFitNULLUCBpvalue)

logNormalFitNULL_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", logNormalFitNULLUCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULL_ppcUCB

###### Mean 
logNormalFitNULLMEANsims <- apply(logNormalFitNULLFit, 
                                  MARGIN = 1,
                                  function(x){
                                    mean(x)
                                  })
logNormalFitNULLMEANpvalueVec <- logNormalFitNULLMEANsims < mean(StormdataTrain3$VMAX)
logNormalFitNULLMEANpvalue <- sum(logNormalFitNULLMEANpvalueVec)
logNormalFitNULLMEANpvalue <- round(logNormalFitNULLMEANpvalue/4000, 3)
logNormalFitNULLMEANpvalue <- min(logNormalFitNULLMEANpvalue, 1 - logNormalFitNULLMEANpvalue)

logNormalFitNULL_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", logNormalFitNULLMEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULL_ppcMEAN

###### Med 
logNormalFitNULLMEDsims <- apply(logNormalFitNULLFit, 
                                 MARGIN = 1,
                                 function(x){
                                   quantile(x, 0.5)
                                 })
logNormalFitNULLMEDpvalueVec <- logNormalFitNULLMEDsims < quantile(StormdataTrain3$VMAX, 0.5)
logNormalFitNULLMEDpvalue <- sum(logNormalFitNULLMEDpvalueVec)
logNormalFitNULLMEDpvalue <- round(logNormalFitNULLMEDpvalue/4000, 3)
logNormalFitNULLMEDpvalue <- min(logNormalFitNULLMEDpvalue, 1 - logNormalFitNULLMEDpvalue)

logNormalFitNULL_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", logNormalFitNULLMEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULL_ppcMED

###### SD 
logNormalFitNULLSDsims <- apply(logNormalFitNULLFit, 
                                MARGIN = 1,
                                function(x){
                                  sd(x)
                                })
logNormalFitNULLSDpvalueVec <- logNormalFitNULLSDsims < sd(StormdataTrain3$VMAX)
logNormalFitNULLSDpvalue <- sum(logNormalFitNULLSDpvalueVec)
logNormalFitNULLSDpvalue <- round(logNormalFitNULLSDpvalue/4000, 3)
logNormalFitNULLSDpvalue <- min(logNormalFitNULLSDpvalue, 1 - logNormalFitNULLSDpvalue)

logNormalFitNULL_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", logNormalFitNULLSDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULL_ppcSD

###### Range 
logNormalFitNULLRANGEsims <- apply(logNormalFitNULLFit, 
                                   MARGIN = 1,
                                   function(x){
                                     max(x)-min(x)
                                   })
logNormalFitNULLRANGEpvalueVec <- logNormalFitNULLRANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
logNormalFitNULLRANGEpvalue <- sum(logNormalFitNULLRANGEpvalueVec)
logNormalFitNULLRANGEpvalue <- round(logNormalFitNULLRANGEpvalue/4000, 3)
logNormalFitNULLRANGEpvalue <- min(logNormalFitNULLRANGEpvalue, 1 - logNormalFitNULLRANGEpvalue)

logNormalFitNULL_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", logNormalFitNULLRANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULL_ppcRANGE

##### Combined Plot ----
logNormalFitNULL_ppcComb <- 
  logNormalFitNULLppcFit /
  (logNormalFitNULL_ppcLCB | logNormalFitNULL_ppcMED | logNormalFitNULL_ppcUCB) /
  (logNormalFitNULL_ppcRANGE | logNormalFitNULL_ppcMEAN | logNormalFitNULL_ppcSD)
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
logNormalFitNULL_ppcComb

##### Bayes p-values ----
logNormalFitNULLpvalues <- tibble(
  Fit = paste0("logNormalFitNULL"),
  LCB = logNormalFitNULLLCBpvalue,
  Median = logNormalFitNULLMEDpvalue,
  UCB = logNormalFitNULLUCBpvalue,
  Range = logNormalFitNULLRANGEpvalue,
  Mean = logNormalFitNULLMEANpvalue,
  SD = logNormalFitNULLSDpvalue
)
logNormalFitNULLpvalues

##### CV ----
kfoldID <- kfold_split_grouped(K = 5, StormdataTrain8$StormID)
logNormalFitNULLkfoldgroup <- kfold(logNormalFitNULL,
                                    folds = kfoldID,
                                    chains = 1,
                                    save_fits = TRUE)
logNormalFitNULLkfoldrand <- kfold(logNormalFitNULL,
                                   K = 5,
                                   chains = 1,
                                   save_fits = TRUE)
logNormalFitNULLkfoldPreds <- kfold_predict(logNormalFitNULLkfoldgroup)
#logNormalFitNULLkfoldPreds <- kfold_predict(logNormalFitNULLkfold)
logNormalFitNULLkfoldPredsDat <- logNormalFitNULLkfoldPreds$yrep
logNormalFitNULLkfoldPredsMean <- colMeans(logNormalFitNULLkfoldPredsDat)
logNormalFitNULLkfoldPredsMed <- apply(logNormalFitNULLkfoldPredsDat, 2, function(x){quantile(x, 0.5)})
logNormalFitNULLkfoldPredsLCB <- apply(logNormalFitNULLkfoldPredsDat, 2, function(x){quantile(x, 0.025)})
logNormalFitNULLkfoldPredsUCB <- apply(logNormalFitNULLkfoldPredsDat, 2, function(x){quantile(x, 0.975)})

logNormalFitNULLkfoldMetrics <- tibble(
  Fit = paste0("logNormalFitNULL"),
  MAE_kfold = mean(abs(logNormalFitNULLkfoldPredsMean - logNormalFitNULLkfoldPreds$y)),
  MAD_kfold = mean(abs(logNormalFitNULLkfoldPredsMed - logNormalFitNULLkfoldPreds$y)),
  COV_kfold = mean(logNormalFitNULLkfoldPredsLCB < logNormalFitNULLkfoldPreds$y & logNormalFitNULLkfoldPreds$y < logNormalFitNULLkfoldPredsUCB)
)
logNormalFitNULLkfoldMetrics

#### Log Sigma ----
logNormalFitNULLlogSigma <- brm(
  bf(
    VMAX ~ HWRF,
    sigma ~ 1
  ),
  data = StormdataTrain8, 
  family = lognormal(link = "identity", link_sigma = "log"),
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)

##### Plot ----
print(logNormalFitNULLlogSigma, digits = 4)
logNormalFitNULLlogSigmappcFit <- pp_check(logNormalFitNULLlogSigma, ndraws = 100) + 
  labs(title = "logNormalFitNULLlogSigma Fit PPC") +
  theme_bw()
logNormalFitNULLlogSigmappcFit

##### LOO ----
logNormalFitNULLlogSigmaloo <- loo(logNormalFitNULLlogSigma)

##### Prediction ----
## Fitted
logNormalFitNULLlogSigmaFit <- posterior_predict(logNormalFitNULLlogSigma)
logNormalFitNULLlogSigmaFitMean <- colMeans(logNormalFitNULLlogSigmaFit)
logNormalFitNULLlogSigmaFitMed <- apply(logNormalFitNULLlogSigmaFit, 2, function(x){quantile(x, 0.5)})
logNormalFitNULLlogSigmaFitLCB <- apply(logNormalFitNULLlogSigmaFit, 2, function(x){quantile(x, 0.025)})
logNormalFitNULLlogSigmaFitUCB <- apply(logNormalFitNULLlogSigmaFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
logNormalFitNULLlogSigmaPreds <- posterior_predict(logNormalFitNULLlogSigma, 
                                                   newdata = StormdataTest8,
                                                   allow_new_levels = TRUE)
logNormalFitNULLlogSigmaPredsMean <- colMeans(logNormalFitNULLlogSigmaPreds)
logNormalFitNULLlogSigmaPredsMed <- apply(logNormalFitNULLlogSigmaPreds, 2, function(x){quantile(x, 0.5)})
logNormalFitNULLlogSigmaPredsLCB <- apply(logNormalFitNULLlogSigmaPreds, 2, function(x){quantile(x, 0.025)})
logNormalFitNULLlogSigmaPredsUCB <- apply(logNormalFitNULLlogSigmaPreds, 2, function(x){quantile(x, 0.975)})

logNormalFitNULLlogSigmapredMetrics <- tibble(
  Fit = "logNormalFitNULLlogSigma",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(logNormalFitNULLlogSigmaFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(logNormalFitNULLlogSigmaFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < logNormalFitNULLlogSigmaFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest3$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(logNormalFitNULLlogSigmaPredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(logNormalFitNULLlogSigmaPredsMed - Actual_Yvec)),
  COV_pred = mean(logNormalFitNULLlogSigmaPredsLCB < Actual_Yvec & Actual_Yvec < logNormalFitNULLlogSigmaPredsUCB)
)
logNormalFitNULLlogSigmapredMetrics

##### PPC ----
###### Quantile 2.5 
logNormalFitNULLlogSigmaLCBsims <- apply(logNormalFitNULLlogSigmaFit, 
                                         MARGIN = 1,
                                         function(x){
                                           quantile(x, 0.025)
                                         })
logNormalFitNULLlogSigmaLCBpvalueVec <- logNormalFitNULLlogSigmaLCBsims < quantile(StormdataTrain3$VMAX, 0.025)
logNormalFitNULLlogSigmaLCBpvalue <- sum(logNormalFitNULLlogSigmaLCBpvalueVec)
logNormalFitNULLlogSigmaLCBpvalue <- round(logNormalFitNULLlogSigmaLCBpvalue/4000, 3)
logNormalFitNULLlogSigmaLCBpvalue <- min(logNormalFitNULLlogSigmaLCBpvalue, 1 - logNormalFitNULLlogSigmaLCBpvalue)

logNormalFitNULLlogSigma_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLlogSigmaFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", logNormalFitNULLlogSigmaLCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLlogSigma_ppcLCB

###### Quantile 97.5 
logNormalFitNULLlogSigmaUCBsims <- apply(logNormalFitNULLlogSigmaFit, 
                                         MARGIN = 1,
                                         function(x){
                                           quantile(x, 0.975)
                                         })
logNormalFitNULLlogSigmaUCBpvalueVec <- logNormalFitNULLlogSigmaUCBsims < quantile(StormdataTrain3$VMAX, 0.975)
logNormalFitNULLlogSigmaUCBpvalue <- as.numeric(sum(logNormalFitNULLlogSigmaUCBpvalueVec))
logNormalFitNULLlogSigmaUCBpvalue <- round(logNormalFitNULLlogSigmaUCBpvalue/4000, 3)
logNormalFitNULLlogSigmaUCBpvalue <- min(logNormalFitNULLlogSigmaUCBpvalue, 1 - logNormalFitNULLlogSigmaUCBpvalue)

logNormalFitNULLlogSigma_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLlogSigmaFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", logNormalFitNULLlogSigmaUCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLlogSigma_ppcUCB

###### Mean 
logNormalFitNULLlogSigmaMEANsims <- apply(logNormalFitNULLlogSigmaFit, 
                                          MARGIN = 1,
                                          function(x){
                                            mean(x)
                                          })
logNormalFitNULLlogSigmaMEANpvalueVec <- logNormalFitNULLlogSigmaMEANsims < mean(StormdataTrain3$VMAX)
logNormalFitNULLlogSigmaMEANpvalue <- sum(logNormalFitNULLlogSigmaMEANpvalueVec)
logNormalFitNULLlogSigmaMEANpvalue <- round(logNormalFitNULLlogSigmaMEANpvalue/4000, 3)
logNormalFitNULLlogSigmaMEANpvalue <- min(logNormalFitNULLlogSigmaMEANpvalue, 1 - logNormalFitNULLlogSigmaMEANpvalue)

logNormalFitNULLlogSigma_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLlogSigmaFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", logNormalFitNULLlogSigmaMEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLlogSigma_ppcMEAN

###### Med 
logNormalFitNULLlogSigmaMEDsims <- apply(logNormalFitNULLlogSigmaFit, 
                                         MARGIN = 1,
                                         function(x){
                                           quantile(x, 0.5)
                                         })
logNormalFitNULLlogSigmaMEDpvalueVec <- logNormalFitNULLlogSigmaMEDsims < quantile(StormdataTrain3$VMAX, 0.5)
logNormalFitNULLlogSigmaMEDpvalue <- sum(logNormalFitNULLlogSigmaMEDpvalueVec)
logNormalFitNULLlogSigmaMEDpvalue <- round(logNormalFitNULLlogSigmaMEDpvalue/4000, 3)
logNormalFitNULLlogSigmaMEDpvalue <- min(logNormalFitNULLlogSigmaMEDpvalue, 1 - logNormalFitNULLlogSigmaMEDpvalue)

logNormalFitNULLlogSigma_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLlogSigmaFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", logNormalFitNULLlogSigmaMEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLlogSigma_ppcMED

###### SD 
logNormalFitNULLlogSigmaSDsims <- apply(logNormalFitNULLlogSigmaFit, 
                                        MARGIN = 1,
                                        function(x){
                                          sd(x)
                                        })
logNormalFitNULLlogSigmaSDpvalueVec <- logNormalFitNULLlogSigmaSDsims < sd(StormdataTrain3$VMAX)
logNormalFitNULLlogSigmaSDpvalue <- sum(logNormalFitNULLlogSigmaSDpvalueVec)
logNormalFitNULLlogSigmaSDpvalue <- round(logNormalFitNULLlogSigmaSDpvalue/4000, 3)
logNormalFitNULLlogSigmaSDpvalue <- min(logNormalFitNULLlogSigmaSDpvalue, 1 - logNormalFitNULLlogSigmaSDpvalue)

logNormalFitNULLlogSigma_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLlogSigmaFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", logNormalFitNULLlogSigmaSDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLlogSigma_ppcSD

###### Range 
logNormalFitNULLlogSigmaRANGEsims <- apply(logNormalFitNULLlogSigmaFit, 
                                           MARGIN = 1,
                                           function(x){
                                             max(x)-min(x)
                                           })
logNormalFitNULLlogSigmaRANGEpvalueVec <- logNormalFitNULLlogSigmaRANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
logNormalFitNULLlogSigmaRANGEpvalue <- sum(logNormalFitNULLlogSigmaRANGEpvalueVec)
logNormalFitNULLlogSigmaRANGEpvalue <- round(logNormalFitNULLlogSigmaRANGEpvalue/4000, 3)
logNormalFitNULLlogSigmaRANGEpvalue <- min(logNormalFitNULLlogSigmaRANGEpvalue, 1 - logNormalFitNULLlogSigmaRANGEpvalue)

logNormalFitNULLlogSigma_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLlogSigmaFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", logNormalFitNULLlogSigmaRANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLlogSigma_ppcRANGE

##### Combined Plot ----
logNormalFitNULLlogSigma_ppcComb <- 
  logNormalFitNULLlogSigmappcFit /
  (logNormalFitNULLlogSigma_ppcLCB | logNormalFitNULLlogSigma_ppcMED | logNormalFitNULLlogSigma_ppcUCB) /
  (logNormalFitNULLlogSigma_ppcRANGE | logNormalFitNULLlogSigma_ppcMEAN | logNormalFitNULLlogSigma_ppcSD)
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
logNormalFitNULLlogSigma_ppcComb

##### Bayes p-values ----
logNormalFitNULLlogSigmapvalues <- tibble(
  Fit = paste0("logNormalFitNULLlogSigma"),
  LCB = logNormalFitNULLlogSigmaLCBpvalue,
  Median = logNormalFitNULLlogSigmaMEDpvalue,
  UCB = logNormalFitNULLlogSigmaUCBpvalue,
  Range = logNormalFitNULLlogSigmaRANGEpvalue,
  Mean = logNormalFitNULLlogSigmaMEANpvalue,
  SD = logNormalFitNULLlogSigmaSDpvalue
)
logNormalFitNULLlogSigmapvalues

##### CV ----
kfoldID <- kfold_split_grouped(K = 5, StormdataTrain8$StormID)
logNormalFitNULLlogSigmakfoldgroup <- kfold(logNormalFitNULLlogSigma,
                                            folds = kfoldID,
                                            chains = 1,
                                            save_fits = TRUE)
# logNormalFitNULLlogSigmakfoldrand <- kfold(logNormalFitNULLlogSigma,
#                                 K = 5,
#                                 chains = 1,
#                                 save_fits = TRUE)
logNormalFitNULLlogSigmakfoldPreds <- kfold_predict(logNormalFitNULLlogSigmakfoldgroup)
#logNormalFitNULLlogSigmakfoldPreds <- kfold_predict(logNormalFitNULLlogSigmakfold)
logNormalFitNULLlogSigmakfoldPredsDat <- logNormalFitNULLlogSigmakfoldPreds$yrep
logNormalFitNULLlogSigmakfoldPredsMean <- colMeans(logNormalFitNULLlogSigmakfoldPredsDat)
logNormalFitNULLlogSigmakfoldPredsMed <- apply(logNormalFitNULLlogSigmakfoldPredsDat, 2, function(x){quantile(x, 0.5)})
logNormalFitNULLlogSigmakfoldPredsLCB <- apply(logNormalFitNULLlogSigmakfoldPredsDat, 2, function(x){quantile(x, 0.025)})
logNormalFitNULLlogSigmakfoldPredsUCB <- apply(logNormalFitNULLlogSigmakfoldPredsDat, 2, function(x){quantile(x, 0.975)})

logNormalFitNULLlogSigmakfoldMetrics <- tibble(
  Fit = paste0("logNormalFitNULLlogSigma"),
  MAE_kfold = mean(abs(logNormalFitNULLlogSigmakfoldPredsMean - logNormalFitNULLlogSigmakfoldPreds$y)),
  MAD_kfold = mean(abs(logNormalFitNULLlogSigmakfoldPredsMed - logNormalFitNULLlogSigmakfoldPreds$y)),
  COV_kfold = mean(logNormalFitNULLlogSigmakfoldPredsLCB < logNormalFitNULLlogSigmakfoldPreds$y & logNormalFitNULLlogSigmakfoldPreds$y < logNormalFitNULLlogSigmakfoldPredsUCB)
)
logNormalFitNULLlogSigmakfoldMetrics

#### Identity Link No Int ----
logNormalFitNULL2 <- brm(
  bf(
    VMAX ~ 0 + Intercept + HWRF
  ),
  data = StormdataTrain8, 
  family = lognormal(link = "identity"),
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)

##### Plot ----
print(logNormalFitNULL2, digits = 4)
logNormalFitNULL2ppcFit <- pp_check(logNormalFitNULL2, ndraws = 100) + 
  labs(title = "logNormalFitNULL2 Fit PPC") +
  theme_bw()
logNormalFitNULL2ppcFit

##### LOO ----
logNormalFitNULL2loo <- loo(logNormalFitNULL2)

##### Prediction ----
## Fitted
logNormalFitNULL2Fit <- posterior_predict(logNormalFitNULL2)
logNormalFitNULL2FitMean <- colMeans(logNormalFitNULL2Fit)
logNormalFitNULL2FitMed <- apply(logNormalFitNULL2Fit, 2, function(x){quantile(x, 0.5)})
logNormalFitNULL2FitLCB <- apply(logNormalFitNULL2Fit, 2, function(x){quantile(x, 0.025)})
logNormalFitNULL2FitUCB <- apply(logNormalFitNULL2Fit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
logNormalFitNULL2Preds <- posterior_predict(logNormalFitNULL2, 
                                            newdata = StormdataTest8,
                                            allow_new_levels = TRUE)
logNormalFitNULL2PredsMean <- colMeans(logNormalFitNULL2Preds)
logNormalFitNULL2PredsMed <- apply(logNormalFitNULL2Preds, 2, function(x){quantile(x, 0.5)})
logNormalFitNULL2PredsLCB <- apply(logNormalFitNULL2Preds, 2, function(x){quantile(x, 0.025)})
logNormalFitNULL2PredsUCB <- apply(logNormalFitNULL2Preds, 2, function(x){quantile(x, 0.975)})

logNormalFitNULL2predMetrics <- tibble(
  Fit = "logNormalFitNULL2",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(logNormalFitNULL2FitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(logNormalFitNULL2FitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < logNormalFitNULL2FitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest3$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(logNormalFitNULL2PredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(logNormalFitNULL2PredsMed - Actual_Yvec)),
  COV_pred = mean(logNormalFitNULL2PredsLCB < Actual_Yvec & Actual_Yvec < logNormalFitNULL2PredsUCB)
)
logNormalFitNULL2predMetrics

##### PPC ----
###### Quantile 2.5 
logNormalFitNULL2LCBsims <- apply(logNormalFitNULL2Fit, 
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
           logNormalFitNULL2Fit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", logNormalFitNULL2LCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULL2_ppcLCB

###### Quantile 97.5 
logNormalFitNULL2UCBsims <- apply(logNormalFitNULL2Fit, 
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
           logNormalFitNULL2Fit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", logNormalFitNULL2UCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULL2_ppcUCB

###### Mean 
logNormalFitNULL2MEANsims <- apply(logNormalFitNULL2Fit, 
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
           logNormalFitNULL2Fit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", logNormalFitNULL2MEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULL2_ppcMEAN

###### Med 
logNormalFitNULL2MEDsims <- apply(logNormalFitNULL2Fit, 
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
           logNormalFitNULL2Fit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", logNormalFitNULL2MEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULL2_ppcMED

###### SD 
logNormalFitNULL2SDsims <- apply(logNormalFitNULL2Fit, 
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
           logNormalFitNULL2Fit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", logNormalFitNULL2SDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULL2_ppcSD

###### Range 
logNormalFitNULL2RANGEsims <- apply(logNormalFitNULL2Fit, 
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
           logNormalFitNULL2Fit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", logNormalFitNULL2RANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULL2_ppcRANGE

##### Combined Plot ----
logNormalFitNULL2_ppcComb <- 
  logNormalFitNULL2ppcFit /
  (logNormalFitNULL2_ppcLCB | logNormalFitNULL2_ppcMED | logNormalFitNULL2_ppcUCB) /
  (logNormalFitNULL2_ppcRANGE | logNormalFitNULL2_ppcMEAN | logNormalFitNULL2_ppcSD)
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
logNormalFitNULL2_ppcComb

##### Bayes p-values ----
logNormalFitNULL2pvalues <- tibble(
  Fit = paste0("logNormalFitNULL2"),
  LCB = logNormalFitNULL2LCBpvalue,
  Median = logNormalFitNULL2MEDpvalue,
  UCB = logNormalFitNULL2UCBpvalue,
  Range = logNormalFitNULL2RANGEpvalue,
  Mean = logNormalFitNULL2MEANpvalue,
  SD = logNormalFitNULL2SDpvalue
)
logNormalFitNULL2pvalues

##### CV ----
kfoldID <- kfold_split_grouped(K = 5, StormdataTrain8$StormID)
logNormalFitNULL2kfoldgroup <- kfold(logNormalFitNULL2,
                                     folds = kfoldID,
                                     chains = 1,
                                     save_fits = TRUE)
# logNormalFitNULL2kfoldrand <- kfold(logNormalFitNULL2,
#                                 K = 5,
#                                 chains = 1,
#                                 save_fits = TRUE)
logNormalFitNULL2kfoldPreds <- kfold_predict(logNormalFitNULL2kfoldgroup)
#logNormalFitNULL2kfoldPreds <- kfold_predict(logNormalFitNULL2kfold)
logNormalFitNULL2kfoldPredsDat <- logNormalFitNULL2kfoldPreds$yrep
logNormalFitNULL2kfoldPredsMean <- colMeans(logNormalFitNULL2kfoldPredsDat)
logNormalFitNULL2kfoldPredsMed <- apply(logNormalFitNULL2kfoldPredsDat, 2, function(x){quantile(x, 0.5)})
logNormalFitNULL2kfoldPredsLCB <- apply(logNormalFitNULL2kfoldPredsDat, 2, function(x){quantile(x, 0.025)})
logNormalFitNULL2kfoldPredsUCB <- apply(logNormalFitNULL2kfoldPredsDat, 2, function(x){quantile(x, 0.975)})

logNormalFitNULL2kfoldMetrics <- tibble(
  Fit = paste0("logNormalFitNULL2"),
  MAE_kfold = mean(abs(logNormalFitNULL2kfoldPredsMean - logNormalFitNULL2kfoldPreds$y)),
  MAD_kfold = mean(abs(logNormalFitNULL2kfoldPredsMed - logNormalFitNULL2kfoldPreds$y)),
  COV_kfold = mean(logNormalFitNULL2kfoldPredsLCB < logNormalFitNULL2kfoldPreds$y & logNormalFitNULL2kfoldPreds$y < logNormalFitNULL2kfoldPredsUCB)
)
logNormalFitNULL2kfoldMetrics

#### Log Sigma ----
logNormalFitNULL2logSigma <- brm(
  bf(
    VMAX ~ 0 + Intercept + HWRF,
    sigma ~ 1
  ),
  data = StormdataTrain8, 
  family = lognormal(link = "identity", link_sigma = "log"),
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)

##### Plot ----
print(logNormalFitNULL2logSigma, digits = 4)
logNormalFitNULL2logSigmappcFit <- pp_check(logNormalFitNULL2logSigma, ndraws = 100) + 
  labs(title = "logNormalFitNULL2logSigma Fit PPC") +
  theme_bw()
logNormalFitNULL2logSigmappcFit

##### LOO ----
logNormalFitNULL2logSigmaloo <- loo(logNormalFitNULL2logSigma)

##### Prediction ----
## Fitted
logNormalFitNULL2logSigmaFit <- posterior_predict(logNormalFitNULL2logSigma)
logNormalFitNULL2logSigmaFitMean <- colMeans(logNormalFitNULL2logSigmaFit)
logNormalFitNULL2logSigmaFitMed <- apply(logNormalFitNULL2logSigmaFit, 2, function(x){quantile(x, 0.5)})
logNormalFitNULL2logSigmaFitLCB <- apply(logNormalFitNULL2logSigmaFit, 2, function(x){quantile(x, 0.025)})
logNormalFitNULL2logSigmaFitUCB <- apply(logNormalFitNULL2logSigmaFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
logNormalFitNULL2logSigmaPreds <- posterior_predict(logNormalFitNULL2logSigma, 
                                                    newdata = StormdataTest8,
                                                    allow_new_levels = TRUE)
logNormalFitNULL2logSigmaPredsMean <- colMeans(logNormalFitNULL2logSigmaPreds)
logNormalFitNULL2logSigmaPredsMed <- apply(logNormalFitNULL2logSigmaPreds, 2, function(x){quantile(x, 0.5)})
logNormalFitNULL2logSigmaPredsLCB <- apply(logNormalFitNULL2logSigmaPreds, 2, function(x){quantile(x, 0.025)})
logNormalFitNULL2logSigmaPredsUCB <- apply(logNormalFitNULL2logSigmaPreds, 2, function(x){quantile(x, 0.975)})

logNormalFitNULL2logSigmapredMetrics <- tibble(
  Fit = "logNormalFitNULL2logSigma",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(logNormalFitNULL2logSigmaFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(logNormalFitNULL2logSigmaFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < logNormalFitNULL2logSigmaFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest3$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(logNormalFitNULL2logSigmaPredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(logNormalFitNULL2logSigmaPredsMed - Actual_Yvec)),
  COV_pred = mean(logNormalFitNULL2logSigmaPredsLCB < Actual_Yvec & Actual_Yvec < logNormalFitNULL2logSigmaPredsUCB)
)
logNormalFitNULL2logSigmapredMetrics

##### PPC ----
###### Quantile 2.5 
logNormalFitNULL2logSigmaLCBsims <- apply(logNormalFitNULL2logSigmaFit, 
                                          MARGIN = 1,
                                          function(x){
                                            quantile(x, 0.025)
                                          })
logNormalFitNULL2logSigmaLCBpvalueVec <- logNormalFitNULL2logSigmaLCBsims < quantile(StormdataTrain3$VMAX, 0.025)
logNormalFitNULL2logSigmaLCBpvalue <- sum(logNormalFitNULL2logSigmaLCBpvalueVec)
logNormalFitNULL2logSigmaLCBpvalue <- round(logNormalFitNULL2logSigmaLCBpvalue/4000, 3)
logNormalFitNULL2logSigmaLCBpvalue <- min(logNormalFitNULL2logSigmaLCBpvalue, 1 - logNormalFitNULL2logSigmaLCBpvalue)

logNormalFitNULL2logSigma_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULL2logSigmaFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", logNormalFitNULL2logSigmaLCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULL2logSigma_ppcLCB

###### Quantile 97.5 
logNormalFitNULL2logSigmaUCBsims <- apply(logNormalFitNULL2logSigmaFit, 
                                          MARGIN = 1,
                                          function(x){
                                            quantile(x, 0.975)
                                          })
logNormalFitNULL2logSigmaUCBpvalueVec <- logNormalFitNULL2logSigmaUCBsims < quantile(StormdataTrain3$VMAX, 0.975)
logNormalFitNULL2logSigmaUCBpvalue <- as.numeric(sum(logNormalFitNULL2logSigmaUCBpvalueVec))
logNormalFitNULL2logSigmaUCBpvalue <- round(logNormalFitNULL2logSigmaUCBpvalue/4000, 3)
logNormalFitNULL2logSigmaUCBpvalue <- min(logNormalFitNULL2logSigmaUCBpvalue, 1 - logNormalFitNULL2logSigmaUCBpvalue)

logNormalFitNULL2logSigma_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULL2logSigmaFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", logNormalFitNULL2logSigmaUCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULL2logSigma_ppcUCB

###### Mean 
logNormalFitNULL2logSigmaMEANsims <- apply(logNormalFitNULL2logSigmaFit, 
                                           MARGIN = 1,
                                           function(x){
                                             mean(x)
                                           })
logNormalFitNULL2logSigmaMEANpvalueVec <- logNormalFitNULL2logSigmaMEANsims < mean(StormdataTrain3$VMAX)
logNormalFitNULL2logSigmaMEANpvalue <- sum(logNormalFitNULL2logSigmaMEANpvalueVec)
logNormalFitNULL2logSigmaMEANpvalue <- round(logNormalFitNULL2logSigmaMEANpvalue/4000, 3)
logNormalFitNULL2logSigmaMEANpvalue <- min(logNormalFitNULL2logSigmaMEANpvalue, 1 - logNormalFitNULL2logSigmaMEANpvalue)

logNormalFitNULL2logSigma_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULL2logSigmaFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", logNormalFitNULL2logSigmaMEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULL2logSigma_ppcMEAN

###### Med 
logNormalFitNULL2logSigmaMEDsims <- apply(logNormalFitNULL2logSigmaFit, 
                                          MARGIN = 1,
                                          function(x){
                                            quantile(x, 0.5)
                                          })
logNormalFitNULL2logSigmaMEDpvalueVec <- logNormalFitNULL2logSigmaMEDsims < quantile(StormdataTrain3$VMAX, 0.5)
logNormalFitNULL2logSigmaMEDpvalue <- sum(logNormalFitNULL2logSigmaMEDpvalueVec)
logNormalFitNULL2logSigmaMEDpvalue <- round(logNormalFitNULL2logSigmaMEDpvalue/4000, 3)
logNormalFitNULL2logSigmaMEDpvalue <- min(logNormalFitNULL2logSigmaMEDpvalue, 1 - logNormalFitNULL2logSigmaMEDpvalue)

logNormalFitNULL2logSigma_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULL2logSigmaFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", logNormalFitNULL2logSigmaMEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULL2logSigma_ppcMED

###### SD 
logNormalFitNULL2logSigmaSDsims <- apply(logNormalFitNULL2logSigmaFit, 
                                         MARGIN = 1,
                                         function(x){
                                           sd(x)
                                         })
logNormalFitNULL2logSigmaSDpvalueVec <- logNormalFitNULL2logSigmaSDsims < sd(StormdataTrain3$VMAX)
logNormalFitNULL2logSigmaSDpvalue <- sum(logNormalFitNULL2logSigmaSDpvalueVec)
logNormalFitNULL2logSigmaSDpvalue <- round(logNormalFitNULL2logSigmaSDpvalue/4000, 3)
logNormalFitNULL2logSigmaSDpvalue <- min(logNormalFitNULL2logSigmaSDpvalue, 1 - logNormalFitNULL2logSigmaSDpvalue)

logNormalFitNULL2logSigma_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULL2logSigmaFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", logNormalFitNULL2logSigmaSDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULL2logSigma_ppcSD

###### Range 
logNormalFitNULL2logSigmaRANGEsims <- apply(logNormalFitNULL2logSigmaFit, 
                                            MARGIN = 1,
                                            function(x){
                                              max(x)-min(x)
                                            })
logNormalFitNULL2logSigmaRANGEpvalueVec <- logNormalFitNULL2logSigmaRANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
logNormalFitNULL2logSigmaRANGEpvalue <- sum(logNormalFitNULL2logSigmaRANGEpvalueVec)
logNormalFitNULL2logSigmaRANGEpvalue <- round(logNormalFitNULL2logSigmaRANGEpvalue/4000, 3)
logNormalFitNULL2logSigmaRANGEpvalue <- min(logNormalFitNULL2logSigmaRANGEpvalue, 1 - logNormalFitNULL2logSigmaRANGEpvalue)

logNormalFitNULL2logSigma_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULL2logSigmaFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", logNormalFitNULL2logSigmaRANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULL2logSigma_ppcRANGE

##### Combined Plot ----
logNormalFitNULL2logSigma_ppcComb <- 
  logNormalFitNULL2logSigmappcFit /
  (logNormalFitNULL2logSigma_ppcLCB | logNormalFitNULL2logSigma_ppcMED | logNormalFitNULL2logSigma_ppcUCB) /
  (logNormalFitNULL2logSigma_ppcRANGE | logNormalFitNULL2logSigma_ppcMEAN | logNormalFitNULL2logSigma_ppcSD)
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
logNormalFitNULL2logSigma_ppcComb

##### Bayes p-values ----
logNormalFitNULL2logSigmapvalues <- tibble(
  Fit = paste0("logNormalFitNULL2logSigma"),
  LCB = logNormalFitNULL2logSigmaLCBpvalue,
  Median = logNormalFitNULL2logSigmaMEDpvalue,
  UCB = logNormalFitNULL2logSigmaUCBpvalue,
  Range = logNormalFitNULL2logSigmaRANGEpvalue,
  Mean = logNormalFitNULL2logSigmaMEANpvalue,
  SD = logNormalFitNULL2logSigmaSDpvalue
)
logNormalFitNULL2logSigmapvalues

##### CV ----
kfoldID <- kfold_split_grouped(K = 5, StormdataTrain8$StormID)
logNormalFitNULL2logSigmakfoldgroup <- kfold(logNormalFitNULL2logSigma,
                                             folds = kfoldID,
                                             chains = 1,
                                             save_fits = TRUE)
# logNormalFitNULL2logSigmakfoldrand <- kfold(logNormalFitNULL2logSigma,
#                                 K = 5,
#                                 chains = 1,
#                                 save_fits = TRUE)
logNormalFitNULL2logSigmakfoldPreds <- kfold_predict(logNormalFitNULL2logSigmakfoldgroup)
#logNormalFitNULL2logSigmakfoldPreds <- kfold_predict(logNormalFitNULL2logSigmakfold)
logNormalFitNULL2logSigmakfoldPredsDat <- logNormalFitNULL2logSigmakfoldPreds$yrep
logNormalFitNULL2logSigmakfoldPredsMean <- colMeans(logNormalFitNULL2logSigmakfoldPredsDat)
logNormalFitNULL2logSigmakfoldPredsMed <- apply(logNormalFitNULL2logSigmakfoldPredsDat, 2, function(x){quantile(x, 0.5)})
logNormalFitNULL2logSigmakfoldPredsLCB <- apply(logNormalFitNULL2logSigmakfoldPredsDat, 2, function(x){quantile(x, 0.025)})
logNormalFitNULL2logSigmakfoldPredsUCB <- apply(logNormalFitNULL2logSigmakfoldPredsDat, 2, function(x){quantile(x, 0.975)})

logNormalFitNULL2logSigmakfoldMetrics <- tibble(
  Fit = paste0("logNormalFitNULL2logSigma"),
  MAE_kfold = mean(abs(logNormalFitNULL2logSigmakfoldPredsMean - logNormalFitNULL2logSigmakfoldPreds$y)),
  MAD_kfold = mean(abs(logNormalFitNULL2logSigmakfoldPredsMed - logNormalFitNULL2logSigmakfoldPreds$y)),
  COV_kfold = mean(logNormalFitNULL2logSigmakfoldPredsLCB < logNormalFitNULL2logSigmakfoldPreds$y & logNormalFitNULL2logSigmakfoldPreds$y < logNormalFitNULL2logSigmakfoldPredsUCB)
)
logNormalFitNULL2logSigmakfoldMetrics



### Scaling ----
#### Identity Link ----
logNormalFitNULLscale <- brm(
  bf(
    VMAX ~ HWRF
  ),
  data = StormdataTrain7scale, 
  family = lognormal(link = "identity"),
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)

##### Plot ----
print(logNormalFitNULLscale, digits = 4)
logNormalFitNULLscaleppcFit <- pp_check(logNormalFitNULLscale, ndraws = 100) + 
  labs(title = "logNormalFitNULLscale Fit PPC") +
  theme_bw()
logNormalFitNULLscaleppcFit

##### LOO ----
logNormalFitNULLscaleloo <- loo(logNormalFitNULLscale)

##### Prediction ----
## Fitted
logNormalFitNULLscaleFit <- posterior_predict(logNormalFitNULLscale)
logNormalFitNULLscaleFitMean <- colMeans(logNormalFitNULLscaleFit)
logNormalFitNULLscaleFitMed <- apply(logNormalFitNULLscaleFit, 2, function(x){quantile(x, 0.5)})
logNormalFitNULLscaleFitLCB <- apply(logNormalFitNULLscaleFit, 2, function(x){quantile(x, 0.025)})
logNormalFitNULLscaleFitUCB <- apply(logNormalFitNULLscaleFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
logNormalFitNULLscalePreds <- posterior_predict(logNormalFitNULLscale, 
                                                newdata = StormdataTest7scale,
                                                allow_new_levels = TRUE)
logNormalFitNULLscalePredsMean <- colMeans(logNormalFitNULLscalePreds)
logNormalFitNULLscalePredsMed <- apply(logNormalFitNULLscalePreds, 2, function(x){quantile(x, 0.5)})
logNormalFitNULLscalePredsLCB <- apply(logNormalFitNULLscalePreds, 2, function(x){quantile(x, 0.025)})
logNormalFitNULLscalePredsUCB <- apply(logNormalFitNULLscalePreds, 2, function(x){quantile(x, 0.975)})

logNormalFitNULLscalepredMetrics <- tibble(
  Fit = "logNormalFitNULLscale",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(logNormalFitNULLscaleFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(logNormalFitNULLscaleFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < logNormalFitNULLscaleFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest3$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(logNormalFitNULLscalePredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(logNormalFitNULLscalePredsMed - Actual_Yvec)),
  COV_pred = mean(logNormalFitNULLscalePredsLCB < Actual_Yvec & Actual_Yvec < logNormalFitNULLscalePredsUCB)
)
logNormalFitNULLscalepredMetrics

##### PPC ----
###### Quantile 2.5 
logNormalFitNULLscaleLCBsims <- apply(logNormalFitNULLscaleFit, 
                                      MARGIN = 1,
                                      function(x){
                                        quantile(x, 0.025)
                                      })
logNormalFitNULLscaleLCBpvalueVec <- logNormalFitNULLscaleLCBsims < quantile(StormdataTrain3$VMAX, 0.025)
logNormalFitNULLscaleLCBpvalue <- sum(logNormalFitNULLscaleLCBpvalueVec)
logNormalFitNULLscaleLCBpvalue <- round(logNormalFitNULLscaleLCBpvalue/4000, 3)
logNormalFitNULLscaleLCBpvalue <- min(logNormalFitNULLscaleLCBpvalue, 1 - logNormalFitNULLscaleLCBpvalue)

logNormalFitNULLscale_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLscaleFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", logNormalFitNULLscaleLCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLscale_ppcLCB

###### Quantile 97.5 
logNormalFitNULLscaleUCBsims <- apply(logNormalFitNULLscaleFit, 
                                      MARGIN = 1,
                                      function(x){
                                        quantile(x, 0.975)
                                      })
logNormalFitNULLscaleUCBpvalueVec <- logNormalFitNULLscaleUCBsims < quantile(StormdataTrain3$VMAX, 0.975)
logNormalFitNULLscaleUCBpvalue <- as.numeric(sum(logNormalFitNULLscaleUCBpvalueVec))
logNormalFitNULLscaleUCBpvalue <- round(logNormalFitNULLscaleUCBpvalue/4000, 3)
logNormalFitNULLscaleUCBpvalue <- min(logNormalFitNULLscaleUCBpvalue, 1 - logNormalFitNULLscaleUCBpvalue)

logNormalFitNULLscale_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLscaleFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", logNormalFitNULLscaleUCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLscale_ppcUCB

###### Mean 
logNormalFitNULLscaleMEANsims <- apply(logNormalFitNULLscaleFit, 
                                       MARGIN = 1,
                                       function(x){
                                         mean(x)
                                       })
logNormalFitNULLscaleMEANpvalueVec <- logNormalFitNULLscaleMEANsims < mean(StormdataTrain3$VMAX)
logNormalFitNULLscaleMEANpvalue <- sum(logNormalFitNULLscaleMEANpvalueVec)
logNormalFitNULLscaleMEANpvalue <- round(logNormalFitNULLscaleMEANpvalue/4000, 3)
logNormalFitNULLscaleMEANpvalue <- min(logNormalFitNULLscaleMEANpvalue, 1 - logNormalFitNULLscaleMEANpvalue)

logNormalFitNULLscale_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLscaleFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", logNormalFitNULLscaleMEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLscale_ppcMEAN

###### Med 
logNormalFitNULLscaleMEDsims <- apply(logNormalFitNULLscaleFit, 
                                      MARGIN = 1,
                                      function(x){
                                        quantile(x, 0.5)
                                      })
logNormalFitNULLscaleMEDpvalueVec <- logNormalFitNULLscaleMEDsims < quantile(StormdataTrain3$VMAX, 0.5)
logNormalFitNULLscaleMEDpvalue <- sum(logNormalFitNULLscaleMEDpvalueVec)
logNormalFitNULLscaleMEDpvalue <- round(logNormalFitNULLscaleMEDpvalue/4000, 3)
logNormalFitNULLscaleMEDpvalue <- min(logNormalFitNULLscaleMEDpvalue, 1 - logNormalFitNULLscaleMEDpvalue)

logNormalFitNULLscale_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLscaleFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", logNormalFitNULLscaleMEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLscale_ppcMED

###### SD 
logNormalFitNULLscaleSDsims <- apply(logNormalFitNULLscaleFit, 
                                     MARGIN = 1,
                                     function(x){
                                       sd(x)
                                     })
logNormalFitNULLscaleSDpvalueVec <- logNormalFitNULLscaleSDsims < sd(StormdataTrain3$VMAX)
logNormalFitNULLscaleSDpvalue <- sum(logNormalFitNULLscaleSDpvalueVec)
logNormalFitNULLscaleSDpvalue <- round(logNormalFitNULLscaleSDpvalue/4000, 3)
logNormalFitNULLscaleSDpvalue <- min(logNormalFitNULLscaleSDpvalue, 1 - logNormalFitNULLscaleSDpvalue)

logNormalFitNULLscale_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLscaleFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", logNormalFitNULLscaleSDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLscale_ppcSD

###### Range 
logNormalFitNULLscaleRANGEsims <- apply(logNormalFitNULLscaleFit, 
                                        MARGIN = 1,
                                        function(x){
                                          max(x)-min(x)
                                        })
logNormalFitNULLscaleRANGEpvalueVec <- logNormalFitNULLscaleRANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
logNormalFitNULLscaleRANGEpvalue <- sum(logNormalFitNULLscaleRANGEpvalueVec)
logNormalFitNULLscaleRANGEpvalue <- round(logNormalFitNULLscaleRANGEpvalue/4000, 3)
logNormalFitNULLscaleRANGEpvalue <- min(logNormalFitNULLscaleRANGEpvalue, 1 - logNormalFitNULLscaleRANGEpvalue)

logNormalFitNULLscale_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLscaleFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", logNormalFitNULLscaleRANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLscale_ppcRANGE

##### Combined Plot ----
logNormalFitNULLscale_ppcComb <- 
  logNormalFitNULLscaleppcFit /
  (logNormalFitNULLscale_ppcLCB | logNormalFitNULLscale_ppcMED | logNormalFitNULLscale_ppcUCB) /
  (logNormalFitNULLscale_ppcRANGE | logNormalFitNULLscale_ppcMEAN | logNormalFitNULLscale_ppcSD)
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
logNormalFitNULLscale_ppcComb

##### Bayes p-values ----
logNormalFitNULLscalepvalues <- tibble(
  Fit = paste0("logNormalFitNULLscale"),
  LCB = logNormalFitNULLscaleLCBpvalue,
  Median = logNormalFitNULLscaleMEDpvalue,
  UCB = logNormalFitNULLscaleUCBpvalue,
  Range = logNormalFitNULLscaleRANGEpvalue,
  Mean = logNormalFitNULLscaleMEANpvalue,
  SD = logNormalFitNULLscaleSDpvalue
)
logNormalFitNULLscalepvalues

##### CV ----
kfoldID <- kfold_split_grouped(K = 5, StormdataTrain8$StormID)
logNormalFitNULLscalekfoldgroup <- kfold(logNormalFitNULLscale,
                                         folds = kfoldID,
                                         chains = 1,
                                         save_fits = TRUE)
# logNormalFitNULLscalekfoldrand <- kfold(logNormalFitNULLscale,
#                                     K = 5,
#                                     chains = 1,
#                                     save_fits = TRUE)
logNormalFitNULLscalekfoldPreds <- kfold_predict(logNormalFitNULLscalekfoldgroup)
#logNormalFitNULLscalekfoldPreds <- kfold_predict(logNormalFitNULLscalekfold)
logNormalFitNULLscalekfoldPredsDat <- logNormalFitNULLscalekfoldPreds$yrep
logNormalFitNULLscalekfoldPredsMean <- colMeans(logNormalFitNULLscalekfoldPredsDat)
logNormalFitNULLscalekfoldPredsMed <- apply(logNormalFitNULLscalekfoldPredsDat, 2, function(x){quantile(x, 0.5)})
logNormalFitNULLscalekfoldPredsLCB <- apply(logNormalFitNULLscalekfoldPredsDat, 2, function(x){quantile(x, 0.025)})
logNormalFitNULLscalekfoldPredsUCB <- apply(logNormalFitNULLscalekfoldPredsDat, 2, function(x){quantile(x, 0.975)})

logNormalFitNULLscalekfoldMetrics <- tibble(
  Fit = paste0("logNormalFitNULLscale"),
  MAE_kfold = mean(abs(logNormalFitNULLscalekfoldPredsMean - logNormalFitNULLscalekfoldPreds$y)),
  MAD_kfold = mean(abs(logNormalFitNULLscalekfoldPredsMed - logNormalFitNULLscalekfoldPreds$y)),
  COV_kfold = mean(logNormalFitNULLscalekfoldPredsLCB < logNormalFitNULLscalekfoldPreds$y & logNormalFitNULLscalekfoldPreds$y < logNormalFitNULLscalekfoldPredsUCB)
)
logNormalFitNULLscalekfoldMetrics

#### Identity Link No Int ----
logNormalFitNULLscale2 <- brm(
  bf(
    VMAX ~ 0 + Intercept + HWRF
  ),
  data = StormdataTrain7scale, 
  family = lognormal(link = "identity"),
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)

##### Plot ----
print(logNormalFitNULLscale2, digits = 4)
logNormalFitNULLscale2ppcFit <- pp_check(logNormalFitNULLscale2, ndraws = 100) + 
  labs(title = "logNormalFitNULLscale2 Fit PPC") +
  theme_bw()
logNormalFitNULLscale2ppcFit

##### LOO ----
logNormalFitNULLscale2loo <- loo(logNormalFitNULLscale2)

##### Prediction ----
## Fitted
logNormalFitNULLscale2Fit <- posterior_predict(logNormalFitNULLscale2)
logNormalFitNULLscale2FitMean <- colMeans(logNormalFitNULLscale2Fit)
logNormalFitNULLscale2FitMed <- apply(logNormalFitNULLscale2Fit, 2, function(x){quantile(x, 0.5)})
logNormalFitNULLscale2FitLCB <- apply(logNormalFitNULLscale2Fit, 2, function(x){quantile(x, 0.025)})
logNormalFitNULLscale2FitUCB <- apply(logNormalFitNULLscale2Fit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
logNormalFitNULLscale2Preds <- posterior_predict(logNormalFitNULLscale2, 
                                                 newdata = StormdataTest7scale,
                                                 allow_new_levels = TRUE)
logNormalFitNULLscale2PredsMean <- colMeans(logNormalFitNULLscale2Preds)
logNormalFitNULLscale2PredsMed <- apply(logNormalFitNULLscale2Preds, 2, function(x){quantile(x, 0.5)})
logNormalFitNULLscale2PredsLCB <- apply(logNormalFitNULLscale2Preds, 2, function(x){quantile(x, 0.025)})
logNormalFitNULLscale2PredsUCB <- apply(logNormalFitNULLscale2Preds, 2, function(x){quantile(x, 0.975)})

logNormalFitNULLscale2predMetrics <- tibble(
  Fit = "logNormalFitNULLscale2",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(logNormalFitNULLscale2FitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(logNormalFitNULLscale2FitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < logNormalFitNULLscale2FitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest3$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(logNormalFitNULLscale2PredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(logNormalFitNULLscale2PredsMed - Actual_Yvec)),
  COV_pred = mean(logNormalFitNULLscale2PredsLCB < Actual_Yvec & Actual_Yvec < logNormalFitNULLscale2PredsUCB)
)
logNormalFitNULLscale2predMetrics

##### PPC ----
###### Quantile 2.5 
logNormalFitNULLscale2LCBsims <- apply(logNormalFitNULLscale2Fit, 
                                       MARGIN = 1,
                                       function(x){
                                         quantile(x, 0.025)
                                       })
logNormalFitNULLscale2LCBpvalueVec <- logNormalFitNULLscale2LCBsims < quantile(StormdataTrain3$VMAX, 0.025)
logNormalFitNULLscale2LCBpvalue <- sum(logNormalFitNULLscale2LCBpvalueVec)
logNormalFitNULLscale2LCBpvalue <- round(logNormalFitNULLscale2LCBpvalue/4000, 3)
logNormalFitNULLscale2LCBpvalue <- min(logNormalFitNULLscale2LCBpvalue, 1 - logNormalFitNULLscale2LCBpvalue)

logNormalFitNULLscale2_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLscale2Fit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", logNormalFitNULLscale2LCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLscale2_ppcLCB

###### Quantile 97.5 
logNormalFitNULLscale2UCBsims <- apply(logNormalFitNULLscale2Fit, 
                                       MARGIN = 1,
                                       function(x){
                                         quantile(x, 0.975)
                                       })
logNormalFitNULLscale2UCBpvalueVec <- logNormalFitNULLscale2UCBsims < quantile(StormdataTrain3$VMAX, 0.975)
logNormalFitNULLscale2UCBpvalue <- as.numeric(sum(logNormalFitNULLscale2UCBpvalueVec))
logNormalFitNULLscale2UCBpvalue <- round(logNormalFitNULLscale2UCBpvalue/4000, 3)
logNormalFitNULLscale2UCBpvalue <- min(logNormalFitNULLscale2UCBpvalue, 1 - logNormalFitNULLscale2UCBpvalue)

logNormalFitNULLscale2_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLscale2Fit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", logNormalFitNULLscale2UCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLscale2_ppcUCB

###### Mean 
logNormalFitNULLscale2MEANsims <- apply(logNormalFitNULLscale2Fit, 
                                        MARGIN = 1,
                                        function(x){
                                          mean(x)
                                        })
logNormalFitNULLscale2MEANpvalueVec <- logNormalFitNULLscale2MEANsims < mean(StormdataTrain3$VMAX)
logNormalFitNULLscale2MEANpvalue <- sum(logNormalFitNULLscale2MEANpvalueVec)
logNormalFitNULLscale2MEANpvalue <- round(logNormalFitNULLscale2MEANpvalue/4000, 3)
logNormalFitNULLscale2MEANpvalue <- min(logNormalFitNULLscale2MEANpvalue, 1 - logNormalFitNULLscale2MEANpvalue)

logNormalFitNULLscale2_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLscale2Fit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", logNormalFitNULLscale2MEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLscale2_ppcMEAN

###### Med 
logNormalFitNULLscale2MEDsims <- apply(logNormalFitNULLscale2Fit, 
                                       MARGIN = 1,
                                       function(x){
                                         quantile(x, 0.5)
                                       })
logNormalFitNULLscale2MEDpvalueVec <- logNormalFitNULLscale2MEDsims < quantile(StormdataTrain3$VMAX, 0.5)
logNormalFitNULLscale2MEDpvalue <- sum(logNormalFitNULLscale2MEDpvalueVec)
logNormalFitNULLscale2MEDpvalue <- round(logNormalFitNULLscale2MEDpvalue/4000, 3)
logNormalFitNULLscale2MEDpvalue <- min(logNormalFitNULLscale2MEDpvalue, 1 - logNormalFitNULLscale2MEDpvalue)

logNormalFitNULLscale2_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLscale2Fit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", logNormalFitNULLscale2MEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLscale2_ppcMED

###### SD 
logNormalFitNULLscale2SDsims <- apply(logNormalFitNULLscale2Fit, 
                                      MARGIN = 1,
                                      function(x){
                                        sd(x)
                                      })
logNormalFitNULLscale2SDpvalueVec <- logNormalFitNULLscale2SDsims < sd(StormdataTrain3$VMAX)
logNormalFitNULLscale2SDpvalue <- sum(logNormalFitNULLscale2SDpvalueVec)
logNormalFitNULLscale2SDpvalue <- round(logNormalFitNULLscale2SDpvalue/4000, 3)
logNormalFitNULLscale2SDpvalue <- min(logNormalFitNULLscale2SDpvalue, 1 - logNormalFitNULLscale2SDpvalue)

logNormalFitNULLscale2_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLscale2Fit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", logNormalFitNULLscale2SDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLscale2_ppcSD

###### Range 
logNormalFitNULLscale2RANGEsims <- apply(logNormalFitNULLscale2Fit, 
                                         MARGIN = 1,
                                         function(x){
                                           max(x)-min(x)
                                         })
logNormalFitNULLscale2RANGEpvalueVec <- logNormalFitNULLscale2RANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
logNormalFitNULLscale2RANGEpvalue <- sum(logNormalFitNULLscale2RANGEpvalueVec)
logNormalFitNULLscale2RANGEpvalue <- round(logNormalFitNULLscale2RANGEpvalue/4000, 3)
logNormalFitNULLscale2RANGEpvalue <- min(logNormalFitNULLscale2RANGEpvalue, 1 - logNormalFitNULLscale2RANGEpvalue)

logNormalFitNULLscale2_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLscale2Fit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", logNormalFitNULLscale2RANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLscale2_ppcRANGE

##### Combined Plot ----
logNormalFitNULLscale2_ppcComb <- 
  logNormalFitNULLscale2ppcFit /
  (logNormalFitNULLscale2_ppcLCB | logNormalFitNULLscale2_ppcMED | logNormalFitNULLscale2_ppcUCB) /
  (logNormalFitNULLscale2_ppcRANGE | logNormalFitNULLscale2_ppcMEAN | logNormalFitNULLscale2_ppcSD)
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
logNormalFitNULLscale2_ppcComb

##### Bayes p-values ----
logNormalFitNULLscale2pvalues <- tibble(
  Fit = paste0("logNormalFitNULLscale2"),
  LCB = logNormalFitNULLscale2LCBpvalue,
  Median = logNormalFitNULLscale2MEDpvalue,
  UCB = logNormalFitNULLscale2UCBpvalue,
  Range = logNormalFitNULLscale2RANGEpvalue,
  Mean = logNormalFitNULLscale2MEANpvalue,
  SD = logNormalFitNULLscale2SDpvalue
)
logNormalFitNULLscale2pvalues

##### CV ----
kfoldID <- kfold_split_grouped(K = 5, StormdataTrain8$StormID)
logNormalFitNULLscale2kfoldgroup <- kfold(logNormalFitNULLscale2,
                                          folds = kfoldID,
                                          chains = 1,
                                          save_fits = TRUE)
logNormalFitNULLscale2kfoldPreds <- kfold_predict(logNormalFitNULLscale2kfoldgroup)
#logNormalFitNULLscale2kfoldPreds <- kfold_predict(logNormalFitNULLscale2kfold)
logNormalFitNULLscale2kfoldPredsDat <- logNormalFitNULLscale2kfoldPreds$yrep
logNormalFitNULLscale2kfoldPredsMean <- colMeans(logNormalFitNULLscale2kfoldPredsDat)
logNormalFitNULLscale2kfoldPredsMed <- apply(logNormalFitNULLscale2kfoldPredsDat, 2, function(x){quantile(x, 0.5)})
logNormalFitNULLscale2kfoldPredsLCB <- apply(logNormalFitNULLscale2kfoldPredsDat, 2, function(x){quantile(x, 0.025)})
logNormalFitNULLscale2kfoldPredsUCB <- apply(logNormalFitNULLscale2kfoldPredsDat, 2, function(x){quantile(x, 0.975)})

logNormalFitNULLscale2kfoldMetrics <- tibble(
  Fit = paste0("logNormalFitNULLscale2"),
  MAE_kfold = mean(abs(logNormalFitNULLscale2kfoldPredsMean - logNormalFitNULLscale2kfoldPreds$y)),
  MAD_kfold = mean(abs(logNormalFitNULLscale2kfoldPredsMed - logNormalFitNULLscale2kfoldPreds$y)),
  COV_kfold = mean(logNormalFitNULLscale2kfoldPredsLCB < logNormalFitNULLscale2kfoldPreds$y & logNormalFitNULLscale2kfoldPreds$y < logNormalFitNULLscale2kfoldPredsUCB)
)
logNormalFitNULLscale2kfoldMetrics

#### Compare Identity Models ----
print(logNormalFitNULL, digits = 4)
print(logNormalFitNULL2, digits = 4)
print(logNormalFitNULLscale, digits = 4)
print(logNormalFitNULLscale2, digits = 4)

logNormalFitNULLpredMetrics
logNormalFitNULL2predMetrics
logNormalFitNULLscalepredMetrics
logNormalFitNULLscale2predMetrics

### Log HWRF ----
#### Identity Link ----
logNormalFitNULLlog <- brm(
  bf(
    VMAX ~ log(HWRF)
  ),
  data = StormdataTrain8, 
  family = lognormal(link = "identity"),
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)

##### Plot ----
print(logNormalFitNULLlog, digits = 4)
logNormalFitNULLlogppcFit <- pp_check(logNormalFitNULLlog, ndraws = 100) + 
  labs(title = "logNormalFitNULLlog Fit PPC") +
  theme_bw()
logNormalFitNULLlogppcFit

##### LOO ----
logNormalFitNULLlogloo <- loo(logNormalFitNULLlog)

##### Prediction ----
## Fitted
logNormalFitNULLlogFit <- posterior_predict(logNormalFitNULLlog)
logNormalFitNULLlogFitMean <- colMeans(logNormalFitNULLlogFit)
logNormalFitNULLlogFitMed <- apply(logNormalFitNULLlogFit, 2, function(x){quantile(x, 0.5)})
logNormalFitNULLlogFitLCB <- apply(logNormalFitNULLlogFit, 2, function(x){quantile(x, 0.025)})
logNormalFitNULLlogFitUCB <- apply(logNormalFitNULLlogFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
logNormalFitNULLlogPreds <- posterior_predict(logNormalFitNULLlog, 
                                              newdata = StormdataTest8,
                                              allow_new_levels = TRUE)
logNormalFitNULLlogPredsMean <- colMeans(logNormalFitNULLlogPreds)
logNormalFitNULLlogPredsMed <- apply(logNormalFitNULLlogPreds, 2, function(x){quantile(x, 0.5)})
logNormalFitNULLlogPredsLCB <- apply(logNormalFitNULLlogPreds, 2, function(x){quantile(x, 0.025)})
logNormalFitNULLlogPredsUCB <- apply(logNormalFitNULLlogPreds, 2, function(x){quantile(x, 0.975)})

logNormalFitNULLlogpredMetrics <- tibble(
  Fit = "logNormalFitNULLlog",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(logNormalFitNULLlogFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(logNormalFitNULLlogFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < logNormalFitNULLlogFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest3$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(logNormalFitNULLlogPredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(logNormalFitNULLlogPredsMed - Actual_Yvec)),
  COV_pred = mean(logNormalFitNULLlogPredsLCB < Actual_Yvec & Actual_Yvec < logNormalFitNULLlogPredsUCB)
)
logNormalFitNULLlogpredMetrics

##### PPC ----
###### Quantile 2.5 
logNormalFitNULLlogLCBsims <- apply(logNormalFitNULLlogFit, 
                                    MARGIN = 1,
                                    function(x){
                                      quantile(x, 0.025)
                                    })
logNormalFitNULLlogLCBpvalueVec <- logNormalFitNULLlogLCBsims < quantile(StormdataTrain3$VMAX, 0.025)
logNormalFitNULLlogLCBpvalue <- sum(logNormalFitNULLlogLCBpvalueVec)
logNormalFitNULLlogLCBpvalue <- round(logNormalFitNULLlogLCBpvalue/4000, 3)
logNormalFitNULLlogLCBpvalue <- min(logNormalFitNULLlogLCBpvalue, 1 - logNormalFitNULLlogLCBpvalue)

logNormalFitNULLlog_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLlogFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", logNormalFitNULLlogLCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLlog_ppcLCB

###### Quantile 97.5 
logNormalFitNULLlogUCBsims <- apply(logNormalFitNULLlogFit, 
                                    MARGIN = 1,
                                    function(x){
                                      quantile(x, 0.975)
                                    })
logNormalFitNULLlogUCBpvalueVec <- logNormalFitNULLlogUCBsims < quantile(StormdataTrain3$VMAX, 0.975)
logNormalFitNULLlogUCBpvalue <- as.numeric(sum(logNormalFitNULLlogUCBpvalueVec))
logNormalFitNULLlogUCBpvalue <- round(logNormalFitNULLlogUCBpvalue/4000, 3)
logNormalFitNULLlogUCBpvalue <- min(logNormalFitNULLlogUCBpvalue, 1 - logNormalFitNULLlogUCBpvalue)

logNormalFitNULLlog_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLlogFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", logNormalFitNULLlogUCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLlog_ppcUCB

###### Mean 
logNormalFitNULLlogMEANsims <- apply(logNormalFitNULLlogFit, 
                                     MARGIN = 1,
                                     function(x){
                                       mean(x)
                                     })
logNormalFitNULLlogMEANpvalueVec <- logNormalFitNULLlogMEANsims < mean(StormdataTrain3$VMAX)
logNormalFitNULLlogMEANpvalue <- sum(logNormalFitNULLlogMEANpvalueVec)
logNormalFitNULLlogMEANpvalue <- round(logNormalFitNULLlogMEANpvalue/4000, 3)
logNormalFitNULLlogMEANpvalue <- min(logNormalFitNULLlogMEANpvalue, 1 - logNormalFitNULLlogMEANpvalue)

logNormalFitNULLlog_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLlogFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", logNormalFitNULLlogMEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLlog_ppcMEAN

###### Med 
logNormalFitNULLlogMEDsims <- apply(logNormalFitNULLlogFit, 
                                    MARGIN = 1,
                                    function(x){
                                      quantile(x, 0.5)
                                    })
logNormalFitNULLlogMEDpvalueVec <- logNormalFitNULLlogMEDsims < quantile(StormdataTrain3$VMAX, 0.5)
logNormalFitNULLlogMEDpvalue <- sum(logNormalFitNULLlogMEDpvalueVec)
logNormalFitNULLlogMEDpvalue <- round(logNormalFitNULLlogMEDpvalue/4000, 3)
logNormalFitNULLlogMEDpvalue <- min(logNormalFitNULLlogMEDpvalue, 1 - logNormalFitNULLlogMEDpvalue)

logNormalFitNULLlog_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLlogFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", logNormalFitNULLlogMEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLlog_ppcMED

###### SD 
logNormalFitNULLlogSDsims <- apply(logNormalFitNULLlogFit, 
                                   MARGIN = 1,
                                   function(x){
                                     sd(x)
                                   })
logNormalFitNULLlogSDpvalueVec <- logNormalFitNULLlogSDsims < sd(StormdataTrain3$VMAX)
logNormalFitNULLlogSDpvalue <- sum(logNormalFitNULLlogSDpvalueVec)
logNormalFitNULLlogSDpvalue <- round(logNormalFitNULLlogSDpvalue/4000, 3)
logNormalFitNULLlogSDpvalue <- min(logNormalFitNULLlogSDpvalue, 1 - logNormalFitNULLlogSDpvalue)

logNormalFitNULLlog_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLlogFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", logNormalFitNULLlogSDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLlog_ppcSD

###### Range 
logNormalFitNULLlogRANGEsims <- apply(logNormalFitNULLlogFit, 
                                      MARGIN = 1,
                                      function(x){
                                        max(x)-min(x)
                                      })
logNormalFitNULLlogRANGEpvalueVec <- logNormalFitNULLlogRANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
logNormalFitNULLlogRANGEpvalue <- sum(logNormalFitNULLlogRANGEpvalueVec)
logNormalFitNULLlogRANGEpvalue <- round(logNormalFitNULLlogRANGEpvalue/4000, 3)
logNormalFitNULLlogRANGEpvalue <- min(logNormalFitNULLlogRANGEpvalue, 1 - logNormalFitNULLlogRANGEpvalue)

logNormalFitNULLlog_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLlogFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", logNormalFitNULLlogRANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLlog_ppcRANGE

##### Combined Plot ----
logNormalFitNULLlog_ppcComb <- 
  logNormalFitNULLlogppcFit /
  (logNormalFitNULLlog_ppcLCB | logNormalFitNULLlog_ppcMED | logNormalFitNULLlog_ppcUCB) /
  (logNormalFitNULLlog_ppcRANGE | logNormalFitNULLlog_ppcMEAN | logNormalFitNULLlog_ppcSD)
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
logNormalFitNULLlog_ppcComb

##### Bayes p-values ----
logNormalFitNULLlogpvalues <- tibble(
  Fit = paste0("logNormalFitNULLlog"),
  LCB = logNormalFitNULLlogLCBpvalue,
  Median = logNormalFitNULLlogMEDpvalue,
  UCB = logNormalFitNULLlogUCBpvalue,
  Range = logNormalFitNULLlogRANGEpvalue,
  Mean = logNormalFitNULLlogMEANpvalue,
  SD = logNormalFitNULLlogSDpvalue
)
logNormalFitNULLlogpvalues

##### CV ----
kfoldID <- kfold_split_grouped(K = 5, StormdataTrain8$StormID)
logNormalFitNULLlogkfoldgroup <- kfold(logNormalFitNULLlog,
                                       folds = kfoldID,
                                       chains = 1,
                                       save_fits = TRUE)
logNormalFitNULLlogkfoldPreds <- kfold_predict(logNormalFitNULLlogkfoldgroup)
#logNormalFitNULLlogkfoldPreds <- kfold_predict(logNormalFitNULLlogkfold)
logNormalFitNULLlogkfoldPredsDat <- logNormalFitNULLlogkfoldPreds$yrep
logNormalFitNULLlogkfoldPredsMean <- colMeans(logNormalFitNULLlogkfoldPredsDat)
logNormalFitNULLlogkfoldPredsMed <- apply(logNormalFitNULLlogkfoldPredsDat, 2, function(x){quantile(x, 0.5)})
logNormalFitNULLlogkfoldPredsLCB <- apply(logNormalFitNULLlogkfoldPredsDat, 2, function(x){quantile(x, 0.025)})
logNormalFitNULLlogkfoldPredsUCB <- apply(logNormalFitNULLlogkfoldPredsDat, 2, function(x){quantile(x, 0.975)})

logNormalFitNULLlogkfoldMetrics <- tibble(
  Fit = paste0("logNormalFitNULLlog"),
  MAE_kfold = mean(abs(logNormalFitNULLlogkfoldPredsMean - logNormalFitNULLlogkfoldPreds$y)),
  MAD_kfold = mean(abs(logNormalFitNULLlogkfoldPredsMed - logNormalFitNULLlogkfoldPreds$y)),
  COV_kfold = mean(logNormalFitNULLlogkfoldPredsLCB < logNormalFitNULLlogkfoldPreds$y & logNormalFitNULLlogkfoldPreds$y < logNormalFitNULLlogkfoldPredsUCB)
)
logNormalFitNULLlogkfoldMetrics

#### Identity Link ----
logNormalFitNULLlog3 <- brm(
  bf(VMAX ~ logHWRF + 
       (1|StormID)
  ),
  data = StormdataTrain9, 
  family = lognormal(link = "identity"),
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)

##### Plot ----
print(logNormalFitNULLlog3, digits = 4)
logNormalFitNULLlog3ppcFit <- pp_check(logNormalFitNULLlog3, ndraws = 100) + 
  labs(title = "logNormalFitNULLlog3 Fit PPC") +
  theme_bw()
logNormalFitNULLlog3ppcFit

##### LOO ----
logNormalFitNULLlog3loo <- loo(logNormalFitNULLlog3)

##### Prediction ----
## Fitted
logNormalFitNULLlog3Fit <- posterior_predict(logNormalFitNULLlog3)
logNormalFitNULLlog3FitMean <- colMeans(logNormalFitNULLlog3Fit)
logNormalFitNULLlog3FitMed <- apply(logNormalFitNULLlog3Fit, 2, function(x){quantile(x, 0.5)})
logNormalFitNULLlog3FitLCB <- apply(logNormalFitNULLlog3Fit, 2, function(x){quantile(x, 0.025)})
logNormalFitNULLlog3FitUCB <- apply(logNormalFitNULLlog3Fit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
logNormalFitNULLlog3Preds <- posterior_predict(logNormalFitNULLlog3, 
                                              newdata = StormdataTest9,
                                              allow_new_levels = TRUE)
logNormalFitNULLlog3PredsMean <- colMeans(logNormalFitNULLlog3Preds)
logNormalFitNULLlog3PredsMed <- apply(logNormalFitNULLlog3Preds, 2, function(x){quantile(x, 0.5)})
logNormalFitNULLlog3PredsLCB <- apply(logNormalFitNULLlog3Preds, 2, function(x){quantile(x, 0.025)})
logNormalFitNULLlog3PredsUCB <- apply(logNormalFitNULLlog3Preds, 2, function(x){quantile(x, 0.975)})

logNormalFitNULLlog3predMetrics <- tibble(
  Fit = "logNormalFitNULLlog3",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(logNormalFitNULLlog3FitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(logNormalFitNULLlog3FitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < logNormalFitNULLlog3FitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest3$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(logNormalFitNULLlog3PredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(logNormalFitNULLlog3PredsMed - Actual_Yvec)),
  COV_pred = mean(logNormalFitNULLlog3PredsLCB < Actual_Yvec & Actual_Yvec < logNormalFitNULLlog3PredsUCB)
)
logNormalFitNULLlog3predMetrics

##### PPC ----
###### Quantile 2.5 
logNormalFitNULLlog3LCBsims <- apply(logNormalFitNULLlog3Fit, 
                                    MARGIN = 1,
                                    function(x){
                                      quantile(x, 0.025)
                                    })
logNormalFitNULLlog3LCBpvalueVec <- logNormalFitNULLlog3LCBsims < quantile(StormdataTrain3$VMAX, 0.025)
logNormalFitNULLlog3LCBpvalue <- sum(logNormalFitNULLlog3LCBpvalueVec)
logNormalFitNULLlog3LCBpvalue <- round(logNormalFitNULLlog3LCBpvalue/4000, 3)
logNormalFitNULLlog3LCBpvalue <- min(logNormalFitNULLlog3LCBpvalue, 1 - logNormalFitNULLlog3LCBpvalue)

logNormalFitNULLlog3_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLlog3Fit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", logNormalFitNULLlog3LCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLlog3_ppcLCB

###### Quantile 97.5 
logNormalFitNULLlog3UCBsims <- apply(logNormalFitNULLlog3Fit, 
                                    MARGIN = 1,
                                    function(x){
                                      quantile(x, 0.975)
                                    })
logNormalFitNULLlog3UCBpvalueVec <- logNormalFitNULLlog3UCBsims < quantile(StormdataTrain3$VMAX, 0.975)
logNormalFitNULLlog3UCBpvalue <- as.numeric(sum(logNormalFitNULLlog3UCBpvalueVec))
logNormalFitNULLlog3UCBpvalue <- round(logNormalFitNULLlog3UCBpvalue/4000, 3)
logNormalFitNULLlog3UCBpvalue <- min(logNormalFitNULLlog3UCBpvalue, 1 - logNormalFitNULLlog3UCBpvalue)

logNormalFitNULLlog3_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLlog3Fit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", logNormalFitNULLlog3UCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLlog3_ppcUCB

###### Mean 
logNormalFitNULLlog3MEANsims <- apply(logNormalFitNULLlog3Fit, 
                                     MARGIN = 1,
                                     function(x){
                                       mean(x)
                                     })
logNormalFitNULLlog3MEANpvalueVec <- logNormalFitNULLlog3MEANsims < mean(StormdataTrain3$VMAX)
logNormalFitNULLlog3MEANpvalue <- sum(logNormalFitNULLlog3MEANpvalueVec)
logNormalFitNULLlog3MEANpvalue <- round(logNormalFitNULLlog3MEANpvalue/4000, 3)
logNormalFitNULLlog3MEANpvalue <- min(logNormalFitNULLlog3MEANpvalue, 1 - logNormalFitNULLlog3MEANpvalue)

logNormalFitNULLlog3_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLlog3Fit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", logNormalFitNULLlog3MEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLlog3_ppcMEAN

###### Med 
logNormalFitNULLlog3MEDsims <- apply(logNormalFitNULLlog3Fit, 
                                    MARGIN = 1,
                                    function(x){
                                      quantile(x, 0.5)
                                    })
logNormalFitNULLlog3MEDpvalueVec <- logNormalFitNULLlog3MEDsims < quantile(StormdataTrain3$VMAX, 0.5)
logNormalFitNULLlog3MEDpvalue <- sum(logNormalFitNULLlog3MEDpvalueVec)
logNormalFitNULLlog3MEDpvalue <- round(logNormalFitNULLlog3MEDpvalue/4000, 3)
logNormalFitNULLlog3MEDpvalue <- min(logNormalFitNULLlog3MEDpvalue, 1 - logNormalFitNULLlog3MEDpvalue)

logNormalFitNULLlog3_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLlog3Fit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", logNormalFitNULLlog3MEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLlog3_ppcMED

###### SD 
logNormalFitNULLlog3SDsims <- apply(logNormalFitNULLlog3Fit, 
                                   MARGIN = 1,
                                   function(x){
                                     sd(x)
                                   })
logNormalFitNULLlog3SDpvalueVec <- logNormalFitNULLlog3SDsims < sd(StormdataTrain3$VMAX)
logNormalFitNULLlog3SDpvalue <- sum(logNormalFitNULLlog3SDpvalueVec)
logNormalFitNULLlog3SDpvalue <- round(logNormalFitNULLlog3SDpvalue/4000, 3)
logNormalFitNULLlog3SDpvalue <- min(logNormalFitNULLlog3SDpvalue, 1 - logNormalFitNULLlog3SDpvalue)

logNormalFitNULLlog3_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLlog3Fit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", logNormalFitNULLlog3SDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLlog3_ppcSD

###### Range 
logNormalFitNULLlog3RANGEsims <- apply(logNormalFitNULLlog3Fit, 
                                      MARGIN = 1,
                                      function(x){
                                        max(x)-min(x)
                                      })
logNormalFitNULLlog3RANGEpvalueVec <- logNormalFitNULLlog3RANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
logNormalFitNULLlog3RANGEpvalue <- sum(logNormalFitNULLlog3RANGEpvalueVec)
logNormalFitNULLlog3RANGEpvalue <- round(logNormalFitNULLlog3RANGEpvalue/4000, 3)
logNormalFitNULLlog3RANGEpvalue <- min(logNormalFitNULLlog3RANGEpvalue, 1 - logNormalFitNULLlog3RANGEpvalue)

logNormalFitNULLlog3_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLlog3Fit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", logNormalFitNULLlog3RANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLlog3_ppcRANGE

##### Combined Plot ----
logNormalFitNULLlog3_ppcComb <- 
  logNormalFitNULLlog3ppcFit /
  (logNormalFitNULLlog3_ppcLCB | logNormalFitNULLlog3_ppcMED | logNormalFitNULLlog3_ppcUCB) /
  (logNormalFitNULLlog3_ppcRANGE | logNormalFitNULLlog3_ppcMEAN | logNormalFitNULLlog3_ppcSD)
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
logNormalFitNULLlog3_ppcComb

##### Bayes p-values ----
logNormalFitNULLlog3pvalues <- tibble(
  Fit = paste0("logNormalFitNULLlog3"),
  LCB = logNormalFitNULLlog3LCBpvalue,
  Median = logNormalFitNULLlog3MEDpvalue,
  UCB = logNormalFitNULLlog3UCBpvalue,
  Range = logNormalFitNULLlog3RANGEpvalue,
  Mean = logNormalFitNULLlog3MEANpvalue,
  SD = logNormalFitNULLlog3SDpvalue
)
logNormalFitNULLlog3pvalues

##### CV ----
kfoldID <- kfold_split_grouped(K = 5, StormdataTrain8$StormID)
logNormalFitNULLlog3kfoldgroup <- kfold(logNormalFitNULLlog3,
                                       folds = kfoldID,
                                       chains = 1,
                                       save_fits = TRUE)
logNormalFitNULLlog3kfoldPreds <- kfold_predict(logNormalFitNULLlog3kfoldgroup)
#logNormalFitNULLlog3kfoldPreds <- kfold_predict(logNormalFitNULLlog3kfold)
logNormalFitNULLlog3kfoldPredsDat <- logNormalFitNULLlog3kfoldPreds$yrep
logNormalFitNULLlog3kfoldPredsMean <- colMeans(logNormalFitNULLlog3kfoldPredsDat)
logNormalFitNULLlog3kfoldPredsMed <- apply(logNormalFitNULLlog3kfoldPredsDat, 2, function(x){quantile(x, 0.5)})
logNormalFitNULLlog3kfoldPredsLCB <- apply(logNormalFitNULLlog3kfoldPredsDat, 2, function(x){quantile(x, 0.025)})
logNormalFitNULLlog3kfoldPredsUCB <- apply(logNormalFitNULLlog3kfoldPredsDat, 2, function(x){quantile(x, 0.975)})

logNormalFitNULLlog3kfoldMetrics <- tibble(
  Fit = paste0("logNormalFitNULLlog3"),
  MAE_kfold = mean(abs(logNormalFitNULLlog3kfoldPredsMean - logNormalFitNULLlog3kfoldPreds$y)),
  MAD_kfold = mean(abs(logNormalFitNULLlog3kfoldPredsMed - logNormalFitNULLlog3kfoldPreds$y)),
  COV_kfold = mean(logNormalFitNULLlog3kfoldPredsLCB < logNormalFitNULLlog3kfoldPreds$y & logNormalFitNULLlog3kfoldPreds$y < logNormalFitNULLlog3kfoldPredsUCB)
)
logNormalFitNULLlog3kfoldMetrics

#### Identity Link Log Sigma ----
logNormalFitNULLlogWsigma <- brm(
  bf(
    VMAX ~ log(HWRF),
    sigma ~ 1
  ),
  data = StormdataTrain8, 
  family = lognormal(link = "identity", link_sigma = "log"),
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)

##### Plot ----
print(logNormalFitNULLlogWsigma, digits = 4)
logNormalFitNULLlogWsigmappcFit <- pp_check(logNormalFitNULLlogWsigma, ndraws = 100) + 
  labs(title = "logNormalFitNULLlogWsigma Fit PPC") +
  theme_bw()
logNormalFitNULLlogWsigmappcFit

##### LOO ----
logNormalFitNULLlogWsigmaloo <- loo(logNormalFitNULLlogWsigma)

##### Prediction ----
## Fitted
logNormalFitNULLlogWsigmaFit <- posterior_predict(logNormalFitNULLlogWsigma)
logNormalFitNULLlogWsigmaFitMean <- colMeans(logNormalFitNULLlogWsigmaFit)
logNormalFitNULLlogWsigmaFitMed <- apply(logNormalFitNULLlogWsigmaFit, 2, function(x){quantile(x, 0.5)})
logNormalFitNULLlogWsigmaFitLCB <- apply(logNormalFitNULLlogWsigmaFit, 2, function(x){quantile(x, 0.025)})
logNormalFitNULLlogWsigmaFitUCB <- apply(logNormalFitNULLlogWsigmaFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
logNormalFitNULLlogWsigmaPreds <- posterior_predict(logNormalFitNULLlogWsigma, 
                                              newdata = StormdataTest8,
                                              allow_new_levels = TRUE)
logNormalFitNULLlogWsigmaPredsMean <- colMeans(logNormalFitNULLlogWsigmaPreds)
logNormalFitNULLlogWsigmaPredsMed <- apply(logNormalFitNULLlogWsigmaPreds, 2, function(x){quantile(x, 0.5)})
logNormalFitNULLlogWsigmaPredsLCB <- apply(logNormalFitNULLlogWsigmaPreds, 2, function(x){quantile(x, 0.025)})
logNormalFitNULLlogWsigmaPredsUCB <- apply(logNormalFitNULLlogWsigmaPreds, 2, function(x){quantile(x, 0.975)})

logNormalFitNULLlogWsigmapredMetrics <- tibble(
  Fit = "logNormalFitNULLlogWsigma",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(logNormalFitNULLlogWsigmaFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(logNormalFitNULLlogWsigmaFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < logNormalFitNULLlogWsigmaFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest3$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(logNormalFitNULLlogWsigmaPredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(logNormalFitNULLlogWsigmaPredsMed - Actual_Yvec)),
  COV_pred = mean(logNormalFitNULLlogWsigmaPredsLCB < Actual_Yvec & Actual_Yvec < logNormalFitNULLlogWsigmaPredsUCB)
)
logNormalFitNULLlogWsigmapredMetrics

##### PPC ----
###### Quantile 2.5 
logNormalFitNULLlogWsigmaLCBsims <- apply(logNormalFitNULLlogWsigmaFit, 
                                    MARGIN = 1,
                                    function(x){
                                      quantile(x, 0.025)
                                    })
logNormalFitNULLlogWsigmaLCBpvalueVec <- logNormalFitNULLlogWsigmaLCBsims < quantile(StormdataTrain3$VMAX, 0.025)
logNormalFitNULLlogWsigmaLCBpvalue <- sum(logNormalFitNULLlogWsigmaLCBpvalueVec)
logNormalFitNULLlogWsigmaLCBpvalue <- round(logNormalFitNULLlogWsigmaLCBpvalue/4000, 3)
logNormalFitNULLlogWsigmaLCBpvalue <- min(logNormalFitNULLlogWsigmaLCBpvalue, 1 - logNormalFitNULLlogWsigmaLCBpvalue)

logNormalFitNULLlogWsigma_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLlogWsigmaFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", logNormalFitNULLlogWsigmaLCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLlogWsigma_ppcLCB

###### Quantile 97.5 
logNormalFitNULLlogWsigmaUCBsims <- apply(logNormalFitNULLlogWsigmaFit, 
                                    MARGIN = 1,
                                    function(x){
                                      quantile(x, 0.975)
                                    })
logNormalFitNULLlogWsigmaUCBpvalueVec <- logNormalFitNULLlogWsigmaUCBsims < quantile(StormdataTrain3$VMAX, 0.975)
logNormalFitNULLlogWsigmaUCBpvalue <- as.numeric(sum(logNormalFitNULLlogWsigmaUCBpvalueVec))
logNormalFitNULLlogWsigmaUCBpvalue <- round(logNormalFitNULLlogWsigmaUCBpvalue/4000, 3)
logNormalFitNULLlogWsigmaUCBpvalue <- min(logNormalFitNULLlogWsigmaUCBpvalue, 1 - logNormalFitNULLlogWsigmaUCBpvalue)

logNormalFitNULLlogWsigma_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLlogWsigmaFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", logNormalFitNULLlogWsigmaUCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLlogWsigma_ppcUCB

###### Mean 
logNormalFitNULLlogWsigmaMEANsims <- apply(logNormalFitNULLlogWsigmaFit, 
                                     MARGIN = 1,
                                     function(x){
                                       mean(x)
                                     })
logNormalFitNULLlogWsigmaMEANpvalueVec <- logNormalFitNULLlogWsigmaMEANsims < mean(StormdataTrain3$VMAX)
logNormalFitNULLlogWsigmaMEANpvalue <- sum(logNormalFitNULLlogWsigmaMEANpvalueVec)
logNormalFitNULLlogWsigmaMEANpvalue <- round(logNormalFitNULLlogWsigmaMEANpvalue/4000, 3)
logNormalFitNULLlogWsigmaMEANpvalue <- min(logNormalFitNULLlogWsigmaMEANpvalue, 1 - logNormalFitNULLlogWsigmaMEANpvalue)

logNormalFitNULLlogWsigma_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLlogWsigmaFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", logNormalFitNULLlogWsigmaMEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLlogWsigma_ppcMEAN

###### Med 
logNormalFitNULLlogWsigmaMEDsims <- apply(logNormalFitNULLlogWsigmaFit, 
                                    MARGIN = 1,
                                    function(x){
                                      quantile(x, 0.5)
                                    })
logNormalFitNULLlogWsigmaMEDpvalueVec <- logNormalFitNULLlogWsigmaMEDsims < quantile(StormdataTrain3$VMAX, 0.5)
logNormalFitNULLlogWsigmaMEDpvalue <- sum(logNormalFitNULLlogWsigmaMEDpvalueVec)
logNormalFitNULLlogWsigmaMEDpvalue <- round(logNormalFitNULLlogWsigmaMEDpvalue/4000, 3)
logNormalFitNULLlogWsigmaMEDpvalue <- min(logNormalFitNULLlogWsigmaMEDpvalue, 1 - logNormalFitNULLlogWsigmaMEDpvalue)

logNormalFitNULLlogWsigma_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLlogWsigmaFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", logNormalFitNULLlogWsigmaMEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLlogWsigma_ppcMED

###### SD 
logNormalFitNULLlogWsigmaSDsims <- apply(logNormalFitNULLlogWsigmaFit, 
                                   MARGIN = 1,
                                   function(x){
                                     sd(x)
                                   })
logNormalFitNULLlogWsigmaSDpvalueVec <- logNormalFitNULLlogWsigmaSDsims < sd(StormdataTrain3$VMAX)
logNormalFitNULLlogWsigmaSDpvalue <- sum(logNormalFitNULLlogWsigmaSDpvalueVec)
logNormalFitNULLlogWsigmaSDpvalue <- round(logNormalFitNULLlogWsigmaSDpvalue/4000, 3)
logNormalFitNULLlogWsigmaSDpvalue <- min(logNormalFitNULLlogWsigmaSDpvalue, 1 - logNormalFitNULLlogWsigmaSDpvalue)

logNormalFitNULLlogWsigma_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLlogWsigmaFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", logNormalFitNULLlogWsigmaSDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLlogWsigma_ppcSD

###### Range 
logNormalFitNULLlogWsigmaRANGEsims <- apply(logNormalFitNULLlogWsigmaFit, 
                                      MARGIN = 1,
                                      function(x){
                                        max(x)-min(x)
                                      })
logNormalFitNULLlogWsigmaRANGEpvalueVec <- logNormalFitNULLlogWsigmaRANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
logNormalFitNULLlogWsigmaRANGEpvalue <- sum(logNormalFitNULLlogWsigmaRANGEpvalueVec)
logNormalFitNULLlogWsigmaRANGEpvalue <- round(logNormalFitNULLlogWsigmaRANGEpvalue/4000, 3)
logNormalFitNULLlogWsigmaRANGEpvalue <- min(logNormalFitNULLlogWsigmaRANGEpvalue, 1 - logNormalFitNULLlogWsigmaRANGEpvalue)

logNormalFitNULLlogWsigma_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLlogWsigmaFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", logNormalFitNULLlogWsigmaRANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLlogWsigma_ppcRANGE

##### Combined Plot ----
logNormalFitNULLlogWsigma_ppcComb <- 
  logNormalFitNULLlogWsigmappcFit /
  (logNormalFitNULLlogWsigma_ppcLCB | logNormalFitNULLlogWsigma_ppcMED | logNormalFitNULLlogWsigma_ppcUCB) /
  (logNormalFitNULLlogWsigma_ppcRANGE | logNormalFitNULLlogWsigma_ppcMEAN | logNormalFitNULLlogWsigma_ppcSD)
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
logNormalFitNULLlogWsigma_ppcComb

##### Bayes p-values ----
logNormalFitNULLlogWsigmapvalues <- tibble(
  Fit = paste0("logNormalFitNULLlogWsigma"),
  LCB = logNormalFitNULLlogWsigmaLCBpvalue,
  Median = logNormalFitNULLlogWsigmaMEDpvalue,
  UCB = logNormalFitNULLlogWsigmaUCBpvalue,
  Range = logNormalFitNULLlogWsigmaRANGEpvalue,
  Mean = logNormalFitNULLlogWsigmaMEANpvalue,
  SD = logNormalFitNULLlogWsigmaSDpvalue
)
logNormalFitNULLlogWsigmapvalues

##### CV ----
kfoldID <- kfold_split_grouped(K = 5, StormdataTrain8$StormID)
logNormalFitNULLlogWsigmakfoldgroup <- kfold(logNormalFitNULLlogWsigma,
                                       folds = kfoldID,
                                       chains = 1,
                                       save_fits = TRUE)
logNormalFitNULLlogWsigmakfoldPreds <- kfold_predict(logNormalFitNULLlogWsigmakfoldgroup)
#logNormalFitNULLlogWsigmakfoldPreds <- kfold_predict(logNormalFitNULLlogWsigmakfold)
logNormalFitNULLlogWsigmakfoldPredsDat <- logNormalFitNULLlogWsigmakfoldPreds$yrep
logNormalFitNULLlogWsigmakfoldPredsMean <- colMeans(logNormalFitNULLlogWsigmakfoldPredsDat)
logNormalFitNULLlogWsigmakfoldPredsMed <- apply(logNormalFitNULLlogWsigmakfoldPredsDat, 2, function(x){quantile(x, 0.5)})
logNormalFitNULLlogWsigmakfoldPredsLCB <- apply(logNormalFitNULLlogWsigmakfoldPredsDat, 2, function(x){quantile(x, 0.025)})
logNormalFitNULLlogWsigmakfoldPredsUCB <- apply(logNormalFitNULLlogWsigmakfoldPredsDat, 2, function(x){quantile(x, 0.975)})

logNormalFitNULLlogWsigmakfoldMetrics <- tibble(
  Fit = paste0("logNormalFitNULLlogWsigma"),
  MAE_kfold = mean(abs(logNormalFitNULLlogWsigmakfoldPredsMean - logNormalFitNULLlogWsigmakfoldPreds$y)),
  MAD_kfold = mean(abs(logNormalFitNULLlogWsigmakfoldPredsMed - logNormalFitNULLlogWsigmakfoldPreds$y)),
  COV_kfold = mean(logNormalFitNULLlogWsigmakfoldPredsLCB < logNormalFitNULLlogWsigmakfoldPreds$y & logNormalFitNULLlogWsigmakfoldPreds$y < logNormalFitNULLlogWsigmakfoldPredsUCB)
)
logNormalFitNULLlogWsigmakfoldMetrics

#### Identity Link No Int ----
logNormalFitNULLlog2 <- brm(
  bf(
    VMAX ~ 0 + Intercept + log(HWRF)
  ),
  data = StormdataTrain8, 
  family = lognormal(link = "identity"),
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)

##### Plot ----
print(logNormalFitNULLlog2, digits = 4)
logNormalFitNULLlog2ppcFit <- pp_check(logNormalFitNULLlog2, ndraws = 100) + 
  labs(title = "logNormalFitNULLlog2 Fit PPC") +
  theme_bw()
logNormalFitNULLlog2ppcFit

##### LOO ----
logNormalFitNULLlog2loo <- loo(logNormalFitNULLlog2)

##### Prediction ----
## Fitted
logNormalFitNULLlog2Fit <- posterior_predict(logNormalFitNULLlog2)
logNormalFitNULLlog2FitMean <- colMeans(logNormalFitNULLlog2Fit)
logNormalFitNULLlog2FitMed <- apply(logNormalFitNULLlog2Fit, 2, function(x){quantile(x, 0.5)})
logNormalFitNULLlog2FitLCB <- apply(logNormalFitNULLlog2Fit, 2, function(x){quantile(x, 0.025)})
logNormalFitNULLlog2FitUCB <- apply(logNormalFitNULLlog2Fit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
logNormalFitNULLlog2Preds <- posterior_predict(logNormalFitNULLlog2, 
                                               newdata = StormdataTest8,
                                               allow_new_levels = TRUE)
logNormalFitNULLlog2PredsMean <- colMeans(logNormalFitNULLlog2Preds)
logNormalFitNULLlog2PredsMed <- apply(logNormalFitNULLlog2Preds, 2, function(x){quantile(x, 0.5)})
logNormalFitNULLlog2PredsLCB <- apply(logNormalFitNULLlog2Preds, 2, function(x){quantile(x, 0.025)})
logNormalFitNULLlog2PredsUCB <- apply(logNormalFitNULLlog2Preds, 2, function(x){quantile(x, 0.975)})

logNormalFitNULLlog2predMetrics <- tibble(
  Fit = "logNormalFitNULLlog2",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(logNormalFitNULLlog2FitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(logNormalFitNULLlog2FitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < logNormalFitNULLlog2FitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest3$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(logNormalFitNULLlog2PredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(logNormalFitNULLlog2PredsMed - Actual_Yvec)),
  COV_pred = mean(logNormalFitNULLlog2PredsLCB < Actual_Yvec & Actual_Yvec < logNormalFitNULLlog2PredsUCB)
)
logNormalFitNULLlog2predMetrics

##### PPC ----
###### Quantile 2.5 
logNormalFitNULLlog2LCBsims <- apply(logNormalFitNULLlog2Fit, 
                                     MARGIN = 1,
                                     function(x){
                                       quantile(x, 0.025)
                                     })
logNormalFitNULLlog2LCBpvalueVec <- logNormalFitNULLlog2LCBsims < quantile(StormdataTrain3$VMAX, 0.025)
logNormalFitNULLlog2LCBpvalue <- sum(logNormalFitNULLlog2LCBpvalueVec)
logNormalFitNULLlog2LCBpvalue <- round(logNormalFitNULLlog2LCBpvalue/4000, 3)
logNormalFitNULLlog2LCBpvalue <- min(logNormalFitNULLlog2LCBpvalue, 1 - logNormalFitNULLlog2LCBpvalue)

logNormalFitNULLlog2_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLlog2Fit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", logNormalFitNULLlog2LCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLlog2_ppcLCB

###### Quantile 97.5 
logNormalFitNULLlog2UCBsims <- apply(logNormalFitNULLlog2Fit, 
                                     MARGIN = 1,
                                     function(x){
                                       quantile(x, 0.975)
                                     })
logNormalFitNULLlog2UCBpvalueVec <- logNormalFitNULLlog2UCBsims < quantile(StormdataTrain3$VMAX, 0.975)
logNormalFitNULLlog2UCBpvalue <- as.numeric(sum(logNormalFitNULLlog2UCBpvalueVec))
logNormalFitNULLlog2UCBpvalue <- round(logNormalFitNULLlog2UCBpvalue/4000, 3)
logNormalFitNULLlog2UCBpvalue <- min(logNormalFitNULLlog2UCBpvalue, 1 - logNormalFitNULLlog2UCBpvalue)

logNormalFitNULLlog2_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLlog2Fit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", logNormalFitNULLlog2UCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLlog2_ppcUCB

###### Mean 
logNormalFitNULLlog2MEANsims <- apply(logNormalFitNULLlog2Fit, 
                                      MARGIN = 1,
                                      function(x){
                                        mean(x)
                                      })
logNormalFitNULLlog2MEANpvalueVec <- logNormalFitNULLlog2MEANsims < mean(StormdataTrain3$VMAX)
logNormalFitNULLlog2MEANpvalue <- sum(logNormalFitNULLlog2MEANpvalueVec)
logNormalFitNULLlog2MEANpvalue <- round(logNormalFitNULLlog2MEANpvalue/4000, 3)
logNormalFitNULLlog2MEANpvalue <- min(logNormalFitNULLlog2MEANpvalue, 1 - logNormalFitNULLlog2MEANpvalue)

logNormalFitNULLlog2_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLlog2Fit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", logNormalFitNULLlog2MEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLlog2_ppcMEAN

###### Med 
logNormalFitNULLlog2MEDsims <- apply(logNormalFitNULLlog2Fit, 
                                     MARGIN = 1,
                                     function(x){
                                       quantile(x, 0.5)
                                     })
logNormalFitNULLlog2MEDpvalueVec <- logNormalFitNULLlog2MEDsims < quantile(StormdataTrain3$VMAX, 0.5)
logNormalFitNULLlog2MEDpvalue <- sum(logNormalFitNULLlog2MEDpvalueVec)
logNormalFitNULLlog2MEDpvalue <- round(logNormalFitNULLlog2MEDpvalue/4000, 3)
logNormalFitNULLlog2MEDpvalue <- min(logNormalFitNULLlog2MEDpvalue, 1 - logNormalFitNULLlog2MEDpvalue)

logNormalFitNULLlog2_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLlog2Fit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", logNormalFitNULLlog2MEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLlog2_ppcMED

###### SD 
logNormalFitNULLlog2SDsims <- apply(logNormalFitNULLlog2Fit, 
                                    MARGIN = 1,
                                    function(x){
                                      sd(x)
                                    })
logNormalFitNULLlog2SDpvalueVec <- logNormalFitNULLlog2SDsims < sd(StormdataTrain3$VMAX)
logNormalFitNULLlog2SDpvalue <- sum(logNormalFitNULLlog2SDpvalueVec)
logNormalFitNULLlog2SDpvalue <- round(logNormalFitNULLlog2SDpvalue/4000, 3)
logNormalFitNULLlog2SDpvalue <- min(logNormalFitNULLlog2SDpvalue, 1 - logNormalFitNULLlog2SDpvalue)

logNormalFitNULLlog2_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLlog2Fit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", logNormalFitNULLlog2SDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLlog2_ppcSD

###### Range 
logNormalFitNULLlog2RANGEsims <- apply(logNormalFitNULLlog2Fit, 
                                       MARGIN = 1,
                                       function(x){
                                         max(x)-min(x)
                                       })
logNormalFitNULLlog2RANGEpvalueVec <- logNormalFitNULLlog2RANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
logNormalFitNULLlog2RANGEpvalue <- sum(logNormalFitNULLlog2RANGEpvalueVec)
logNormalFitNULLlog2RANGEpvalue <- round(logNormalFitNULLlog2RANGEpvalue/4000, 3)
logNormalFitNULLlog2RANGEpvalue <- min(logNormalFitNULLlog2RANGEpvalue, 1 - logNormalFitNULLlog2RANGEpvalue)

logNormalFitNULLlog2_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLlog2Fit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", logNormalFitNULLlog2RANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLlog2_ppcRANGE

##### Combined Plot ----
logNormalFitNULLlog2_ppcComb <- 
  logNormalFitNULLlog2ppcFit /
  (logNormalFitNULLlog2_ppcLCB | logNormalFitNULLlog2_ppcMED | logNormalFitNULLlog2_ppcUCB) /
  (logNormalFitNULLlog2_ppcRANGE | logNormalFitNULLlog2_ppcMEAN | logNormalFitNULLlog2_ppcSD)
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
logNormalFitNULLlog2_ppcComb

##### Bayes p-values ----
logNormalFitNULLlog2pvalues <- tibble(
  Fit = paste0("logNormalFitNULLlog2"),
  LCB = logNormalFitNULLlog2LCBpvalue,
  Median = logNormalFitNULLlog2MEDpvalue,
  UCB = logNormalFitNULLlog2UCBpvalue,
  Range = logNormalFitNULLlog2RANGEpvalue,
  Mean = logNormalFitNULLlog2MEANpvalue,
  SD = logNormalFitNULLlog2SDpvalue
)
logNormalFitNULLlog2pvalues

##### CV ----
kfoldID <- kfold_split_grouped(K = 5, StormdataTrain8$StormID)
logNormalFitNULLlog2kfoldgroup <- kfold(logNormalFitNULLlog2,
                                        folds = kfoldID,
                                        chains = 1,
                                        save_fits = TRUE)
logNormalFitNULLlog2kfoldPreds <- kfold_predict(logNormalFitNULLlog2kfoldgroup)
#logNormalFitNULLlog2kfoldPreds <- kfold_predict(logNormalFitNULLlog2kfold)
logNormalFitNULLlog2kfoldPredsDat <- logNormalFitNULLlog2kfoldPreds$yrep
logNormalFitNULLlog2kfoldPredsMean <- colMeans(logNormalFitNULLlog2kfoldPredsDat)
logNormalFitNULLlog2kfoldPredsMed <- apply(logNormalFitNULLlog2kfoldPredsDat, 2, function(x){quantile(x, 0.5)})
logNormalFitNULLlog2kfoldPredsLCB <- apply(logNormalFitNULLlog2kfoldPredsDat, 2, function(x){quantile(x, 0.025)})
logNormalFitNULLlog2kfoldPredsUCB <- apply(logNormalFitNULLlog2kfoldPredsDat, 2, function(x){quantile(x, 0.975)})

logNormalFitNULLlog2kfoldMetrics <- tibble(
  Fit = paste0("logNormalFitNULLlog2"),
  MAE_kfold = mean(abs(logNormalFitNULLlog2kfoldPredsMean - logNormalFitNULLlog2kfoldPreds$y)),
  MAD_kfold = mean(abs(logNormalFitNULLlog2kfoldPredsMed - logNormalFitNULLlog2kfoldPreds$y)),
  COV_kfold = mean(logNormalFitNULLlog2kfoldPredsLCB < logNormalFitNULLlog2kfoldPreds$y & logNormalFitNULLlog2kfoldPreds$y < logNormalFitNULLlog2kfoldPredsUCB)
)
logNormalFitNULLlog2kfoldMetrics

#### Identity Link Log Sigma ----
logNormalFitNULLlog2Wsigma <- brm(
  bf(
    VMAX ~ 0 + Intercept + log(HWRF),
    sigma ~ 1
  ),
  data = StormdataTrain8, 
  family = lognormal(link = "identity", link_sigma = "log"),
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)

##### Plot ----
print(logNormalFitNULLlog2Wsigma, digits = 4)
logNormalFitNULLlog2WsigmappcFit <- pp_check(logNormalFitNULLlog2Wsigma, ndraws = 100) + 
  labs(title = "logNormalFitNULLlog2Wsigma Fit PPC") +
  theme_bw()
logNormalFitNULLlog2WsigmappcFit

##### LOO ----
logNormalFitNULLlog2Wsigmaloo <- loo(logNormalFitNULLlog2Wsigma)

##### Prediction ----
## Fitted
logNormalFitNULLlog2WsigmaFit <- posterior_predict(logNormalFitNULLlog2Wsigma)
logNormalFitNULLlog2WsigmaFitMean <- colMeans(logNormalFitNULLlog2WsigmaFit)
logNormalFitNULLlog2WsigmaFitMed <- apply(logNormalFitNULLlog2WsigmaFit, 2, function(x){quantile(x, 0.5)})
logNormalFitNULLlog2WsigmaFitLCB <- apply(logNormalFitNULLlog2WsigmaFit, 2, function(x){quantile(x, 0.025)})
logNormalFitNULLlog2WsigmaFitUCB <- apply(logNormalFitNULLlog2WsigmaFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
logNormalFitNULLlog2WsigmaPreds <- posterior_predict(logNormalFitNULLlog2Wsigma, 
                                                    newdata = StormdataTest8,
                                                    allow_new_levels = TRUE)
logNormalFitNULLlog2WsigmaPredsMean <- colMeans(logNormalFitNULLlog2WsigmaPreds)
logNormalFitNULLlog2WsigmaPredsMed <- apply(logNormalFitNULLlog2WsigmaPreds, 2, function(x){quantile(x, 0.5)})
logNormalFitNULLlog2WsigmaPredsLCB <- apply(logNormalFitNULLlog2WsigmaPreds, 2, function(x){quantile(x, 0.025)})
logNormalFitNULLlog2WsigmaPredsUCB <- apply(logNormalFitNULLlog2WsigmaPreds, 2, function(x){quantile(x, 0.975)})

logNormalFitNULLlog2WsigmapredMetrics <- tibble(
  Fit = "logNormalFitNULLlog2Wsigma",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(logNormalFitNULLlog2WsigmaFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(logNormalFitNULLlog2WsigmaFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < logNormalFitNULLlog2WsigmaFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest3$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(logNormalFitNULLlog2WsigmaPredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(logNormalFitNULLlog2WsigmaPredsMed - Actual_Yvec)),
  COV_pred = mean(logNormalFitNULLlog2WsigmaPredsLCB < Actual_Yvec & Actual_Yvec < logNormalFitNULLlog2WsigmaPredsUCB)
)
logNormalFitNULLlog2WsigmapredMetrics

##### PPC ----
###### Quantile 2.5 
logNormalFitNULLlog2WsigmaLCBsims <- apply(logNormalFitNULLlog2WsigmaFit, 
                                          MARGIN = 1,
                                          function(x){
                                            quantile(x, 0.025)
                                          })
logNormalFitNULLlog2WsigmaLCBpvalueVec <- logNormalFitNULLlog2WsigmaLCBsims < quantile(StormdataTrain3$VMAX, 0.025)
logNormalFitNULLlog2WsigmaLCBpvalue <- sum(logNormalFitNULLlog2WsigmaLCBpvalueVec)
logNormalFitNULLlog2WsigmaLCBpvalue <- round(logNormalFitNULLlog2WsigmaLCBpvalue/4000, 3)
logNormalFitNULLlog2WsigmaLCBpvalue <- min(logNormalFitNULLlog2WsigmaLCBpvalue, 1 - logNormalFitNULLlog2WsigmaLCBpvalue)

logNormalFitNULLlog2Wsigma_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLlog2WsigmaFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", logNormalFitNULLlog2WsigmaLCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLlog2Wsigma_ppcLCB

###### Quantile 97.5 
logNormalFitNULLlog2WsigmaUCBsims <- apply(logNormalFitNULLlog2WsigmaFit, 
                                          MARGIN = 1,
                                          function(x){
                                            quantile(x, 0.975)
                                          })
logNormalFitNULLlog2WsigmaUCBpvalueVec <- logNormalFitNULLlog2WsigmaUCBsims < quantile(StormdataTrain3$VMAX, 0.975)
logNormalFitNULLlog2WsigmaUCBpvalue <- as.numeric(sum(logNormalFitNULLlog2WsigmaUCBpvalueVec))
logNormalFitNULLlog2WsigmaUCBpvalue <- round(logNormalFitNULLlog2WsigmaUCBpvalue/4000, 3)
logNormalFitNULLlog2WsigmaUCBpvalue <- min(logNormalFitNULLlog2WsigmaUCBpvalue, 1 - logNormalFitNULLlog2WsigmaUCBpvalue)

logNormalFitNULLlog2Wsigma_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLlog2WsigmaFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", logNormalFitNULLlog2WsigmaUCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLlog2Wsigma_ppcUCB

###### Mean 
logNormalFitNULLlog2WsigmaMEANsims <- apply(logNormalFitNULLlog2WsigmaFit, 
                                           MARGIN = 1,
                                           function(x){
                                             mean(x)
                                           })
logNormalFitNULLlog2WsigmaMEANpvalueVec <- logNormalFitNULLlog2WsigmaMEANsims < mean(StormdataTrain3$VMAX)
logNormalFitNULLlog2WsigmaMEANpvalue <- sum(logNormalFitNULLlog2WsigmaMEANpvalueVec)
logNormalFitNULLlog2WsigmaMEANpvalue <- round(logNormalFitNULLlog2WsigmaMEANpvalue/4000, 3)
logNormalFitNULLlog2WsigmaMEANpvalue <- min(logNormalFitNULLlog2WsigmaMEANpvalue, 1 - logNormalFitNULLlog2WsigmaMEANpvalue)

logNormalFitNULLlog2Wsigma_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLlog2WsigmaFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", logNormalFitNULLlog2WsigmaMEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLlog2Wsigma_ppcMEAN

###### Med 
logNormalFitNULLlog2WsigmaMEDsims <- apply(logNormalFitNULLlog2WsigmaFit, 
                                          MARGIN = 1,
                                          function(x){
                                            quantile(x, 0.5)
                                          })
logNormalFitNULLlog2WsigmaMEDpvalueVec <- logNormalFitNULLlog2WsigmaMEDsims < quantile(StormdataTrain3$VMAX, 0.5)
logNormalFitNULLlog2WsigmaMEDpvalue <- sum(logNormalFitNULLlog2WsigmaMEDpvalueVec)
logNormalFitNULLlog2WsigmaMEDpvalue <- round(logNormalFitNULLlog2WsigmaMEDpvalue/4000, 3)
logNormalFitNULLlog2WsigmaMEDpvalue <- min(logNormalFitNULLlog2WsigmaMEDpvalue, 1 - logNormalFitNULLlog2WsigmaMEDpvalue)

logNormalFitNULLlog2Wsigma_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLlog2WsigmaFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", logNormalFitNULLlog2WsigmaMEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLlog2Wsigma_ppcMED

###### SD 
logNormalFitNULLlog2WsigmaSDsims <- apply(logNormalFitNULLlog2WsigmaFit, 
                                         MARGIN = 1,
                                         function(x){
                                           sd(x)
                                         })
logNormalFitNULLlog2WsigmaSDpvalueVec <- logNormalFitNULLlog2WsigmaSDsims < sd(StormdataTrain3$VMAX)
logNormalFitNULLlog2WsigmaSDpvalue <- sum(logNormalFitNULLlog2WsigmaSDpvalueVec)
logNormalFitNULLlog2WsigmaSDpvalue <- round(logNormalFitNULLlog2WsigmaSDpvalue/4000, 3)
logNormalFitNULLlog2WsigmaSDpvalue <- min(logNormalFitNULLlog2WsigmaSDpvalue, 1 - logNormalFitNULLlog2WsigmaSDpvalue)

logNormalFitNULLlog2Wsigma_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLlog2WsigmaFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", logNormalFitNULLlog2WsigmaSDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLlog2Wsigma_ppcSD

###### Range 
logNormalFitNULLlog2WsigmaRANGEsims <- apply(logNormalFitNULLlog2WsigmaFit, 
                                            MARGIN = 1,
                                            function(x){
                                              max(x)-min(x)
                                            })
logNormalFitNULLlog2WsigmaRANGEpvalueVec <- logNormalFitNULLlog2WsigmaRANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
logNormalFitNULLlog2WsigmaRANGEpvalue <- sum(logNormalFitNULLlog2WsigmaRANGEpvalueVec)
logNormalFitNULLlog2WsigmaRANGEpvalue <- round(logNormalFitNULLlog2WsigmaRANGEpvalue/4000, 3)
logNormalFitNULLlog2WsigmaRANGEpvalue <- min(logNormalFitNULLlog2WsigmaRANGEpvalue, 1 - logNormalFitNULLlog2WsigmaRANGEpvalue)

logNormalFitNULLlog2Wsigma_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitNULLlog2WsigmaFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", logNormalFitNULLlog2WsigmaRANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULLlog2Wsigma_ppcRANGE

##### Combined Plot ----
logNormalFitNULLlog2Wsigma_ppcComb <- 
  logNormalFitNULLlog2WsigmappcFit /
  (logNormalFitNULLlog2Wsigma_ppcLCB | logNormalFitNULLlog2Wsigma_ppcMED | logNormalFitNULLlog2Wsigma_ppcUCB) /
  (logNormalFitNULLlog2Wsigma_ppcRANGE | logNormalFitNULLlog2Wsigma_ppcMEAN | logNormalFitNULLlog2Wsigma_ppcSD)
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
logNormalFitNULLlog2Wsigma_ppcComb

##### Bayes p-values ----
logNormalFitNULLlog2Wsigmapvalues <- tibble(
  Fit = paste0("logNormalFitNULLlog2Wsigma"),
  LCB = logNormalFitNULLlog2WsigmaLCBpvalue,
  Median = logNormalFitNULLlog2WsigmaMEDpvalue,
  UCB = logNormalFitNULLlog2WsigmaUCBpvalue,
  Range = logNormalFitNULLlog2WsigmaRANGEpvalue,
  Mean = logNormalFitNULLlog2WsigmaMEANpvalue,
  SD = logNormalFitNULLlog2WsigmaSDpvalue
)
logNormalFitNULLlog2Wsigmapvalues

##### CV ----
kfoldID <- kfold_split_grouped(K = 5, StormdataTrain8$StormID)
logNormalFitNULLlog2Wsigmakfoldgroup <- kfold(logNormalFitNULLlog2Wsigma,
                                             folds = kfoldID,
                                             chains = 1,
                                             save_fits = TRUE)
logNormalFitNULLlog2WsigmakfoldPreds <- kfold_predict(logNormalFitNULLlog2Wsigmakfoldgroup)
#logNormalFitNULLlog2WsigmakfoldPreds <- kfold_predict(logNormalFitNULLlog2Wsigmakfold)
logNormalFitNULLlog2WsigmakfoldPredsDat <- logNormalFitNULLlog2WsigmakfoldPreds$yrep
logNormalFitNULLlog2WsigmakfoldPredsMean <- colMeans(logNormalFitNULLlog2WsigmakfoldPredsDat)
logNormalFitNULLlog2WsigmakfoldPredsMed <- apply(logNormalFitNULLlog2WsigmakfoldPredsDat, 2, function(x){quantile(x, 0.5)})
logNormalFitNULLlog2WsigmakfoldPredsLCB <- apply(logNormalFitNULLlog2WsigmakfoldPredsDat, 2, function(x){quantile(x, 0.025)})
logNormalFitNULLlog2WsigmakfoldPredsUCB <- apply(logNormalFitNULLlog2WsigmakfoldPredsDat, 2, function(x){quantile(x, 0.975)})

logNormalFitNULLlog2WsigmakfoldMetrics <- tibble(
  Fit = paste0("logNormalFitNULLlog2Wsigma"),
  MAE_kfold = mean(abs(logNormalFitNULLlog2WsigmakfoldPredsMean - logNormalFitNULLlog2WsigmakfoldPreds$y)),
  MAD_kfold = mean(abs(logNormalFitNULLlog2WsigmakfoldPredsMed - logNormalFitNULLlog2WsigmakfoldPreds$y)),
  COV_kfold = mean(logNormalFitNULLlog2WsigmakfoldPredsLCB < logNormalFitNULLlog2WsigmakfoldPreds$y & logNormalFitNULLlog2WsigmakfoldPreds$y < logNormalFitNULLlog2WsigmakfoldPredsUCB)
)
logNormalFitNULLlog2WsigmakfoldMetrics

#### Identity Link ----
shiftedLogNormalFitNULLlog <- brm(
  bf(
    VMAX ~ log(HWRF)
  ),
  data = StormdataTrain8, 
  family = shifted_lognormal(link = "identity"),
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)

##### Plot ----
print(shiftedLogNormalFitNULLlog, digits = 4)
shiftedLogNormalFitNULLlogppcFit <- pp_check(shiftedLogNormalFitNULLlog, ndraws = 100) + 
  labs(title = "shiftedLogNormalFitNULLlog Fit PPC") +
  theme_bw()
shiftedLogNormalFitNULLlogppcFit

##### LOO ----
shiftedLogNormalFitNULLlogloo <- loo(shiftedLogNormalFitNULLlog)

##### Prediction ----
## Fitted
shiftedLogNormalFitNULLlogFit <- posterior_predict(shiftedLogNormalFitNULLlog)
shiftedLogNormalFitNULLlogFitMean <- colMeans(shiftedLogNormalFitNULLlogFit)
shiftedLogNormalFitNULLlogFitMed <- apply(shiftedLogNormalFitNULLlogFit, 2, function(x){quantile(x, 0.5)})
shiftedLogNormalFitNULLlogFitLCB <- apply(shiftedLogNormalFitNULLlogFit, 2, function(x){quantile(x, 0.025)})
shiftedLogNormalFitNULLlogFitUCB <- apply(shiftedLogNormalFitNULLlogFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
shiftedLogNormalFitNULLlogPreds <- posterior_predict(shiftedLogNormalFitNULLlog, 
                                              newdata = StormdataTest8,
                                              allow_new_levels = TRUE)
shiftedLogNormalFitNULLlogPredsMean <- colMeans(shiftedLogNormalFitNULLlogPreds)
shiftedLogNormalFitNULLlogPredsMed <- apply(shiftedLogNormalFitNULLlogPreds, 2, function(x){quantile(x, 0.5)})
shiftedLogNormalFitNULLlogPredsLCB <- apply(shiftedLogNormalFitNULLlogPreds, 2, function(x){quantile(x, 0.025)})
shiftedLogNormalFitNULLlogPredsUCB <- apply(shiftedLogNormalFitNULLlogPreds, 2, function(x){quantile(x, 0.975)})

shiftedLogNormalFitNULLlogpredMetrics <- tibble(
  Fit = "shiftedLogNormalFitNULLlog",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(shiftedLogNormalFitNULLlogFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(shiftedLogNormalFitNULLlogFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < shiftedLogNormalFitNULLlogFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest3$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(shiftedLogNormalFitNULLlogPredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(shiftedLogNormalFitNULLlogPredsMed - Actual_Yvec)),
  COV_pred = mean(shiftedLogNormalFitNULLlogPredsLCB < Actual_Yvec & Actual_Yvec < shiftedLogNormalFitNULLlogPredsUCB)
)
shiftedLogNormalFitNULLlogpredMetrics

##### PPC ----
###### Quantile 2.5 
shiftedLogNormalFitNULLlogLCBsims <- apply(shiftedLogNormalFitNULLlogFit, 
                                    MARGIN = 1,
                                    function(x){
                                      quantile(x, 0.025)
                                    })
shiftedLogNormalFitNULLlogLCBpvalueVec <- shiftedLogNormalFitNULLlogLCBsims < quantile(StormdataTrain3$VMAX, 0.025)
shiftedLogNormalFitNULLlogLCBpvalue <- sum(shiftedLogNormalFitNULLlogLCBpvalueVec)
shiftedLogNormalFitNULLlogLCBpvalue <- round(shiftedLogNormalFitNULLlogLCBpvalue/4000, 3)
shiftedLogNormalFitNULLlogLCBpvalue <- min(shiftedLogNormalFitNULLlogLCBpvalue, 1 - shiftedLogNormalFitNULLlogLCBpvalue)

shiftedLogNormalFitNULLlog_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           shiftedLogNormalFitNULLlogFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", shiftedLogNormalFitNULLlogLCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#shiftedLogNormalFitNULLlog_ppcLCB

###### Quantile 97.5 
shiftedLogNormalFitNULLlogUCBsims <- apply(shiftedLogNormalFitNULLlogFit, 
                                    MARGIN = 1,
                                    function(x){
                                      quantile(x, 0.975)
                                    })
shiftedLogNormalFitNULLlogUCBpvalueVec <- shiftedLogNormalFitNULLlogUCBsims < quantile(StormdataTrain3$VMAX, 0.975)
shiftedLogNormalFitNULLlogUCBpvalue <- as.numeric(sum(shiftedLogNormalFitNULLlogUCBpvalueVec))
shiftedLogNormalFitNULLlogUCBpvalue <- round(shiftedLogNormalFitNULLlogUCBpvalue/4000, 3)
shiftedLogNormalFitNULLlogUCBpvalue <- min(shiftedLogNormalFitNULLlogUCBpvalue, 1 - shiftedLogNormalFitNULLlogUCBpvalue)

shiftedLogNormalFitNULLlog_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           shiftedLogNormalFitNULLlogFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", shiftedLogNormalFitNULLlogUCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#shiftedLogNormalFitNULLlog_ppcUCB

###### Mean 
shiftedLogNormalFitNULLlogMEANsims <- apply(shiftedLogNormalFitNULLlogFit, 
                                     MARGIN = 1,
                                     function(x){
                                       mean(x)
                                     })
shiftedLogNormalFitNULLlogMEANpvalueVec <- shiftedLogNormalFitNULLlogMEANsims < mean(StormdataTrain3$VMAX)
shiftedLogNormalFitNULLlogMEANpvalue <- sum(shiftedLogNormalFitNULLlogMEANpvalueVec)
shiftedLogNormalFitNULLlogMEANpvalue <- round(shiftedLogNormalFitNULLlogMEANpvalue/4000, 3)
shiftedLogNormalFitNULLlogMEANpvalue <- min(shiftedLogNormalFitNULLlogMEANpvalue, 1 - shiftedLogNormalFitNULLlogMEANpvalue)

shiftedLogNormalFitNULLlog_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           shiftedLogNormalFitNULLlogFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", shiftedLogNormalFitNULLlogMEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#shiftedLogNormalFitNULLlog_ppcMEAN

###### Med 
shiftedLogNormalFitNULLlogMEDsims <- apply(shiftedLogNormalFitNULLlogFit, 
                                    MARGIN = 1,
                                    function(x){
                                      quantile(x, 0.5)
                                    })
shiftedLogNormalFitNULLlogMEDpvalueVec <- shiftedLogNormalFitNULLlogMEDsims < quantile(StormdataTrain3$VMAX, 0.5)
shiftedLogNormalFitNULLlogMEDpvalue <- sum(shiftedLogNormalFitNULLlogMEDpvalueVec)
shiftedLogNormalFitNULLlogMEDpvalue <- round(shiftedLogNormalFitNULLlogMEDpvalue/4000, 3)
shiftedLogNormalFitNULLlogMEDpvalue <- min(shiftedLogNormalFitNULLlogMEDpvalue, 1 - shiftedLogNormalFitNULLlogMEDpvalue)

shiftedLogNormalFitNULLlog_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           shiftedLogNormalFitNULLlogFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", shiftedLogNormalFitNULLlogMEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#shiftedLogNormalFitNULLlog_ppcMED

###### SD 
shiftedLogNormalFitNULLlogSDsims <- apply(shiftedLogNormalFitNULLlogFit, 
                                   MARGIN = 1,
                                   function(x){
                                     sd(x)
                                   })
shiftedLogNormalFitNULLlogSDpvalueVec <- shiftedLogNormalFitNULLlogSDsims < sd(StormdataTrain3$VMAX)
shiftedLogNormalFitNULLlogSDpvalue <- sum(shiftedLogNormalFitNULLlogSDpvalueVec)
shiftedLogNormalFitNULLlogSDpvalue <- round(shiftedLogNormalFitNULLlogSDpvalue/4000, 3)
shiftedLogNormalFitNULLlogSDpvalue <- min(shiftedLogNormalFitNULLlogSDpvalue, 1 - shiftedLogNormalFitNULLlogSDpvalue)

shiftedLogNormalFitNULLlog_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           shiftedLogNormalFitNULLlogFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", shiftedLogNormalFitNULLlogSDpvalue, ")")) +
  theme_bw() +
  legend_none()
#shiftedLogNormalFitNULLlog_ppcSD

###### Range 
shiftedLogNormalFitNULLlogRANGEsims <- apply(shiftedLogNormalFitNULLlogFit, 
                                      MARGIN = 1,
                                      function(x){
                                        max(x)-min(x)
                                      })
shiftedLogNormalFitNULLlogRANGEpvalueVec <- shiftedLogNormalFitNULLlogRANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
shiftedLogNormalFitNULLlogRANGEpvalue <- sum(shiftedLogNormalFitNULLlogRANGEpvalueVec)
shiftedLogNormalFitNULLlogRANGEpvalue <- round(shiftedLogNormalFitNULLlogRANGEpvalue/4000, 3)
shiftedLogNormalFitNULLlogRANGEpvalue <- min(shiftedLogNormalFitNULLlogRANGEpvalue, 1 - shiftedLogNormalFitNULLlogRANGEpvalue)

shiftedLogNormalFitNULLlog_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           shiftedLogNormalFitNULLlogFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", shiftedLogNormalFitNULLlogRANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#shiftedLogNormalFitNULLlog_ppcRANGE

##### Combined Plot ----
shiftedLogNormalFitNULLlog_ppcComb <- 
  shiftedLogNormalFitNULLlogppcFit /
  (shiftedLogNormalFitNULLlog_ppcLCB | shiftedLogNormalFitNULLlog_ppcMED | shiftedLogNormalFitNULLlog_ppcUCB) /
  (shiftedLogNormalFitNULLlog_ppcRANGE | shiftedLogNormalFitNULLlog_ppcMEAN | shiftedLogNormalFitNULLlog_ppcSD)
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
shiftedLogNormalFitNULLlog_ppcComb

##### Bayes p-values ----
shiftedLogNormalFitNULLlogpvalues <- tibble(
  Fit = paste0("shiftedLogNormalFitNULLlog"),
  LCB = shiftedLogNormalFitNULLlogLCBpvalue,
  Median = shiftedLogNormalFitNULLlogMEDpvalue,
  UCB = shiftedLogNormalFitNULLlogUCBpvalue,
  Range = shiftedLogNormalFitNULLlogRANGEpvalue,
  Mean = shiftedLogNormalFitNULLlogMEANpvalue,
  SD = shiftedLogNormalFitNULLlogSDpvalue
)
shiftedLogNormalFitNULLlogpvalues

##### CV ----
kfoldID <- kfold_split_grouped(K = 5, StormdataTrain8$StormID)
shiftedLogNormalFitNULLlogkfoldgroup <- kfold(shiftedLogNormalFitNULLlog,
                                       folds = kfoldID,
                                       chains = 1,
                                       save_fits = TRUE)
shiftedLogNormalFitNULLlogkfoldPreds <- kfold_predict(shiftedLogNormalFitNULLlogkfoldgroup)
#shiftedLogNormalFitNULLlogkfoldPreds <- kfold_predict(shiftedLogNormalFitNULLlogkfold)
shiftedLogNormalFitNULLlogkfoldPredsDat <- shiftedLogNormalFitNULLlogkfoldPreds$yrep
shiftedLogNormalFitNULLlogkfoldPredsMean <- colMeans(shiftedLogNormalFitNULLlogkfoldPredsDat)
shiftedLogNormalFitNULLlogkfoldPredsMed <- apply(shiftedLogNormalFitNULLlogkfoldPredsDat, 2, function(x){quantile(x, 0.5)})
shiftedLogNormalFitNULLlogkfoldPredsLCB <- apply(shiftedLogNormalFitNULLlogkfoldPredsDat, 2, function(x){quantile(x, 0.025)})
shiftedLogNormalFitNULLlogkfoldPredsUCB <- apply(shiftedLogNormalFitNULLlogkfoldPredsDat, 2, function(x){quantile(x, 0.975)})

shiftedLogNormalFitNULLlogkfoldMetrics <- tibble(
  Fit = paste0("shiftedLogNormalFitNULLlog"),
  MAE_kfold = mean(abs(shiftedLogNormalFitNULLlogkfoldPredsMean - shiftedLogNormalFitNULLlogkfoldPreds$y)),
  MAD_kfold = mean(abs(shiftedLogNormalFitNULLlogkfoldPredsMed - shiftedLogNormalFitNULLlogkfoldPreds$y)),
  COV_kfold = mean(shiftedLogNormalFitNULLlogkfoldPredsLCB < shiftedLogNormalFitNULLlogkfoldPreds$y & shiftedLogNormalFitNULLlogkfoldPreds$y < shiftedLogNormalFitNULLlogkfoldPredsUCB)
)
shiftedLogNormalFitNULLlogkfoldMetrics

### Compare NULLs ----
logNormalFitNULLppcFitComb <- 
  logNormalFitNULLppcFit /
  logNormalFitNULLppcFit /
  logNormalFitNULLppcFit /
  logNormalFitNULL2ppcFit /
  logNormalFitNULLscaleppcFit /
  logNormalFitNULLscale2ppcFit /
  logNormalFitNULLlogppcFit /
  logNormalFitNULLlog2ppcFit

(logNormalFitNULLidentity_ppcLCB | 
    logNormalFitNULLidentity_ppcMED | 
    logNormalFitNULLidentity_ppcUCB | 
    logNormalFitNULLidentity_ppcRANGE |
    logNormalFitNULLidentity_ppcMEAN | 
    logNormalFitNULLidentity_ppcSD) /
  (logNormalFitNULLlog2_ppcLCB | 
     logNormalFitNULLlog2_ppcMED | 
     logNormalFitNULLlog2_ppcUCB | 
     logNormalFitNULLlog2_ppcRANGE | 
     logNormalFitNULLlog2_ppcMEAN | 
     logNormalFitNULLlog2_ppcSD) /
  (logNormalFitNULLscaleidentity_ppcLCB | 
     logNormalFitNULLscaleidentity_ppcMED | 
     logNormalFitNULLscaleidentity_ppcUCB | 
     logNormalFitNULLscaleidentity_ppcRANGE | 
     logNormalFitNULLscaleidentity_ppcMEAN | 
     logNormalFitNULLscaleidentity_ppcSD) /
  (logNormalFitNULLlog_ppcLCB | 
     logNormalFitNULLlog_ppcMED | 
     logNormalFitNULLlog_ppcUCB | 
     logNormalFitNULLlog_ppcRANGE | 
     logNormalFitNULLlog_ppcMEAN | 
     logNormalFitNULLlog_ppcSD) 

logNormalFitNULLpvaluesComb <- bind_rows(
  logNormalFitNULLpvalues,
  logNormalFitNULLlogSigmapvalues,
  logNormalFitNULL2pvalues,
  logNormalFitNULL2logSigmapvalues,
  #logNormalFitNULLidentitypvalues,
  logNormalFitNULLscalepvalues,
  logNormalFitNULLscale2pvalues,
  #logNormalFitNULLscaleidentitypvalues,
  logNormalFitNULLlogpvalues,
  logNormalFitNULLlogWsigmapvalues,
  logNormalFitNULLlog2pvalues,
  logNormalFitNULLlog2Wsigmapvalues,
  shiftedLogNormalFitNULLlogpvalues
)
logNormalFitNULLpvaluesComb

logNormalFitNULLlooComp <- loo_compare(
  logNormalFitNULLloo,
  logNormalFitNULLlogSigmaloo,
  logNormalFitNULL2loo,
  logNormalFitNULL2logSigmaloo,
  #logNormalFitNULLidentityloo,
  logNormalFitNULLscaleloo,
  logNormalFitNULLscale2loo,
  #logNormalFitNULLscaleidentityloo,
  logNormalFitNULLlogloo,
  logNormalFitNULLlogWsigmaloo,
  logNormalFitNULLlog2loo,
  logNormalFitNULLlog2Wsigmaloo,
  shiftedLogNormalFitNULLlogloo
)
logNormalFitNULLlooComp

logNormalFitNULLpredComp <- bind_rows(
  logNormalFitNULLpredMetrics,
  logNormalFitNULLlogSigmapredMetrics,
  logNormalFitNULL2predMetrics,
  logNormalFitNULL2logSigmapredMetrics,
  logNormalFitNULLscalepredMetrics,
  logNormalFitNULLscale2predMetrics,
  logNormalFitNULLlogpredMetrics,
  logNormalFitNULLlogWsigmapredMetrics,
  logNormalFitNULLlog2predMetrics,
  logNormalFitNULLlog2WsigmapredMetrics,
  shiftedLogNormalFitNULLlogpredMetrics
) |>
  arrange(MAE_pred)
logNormalFitNULLpredComp

logNormalFitNULLkfoldComp <- bind_rows(
  logNormalFitNULLkfoldMetrics,
  logNormalFitNULLlogSigmakfoldMetrics,
  logNormalFitNULL2kfoldMetrics,
  logNormalFitNULL2logSigmakfoldMetrics,
  logNormalFitNULLscalekfoldMetrics,
  logNormalFitNULLscale2kfoldMetrics,
  logNormalFitNULLlogkfoldMetrics,
  logNormalFitNULLlogWsigmakfoldMetrics,
  logNormalFitNULLlog2kfoldMetrics,
  logNormalFitNULLlog2WsigmakfoldMetrics,
  shiftedLogNormalFitNULLlogkfoldMetrics
) |>
  arrange(MAE_kfold)

logNormalFitNULLkfoldComp

logNormalFitNULLlog
logNormalFitNULLlogWsigma
logNormalFitNULLlog2

### Clear Environment ----
NUllenv <- ls()
exNULL <- str_detect(NUllenv, "NULL")
keepNULL <- NUllenv %in% c("logNormalFitNULLppcFitComb",
                           "logNormalFitNULLpvaluesComb",
                           "logNormalFitNULLlooComp",
                           "logNormalFitNULLpredComp",
                           "logNormalFitNULLkfoldComp",
                           "logNormalFitNULLidentity",
                           "logNormalFitNULLlog2",
                           "logNormalFitNULLscaleidentity",
                           "logNormalFitNULLlog",
                           "logNormalFitNULLlog2",
                           "logNormalFitNULLlog2Wsigma",
                           "logNormalFitNULLlogWsigma",
                           "shiftedLogNormalFitNULLlog")
cutNULL <- exNULL & !keepNULL
NULLrm <- NUllenv[cutNULL]
rm(list = NULLrm)


## logNormalFit1 <- logNormalFit
## logNormalFit2 <- logNormalFit

# Fit Model ----
# Fit1: VMAX ~ offset(log(HWRF)) + (1|StormID)
##      data = StormdataTrain8
# Fit2: VMAX ~ offset(log(HWRF)) + (1|StormID) + s(StormElapsedTime)
##      data = StormdataTrain8
# Fit3: VMAX ~ offset(log(HWRF)) + (1|StormID) + s(StormElapsedTime, m = 1)
##      data = StormdataTrain8
# Fit4: VMAX ~ offset(log(HWRF)) + s(StormElapsedTime) + s(StormID, bs = "re")
##      data = StormdataTrain8
# Fit5: VMAX ~ offset(log(HWRF)) + (1|StormID) + s(StormElapsedTime, bs = "cr") 
##      data = StormdataTrain8
# Fit6: VMAX ~ offset(log(HWRF)) + StormElapsedTime + I(StormElapsedTime^2) + (1|StormID)
##      data = StormdataTrain8
# Fit7: VMAX ~ offset(log(HWRF)) + StormElapsedTime + I(StormElapsedTime^2) +
##            (1 + StormElapsedTime + I(StormElapsedTime^2)|StormID)
##      data = StormdataTrain8
##Fit8: bf(VMAX ~ offset(log(HWRF)) +
# t2(LON, LAT, StormElapsedTime, d = c(2,1), bs = c("tp", "cr")) +
#   MINSLP +
#   SHR_MAG +
#   STM_SPD +
#   SST +
#   RHLO +
#   CAPE1 +
#   CAPE3 +
#   SHTFL2 +
#   TCOND7002 +
#   INST2 +
#   CP1 +
#   TCONDSYM2 +
#   COUPLSYM3 +
#   HWFI +
#   VMAX_OP_T0 +
#   (1|StormID)
# ),
# data = StormdataTrain8

### Model ----
fit <- 9
tic()
logNormalFit <- brm(
  bf(VMAX ~ 
       #0 + Intercept +
       s(Day, bs = "cc") +
       s(StormElapsedTime, bs = "tp") +
       s(LON, bs = "tp") +
       s(LAT, bs = "tp") +
       basin + 
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
  data = StormdataTrain3,
  family = lognormal(link = "identity"),
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000,
  # prior = c(
  #   prior(normal(62, 30), class = b, coef = Intercept)
  # ),
  knots = list(
    Day = c(0, 365)),
  # Day = c(scale(0,
  #               center = attr(StormdataTrain7scale |> pull(Day), "scaled:center"),
  #               scale = attr(StormdataTrain7scale |> pull(Day), "scaled:scale")),
  #         scale(365,
  #               center = attr(StormdataTrain7scale |> pull(Day), "scaled:center"),
  #               scale = attr(StormdataTrain7scale |> pull(Day), "scaled:scale"))
  # )),
  control = list(adapt_delta = 0.95)
)
toc()

logNormalFit2 <- logNormalFit
# save(logNormalFit1,
#      file = "~/Desktop/Temp Hurricane Model Data/logNormalFit1.RData")
prior_summary(logNormalFit)
posterior_summary(logNormalFit)
logNormalFit
launch_shinystan(logNormalFit)

print(logNormalFit, digits = 4)
print(logNormalFit1, digits = 4)
print(logNormalFit2, digits = 4)
print(logNormalFit3, digits = 4)
print(logNormalFit4, digits = 4)
print(logNormalFit5, digits = 4)
print(logNormalFit6, digits = 4)
print(logNormalFit7, digits = 4)
print(logNormalFit8, digits = 4)
print(logNormalFit9, digits = 4)
print(logNormalFit10, digits = 4)
plot(logNormalFit)
logNormalFitppcFit <- pp_check(logNormalFit, ndraws = 100) + 
  labs(title = paste0("logNormalFit", fit, " Fit PPC")) +
  theme_bw()
logNormalFitppcFit
logNormalFitloo <- loo(logNormalFit)
waic(logNormalFit)
performance::check_distribution(logNormalFit)
performance::check_outliers(logNormalFit)
performance::check_heteroskedasticity(logNormalFit)
performance_rmse(logNormalFit)
performance_mae(logNormalFit)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(logNormalFit)

variance_decomposition(logNormalFit)
fixef(logNormalFit)
ranef(logNormalFit)

bayes_R2(logNormalFit)

logNormalFitC <- logNormalFit
logNormalFit <- logNormalFitB
bayes_factor(logNormalFit, logNormalFit1)
bayes_factor(logNormalFit, logNormalFit2)
bayes_factor(logNormalFit, logNormalFit3)
bayes_factor(logNormalFit, logNormalFit4)
bayes_factor(logNormalFit, logNormalFit5)
bayes_factor(logNormalFit, logNormalFit6)
bayes_factor(logNormalFit, logNormalFit7)
bayes_factor(logNormalFit, logNormalFit10)

logNormalFitsmooths <- conditional_smooths(logNormalFit)
plot(logNormalFitsmooths, 
     stype = "contour", 
     ask = FALSE,
     theme = theme(legend.position = "bottom"))

logNormalFiteffects <- conditional_effects(logNormalFit, 
                                           surface = TRUE)
plot(logNormalFiteffects, 
     points = TRUE, 
     ask = FALSE, 
     stype = "contour")

##### Prediction ----
## Fitted
pred1 <- predict(logNormalFit, type = "link")

logNormalFitfinalFit <- posterior_predict(logNormalFit)
logNormalFitfinalResiduals <- t(StormdataTrain3$VMAX - t(logNormalFitfinalFit))
logNormalFitfinalResidualsMean <- colMeans(logNormalFitfinalResiduals)
logNormalFitfinalFitMean <- colMeans(logNormalFitfinalFit)
logNormalFitfinalFitMed <- apply(logNormalFitfinalFit, 2, function(x){quantile(x, 0.5)})
logNormalFitfinalFitLCB <- apply(logNormalFitfinalFit, 2, function(x){quantile(x, 0.025)})
logNormalFitfinalFitUCB <- apply(logNormalFitfinalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
logNormalFitfinalPreds <- posterior_predict(logNormalFit, 
                                            newdata = StormdataTest3,
                                            allow_new_levels = TRUE, 
                                            re_formula = NULL)
logNormalFitfinalPredsMean <- colMeans(logNormalFitfinalPreds)
logNormalFitfinalPredsMed <- apply(logNormalFitfinalPreds, 2, function(x){quantile(x, 0.5, na.rm = TRUE)})
logNormalFitfinalPredsLCB <- apply(logNormalFitfinalPreds, 2, function(x){quantile(x, 0.025, na.rm = TRUE)})
logNormalFitfinalPredsUCB <- apply(logNormalFitfinalPreds, 2, function(x){quantile(x, 0.975, na.rm = TRUE)})

logNormalFitpredMetrics <- tibble(
  Fit = paste0("logNormalFit", fit),
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(logNormalFitfinalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(logNormalFitfinalFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < logNormalFitfinalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest3$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(logNormalFitfinalPredsMean - Actual_Yvec), na.rm = TRUE),
  MAD_pred = mean(abs(logNormalFitfinalPredsMed - Actual_Yvec), na.rm = TRUE),
  COV_pred = mean(logNormalFitfinalPredsLCB < Actual_Yvec & Actual_Yvec < logNormalFitfinalPredsUCB)
)
logNormalFitpredMetrics


##### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = logNormalFitfinalPreds) +
  labs(title = "logNormalFit Predict") +
  theme_bw()

logNormalFitFitDF <- bind_cols(
  StormdataTrain3,
  LCB = logNormalFitfinalFitLCB,
  Mean = logNormalFitfinalFitMean,
  Med = logNormalFitfinalFitMed,
  UCB = logNormalFitfinalFitUCB
)

logNormalFitstormsFitplot <- ggplot(data = logNormalFitFitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX, color = "Observed")) +
  geom_line(aes(y = Mean, color = "PPD Mean")) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50)) +
  labs(title = "logNormalFit PPD Mean vs Observed VMAX",
       subtitle = "95% Credible Interval about PPD Mean") +
  scale_color_manual(name = NULL, values = c("black","red")) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1))
  ) +
  theme_bw()
logNormalFitstormsFitplot

## Residuals
logNormalFitResiduals <- residuals(logNormalFit, method = "posterior_predict")
ppc_error_hist(StormdataTrain8$VMAX, logNormalFitfinalFit[c(1,20,100), ])
ppc_error_scatter_avg(StormdataTrain9$VMAX, logNormalFitfinalFit)
ppc_error_scatter_avg_vs_x(StormdataTrain8$VMAX, 
                           logNormalFitfinalFit,
                           StormdataTrain3$StormElapsedTime)
ppc_error_scatter_avg_vs_x(StormdataTrain8$VMAX, 
                           logNormalFitfinalFit,
                           StormdataTrain3$HWRF)
ppc_scatter_avg_grouped(StormdataTrain9$VMAX, 
                        logNormalFitfinalFit,
                        StormdataTrain3$StormID,
                        facet_args = list(scales = "free_x"))

## Prediction
logNormalFitPredDF <- bind_cols(
  StormdataTest3,
  LCB = logNormalFitfinalPredsLCB,
  Mean = logNormalFitfinalPredsMean,
  Med = logNormalFitfinalPredsMed,
  UCB = logNormalFitfinalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

logNormalFitstormsPredplot <- ggplot(data = logNormalFitPredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX, color = "Observed")) +
  geom_line(aes(y = Mean, color = "PPD Mean")) +
  facet_wrap(vars(StormID), ncol = 6)+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50)) +
  labs(title = "logNormalFit PPD Mean vs Observed VMAX",
       subtitle = "95% Credible Interval about PPD Mean") +
  scale_color_manual(name = NULL, values = c("black","red")) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1))
  ) +
  theme_bw()
logNormalFitstormsPredplot

##### PPC ----
###### Quantile 2.5 
logNormalFitLCBsims <- apply(logNormalFitfinalFit, 
                             MARGIN = 1,
                             function(x){
                               quantile(x, 0.025)
                             })
logNormalFitLCBpvalueVec <- logNormalFitLCBsims < quantile(StormdataTrain3$VMAX, 0.025)
logNormalFitLCBpvalue <- sum(logNormalFitLCBpvalueVec)
logNormalFitLCBpvalue <- round(logNormalFitLCBpvalue/4000, 3)
logNormalFitLCBpvalue <- min(logNormalFitLCBpvalue, 1 - logNormalFitLCBpvalue)

logNormalFit_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitfinalFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", logNormalFitLCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFit_ppcLCB

###### Quantile 97.5 
logNormalFitUCBsims <- apply(logNormalFitfinalFit, 
                             MARGIN = 1,
                             function(x){
                               quantile(x, 0.975)
                             })
logNormalFitUCBpvalueVec <- logNormalFitUCBsims < quantile(StormdataTrain3$VMAX, 0.975)
logNormalFitUCBpvalue <- as.numeric(sum(logNormalFitUCBpvalueVec))
logNormalFitUCBpvalue <- round(logNormalFitUCBpvalue/4000, 3)
logNormalFitUCBpvalue <- min(logNormalFitUCBpvalue, 1 - logNormalFitUCBpvalue)

logNormalFit_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitfinalFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", logNormalFitUCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFit_ppcUCB

###### Mean 
logNormalFitMEANsims <- apply(logNormalFitfinalFit, 
                              MARGIN = 1,
                              function(x){
                                mean(x)
                              })
logNormalFitMEANpvalueVec <- logNormalFitMEANsims < mean(StormdataTrain3$VMAX)
logNormalFitMEANpvalue <- sum(logNormalFitMEANpvalueVec)
logNormalFitMEANpvalue <- round(logNormalFitMEANpvalue/4000, 3)
logNormalFitMEANpvalue <- min(logNormalFitMEANpvalue, 1 - logNormalFitMEANpvalue)

logNormalFit_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitfinalFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", logNormalFitMEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFit_ppcMEAN

###### Med 
logNormalFitMEDsims <- apply(logNormalFitfinalFit, 
                             MARGIN = 1,
                             function(x){
                               quantile(x, 0.5)
                             })
logNormalFitMEDpvalueVec <- logNormalFitMEDsims < quantile(StormdataTrain3$VMAX, 0.5)
logNormalFitMEDpvalue <- sum(logNormalFitMEDpvalueVec)
logNormalFitMEDpvalue <- round(logNormalFitMEDpvalue/4000, 3)
logNormalFitMEDpvalue <- min(logNormalFitMEDpvalue, 1 - logNormalFitMEDpvalue)

logNormalFit_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitfinalFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", logNormalFitMEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFit_ppcMED

###### SD 
logNormalFitSDsims <- apply(logNormalFitfinalFit, 
                            MARGIN = 1,
                            function(x){
                              sd(x)
                            })
logNormalFitSDpvalueVec <- logNormalFitSDsims < sd(StormdataTrain3$VMAX)
logNormalFitSDpvalue <- sum(logNormalFitSDpvalueVec)
logNormalFitSDpvalue <- round(logNormalFitSDpvalue/4000, 3)
logNormalFitSDpvalue <- min(logNormalFitSDpvalue, 1 - logNormalFitSDpvalue)

logNormalFit_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitfinalFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", logNormalFitSDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFit_ppcSD

###### Range 
logNormalFitRANGEsims <- apply(logNormalFitfinalFit, 
                               MARGIN = 1,
                               function(x){
                                 max(x)-min(x)
                               })
logNormalFitRANGEpvalueVec <- logNormalFitRANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
logNormalFitRANGEpvalue <- sum(logNormalFitRANGEpvalueVec)
logNormalFitRANGEpvalue <- round(logNormalFitRANGEpvalue/4000, 3)
logNormalFitRANGEpvalue <- min(logNormalFitRANGEpvalue, 1 - logNormalFitRANGEpvalue)

logNormalFit_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           logNormalFitfinalFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", logNormalFitRANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFit_ppcRANGE

##### Combined Plot ----
logNormalFit_ppcComb <- 
  logNormalFitppcFit /
  (logNormalFit_ppcLCB | logNormalFit_ppcMED | logNormalFit_ppcUCB) /
  (logNormalFit_ppcRANGE | logNormalFit_ppcMEAN | logNormalFit_ppcSD)
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
logNormalFit_ppcComb

##### Bayes p-values ----
logNormalFitpvalues <- tibble(
  Fit = paste0("logNormalFit", fit),
  LCB = logNormalFitLCBpvalue,
  Median = logNormalFitMEDpvalue,
  UCB = logNormalFitUCBpvalue,
  Range = logNormalFitRANGEpvalue,
  Mean = logNormalFitMEANpvalue,
  SD = logNormalFitSDpvalue
)
logNormalFitpvalues

##### CV ----
#kfoldID <- kfold_split_grouped(K = 5, StormdataTrain8$StormID)
logNormalFitkfoldgroup <- kfold(logNormalFit,
                                folds = kfoldID,
                                chains = 1,
                                save_fits = TRUE)
logNormalFitkfoldrand <- kfold(logNormalFit,
                               K = 5,
                               chains = 1,
                               save_fits = TRUE)
logNormalFitkfoldPreds <- kfold_predict(logNormalFitkfoldgroup)
#logNormalFitkfoldPreds <- kfold_predict(logNormalFitkfold)
logNormalFitkfoldPredsDat <- logNormalFitkfoldPreds$yrep
logNormalFitkfoldPredsMean <- colMeans(logNormalFitkfoldPredsDat)
logNormalFitkfoldPredsMed <- apply(logNormalFitkfoldPredsDat, 2, function(x){quantile(x, 0.5)})
logNormalFitkfoldPredsLCB <- apply(logNormalFitkfoldPredsDat, 2, function(x){quantile(x, 0.025)})
logNormalFitkfoldPredsUCB <- apply(logNormalFitkfoldPredsDat, 2, function(x){quantile(x, 0.975)})

logNormalFitkfoldMetrics <- tibble(
  Fit = paste0("logNormalFit", fit),
  MAE_kfold = mean(abs(logNormalFitkfoldPredsMean - logNormalFitkfoldPreds$y)),
  MAD_kfold = mean(abs(logNormalFitkfoldPredsMed - logNormalFitkfoldPreds$y)),
  COV_kfold = mean(logNormalFitkfoldPredsLCB < logNormalFitkfoldPreds$y & logNormalFitkfoldPreds$y < logNormalFitkfoldPredsUCB)
)
logNormalFitkfoldMetrics


# Compare Models ----
## Bayes pvalue ----
pvalueComp_logNormal <- logNormalFitNULLidentitypvalues
pvalueComp_logNormaltemp <- logNormalFitpvalues
pvalueComp_logNormal <- bind_rows(
  pvalueComp_logNormal,
  pvalueComp_logNormaltemp
)
pvalueComp_logNormal
save(pvalueComp_logNormal, 
     file = "~/Desktop/Temp Hurricane Model Data/pvalueComp_logNormal.RData")

## LOO ----
# logNormalFitloo1 <- logNormalFitloo
# attributes(logNormalFitloo1)$model_name <- "logNormalFit1"
logNormalFitloo2 <- logNormalFitloo
attributes(logNormalFitloo2)$model_name <- "logNormalFit2"
# logNormalFitloo3 <- logNormalFitloo
# attributes(logNormalFitloo3)$model_name <- "logNormalFit3"
#logNormalFitloo5 <- loo(logNormalFit5)
#attributes(logNormalFitloo5)$model_name <- "logNormalFit5"
#logNormalFitloo6 <- loo(logNormalFit)
#attributes(logNormalFitloo6)$model_name <- "logNormalFit6"
# logNormalFitloo8 <- loo(logNormalFit)
# attributes(logNormalFitloo8)$model_name <- "logNormalFit8"
# logNormalFitloo9 <- loo(logNormalFit)
# attributes(logNormalFitloo9)$model_name <- "logNormalFit9"
# logNormalFitloo11 <- loo(logNormalFit)
# attributes(logNormalFitloo11)$model_name <- "logNormalFit11"
# logNormalFitloo12 <- loo(logNormalFit)
# attributes(logNormalFitloo12)$model_name <- "logNormalFit12"
looComp_logNormal <- loo_compare(logNormalFitNULLidentityloo,
                                 logNormalFitloo1,
                                 logNormalFitloo2
                                 # logNormalFitloo3,
                                 # logNormalFitloo5,
                                 # logNormalFitloo6,
                                 # logNormalFitloo8,
                                 # logNormalFitloo9,
                                 # logNormalFitloo11,
                                 # logNormalFitloo12
)

looComp_logNormal
save(looComp_logNormal, 
     file = "~/Desktop/Temp Hurricane Model Data/looComp_logNormal.RData")

## CV ----
#cvComp_logNormal <- logNormalFitNULL2kfoldMetrics
cvComp_logNormaltemp <- logNormalFitkfoldMetrics
cvComp_logNormal <- bind_rows(
  cvComp_logNormal,
  cvComp_logNormaltemp
)
cvComp_logNormal <- cvComp_logNormal |> arrange(MAE_kfold)
cvComp_logNormal
save(cvComp_logNormal, 
     file = "~/Desktop/Temp Hurricane Model Data/cvComp_logNormal.RData")


## Preds ----
predComp_logNormal <- logNormalFitNULLidentitypredMetrics
predComp_logNormaltemp <- logNormalFitpredMetrics
predComp_logNormal <- bind_rows(
  predComp_logNormal,
  predComp_logNormaltemp
)
predComp_logNormal <- predComp_logNormal |> arrange(MAE_pred)
predComp_logNormal
save(predComp_logNormal, 
     file = "~/Desktop/Temp Hurricane Model Data/predComp_logNormal.RData")


### Bayes R2
# bayesR2_logNormal <- bayes_R2(logNormalFitNULL2) |> 
#   bind_cols(Fit = "NULL")
bayesR2_logNormaltemp <- bayes_R2(logNormalFit) |> 
  bind_cols(Fit = paste0("logNormalFit", fit))
bayesR2_logNormal <- bind_rows(
  bayesR2_logNormal,
  bayesR2_logNormaltemp
)
bayesR2_logNormal <- bayesR2_logNormal |> arrange(desc(Estimate))
bayesR2_logNormal

logNormalFit6 <- logNormalFit
save(logNormalFit8,
     file = "~/Desktop/Temp Hurricane Model Data/logNormalFit8.RData")


logNormalFitNULL2files <- ls()[str_detect(ls(), pattern = "logNormalFitNULL2")]
logNormalFitNULL2filesRM <- logNormalFitNULL2files[(logNormalFitNULL2files %in% c(
  "logNormalFitNULL2",
  "logNormalFitNULL2_ppcComb",
  "logNormalFitNULL2kfold",
  "logNormalFitNULL2kfoldMetrics",
  "logNormalFitNULL2loo",
  "logNormalFitNULL2predMetrics",
  "logNormalFitNULL2pvalues"))]

rm(list = ls()[!(ls() %in% logNormalFitNULL2filesRM)])
rm(logNormalFitNULLfilesRM)

