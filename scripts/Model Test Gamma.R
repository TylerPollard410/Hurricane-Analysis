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
  ) |>
  mutate(
    across(where(is.numeric) & !c(VMAX, 
                                  HWRF, 
                                  StormElapsedTime,StormElapsedTime2,
                                  LAT, LON),
           function(x){scale(x)})
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
  ) |>
  mutate(
    across(where(is.numeric) & !c(VMAX, HWRF, StormElapsedTime,StormElapsedTime2),
           function(x){scale(x)})
  ) |>
  mutate(
    propVMAX = VMAX/HWRF
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

ddistLog <- descdist(log(StormdataTrain3$VMAX), discrete = FALSE)
summary(ddist)

ddistLogProp <- descdist(log(StormdataTrain3$VMAX/StormdataTrain3$HWRF), discrete = FALSE)
summary(ddist)

fitnorm1 <- fitdist(StormdataTrain3$VMAX, distr = "norm", method = "mle")
summary(fitnorm1)
fitdist(StormdataTrain3$VMAX, distr = "norm", method = "mme")
fitdist(StormdataTrain3$VMAX, distr = "norm", method = "mle")
fitdist(StormdataTrain3$VMAX, distr = "norm", method = "mle")

# GAMMA ----
## NULL Models ----
### No Scaling ----
#### Log Link ----
gammaFitNULL <- brm(
  bf(
    VMAX ~ HWRF
  ),
  data = StormdataTrain8, 
  family = Gamma(link = "log"),
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)

##### Plot ----
print(gammaFitNULL, digits = 4)
gammaFitNULLppcFit <- pp_check(gammaFitNULL, ndraws = 100) + 
  labs(title = "gammaFitNULL Fit PPC") +
  theme_bw()
gammaFitNULLppcFit

##### LOO ----
gammaFitNULLloo <- loo(gammaFitNULL)

##### Prediction ----
## Fitted
gammaFitNULLFit <- posterior_predict(gammaFitNULL)
gammaFitNULLFitMean <- colMeans(gammaFitNULLFit)
gammaFitNULLFitMed <- apply(gammaFitNULLFit, 2, function(x){quantile(x, 0.5)})
gammaFitNULLFitLCB <- apply(gammaFitNULLFit, 2, function(x){quantile(x, 0.025)})
gammaFitNULLFitUCB <- apply(gammaFitNULLFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
gammaFitNULLPreds <- posterior_predict(gammaFitNULL, 
                                       newdata = StormdataTest8,
                                       allow_new_levels = TRUE)
gammaFitNULLPredsMean <- colMeans(gammaFitNULLPreds)
gammaFitNULLPredsMed <- apply(gammaFitNULLPreds, 2, function(x){quantile(x, 0.5)})
gammaFitNULLPredsLCB <- apply(gammaFitNULLPreds, 2, function(x){quantile(x, 0.025)})
gammaFitNULLPredsUCB <- apply(gammaFitNULLPreds, 2, function(x){quantile(x, 0.975)})

gammaFitNULLpredMetrics <- tibble(
  Fit = "gammaFitNULL",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(gammaFitNULLFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(gammaFitNULLFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < gammaFitNULLFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(gammaFitNULLPredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(gammaFitNULLPredsMed - Actual_Yvec)),
  COV_pred = mean(gammaFitNULLPredsLCB < Actual_Yvec & Actual_Yvec < gammaFitNULLPredsUCB)
)
gammaFitNULLpredMetrics

##### PPC ----
###### Quantile 2.5 
gammaFitNULLLCBsims <- apply(gammaFitNULLFit, 
                             MARGIN = 1,
                             function(x){
                               quantile(x, 0.025)
                             })
gammaFitNULLLCBpvalueVec <- gammaFitNULLLCBsims < quantile(StormdataTrain3$VMAX, 0.025)
gammaFitNULLLCBpvalue <- sum(gammaFitNULLLCBpvalueVec)
gammaFitNULLLCBpvalue <- round(gammaFitNULLLCBpvalue/4000, 3)
gammaFitNULLLCBpvalue <- min(gammaFitNULLLCBpvalue, 1 - gammaFitNULLLCBpvalue)

gammaFitNULL_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", gammaFitNULLLCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULL_ppcLCB

###### Quantile 97.5 
gammaFitNULLUCBsims <- apply(gammaFitNULLFit, 
                             MARGIN = 1,
                             function(x){
                               quantile(x, 0.975)
                             })
gammaFitNULLUCBpvalueVec <- gammaFitNULLUCBsims < quantile(StormdataTrain3$VMAX, 0.975)
gammaFitNULLUCBpvalue <- as.numeric(sum(gammaFitNULLUCBpvalueVec))
gammaFitNULLUCBpvalue <- round(gammaFitNULLUCBpvalue/4000, 3)
gammaFitNULLUCBpvalue <- min(gammaFitNULLUCBpvalue, 1 - gammaFitNULLUCBpvalue)

gammaFitNULL_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", gammaFitNULLUCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULL_ppcUCB

###### Mean 
gammaFitNULLMEANsims <- apply(gammaFitNULLFit, 
                              MARGIN = 1,
                              function(x){
                                mean(x)
                              })
gammaFitNULLMEANpvalueVec <- gammaFitNULLMEANsims < mean(StormdataTrain3$VMAX)
gammaFitNULLMEANpvalue <- sum(gammaFitNULLMEANpvalueVec)
gammaFitNULLMEANpvalue <- round(gammaFitNULLMEANpvalue/4000, 3)
gammaFitNULLMEANpvalue <- min(gammaFitNULLMEANpvalue, 1 - gammaFitNULLMEANpvalue)

gammaFitNULL_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", gammaFitNULLMEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULL_ppcMEAN

###### Med 
gammaFitNULLMEDsims <- apply(gammaFitNULLFit, 
                             MARGIN = 1,
                             function(x){
                               quantile(x, 0.5)
                             })
gammaFitNULLMEDpvalueVec <- gammaFitNULLMEDsims < quantile(StormdataTrain3$VMAX, 0.5)
gammaFitNULLMEDpvalue <- sum(gammaFitNULLMEDpvalueVec)
gammaFitNULLMEDpvalue <- round(gammaFitNULLMEDpvalue/4000, 3)
gammaFitNULLMEDpvalue <- min(gammaFitNULLMEDpvalue, 1 - gammaFitNULLMEDpvalue)

gammaFitNULL_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", gammaFitNULLMEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULL_ppcMED

###### SD 
gammaFitNULLSDsims <- apply(gammaFitNULLFit, 
                            MARGIN = 1,
                            function(x){
                              sd(x)
                            })
gammaFitNULLSDpvalueVec <- gammaFitNULLSDsims < sd(StormdataTrain3$VMAX)
gammaFitNULLSDpvalue <- sum(gammaFitNULLSDpvalueVec)
gammaFitNULLSDpvalue <- round(gammaFitNULLSDpvalue/4000, 3)
gammaFitNULLSDpvalue <- min(gammaFitNULLSDpvalue, 1 - gammaFitNULLSDpvalue)

gammaFitNULL_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", gammaFitNULLSDpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULL_ppcSD

###### Range 
gammaFitNULLRANGEsims <- apply(gammaFitNULLFit, 
                               MARGIN = 1,
                               function(x){
                                 max(x)-min(x)
                               })
gammaFitNULLRANGEpvalueVec <- gammaFitNULLRANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
gammaFitNULLRANGEpvalue <- sum(gammaFitNULLRANGEpvalueVec)
gammaFitNULLRANGEpvalue <- round(gammaFitNULLRANGEpvalue/4000, 3)
gammaFitNULLRANGEpvalue <- min(gammaFitNULLRANGEpvalue, 1 - gammaFitNULLRANGEpvalue)

gammaFitNULL_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", gammaFitNULLRANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULL_ppcRANGE

##### Combined Plot ----
gammaFitNULL_ppcComb <- 
  gammaFitNULLppcFit /
  (gammaFitNULL_ppcLCB | gammaFitNULL_ppcMED | gammaFitNULL_ppcUCB) /
  (gammaFitNULL_ppcRANGE | gammaFitNULL_ppcMEAN | gammaFitNULL_ppcSD)
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
gammaFitNULL_ppcComb

##### Bayes p-values ----
gammaFitNULLpvalues <- tibble(
  Fit = paste0("gammaFitNULL"),
  LCB = gammaFitNULLLCBpvalue,
  Median = gammaFitNULLMEDpvalue,
  UCB = gammaFitNULLUCBpvalue,
  Range = gammaFitNULLRANGEpvalue,
  Mean = gammaFitNULLMEANpvalue,
  SD = gammaFitNULLSDpvalue
)
gammaFitNULLpvalues

##### CV ----
kfoldID <- kfold_split_grouped(K = 5, StormdataTrain8$StormID)
gammaFitNULLkfoldgroup <- kfold(gammaFitNULL,
                                folds = kfoldID,
                                chains = 1,
                                save_fits = TRUE)
gammaFitNULLkfoldrand <- kfold(gammaFitNULL,
                               K = 5,
                               chains = 1,
                               save_fits = TRUE)
gammaFitNULLkfoldPreds <- kfold_predict(gammaFitNULLkfoldgroup)
#gammaFitNULLkfoldPreds <- kfold_predict(gammaFitNULLkfold)
gammaFitNULLkfoldPredsDat <- gammaFitNULLkfoldPreds$yrep
gammaFitNULLkfoldPredsMean <- colMeans(gammaFitNULLkfoldPredsDat)
gammaFitNULLkfoldPredsMed <- apply(gammaFitNULLkfoldPredsDat, 2, function(x){quantile(x, 0.5)})
gammaFitNULLkfoldPredsLCB <- apply(gammaFitNULLkfoldPredsDat, 2, function(x){quantile(x, 0.025)})
gammaFitNULLkfoldPredsUCB <- apply(gammaFitNULLkfoldPredsDat, 2, function(x){quantile(x, 0.975)})

gammaFitNULLkfoldMetrics <- tibble(
  Fit = paste0("gammaFitNULL"),
  MAE_kfold = mean(abs(gammaFitNULLkfoldPredsMean - gammaFitNULLkfoldPreds$y)),
  MAD_kfold = mean(abs(gammaFitNULLkfoldPredsMed - gammaFitNULLkfoldPreds$y)),
  COV_kfold = mean(gammaFitNULLkfoldPredsLCB < gammaFitNULLkfoldPreds$y & gammaFitNULLkfoldPreds$y < gammaFitNULLkfoldPredsUCB)
)
gammaFitNULLkfoldMetrics

#### Log Link No Int ----
gammaFitNULL2 <- brm(
  bf(
    VMAX ~ 0 + Intercept + HWRF
  ),
  data = StormdataTrain8, 
  family = Gamma(link = "log"),
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)

##### Plot ----
print(gammaFitNULL2, digits = 4)
gammaFitNULL2ppcFit <- pp_check(gammaFitNULL2, ndraws = 100) + 
  labs(title = "gammaFitNULL2 Fit PPC") +
  theme_bw()
gammaFitNULL2ppcFit

##### LOO ----
gammaFitNULL2loo <- loo(gammaFitNULL2)

##### Prediction ----
## Fitted
gammaFitNULL2Fit <- posterior_predict(gammaFitNULL2)
gammaFitNULL2FitMean <- colMeans(gammaFitNULL2Fit)
gammaFitNULL2FitMed <- apply(gammaFitNULL2Fit, 2, function(x){quantile(x, 0.5)})
gammaFitNULL2FitLCB <- apply(gammaFitNULL2Fit, 2, function(x){quantile(x, 0.025)})
gammaFitNULL2FitUCB <- apply(gammaFitNULL2Fit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
gammaFitNULL2Preds <- posterior_predict(gammaFitNULL2, 
                                        newdata = StormdataTest8,
                                        allow_new_levels = TRUE)
gammaFitNULL2PredsMean <- colMeans(gammaFitNULL2Preds)
gammaFitNULL2PredsMed <- apply(gammaFitNULL2Preds, 2, function(x){quantile(x, 0.5)})
gammaFitNULL2PredsLCB <- apply(gammaFitNULL2Preds, 2, function(x){quantile(x, 0.025)})
gammaFitNULL2PredsUCB <- apply(gammaFitNULL2Preds, 2, function(x){quantile(x, 0.975)})

gammaFitNULL2predMetrics <- tibble(
  Fit = "gammaFitNULL2",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(gammaFitNULL2FitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(gammaFitNULL2FitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < gammaFitNULL2FitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(gammaFitNULL2PredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(gammaFitNULL2PredsMed - Actual_Yvec)),
  COV_pred = mean(gammaFitNULL2PredsLCB < Actual_Yvec & Actual_Yvec < gammaFitNULL2PredsUCB)
)
gammaFitNULL2predMetrics

##### PPC ----
###### Quantile 2.5 
gammaFitNULL2LCBsims <- apply(gammaFitNULL2Fit, 
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
           gammaFitNULL2Fit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", gammaFitNULL2LCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULL2_ppcLCB

###### Quantile 97.5 
gammaFitNULL2UCBsims <- apply(gammaFitNULL2Fit, 
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
           gammaFitNULL2Fit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", gammaFitNULL2UCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULL2_ppcUCB

###### Mean 
gammaFitNULL2MEANsims <- apply(gammaFitNULL2Fit, 
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
           gammaFitNULL2Fit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", gammaFitNULL2MEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULL2_ppcMEAN

###### Med 
gammaFitNULL2MEDsims <- apply(gammaFitNULL2Fit, 
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
           gammaFitNULL2Fit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", gammaFitNULL2MEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULL2_ppcMED

###### SD 
gammaFitNULL2SDsims <- apply(gammaFitNULL2Fit, 
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
           gammaFitNULL2Fit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", gammaFitNULL2SDpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULL2_ppcSD

###### Range 
gammaFitNULL2RANGEsims <- apply(gammaFitNULL2Fit, 
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
           gammaFitNULL2Fit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", gammaFitNULL2RANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULL2_ppcRANGE

##### Combined Plot ----
gammaFitNULL2_ppcComb <- 
  gammaFitNULL2ppcFit /
  (gammaFitNULL2_ppcLCB | gammaFitNULL2_ppcMED | gammaFitNULL2_ppcUCB) /
  (gammaFitNULL2_ppcRANGE | gammaFitNULL2_ppcMEAN | gammaFitNULL2_ppcSD)
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
gammaFitNULL2_ppcComb

##### Bayes p-values ----
gammaFitNULL2pvalues <- tibble(
  Fit = paste0("gammaFitNULL2"),
  LCB = gammaFitNULL2LCBpvalue,
  Median = gammaFitNULL2MEDpvalue,
  UCB = gammaFitNULL2UCBpvalue,
  Range = gammaFitNULL2RANGEpvalue,
  Mean = gammaFitNULL2MEANpvalue,
  SD = gammaFitNULL2SDpvalue
)
gammaFitNULL2pvalues

##### CV ----
kfoldID <- kfold_split_grouped(K = 5, StormdataTrain8$StormID)
gammaFitNULL2kfoldgroup <- kfold(gammaFitNULL2,
                                 folds = kfoldID,
                                 chains = 1,
                                 save_fits = TRUE)
gammaFitNULL2kfoldrand <- kfold(gammaFitNULL2,
                                K = 5,
                                chains = 1,
                                save_fits = TRUE)
gammaFitNULL2kfoldPreds <- kfold_predict(gammaFitNULL2kfoldgroup)
#gammaFitNULL2kfoldPreds <- kfold_predict(gammaFitNULL2kfold)
gammaFitNULL2kfoldPredsDat <- gammaFitNULL2kfoldPreds$yrep
gammaFitNULL2kfoldPredsMean <- colMeans(gammaFitNULL2kfoldPredsDat)
gammaFitNULL2kfoldPredsMed <- apply(gammaFitNULL2kfoldPredsDat, 2, function(x){quantile(x, 0.5)})
gammaFitNULL2kfoldPredsLCB <- apply(gammaFitNULL2kfoldPredsDat, 2, function(x){quantile(x, 0.025)})
gammaFitNULL2kfoldPredsUCB <- apply(gammaFitNULL2kfoldPredsDat, 2, function(x){quantile(x, 0.975)})

gammaFitNULL2kfoldMetrics <- tibble(
  Fit = paste0("gammaFitNULL2"),
  MAE_kfold = mean(abs(gammaFitNULL2kfoldPredsMean - gammaFitNULL2kfoldPreds$y)),
  MAD_kfold = mean(abs(gammaFitNULL2kfoldPredsMed - gammaFitNULL2kfoldPreds$y)),
  COV_kfold = mean(gammaFitNULL2kfoldPredsLCB < gammaFitNULL2kfoldPreds$y & gammaFitNULL2kfoldPreds$y < gammaFitNULL2kfoldPredsUCB)
)
gammaFitNULL2kfoldMetrics

#### Identity link No Int ----
gammaFitNULLidentity <- brm(
  bf(
    VMAX ~ 0 + Intercept + HWRF
  ),
  data = StormdataTrain8, 
  family = brmsfamily(family = "Gamma", link = "identity"),
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)

##### Plot ----
print(gammaFitNULLidentity, digits = 4)
gammaFitNULLidentityppcFit <- pp_check(gammaFitNULLidentity, ndraws = 100) + 
  labs(title = "gammaFitNULLidentity Fit PPC") +
  theme_bw()
gammaFitNULLidentityppcFit

##### LOO ----
gammaFitNULLidentityloo <- loo(gammaFitNULLidentity)

##### Prediction ----
## Fitted
gammaFitNULLidentityFit <- posterior_predict(gammaFitNULLidentity)
gammaFitNULLidentityFitMean <- colMeans(gammaFitNULLidentityFit)
gammaFitNULLidentityFitMed <- apply(gammaFitNULLidentityFit, 2, function(x){quantile(x, 0.5)})
gammaFitNULLidentityFitLCB <- apply(gammaFitNULLidentityFit, 2, function(x){quantile(x, 0.025)})
gammaFitNULLidentityFitUCB <- apply(gammaFitNULLidentityFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
gammaFitNULLidentityPreds <- posterior_predict(gammaFitNULLidentity, 
                                               newdata = StormdataTest8,
                                               allow_new_levels = TRUE)
gammaFitNULLidentityPredsMean <- colMeans(gammaFitNULLidentityPreds)
gammaFitNULLidentityPredsMed <- apply(gammaFitNULLidentityPreds, 2, function(x){quantile(x, 0.5)})
gammaFitNULLidentityPredsLCB <- apply(gammaFitNULLidentityPreds, 2, function(x){quantile(x, 0.025)})
gammaFitNULLidentityPredsUCB <- apply(gammaFitNULLidentityPreds, 2, function(x){quantile(x, 0.975)})

gammaFitNULLidentitypredMetrics <- tibble(
  Fit = "gammaFitNULLidentity",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(gammaFitNULLidentityFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(gammaFitNULLidentityFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < gammaFitNULLidentityFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(gammaFitNULLidentityPredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(gammaFitNULLidentityPredsMed - Actual_Yvec)),
  COV_pred = mean(gammaFitNULLidentityPredsLCB < Actual_Yvec & Actual_Yvec < gammaFitNULLidentityPredsUCB)
)
gammaFitNULLidentitypredMetrics

##### PPC ----
###### Quantile 2.5 
gammaFitNULLidentityLCBsims <- apply(gammaFitNULLidentityFit, 
                                     MARGIN = 1,
                                     function(x){
                                       quantile(x, 0.025)
                                     })
gammaFitNULLidentityLCBpvalueVec <- gammaFitNULLidentityLCBsims < quantile(StormdataTrain3$VMAX, 0.025)
gammaFitNULLidentityLCBpvalue <- sum(gammaFitNULLidentityLCBpvalueVec)
gammaFitNULLidentityLCBpvalue <- round(gammaFitNULLidentityLCBpvalue/4000, 3)
gammaFitNULLidentityLCBpvalue <- min(gammaFitNULLidentityLCBpvalue, 1 - gammaFitNULLidentityLCBpvalue)

gammaFitNULLidentity_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLidentityFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", gammaFitNULLidentityLCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLidentity_ppcLCB

###### Quantile 97.5 
gammaFitNULLidentityUCBsims <- apply(gammaFitNULLidentityFit, 
                                     MARGIN = 1,
                                     function(x){
                                       quantile(x, 0.975)
                                     })
gammaFitNULLidentityUCBpvalueVec <- gammaFitNULLidentityUCBsims < quantile(StormdataTrain3$VMAX, 0.975)
gammaFitNULLidentityUCBpvalue <- as.numeric(sum(gammaFitNULLidentityUCBpvalueVec))
gammaFitNULLidentityUCBpvalue <- round(gammaFitNULLidentityUCBpvalue/4000, 3)
gammaFitNULLidentityUCBpvalue <- min(gammaFitNULLidentityUCBpvalue, 1 - gammaFitNULLidentityUCBpvalue)

gammaFitNULLidentity_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLidentityFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", gammaFitNULLidentityUCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLidentity_ppcUCB

###### Mean 
gammaFitNULLidentityMEANsims <- apply(gammaFitNULLidentityFit, 
                                      MARGIN = 1,
                                      function(x){
                                        mean(x)
                                      })
gammaFitNULLidentityMEANpvalueVec <- gammaFitNULLidentityMEANsims < mean(StormdataTrain3$VMAX)
gammaFitNULLidentityMEANpvalue <- sum(gammaFitNULLidentityMEANpvalueVec)
gammaFitNULLidentityMEANpvalue <- round(gammaFitNULLidentityMEANpvalue/4000, 3)
gammaFitNULLidentityMEANpvalue <- min(gammaFitNULLidentityMEANpvalue, 1 - gammaFitNULLidentityMEANpvalue)

gammaFitNULLidentity_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLidentityFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", gammaFitNULLidentityMEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLidentity_ppcMEAN

###### Med 
gammaFitNULLidentityMEDsims <- apply(gammaFitNULLidentityFit, 
                                     MARGIN = 1,
                                     function(x){
                                       quantile(x, 0.5)
                                     })
gammaFitNULLidentityMEDpvalueVec <- gammaFitNULLidentityMEDsims < quantile(StormdataTrain3$VMAX, 0.5)
gammaFitNULLidentityMEDpvalue <- sum(gammaFitNULLidentityMEDpvalueVec)
gammaFitNULLidentityMEDpvalue <- round(gammaFitNULLidentityMEDpvalue/4000, 3)
gammaFitNULLidentityMEDpvalue <- min(gammaFitNULLidentityMEDpvalue, 1 - gammaFitNULLidentityMEDpvalue)

gammaFitNULLidentity_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLidentityFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", gammaFitNULLidentityMEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLidentity_ppcMED

###### SD 
gammaFitNULLidentitySDsims <- apply(gammaFitNULLidentityFit, 
                                    MARGIN = 1,
                                    function(x){
                                      sd(x)
                                    })
gammaFitNULLidentitySDpvalueVec <- gammaFitNULLidentitySDsims < sd(StormdataTrain3$VMAX)
gammaFitNULLidentitySDpvalue <- sum(gammaFitNULLidentitySDpvalueVec)
gammaFitNULLidentitySDpvalue <- round(gammaFitNULLidentitySDpvalue/4000, 3)
gammaFitNULLidentitySDpvalue <- min(gammaFitNULLidentitySDpvalue, 1 - gammaFitNULLidentitySDpvalue)

gammaFitNULLidentity_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLidentityFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", gammaFitNULLidentitySDpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLidentity_ppcSD

###### Range 
gammaFitNULLidentityRANGEsims <- apply(gammaFitNULLidentityFit, 
                                       MARGIN = 1,
                                       function(x){
                                         max(x)-min(x)
                                       })
gammaFitNULLidentityRANGEpvalueVec <- gammaFitNULLidentityRANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
gammaFitNULLidentityRANGEpvalue <- sum(gammaFitNULLidentityRANGEpvalueVec)
gammaFitNULLidentityRANGEpvalue <- round(gammaFitNULLidentityRANGEpvalue/4000, 3)
gammaFitNULLidentityRANGEpvalue <- min(gammaFitNULLidentityRANGEpvalue, 1 - gammaFitNULLidentityRANGEpvalue)

gammaFitNULLidentity_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLidentityFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", gammaFitNULLidentityRANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLidentity_ppcRANGE

##### Combined Plot ----
gammaFitNULLidentity_ppcComb <- 
  gammaFitNULLidentityppcFit /
  (gammaFitNULLidentity_ppcLCB | gammaFitNULLidentity_ppcMED | gammaFitNULLidentity_ppcUCB) /
  (gammaFitNULLidentity_ppcRANGE | gammaFitNULLidentity_ppcMEAN | gammaFitNULLidentity_ppcSD)
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
gammaFitNULLidentity_ppcComb

##### Bayes p-values ----
gammaFitNULLidentitypvalues <- tibble(
  Fit = paste0("gammaFitNULLidentity"),
  LCB = gammaFitNULLidentityLCBpvalue,
  Median = gammaFitNULLidentityMEDpvalue,
  UCB = gammaFitNULLidentityUCBpvalue,
  Range = gammaFitNULLidentityRANGEpvalue,
  Mean = gammaFitNULLidentityMEANpvalue,
  SD = gammaFitNULLidentitySDpvalue
)
gammaFitNULLidentitypvalues

##### CV ----
kfoldID <- kfold_split_grouped(K = 5, StormdataTrain8$StormID)
gammaFitNULLidentitykfoldgroup <- kfold(gammaFitNULLidentity,
                                        folds = kfoldID,
                                        chains = 1,
                                        save_fits = TRUE)
gammaFitNULLidentitykfoldrand <- kfold(gammaFitNULLidentity,
                                       K = 5,
                                       chains = 1,
                                       save_fits = TRUE)
gammaFitNULLidentitykfoldPreds <- kfold_predict(gammaFitNULLidentitykfoldgroup)
#gammaFitNULLidentitykfoldPreds <- kfold_predict(gammaFitNULLidentitykfold)
gammaFitNULLidentitykfoldPredsDat <- gammaFitNULLidentitykfoldPreds$yrep
gammaFitNULLidentitykfoldPredsMean <- colMeans(gammaFitNULLidentitykfoldPredsDat)
gammaFitNULLidentitykfoldPredsMed <- apply(gammaFitNULLidentitykfoldPredsDat, 2, function(x){quantile(x, 0.5)})
gammaFitNULLidentitykfoldPredsLCB <- apply(gammaFitNULLidentitykfoldPredsDat, 2, function(x){quantile(x, 0.025)})
gammaFitNULLidentitykfoldPredsUCB <- apply(gammaFitNULLidentitykfoldPredsDat, 2, function(x){quantile(x, 0.975)})

gammaFitNULLidentitykfoldMetrics <- tibble(
  Fit = paste0("gammaFitNULLidentity"),
  MAE_kfold = mean(abs(gammaFitNULLidentitykfoldPredsMean - gammaFitNULLidentitykfoldPreds$y)),
  MAD_kfold = mean(abs(gammaFitNULLidentitykfoldPredsMed - gammaFitNULLidentitykfoldPreds$y)),
  COV_kfold = mean(gammaFitNULLidentitykfoldPredsLCB < gammaFitNULLidentitykfoldPreds$y & gammaFitNULLidentitykfoldPreds$y < gammaFitNULLidentitykfoldPredsUCB)
)
gammaFitNULLidentitykfoldMetrics

### Scaling ----
#### Log Link ----
gammaFitNULLscale <- brm(
  bf(
    VMAX ~ HWRF
  ),
  data = StormdataTrain7scale, 
  family = Gamma(link = "log"),
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)

##### Plot ----
print(gammaFitNULLscale, digits = 4)
gammaFitNULLscaleppcFit <- pp_check(gammaFitNULLscale, ndraws = 100) + 
  labs(title = "gammaFitNULLscale Fit PPC") +
  theme_bw()
gammaFitNULLscaleppcFit

##### LOO ----
gammaFitNULLscaleloo <- loo(gammaFitNULLscale)

##### Prediction ----
## Fitted
gammaFitNULLscaleFit <- posterior_predict(gammaFitNULLscale)
gammaFitNULLscaleFitMean <- colMeans(gammaFitNULLscaleFit)
gammaFitNULLscaleFitMed <- apply(gammaFitNULLscaleFit, 2, function(x){quantile(x, 0.5)})
gammaFitNULLscaleFitLCB <- apply(gammaFitNULLscaleFit, 2, function(x){quantile(x, 0.025)})
gammaFitNULLscaleFitUCB <- apply(gammaFitNULLscaleFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
gammaFitNULLscalePreds <- posterior_predict(gammaFitNULLscale, 
                                            newdata = StormdataTest7scale,
                                            allow_new_levels = TRUE)
gammaFitNULLscalePredsMean <- colMeans(gammaFitNULLscalePreds)
gammaFitNULLscalePredsMed <- apply(gammaFitNULLscalePreds, 2, function(x){quantile(x, 0.5)})
gammaFitNULLscalePredsLCB <- apply(gammaFitNULLscalePreds, 2, function(x){quantile(x, 0.025)})
gammaFitNULLscalePredsUCB <- apply(gammaFitNULLscalePreds, 2, function(x){quantile(x, 0.975)})

gammaFitNULLscalepredMetrics <- tibble(
  Fit = "gammaFitNULLscale",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(gammaFitNULLscaleFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(gammaFitNULLscaleFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < gammaFitNULLscaleFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(gammaFitNULLscalePredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(gammaFitNULLscalePredsMed - Actual_Yvec)),
  COV_pred = mean(gammaFitNULLscalePredsLCB < Actual_Yvec & Actual_Yvec < gammaFitNULLscalePredsUCB)
)
gammaFitNULLscalepredMetrics

##### PPC ----
###### Quantile 2.5 
gammaFitNULLscaleLCBsims <- apply(gammaFitNULLscaleFit, 
                                  MARGIN = 1,
                                  function(x){
                                    quantile(x, 0.025)
                                  })
gammaFitNULLscaleLCBpvalueVec <- gammaFitNULLscaleLCBsims < quantile(StormdataTrain3$VMAX, 0.025)
gammaFitNULLscaleLCBpvalue <- sum(gammaFitNULLscaleLCBpvalueVec)
gammaFitNULLscaleLCBpvalue <- round(gammaFitNULLscaleLCBpvalue/4000, 3)
gammaFitNULLscaleLCBpvalue <- min(gammaFitNULLscaleLCBpvalue, 1 - gammaFitNULLscaleLCBpvalue)

gammaFitNULLscale_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLscaleFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", gammaFitNULLscaleLCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLscale_ppcLCB

###### Quantile 97.5 
gammaFitNULLscaleUCBsims <- apply(gammaFitNULLscaleFit, 
                                  MARGIN = 1,
                                  function(x){
                                    quantile(x, 0.975)
                                  })
gammaFitNULLscaleUCBpvalueVec <- gammaFitNULLscaleUCBsims < quantile(StormdataTrain3$VMAX, 0.975)
gammaFitNULLscaleUCBpvalue <- as.numeric(sum(gammaFitNULLscaleUCBpvalueVec))
gammaFitNULLscaleUCBpvalue <- round(gammaFitNULLscaleUCBpvalue/4000, 3)
gammaFitNULLscaleUCBpvalue <- min(gammaFitNULLscaleUCBpvalue, 1 - gammaFitNULLscaleUCBpvalue)

gammaFitNULLscale_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLscaleFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", gammaFitNULLscaleUCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLscale_ppcUCB

###### Mean 
gammaFitNULLscaleMEANsims <- apply(gammaFitNULLscaleFit, 
                                   MARGIN = 1,
                                   function(x){
                                     mean(x)
                                   })
gammaFitNULLscaleMEANpvalueVec <- gammaFitNULLscaleMEANsims < mean(StormdataTrain3$VMAX)
gammaFitNULLscaleMEANpvalue <- sum(gammaFitNULLscaleMEANpvalueVec)
gammaFitNULLscaleMEANpvalue <- round(gammaFitNULLscaleMEANpvalue/4000, 3)
gammaFitNULLscaleMEANpvalue <- min(gammaFitNULLscaleMEANpvalue, 1 - gammaFitNULLscaleMEANpvalue)

gammaFitNULLscale_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLscaleFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", gammaFitNULLscaleMEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLscale_ppcMEAN

###### Med 
gammaFitNULLscaleMEDsims <- apply(gammaFitNULLscaleFit, 
                                  MARGIN = 1,
                                  function(x){
                                    quantile(x, 0.5)
                                  })
gammaFitNULLscaleMEDpvalueVec <- gammaFitNULLscaleMEDsims < quantile(StormdataTrain3$VMAX, 0.5)
gammaFitNULLscaleMEDpvalue <- sum(gammaFitNULLscaleMEDpvalueVec)
gammaFitNULLscaleMEDpvalue <- round(gammaFitNULLscaleMEDpvalue/4000, 3)
gammaFitNULLscaleMEDpvalue <- min(gammaFitNULLscaleMEDpvalue, 1 - gammaFitNULLscaleMEDpvalue)

gammaFitNULLscale_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLscaleFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", gammaFitNULLscaleMEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLscale_ppcMED

###### SD 
gammaFitNULLscaleSDsims <- apply(gammaFitNULLscaleFit, 
                                 MARGIN = 1,
                                 function(x){
                                   sd(x)
                                 })
gammaFitNULLscaleSDpvalueVec <- gammaFitNULLscaleSDsims < sd(StormdataTrain3$VMAX)
gammaFitNULLscaleSDpvalue <- sum(gammaFitNULLscaleSDpvalueVec)
gammaFitNULLscaleSDpvalue <- round(gammaFitNULLscaleSDpvalue/4000, 3)
gammaFitNULLscaleSDpvalue <- min(gammaFitNULLscaleSDpvalue, 1 - gammaFitNULLscaleSDpvalue)

gammaFitNULLscale_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLscaleFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", gammaFitNULLscaleSDpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLscale_ppcSD

###### Range 
gammaFitNULLscaleRANGEsims <- apply(gammaFitNULLscaleFit, 
                                    MARGIN = 1,
                                    function(x){
                                      max(x)-min(x)
                                    })
gammaFitNULLscaleRANGEpvalueVec <- gammaFitNULLscaleRANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
gammaFitNULLscaleRANGEpvalue <- sum(gammaFitNULLscaleRANGEpvalueVec)
gammaFitNULLscaleRANGEpvalue <- round(gammaFitNULLscaleRANGEpvalue/4000, 3)
gammaFitNULLscaleRANGEpvalue <- min(gammaFitNULLscaleRANGEpvalue, 1 - gammaFitNULLscaleRANGEpvalue)

gammaFitNULLscale_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLscaleFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", gammaFitNULLscaleRANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLscale_ppcRANGE

##### Combined Plot ----
gammaFitNULLscale_ppcComb <- 
  gammaFitNULLscaleppcFit /
  (gammaFitNULLscale_ppcLCB | gammaFitNULLscale_ppcMED | gammaFitNULLscale_ppcUCB) /
  (gammaFitNULLscale_ppcRANGE | gammaFitNULLscale_ppcMEAN | gammaFitNULLscale_ppcSD)
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
gammaFitNULLscale_ppcComb

##### Bayes p-values ----
gammaFitNULLscalepvalues <- tibble(
  Fit = paste0("gammaFitNULLscale"),
  LCB = gammaFitNULLscaleLCBpvalue,
  Median = gammaFitNULLscaleMEDpvalue,
  UCB = gammaFitNULLscaleUCBpvalue,
  Range = gammaFitNULLscaleRANGEpvalue,
  Mean = gammaFitNULLscaleMEANpvalue,
  SD = gammaFitNULLscaleSDpvalue
)
gammaFitNULLscalepvalues

##### CV ----
kfoldID <- kfold_split_grouped(K = 5, StormdataTrain8$StormID)
gammaFitNULLscalekfoldgroup <- kfold(gammaFitNULLscale,
                                     folds = kfoldID,
                                     chains = 1,
                                     save_fits = TRUE)
gammaFitNULLscalekfoldrand <- kfold(gammaFitNULLscale,
                                    K = 5,
                                    chains = 1,
                                    save_fits = TRUE)
gammaFitNULLscalekfoldPreds <- kfold_predict(gammaFitNULLscalekfoldgroup)
#gammaFitNULLscalekfoldPreds <- kfold_predict(gammaFitNULLscalekfold)
gammaFitNULLscalekfoldPredsDat <- gammaFitNULLscalekfoldPreds$yrep
gammaFitNULLscalekfoldPredsMean <- colMeans(gammaFitNULLscalekfoldPredsDat)
gammaFitNULLscalekfoldPredsMed <- apply(gammaFitNULLscalekfoldPredsDat, 2, function(x){quantile(x, 0.5)})
gammaFitNULLscalekfoldPredsLCB <- apply(gammaFitNULLscalekfoldPredsDat, 2, function(x){quantile(x, 0.025)})
gammaFitNULLscalekfoldPredsUCB <- apply(gammaFitNULLscalekfoldPredsDat, 2, function(x){quantile(x, 0.975)})

gammaFitNULLscalekfoldMetrics <- tibble(
  Fit = paste0("gammaFitNULLscale"),
  MAE_kfold = mean(abs(gammaFitNULLscalekfoldPredsMean - gammaFitNULLscalekfoldPreds$y)),
  MAD_kfold = mean(abs(gammaFitNULLscalekfoldPredsMed - gammaFitNULLscalekfoldPreds$y)),
  COV_kfold = mean(gammaFitNULLscalekfoldPredsLCB < gammaFitNULLscalekfoldPreds$y & gammaFitNULLscalekfoldPreds$y < gammaFitNULLscalekfoldPredsUCB)
)
gammaFitNULLscalekfoldMetrics

#### Log Link No Int ----
gammaFitNULLscale2 <- brm(
  bf(
    VMAX ~ 0 + Intercept + HWRF
  ),
  data = StormdataTrain7scale, 
  family = Gamma(link = "log"),
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)

##### Plot ----
print(gammaFitNULLscale2, digits = 4)
gammaFitNULLscale2ppcFit <- pp_check(gammaFitNULLscale2, ndraws = 100) + 
  labs(title = "gammaFitNULLscale2 Fit PPC") +
  theme_bw()
gammaFitNULLscale2ppcFit

##### LOO ----
gammaFitNULLscale2loo <- loo(gammaFitNULLscale2)

##### Prediction ----
## Fitted
gammaFitNULLscale2Fit <- posterior_predict(gammaFitNULLscale2)
gammaFitNULLscale2FitMean <- colMeans(gammaFitNULLscale2Fit)
gammaFitNULLscale2FitMed <- apply(gammaFitNULLscale2Fit, 2, function(x){quantile(x, 0.5)})
gammaFitNULLscale2FitLCB <- apply(gammaFitNULLscale2Fit, 2, function(x){quantile(x, 0.025)})
gammaFitNULLscale2FitUCB <- apply(gammaFitNULLscale2Fit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
gammaFitNULLscale2Preds <- posterior_predict(gammaFitNULLscale2, 
                                             newdata = StormdataTest7scale,
                                             allow_new_levels = TRUE)
gammaFitNULLscale2PredsMean <- colMeans(gammaFitNULLscale2Preds)
gammaFitNULLscale2PredsMed <- apply(gammaFitNULLscale2Preds, 2, function(x){quantile(x, 0.5)})
gammaFitNULLscale2PredsLCB <- apply(gammaFitNULLscale2Preds, 2, function(x){quantile(x, 0.025)})
gammaFitNULLscale2PredsUCB <- apply(gammaFitNULLscale2Preds, 2, function(x){quantile(x, 0.975)})

gammaFitNULLscale2predMetrics <- tibble(
  Fit = "gammaFitNULLscale2",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(gammaFitNULLscale2FitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(gammaFitNULLscale2FitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < gammaFitNULLscale2FitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(gammaFitNULLscale2PredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(gammaFitNULLscale2PredsMed - Actual_Yvec)),
  COV_pred = mean(gammaFitNULLscale2PredsLCB < Actual_Yvec & Actual_Yvec < gammaFitNULLscale2PredsUCB)
)
gammaFitNULLscale2predMetrics

##### PPC ----
###### Quantile 2.5 
gammaFitNULLscale2LCBsims <- apply(gammaFitNULLscale2Fit, 
                                   MARGIN = 1,
                                   function(x){
                                     quantile(x, 0.025)
                                   })
gammaFitNULLscale2LCBpvalueVec <- gammaFitNULLscale2LCBsims < quantile(StormdataTrain3$VMAX, 0.025)
gammaFitNULLscale2LCBpvalue <- sum(gammaFitNULLscale2LCBpvalueVec)
gammaFitNULLscale2LCBpvalue <- round(gammaFitNULLscale2LCBpvalue/4000, 3)
gammaFitNULLscale2LCBpvalue <- min(gammaFitNULLscale2LCBpvalue, 1 - gammaFitNULLscale2LCBpvalue)

gammaFitNULLscale2_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLscale2Fit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", gammaFitNULLscale2LCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLscale2_ppcLCB

###### Quantile 97.5 
gammaFitNULLscale2UCBsims <- apply(gammaFitNULLscale2Fit, 
                                   MARGIN = 1,
                                   function(x){
                                     quantile(x, 0.975)
                                   })
gammaFitNULLscale2UCBpvalueVec <- gammaFitNULLscale2UCBsims < quantile(StormdataTrain3$VMAX, 0.975)
gammaFitNULLscale2UCBpvalue <- as.numeric(sum(gammaFitNULLscale2UCBpvalueVec))
gammaFitNULLscale2UCBpvalue <- round(gammaFitNULLscale2UCBpvalue/4000, 3)
gammaFitNULLscale2UCBpvalue <- min(gammaFitNULLscale2UCBpvalue, 1 - gammaFitNULLscale2UCBpvalue)

gammaFitNULLscale2_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLscale2Fit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", gammaFitNULLscale2UCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLscale2_ppcUCB

###### Mean 
gammaFitNULLscale2MEANsims <- apply(gammaFitNULLscale2Fit, 
                                    MARGIN = 1,
                                    function(x){
                                      mean(x)
                                    })
gammaFitNULLscale2MEANpvalueVec <- gammaFitNULLscale2MEANsims < mean(StormdataTrain3$VMAX)
gammaFitNULLscale2MEANpvalue <- sum(gammaFitNULLscale2MEANpvalueVec)
gammaFitNULLscale2MEANpvalue <- round(gammaFitNULLscale2MEANpvalue/4000, 3)
gammaFitNULLscale2MEANpvalue <- min(gammaFitNULLscale2MEANpvalue, 1 - gammaFitNULLscale2MEANpvalue)

gammaFitNULLscale2_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLscale2Fit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", gammaFitNULLscale2MEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLscale2_ppcMEAN

###### Med 
gammaFitNULLscale2MEDsims <- apply(gammaFitNULLscale2Fit, 
                                   MARGIN = 1,
                                   function(x){
                                     quantile(x, 0.5)
                                   })
gammaFitNULLscale2MEDpvalueVec <- gammaFitNULLscale2MEDsims < quantile(StormdataTrain3$VMAX, 0.5)
gammaFitNULLscale2MEDpvalue <- sum(gammaFitNULLscale2MEDpvalueVec)
gammaFitNULLscale2MEDpvalue <- round(gammaFitNULLscale2MEDpvalue/4000, 3)
gammaFitNULLscale2MEDpvalue <- min(gammaFitNULLscale2MEDpvalue, 1 - gammaFitNULLscale2MEDpvalue)

gammaFitNULLscale2_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLscale2Fit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", gammaFitNULLscale2MEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLscale2_ppcMED

###### SD 
gammaFitNULLscale2SDsims <- apply(gammaFitNULLscale2Fit, 
                                  MARGIN = 1,
                                  function(x){
                                    sd(x)
                                  })
gammaFitNULLscale2SDpvalueVec <- gammaFitNULLscale2SDsims < sd(StormdataTrain3$VMAX)
gammaFitNULLscale2SDpvalue <- sum(gammaFitNULLscale2SDpvalueVec)
gammaFitNULLscale2SDpvalue <- round(gammaFitNULLscale2SDpvalue/4000, 3)
gammaFitNULLscale2SDpvalue <- min(gammaFitNULLscale2SDpvalue, 1 - gammaFitNULLscale2SDpvalue)

gammaFitNULLscale2_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLscale2Fit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", gammaFitNULLscale2SDpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLscale2_ppcSD

###### Range 
gammaFitNULLscale2RANGEsims <- apply(gammaFitNULLscale2Fit, 
                                     MARGIN = 1,
                                     function(x){
                                       max(x)-min(x)
                                     })
gammaFitNULLscale2RANGEpvalueVec <- gammaFitNULLscale2RANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
gammaFitNULLscale2RANGEpvalue <- sum(gammaFitNULLscale2RANGEpvalueVec)
gammaFitNULLscale2RANGEpvalue <- round(gammaFitNULLscale2RANGEpvalue/4000, 3)
gammaFitNULLscale2RANGEpvalue <- min(gammaFitNULLscale2RANGEpvalue, 1 - gammaFitNULLscale2RANGEpvalue)

gammaFitNULLscale2_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLscale2Fit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", gammaFitNULLscale2RANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLscale2_ppcRANGE

##### Combined Plot ----
gammaFitNULLscale2_ppcComb <- 
  gammaFitNULLscale2ppcFit /
  (gammaFitNULLscale2_ppcLCB | gammaFitNULLscale2_ppcMED | gammaFitNULLscale2_ppcUCB) /
  (gammaFitNULLscale2_ppcRANGE | gammaFitNULLscale2_ppcMEAN | gammaFitNULLscale2_ppcSD)
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
gammaFitNULLscale2_ppcComb

##### Bayes p-values ----
gammaFitNULLscale2pvalues <- tibble(
  Fit = paste0("gammaFitNULLscale2"),
  LCB = gammaFitNULLscale2LCBpvalue,
  Median = gammaFitNULLscale2MEDpvalue,
  UCB = gammaFitNULLscale2UCBpvalue,
  Range = gammaFitNULLscale2RANGEpvalue,
  Mean = gammaFitNULLscale2MEANpvalue,
  SD = gammaFitNULLscale2SDpvalue
)
gammaFitNULLscale2pvalues

##### CV ----
kfoldID <- kfold_split_grouped(K = 5, StormdataTrain8$StormID)
gammaFitNULLscale2kfoldgroup <- kfold(gammaFitNULLscale2,
                                      folds = kfoldID,
                                      chains = 1,
                                      save_fits = TRUE)
gammaFitNULLscale2kfoldrand <- kfold(gammaFitNULLscale2,
                                     K = 5,
                                     chains = 1,
                                     save_fits = TRUE)
gammaFitNULLscale2kfoldPreds <- kfold_predict(gammaFitNULLscale2kfoldgroup)
#gammaFitNULLscale2kfoldPreds <- kfold_predict(gammaFitNULLscale2kfold)
gammaFitNULLscale2kfoldPredsDat <- gammaFitNULLscale2kfoldPreds$yrep
gammaFitNULLscale2kfoldPredsMean <- colMeans(gammaFitNULLscale2kfoldPredsDat)
gammaFitNULLscale2kfoldPredsMed <- apply(gammaFitNULLscale2kfoldPredsDat, 2, function(x){quantile(x, 0.5)})
gammaFitNULLscale2kfoldPredsLCB <- apply(gammaFitNULLscale2kfoldPredsDat, 2, function(x){quantile(x, 0.025)})
gammaFitNULLscale2kfoldPredsUCB <- apply(gammaFitNULLscale2kfoldPredsDat, 2, function(x){quantile(x, 0.975)})

gammaFitNULLscale2kfoldMetrics <- tibble(
  Fit = paste0("gammaFitNULLscale2"),
  MAE_kfold = mean(abs(gammaFitNULLscale2kfoldPredsMean - gammaFitNULLscale2kfoldPreds$y)),
  MAD_kfold = mean(abs(gammaFitNULLscale2kfoldPredsMed - gammaFitNULLscale2kfoldPreds$y)),
  COV_kfold = mean(gammaFitNULLscale2kfoldPredsLCB < gammaFitNULLscale2kfoldPreds$y & gammaFitNULLscale2kfoldPreds$y < gammaFitNULLscale2kfoldPredsUCB)
)
gammaFitNULLscale2kfoldMetrics

#### Identity link No Int ----
gammaFitNULLscaleidentity <- brm(
  bf(
    VMAX ~ 0 + Intercept + HWRF
  ),
  data = StormdataTrain7scale, 
  family = brmsfamily(family = "Gamma", link = "identity"),
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)

##### Plot ----
print(gammaFitNULLscaleidentity, digits = 4)
gammaFitNULLscaleidentityppcFit <- pp_check(gammaFitNULLscaleidentity, ndraws = 100) + 
  labs(title = "gammaFitNULLscaleidentity Fit PPC") +
  theme_bw()
gammaFitNULLscaleidentityppcFit

##### LOO ----
gammaFitNULLscaleidentityloo <- loo(gammaFitNULLscaleidentity)

##### Prediction ----
## Fitted
gammaFitNULLscaleidentityFit <- posterior_predict(gammaFitNULLscaleidentity)
gammaFitNULLscaleidentityFitMean <- colMeans(gammaFitNULLscaleidentityFit)
gammaFitNULLscaleidentityFitMed <- apply(gammaFitNULLscaleidentityFit, 2, function(x){quantile(x, 0.5)})
gammaFitNULLscaleidentityFitLCB <- apply(gammaFitNULLscaleidentityFit, 2, function(x){quantile(x, 0.025)})
gammaFitNULLscaleidentityFitUCB <- apply(gammaFitNULLscaleidentityFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
gammaFitNULLscaleidentityPreds <- posterior_predict(gammaFitNULLscaleidentity, 
                                                    newdata = StormdataTest7scale,
                                                    allow_new_levels = TRUE)
gammaFitNULLscaleidentityPredsMean <- colMeans(gammaFitNULLscaleidentityPreds)
gammaFitNULLscaleidentityPredsMed <- apply(gammaFitNULLscaleidentityPreds, 2, function(x){quantile(x, 0.5)})
gammaFitNULLscaleidentityPredsLCB <- apply(gammaFitNULLscaleidentityPreds, 2, function(x){quantile(x, 0.025)})
gammaFitNULLscaleidentityPredsUCB <- apply(gammaFitNULLscaleidentityPreds, 2, function(x){quantile(x, 0.975)})

gammaFitNULLscaleidentitypredMetrics <- tibble(
  Fit = "gammaFitNULLscaleidentity",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(gammaFitNULLscaleidentityFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(gammaFitNULLscaleidentityFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < gammaFitNULLscaleidentityFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(gammaFitNULLscaleidentityPredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(gammaFitNULLscaleidentityPredsMed - Actual_Yvec)),
  COV_pred = mean(gammaFitNULLscaleidentityPredsLCB < Actual_Yvec & Actual_Yvec < gammaFitNULLscaleidentityPredsUCB)
)
gammaFitNULLscaleidentitypredMetrics

##### PPC ----
###### Quantile 2.5 
gammaFitNULLscaleidentityLCBsims <- apply(gammaFitNULLscaleidentityFit, 
                                          MARGIN = 1,
                                          function(x){
                                            quantile(x, 0.025)
                                          })
gammaFitNULLscaleidentityLCBpvalueVec <- gammaFitNULLscaleidentityLCBsims < quantile(StormdataTrain3$VMAX, 0.025)
gammaFitNULLscaleidentityLCBpvalue <- sum(gammaFitNULLscaleidentityLCBpvalueVec)
gammaFitNULLscaleidentityLCBpvalue <- round(gammaFitNULLscaleidentityLCBpvalue/4000, 3)
gammaFitNULLscaleidentityLCBpvalue <- min(gammaFitNULLscaleidentityLCBpvalue, 1 - gammaFitNULLscaleidentityLCBpvalue)

gammaFitNULLscaleidentity_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLscaleidentityFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", gammaFitNULLscaleidentityLCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLscaleidentity_ppcLCB

###### Quantile 97.5 
gammaFitNULLscaleidentityUCBsims <- apply(gammaFitNULLscaleidentityFit, 
                                          MARGIN = 1,
                                          function(x){
                                            quantile(x, 0.975)
                                          })
gammaFitNULLscaleidentityUCBpvalueVec <- gammaFitNULLscaleidentityUCBsims < quantile(StormdataTrain3$VMAX, 0.975)
gammaFitNULLscaleidentityUCBpvalue <- as.numeric(sum(gammaFitNULLscaleidentityUCBpvalueVec))
gammaFitNULLscaleidentityUCBpvalue <- round(gammaFitNULLscaleidentityUCBpvalue/4000, 3)
gammaFitNULLscaleidentityUCBpvalue <- min(gammaFitNULLscaleidentityUCBpvalue, 1 - gammaFitNULLscaleidentityUCBpvalue)

gammaFitNULLscaleidentity_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLscaleidentityFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", gammaFitNULLscaleidentityUCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLscaleidentity_ppcUCB

###### Mean 
gammaFitNULLscaleidentityMEANsims <- apply(gammaFitNULLscaleidentityFit, 
                                           MARGIN = 1,
                                           function(x){
                                             mean(x)
                                           })
gammaFitNULLscaleidentityMEANpvalueVec <- gammaFitNULLscaleidentityMEANsims < mean(StormdataTrain3$VMAX)
gammaFitNULLscaleidentityMEANpvalue <- sum(gammaFitNULLscaleidentityMEANpvalueVec)
gammaFitNULLscaleidentityMEANpvalue <- round(gammaFitNULLscaleidentityMEANpvalue/4000, 3)
gammaFitNULLscaleidentityMEANpvalue <- min(gammaFitNULLscaleidentityMEANpvalue, 1 - gammaFitNULLscaleidentityMEANpvalue)

gammaFitNULLscaleidentity_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLscaleidentityFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", gammaFitNULLscaleidentityMEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLscaleidentity_ppcMEAN

###### Med 
gammaFitNULLscaleidentityMEDsims <- apply(gammaFitNULLscaleidentityFit, 
                                          MARGIN = 1,
                                          function(x){
                                            quantile(x, 0.5)
                                          })
gammaFitNULLscaleidentityMEDpvalueVec <- gammaFitNULLscaleidentityMEDsims < quantile(StormdataTrain3$VMAX, 0.5)
gammaFitNULLscaleidentityMEDpvalue <- sum(gammaFitNULLscaleidentityMEDpvalueVec)
gammaFitNULLscaleidentityMEDpvalue <- round(gammaFitNULLscaleidentityMEDpvalue/4000, 3)
gammaFitNULLscaleidentityMEDpvalue <- min(gammaFitNULLscaleidentityMEDpvalue, 1 - gammaFitNULLscaleidentityMEDpvalue)

gammaFitNULLscaleidentity_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLscaleidentityFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", gammaFitNULLscaleidentityMEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLscaleidentity_ppcMED

###### SD 
gammaFitNULLscaleidentitySDsims <- apply(gammaFitNULLscaleidentityFit, 
                                         MARGIN = 1,
                                         function(x){
                                           sd(x)
                                         })
gammaFitNULLscaleidentitySDpvalueVec <- gammaFitNULLscaleidentitySDsims < sd(StormdataTrain3$VMAX)
gammaFitNULLscaleidentitySDpvalue <- sum(gammaFitNULLscaleidentitySDpvalueVec)
gammaFitNULLscaleidentitySDpvalue <- round(gammaFitNULLscaleidentitySDpvalue/4000, 3)
gammaFitNULLscaleidentitySDpvalue <- min(gammaFitNULLscaleidentitySDpvalue, 1 - gammaFitNULLscaleidentitySDpvalue)

gammaFitNULLscaleidentity_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLscaleidentityFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", gammaFitNULLscaleidentitySDpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLscaleidentity_ppcSD

###### Range 
gammaFitNULLscaleidentityRANGEsims <- apply(gammaFitNULLscaleidentityFit, 
                                            MARGIN = 1,
                                            function(x){
                                              max(x)-min(x)
                                            })
gammaFitNULLscaleidentityRANGEpvalueVec <- gammaFitNULLscaleidentityRANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
gammaFitNULLscaleidentityRANGEpvalue <- sum(gammaFitNULLscaleidentityRANGEpvalueVec)
gammaFitNULLscaleidentityRANGEpvalue <- round(gammaFitNULLscaleidentityRANGEpvalue/4000, 3)
gammaFitNULLscaleidentityRANGEpvalue <- min(gammaFitNULLscaleidentityRANGEpvalue, 1 - gammaFitNULLscaleidentityRANGEpvalue)

gammaFitNULLscaleidentity_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLscaleidentityFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", gammaFitNULLscaleidentityRANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLscaleidentity_ppcRANGE

##### Combined Plot ----
gammaFitNULLscaleidentity_ppcComb <- 
  gammaFitNULLscaleidentityppcFit /
  (gammaFitNULLscaleidentity_ppcLCB | gammaFitNULLscaleidentity_ppcMED | gammaFitNULLscaleidentity_ppcUCB) /
  (gammaFitNULLscaleidentity_ppcRANGE | gammaFitNULLscaleidentity_ppcMEAN | gammaFitNULLscaleidentity_ppcSD)
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
gammaFitNULLscaleidentity_ppcComb

##### Bayes p-values ----
gammaFitNULLscaleidentitypvalues <- tibble(
  Fit = paste0("gammaFitNULLscaleidentity"),
  LCB = gammaFitNULLscaleidentityLCBpvalue,
  Median = gammaFitNULLscaleidentityMEDpvalue,
  UCB = gammaFitNULLscaleidentityUCBpvalue,
  Range = gammaFitNULLscaleidentityRANGEpvalue,
  Mean = gammaFitNULLscaleidentityMEANpvalue,
  SD = gammaFitNULLscaleidentitySDpvalue
)
gammaFitNULLscaleidentitypvalues

##### CV ----
kfoldID <- kfold_split_grouped(K = 5, StormdataTrain8$StormID)
gammaFitNULLscaleidentitykfoldgroup <- kfold(gammaFitNULLscaleidentity,
                                             folds = kfoldID,
                                             chains = 1,
                                             save_fits = TRUE)
gammaFitNULLscaleidentitykfoldrand <- kfold(gammaFitNULLscaleidentity,
                                            K = 5,
                                            chains = 1,
                                            save_fits = TRUE)
gammaFitNULLscaleidentitykfoldPreds <- kfold_predict(gammaFitNULLscaleidentitykfoldgroup)
#gammaFitNULLscaleidentitykfoldPreds <- kfold_predict(gammaFitNULLscaleidentitykfold)
gammaFitNULLscaleidentitykfoldPredsDat <- gammaFitNULLscaleidentitykfoldPreds$yrep
gammaFitNULLscaleidentitykfoldPredsMean <- colMeans(gammaFitNULLscaleidentitykfoldPredsDat)
gammaFitNULLscaleidentitykfoldPredsMed <- apply(gammaFitNULLscaleidentitykfoldPredsDat, 2, function(x){quantile(x, 0.5)})
gammaFitNULLscaleidentitykfoldPredsLCB <- apply(gammaFitNULLscaleidentitykfoldPredsDat, 2, function(x){quantile(x, 0.025)})
gammaFitNULLscaleidentitykfoldPredsUCB <- apply(gammaFitNULLscaleidentitykfoldPredsDat, 2, function(x){quantile(x, 0.975)})

gammaFitNULLscaleidentitykfoldMetrics <- tibble(
  Fit = paste0("gammaFitNULLscaleidentity"),
  MAE_kfold = mean(abs(gammaFitNULLscaleidentitykfoldPredsMean - gammaFitNULLscaleidentitykfoldPreds$y)),
  MAD_kfold = mean(abs(gammaFitNULLscaleidentitykfoldPredsMed - gammaFitNULLscaleidentitykfoldPreds$y)),
  COV_kfold = mean(gammaFitNULLscaleidentitykfoldPredsLCB < gammaFitNULLscaleidentitykfoldPreds$y & gammaFitNULLscaleidentitykfoldPreds$y < gammaFitNULLscaleidentitykfoldPredsUCB)
)
gammaFitNULLscaleidentitykfoldMetrics

### Log ----
#### Log Link ----
gammaFitNULLlog <- brm(
  bf(
    VMAX ~ log(HWRF)
  ),
  data = StormdataTrain8, 
  family = Gamma(link = "log"),
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)

##### Plot ----
print(gammaFitNULLlog, digits = 4)
gammaFitNULLlogppcFit <- pp_check(gammaFitNULLlog, ndraws = 100) + 
  labs(title = "gammaFitNULLlog Fit PPC") +
  theme_bw()
gammaFitNULLlogppcFit

##### LOO ----
gammaFitNULLlogloo <- loo(gammaFitNULLlog)

##### Prediction ----
## Fitted
gammaFitNULLlogFit <- posterior_predict(gammaFitNULLlog)
gammaFitNULLlogFitMean <- colMeans(gammaFitNULLlogFit)
gammaFitNULLlogFitMed <- apply(gammaFitNULLlogFit, 2, function(x){quantile(x, 0.5)})
gammaFitNULLlogFitLCB <- apply(gammaFitNULLlogFit, 2, function(x){quantile(x, 0.025)})
gammaFitNULLlogFitUCB <- apply(gammaFitNULLlogFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
gammaFitNULLlogPreds <- posterior_predict(gammaFitNULLlog, 
                                          newdata = StormdataTest8,
                                          allow_new_levels = TRUE)
gammaFitNULLlogPredsMean <- colMeans(gammaFitNULLlogPreds)
gammaFitNULLlogPredsMed <- apply(gammaFitNULLlogPreds, 2, function(x){quantile(x, 0.5)})
gammaFitNULLlogPredsLCB <- apply(gammaFitNULLlogPreds, 2, function(x){quantile(x, 0.025)})
gammaFitNULLlogPredsUCB <- apply(gammaFitNULLlogPreds, 2, function(x){quantile(x, 0.975)})

gammaFitNULLlogpredMetrics <- tibble(
  Fit = "gammaFitNULLlog",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(gammaFitNULLlogFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(gammaFitNULLlogFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < gammaFitNULLlogFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(gammaFitNULLlogPredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(gammaFitNULLlogPredsMed - Actual_Yvec)),
  COV_pred = mean(gammaFitNULLlogPredsLCB < Actual_Yvec & Actual_Yvec < gammaFitNULLlogPredsUCB)
)
gammaFitNULLlogpredMetrics

##### PPC ----
###### Quantile 2.5 
gammaFitNULLlogLCBsims <- apply(gammaFitNULLlogFit, 
                                MARGIN = 1,
                                function(x){
                                  quantile(x, 0.025)
                                })
gammaFitNULLlogLCBpvalueVec <- gammaFitNULLlogLCBsims < quantile(StormdataTrain3$VMAX, 0.025)
gammaFitNULLlogLCBpvalue <- sum(gammaFitNULLlogLCBpvalueVec)
gammaFitNULLlogLCBpvalue <- round(gammaFitNULLlogLCBpvalue/4000, 3)
gammaFitNULLlogLCBpvalue <- min(gammaFitNULLlogLCBpvalue, 1 - gammaFitNULLlogLCBpvalue)

gammaFitNULLlog_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLlogFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", gammaFitNULLlogLCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLlog_ppcLCB

###### Quantile 97.5 
gammaFitNULLlogUCBsims <- apply(gammaFitNULLlogFit, 
                                MARGIN = 1,
                                function(x){
                                  quantile(x, 0.975)
                                })
gammaFitNULLlogUCBpvalueVec <- gammaFitNULLlogUCBsims < quantile(StormdataTrain3$VMAX, 0.975)
gammaFitNULLlogUCBpvalue <- as.numeric(sum(gammaFitNULLlogUCBpvalueVec))
gammaFitNULLlogUCBpvalue <- round(gammaFitNULLlogUCBpvalue/4000, 3)
gammaFitNULLlogUCBpvalue <- min(gammaFitNULLlogUCBpvalue, 1 - gammaFitNULLlogUCBpvalue)

gammaFitNULLlog_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLlogFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", gammaFitNULLlogUCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLlog_ppcUCB

###### Mean 
gammaFitNULLlogMEANsims <- apply(gammaFitNULLlogFit, 
                                 MARGIN = 1,
                                 function(x){
                                   mean(x)
                                 })
gammaFitNULLlogMEANpvalueVec <- gammaFitNULLlogMEANsims < mean(StormdataTrain3$VMAX)
gammaFitNULLlogMEANpvalue <- sum(gammaFitNULLlogMEANpvalueVec)
gammaFitNULLlogMEANpvalue <- round(gammaFitNULLlogMEANpvalue/4000, 3)
gammaFitNULLlogMEANpvalue <- min(gammaFitNULLlogMEANpvalue, 1 - gammaFitNULLlogMEANpvalue)

gammaFitNULLlog_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLlogFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", gammaFitNULLlogMEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLlog_ppcMEAN

###### Med 
gammaFitNULLlogMEDsims <- apply(gammaFitNULLlogFit, 
                                MARGIN = 1,
                                function(x){
                                  quantile(x, 0.5)
                                })
gammaFitNULLlogMEDpvalueVec <- gammaFitNULLlogMEDsims < quantile(StormdataTrain3$VMAX, 0.5)
gammaFitNULLlogMEDpvalue <- sum(gammaFitNULLlogMEDpvalueVec)
gammaFitNULLlogMEDpvalue <- round(gammaFitNULLlogMEDpvalue/4000, 3)
gammaFitNULLlogMEDpvalue <- min(gammaFitNULLlogMEDpvalue, 1 - gammaFitNULLlogMEDpvalue)

gammaFitNULLlog_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLlogFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", gammaFitNULLlogMEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLlog_ppcMED

###### SD 
gammaFitNULLlogSDsims <- apply(gammaFitNULLlogFit, 
                               MARGIN = 1,
                               function(x){
                                 sd(x)
                               })
gammaFitNULLlogSDpvalueVec <- gammaFitNULLlogSDsims < sd(StormdataTrain3$VMAX)
gammaFitNULLlogSDpvalue <- sum(gammaFitNULLlogSDpvalueVec)
gammaFitNULLlogSDpvalue <- round(gammaFitNULLlogSDpvalue/4000, 3)
gammaFitNULLlogSDpvalue <- min(gammaFitNULLlogSDpvalue, 1 - gammaFitNULLlogSDpvalue)

gammaFitNULLlog_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLlogFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", gammaFitNULLlogSDpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLlog_ppcSD

###### Range 
gammaFitNULLlogRANGEsims <- apply(gammaFitNULLlogFit, 
                                  MARGIN = 1,
                                  function(x){
                                    max(x)-min(x)
                                  })
gammaFitNULLlogRANGEpvalueVec <- gammaFitNULLlogRANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
gammaFitNULLlogRANGEpvalue <- sum(gammaFitNULLlogRANGEpvalueVec)
gammaFitNULLlogRANGEpvalue <- round(gammaFitNULLlogRANGEpvalue/4000, 3)
gammaFitNULLlogRANGEpvalue <- min(gammaFitNULLlogRANGEpvalue, 1 - gammaFitNULLlogRANGEpvalue)

gammaFitNULLlog_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLlogFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", gammaFitNULLlogRANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLlog_ppcRANGE

##### Combined Plot ----
gammaFitNULLlog_ppcComb <- 
  gammaFitNULLlogppcFit /
  (gammaFitNULLlog_ppcLCB | gammaFitNULLlog_ppcMED | gammaFitNULLlog_ppcUCB) /
  (gammaFitNULLlog_ppcRANGE | gammaFitNULLlog_ppcMEAN | gammaFitNULLlog_ppcSD)
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
gammaFitNULLlog_ppcComb

##### Bayes p-values ----
gammaFitNULLlogpvalues <- tibble(
  Fit = paste0("gammaFitNULLlog"),
  LCB = gammaFitNULLlogLCBpvalue,
  Median = gammaFitNULLlogMEDpvalue,
  UCB = gammaFitNULLlogUCBpvalue,
  Range = gammaFitNULLlogRANGEpvalue,
  Mean = gammaFitNULLlogMEANpvalue,
  SD = gammaFitNULLlogSDpvalue
)
gammaFitNULLlogpvalues

##### CV ----
kfoldID <- kfold_split_grouped(K = 5, StormdataTrain8$StormID)
gammaFitNULLlogkfoldgroup <- kfold(gammaFitNULLlog,
                                   folds = kfoldID,
                                   chains = 1,
                                   save_fits = TRUE)
gammaFitNULLlogkfoldrand <- kfold(gammaFitNULLlog,
                                  K = 5,
                                  chains = 1,
                                  save_fits = TRUE)
gammaFitNULLlogkfoldPreds <- kfold_predict(gammaFitNULLlogkfoldgroup)
#gammaFitNULLlogkfoldPreds <- kfold_predict(gammaFitNULLlogkfold)
gammaFitNULLlogkfoldPredsDat <- gammaFitNULLlogkfoldPreds$yrep
gammaFitNULLlogkfoldPredsMean <- colMeans(gammaFitNULLlogkfoldPredsDat)
gammaFitNULLlogkfoldPredsMed <- apply(gammaFitNULLlogkfoldPredsDat, 2, function(x){quantile(x, 0.5)})
gammaFitNULLlogkfoldPredsLCB <- apply(gammaFitNULLlogkfoldPredsDat, 2, function(x){quantile(x, 0.025)})
gammaFitNULLlogkfoldPredsUCB <- apply(gammaFitNULLlogkfoldPredsDat, 2, function(x){quantile(x, 0.975)})

gammaFitNULLlogkfoldMetrics <- tibble(
  Fit = paste0("gammaFitNULLlog"),
  MAE_kfold = mean(abs(gammaFitNULLlogkfoldPredsMean - gammaFitNULLlogkfoldPreds$y)),
  MAD_kfold = mean(abs(gammaFitNULLlogkfoldPredsMed - gammaFitNULLlogkfoldPreds$y)),
  COV_kfold = mean(gammaFitNULLlogkfoldPredsLCB < gammaFitNULLlogkfoldPreds$y & gammaFitNULLlogkfoldPreds$y < gammaFitNULLlogkfoldPredsUCB)
)
gammaFitNULLlogkfoldMetrics

#### Log Link No Int ----
gammaFitNULLlog2 <- brm(
  bf(
    VMAX ~ 0 + Intercept + log(HWRF)
  ),
  data = StormdataTrain8, 
  family = Gamma(link = "log"),
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)

##### Plot ----
print(gammaFitNULLlog2, digits = 4)
gammaFitNULLlog2ppcFit <- pp_check(gammaFitNULLlog2, ndraws = 100) + 
  labs(title = "gammaFitNULLlog2 Fit PPC") +
  theme_bw()
#gammaFitNULLlog2ppcFit

##### LOO ----
gammaFitNULLlog2loo <- loo(gammaFitNULLlog2)

##### Prediction ----
## Fitted
gammaFitNULLlog2Fit <- posterior_predict(gammaFitNULLlog2)
gammaFitNULLlog2FitMean <- colMeans(gammaFitNULLlog2Fit)
gammaFitNULLlog2FitMed <- apply(gammaFitNULLlog2Fit, 2, function(x){quantile(x, 0.5)})
gammaFitNULLlog2FitLCB <- apply(gammaFitNULLlog2Fit, 2, function(x){quantile(x, 0.025)})
gammaFitNULLlog2FitUCB <- apply(gammaFitNULLlog2Fit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
gammaFitNULLlog2Preds <- posterior_predict(gammaFitNULLlog2, 
                                           newdata = StormdataTest8,
                                           allow_new_levels = TRUE)
gammaFitNULLlog2PredsMean <- colMeans(gammaFitNULLlog2Preds)
gammaFitNULLlog2PredsMed <- apply(gammaFitNULLlog2Preds, 2, function(x){quantile(x, 0.5)})
gammaFitNULLlog2PredsLCB <- apply(gammaFitNULLlog2Preds, 2, function(x){quantile(x, 0.025)})
gammaFitNULLlog2PredsUCB <- apply(gammaFitNULLlog2Preds, 2, function(x){quantile(x, 0.975)})

gammaFitNULLlog2predMetrics <- tibble(
  Fit = "gammaFitNULLlog2",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(gammaFitNULLlog2FitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(gammaFitNULLlog2FitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < gammaFitNULLlog2FitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(gammaFitNULLlog2PredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(gammaFitNULLlog2PredsMed - Actual_Yvec)),
  COV_pred = mean(gammaFitNULLlog2PredsLCB < Actual_Yvec & Actual_Yvec < gammaFitNULLlog2PredsUCB)
)
gammaFitNULLlog2predMetrics

##### PPC ----
###### Quantile 2.5 
gammaFitNULLlog2LCBsims <- apply(gammaFitNULLlog2Fit, 
                                 MARGIN = 1,
                                 function(x){
                                   quantile(x, 0.025)
                                 })
gammaFitNULLlog2LCBpvalueVec <- gammaFitNULLlog2LCBsims < quantile(StormdataTrain3$VMAX, 0.025)
gammaFitNULLlog2LCBpvalue <- sum(gammaFitNULLlog2LCBpvalueVec)
gammaFitNULLlog2LCBpvalue <- round(gammaFitNULLlog2LCBpvalue/4000, 3)
gammaFitNULLlog2LCBpvalue <- min(gammaFitNULLlog2LCBpvalue, 1 - gammaFitNULLlog2LCBpvalue)

gammaFitNULLlog2_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLlog2Fit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", gammaFitNULLlog2LCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLlog2_ppcLCB

###### Quantile 97.5 
gammaFitNULLlog2UCBsims <- apply(gammaFitNULLlog2Fit, 
                                 MARGIN = 1,
                                 function(x){
                                   quantile(x, 0.975)
                                 })
gammaFitNULLlog2UCBpvalueVec <- gammaFitNULLlog2UCBsims < quantile(StormdataTrain3$VMAX, 0.975)
gammaFitNULLlog2UCBpvalue <- as.numeric(sum(gammaFitNULLlog2UCBpvalueVec))
gammaFitNULLlog2UCBpvalue <- round(gammaFitNULLlog2UCBpvalue/4000, 3)
gammaFitNULLlog2UCBpvalue <- min(gammaFitNULLlog2UCBpvalue, 1 - gammaFitNULLlog2UCBpvalue)

gammaFitNULLlog2_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLlog2Fit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", gammaFitNULLlog2UCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLlog2_ppcUCB

###### Mean 
gammaFitNULLlog2MEANsims <- apply(gammaFitNULLlog2Fit, 
                                  MARGIN = 1,
                                  function(x){
                                    mean(x)
                                  })
gammaFitNULLlog2MEANpvalueVec <- gammaFitNULLlog2MEANsims < mean(StormdataTrain3$VMAX)
gammaFitNULLlog2MEANpvalue <- sum(gammaFitNULLlog2MEANpvalueVec)
gammaFitNULLlog2MEANpvalue <- round(gammaFitNULLlog2MEANpvalue/4000, 3)
gammaFitNULLlog2MEANpvalue <- min(gammaFitNULLlog2MEANpvalue, 1 - gammaFitNULLlog2MEANpvalue)

gammaFitNULLlog2_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLlog2Fit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", gammaFitNULLlog2MEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLlog2_ppcMEAN

###### Med 
gammaFitNULLlog2MEDsims <- apply(gammaFitNULLlog2Fit, 
                                 MARGIN = 1,
                                 function(x){
                                   quantile(x, 0.5)
                                 })
gammaFitNULLlog2MEDpvalueVec <- gammaFitNULLlog2MEDsims < quantile(StormdataTrain3$VMAX, 0.5)
gammaFitNULLlog2MEDpvalue <- sum(gammaFitNULLlog2MEDpvalueVec)
gammaFitNULLlog2MEDpvalue <- round(gammaFitNULLlog2MEDpvalue/4000, 3)
gammaFitNULLlog2MEDpvalue <- min(gammaFitNULLlog2MEDpvalue, 1 - gammaFitNULLlog2MEDpvalue)

gammaFitNULLlog2_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLlog2Fit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", gammaFitNULLlog2MEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLlog2_ppcMED

###### SD 
gammaFitNULLlog2SDsims <- apply(gammaFitNULLlog2Fit, 
                                MARGIN = 1,
                                function(x){
                                  sd(x)
                                })
gammaFitNULLlog2SDpvalueVec <- gammaFitNULLlog2SDsims < sd(StormdataTrain3$VMAX)
gammaFitNULLlog2SDpvalue <- sum(gammaFitNULLlog2SDpvalueVec)
gammaFitNULLlog2SDpvalue <- round(gammaFitNULLlog2SDpvalue/4000, 3)
gammaFitNULLlog2SDpvalue <- min(gammaFitNULLlog2SDpvalue, 1 - gammaFitNULLlog2SDpvalue)

gammaFitNULLlog2_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLlog2Fit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", gammaFitNULLlog2SDpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLlog2_ppcSD

###### Range 
gammaFitNULLlog2RANGEsims <- apply(gammaFitNULLlog2Fit, 
                                   MARGIN = 1,
                                   function(x){
                                     max(x)-min(x)
                                   })
gammaFitNULLlog2RANGEpvalueVec <- gammaFitNULLlog2RANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
gammaFitNULLlog2RANGEpvalue <- sum(gammaFitNULLlog2RANGEpvalueVec)
gammaFitNULLlog2RANGEpvalue <- round(gammaFitNULLlog2RANGEpvalue/4000, 3)
gammaFitNULLlog2RANGEpvalue <- min(gammaFitNULLlog2RANGEpvalue, 1 - gammaFitNULLlog2RANGEpvalue)

gammaFitNULLlog2_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLlog2Fit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", gammaFitNULLlog2RANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLlog2_ppcRANGE

##### Combined Plot ----
gammaFitNULLlog2_ppcComb <- 
  gammaFitNULLlog2ppcFit /
  (gammaFitNULLlog2_ppcLCB | gammaFitNULLlog2_ppcMED | gammaFitNULLlog2_ppcUCB) /
  (gammaFitNULLlog2_ppcRANGE | gammaFitNULLlog2_ppcMEAN | gammaFitNULLlog2_ppcSD)
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
gammaFitNULLlog2_ppcComb

##### Bayes p-values ----
gammaFitNULLlog2pvalues <- tibble(
  Fit = paste0("gammaFitNULLlog2"),
  LCB = gammaFitNULLlog2LCBpvalue,
  Median = gammaFitNULLlog2MEDpvalue,
  UCB = gammaFitNULLlog2UCBpvalue,
  Range = gammaFitNULLlog2RANGEpvalue,
  Mean = gammaFitNULLlog2MEANpvalue,
  SD = gammaFitNULLlog2SDpvalue
)
gammaFitNULLlog2pvalues

##### CV ----
kfoldID <- kfold_split_grouped(K = 5, StormdataTrain8$StormID)
gammaFitNULLlog2kfoldgroup <- kfold(gammaFitNULLlog2,
                                    folds = kfoldID,
                                    chains = 1,
                                    save_fits = TRUE)
gammaFitNULLlog2kfoldrand <- kfold(gammaFitNULLlog2,
                                   K = 5,
                                   chains = 1,
                                   save_fits = TRUE)
gammaFitNULLlog2kfoldPreds <- kfold_predict(gammaFitNULLlog2kfoldgroup)
#gammaFitNULLlog2kfoldPreds <- kfold_predict(gammaFitNULLlog2kfold)
gammaFitNULLlog2kfoldPredsDat <- gammaFitNULLlog2kfoldPreds$yrep
gammaFitNULLlog2kfoldPredsMean <- colMeans(gammaFitNULLlog2kfoldPredsDat)
gammaFitNULLlog2kfoldPredsMed <- apply(gammaFitNULLlog2kfoldPredsDat, 2, function(x){quantile(x, 0.5)})
gammaFitNULLlog2kfoldPredsLCB <- apply(gammaFitNULLlog2kfoldPredsDat, 2, function(x){quantile(x, 0.025)})
gammaFitNULLlog2kfoldPredsUCB <- apply(gammaFitNULLlog2kfoldPredsDat, 2, function(x){quantile(x, 0.975)})

gammaFitNULLlog2kfoldMetrics <- tibble(
  Fit = paste0("gammaFitNULLlog2"),
  MAE_kfold = mean(abs(gammaFitNULLlog2kfoldPredsMean - gammaFitNULLlog2kfoldPreds$y)),
  MAD_kfold = mean(abs(gammaFitNULLlog2kfoldPredsMed - gammaFitNULLlog2kfoldPreds$y)),
  COV_kfold = mean(gammaFitNULLlog2kfoldPredsLCB < gammaFitNULLlog2kfoldPreds$y & gammaFitNULLlog2kfoldPreds$y < gammaFitNULLlog2kfoldPredsUCB)
)
gammaFitNULLlog2kfoldMetrics

#### Identity link No Int ----
gammaFitNULLlogidentity <- brm(
  bf(
    VMAX ~ 0 + Intercept + log(HWRF)
  ),
  data = StormdataTrain7scale, 
  family = brmsfamily(family = "Gamma", link = "identity"),
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)

##### Plot ----
print(gammaFitNULLlogidentity, digits = 4)
gammaFitNULLlogidentityppcFit <- pp_check(gammaFitNULLlogidentity, ndraws = 100) + 
  labs(title = "gammaFitNULLlogidentity Fit PPC") +
  theme_bw()
gammaFitNULLlogidentityppcFit

##### LOO ----
gammaFitNULLlogidentityloo <- loo(gammaFitNULLlogidentity)

##### Prediction ----
## Fitted
gammaFitNULLlogidentityFit <- posterior_predict(gammaFitNULLlogidentity)
gammaFitNULLlogidentityFitMean <- colMeans(gammaFitNULLlogidentityFit)
gammaFitNULLlogidentityFitMed <- apply(gammaFitNULLlogidentityFit, 2, function(x){quantile(x, 0.5)})
gammaFitNULLlogidentityFitLCB <- apply(gammaFitNULLlogidentityFit, 2, function(x){quantile(x, 0.025)})
gammaFitNULLlogidentityFitUCB <- apply(gammaFitNULLlogidentityFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
gammaFitNULLlogidentityPreds <- posterior_predict(gammaFitNULLlogidentity, 
                                                  newdata = StormdataTest8,
                                                  allow_new_levels = TRUE)
gammaFitNULLlogidentityPredsMean <- colMeans(gammaFitNULLlogidentityPreds)
gammaFitNULLlogidentityPredsMed <- apply(gammaFitNULLlogidentityPreds, 2, function(x){quantile(x, 0.5)})
gammaFitNULLlogidentityPredsLCB <- apply(gammaFitNULLlogidentityPreds, 2, function(x){quantile(x, 0.025)})
gammaFitNULLlogidentityPredsUCB <- apply(gammaFitNULLlogidentityPreds, 2, function(x){quantile(x, 0.975)})

gammaFitNULLlogidentitypredMetrics <- tibble(
  Fit = "gammaFitNULLlogidentity",
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(gammaFitNULLlogidentityFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(gammaFitNULLlogidentityFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < gammaFitNULLlogidentityFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(gammaFitNULLlogidentityPredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(gammaFitNULLlogidentityPredsMed - Actual_Yvec)),
  COV_pred = mean(gammaFitNULLlogidentityPredsLCB < Actual_Yvec & Actual_Yvec < gammaFitNULLlogidentityPredsUCB)
)
gammaFitNULLlogidentitypredMetrics

##### PPC ----
###### Quantile 2.5 
gammaFitNULLlogidentityLCBsims <- apply(gammaFitNULLlogidentityFit, 
                                        MARGIN = 1,
                                        function(x){
                                          quantile(x, 0.025)
                                        })
gammaFitNULLlogidentityLCBpvalueVec <- gammaFitNULLlogidentityLCBsims < quantile(StormdataTrain3$VMAX, 0.025)
gammaFitNULLlogidentityLCBpvalue <- sum(gammaFitNULLlogidentityLCBpvalueVec)
gammaFitNULLlogidentityLCBpvalue <- round(gammaFitNULLlogidentityLCBpvalue/4000, 3)
gammaFitNULLlogidentityLCBpvalue <- min(gammaFitNULLlogidentityLCBpvalue, 1 - gammaFitNULLlogidentityLCBpvalue)

gammaFitNULLlogidentity_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLlogidentityFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", gammaFitNULLlogidentityLCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLlogidentity_ppcLCB

###### Quantile 97.5 
gammaFitNULLlogidentityUCBsims <- apply(gammaFitNULLlogidentityFit, 
                                        MARGIN = 1,
                                        function(x){
                                          quantile(x, 0.975)
                                        })
gammaFitNULLlogidentityUCBpvalueVec <- gammaFitNULLlogidentityUCBsims < quantile(StormdataTrain3$VMAX, 0.975)
gammaFitNULLlogidentityUCBpvalue <- as.numeric(sum(gammaFitNULLlogidentityUCBpvalueVec))
gammaFitNULLlogidentityUCBpvalue <- round(gammaFitNULLlogidentityUCBpvalue/4000, 3)
gammaFitNULLlogidentityUCBpvalue <- min(gammaFitNULLlogidentityUCBpvalue, 1 - gammaFitNULLlogidentityUCBpvalue)

gammaFitNULLlogidentity_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLlogidentityFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", gammaFitNULLlogidentityUCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLlogidentity_ppcUCB

###### Mean 
gammaFitNULLlogidentityMEANsims <- apply(gammaFitNULLlogidentityFit, 
                                         MARGIN = 1,
                                         function(x){
                                           mean(x)
                                         })
gammaFitNULLlogidentityMEANpvalueVec <- gammaFitNULLlogidentityMEANsims < mean(StormdataTrain3$VMAX)
gammaFitNULLlogidentityMEANpvalue <- sum(gammaFitNULLlogidentityMEANpvalueVec)
gammaFitNULLlogidentityMEANpvalue <- round(gammaFitNULLlogidentityMEANpvalue/4000, 3)
gammaFitNULLlogidentityMEANpvalue <- min(gammaFitNULLlogidentityMEANpvalue, 1 - gammaFitNULLlogidentityMEANpvalue)

gammaFitNULLlogidentity_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLlogidentityFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", gammaFitNULLlogidentityMEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLlogidentity_ppcMEAN

###### Med 
gammaFitNULLlogidentityMEDsims <- apply(gammaFitNULLlogidentityFit, 
                                        MARGIN = 1,
                                        function(x){
                                          quantile(x, 0.5)
                                        })
gammaFitNULLlogidentityMEDpvalueVec <- gammaFitNULLlogidentityMEDsims < quantile(StormdataTrain3$VMAX, 0.5)
gammaFitNULLlogidentityMEDpvalue <- sum(gammaFitNULLlogidentityMEDpvalueVec)
gammaFitNULLlogidentityMEDpvalue <- round(gammaFitNULLlogidentityMEDpvalue/4000, 3)
gammaFitNULLlogidentityMEDpvalue <- min(gammaFitNULLlogidentityMEDpvalue, 1 - gammaFitNULLlogidentityMEDpvalue)

gammaFitNULLlogidentity_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLlogidentityFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", gammaFitNULLlogidentityMEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLlogidentity_ppcMED

###### SD 
gammaFitNULLlogidentitySDsims <- apply(gammaFitNULLlogidentityFit, 
                                       MARGIN = 1,
                                       function(x){
                                         sd(x)
                                       })
gammaFitNULLlogidentitySDpvalueVec <- gammaFitNULLlogidentitySDsims < sd(StormdataTrain3$VMAX)
gammaFitNULLlogidentitySDpvalue <- sum(gammaFitNULLlogidentitySDpvalueVec)
gammaFitNULLlogidentitySDpvalue <- round(gammaFitNULLlogidentitySDpvalue/4000, 3)
gammaFitNULLlogidentitySDpvalue <- min(gammaFitNULLlogidentitySDpvalue, 1 - gammaFitNULLlogidentitySDpvalue)

gammaFitNULLlogidentity_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLlogidentityFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", gammaFitNULLlogidentitySDpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLlogidentity_ppcSD

###### Range 
gammaFitNULLlogidentityRANGEsims <- apply(gammaFitNULLlogidentityFit, 
                                          MARGIN = 1,
                                          function(x){
                                            max(x)-min(x)
                                          })
gammaFitNULLlogidentityRANGEpvalueVec <- gammaFitNULLlogidentityRANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
gammaFitNULLlogidentityRANGEpvalue <- sum(gammaFitNULLlogidentityRANGEpvalueVec)
gammaFitNULLlogidentityRANGEpvalue <- round(gammaFitNULLlogidentityRANGEpvalue/4000, 3)
gammaFitNULLlogidentityRANGEpvalue <- min(gammaFitNULLlogidentityRANGEpvalue, 1 - gammaFitNULLlogidentityRANGEpvalue)

gammaFitNULLlogidentity_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitNULLlogidentityFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", gammaFitNULLlogidentityRANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFitNULLlogidentity_ppcRANGE

##### Combined Plot ----
gammaFitNULLlogidentity_ppcComb <- 
  gammaFitNULLlogidentityppcFit /
  (gammaFitNULLlogidentity_ppcLCB | gammaFitNULLlogidentity_ppcMED | gammaFitNULLlogidentity_ppcUCB) /
  (gammaFitNULLlogidentity_ppcRANGE | gammaFitNULLlogidentity_ppcMEAN | gammaFitNULLlogidentity_ppcSD)
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
gammaFitNULLlogidentity_ppcComb

##### Bayes p-values ----
gammaFitNULLlogidentitypvalues <- tibble(
  Fit = paste0("gammaFitNULLlogidentity"),
  LCB = gammaFitNULLlogidentityLCBpvalue,
  Median = gammaFitNULLlogidentityMEDpvalue,
  UCB = gammaFitNULLlogidentityUCBpvalue,
  Range = gammaFitNULLlogidentityRANGEpvalue,
  Mean = gammaFitNULLlogidentityMEANpvalue,
  SD = gammaFitNULLlogidentitySDpvalue
)
gammaFitNULLlogidentitypvalues

##### CV ----
kfoldID <- kfold_split_grouped(K = 5, StormdataTrain8$StormID)
gammaFitNULLlogidentitykfoldgroup <- kfold(gammaFitNULLlogidentity,
                                           folds = kfoldID,
                                           chains = 1,
                                           save_fits = TRUE)
gammaFitNULLlogidentitykfoldrand <- kfold(gammaFitNULLlogidentity,
                                          K = 5,
                                          chains = 1,
                                          save_fits = TRUE)
gammaFitNULLlogidentitykfoldPreds <- kfold_predict(gammaFitNULLlogidentitykfoldgroup)
#gammaFitNULLlogidentitykfoldPreds <- kfold_predict(gammaFitNULLlogidentitykfold)
gammaFitNULLlogidentitykfoldPredsDat <- gammaFitNULLlogidentitykfoldPreds$yrep
gammaFitNULLlogidentitykfoldPredsMean <- colMeans(gammaFitNULLlogidentitykfoldPredsDat)
gammaFitNULLlogidentitykfoldPredsMed <- apply(gammaFitNULLlogidentitykfoldPredsDat, 2, function(x){quantile(x, 0.5)})
gammaFitNULLlogidentitykfoldPredsLCB <- apply(gammaFitNULLlogidentitykfoldPredsDat, 2, function(x){quantile(x, 0.025)})
gammaFitNULLlogidentitykfoldPredsUCB <- apply(gammaFitNULLlogidentitykfoldPredsDat, 2, function(x){quantile(x, 0.975)})

gammaFitNULLlogidentitykfoldMetrics <- tibble(
  Fit = paste0("gammaFitNULLlogidentity"),
  MAE_kfold = mean(abs(gammaFitNULLlogidentitykfoldPredsMean - gammaFitNULLlogidentitykfoldPreds$y)),
  MAD_kfold = mean(abs(gammaFitNULLlogidentitykfoldPredsMed - gammaFitNULLlogidentitykfoldPreds$y)),
  COV_kfold = mean(gammaFitNULLlogidentitykfoldPredsLCB < gammaFitNULLlogidentitykfoldPreds$y & gammaFitNULLlogidentitykfoldPreds$y < gammaFitNULLlogidentitykfoldPredsUCB)
)
gammaFitNULLlogidentitykfoldMetrics

### Compare NULLs ----
gammaFitNULLppcFitComb <- 
  gammaFitNULLppcFit /
  gammaFitNULL2ppcFit /
  gammaFitNULLidentityppcFit /
  gammaFitNULLscaleppcFit /
  gammaFitNULLscale2ppcFit /
  gammaFitNULLscaleidentityppcFit /
  gammaFitNULLlogppcFit /
  gammaFitNULLlog2ppcFit

(gammaFitNULLidentity_ppcLCB | 
    gammaFitNULLidentity_ppcMED | 
    gammaFitNULLidentity_ppcUCB | 
    gammaFitNULLidentity_ppcRANGE |
    gammaFitNULLidentity_ppcMEAN | 
    gammaFitNULLidentity_ppcSD) /
  (gammaFitNULLlog2_ppcLCB | 
     gammaFitNULLlog2_ppcMED | 
     gammaFitNULLlog2_ppcUCB | 
     gammaFitNULLlog2_ppcRANGE | 
     gammaFitNULLlog2_ppcMEAN | 
     gammaFitNULLlog2_ppcSD) /
  (gammaFitNULLscaleidentity_ppcLCB | 
     gammaFitNULLscaleidentity_ppcMED | 
     gammaFitNULLscaleidentity_ppcUCB | 
     gammaFitNULLscaleidentity_ppcRANGE | 
     gammaFitNULLscaleidentity_ppcMEAN | 
     gammaFitNULLscaleidentity_ppcSD) /
  (gammaFitNULLlog_ppcLCB | 
     gammaFitNULLlog_ppcMED | 
     gammaFitNULLlog_ppcUCB | 
     gammaFitNULLlog_ppcRANGE | 
     gammaFitNULLlog_ppcMEAN | 
     gammaFitNULLlog_ppcSD) 

gammaFitNULLpvaluesComb <- bind_rows(
  gammaFitNULLpvalues,
  gammaFitNULL2pvalues,
  gammaFitNULLidentitypvalues,
  gammaFitNULLscalepvalues,
  gammaFitNULLscale2pvalues,
  gammaFitNULLscaleidentitypvalues,
  gammaFitNULLlogpvalues,
  gammaFitNULLlog2pvalues
)

gammaFitNULLlooComp <- loo_compare(gammaFitNULLloo, 
                                   gammaFitNULL2loo,
                                   gammaFitNULLidentityloo,
                                   gammaFitNULLscaleloo, 
                                   gammaFitNULLscale2loo,
                                   gammaFitNULLscaleidentityloo,
                                   gammaFitNULLlogloo,
                                   gammaFitNULLlog2loo
)

gammaFitNULLpredComp <- bind_rows(
  gammaFitNULLpredMetrics,
  gammaFitNULL2predMetrics,
  gammaFitNULLidentitypredMetrics,
  gammaFitNULLscalepredMetrics,
  gammaFitNULLscale2predMetrics,
  gammaFitNULLscaleidentitypredMetrics,
  gammaFitNULLlogpredMetrics,
  gammaFitNULLlog2predMetrics
) |>
  arrange(MAE_pred)

gammaFitNULLkfoldComp <- bind_rows(
  gammaFitNULLkfoldMetrics,
  gammaFitNULL2kfoldMetrics,
  gammaFitNULLidentitykfoldMetrics,
  gammaFitNULLscalekfoldMetrics,
  gammaFitNULLscale2kfoldMetrics,
  gammaFitNULLscaleidentitykfoldMetrics,
  gammaFitNULLlogkfoldMetrics,
  gammaFitNULLlog2kfoldMetrics
) |>
  arrange(MAE_kfold)

gammaFitNULLidentity
gammaFitNULLlog2
gammaFitNULLscaleidentity
gammaFitNULLlog

### Clear Environment ----
NUllenv <- ls()
exNULL <- str_detect(NUllenv, "NULL")
keepNULL <- NUllenv %in% c("gammaFitNULLppcFitComb",
                        "gammaFitNULLpvaluesComb",
                        "gammaFitNULLlooComp",
                        "gammaFitNULLpredComp",
                        "gammaFitNULLkfoldComp",
                        "gammaFitNULLidentity",
                        "gammaFitNULLlog2",
                        "gammaFitNULLscaleidentity",
                        "gammaFitNULLlog")
cutNULL <- exNULL & !keepNULL
NULLrm <- NUllenv[cutNULL]
rm(list = NULLrm)

## gammaFit1 <- gammaFit
## gammaFit2 <- gammaFit

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
fit <- 3
tic()
gammaFit <- brm(
  bf(VMAX ~ 
       #0 + Intercept +
       s(Day, bs = "cc") +
       s(StormElapsedTime, bs = "tp") +
       s(LON, bs = "tp") +
       s(LAT, bs = "tp") +
       # Day +
       # StormElapsedTime +
       # LON +
       # LAT +
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
  family = Gamma(link = "log"),
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

gammaFit2 <- gammaFit
# save(gammaFit1,
#      file = "~/Desktop/Temp Hurricane Model Data/gammaFit1.RData")
prior_summary(gammaFit)
posterior_summary(gammaFit)
gammaFit
launch_shinystan(gammaFit)

print(gammaFit, digits = 4)
print(gammaFit1, digits = 4)
print(gammaFit2, digits = 4)
print(gammaFit3, digits = 4)
print(gammaFit4, digits = 4)
print(gammaFit5, digits = 4)
print(gammaFit6, digits = 4)
print(gammaFit7, digits = 4)
print(gammaFit8, digits = 4)
print(gammaFit9, digits = 4)
print(gammaFit10, digits = 4)
plot(gammaFit)
gammaFitppcFit <- pp_check(gammaFit, ndraws = 100) + 
  labs(title = paste0("gammaFit", fit, " Fit PPC")) +
  theme_bw()
gammaFitppcFit
gammaFitloo <- loo(gammaFit)
waic(gammaFit)
performance::check_distribution(gammaFit)
performance::check_outliers(gammaFit)
performance::check_heteroskedasticity(gammaFit)
performance_rmse(gammaFit)
performance_mae(gammaFit)
mean(abs(StormdataTrain3$VMAX - StormdataTrain3$HWRF))
model_performance(gammaFit)

variance_decomposition(gammaFit)
fixef(gammaFit)
ranef(gammaFit)

bayes_R2(gammaFit)

gammaFitC <- gammaFit
gammaFit <- gammaFitB
bayes_factor(gammaFit, gammaFit1)
bayes_factor(gammaFit, gammaFit2)
bayes_factor(gammaFit, gammaFit3)
bayes_factor(gammaFit, gammaFit4)
bayes_factor(gammaFit, gammaFit5)
bayes_factor(gammaFit, gammaFit6)
bayes_factor(gammaFit, gammaFit7)
bayes_factor(gammaFit, gammaFit10)

gammaFitsmooths <- conditional_smooths(gammaFit)
plot(gammaFitsmooths, 
     stype = "contour", 
     ask = FALSE,
     theme = theme(legend.position = "bottom"))

gammaFiteffects <- conditional_effects(gammaFit, 
                                       surface = TRUE)
plot(gammaFiteffects, 
     points = TRUE, 
     ask = FALSE, 
     stype = "contour")

##### Prediction ----
## Fitted
pred1 <- predict(gammaFit, type = "link")

gammaFitfinalFit <- posterior_predict(gammaFit)
gammaFitfinalResiduals <- t(StormdataTrain3$VMAX - t(gammaFitfinalFit))
gammaFitfinalResidualsMean <- colMeans(gammaFitfinalResiduals)
gammaFitfinalFitMean <- colMeans(gammaFitfinalFit)
gammaFitfinalFitMed <- apply(gammaFitfinalFit, 2, function(x){quantile(x, 0.5)})
gammaFitfinalFitLCB <- apply(gammaFitfinalFit, 2, function(x){quantile(x, 0.025)})
gammaFitfinalFitUCB <- apply(gammaFitfinalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
gammaFitfinalPreds <- posterior_predict(gammaFit, 
                                        newdata = StormdataTest3,
                                        allow_new_levels = TRUE, 
                                        re_formula = NULL)
gammaFitfinalPredsMean <- colMeans(gammaFitfinalPreds)
gammaFitfinalPredsMed <- apply(gammaFitfinalPreds, 2, function(x){quantile(x, 0.5, na.rm = TRUE)})
gammaFitfinalPredsLCB <- apply(gammaFitfinalPreds, 2, function(x){quantile(x, 0.025, na.rm = TRUE)})
gammaFitfinalPredsUCB <- apply(gammaFitfinalPreds, 2, function(x){quantile(x, 0.975, na.rm = TRUE)})

gammaFitpredMetrics <- tibble(
  Fit = paste0("gammaFit", fit),
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(gammaFitfinalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(gammaFitfinalFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < gammaFitfinalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(gammaFitfinalPredsMean - Actual_Yvec), na.rm = TRUE),
  MAD_pred = mean(abs(gammaFitfinalPredsMed - Actual_Yvec), na.rm = TRUE),
  COV_pred = mean(gammaFitfinalPredsLCB < Actual_Yvec & Actual_Yvec < gammaFitfinalPredsUCB)
)
gammaFitpredMetrics


##### Plotting ----
## Fit
ppc_dens_overlay(y = Actual_Yvec, yrep = gammaFitfinalPreds) +
  labs(title = "gammaFit Predict") +
  theme_bw()

gammaFitFitDF <- bind_cols(
  StormdataTrain3,
  LCB = gammaFitfinalFitLCB,
  Mean = gammaFitfinalFitMean,
  Med = gammaFitfinalFitMed,
  UCB = gammaFitfinalFitUCB
)

gammaFitstormsFitplot <- ggplot(data = gammaFitFitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX, color = "Observed")) +
  geom_line(aes(y = Mean, color = "PPD Mean")) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50)) +
  labs(title = "gammaFit PPD Mean vs Observed VMAX",
       subtitle = "95% Credible Interval about PPD Mean") +
  scale_color_manual(name = NULL, values = c("black","red")) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1))
  ) +
  theme_bw()
gammaFitstormsFitplot

## Residuals
gammaFitResiduals <- residuals(gammaFit, method = "posterior_predict")
ppc_error_hist(StormdataTrain8$VMAX, gammaFitfinalFit[c(1,20,100), ])
ppc_error_scatter_avg(StormdataTrain9$VMAX, gammaFitfinalFit)
ppc_error_scatter_avg_vs_x(StormdataTrain8$VMAX, 
                           gammaFitfinalFit,
                           StormdataTrain3$StormElapsedTime)
ppc_error_scatter_avg_vs_x(StormdataTrain8$VMAX, 
                           gammaFitfinalFit,
                           StormdataTrain3$HWRF)
ppc_scatter_avg_grouped(StormdataTrain9$VMAX, 
                        gammaFitfinalFit,
                        StormdataTrain3$StormID,
                        facet_args = list(scales = "free_x"))

## Prediction
gammaFitPredDF <- bind_cols(
  StormdataTest3,
  LCB = gammaFitfinalPredsLCB,
  Mean = gammaFitfinalPredsMean,
  Med = gammaFitfinalPredsMed,
  UCB = gammaFitfinalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) #|>
# filter(StormID %in% c(1812014))

gammaFitstormsPredplot <- ggplot(data = gammaFitPredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX, color = "Observed")) +
  geom_line(aes(y = Mean, color = "PPD Mean")) +
  facet_wrap(vars(StormID), ncol = 6)+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50)) +
  labs(title = "gammaFit PPD Mean vs Observed VMAX",
       subtitle = "95% Credible Interval about PPD Mean") +
  scale_color_manual(name = NULL, values = c("black","red")) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1))
  ) +
  theme_bw()
gammaFitstormsPredplot

##### PPC ----
###### Quantile 2.5 
gammaFitLCBsims <- apply(gammaFitfinalFit, 
                         MARGIN = 1,
                         function(x){
                           quantile(x, 0.025)
                         })
gammaFitLCBpvalueVec <- gammaFitLCBsims < quantile(StormdataTrain3$VMAX, 0.025)
gammaFitLCBpvalue <- sum(gammaFitLCBpvalueVec)
gammaFitLCBpvalue <- round(gammaFitLCBpvalue/4000, 3)
gammaFitLCBpvalue <- min(gammaFitLCBpvalue, 1 - gammaFitLCBpvalue)

gammaFit_ppcLCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitfinalFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", gammaFitLCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFit_ppcLCB

###### Quantile 97.5 
gammaFitUCBsims <- apply(gammaFitfinalFit, 
                         MARGIN = 1,
                         function(x){
                           quantile(x, 0.975)
                         })
gammaFitUCBpvalueVec <- gammaFitUCBsims < quantile(StormdataTrain3$VMAX, 0.975)
gammaFitUCBpvalue <- as.numeric(sum(gammaFitUCBpvalueVec))
gammaFitUCBpvalue <- round(gammaFitUCBpvalue/4000, 3)
gammaFitUCBpvalue <- min(gammaFitUCBpvalue, 1 - gammaFitUCBpvalue)

gammaFit_ppcUCB <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitfinalFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", gammaFitUCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFit_ppcUCB

###### Mean 
gammaFitMEANsims <- apply(gammaFitfinalFit, 
                          MARGIN = 1,
                          function(x){
                            mean(x)
                          })
gammaFitMEANpvalueVec <- gammaFitMEANsims < mean(StormdataTrain3$VMAX)
gammaFitMEANpvalue <- sum(gammaFitMEANpvalueVec)
gammaFitMEANpvalue <- round(gammaFitMEANpvalue/4000, 3)
gammaFitMEANpvalue <- min(gammaFitMEANpvalue, 1 - gammaFitMEANpvalue)

gammaFit_ppcMEAN <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitfinalFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", gammaFitMEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFit_ppcMEAN

###### Med 
gammaFitMEDsims <- apply(gammaFitfinalFit, 
                         MARGIN = 1,
                         function(x){
                           quantile(x, 0.5)
                         })
gammaFitMEDpvalueVec <- gammaFitMEDsims < quantile(StormdataTrain3$VMAX, 0.5)
gammaFitMEDpvalue <- sum(gammaFitMEDpvalueVec)
gammaFitMEDpvalue <- round(gammaFitMEDpvalue/4000, 3)
gammaFitMEDpvalue <- min(gammaFitMEDpvalue, 1 - gammaFitMEDpvalue)

gammaFit_ppcMED <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitfinalFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", gammaFitMEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFit_ppcMED

###### SD 
gammaFitSDsims <- apply(gammaFitfinalFit, 
                        MARGIN = 1,
                        function(x){
                          sd(x)
                        })
gammaFitSDpvalueVec <- gammaFitSDsims < sd(StormdataTrain3$VMAX)
gammaFitSDpvalue <- sum(gammaFitSDpvalueVec)
gammaFitSDpvalue <- round(gammaFitSDpvalue/4000, 3)
gammaFitSDpvalue <- min(gammaFitSDpvalue, 1 - gammaFitSDpvalue)

gammaFit_ppcSD <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitfinalFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", gammaFitSDpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFit_ppcSD

###### Range 
gammaFitRANGEsims <- apply(gammaFitfinalFit, 
                           MARGIN = 1,
                           function(x){
                             max(x)-min(x)
                           })
gammaFitRANGEpvalueVec <- gammaFitRANGEsims < (max(StormdataTrain3$VMAX)-min(StormdataTrain3$VMAX))
gammaFitRANGEpvalue <- sum(gammaFitRANGEpvalueVec)
gammaFitRANGEpvalue <- round(gammaFitRANGEpvalue/4000, 3)
gammaFitRANGEpvalue <- min(gammaFitRANGEpvalue, 1 - gammaFitRANGEpvalue)

gammaFit_ppcRANGE <- 
  ppc_stat(StormdataTrain3$VMAX,
           gammaFitfinalFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", gammaFitRANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#gammaFit_ppcRANGE

##### Combined Plot ----
gammaFit_ppcComb <- 
  gammaFitppcFit /
  (gammaFit_ppcLCB | gammaFit_ppcMED | gammaFit_ppcUCB) /
  (gammaFit_ppcRANGE | gammaFit_ppcMEAN | gammaFit_ppcSD)
# plot_layout(ncol = 3, 
#             widths = c(3,1,1,1,1,1,1), 
#             heights = c(1,1,1,1,1,1,1),
#             byrow = FALSE)
gammaFit_ppcComb

##### Bayes p-values ----
gammaFitpvalues <- tibble(
  Fit = paste0("gammaFit", fit),
  LCB = gammaFitLCBpvalue,
  Median = gammaFitMEDpvalue,
  UCB = gammaFitUCBpvalue,
  Range = gammaFitRANGEpvalue,
  Mean = gammaFitMEANpvalue,
  SD = gammaFitSDpvalue
)
gammaFitpvalues

##### CV ----
#kfoldID <- kfold_split_grouped(K = 5, StormdataTrain8$StormID)
gammaFitkfoldgroup <- kfold(gammaFit,
                            folds = kfoldID,
                            chains = 1,
                            save_fits = TRUE)
gammaFitkfoldrand <- kfold(gammaFit,
                           K = 5,
                           chains = 1,
                           save_fits = TRUE)
gammaFitkfoldPreds <- kfold_predict(gammaFitkfoldgroup)
#gammaFitkfoldPreds <- kfold_predict(gammaFitkfold)
gammaFitkfoldPredsDat <- gammaFitkfoldPreds$yrep
gammaFitkfoldPredsMean <- colMeans(gammaFitkfoldPredsDat)
gammaFitkfoldPredsMed <- apply(gammaFitkfoldPredsDat, 2, function(x){quantile(x, 0.5)})
gammaFitkfoldPredsLCB <- apply(gammaFitkfoldPredsDat, 2, function(x){quantile(x, 0.025)})
gammaFitkfoldPredsUCB <- apply(gammaFitkfoldPredsDat, 2, function(x){quantile(x, 0.975)})

gammaFitkfoldMetrics <- tibble(
  Fit = paste0("gammaFit", fit),
  MAE_kfold = mean(abs(gammaFitkfoldPredsMean - gammaFitkfoldPreds$y)),
  MAD_kfold = mean(abs(gammaFitkfoldPredsMed - gammaFitkfoldPreds$y)),
  COV_kfold = mean(gammaFitkfoldPredsLCB < gammaFitkfoldPreds$y & gammaFitkfoldPreds$y < gammaFitkfoldPredsUCB)
)
gammaFitkfoldMetrics


# Compare Models ----
## Bayes pvalue ----
pvalueComp_gamma <- gammaFitNULLidentitypvalues
pvalueComp_gammatemp <- gammaFitpvalues
pvalueComp_gamma <- bind_rows(
  pvalueComp_gamma,
  pvalueComp_gammatemp
)
pvalueComp_gamma
save(pvalueComp_gamma, 
     file = "~/Desktop/Temp Hurricane Model Data/pvalueComp_gamma.RData")

## LOO ----
# gammaFitloo1 <- gammaFitloo
# attributes(gammaFitloo1)$model_name <- "gammaFit1"
gammaFitloo2 <- gammaFitloo
attributes(gammaFitloo2)$model_name <- "gammaFit2"
# gammaFitloo3 <- gammaFitloo
# attributes(gammaFitloo3)$model_name <- "gammaFit3"
#gammaFitloo5 <- loo(gammaFit5)
#attributes(gammaFitloo5)$model_name <- "gammaFit5"
#gammaFitloo6 <- loo(gammaFit)
#attributes(gammaFitloo6)$model_name <- "gammaFit6"
# gammaFitloo8 <- loo(gammaFit)
# attributes(gammaFitloo8)$model_name <- "gammaFit8"
# gammaFitloo9 <- loo(gammaFit)
# attributes(gammaFitloo9)$model_name <- "gammaFit9"
# gammaFitloo11 <- loo(gammaFit)
# attributes(gammaFitloo11)$model_name <- "gammaFit11"
# gammaFitloo12 <- loo(gammaFit)
# attributes(gammaFitloo12)$model_name <- "gammaFit12"
looComp_gamma <- loo_compare(gammaFitNULLidentityloo,
                             gammaFitloo1,
                             gammaFitloo2
                             # gammaFitloo3,
                             # gammaFitloo5,
                             # gammaFitloo6,
                             # gammaFitloo8,
                             # gammaFitloo9,
                             # gammaFitloo11,
                             # gammaFitloo12
)

looComp_gamma
save(looComp_gamma, 
     file = "~/Desktop/Temp Hurricane Model Data/looComp_gamma.RData")

## CV ----
#cvComp_gamma <- gammaFitNULL2kfoldMetrics
cvComp_gammatemp <- gammaFitkfoldMetrics
cvComp_gamma <- bind_rows(
  cvComp_gamma,
  cvComp_gammatemp
)
cvComp_gamma <- cvComp_gamma |> arrange(MAE_kfold)
cvComp_gamma
save(cvComp_gamma, 
     file = "~/Desktop/Temp Hurricane Model Data/cvComp_gamma.RData")


## Preds ----
predComp_gamma <- gammaFitNULLidentitypredMetrics
predComp_gammatemp <- gammaFitpredMetrics
predComp_gamma <- bind_rows(
  predComp_gamma,
  predComp_gammatemp
)
predComp_gamma <- predComp_gamma |> arrange(MAE_pred)
predComp_gamma
save(predComp_gamma, 
     file = "~/Desktop/Temp Hurricane Model Data/predComp_gamma.RData")


### Bayes R2
# bayesR2_gamma <- bayes_R2(gammaFitNULL2) |> 
#   bind_cols(Fit = "NULL")
bayesR2_gammatemp <- bayes_R2(gammaFit) |> 
  bind_cols(Fit = paste0("gammaFit", fit))
bayesR2_gamma <- bind_rows(
  bayesR2_gamma,
  bayesR2_gammatemp
)
bayesR2_gamma <- bayesR2_gamma |> arrange(desc(Estimate))
bayesR2_gamma

gammaFit6 <- gammaFit
save(gammaFit8,
     file = "~/Desktop/Temp Hurricane Model Data/gammaFit8.RData")


gammaFitNULL2files <- ls()[str_detect(ls(), pattern = "gammaFitNULL2")]
gammaFitNULL2filesRM <- gammaFitNULL2files[(gammaFitNULL2files %in% c(
  "gammaFitNULL2",
  "gammaFitNULL2_ppcComb",
  "gammaFitNULL2kfold",
  "gammaFitNULL2kfoldMetrics",
  "gammaFitNULL2loo",
  "gammaFitNULL2predMetrics",
  "gammaFitNULL2pvalues"))]

rm(list = ls()[!(ls() %in% gammaFitNULL2filesRM)])
rm(gammaFitNULLfilesRM)

