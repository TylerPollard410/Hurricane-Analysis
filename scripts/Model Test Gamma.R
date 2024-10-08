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
StormdataTrain1 <- Stormdata |> filter(complete.cases(VMAX))

# Remove not varying 
StormdataTrain2 <- StormdataTrain1 |>
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
    across(where(is.numeric) & !c(VMAX, 
                                  HWRF, 
                                  StormElapsedTime,StormElapsedTime2,
                                  LAT, LON),
           function(x){scale(x)})
  )
str(StormdataTrain9)

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
## LOGNORMAL ----
gammaFitNULL2 <- brm(
  bf(
    VMAX ~ offset(log(HWRF)), family = lognormal()
  ),
  data = StormdataTrain8, 
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000
)
print(gammaFitNULL2, digits = 4)

## gammaFit1 <- gammaFit
## gammaFit2 <- gammaFit

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

gammaFitcode <- stancode(
  bf(VMAX ~ offset(log(HWRF)) +
       t2(LON, LAT, StormElapsedTime, d = c(2,1), bs = c("tp", "cr")) +
       # Year +
       # Month +
       # basin + 
       # LON +
       # LAT +
       # s(Day) +
       # t2(LON, LAT) +
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
gammaFitcode

### Model ----
fit <- 12
gammaFit <- brm(
  bf(VMAX ~ 
       s(StormElapsedTime, bs = "cr") +
       s(LON, bs = "cr") +
       s(LAT, bs = "cr") +
       basin + 
       #s(LON, basin, bs = "fs") +
       #s(LAT, basin, bs = "fs") +
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
       log(HWFI) +
       VMAX_OP_T0 +
       log(HWRF) +
       (1|StormID)
  ),
  data = StormdataTest3,
  family = lognormal(),
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 2000,
  seed = 52,
  warmup = 1000,
  control = list(adapt_delta = 0.95)
)

gammaFit12 <- gammaFit
save(gammaFit8,
     file = "~/Desktop/Temp Hurricane Model Data/gammaFit8.RData")
prior_summary(gammaFit)
posterior_summary(gammaFit)
gammaFit

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
gammaFitfinalResiduals <- t(StormdataTrain9$VMAX - t(gammaFitfinalFit))
gammaFitfinalResidualsMean <- colMeans(gammaFitfinalResiduals)
gammaFitfinalFitMean <- colMeans(gammaFitfinalFit)
gammaFitfinalFitMed <- apply(gammaFitfinalFit, 2, function(x){quantile(x, 0.5)})
gammaFitfinalFitLCB <- apply(gammaFitfinalFit, 2, function(x){quantile(x, 0.025)})
gammaFitfinalFitUCB <- apply(gammaFitfinalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
gammaFitfinalPreds <- posterior_predict(gammaFit, 
                                            newdata = StormdataTest7scale,
                                            allow_new_levels = TRUE, 
                                            re_formula = NULL)
gammaFitfinalPredsMean <- colMeans(gammaFitfinalPreds)
gammaFitfinalPredsMed <- apply(gammaFitfinalPreds, 2, function(x){quantile(x, 0.5)})
gammaFitfinalPredsLCB <- apply(gammaFitfinalPreds, 2, function(x){quantile(x, 0.025)})
gammaFitfinalPredsUCB <- apply(gammaFitfinalPreds, 2, function(x){quantile(x, 0.975)})

gammaFitpredMetrics <- tibble(
  Fit = paste0("gammaFit", fit),
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(gammaFitfinalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(gammaFitfinalFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < gammaFitfinalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(gammaFitfinalPredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(gammaFitfinalPredsMed - Actual_Yvec)),
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

gammaFitstormsFitplot <- 
  ggplot(data = gammaFitFitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))
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
#pvalueComp_gamma <- pvalueCompNULL |> filter(Fit == "gammaFitNULL2")
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
# gammaFitloo2 <- gammaFitloo
# attributes(gammaFitloo2)$model_name <- "gammaFit2"
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
gammaFitloo12 <- loo(gammaFit)
attributes(gammaFitloo12)$model_name <- "gammaFit12"
looComp_gamma <- loo_compare(gammaFitNULL2loo,
                                 gammaFitloo1,
                                 gammaFitloo2,
                                 gammaFitloo3,
                                 gammaFitloo5,
                                 gammaFitloo6,
                                 gammaFitloo8,
                                 gammaFitloo9,
                                 gammaFitloo11,
                                 gammaFitloo12)

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
#predComp_gamma <- gammaFitNULL2predMetrics
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



