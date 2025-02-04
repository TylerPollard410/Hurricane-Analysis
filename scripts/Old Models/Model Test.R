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
    across(where(is.numeric) & !c(VMAX, HWRF, StormElapsedTime,StormElapsedTime2),
           function(x){scale(x)})
  ) |>
  mutate(
    propVMAX = VMAX/HWRF
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
logNormalFitNULL2 <- brm(
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
print(logNormalFitNULL2, digits = 4)

## logNormalFit1 <- logNormalFit
## logNormalFit2 <- logNormalFit

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

logNormalFitcode <- stancode(
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
logNormalFitcode

### Model ----
fit <- 12
logNormalFit <- brm(
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

logNormalFit12 <- logNormalFit
save(logNormalFit8,
     file = "~/Desktop/Temp Hurricane Model Data/logNormalFit8.RData")
prior_summary(logNormalFit)
posterior_summary(logNormalFit)
logNormalFit

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
logNormalFitfinalResiduals <- t(StormdataTrain9$VMAX - t(logNormalFitfinalFit))
logNormalFitfinalResidualsMean <- colMeans(logNormalFitfinalResiduals)
logNormalFitfinalFitMean <- colMeans(logNormalFitfinalFit)
logNormalFitfinalFitMed <- apply(logNormalFitfinalFit, 2, function(x){quantile(x, 0.5)})
logNormalFitfinalFitLCB <- apply(logNormalFitfinalFit, 2, function(x){quantile(x, 0.025)})
logNormalFitfinalFitUCB <- apply(logNormalFitfinalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
logNormalFitfinalPreds <- posterior_predict(logNormalFit, 
                                            newdata = StormdataTest7scale,
                                            allow_new_levels = TRUE, 
                                            re_formula = NULL)
logNormalFitfinalPredsMean <- colMeans(logNormalFitfinalPreds)
logNormalFitfinalPredsMed <- apply(logNormalFitfinalPreds, 2, function(x){quantile(x, 0.5)})
logNormalFitfinalPredsLCB <- apply(logNormalFitfinalPreds, 2, function(x){quantile(x, 0.025)})
logNormalFitfinalPredsUCB <- apply(logNormalFitfinalPreds, 2, function(x){quantile(x, 0.975)})

logNormalFitpredMetrics <- tibble(
  Fit = paste0("logNormalFit", fit),
  MAE_HWRF_fit = mean(abs(StormdataTrain3$HWRF - StormdataTrain3$VMAX)),
  MAE_fit = mean(abs(logNormalFitfinalFitMean - StormdataTrain3$VMAX)),
  COV_fit = mean(logNormalFitfinalFitLCB < StormdataTrain3$VMAX & StormdataTrain3$VMAX < logNormalFitfinalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest2$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(logNormalFitfinalPredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(logNormalFitfinalPredsMed - Actual_Yvec)),
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

logNormalFitstormsFitplot <- 
  ggplot(data = logNormalFitFitDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX), color = "red") +
  geom_line(aes(y = Mean)) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50))
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
#pvalueComp_logNormal <- pvalueCompNULL |> filter(Fit == "logNormalFitNULL2")
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
# logNormalFitloo2 <- logNormalFitloo
# attributes(logNormalFitloo2)$model_name <- "logNormalFit2"
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
logNormalFitloo12 <- loo(logNormalFit)
attributes(logNormalFitloo12)$model_name <- "logNormalFit12"
looComp_logNormal <- loo_compare(logNormalFitNULL2loo,
                                 logNormalFitloo1,
                                 logNormalFitloo2,
                                 logNormalFitloo3,
                                 logNormalFitloo5,
                                 logNormalFitloo6,
                                 logNormalFitloo8,
                                 logNormalFitloo9,
                                 logNormalFitloo11,
                                 logNormalFitloo12)

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
#predComp_logNormal <- logNormalFitNULL2predMetrics
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



