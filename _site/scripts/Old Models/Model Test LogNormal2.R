# Load Libraries ----
library(knitr)
library(data.table)
library(MASS)
library(bestNormalize)
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

### Summary of data ----
dataSum <- Stormdata |>
  group_by(DataType) |>
  reframe(
    across(where(is.numeric), function(x){as.numeric(summary(x))})
  )

## Training ----
StormdataTrain <- Stormdata |> 
  filter(DataType == "Train") |>
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
    everything(),
    -DataType
  )

#### Train 1 ----
StormdataTrain1 <- StormdataTrain |>
  select(
    "StormID",
    #Date,
    Year,
    Month,
    Day,
    StormElapsedTime,
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
                                  StormElapsedTime,
                                  LAT, LON),
           function(x){scale(x)})
  )
str(StormdataTrain1)

#### Train 2 ----
StormdataTrain2 <- StormdataTrain |>
  select(
    "StormID",
    #Date,
    Year,
    Month,
    Day,
    StormElapsedTime,
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
    across(where(is.numeric) & !c(VMAX),
           function(x){scale(x)})
  )
str(StormdataTrain2)

#### Train 3 ----
StormdataTrain3 <- StormdataTrain |>
  select(
    "StormID",
    #Date,
    Year,
    Month,
    Day,
    StormElapsedTime,
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
    across(where(is.numeric) & !c(VMAX, HWRF),
           function(x){scale(x)})
  ) |>
  mutate(
    propVMAX = VMAX/HWRF
  )
str(StormdataTrain3)

#### Train 4 ----
StormdataTrain4 <- StormdataTrain |>
  select(
    "StormID",
    #Date,
    Year,
    Month,
    Day,
    StormElapsedTime,
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
    LON2 = LON,
    LAT2 = LAT,
    StormElapsedTime2 = StormElapsedTime,
    logMINSLP = log(MINSLP),
    logSHR_MAG = log(SHR_MAG),
    logSTM_SPD = log(STM_SPD),
    logSST = log(SST + 1 - min(SST)),
    logRHLO = log(RHLO),
    logCAPE1 = log(CAPE1),
    logCAPE3 = log(CAPE3),
    logSHTFL2 = log(SHTFL2 + 1 - min(SHTFL2)),
    logTCOND7002 = log(TCOND7002 + 1 - min(TCOND7002)),
    logINST2 = log(INST2 + 1 - min(INST2)),
    logCP1 = log(CP1 + 1 - min(CP1)),
    logTCONDSYM2 = log(TCONDSYM2),
    logCOUPLSYM3 = log(COUPLSYM3 + 1 - min(COUPLSYM3)),
    logHWFI = log(HWFI),
    logVMAX_OP_T0 = log(VMAX_OP_T0),
    logHWRF = log(HWRF)
  ) |>
  mutate(
    across(where(is.numeric) & !c(VMAX, 
                                  # HWRF, 
                                  logHWRF,
                                  StormElapsedTime2,
                                  LAT2, LON2
    ),
    function(x){scale(x)}),
    # logMINSLPscale = scale(logMINSLP, center = TRUE, scale = TRUE),
    # logSHR_MAGscale = scale(logSHR_MAG, center = TRUE, scale = TRUE),
    # logSTM_SPDscale = scale(logSTM_SPD, center = TRUE, scale = TRUE),
    # logSSTscale = scale(logSST, center = TRUE, scale = TRUE),
    # logRHLOscale = scale(logRHLO, center = TRUE, scale = TRUE),
    # logCAPE1scale = scale(logCAPE1, center = TRUE, scale = TRUE),
    # logCAPE3scale = scale(logCAPE3, center = TRUE, scale = TRUE),
    # logSHTFL2scale = scale(logSHTFL2, center = TRUE, scale = TRUE),
    # logTCOND7002scale = scale(logTCOND7002, center = TRUE, scale = TRUE),
    # logINST2scale = scale(logINST2, center = TRUE, scale = TRUE),
    # logCP1scale = scale(logCP1, center = TRUE, scale = TRUE),
    # logCOUPLSYM3scale = scale(logCOUPLSYM3, center = TRUE, scale = TRUE),
    # logHWFIscale = scale(logHWFI, center = TRUE, scale = TRUE),
    # logVMAX_OP_T0scale = scale(logVMAX_OP_T0, center = TRUE, scale = TRUE),
    logHWRFcenter = scale(logHWRF, center = TRUE, scale = FALSE),
    logHWRFscale = scale(logHWRF, center = TRUE, scale = TRUE),
    logHWRF = scale(logHWRF)
  )
str(StormdataTrain4)

#### Train 4 Unscale ----
StormdataTrain4unscale <- StormdataTrain |>
  select(
    "StormID",
    #Date,
    Year,
    Month,
    Day,
    StormElapsedTime,
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
    logMINSLP = log(MINSLP),
    logSHR_MAG = log(SHR_MAG),
    logSTM_SPD = log(STM_SPD),
    logSST = log(SST + 1 - min(SST)),
    logRHLO = log(RHLO),
    logCAPE1 = log(CAPE1),
    logCAPE3 = log(CAPE3),
    logSHTFL2 = log(SHTFL2 + 1 - min(SHTFL2)),
    logTCOND7002 = log(TCOND7002 + 1 - min(TCOND7002)),
    logINST2 = log(INST2 + 1 - min(INST2)),
    logCP1 = log(CP1 + 1 - min(CP1)),
    logTCONDSYM2 = log(TCONDSYM2),
    logCOUPLSYM3 = log(COUPLSYM3 + 1 - min(COUPLSYM3)),
    logHWFI = log(HWFI),
    logVMAX_OP_T0 = log(VMAX_OP_T0),
    logHWRF = log(HWRF)
  )
str(StormdataTrain4unscale)

#### Train 5 ----
StormdataTrain5 <- StormdataTrain |>
  select(
    "StormID",
    #Date,
    Year,
    Month,
    Day,
    StormElapsedTime,
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
    LON2 = LON,
    LAT2 = LAT,
    StormElapsedTime2 = StormElapsedTime,
    arcsinhMINSLP = log(MINSLP + sqrt(MINSLP^2 + 1)),
    arcsinhSHR_MAG = log(SHR_MAG + sqrt(SHR_MAG^2 + 1)),
    arcsinhSTM_SPD = log(STM_SPD + sqrt(STM_SPD^2 + 1)),
    arcsinhSST = log(SST + sqrt(SST^2 + 1)),
    arcsinhRHLO = log(RHLO + sqrt(RHLO^2 + 1)),
    arcsinhCAPE1 = log(CAPE1 + sqrt(CAPE1^2 + 1)),
    arcsinhCAPE3 = log(CAPE3 + sqrt(CAPE3^2 + 1)),
    arcsinhSHTFL2 = log(SHTFL2 + sqrt(SHTFL2^2 + 1)),
    arcsinhTCOND7002 = log(TCOND7002 + sqrt(TCOND7002^2 + 1)),
    arcsinhINST2 = log(INST2 + sqrt(INST2^2 + 1)),
    arcsinhCP1 = log(CP1 + sqrt(CP1^2 + 1)),
    arcsinhTCONDSYM2 = log(TCONDSYM2 + sqrt(TCONDSYM2^2 + 1)),
    arcsinhCOUPLSYM3 = log(COUPLSYM3 + sqrt(COUPLSYM3^2 + 1)),
    arcsinhHWFI = log(HWFI + sqrt(HWFI^2 + 1)),
    arcsinhVMAX_OP_T0 = log(VMAX_OP_T0 + sqrt(VMAX_OP_T0^2 + 1)),
    arcsinhHWRF = log(HWRF + sqrt(HWRF^2 + 1))
  ) |>
  mutate(
    across(where(is.numeric) & !c(VMAX, StormElapsedTime2, LAT2, LON2
    ),
    function(x){scale(x)})
  )
str(StormdataTrain5)

#### Train 5 Unscale ----
StormdataTrain5unscale <- StormdataTrain |>
  select(
    "StormID",
    #Date,
    Year,
    Month,
    Day,
    StormElapsedTime,
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
    arcsinhMINSLP = log(MINSLP + sqrt(MINSLP^2 + 1)),
    arcsinhSHR_MAG = log(SHR_MAG + sqrt(SHR_MAG^2 + 1)),
    arcsinhSTM_SPD = log(STM_SPD + sqrt(STM_SPD^2 + 1)),
    arcsinhSST = log(SST + sqrt(SST^2 + 1)),
    arcsinhRHLO = log(RHLO + sqrt(RHLO^2 + 1)),
    arcsinhCAPE1 = log(CAPE1 + sqrt(CAPE1^2 + 1)),
    arcsinhCAPE3 = log(CAPE3 + sqrt(CAPE3^2 + 1)),
    arcsinhSHTFL2 = log(SHTFL2 + sqrt(SHTFL2^2 + 1)),
    arcsinhTCOND7002 = log(TCOND7002 + sqrt(TCOND7002^2 + 1)),
    arcsinhINST2 = log(INST2 + sqrt(INST2^2 + 1)),
    arcsinhCP1 = log(CP1 + sqrt(CP1^2 + 1)),
    arcsinhTCONDSYM2 = log(TCONDSYM2 + sqrt(TCONDSYM2^2 + 1)),
    arcsinhCOUPLSYM3 = log(COUPLSYM3 + sqrt(COUPLSYM3^2 + 1)),
    arcsinhHWFI = log(HWFI + sqrt(HWFI^2 + 1)),
    arcsinhVMAX_OP_T0 = log(VMAX_OP_T0 + sqrt(VMAX_OP_T0^2 + 1)),
    arcsinhHWRF = log(HWRF + sqrt(HWRF^2 + 1))
  )
str(StormdataTrain5unscale)

#### Train 6 ----
StormdataTrain6 <- StormdataTrain |>
  select(
    "StormID",
    #Date,
    Year,
    Month,
    Day,
    StormElapsedTime,
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
    LON2 = LON,
    LAT2 = LAT,
    StormElapsedTime2 = StormElapsedTime,
    arcsinhLON = log(LON + sqrt(LON^2 + 1)),
    arcsinhLAT = log(LAT + sqrt(LAT^2 + 1)),
    arcsinhStormElapsedTime = log(StormElapsedTime + sqrt(StormElapsedTime^2 + 1)),
    arcsinhMINSLP = log(MINSLP + sqrt(MINSLP^2 + 1)),
    arcsinhSHR_MAG = log(SHR_MAG + sqrt(SHR_MAG^2 + 1)),
    arcsinhSTM_SPD = log(STM_SPD + sqrt(STM_SPD^2 + 1)),
    arcsinhSST = log(SST + sqrt(SST^2 + 1)),
    arcsinhRHLO = log(RHLO + sqrt(RHLO^2 + 1)),
    arcsinhCAPE1 = log(CAPE1 + sqrt(CAPE1^2 + 1)),
    arcsinhCAPE3 = log(CAPE3 + sqrt(CAPE3^2 + 1)),
    arcsinhSHTFL2 = log(SHTFL2 + sqrt(SHTFL2^2 + 1)),
    arcsinhTCOND7002 = log(TCOND7002 + sqrt(TCOND7002^2 + 1)),
    arcsinhINST2 = log(INST2 + sqrt(INST2^2 + 1)),
    arcsinhCP1 = log(CP1 + sqrt(CP1^2 + 1)),
    arcsinhTCONDSYM2 = log(TCONDSYM2 + sqrt(TCONDSYM2^2 + 1)),
    arcsinhCOUPLSYM3 = log(COUPLSYM3 + sqrt(COUPLSYM3^2 + 1)),
    arcsinhHWFI = log(HWFI + sqrt(HWFI^2 + 1)),
    arcsinhVMAX_OP_T0 = log(VMAX_OP_T0 + sqrt(VMAX_OP_T0^2 + 1)),
    arcsinhHWRF = log(HWRF + sqrt(HWRF^2 + 1))
  ) |>
  mutate(
    across(where(is.numeric) & !c(VMAX, StormElapsedTime2, LAT2, LON2
    ),
    function(x){scale(x)})
  )
str(StormdataTrain6)

#### Train 6 Unscale ----
StormdataTrain5unscale <- StormdataTrain |>
  select(
    "StormID",
    #Date,
    Year,
    Month,
    Day,
    StormElapsedTime,
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
    arcsinhMINSLP = log(MINSLP + sqrt(MINSLP^2 + 1)),
    arcsinhSHR_MAG = log(SHR_MAG + sqrt(SHR_MAG^2 + 1)),
    arcsinhSTM_SPD = log(STM_SPD + sqrt(STM_SPD^2 + 1)),
    arcsinhSST = log(SST + sqrt(SST^2 + 1)),
    arcsinhRHLO = log(RHLO + sqrt(RHLO^2 + 1)),
    arcsinhCAPE1 = log(CAPE1 + sqrt(CAPE1^2 + 1)),
    arcsinhCAPE3 = log(CAPE3 + sqrt(CAPE3^2 + 1)),
    arcsinhSHTFL2 = log(SHTFL2 + sqrt(SHTFL2^2 + 1)),
    arcsinhTCOND7002 = log(TCOND7002 + sqrt(TCOND7002^2 + 1)),
    arcsinhINST2 = log(INST2 + sqrt(INST2^2 + 1)),
    arcsinhCP1 = log(CP1 + sqrt(CP1^2 + 1)),
    arcsinhTCONDSYM2 = log(TCONDSYM2 + sqrt(TCONDSYM2^2 + 1)),
    arcsinhCOUPLSYM3 = log(COUPLSYM3 + sqrt(COUPLSYM3^2 + 1)),
    arcsinhHWFI = log(HWFI + sqrt(HWFI^2 + 1)),
    arcsinhVMAX_OP_T0 = log(VMAX_OP_T0 + sqrt(VMAX_OP_T0^2 + 1)),
    arcsinhHWRF = log(HWRF + sqrt(HWRF^2 + 1))
  )
str(StormdataTrain5unscale)

## Test ----
StormdataTest <- Stormdata |> 
  filter(DataType == "Test") |>
  mutate(
    StormID = droplevels(StormID),
    Year = factor(Year, ordered = TRUE)
    #Month = dataTestMonths,
    #Day = dataTestYearDay
  ) |>
  select(
    StormID,
    Date,
    Year,
    Month,
    Day,
    StormElapsedTime,
    StormElapsedTime2,
    everything(),
    -DataType
  )

#### Test 1 ----
StormdataTest1 <- StormdataTest |>
  select(
    "StormID",
    #Date,
    Year,
    Month,
    Day,
    StormElapsedTime,
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
                                  StormElapsedTime,
                                  LAT, LON),
           function(x){scale(x,
                             center = attr(StormdataTrain1 |> pull(x), "scaled:center"),
                             scale = attr(StormdataTrain1 |> pull(x), "scaled:scale"))
           })
  )
str(StormdataTest1)

#### Test 2 ----
StormdataTest2 <- StormdataTest |>
  select(
    "StormID",
    #Date,
    Year,
    Month,
    Day,
    StormElapsedTime,
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
    across(where(is.numeric) & !c(VMAX),
           function(x){scale(x,
                             center = attr(StormdataTrain2 |> pull(x), "scaled:center"),
                             scale = attr(StormdataTrain2 |> pull(x), "scaled:scale"))
           })
  )
str(StormdataTest2)

#### Test 3 ----
StormdataTest3 <- StormdataTest |>
  select(
    "StormID",
    #Date,
    Year,
    Month,
    Day,
    StormElapsedTime,
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
    across(where(is.numeric) & !c(VMAX, HWRF),
           function(x){scale(x,
                             center = attr(StormdataTrain3 |> pull(x), "scaled:center"),
                             scale = attr(StormdataTrain3 |> pull(x), "scaled:scale"))
           })
  ) |>
  mutate(
    propVMAX = VMAX/HWRF
  )
str(StormdataTest3)

#### Test 4 ----
StormdataTest4 <- StormdataTest |>
  select(
    "StormID",
    #Date,
    Year,
    Month,
    Day,
    StormElapsedTime,
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
    LON2 = LON,
    LAT2 = LAT,
    StormElapsedTime2 = StormElapsedTime,
    logMINSLP = log(MINSLP),
    logSHR_MAG = log(SHR_MAG),
    logSTM_SPD = log(STM_SPD),
    #logSST = log(SST + 1 - min(StormdataTrain4$SST)),
    logRHLO = log(RHLO),
    logCAPE1 = log(CAPE1),
    logCAPE3 = log(CAPE3),
    logSHTFL2 = log(SHTFL2 + 1 - min(StormdataTrain4$SHTFL2)),
    logTCOND7002 = log(TCOND7002 + 1 - min(StormdataTrain4$TCOND7002)),
    #logINST2 = log(INST2 + 1 - min(StormdataTrain4$INST2)),
    logCP1 = log(CP1 + 1 - min(StormdataTrain4$CP1)),
    logTCONDSYM2 = log(TCONDSYM2),
    #logCOUPLSYM3 = log(COUPLSYM3 + 1 - min(StormdataTrain4$COUPLSYM3)),
    logHWFI = log(HWFI),
    logVMAX_OP_T0 = log(VMAX_OP_T0),
    logHWRF = log(HWRF)
  ) |>
  mutate(
    across(where(is.numeric) & !c(VMAX, 
                                  # HWRF, 
                                  #logHWRF
                                  StormElapsedTime2,
                                  LAT2, LON2
    ),
    function(x){scale(x,
                      center = attr(StormdataTrain4 |> pull(x), "scaled:center"),
                      scale = attr(StormdataTrain4 |> pull(x), "scaled:scale"))
    })
    #   logMINSLPscale = scale(logMINSLP,
    #                          center = attr(StormdataTrain4 |> pull(logMINSLPscale), "scaled:center"),
    #                          scale = attr(StormdataTrain4 |> pull(logMINSLPscale), "scaled:scale")),
    #   logSHR_MAGscale = scale(logSHR_MAG,
    #                           center = attr(StormdataTrain4 |> pull(logSHR_MAGscale), "scaled:center"),
    #                           scale = attr(StormdataTrain4 |> pull(logSHR_MAGscale), "scaled:scale")),
    #   logSTM_SPDscale = scale(logSTM_SPD,
    #                           center = attr(StormdataTrain4 |> pull(logSTM_SPDscale), "scaled:center"),
    #                           scale = attr(StormdataTrain4 |> pull(logSTM_SPDscale), "scaled:scale")),
    #   logSSTscale = scale(logSST,
    #                       center = attr(StormdataTrain4 |> pull(logSSTscale), "scaled:center"),
    #                       scale = attr(StormdataTrain4 |> pull(logSSTscale), "scaled:scale")),
    #   logRHLOscale = scale(logRHLO,
    #                        center = attr(StormdataTrain4 |> pull(logRHLOscale), "scaled:center"),
    #                        scale = attr(StormdataTrain4 |> pull(logRHLOscale), "scaled:scale")),
    #   logCAPE1scale = scale(logCAPE1,
    #                         center = attr(StormdataTrain4 |> pull(logCAPE1scale), "scaled:center"),
    #                         scale = attr(StormdataTrain4 |> pull(logCAPE1scale), "scaled:scale")),
    #   logCAPE3scale = scale(logCAPE3,
    #                         center = attr(StormdataTrain4 |> pull(logCAPE3scale), "scaled:center"),
    #                         scale = attr(StormdataTrain4 |> pull(logCAPE3scale), "scaled:scale")),
    #   logTCOND7002scale = scale(logTCOND7002,
    #                             center = attr(StormdataTrain4 |> pull(logTCOND7002scale), "scaled:center"),
    #                             scale = attr(StormdataTrain4 |> pull(logTCOND7002scale), "scaled:scale")),
    #   logINST2scale = scale(logINST2,
    #                         center = attr(StormdataTrain4 |> pull(logINST2scale), "scaled:center"),
    #                         scale = attr(StormdataTrain4 |> pull(logINST2scale), "scaled:scale")),
    #   logCP1scale = scale(logCP1,
    #                       center = attr(StormdataTrain4 |> pull(logCP1scale), "scaled:center"),
    #                       scale = attr(StormdataTrain4 |> pull(logCP1scale), "scaled:scale")),
    #   logCOUPLSYM3scale = scale(logCOUPLSYM3,
    #                             center = attr(StormdataTrain4 |> pull(logCOUPLSYM3scale), "scaled:center"),
    #                             scale = attr(StormdataTrain4 |> pull(logCOUPLSYM3scale), "scaled:scale")),
    #   logHWFIscale = scale(logHWFI,
    #                        center = attr(StormdataTrain4 |> pull(logHWFIscale), "scaled:center"),
    #                        scale = attr(StormdataTrain4 |> pull(logHWFIscale), "scaled:scale")),
    #   logVMAX_OP_T0scale = scale(logVMAX_OP_T0,
    #                              center = attr(StormdataTrain4 |> pull(logVMAX_OP_T0scale), "scaled:center"),
    #                              scale = attr(StormdataTrain4 |> pull(logVMAX_OP_T0scale), "scaled:scale")),
    #   logHWRFcenter = scale(logHWRF, 
    #                         center = attr(StormdataTrain4 |> pull(logHWRFcenter), "scaled:center"),
    #                         scale = FALSE),
    #   logHWRFscale = scale(logHWRF, 
    #                        center = attr(StormdataTrain4 |> pull(logHWRFscale), "scaled:center"),
    #                        scale = attr(StormdataTrain4 |> pull(logHWRFscale), "scaled:scale"))
  )
str(StormdataTest4)

#### Test 4 ----
StormdataTest4unscale <- StormdataTest |>
  select(
    "StormID",
    #Date,
    Year,
    Month,
    Day,
    StormElapsedTime,
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
    logMINSLP = log(MINSLP),
    logSHR_MAG = log(SHR_MAG),
    logSTM_SPD = log(STM_SPD),
    #logSST = log(SST + 1 - min(StormdataTrain4$SST)),
    logRHLO = log(RHLO),
    logCAPE1 = log(CAPE1),
    logCAPE3 = log(CAPE3),
    logSHTFL2 = log(SHTFL2 + 1 - min(StormdataTrain4$SHTFL2)),
    logTCOND7002 = log(TCOND7002 + 1 - min(StormdataTrain4$TCOND7002)),
    #logINST2 = log(INST2 + 1 - min(StormdataTrain4$INST2)),
    logCP1 = log(CP1 + 1 - min(StormdataTrain4$CP1)),
    logTCONDSYM2 = log(TCONDSYM2),
    #logCOUPLSYM3 = log(COUPLSYM3 + 1 - min(StormdataTrain4$COUPLSYM3)),
    logHWFI = log(HWFI),
    logVMAX_OP_T0 = log(VMAX_OP_T0),
    logHWRF = log(HWRF)
  ) |>
  mutate(
    across(where(is.numeric) & !c(VMAX, 
                                  # HWRF, 
                                  #logHWRF
                                  # StormElapsedTime,
                                  # LAT, LON
    )
    )
  )
str(StormdataTest4unscale)

#### Test 5 ----
StormdataTest5 <- StormdataTest |>
  select(
    "StormID",
    #Date,
    Year,
    Month,
    Day,
    StormElapsedTime,
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
    LON2 = LON,
    LAT2 = LAT,
    StormElapsedTime2 = StormElapsedTime,
    arcsinhMINSLP = log(MINSLP + sqrt(MINSLP^2 + 1)),
    arcsinhSHR_MAG = log(SHR_MAG + sqrt(SHR_MAG^2 + 1)),
    arcsinhSTM_SPD = log(STM_SPD + sqrt(STM_SPD^2 + 1)),
    arcsinhSST = log(SST + sqrt(SST^2 + 1)),
    arcsinhRHLO = log(RHLO + sqrt(RHLO^2 + 1)),
    arcsinhCAPE1 = log(CAPE1 + sqrt(CAPE1^2 + 1)),
    arcsinhCAPE3 = log(CAPE3 + sqrt(CAPE3^2 + 1)),
    arcsinhSHTFL2 = log(SHTFL2 + sqrt(SHTFL2^2 + 1)),
    arcsinhTCOND7002 = log(TCOND7002 + sqrt(TCOND7002^2 + 1)),
    arcsinhINST2 = log(INST2 + sqrt(INST2^2 + 1)),
    arcsinhCP1 = log(CP1 + sqrt(CP1^2 + 1)),
    arcsinhTCONDSYM2 = log(TCONDSYM2 + sqrt(TCONDSYM2^2 + 1)),
    arcsinhCOUPLSYM3 = log(COUPLSYM3 + sqrt(COUPLSYM3^2 + 1)),
    arcsinhHWFI = log(HWFI + sqrt(HWFI^2 + 1)),
    arcsinhVMAX_OP_T0 = log(VMAX_OP_T0 + sqrt(VMAX_OP_T0^2 + 1)),
    arcsinhHWRF = log(HWRF + sqrt(HWRF^2 + 1))
  ) |>
  mutate(
    across(where(is.numeric) & !c(VMAX, StormElapsedTime2, LAT2, LON2
    ),
    function(x){scale(x,
                      center = attr(StormdataTrain5 |> pull(x), "scaled:center"),
                      scale = attr(StormdataTrain5 |> pull(x), "scaled:scale"))
    })
  )
str(StormdataTest5)

#### Test 6 ----
StormdataTest6 <- StormdataTest |>
  select(
    "StormID",
    #Date,
    Year,
    Month,
    Day,
    StormElapsedTime,
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
    LON2 = LON,
    LAT2 = LAT,
    StormElapsedTime2 = StormElapsedTime,
    arcsinhLON = log(LON + sqrt(LON^2 + 1)),
    arcsinhLAT = log(LAT + sqrt(LAT^2 + 1)),
    arcsinhStormElapsedTime = log(StormElapsedTime + sqrt(StormElapsedTime)),
    arcsinhMINSLP = log(MINSLP + sqrt(MINSLP^2 + 1)),
    arcsinhSHR_MAG = log(SHR_MAG + sqrt(SHR_MAG^2 + 1)),
    arcsinhSTM_SPD = log(STM_SPD + sqrt(STM_SPD^2 + 1)),
    arcsinhSST = log(SST + sqrt(SST^2 + 1)),
    arcsinhRHLO = log(RHLO + sqrt(RHLO^2 + 1)),
    arcsinhCAPE1 = log(CAPE1 + sqrt(CAPE1^2 + 1)),
    arcsinhCAPE3 = log(CAPE3 + sqrt(CAPE3^2 + 1)),
    arcsinhSHTFL2 = log(SHTFL2 + sqrt(SHTFL2^2 + 1)),
    arcsinhTCOND7002 = log(TCOND7002 + sqrt(TCOND7002^2 + 1)),
    arcsinhINST2 = log(INST2 + sqrt(INST2^2 + 1)),
    arcsinhCP1 = log(CP1 + sqrt(CP1^2 + 1)),
    arcsinhTCONDSYM2 = log(TCONDSYM2 + sqrt(TCONDSYM2^2 + 1)),
    arcsinhCOUPLSYM3 = log(COUPLSYM3 + sqrt(COUPLSYM3^2 + 1)),
    arcsinhHWFI = log(HWFI + sqrt(HWFI^2 + 1)),
    arcsinhVMAX_OP_T0 = log(VMAX_OP_T0 + sqrt(VMAX_OP_T0^2 + 1)),
    arcsinhHWRF = log(HWRF + sqrt(HWRF^2 + 1))
  ) |>
  mutate(
    across(where(is.numeric) & !c(VMAX, StormElapsedTime2, LAT2, LON2
    ),
    function(x){scale(x,
                      center = attr(StormdataTrain6 |> pull(x), "scaled:center"),
                      scale = attr(StormdataTrain6 |> pull(x), "scaled:scale"))
    })
  )
str(StormdataTest6)

## Actual ----
Actual_Y <- fread("_data/Actual Y.csv")
Actual_Yvec <- Actual_Y |> filter(complete.cases(x)) |> pull(x)

### Total Data ----
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
## Correlations ----
### Normal vs Normal ----
ggPairPlot <- ggpairs(StormdataTrain |>
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
                          "VMAX")) +
  labs(title = "Normal VMAX vs Normal Predictors")
ggPairPlot

### Log vs Normal ----
ggPairPlotlogVMAX <- ggpairs(StormdataTrain |>
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
                                 "VMAX") |>
                               mutate(
                                 VMAX = log(VMAX)
                               )) +
  labs(title = "Log VMAX vs Normal Predictors")
ggPairPlotlogVMAX

### Normal vs Log ----
ggPairPlotlogPreds <- ggpairs(StormdataTrain4unscale |>
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
                                  "Land",
                                  "logMINSLP",
                                  "logSHR_MAG",
                                  "logSTM_SPD",
                                  "logSST",
                                  "logRHLO",
                                  "logCAPE1",
                                  "logCAPE3",
                                  "logSHTFL2",
                                  "logTCOND7002",
                                  "logINST2",
                                  "logCP1",
                                  "logTCONDSYM2",
                                  "logCOUPLSYM3",
                                  "logHWFI",
                                  "logVMAX_OP_T0",
                                  "logHWRF",
                                  "VMAX")) +
  labs(title = "Normal VMAX vs Log Predictors")
ggPairPlotlogPreds

### Log vs Log ----
ggPairPlotlogBoth <- ggpairs(StormdataTrain4unscale |>
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
                                 "Land",
                                 "logMINSLP",
                                 "logSHR_MAG",
                                 "logSTM_SPD",
                                 "logSST",
                                 "logRHLO",
                                 "logCAPE1",
                                 "logCAPE3",
                                 "logSHTFL2",
                                 "logTCOND7002",
                                 "logINST2",
                                 "logCP1",
                                 "logTCONDSYM2",
                                 "logCOUPLSYM3",
                                 "logHWFI",
                                 "logVMAX_OP_T0",
                                 "logHWRF",
                                 "VMAX") |>
                               mutate(
                                 VMAX = log(VMAX)
                               )) +
  labs(title = "Log VMAX vs Log Predictors")
ggPairPlotlogBoth

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

mean_var_plot <- StormdataTrain4 %>%
  group_by(StormID) %>%
  summarise(mean_wind = mean(VMAX),
            var_wind = var(VMAX))

plot(mean_var_plot$mean_wind, mean_var_plot$var_wind, 
     main = "Variance vs Mean Wind Speed", 
     xlab = "Mean Wind Speed", ylab = "Variance")

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

ggplot() + 
  geom_point(
    data = StormdataTrain,
    aes(x = LON, y = LAT, 
        color = VMAX,
        shape = Land
    )
  ) +
  # xlim(c(-190,0)) +
  # ylim(c(0,60)) +
  scale_color_continuous(low = "green", high = "red") +
  facet_wrap(vars(StormElapsedTime)) +
  theme_bw()

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

# Load prior models ----
modelFiles <- list.files("~/Desktop/Temp Hurricane Model Data")
for(i in modelFiles){
  load(file = paste0("~/Desktop/Temp Hurricane Model Data/", i))
}

# LOGNORMAL ----
## NULL Model ----
logNormalFitNULL <- brm(
  bf(
    VMAX ~ 1
  ),
  data = StormdataTrain4, 
  family = lognormal(link = "identity"),
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 5000,
  seed = 52,
  warmup = 2500
)

### Plot ----
print(logNormalFitNULL, digits = 4)
logNormalFitNULLppcFit <- pp_check(logNormalFitNULL, ndraws = 100) + 
  labs(title = "logNormalFitNULL Fit PPC") +
  theme_bw()
logNormalFitNULLppcFit

### Bayes R2 ----
logNormalFitNULLR2 <- bayes_R2(logNormalFitNULL) |>
  bind_cols(Fit = "logNormalFitNULL") |>
  select(Fit, everything())

### LOO ----
logNormalFitNULLloo <- loo(logNormalFitNULL)

### Prediction ----
## Fitted
logNormalFitNULLFit <- posterior_predict(logNormalFitNULL)
logNormalFitNULLFitMean <- colMeans(logNormalFitNULLFit)
logNormalFitNULLFitMed <- apply(logNormalFitNULLFit, 2, function(x){quantile(x, 0.5)})
logNormalFitNULLFitLCB <- apply(logNormalFitNULLFit, 2, function(x){quantile(x, 0.025)})
logNormalFitNULLFitUCB <- apply(logNormalFitNULLFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
logNormalFitNULLPreds <- posterior_predict(logNormalFitNULL, 
                                           newdata = StormdataTest4,
                                           allow_new_levels = TRUE)
logNormalFitNULLPredsMean <- colMeans(logNormalFitNULLPreds)
logNormalFitNULLPredsMed <- apply(logNormalFitNULLPreds, 2, function(x){quantile(x, 0.5)})
logNormalFitNULLPredsLCB <- apply(logNormalFitNULLPreds, 2, function(x){quantile(x, 0.025)})
logNormalFitNULLPredsUCB <- apply(logNormalFitNULLPreds, 2, function(x){quantile(x, 0.975)})

logNormalFitNULLpredMetrics <- tibble(
  Fit = "logNormalFitNULL",
  MAE_HWRF_fit = mean(abs(StormdataTrain$HWRF - StormdataTrain$VMAX)),
  MAE_fit = mean(abs(logNormalFitNULLFitMean - StormdataTrain$VMAX)),
  COV_fit = mean(logNormalFitNULLFitLCB < StormdataTrain$VMAX & StormdataTrain$VMAX < logNormalFitNULLFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(logNormalFitNULLPredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(logNormalFitNULLPredsMed - Actual_Yvec)),
  COV_pred = mean(logNormalFitNULLPredsLCB < Actual_Yvec & Actual_Yvec < logNormalFitNULLPredsUCB)
)
logNormalFitNULLpredMetrics

### PPC ----
###### Quantile 2.5 
logNormalFitNULLLCBsims <- apply(logNormalFitNULLFit, 
                                 MARGIN = 1,
                                 function(x){
                                   quantile(x, 0.025)
                                 })
logNormalFitNULLLCBpvalueVec <- logNormalFitNULLLCBsims < quantile(StormdataTrain$VMAX, 0.025)
logNormalFitNULLLCBpvalue <- sum(logNormalFitNULLLCBpvalueVec)
logNormalFitNULLLCBpvalue <- round(logNormalFitNULLLCBpvalue/10000, 3)
logNormalFitNULLLCBpvalue <- min(logNormalFitNULLLCBpvalue, 1 - logNormalFitNULLLCBpvalue)

logNormalFitNULL_ppcLCB <- 
  ppc_stat(StormdataTrain$VMAX,
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
logNormalFitNULLUCBpvalueVec <- logNormalFitNULLUCBsims < quantile(StormdataTrain$VMAX, 0.975)
logNormalFitNULLUCBpvalue <- as.numeric(sum(logNormalFitNULLUCBpvalueVec))
logNormalFitNULLUCBpvalue <- round(logNormalFitNULLUCBpvalue/10000, 3)
logNormalFitNULLUCBpvalue <- min(logNormalFitNULLUCBpvalue, 1 - logNormalFitNULLUCBpvalue)

logNormalFitNULL_ppcUCB <- 
  ppc_stat(StormdataTrain$VMAX,
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
logNormalFitNULLMEANpvalueVec <- logNormalFitNULLMEANsims < mean(StormdataTrain$VMAX)
logNormalFitNULLMEANpvalue <- sum(logNormalFitNULLMEANpvalueVec)
logNormalFitNULLMEANpvalue <- round(logNormalFitNULLMEANpvalue/10000, 3)
logNormalFitNULLMEANpvalue <- min(logNormalFitNULLMEANpvalue, 1 - logNormalFitNULLMEANpvalue)

logNormalFitNULL_ppcMEAN <- 
  ppc_stat(StormdataTrain$VMAX,
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
logNormalFitNULLMEDpvalueVec <- logNormalFitNULLMEDsims < quantile(StormdataTrain$VMAX, 0.5)
logNormalFitNULLMEDpvalue <- sum(logNormalFitNULLMEDpvalueVec)
logNormalFitNULLMEDpvalue <- round(logNormalFitNULLMEDpvalue/10000, 3)
logNormalFitNULLMEDpvalue <- min(logNormalFitNULLMEDpvalue, 1 - logNormalFitNULLMEDpvalue)

logNormalFitNULL_ppcMED <- 
  ppc_stat(StormdataTrain$VMAX,
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
logNormalFitNULLSDpvalueVec <- logNormalFitNULLSDsims < sd(StormdataTrain$VMAX)
logNormalFitNULLSDpvalue <- sum(logNormalFitNULLSDpvalueVec)
logNormalFitNULLSDpvalue <- round(logNormalFitNULLSDpvalue/10000, 3)
logNormalFitNULLSDpvalue <- min(logNormalFitNULLSDpvalue, 1 - logNormalFitNULLSDpvalue)

logNormalFitNULL_ppcSD <- 
  ppc_stat(StormdataTrain$VMAX,
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
logNormalFitNULLRANGEpvalueVec <- logNormalFitNULLRANGEsims < (max(StormdataTrain$VMAX)-min(StormdataTrain3$VMAX))
logNormalFitNULLRANGEpvalue <- sum(logNormalFitNULLRANGEpvalueVec)
logNormalFitNULLRANGEpvalue <- round(logNormalFitNULLRANGEpvalue/10000, 3)
logNormalFitNULLRANGEpvalue <- min(logNormalFitNULLRANGEpvalue, 1 - logNormalFitNULLRANGEpvalue)

logNormalFitNULL_ppcRANGE <- 
  ppc_stat(StormdataTrain$VMAX,
           logNormalFitNULLFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", logNormalFitNULLRANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitNULL_ppcRANGE

### Combined Plot ----
logNormalFitNULL_ppcComb <- 
  logNormalFitNULLppcFit /
  (logNormalFitNULL_ppcLCB | logNormalFitNULL_ppcMED | logNormalFitNULL_ppcUCB) /
  (logNormalFitNULL_ppcRANGE | logNormalFitNULL_ppcMEAN | logNormalFitNULL_ppcSD)
logNormalFitNULL_ppcComb

### Bayes p-values ----
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

### CV ----
set.seed(52)
kfoldID <- kfold_split_grouped(K = 5, StormdataTrain$StormID)
logNormalFitNULLkfoldgroup <- kfold(logNormalFitNULL,
                                    folds = kfoldID,
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


### Clear Environment ----
NUllenv <- ls()
exNULL <- str_detect(NUllenv, "NULL")
keepNULL <- NUllenv %in% c("logNormalFitNULLppcFitComb",
                           "logNormalFitNULLpvalues",
                           "logNormalFitNULLloo",
                           "logNormalFitNULLpredMetrics",
                           "logNormalFitNULLkfoldMetrics",
                           "logNormalFitNULL")
cutNULL <- exNULL & !keepNULL
NULLrm <- NUllenv[cutNULL]
rm(list = NULLrm)

## HWRF Model ----
logNormalFitHWRF2 <- brm(
  bf(
    VMAX ~ logHWRF
  ),
  data = StormdataTrain4, 
  family = lognormal(link = "identity"),
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 5000,
  seed = 52,
  warmup = 2500
)

### Plot ----
print(logNormalFitHWRF, digits = 4)
logNormalFitHWRFppcFit <- pp_check(logNormalFitHWRF, ndraws = 100) + 
  labs(title = "logNormalFitHWRF Fit PPC") +
  theme_bw()
logNormalFitHWRFppcFit

### Bayes R2 ----
logNormalFitHWRFR2 <- bayes_R2(logNormalFitHWRF) |>
  bind_cols(Fit = "logNormalFitHWRF") |>
  select(Fit, everything())

### LOO ----
logNormalFitHWRFloo <- loo(logNormalFitHWRF)

### Prediction ----
## Fitted
logNormalFitHWRFFit <- posterior_predict(logNormalFitHWRF)
logNormalFitHWRFFitMean <- colMeans(logNormalFitHWRFFit)
logNormalFitHWRFFitMed <- apply(logNormalFitHWRFFit, 2, function(x){quantile(x, 0.5)})
logNormalFitHWRFFitLCB <- apply(logNormalFitHWRFFit, 2, function(x){quantile(x, 0.025)})
logNormalFitHWRFFitUCB <- apply(logNormalFitHWRFFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
logNormalFitHWRFPreds <- posterior_predict(logNormalFitHWRF, 
                                           newdata = StormdataTest4,
                                           allow_new_levels = TRUE)
logNormalFitHWRFPredsMean <- colMeans(logNormalFitHWRFPreds)
logNormalFitHWRFPredsMed <- apply(logNormalFitHWRFPreds, 2, function(x){quantile(x, 0.5)})
logNormalFitHWRFPredsLCB <- apply(logNormalFitHWRFPreds, 2, function(x){quantile(x, 0.025)})
logNormalFitHWRFPredsUCB <- apply(logNormalFitHWRFPreds, 2, function(x){quantile(x, 0.975)})

logNormalFitHWRFpredMetrics <- tibble(
  Fit = "logNormalFitHWRF",
  MAE_HWRF_fit = mean(abs(StormdataTrain$HWRF - StormdataTrain$VMAX)),
  MAE_fit = mean(abs(logNormalFitHWRFFitMean - StormdataTrain$VMAX)),
  COV_fit = mean(logNormalFitHWRFFitLCB < StormdataTrain$VMAX & StormdataTrain$VMAX < logNormalFitHWRFFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(logNormalFitHWRFPredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(logNormalFitHWRFPredsMed - Actual_Yvec)),
  COV_pred = mean(logNormalFitHWRFPredsLCB < Actual_Yvec & Actual_Yvec < logNormalFitHWRFPredsUCB)
)
logNormalFitHWRFpredMetrics

### PPC ----
###### Quantile 2.5 
logNormalFitHWRFLCBsims <- apply(logNormalFitHWRFFit, 
                                 MARGIN = 1,
                                 function(x){
                                   quantile(x, 0.025)
                                 })
logNormalFitHWRFLCBpvalueVec <- logNormalFitHWRFLCBsims < quantile(StormdataTrain$VMAX, 0.025)
logNormalFitHWRFLCBpvalue <- sum(logNormalFitHWRFLCBpvalueVec)
logNormalFitHWRFLCBpvalue <- round(logNormalFitHWRFLCBpvalue/10000, 3)
logNormalFitHWRFLCBpvalue <- min(logNormalFitHWRFLCBpvalue, 1 - logNormalFitHWRFLCBpvalue)

logNormalFitHWRF_ppcLCB <- 
  ppc_stat(StormdataTrain$VMAX,
           logNormalFitHWRFFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", logNormalFitHWRFLCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitHWRF_ppcLCB

###### Quantile 97.5 
logNormalFitHWRFUCBsims <- apply(logNormalFitHWRFFit, 
                                 MARGIN = 1,
                                 function(x){
                                   quantile(x, 0.975)
                                 })
logNormalFitHWRFUCBpvalueVec <- logNormalFitHWRFUCBsims < quantile(StormdataTrain$VMAX, 0.975)
logNormalFitHWRFUCBpvalue <- as.numeric(sum(logNormalFitHWRFUCBpvalueVec))
logNormalFitHWRFUCBpvalue <- round(logNormalFitHWRFUCBpvalue/10000, 3)
logNormalFitHWRFUCBpvalue <- min(logNormalFitHWRFUCBpvalue, 1 - logNormalFitHWRFUCBpvalue)

logNormalFitHWRF_ppcUCB <- 
  ppc_stat(StormdataTrain$VMAX,
           logNormalFitHWRFFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", logNormalFitHWRFUCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitHWRF_ppcUCB

###### Mean 
logNormalFitHWRFMEANsims <- apply(logNormalFitHWRFFit, 
                                  MARGIN = 1,
                                  function(x){
                                    mean(x)
                                  })
logNormalFitHWRFMEANpvalueVec <- logNormalFitHWRFMEANsims < mean(StormdataTrain$VMAX)
logNormalFitHWRFMEANpvalue <- sum(logNormalFitHWRFMEANpvalueVec)
logNormalFitHWRFMEANpvalue <- round(logNormalFitHWRFMEANpvalue/10000, 3)
logNormalFitHWRFMEANpvalue <- min(logNormalFitHWRFMEANpvalue, 1 - logNormalFitHWRFMEANpvalue)

logNormalFitHWRF_ppcMEAN <- 
  ppc_stat(StormdataTrain$VMAX,
           logNormalFitHWRFFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", logNormalFitHWRFMEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitHWRF_ppcMEAN

###### Med 
logNormalFitHWRFMEDsims <- apply(logNormalFitHWRFFit, 
                                 MARGIN = 1,
                                 function(x){
                                   quantile(x, 0.5)
                                 })
logNormalFitHWRFMEDpvalueVec <- logNormalFitHWRFMEDsims < quantile(StormdataTrain$VMAX, 0.5)
logNormalFitHWRFMEDpvalue <- sum(logNormalFitHWRFMEDpvalueVec)
logNormalFitHWRFMEDpvalue <- round(logNormalFitHWRFMEDpvalue/10000, 3)
logNormalFitHWRFMEDpvalue <- min(logNormalFitHWRFMEDpvalue, 1 - logNormalFitHWRFMEDpvalue)

logNormalFitHWRF_ppcMED <- 
  ppc_stat(StormdataTrain$VMAX,
           logNormalFitHWRFFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", logNormalFitHWRFMEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitHWRF_ppcMED

###### SD 
logNormalFitHWRFSDsims <- apply(logNormalFitHWRFFit, 
                                MARGIN = 1,
                                function(x){
                                  sd(x)
                                })
logNormalFitHWRFSDpvalueVec <- logNormalFitHWRFSDsims < sd(StormdataTrain$VMAX)
logNormalFitHWRFSDpvalue <- sum(logNormalFitHWRFSDpvalueVec)
logNormalFitHWRFSDpvalue <- round(logNormalFitHWRFSDpvalue/10000, 3)
logNormalFitHWRFSDpvalue <- min(logNormalFitHWRFSDpvalue, 1 - logNormalFitHWRFSDpvalue)

logNormalFitHWRF_ppcSD <- 
  ppc_stat(StormdataTrain$VMAX,
           logNormalFitHWRFFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", logNormalFitHWRFSDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitHWRF_ppcSD

###### Range 
logNormalFitHWRFRANGEsims <- apply(logNormalFitHWRFFit, 
                                   MARGIN = 1,
                                   function(x){
                                     max(x)-min(x)
                                   })
logNormalFitHWRFRANGEpvalueVec <- logNormalFitHWRFRANGEsims < (max(StormdataTrain$VMAX)-min(StormdataTrain3$VMAX))
logNormalFitHWRFRANGEpvalue <- sum(logNormalFitHWRFRANGEpvalueVec)
logNormalFitHWRFRANGEpvalue <- round(logNormalFitHWRFRANGEpvalue/10000, 3)
logNormalFitHWRFRANGEpvalue <- min(logNormalFitHWRFRANGEpvalue, 1 - logNormalFitHWRFRANGEpvalue)

logNormalFitHWRF_ppcRANGE <- 
  ppc_stat(StormdataTrain$VMAX,
           logNormalFitHWRFFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", logNormalFitHWRFRANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitHWRF_ppcRANGE

### Combined Plot ----
logNormalFitHWRF_ppcComb <- 
  logNormalFitHWRFppcFit /
  (logNormalFitHWRF_ppcLCB | logNormalFitHWRF_ppcMED | logNormalFitHWRF_ppcUCB) /
  (logNormalFitHWRF_ppcRANGE | logNormalFitHWRF_ppcMEAN | logNormalFitHWRF_ppcSD)
logNormalFitHWRF_ppcComb

### Bayes p-values ----
logNormalFitHWRFpvalues <- tibble(
  Fit = paste0("logNormalFitHWRF"),
  LCB = logNormalFitHWRFLCBpvalue,
  Median = logNormalFitHWRFMEDpvalue,
  UCB = logNormalFitHWRFUCBpvalue,
  Range = logNormalFitHWRFRANGEpvalue,
  Mean = logNormalFitHWRFMEANpvalue,
  SD = logNormalFitHWRFSDpvalue
)
logNormalFitHWRFpvalues

### CV ----
set.seed(52)
kfoldID <- kfold_split_grouped(K = 5, StormdataTrain$StormID)
logNormalFitHWRFkfoldgroup <- kfold(logNormalFitHWRF,
                                    folds = kfoldID,
                                    chains = 1,
                                    save_fits = TRUE)
logNormalFitHWRFkfoldPreds <- kfold_predict(logNormalFitHWRFkfoldgroup)
logNormalFitHWRFkfoldPredsDat <- logNormalFitHWRFkfoldPreds$yrep
logNormalFitHWRFkfoldPredsMean <- colMeans(logNormalFitHWRFkfoldPredsDat)
logNormalFitHWRFkfoldPredsMed <- apply(logNormalFitHWRFkfoldPredsDat, 2, function(x){quantile(x, 0.5)})
logNormalFitHWRFkfoldPredsLCB <- apply(logNormalFitHWRFkfoldPredsDat, 2, function(x){quantile(x, 0.025)})
logNormalFitHWRFkfoldPredsUCB <- apply(logNormalFitHWRFkfoldPredsDat, 2, function(x){quantile(x, 0.975)})

logNormalFitHWRFkfoldMetrics <- tibble(
  Fit = paste0("logNormalFitHWRF"),
  MAE_kfold = mean(abs(logNormalFitHWRFkfoldPredsMean - logNormalFitHWRFkfoldPreds$y)),
  MAD_kfold = mean(abs(logNormalFitHWRFkfoldPredsMed - logNormalFitHWRFkfoldPreds$y)),
  COV_kfold = mean(logNormalFitHWRFkfoldPredsLCB < logNormalFitHWRFkfoldPreds$y & logNormalFitHWRFkfoldPreds$y < logNormalFitHWRFkfoldPredsUCB)
)
logNormalFitHWRFkfoldMetrics

### Clear Environment ----
NUllenv <- ls()
exHWRF <- str_detect(NUllenv, "HWRF")
keepHWRF <- NUllenv %in% c("logNormalFitHWRFppcFitComb",
                           "logNormalFitHWRFpvalues",
                           "logNormalFitHWRFloo",
                           "logNormalFitHWRFpredMetrics",
                           "logNormalFitHWRFkfoldMetrics",
                           "logNormalFitHWRF")
cutHWRF <- exHWRF & !keepHWRF
HWRFrm <- NUllenv[cutHWRF]
rm(list = HWRFrm)

## HWRF StormID Model ----
logNormalFitHWRFstormID2 <- brm(
  bf(
    VMAX ~ logHWRF +
      (1|StormID)
  ),
  data = StormdataTrain4, 
  family = lognormal(link = "identity"),
  save_pars = save_pars(all = TRUE), 
  chains = 4,
  iter = 5000,
  seed = 52,
  warmup = 2500
)

### Plot ----
print(logNormalFitHWRFstormID, digits = 4)
logNormalFitHWRFstormIDppcFit <- pp_check(logNormalFitHWRFstormID, ndraws = 100) + 
  labs(title = "logNormalFitHWRFstormID Fit PPC") +
  theme_bw()
logNormalFitHWRFstormIDppcFit

### Bayes R2 ----
logNormalFitHWRFstormIDR2 <- bayes_R2(logNormalFitHWRFstormID) |>
  bind_cols(Fit = "logNormalFitHWRFstormID") |>
  select(Fit, everything())

### LOO ----
logNormalFitHWRFstormIDloo <- loo(logNormalFitHWRFstormID)

### Prediction ----
## Fitted
logNormalFitHWRFstormIDFit <- posterior_predict(logNormalFitHWRFstormID)
logNormalFitHWRFstormIDFitMean <- colMeans(logNormalFitHWRFstormIDFit)
logNormalFitHWRFstormIDFitMed <- apply(logNormalFitHWRFstormIDFit, 2, function(x){quantile(x, 0.5)})
logNormalFitHWRFstormIDFitLCB <- apply(logNormalFitHWRFstormIDFit, 2, function(x){quantile(x, 0.025)})
logNormalFitHWRFstormIDFitUCB <- apply(logNormalFitHWRFstormIDFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
logNormalFitHWRFstormIDPreds <- posterior_predict(logNormalFitHWRFstormID, 
                                                  newdata = StormdataTest4,
                                                  allow_new_levels = TRUE,
                                                  re_formula = NULL)
logNormalFitHWRFstormIDPredsMean <- colMeans(logNormalFitHWRFstormIDPreds)
logNormalFitHWRFstormIDPredsMed <- apply(logNormalFitHWRFstormIDPreds, 2, function(x){quantile(x, 0.5)})
logNormalFitHWRFstormIDPredsLCB <- apply(logNormalFitHWRFstormIDPreds, 2, function(x){quantile(x, 0.025)})
logNormalFitHWRFstormIDPredsUCB <- apply(logNormalFitHWRFstormIDPreds, 2, function(x){quantile(x, 0.975)})

logNormalFitHWRFstormIDpredMetrics <- tibble(
  Fit = "logNormalFitHWRFstormID",
  MAE_HWRF_fit = mean(abs(StormdataTrain$HWRF - StormdataTrain$VMAX)),
  MAE_fit = mean(abs(logNormalFitHWRFstormIDFitMean - StormdataTrain$VMAX)),
  COV_fit = mean(logNormalFitHWRFstormIDFitLCB < StormdataTrain$VMAX & StormdataTrain$VMAX < logNormalFitHWRFstormIDFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(logNormalFitHWRFstormIDPredsMean - Actual_Yvec)),
  MAD_pred = mean(abs(logNormalFitHWRFstormIDPredsMed - Actual_Yvec)),
  COV_pred = mean(logNormalFitHWRFstormIDPredsLCB < Actual_Yvec & Actual_Yvec < logNormalFitHWRFstormIDPredsUCB)
)
logNormalFitHWRFstormIDpredMetrics

### PPC ----
###### Quantile 2.5 
logNormalFitHWRFstormIDLCBsims <- apply(logNormalFitHWRFstormIDFit, 
                                        MARGIN = 1,
                                        function(x){
                                          quantile(x, 0.025)
                                        })
logNormalFitHWRFstormIDLCBpvalueVec <- logNormalFitHWRFstormIDLCBsims < quantile(StormdataTrain$VMAX, 0.025)
logNormalFitHWRFstormIDLCBpvalue <- sum(logNormalFitHWRFstormIDLCBpvalueVec)
logNormalFitHWRFstormIDLCBpvalue <- round(logNormalFitHWRFstormIDLCBpvalue/10000, 3)
logNormalFitHWRFstormIDLCBpvalue <- min(logNormalFitHWRFstormIDLCBpvalue, 1 - logNormalFitHWRFstormIDLCBpvalue)

logNormalFitHWRFstormID_ppcLCB <- 
  ppc_stat(StormdataTrain$VMAX,
           logNormalFitHWRFstormIDFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", logNormalFitHWRFstormIDLCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitHWRFstormID_ppcLCB

###### Quantile 97.5 
logNormalFitHWRFstormIDUCBsims <- apply(logNormalFitHWRFstormIDFit, 
                                        MARGIN = 1,
                                        function(x){
                                          quantile(x, 0.975)
                                        })
logNormalFitHWRFstormIDUCBpvalueVec <- logNormalFitHWRFstormIDUCBsims < quantile(StormdataTrain$VMAX, 0.975)
logNormalFitHWRFstormIDUCBpvalue <- as.numeric(sum(logNormalFitHWRFstormIDUCBpvalueVec))
logNormalFitHWRFstormIDUCBpvalue <- round(logNormalFitHWRFstormIDUCBpvalue/10000, 3)
logNormalFitHWRFstormIDUCBpvalue <- min(logNormalFitHWRFstormIDUCBpvalue, 1 - logNormalFitHWRFstormIDUCBpvalue)

logNormalFitHWRFstormID_ppcUCB <- 
  ppc_stat(StormdataTrain$VMAX,
           logNormalFitHWRFstormIDFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", logNormalFitHWRFstormIDUCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitHWRFstormID_ppcUCB

###### Mean 
logNormalFitHWRFstormIDMEANsims <- apply(logNormalFitHWRFstormIDFit, 
                                         MARGIN = 1,
                                         function(x){
                                           mean(x)
                                         })
logNormalFitHWRFstormIDMEANpvalueVec <- logNormalFitHWRFstormIDMEANsims < mean(StormdataTrain$VMAX)
logNormalFitHWRFstormIDMEANpvalue <- sum(logNormalFitHWRFstormIDMEANpvalueVec)
logNormalFitHWRFstormIDMEANpvalue <- round(logNormalFitHWRFstormIDMEANpvalue/10000, 3)
logNormalFitHWRFstormIDMEANpvalue <- min(logNormalFitHWRFstormIDMEANpvalue, 1 - logNormalFitHWRFstormIDMEANpvalue)

logNormalFitHWRFstormID_ppcMEAN <- 
  ppc_stat(StormdataTrain$VMAX,
           logNormalFitHWRFstormIDFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", logNormalFitHWRFstormIDMEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitHWRFstormID_ppcMEAN

###### Med 
logNormalFitHWRFstormIDMEDsims <- apply(logNormalFitHWRFstormIDFit, 
                                        MARGIN = 1,
                                        function(x){
                                          quantile(x, 0.5)
                                        })
logNormalFitHWRFstormIDMEDpvalueVec <- logNormalFitHWRFstormIDMEDsims < quantile(StormdataTrain$VMAX, 0.5)
logNormalFitHWRFstormIDMEDpvalue <- sum(logNormalFitHWRFstormIDMEDpvalueVec)
logNormalFitHWRFstormIDMEDpvalue <- round(logNormalFitHWRFstormIDMEDpvalue/10000, 3)
logNormalFitHWRFstormIDMEDpvalue <- min(logNormalFitHWRFstormIDMEDpvalue, 1 - logNormalFitHWRFstormIDMEDpvalue)

logNormalFitHWRFstormID_ppcMED <- 
  ppc_stat(StormdataTrain$VMAX,
           logNormalFitHWRFstormIDFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", logNormalFitHWRFstormIDMEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitHWRFstormID_ppcMED

###### SD 
logNormalFitHWRFstormIDSDsims <- apply(logNormalFitHWRFstormIDFit, 
                                       MARGIN = 1,
                                       function(x){
                                         sd(x)
                                       })
logNormalFitHWRFstormIDSDpvalueVec <- logNormalFitHWRFstormIDSDsims < sd(StormdataTrain$VMAX)
logNormalFitHWRFstormIDSDpvalue <- sum(logNormalFitHWRFstormIDSDpvalueVec)
logNormalFitHWRFstormIDSDpvalue <- round(logNormalFitHWRFstormIDSDpvalue/10000, 3)
logNormalFitHWRFstormIDSDpvalue <- min(logNormalFitHWRFstormIDSDpvalue, 1 - logNormalFitHWRFstormIDSDpvalue)

logNormalFitHWRFstormID_ppcSD <- 
  ppc_stat(StormdataTrain$VMAX,
           logNormalFitHWRFstormIDFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", logNormalFitHWRFstormIDSDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitHWRFstormID_ppcSD

###### Range 
logNormalFitHWRFstormIDRANGEsims <- apply(logNormalFitHWRFstormIDFit, 
                                          MARGIN = 1,
                                          function(x){
                                            max(x)-min(x)
                                          })
logNormalFitHWRFstormIDRANGEpvalueVec <- logNormalFitHWRFstormIDRANGEsims < (max(StormdataTrain$VMAX)-min(StormdataTrain3$VMAX))
logNormalFitHWRFstormIDRANGEpvalue <- sum(logNormalFitHWRFstormIDRANGEpvalueVec)
logNormalFitHWRFstormIDRANGEpvalue <- round(logNormalFitHWRFstormIDRANGEpvalue/10000, 3)
logNormalFitHWRFstormIDRANGEpvalue <- min(logNormalFitHWRFstormIDRANGEpvalue, 1 - logNormalFitHWRFstormIDRANGEpvalue)

logNormalFitHWRFstormID_ppcRANGE <- 
  ppc_stat(StormdataTrain$VMAX,
           logNormalFitHWRFstormIDFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", logNormalFitHWRFstormIDRANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitHWRFstormID_ppcRANGE

### Combined Plot ----
logNormalFitHWRFstormID_ppcComb <- 
  logNormalFitHWRFstormIDppcFit /
  (logNormalFitHWRFstormID_ppcLCB | logNormalFitHWRFstormID_ppcMED | logNormalFitHWRFstormID_ppcUCB) /
  (logNormalFitHWRFstormID_ppcRANGE | logNormalFitHWRFstormID_ppcMEAN | logNormalFitHWRFstormID_ppcSD)
logNormalFitHWRFstormID_ppcComb

### Bayes p-values ----
logNormalFitHWRFstormIDpvalues <- tibble(
  Fit = paste0("logNormalFitHWRFstormID"),
  LCB = logNormalFitHWRFstormIDLCBpvalue,
  Median = logNormalFitHWRFstormIDMEDpvalue,
  UCB = logNormalFitHWRFstormIDUCBpvalue,
  Range = logNormalFitHWRFstormIDRANGEpvalue,
  Mean = logNormalFitHWRFstormIDMEANpvalue,
  SD = logNormalFitHWRFstormIDSDpvalue
)
logNormalFitHWRFstormIDpvalues

### CV ----
set.seed(52)
kfoldID <- kfold_split_grouped(K = 5, StormdataTrain$StormID)
logNormalFitHWRFstormIDkfoldgroup <- kfold(logNormalFitHWRFstormID,
                                           folds = kfoldID,
                                           chains = 1,
                                           save_fits = TRUE)
logNormalFitHWRFstormIDkfoldPreds <- kfold_predict(logNormalFitHWRFstormIDkfoldgroup)
logNormalFitHWRFstormIDkfoldPredsDat <- logNormalFitHWRFstormIDkfoldPreds$yrep
logNormalFitHWRFstormIDkfoldPredsMean <- colMeans(logNormalFitHWRFstormIDkfoldPredsDat)
logNormalFitHWRFstormIDkfoldPredsMed <- apply(logNormalFitHWRFstormIDkfoldPredsDat, 2, function(x){quantile(x, 0.5)})
logNormalFitHWRFstormIDkfoldPredsLCB <- apply(logNormalFitHWRFstormIDkfoldPredsDat, 2, function(x){quantile(x, 0.025)})
logNormalFitHWRFstormIDkfoldPredsUCB <- apply(logNormalFitHWRFstormIDkfoldPredsDat, 2, function(x){quantile(x, 0.975)})

logNormalFitHWRFstormIDkfoldMetrics <- tibble(
  Fit = paste0("logNormalFitHWRFstormID"),
  MAE_kfold = mean(abs(logNormalFitHWRFstormIDkfoldPredsMean - logNormalFitHWRFstormIDkfoldPreds$y)),
  MAD_kfold = mean(abs(logNormalFitHWRFstormIDkfoldPredsMed - logNormalFitHWRFstormIDkfoldPreds$y)),
  COV_kfold = mean(logNormalFitHWRFstormIDkfoldPredsLCB < logNormalFitHWRFstormIDkfoldPreds$y & logNormalFitHWRFstormIDkfoldPreds$y < logNormalFitHWRFstormIDkfoldPredsUCB)
)
logNormalFitHWRFstormIDkfoldMetrics

### Clear Environment ----
NUllenv <- ls()
exHWRF <- str_detect(NUllenv, "HWRFstormID")
keepHWRF <- NUllenv %in% c("logNormalFitHWRFstormIDppcFitComb",
                           "logNormalFitHWRFstormIDpvalues",
                           "logNormalFitHWRFstormIDloo",
                           "logNormalFitHWRFstormIDpredMetrics",
                           "logNormalFitHWRFstormIDkfoldMetrics",
                           "logNormalFitHWRFstormID")
cutHWRF <- exHWRF & !keepHWRF
HWRFrm <- NUllenv[cutHWRF]
rm(list = HWRFrm)

## Compare Base Models ----
### Bayes pvalue ----
pvalueComp_logNormalBase <- bind_rows(
  logNormalFitNULLpvalues,
  logNormalFitHWRFpvalues,
  logNormalFitHWRFstormIDpvalues
)
pvalueComp_logNormalBase
save(pvalueComp_logNormalBase, 
     file = "~/Desktop/Temp Hurricane Model Data/pvalueComp_logNormalBase.RData")

### LOO ----
looComp_logNormalBase <- loo_compare(
  logNormalFitNULLloo,
  logNormalFitHWRFloo,
  logNormalFitHWRFstormIDloo
)
looComp_logNormalBase
save(looComp_logNormalBase, 
     file = "~/Desktop/Temp Hurricane Model Data/looComp_logNormalBase.RData")

### CV ----
cvComp_logNormalBase <- bind_rows(
  logNormalFitNULLkfoldMetrics,
  logNormalFitHWRFkfoldMetrics,
  logNormalFitHWRFstormIDkfoldMetrics
)
cvComp_logNormalBase <- cvComp_logNormalBase |> arrange(MAE_kfold)
cvComp_logNormalBase
save(cvComp_logNormalBase, 
     file = "~/Desktop/Temp Hurricane Model Data/cvComp_logNormalBase.RData")

### Preds ----
predComp_logNormalBase <- bind_rows(
  logNormalFitNULLpredMetrics,
  logNormalFitHWRFpredMetrics,
  logNormalFitHWRFstormIDpredMetrics
)
predComp_logNormalBase <- predComp_logNormalBase |> arrange(MAE_pred)
predComp_logNormalBase
save(predComp_logNormalBase, 
     file = "~/Desktop/Temp Hurricane Model Data/predComp_logNormalBase.RData")

### Bayes R2 ----
bayesR2_logNormalBase <- bind_rows(
  logNormalFitNULLR2,
  logNormalFitHWRFR2,
  logNormalFitHWRFstormIDR2
) |>
  arrange(desc(Estimate))
bayesR2_logNormalBase
save(bayesR2_logNormalBase, 
     file = "~/Desktop/Temp Hurricane Model Data/bayesR2_logNormalBase.RData")

### Clear Environment ----
BASEenv <- ls()
save(list = BASEenv,
     file = "~/Desktop/Temp Hurricane Model Data/logNormalBaseEnv.RData")

# Fit Model ----

### POOLED
# s(Day, bs = "cc") +
# s(StormElapsedTime, bs = "tp") +
# s(LON, bs = "tp") +
# s(LAT, bs = "tp") +
# Day + 
# StormElapsedTime +
# LAT +
# LON +
# basin +
# MINSLP +
# SHR_MAG +
# STM_SPD +
# SST +
# RHLO +
# CAPE1 +
# CAPE3 +
# SHTFL2 +
# TCOND7002 +
# INST2 +
# CP1 +
# TCONDSYM2 +
# COUPLSYM3 +
# HWFI +
# VMAX_OP_T0 +
# HWRF

### UNPOOLED

# Fit1: VMAX ~ logHWRFcenter + (1|StormID)
##      data = StormdataTrain4
# Fit2: VMAX ~ logHWRFscale + (1|StormID)
##      data = StormdataTrain4
# Fit3: VMAX ~ HWRF + (1|StormID)
##      data = StormdataTrain4
# Fit4: VMAX ~   bf(VMAX ~ 
#                   s(Day, bs = "cc") +
#                   s(StormElapsedTime, bs = "tp") +
#                   s(LON, bs = "tp") +
#                   s(LAT, bs = "tp") +
#                   logHWRF +
#                   (1|StormID)
#                   )
##      data = StormdataTrain4
# Fit5: VMAX ~   bf(VMAX ~ 
#                   Day + I(Day^2) +
#                   s(StormElapsedTime, bs = "tp") +
#                   s(LON, bs = "tp") +
#                   s(LAT, bs = "tp") +
#                   logHWRF +
#                   (1|StormID)
#                   )
##      data = StormdataTrain4
# Fit6: VMAX ~   bf(VMAX ~ 
#                   Day + I(Day^2) +
#                   s(StormElapsedTime, bs = "tp") +
#                   s(LON, bs = "tp") +
#                   s(LAT, bs = "tp") +
#                   logHWRF +
#                   (1|StormID)
#                   ),
#                   sigma ~ logHWRF
##      data = StormdataTrain4
# Fit7: VMAX ~   bf(VMAX ~ 
#                   s(Day, bs = "cc") +
#                   s(StormElapsedTime, bs = "tp") +
#                   s(LON, LAT) +
#                   logHWRF +
#                   (1 + StormElapsedTime|StormID)
#                   ),
#                   sigma ~ logHWRF
##      data = StormdataTrain4
# Fit8: VMAX ~   bf(VMAX ~ 
#                   s(StormElapsedTime, bs = "tp") +
#                   s(LON, LAT) +
#                   logHWRF +
#                   (1 + StormElapsedTime|StormID)
#                   ),
#                   sigma ~ logHWRF
##      data = StormdataTrain4
# Fit9: VMAX ~   bf(VMAX ~ 
#                   s(StormElapsedTime, bs = "tp") +
#                   t2(LON, LAT) +
#                   logHWRF +
#                   (1 + StormElapsedTime|StormID)
#                   ),
#                   sigma ~ logHWRF
##      data = StormdataTrain4
# Fit10: VMAX ~   bf(VMAX ~ 
#                   s(StormElapsedTime, bs = "tp") +
#                   s(LON, LAT, by = basin) +
#                   logHWRF +
#                   (1|StormID)
#                   ),
#                   sigma ~ logHWRF
##      data = StormdataTrain4
# Fit11: VMAX ~   bf(VMAX ~ 
#                   basin +
#                   s(LON) +
#                   s(LAT) +
#                   Land +
#                   logHWRF +
#                   (1|StormID)
#                   ),
#                   sigma ~ logHWRF
##      data = StormdataTrain4
# Fit12: VMAX ~   bf(VMAX ~ 
#                   basin +
#                   LON +
#                   LAT +
#                   Land +
#                   logMINSLP +
#                   logSTM_SPD +
#                   logHWFI +
#                   logVMAX_OP_T0 +
#                   logHWRF +
#                   (1|StormID)
#                   ),
#                   sigma ~ logHWRF
##      data = StormdataTrain4
# Fit13: VMAX ~   bf(VMAX ~ 
#                   basin +
#                   LON +
#                   LAT +
#                   LON:basin +
#                   LAT:basin +
#                   Land +
#                   logMINSLP +
#                   logSTM_SPD +
#                   logHWFI +
#                   logVMAX_OP_T0 +
#                   logHWRF +
#                   (1|StormID)
#                   ),
#                   sigma ~ logHWRF
##      data = StormdataTrain4
# Fit14: bf(VMAX ~ 
#             basin +
#             LON +
#             LAT +
#             Land +
#             arcsinhMINSLP +
#             arcsinhSHR_MAG +
#             arcsinhSTM_SPD +
#      arcsinhSST +
#      arcsinhRHLO +
#      arcsinhCAPE1 +
#      arcsinhCAPE3 +
#      arcsinhSHTFL2 +
#      arcsinhTCOND7002 +
#      arcsinhINST2 +
#      arcsinhCP1 +
#      arcsinhTCONDSYM2 +
#      arcsinhCOUPLSYM3 +
#      arcsinhHWFI +
#      arcsinhVMAX_OP_T0 +
#      arcsinhHWRF +
#      (1|StormID),
#    sigma ~ arcsinhHWRF
# ),
# data = StormdataTrain5,

## Model ----
iters <- 5000
burn <- 1000
chains <- 2
sims <- (iters-burn)*chains

formulaVMAX <- 
  bf(VMAX ~
       #basin +
       #LON + 
       #LAT +
       Land +
       #arcsinhMINSLP +
       arcsinhSHR_MAG +
       arcsinhSTM_SPD +
       #arcsinhSST +
       arcsinhRHLO +
       #arcsinhCAPE1 +
       arcsinhCAPE3 +
       arcsinhSHTFL2 +
       #arcsinhTCOND7002 +
       #arcsinhINST2 +
       #arcsinhCP1 +
       arcsinhTCONDSYM2 +
       arcsinhCOUPLSYM3 +
       arcsinhHWFI +
       #arcsinhVMAX_OP_T0 +
       arcsinhHWRF +
       (1|StormID),
     sigma ~ arcsinhHWRF
  ) + brmsfamily(family = "lognormal", link = "identity")

default_prior(formulaVMAX, data = StormdataTrain6)

priorsVMAX <- c(
  #prior(horseshoe(1), class = "b"),
  prior(normal(0,10), class = "b")
)

system.time(
  logNormalFit <- brm(
    formulaVMAX,
    data = StormdataTrain6,
    prior = priorsVMAX,
    save_pars = save_pars(all = TRUE), 
    chains = chains,
    iter = iters,
    cores = parallel::detectCores(),
    seed = 52,
    warmup = burn,
    normalize = FALSE,
    # prior = c(
    #   prior(normal(62, 30), class = b, coef = Intercept)
    # ),
    # knots = list(
    #   # Day = c(0, 365)),
    # Day = c(scale(0,
    #               center = attr(StormdataTrain4 |> pull(Day), "scaled:center"),
    #               scale = attr(StormdataTrain4 |> pull(Day), "scaled:scale")),
    #         scale(365,
    #               center = attr(StormdataTrain4 |> pull(Day), "scaled:center"),
    #               scale = attr(StormdataTrain4 |> pull(Day), "scaled:scale"))
    # )),
    control = list(adapt_delta = 0.95),
    backend = "cmdstanr"
  )
)

fit <- 37

### Save Fit ----
logNormalFit37 <- logNormalFit
fileLoc <- "~/Desktop/Temp Hurricane Model Data/logNormal Fits/" 
save(logNormalFit37,
     file = paste0(fileLoc, "logNormalFit", fit, ".RData"))

#load(file = "~/Desktop/Temp Hurricane Model Data/logNormal Fits/logNormalFit37.RData" )

### Diagnostics ----
prior_summary(logNormalFit)
posterior_summary(logNormalFit)
#launch_shinystan(logNormalFit)

print(logNormalFitNULL, digits = 4)
print(logNormalFitHWRF, digits = 4)
print(logNormalFitHWRFstormID, digits = 4)
print(logNormalFit12, digits = 4)
print(logNormalFit12B, digits = 4)
print(logNormalFit12C, digits = 4)
print(logNormalFit12D, digits = 4)
print(logNormalFit15, digits = 4)
print(logNormalFit16, digits = 4)
print(logNormalFit17, digits = 4)
print(logNormalFit32, digits = 4)

print(logNormalFit, digits = 4)
fixedEff <- fixef(logNormalFit)
fixedEff <- data.frame(fixedEff) |>
  mutate(
    p_val = dnorm(Estimate/Est.Error)
  ) |>
  mutate(
    across(everything(), function(x){round(x, 4)})
  ) |>
  mutate(
    Sig = ifelse(p_val < 0.01, "***",
                 ifelse(p_val < 0.05, "**",
                        ifelse(p_val < 0.1, "*", "")))
  )
#print(fixedEff, digits = 4)
fixedSigEff <- fixedEff |> filter(p_val < 0.2)
# fixedSigEff <- fixedSigEff |> 
#   rownames_to_column() |>
#   mutate(
#     response = str_split_i(rowname, "_", i = 1),
#     param = str_remove(rowname, paste0(response,"_"))
#   ) |> 
#   relocate(c(response, param), .after = "rowname") |>
#   select(-rowname)
print(fixedSigEff)



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
mean(abs(StormdataTrain$VMAX - StormdataTrain$HWRF))
model_performance(logNormalFit)

variance_decomposition(logNormalFit)
VarCorr(logNormalFit)
fixedEff <- fixef(logNormalFit)
fixedEff2 <- data.frame(fixedEff) |>
  mutate(
    p_val = dnorm(Estimate/Est.Error)
  ) |>
  mutate(
    across(everything(), function(x){round(x, 4)})
  ) |>
  mutate(
    Sig = ifelse(p_val < 0.01, "***",
                 ifelse(p_val < 0.05, "**",
                        ifelse(p_val < 0.1, "*", "")))
  )
fixedEff2
ranef(logNormalFit)



logNormalFitR2 <- bayes_R2(logNormalFit) |>
  bind_cols(Fit = paste0("logNormalFit", fit)) |>
  select(Fit, everything())
logNormalFitR2

bayes_factor(logNormalFit, logNormalFitHWRFstormID)
bayes_factor(logNormalFit, logNormalFit8)
bayes_factor(logNormalFit, logNormalFit15)

logNormalFitsmooths <- conditional_smooths(logNormalFit,
                                           method = "posterior_predict")
plot(logNormalFitsmooths, 
     stype = "raster", 
     ask = FALSE,
     theme = theme(legend.position = "bottom"))
plot(logNormalFitsmooths, 
     stype = "contour", 
     ask = FALSE,
     theme = theme(legend.position = "bottom"))

logNormalFiteffects <- conditional_effects(logNormalFit, 
                                           method = "posterior_predict",
                                           robust = FALSE)
plot(logNormalFiteffects, 
     points = TRUE, 
     ask = FALSE)

logNormalFiteffects2 <- conditional_effects(logNormalFit12B, 
                                            method = "posterior_predict",
                                            robust = FALSE,
                                            effects = "logHWRF")
plot(logNormalFiteffects2, 
     points = TRUE, 
     ask = FALSE)

logNormalFiteffects3 <- logNormalFiteffects2
logNormalFiteffects3$logHWRF <- logNormalFiteffects3$logHWRF |>
  mutate(
    logHWRF = exp(logHWRF),
    effect1__ = exp(effect1__)
  )
plot(logNormalFiteffects3, 
     points = TRUE, 
     ask = FALSE)

### Prediction ----
logNormalFit <- logNormalFit35
fit <- 35
## Fitted
logNormalFitfinalFit <- posterior_predict(logNormalFit)
# logNormalFitfinalResiduals <- t(StormdataTrain3$VMAX - t(logNormalFitfinalFit))
# logNormalFitfinalResidualsMean <- colMeans(logNormalFitfinalResiduals)
logNormalFitfinalFitMean <- colMeans(logNormalFitfinalFit)
logNormalFitfinalFitMed <- apply(logNormalFitfinalFit, 2, function(x){quantile(x, 0.5)})
logNormalFitfinalFitLCB <- apply(logNormalFitfinalFit, 2, function(x){quantile(x, 0.025)})
logNormalFitfinalFitUCB <- apply(logNormalFitfinalFit, 2, function(x){quantile(x, 0.975)})

## Prediction on new data
logNormalFitfinalPreds <- posterior_predict(logNormalFit, 
                                            newdata = StormdataTest6,
                                            allow_new_levels = TRUE, 
                                            re_formula = NULL)
logNormalFitfinalPredsMean <- colMeans(logNormalFitfinalPreds)
logNormalFitfinalPredsMed <- apply(logNormalFitfinalPreds, 2, function(x){quantile(x, 0.5, na.rm = TRUE)})
logNormalFitfinalPredsLCB <- apply(logNormalFitfinalPreds, 2, function(x){quantile(x, 0.025, na.rm = TRUE)})
logNormalFitfinalPredsUCB <- apply(logNormalFitfinalPreds, 2, function(x){quantile(x, 0.975, na.rm = TRUE)})

logNormalFitpredMetrics <- tibble(
  Fit = paste0("logNormalFit", fit),
  MAE_HWRF_fit = mean(abs(StormdataTrain$HWRF - StormdataTrain$VMAX)),
  MAE_fit = mean(abs(logNormalFitfinalFitMean - StormdataTrain$VMAX)),
  COV_fit = mean(logNormalFitfinalFitLCB < StormdataTrain$VMAX & StormdataTrain$VMAX < logNormalFitfinalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(logNormalFitfinalPredsMean - Actual_Yvec), na.rm = TRUE),
  MAD_pred = mean(abs(logNormalFitfinalPredsMed - Actual_Yvec), na.rm = TRUE),
  COV_pred = mean(logNormalFitfinalPredsLCB < Actual_Yvec & Actual_Yvec < logNormalFitfinalPredsUCB)
)
logNormalFitpredMetrics

### Plotting ----
#### Fit ----
ppc_dens_overlay(y = Actual_Yvec, 
                 yrep = logNormalFitfinalPreds[sample(1:sims, 200, replace = FALSE),]) +
  labs(title = "logNormalFit Predict") +
  theme_bw()
pp_check(logNormalFit,
         ndraws = 200,
         type = "dens_overlay") +
  labs(title = "logNormalFit Predict") +
  theme_bw()

logNormalFitFitDF <- bind_cols(
  StormdataTrain,
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

#### Residuals ----
logNormalFitResiduals <- residuals(logNormalFit, 
                                   method = "posterior_predict",
                                   summary = FALSE)

logNormalFitResidualsSum <- colMeans(logNormalFitResiduals)

# Extract residuals and fitted values from the baseline model
residuals_baseline <- residuals(logNormalFit, summary = TRUE)
fitted_vals <- fitted(fit_baseline)

# Plot residuals vs Sea_Surface_Temp to check for heteroscedasticity
plot(StormdataTrain5$arcsinhMINSLP, logNormalFitResidualsSum, 
     xlab = "arcsinhMINSLP", ylab = "Residuals", 
     main = "Residuals vs arcsinhMINSLP")

# Similarly, plot residuals against other predictors
plot(hurricane_data$Pressure, residuals_baseline, 
     xlab = "Pressure", ylab = "Residuals", 
     main = "Residuals vs Pressure")

modelParams <- row.names(fixef(logNormalFit15))
modelParams <- str_subset(modelParams, "sigma", negate = TRUE)
modelParams <- str_subset(modelParams, "Intercept", negate = TRUE)
modelParams[1] <- "basin"
modelParams[4] <- "Land"
modelParams
resids_list <- list()
for(i in modelParams){
  resids_list[[i]] <- ppc_error_scatter_avg_vs_x(StormdataTrain5$VMAX, 
                                                 logNormalFitfinalFit,
                                                 as.numeric(StormdataTrain5[[i]])) +
    geom_smooth(method = "lm", orientation = "y", level = 0.95) +
    labs(y = i)
}
resids_list

ppc_error_scatter_avg(StormdataTrain5$VMAX, 
                      logNormalFitfinalFit)

ppc_error_scatter_avg_vs_x(StormdataTrain5$VMAX, 
                           logNormalFitfinalFit,
                           as.numeric(StormdataTrain5[["arcsinhMINSLP"]])) +
  geom_smooth(orientation = "y", level = 0.95) +
  labs(y = "arcsinhMINSLP")



#### Prediction ----
logNormalFitPredDF <- bind_cols(
  StormdataTest,
  LCB = logNormalFitfinalPredsLCB,
  Mean = logNormalFitfinalPredsMean,
  Med = logNormalFitfinalPredsMed,
  UCB = logNormalFitfinalPredsUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) 

logNormalFitstormsPredplot <- ggplot(data = logNormalFitPredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX, color = "Observed")) +
  geom_line(aes(y = Mean, color = "PPD Mean")) +
  facet_wrap(vars(StormID))+#, ncol = 6)+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50)) +
  labs(title = "logNormalFit PPD Mean vs Observed VMAX",
       subtitle = "95% Credible Interval about PPD Mean") +
  scale_color_manual(name = NULL, values = c("black","red")) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1))
  ) +
  theme_bw()
logNormalFitstormsPredplot

stormPlot <- c(
  "3612015",
  
)

logNormalFitstormsPredplot2 <- 
  ggplot(data = logNormalFitPredDF |> 
           filter(StormID == "1002017"), 
         aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue") +
  geom_line(aes(y = VMAX, color = "Observed")) +
  geom_line(aes(y = Mean, color = "PPD Mean")) +
  #facet_wrap(vars(StormID))+#, ncol = 6)+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50)) +
  labs(title = "logNormalFit PPD Mean vs Observed VMAX",
       subtitle = "95% Credible Interval about PPD Mean") +
  scale_color_manual(name = NULL, values = c("black","red")) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1))
  ) +
  theme_bw()
logNormalFitstormsPredplot2

### PPC ----
###### Quantile 2.5 
logNormalFitLCBsims <- apply(logNormalFitfinalFit, 
                             MARGIN = 1,
                             function(x){
                               quantile(x, 0.025)
                             })
logNormalFitLCBpvalueVec <- logNormalFitLCBsims < quantile(StormdataTrain$VMAX, 0.025)
logNormalFitLCBpvalue <- sum(logNormalFitLCBpvalueVec)
logNormalFitLCBpvalue <- round(logNormalFitLCBpvalue/(sims), 3)
logNormalFitLCBpvalue <- min(logNormalFitLCBpvalue, 1 - logNormalFitLCBpvalue)

logNormalFit_ppcLCB <- 
  ppc_stat(StormdataTrain$VMAX,
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
logNormalFitUCBpvalueVec <- logNormalFitUCBsims < quantile(StormdataTrain$VMAX, 0.975)
logNormalFitUCBpvalue <- as.numeric(sum(logNormalFitUCBpvalueVec))
logNormalFitUCBpvalue <- round(logNormalFitUCBpvalue/sims, 3)
logNormalFitUCBpvalue <- min(logNormalFitUCBpvalue, 1 - logNormalFitUCBpvalue)

logNormalFit_ppcUCB <- 
  ppc_stat(StormdataTrain$VMAX,
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
logNormalFitMEANpvalueVec <- logNormalFitMEANsims < mean(StormdataTrain$VMAX)
logNormalFitMEANpvalue <- sum(logNormalFitMEANpvalueVec)
logNormalFitMEANpvalue <- round(logNormalFitMEANpvalue/sims, 3)
logNormalFitMEANpvalue <- min(logNormalFitMEANpvalue, 1 - logNormalFitMEANpvalue)

logNormalFit_ppcMEAN <- 
  ppc_stat(StormdataTrain$VMAX,
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
logNormalFitMEDpvalueVec <- logNormalFitMEDsims < quantile(StormdataTrain$VMAX, 0.5)
logNormalFitMEDpvalue <- sum(logNormalFitMEDpvalueVec)
logNormalFitMEDpvalue <- round(logNormalFitMEDpvalue/sims, 3)
logNormalFitMEDpvalue <- min(logNormalFitMEDpvalue, 1 - logNormalFitMEDpvalue)

logNormalFit_ppcMED <- 
  ppc_stat(StormdataTrain$VMAX,
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
logNormalFitSDpvalueVec <- logNormalFitSDsims < sd(StormdataTrain$VMAX)
logNormalFitSDpvalue <- sum(logNormalFitSDpvalueVec)
logNormalFitSDpvalue <- round(logNormalFitSDpvalue/sims, 3)
logNormalFitSDpvalue <- min(logNormalFitSDpvalue, 1 - logNormalFitSDpvalue)

logNormalFit_ppcSD <- 
  ppc_stat(StormdataTrain$VMAX,
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
logNormalFitRANGEpvalueVec <- logNormalFitRANGEsims < (max(StormdataTrain$VMAX)-min(StormdataTrain3$VMAX))
logNormalFitRANGEpvalue <- sum(logNormalFitRANGEpvalueVec)
logNormalFitRANGEpvalue <- round(logNormalFitRANGEpvalue/sims, 3)
logNormalFitRANGEpvalue <- min(logNormalFitRANGEpvalue, 1 - logNormalFitRANGEpvalue)

logNormalFit_ppcRANGE <- 
  ppc_stat(StormdataTrain$VMAX,
           logNormalFitfinalFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", logNormalFitRANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFit_ppcRANGE

### Combined Plot ----
logNormalFit_ppcComb <- 
  logNormalFitppcFit /
  (logNormalFit_ppcLCB | logNormalFit_ppcMED | logNormalFit_ppcUCB) /
  (logNormalFit_ppcRANGE | logNormalFit_ppcMEAN | logNormalFit_ppcSD)
logNormalFit_ppcComb

### Bayes p-values ----
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

### CV ----
set.seed(52)
kfoldID <- kfold_split_grouped(K = 5, StormdataTrain$StormID)
logNormalFitkfoldgroup <- kfold(logNormalFit,
                                folds = kfoldID,
                                chains = 1,
                                save_fits = TRUE)
save(logNormalFitkfoldgroup,
     file = "~/Desktop/Temp Hurricane Model Data/logNormalFit12kfold.RData")
logNormalFitkfoldPreds <- kfold_predict(logNormalFitkfoldgroup)
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
#pvalueComp_logNormalAll <- pvalueComp_logNormalBase
pvalueComp_logNormalAll <- bind_rows(
  pvalueComp_logNormalAll,
  logNormalFitpvalues
) |>
  distinct(Fit, .keep_all = TRUE)
pvalueComp_logNormalAll
save(pvalueComp_logNormalAll, 
     file = "~/Desktop/Temp Hurricane Model Data/pvalueComp_logNormalAll.RData")

## LOO ----
logNormalFitloo1 <- loo(logNormalFit1)
#attributes(logNormalFitloo1)$model_name <- "logNormalFit1"
logNormalFitloo2 <- loo(logNormalFit2)
logNormalFitloo3 <- loo(logNormalFit3)
logNormalFitloo4 <- loo(logNormalFit4)
NormalFitloo4 <- loo(NormalFit4)
logNormalFitloo5 <- loo(logNormalFit5)
logNormalFitloo6 <- loo(logNormalFit6)
logNormalFitloo7 <- loo(logNormalFit7)
logNormalFitloo8 <- loo(logNormalFit8)
logNormalFitloo9 <- loo(logNormalFit9)
logNormalFitloo10 <- loo(logNormalFit10)
logNormalFitloo11 <- loo(logNormalFit11)
logNormalFitloo12 <- loo(logNormalFit12)
logNormalFitloo12B <- loo(logNormalFit12B)
logNormalFitloo12C <- loo(logNormalFit12C)
logNormalFitloo12D <- loo(logNormalFit12D)
logNormalFitloo13 <- loo(logNormalFit13)
logNormalFitloo14 <- loo(logNormalFit14)
logNormalFitloo15 <- loo(logNormalFit15)
logNormalFitloo16 <- loo(logNormalFit16)
logNormalFitloo17 <- loo(logNormalFit17)
logNormalFitloo18 <- loo(logNormalFit18)
logNormalFitloo19 <- loo(logNormalFit19)
logNormalFitloo22 <- loo(logNormalFit22)
logNormalFitloo24 <- loo(logNormalFit24)
logNormalFitloo25 <- loo(logNormalFit25)
logNormalFitloo26 <- loo(logNormalFit26)
logNormalFitloo27 <- loo(logNormalFit27)
logNormalFitloo32 <- loo(logNormalFit32)
logNormalFitloo33 <- loo(logNormalFit33)
logNormalFitloo34 <- loo(logNormalFit34)
logNormalFitloo35 <- loo(logNormalFit35)
logNormalFitloo36 <- loo(logNormalFit36)
logNormalFitloo37 <- loo(logNormalFit37)

looComp_logNormalAll <- loo_compare(
  logNormalFitNULLloo,
  logNormalFitHWRFloo,
  logNormalFitHWRFstormIDloo,
  logNormalFitloo1,
  logNormalFitloo2,
  logNormalFitloo3,
  logNormalFitloo4,
  NormalFitloo4,
  logNormalFitloo5,
  logNormalFitloo6,
  logNormalFitloo7,
  logNormalFitloo8,
  logNormalFitloo9,
  logNormalFitloo10,
  logNormalFitloo11,
  logNormalFitloo12,
  logNormalFitloo12B,
  logNormalFitloo12C,
  logNormalFitloo12D,
  logNormalFitloo13,
  logNormalFitloo14,
  logNormalFitloo15,
  logNormalFitloo16,
  logNormalFitloo17,
  logNormalFitloo18,
  logNormalFitloo19,
  logNormalFitloo22,
  logNormalFitloo24,
  logNormalFitloo25,
  logNormalFitloo26,
  logNormalFitloo27,
  logNormalFitloo30,
  logNormalFitloo31,
  logNormalFitloo32,
  logNormalFitloo33,
  logNormalFitloo34,
  logNormalFitloo35,
  logNormalFitloo37
)
looComp_logNormalAll

loo_compare(
  logNormalFitloo33,
  logNormalFitloo34,
  logNormalFitloo35,
  logNormalFitloo36,
  logNormalFitloo37
)

save(looComp_logNormalAll, 
     file = "~/Desktop/Temp Hurricane Model Data/looComp_logNormalAll.RData")

## CV ----
#cvComp_logNormalAll <- cvComp_logNormalBase
cvComp_logNormalAll <- bind_rows(
  cvComp_logNormalAll,
  logNormalFitkfoldMetrics
) |>
  distinct(Fit, .keep_all = TRUE) |>
  arrange(MAE_kfold)

cvComp_logNormalAll
save(cvComp_logNormalAll, 
     file = "~/Desktop/Temp Hurricane Model Data/cvComp_logNormalAll.RData")

## Preds ----
#predComp_logNormalAll <- predComp_logNormalBase
predComp_logNormalAll <- bind_rows(
  predComp_logNormalAll,
  logNormalFitpredMetrics
) |>
  distinct(Fit, .keep_all = TRUE) |>
  arrange(MAE_pred)

predComp_logNormalAll
save(predComp_logNormalAll, 
     file = "~/Desktop/Temp Hurricane Model Data/predComp_logNormalAll.RData")


## Bayes R2 ----
#bayesR2_logNormalAll <- bayesR2_logNormalBase
bayesR2_logNormalAll <- bind_rows(
  bayesR2_logNormalAll,
  logNormalFitR2
) |> 
  distinct(Fit, .keep_all = TRUE) |>
  arrange(desc(Estimate))

bayesR2_logNormalAll
save(bayesR2_logNormalAll, 
     file = "~/Desktop/Temp Hurricane Model Data/bayesR2_logNormalAll.RData")


pvalueComp_logNormalAll
looComp_logNormalAll
cvComp_logNormalAll
predComp_logNormalAll
bayesR2_logNormalAll

## Print Models ----
logNormalFitModels <- list(
  logNormalFitNULL$formula,
  logNormalFitHWRF$formula,
  logNormalFitHWRFstormID$formula,
  logNormalFit1$formula,
  logNormalFit2$formula,
  logNormalFit3$formula,
  logNormalFit4$formula,
  logNormalFit5$formula,
  logNormalFit6$formula,
  logNormalFit7$formula,
  logNormalFit8$formula,
  logNormalFit9$formula,
  logNormalFit10$formula,
  logNormalFit11$formula,
  logNormalFit12$formula,
  logNormalFit13$formula,
  logNormalFit14$formula,
  logNormalFit15$formula,
  logNormalFit16$formula,
  logNormalFit17$formula,
  logNormalFit18$formula,
  logNormalFit19$formula,
  logNormalFit20$formula,
  logNormalFit21$formula,
  logNormalFit22$formula,
  logNormalFit23$formula,
  logNormalFit24$formula,
  logNormalFit25$formula,
  logNormalFit26$formula,
  logNormalFit27$formula,
  logNormalFit30$formula,
  logNormalFit31$formula)

save(logNormalFitModels,
     file = "~/Desktop/Temp Hurricane Model Data/logNormalModelFits.RData")

print(logNormalFitNULL, digits = 4)
print(logNormalFitHWRF, digits = 4)
print(logNormalFitHWRFstormID, digits = 4)
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
print(logNormalFit11, digits = 4)
print(logNormalFit12, digits = 4)
print(logNormalFit13, digits = 4)
print(logNormalFit31, digits = 4)
