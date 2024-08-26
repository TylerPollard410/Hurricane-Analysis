# Load Libraries ----
library(knitr)
library(data.table)
library(MASS)
library(rjags)
library(plyr)
library(stringr)
library(tidyverse)
library(tictoc)
library(cowplot)
library(caret)
library(splines)
library(DescTools)
library(bayesplot)
library(BayesFactor)
library(brms)
library(rstanarm)
library(tidybayes)
library(lubridate)
library(tidyverse)

# Read in data ----
## Clean data ----
Stormdata_raw <- fread("~/Desktop/Hurricane Analysis/_data/E2_data.csv")
Stormdata <- Stormdata_raw |>
  mutate(
    StormID = factor(StormID),
    basin = factor(basin),
    Date = as_datetime(Date, tz = "UTC")
    #Date = as_datetime(format(Date, "%Y-%m-%d %H:%M:%S"),  tz = "UTC")
  )
#as.POSIXct(Stormdata$Date[1])
# str(Stormdata)

## Create training and test data sets ----
StormdataTrain <- Stormdata |> filter(complete.cases(VMAX))
StormdataTest <- Stormdata |> filter(!complete.cases(VMAX))

## Transform data ----
# Remove not varying 
StormdataTrain2 <- StormdataTrain |>
  select(-lead_time)

# Create Date vars 
dataYears <- year(StormdataTrain2$Date)
dataMonths <- month(StormdataTrain2$Date, label = TRUE)
dataDays <- day(StormdataTrain2$Date)

StormdataTrain3 <- StormdataTrain2 |>
  mutate(
    Year = factor(dataYears, ordered = TRUE),
    Month = dataMonths
  ) |>
  group_by(StormID) |>
  mutate(
    StormElapsedTime = as.numeric(difftime(Date, min(Date), units = "hours"))
  ) |>
  select(
    StormID,
    Date,
    Year,
    Month,
    StormElapsedTime,
    everything()
  ) |>
  ungroup()

summary(StormdataTrain3$VMAX)

## Plot VMAX ----
### Histogram ----
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

### Time ----
ggplot(data = StormdataTrain3) +
  geom_point(aes(x = Date, y = VMAX)) +
  scale_x_datetime(date_breaks = "month", 
                   date_minor_breaks = "day", 
                   date_labels = "%b-%Y") + 
  theme(
    axis.text.x = element_text(angle = 90)
  )

ggplot(data = StormdataTrain3) +
  geom_point(aes(y = StormID, x = VMAX, color = StormID))

### Map ----
# Do by time next
world_coordinates <- map_data("world", ) 
ggplot() + 
  # geom_map() function takes world coordinates  
  # as input to plot world map 
  geom_map( 
    data = world_coordinates, map = world_coordinates, 
    aes(x = long, y = lat, map_id = region) 
  ) + 
  geom_point(
    data = StormdataTrain3,
    aes(x = LON-360, y = LAT, 
        color = VMAX)
  ) +
  xlim(c(-180,0)) +
  ylim(c(0,60)) +
  scale_color_continuous(low = "green", high = "red") +
  theme_bw()

### Plots by Factor ----
num_cols <- colnames(StormdataTrain3)[sapply(StormdataTrain3, function(x){is.numeric(x)})]
num_cols

# Fit model ----
generalBayes1 <- generalTestBF(
  formula = VMAX ~ .,
  data = modelData,
  whichRandom = "StormID",
  multicore = TRUE,
  iterations = 500
)

# Fit prelim models
Y <- StormdataTrain3 |> select(VMAX)
X <- StormdataTrain3 |> 
  select(
    #"StormID",
    #Date,
    #Year,
    #Month,
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
    "HWRF"
  )
modelData <- StormdataTrain3 |>
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
str(modelData)

L <- 10 # number of knots
B1   <- bs(modelData$LON, df=2*L, intercept=TRUE) # Longitude basis functions
B2   <- bs(modelData$LON, df=L, intercept=TRUE)   # Latitude basis functions
X    <- NULL
for(j in 1:ncol(B1)){
  for(k in 1:ncol(B2)){
    X <- cbind(X,B1[,j]*B2[,k])  # Products
  }
}
X    <- X[,apply(X,2,max)>0.1]  # Remove basis function that are near zero for all sites
X    <- ifelse(X>0.001,X,0)
p    <- ncol(X)

fields::BN

index <- sample(1:p, 1)
spline_data <- modelData |>
  select(
    LON, 
    LAT
  ) |>
  add_column(spline = X[spline])
#world_coordinates <- map_data("world") 
ggplot() + 
  # geom_map() function takes world coordinates  
  # as input to plot world map 
  geom_map( 
    data = world_coordinates, map = world_coordinates, 
    aes(x = long, y = lat, map_id = region) 
  ) + 
  geom_point(
    data = spline_data,
    aes(x = LON-360, y = LAT, 
        color = spline)
  ) +
  xlim(c(-180,0)) +
  ylim(c(0,60)) +
  scale_color_continuous(low = "green", high = "red") +
  theme_bw()

















## Create training and test data sets ----
training_data <- Stormdata |> filter(complete.cases(VMAX))
test_data <- Stormdata |> filter(!complete.cases(VMAX))


# 1. Model description ----
## Load data ----
Y <- training_data$VMAX
Z <- log(Y)
X <- training_data |>
  select(-c(
    StormID,
    Date,
    basin,
    lead_time,
    LAT,
    LON,
    MINSLP,
    SHR_MAG,
    STM_SPD,
    SST,
    RHLO,
    CAPE1,
    CAPE3,
    SHTFL2,
    TCOND7002,
    INST2,
    CP1,
    TCONDSYM2,
    HWFI,
    VMAX_OP_T0,
    VMAX
  ))
Xscale <- scale(X)
Xscale <- cbind("Intercept" = rep(1,length(Y)), Xscale)
HWRF <- training_data$HWRF

### Create StormID indicator for random effects model ----
StormIDs <- training_data |> select(StormID)
Storms <- unique(StormIDs$StormID)
IDs <- 1:length(Storms)
Storms <- data.frame(
  StormID = Storms,
  ID = IDs
)
StormIDs <- left_join(StormIDs, Storms)
StormID <- StormIDs$ID

n <- length(Y)
p <- ncol(Xscale)
N <- length(IDs)

burn <- 10000
iters <- 20000
chains <- 2
thin <- 5

# 2. Old Analysis ----
outputDF <- fread(file = "PollardTyler.csv")
colnames(outputDF) <- c("HWRF", "VMAXtlyer", "Ltyler", "Utyler")

ActualY <- read.csv("Actual Y.csv")

TotalY <- cbind(outputDF, ActualY)
TotalYpreds <- TotalY |> filter(!is.na(x))

MAD1 <- mean(abs(TotalYpreds$Y - TotalYpreds$x))
COV1 <- mean( (TotalYpreds$L <= TotalYpreds$x) & (TotalYpreds$x <= TotalYpreds$U))

MADtyler <- mean(abs(TotalYpreds$VMAXtlyer - TotalYpreds$x))
COVtyler <- mean( (TotalYpreds$Ltyler <= TotalYpreds$x) & (TotalYpreds$x <= TotalYpreds$Utyler))

MAD1
MADtyler

COV1
COVtyler


# 3. New Analysis ----
Y <- Stormdata$VMAX[complete.cases(Stormdata$VMAX)]
newDF <- cbind(Xscale, Y) |> as.data.frame()
fit1 <- generalTestBF(Y ~ .,
                      data = newDF)
fit2 <- stan_glm(Y ~ .,
                 data = newDF)
fit2
summary(fit2, digits = 3, probs = c(0.025, 0.975))
car::vif(fit2)

pp_check(fit2, nreps = 1000)
ppc_bars(Y, posterior_predict(fit2))

fit2 <- stan_glm(log(Y) ~ .,
                 data = newDF,
                 family = gaussian())
fit2
summary(fit2, digits = 3, probs = c(0.025, 0.975))
car::vif(fit2)

