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
Stormdata <- fread("_data/E2_data.csv")
Stormdata <- Stormdata |>
  mutate(
    StormID = factor(StormID),
    basin = factor(basin),
    Date = as_datetime(Date)
  )
str(Stormdata)

## Create training and test data sets ----
StormdataTrain <- Stormdata |> filter(complete.cases(VMAX))
StormdataTest <- Stormdata |> filter(!complete.cases(VMAX))

lapply(StormdataTrain, function(x){length(unique(x))})

# Remove not varying 
StormdataTrain2 <- StormdataTrain |>
  select(-lead_time)

summary(StormdataTrain2$VMAX)

## Plot VMAX ----
### Histogram ----
ggplot(data = StormdataTrain2) +
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
ggplot(data = StormdataTrain2) +
  geom_point(aes(x = Date, y = VMAX)) +
  scale_x_datetime(date_breaks = "month", 
                   date_minor_breaks = "day", 
                   date_labels = "%m-%Y") + 
  theme(
    axis.text.x = element_text(angle = 90)
  )

ggplot(data = StormdataTrain2) +
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
    data = StormdataTrain2,
    aes(x = LON-360, y = LAT, 
        color = VMAX)
  ) +
  xlim(c(-180,0)) +
  ylim(c(0,60)) +
  scale_color_continuous(low = "green", high = "red") +
  theme_bw()

### Plots by Factor ----
num_cols <- colnames(StormdataTrain2)[sapply(StormdataTrain2, function(x){is.numeric(x)})]
num_cols


# Fit prelim models




parse_fact <- c(
  "YearSemester",
  "School",
  "Grade",
  "Gender2",
  "Race2"
)

fact_plots <- list()
for(i in 1:length(parse_fact)){
  fact <- parse_fact[i]
  fact_plot <- ggplot(data = Sciencesurvey_data) +
    geom_histogram(aes(x = ScienceScore, after_stat(density), fill = !!sym(fact)),
                   binwidth = 0.25,
                   position = position_dodge()) +
    geom_density(aes(x = ScienceScore)) +
    facet_wrap(vars(!!sym(fact))) +
    theme_bw() +
    theme(
      legend.position = "bottom"
    )
  fact_plots[[i]] <- fact_plot
}

cowplot::plot_grid(plotlist = fact_plots, ncol = 2)























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

