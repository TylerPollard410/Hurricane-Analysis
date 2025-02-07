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
library(bayesplot)
library(performance)
library(tidyverse)

# Read in data ----
Stormdata_raw <- fread("_data/E2_data.csv")
Actual_Y <- fread("_data/Actual Y.csv")
Actual_Yvec <- Actual_Y |> filter(complete.cases(x)) |> pull(x)

# Data ----
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
  ) |>
  mutate(
    across(c(StormID, Year, Month, basin),
           ~factor(.x)
    )
  )

### Summary of data ----
dataSum <- Stormdata |>
  group_by(DataType) |>
  summarise(
    across(where(is.numeric), 
           list(
             Min = ~min(.x, na.rm = TRUE),
             QLow = ~quantile(.x, 0.025, na.rm = TRUE),
             Median = ~quantile(.x, 0.5, na.rm = TRUE),
             Mean = ~mean(.x, na.rm = TRUE),
             QHigh = ~quantile(.x, 0.975, na.rm = TRUE),
             Max = ~max(.x, na.rm = TRUE)
           )
    )
  ) |>
  pivot_longer(
    -DataType,
    names_to = "name"
  ) |>
  mutate(
    Predictor = str_split_i(name, "_(?!.*_)", 1),
    Stat = str_split_i(name, "_(?!.*_)", 2),
    .after = name
  ) |>
  select(-name) |>
  pivot_wider(
    names_from = Predictor
  )

## Split data -----
StormdataTrain <- Stormdata |>
  filter(DataType == "Train") |>
  select(-c(
    DataType,
    Date,
    Day,
    LON2,
    StormElapsedTime2,
    Obs
  ))|>
  mutate(
    across(c(StormID, Year, Month, basin),
           ~droplevels(.x)
    )
  )

StormdataTest <- Stormdata |>
  filter(DataType == "Test") |>
  select(-c(
    DataType,
    Date,
    Day,
    LON2,
    StormElapsedTime2,
    Obs
  ))|>
  mutate(
    across(c(StormID, Year, Month, basin),
           ~droplevels(.x)
    )
  ) |>
  mutate(
    VMAX = Actual_Yvec
  )

## PreProcess data ----
corData <- cor(StormdataTrain |> select(where(is.numeric)))

### Regular ----
preProc <- preProcess(StormdataTrain |> 
                           select(
                             where(is.numeric),
                             -VMAX
                             #-StormElapsedTime,
                             #-LAT,
                             #-LON
                           ),
                         method = c("scale", "center"))
preProc
StormdataTrain2 <- predict(preProc, StormdataTrain)
StormdataTes2 <- predict(preProc, StormdataTest)
colnames(StormdataTrain)

### YeoJohnson ----
preProcYeo <- preProcess(StormdataTrain |> 
                           select(
                             where(is.numeric),
                             -VMAX,
                             #-HWRF,
                             -StormElapsedTime
                             #-LAT,
                             #-LON
                           ),
                         method = c("scale", "center", "YeoJohnson"))
preProcYeo
StormdataTrainYeo <- predict(preProcYeo, StormdataTrain)
StormdataTestYeo <- predict(preProcYeo, StormdataTest)
colnames(StormdataTrainYeo)

preProcYeoCorr <- preProcess(StormdataTrain |> 
                               select(
                                 where(is.numeric),
                                 -VMAX,
                                 -StormElapsedTime,
                                 -LAT,
                                 -LON
                               ),
                             method = c("scale", "center", "corr", "YeoJohnson"))
preProcYeoCorr
StormdataTrainYeoB <- predict(preProcYeoCorr, StormdataTrain)
StormdataTestYeoB <- predict(preProcYeoCorr, StormdataTest)
colnames(StormdataTrainYeoB)


### Arcsinh ----
StormdataTrainArcsinhPre <- StormdataTrain |>
  mutate(
    across(c(where(is.numeric), 
             -VMAX, #-HWRF,
             -StormElapsedTime),# -LAT, -LON),
           ~predict(arcsinh_x(.x, standardize = FALSE), newdata = .x)
    ))

StormdataTestArcsinhPre <- StormdataTest |>
  mutate(
    across(c(where(is.numeric), 
             -VMAX, #-HWRF,
             -StormElapsedTime),# -LAT, -LON),
           ~predict(arcsinh_x(.x, standardize = FALSE), newdata = .x)
    ))

preProcArcsinh <- preProcess(StormdataTrainArcsinhPre |> 
                               select(
                                 where(is.numeric),
                                 -VMAX, #-HWRF,
                                 -StormElapsedTime
                                 #-LAT,
                                 #-LON
                               ),
                             method = c("scale", "center"))
preProcArcsinh
StormdataTrainArcsinh <- predict(preProcArcsinh, StormdataTrainArcsinhPre)
StormdataTestArcsinh <- predict(preProcArcsinh, StormdataTestArcsinhPre)


preProcArcsinhCorr <- preProcess(StormdataTrainArcsinhPre |> 
                                   select(
                                     where(is.numeric),
                                     -VMAX,
                                     -StormElapsedTime,
                                     -LAT,
                                     -LON
                                   ),
                                 method = c("scale", "center", "corr"))
preProcArcsinhCorr
StormdataTrainArcsinhB <- predict(preProcArcsinhCorr, StormdataTrainArcsinhPre)
StormdataTestArcsinhB <- predict(preProcArcsinhCorr, StormdataTestArcsinhPre)


# Model ##################################################
# bf(VMAX ~ 
#      StormElapsedTime + 
#      LAT + 
#      LON + 
#      basin + 
#      Land +
#      MINSLP +
#      SHR_MAG +
#      STM_SPD +
#      SST +
#      RHLO +
#      CAPE1 +
#      CAPE3 +
#      SHTFL2 +
#      TCOND7002 +
#      INST2 +
#      CP1 +
#      TCONDSYM2 +
#      COUPLSYM3 +
#      HWFI +
#      VMAX_OP_T0 +
#      HWRF
# )

formulaVMAX <- 
  bf(VMAX ~ 
       #StormElapsedTime +
       LAT + 
       LON +  # Baseline spatial effects
       basin +
       Land +  # Categorical spatial context
       MINSLP + 
       SHR_MAG +
       STM_SPD + 
       SST + 
       RHLO + 
       CAPE1 + 
       CAPE3 + 
       SHTFL2 + 
       TCOND7002 + INST2 + 
       CP1 + 
       TCONDSYM2 +
       COUPLSYM3 +  # Physical storm predictors
       HWFI +
       VMAX_OP_T0 +  # Operational estimates
       HWRF +  # Benchmark
       (1 | StormID)#,  # Random effect for storm-specific variation
     #sigma ~ HWRF #+ HWFI + STM_SPD
     #nl = TRUE
  ) + brmsfamily(family = "Gamma", link = "log")

default_prior(formulaVMAX, data = StormdataTrainYeo)

priorsVMAX <- c(
  #prior(horseshoe(1), class = "b")
  prior(normal(0, 5), class = "b"),
  #prior(inv_gamma(0.1, 0.1), class = "sigma"),
  prior(inv_gamma(0.1, 0.1), class = "shape"),
  prior(inv_gamma(0.1, 0.1), class = "sd")
)

## Fit brms ----
iters <- 3000
burn <- 1000
chains <- 2
sims <- (iters-burn)*chains

system.time(
  logNormalFit <- brm(
    formulaVMAX,
    #data = StormdataTrainArcsinh,
    data = StormdataTrainYeo,
    #data = StormdataTrain2,
    prior = priorsVMAX,
    save_pars = save_pars(all = TRUE), 
    chains = chains,
    iter = iters,
    cores = parallel::detectCores(),
    seed = 52,
    warmup = burn,
    #init = 0,
    normalize = FALSE,
    control = list(adapt_delta = 0.95),
    backend = "cmdstanr"
  )
)

## Diagnostics ----
fit <- 4
assign(paste0("Fit", fit), logNormalFit)
#logNormalFit <- Fit1

plot(logNormalFit, ask = FALSE)
#prior_summary(logNormalFit)

print(logNormalFit, digits = 4)

waicList <- list()
waic <- waic(logNormalFit)
attributes(waic)$model_name <- paste0("logNormalFit", fit)
waicList[[paste0("fit", 4)]] <- waic4

### Compare Candidate Models ----
waicList

loo_compare(
  waicList 
)

save(waicList, file = "_data/waicList.RData")

### Multicollinearity ----
check_collinearity(logNormalFit)

### Check heteroskedasity ----
check_heteroscedasticity(logNormalFit)

### Fixed Effects ----
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
print(fixedEff, digits = 4)
fixedSigEff <- fixedEff |> filter(p_val < 0.2)
print(fixedSigEff)

pp_check(logNormalFit, ndraws = 100)

### Hypothesis Tests -----
posterior_summary(logNormalFit)

xVars <- str_subset(variables(logNormalFit), pattern = "b_")
hypothesis(logNormalFit, paste(xVars, "= 0"), 
           class = NULL, 
           alpha = 0.1)

hypothesis(logNormalFit, "sHWRF_1 = 0", class = "bs")
hypID <- hypothesis(logNormalFit, 
           "Intercept = 0", 
           group = "StormID", 
           scope = "coef")
plot(hypID)

#variance_decomposition(logNormalFit)
#VarCorr(logNormalFit)

### Residuals ----
fitResiduals <- 
  residuals(
    logNormalFit,
    method = "posterior_predict",
    re_formula = NULL,
    robust = FALSE,
    probs = c(0.025, 0.975)) |>
  data.frame()
mean(abs(fitResiduals$Estimate))

predResiduals <- 
  residuals(
    logNormalFit, 
    newdata = StormdataTestArcsinh,
    #newdata = StormdataTestYeo,
    #newdata = StormdataTest2,
    method = "posterior_predict",
    allow_new_levels = TRUE,
    re_formula = NULL,
    robust = FALSE,
    probs = c(0.025, 0.975)) |>
  data.frame()
mean(abs(predResiduals$Estimate))

## Posteriors ----
#logNormalFit <- Fit8
### Training ----
logNormalFitfinalFit <- posterior_predict(logNormalFit)
# logNormalFitfinalResiduals <- t(StormdataTrain3$VMAX - t(logNormalFitfinalFit))
# logNormalFitfinalResidualsMean <- colMeans(logNormalFitfinalResiduals)
logNormalFitfinalFitMean <- colMeans(logNormalFitfinalFit)
logNormalFitfinalFitMed <- apply(logNormalFitfinalFit, 2, function(x){quantile(x, 0.5)})
logNormalFitfinalFitLCB <- apply(logNormalFitfinalFit, 2, function(x){quantile(x, 0.025)})
logNormalFitfinalFitUCB <- apply(logNormalFitfinalFit, 2, function(x){quantile(x, 0.975)})

### Prediction ----
logNormalFitfinalPreds <- posterior_predict(logNormalFit, 
                                            newdata = StormdataTestArcsinh,
                                            #newdata = StormdataTestYeo,
                                            allow_new_levels = TRUE, 
                                            re_formula = NULL)
logNormalFitfinalPredsMean <- colMeans(logNormalFitfinalPreds)
logNormalFitfinalPredsMed <- apply(logNormalFitfinalPreds, 2, function(x){quantile(x, 0.5, na.rm = TRUE)})
logNormalFitfinalPredsLCB <- apply(logNormalFitfinalPreds, 2, function(x){quantile(x, 0.025, na.rm = TRUE)})
logNormalFitfinalPredsUCB <- apply(logNormalFitfinalPreds, 2, function(x){quantile(x, 0.975, na.rm = TRUE)})


## Plot Fit ----
### PPC ----
trainVMAX <- StormdataTrain$VMAX
testVMAX <- Actual_Yvec

ppc_dens_overlay(y = trainVMAX,
                 yrep = logNormalFitfinalFit[sample(1:sims, 100, replace = FALSE), ]) +
  labs(title = paste0("Fit", fit, " PPC"))

ppc_dens_overlay(y = testVMAX,
                 yrep = logNormalFitfinalPreds[sample(1:sims, 100, replace = FALSE), ]) +
  labs(title = paste0("Fit", fit, " PPD"))

### Residuals ----
xVarsCoef <- str_subset(variables(logNormalFit), pattern = "b_")
xVars <- str_split_i(xVarsCoef, "^[^_]+_", i = 2)
xVars <- ifelse(str_detect(xVars, "basin"), "basin", 
                ifelse(str_detect(xVars, "Land"), "Land", xVars))
xVars <- xVars[-1]

ppcErrorPlotXlist <- list()
for(i in xVars){
  if(is.numeric(StormdataTrainArcsinh[[i]])){
    ppcErrorPlotX <- ppc_error_scatter_avg_vs_x(
      y = trainVMAX,
      yrep = logNormalFitfinalFit[sample(1:sims, 100, replace = FALSE), ],
      x = StormdataTrainArcsinh |> pull(contains(i))
    ) +
      labs(title = paste0("Avg Error ", i))
    ppcErrorPlotXlist[[i]] <- ppcErrorPlotX
  }
}

# for(i in xVars){
#   if(is.numeric(StormdataTrainYeo[[i]])){
#     ppcErrorPlotX <- ppc_error_scatter_avg_vs_x(
#       y = trainVMAX,
#       yrep = logNormalFitfinalFit[sample(1:sims, 100, replace = FALSE), ],
#       x = StormdataTrainYeo[[i]]
#     ) +
#       labs(title = paste0("Avg Error ", i))
#     ppcErrorPlotXlist[[i]] <- ppcErrorPlotX
#   }
# }

wrap_plots(ppcErrorPlotXlist)

### Smooth Effects ----
#logNormalFitsmooths <- conditional_smooths(logNormalFit,
                                          # method = "posterior_predict")
# plot(logNormalFitsmooths, 
#      stype = "raster", 
#      ask = FALSE,
#      theme = theme(legend.position = "bottom"))
# plot(logNormalFitsmooths,
#      stype = "contour",
#      ask = FALSE,
#      theme = theme(legend.position = "bottom"))

### Fixed Effects ----
logNormalFiteffects <- conditional_effects(logNormalFit, 
                                           method = "posterior_predict",
                                           robust = FALSE)
plot(logNormalFiteffects, 
     points = TRUE, 
     ask = FALSE)

bayesR2 <- bayes_R2(logNormalFit) |>
  bind_cols(Fit = paste0("logNormalFit", fit)) |>
  select(Fit, everything())

bayesR2Comb <- bind_rows(
  bayesR2,
  bayesR2Comb
)
bayesR2Comb #<- bayesR2



looFit <- loo(logNormalFit)
waicFit <- waic(logNormalFit)



## Performance Metrics ----
logNormalFitpredMetrics <- tibble(
  Fit = paste0("Fit", fit),
  MAE_HWRF_fit = mean(abs(StormdataTrain$HWRF - StormdataTrain$VMAX)),
  MAE_fit = mean(abs(logNormalFitfinalFitMean - StormdataTrain$VMAX)),
  COV_fit = mean(logNormalFitfinalFitLCB < StormdataTrain$VMAX & StormdataTrain$VMAX < logNormalFitfinalFitUCB),
  MAE_HWRF_pred = mean(abs(StormdataTest$HWRF - Actual_Yvec)),
  MAE_pred = mean(abs(logNormalFitfinalPredsMean - Actual_Yvec), na.rm = TRUE),
  MAD_pred = mean(abs(logNormalFitfinalPredsMed - Actual_Yvec), na.rm = TRUE),
  COV_pred = mean(logNormalFitfinalPredsLCB < Actual_Yvec & Actual_Yvec < logNormalFitfinalPredsUCB)
)
logNormalFitpredMetrics

predMetricsComb <- bind_rows(
  predMetricsComb,
  logNormalFitpredMetrics
)
predMetricsComb #<- logNormalFitpredMetrics

### Loo ----
loo6 <- loo(Fit6)
loo14 <- loo(Fit14)
loo15 <- loo(Fit15)
loo16 <- loo(Fit16)
loo18 <- loo(Fit18)
loo19 <- loo(Fit19)

looComp <- loo_compare(
  loo6,
  loo14,
  loo15,
  loo16,
  loo18,
  loo19
)
looComp

### CV ----
set.seed(52)
kfoldID <- kfold_split_grouped(K = 5, StormdataTrain$StormID)
logNormalFitkfoldgroup <- kfold(Fit14,
                                folds = kfoldID,
                                chains = 1,
                                save_fits = TRUE)
save(logNormalFitkfoldgroup,
     file = "~/Desktop/Temp Hurricane Model Data/logNormalFit14kfold.RData")
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

## State-Space Plots ----  
### Training ----
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


### Prediction ----
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
