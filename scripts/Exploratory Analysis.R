# Load Libraries ----
library(knitr)
library(data.table)
library(MASS)
library(bestNormalize)
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
#library(rstanarm)
library(tidybayes)
library(loo)
library(brms)
library(bayesplot)
library(performance)
library(gt)
library(gtsummary)
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
    LON180 = LON - 360,
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
    LON180,
    everything(),
    -lead_time
  ) |>
  mutate(
    LAT2 = LAT,
    .after = LAT
  ) |>
  mutate(
    LON2 = LON,
    .after = LON
  )

### Land Indicator ----
pts <- st_as_sf(Stormdata1, # |> select(LON, LAT), 
                coords = c("LON180", "LAT"),
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
    LAT2,
    LON,
    LON2,
    LON180,
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
    LON180,
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
    LON180,
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
#preProc
StormdataTrain2 <- predict(preProc, StormdataTrain)
StormdataTes2 <- predict(preProc, StormdataTest)
#colnames(StormdataTrain)

### YeoJohnson ----
preProcYeo <- preProcess(StormdataTrain |> 
                           select(
                             where(is.numeric),
                             -VMAX,
                             #-HWRF,
                             -StormElapsedTime,
                             -LAT2,
                             -LON2
                           ),
                         method = c("scale", "center", "YeoJohnson"))
#preProcYeo
StormdataTrainYeo <- predict(preProcYeo, StormdataTrain)
StormdataTestYeo <- predict(preProcYeo, StormdataTest)
#colnames(StormdataTrainYeo)

preProcYeoCorr <- preProcess(StormdataTrain |> 
                               select(
                                 where(is.numeric),
                                 -VMAX,
                                 -StormElapsedTime,
                                 -LAT,
                                 -LON
                               ),
                             method = c("scale", "center", "corr", "YeoJohnson"))
#preProcYeoCorr
StormdataTrainYeoB <- predict(preProcYeoCorr, StormdataTrain)
StormdataTestYeoB <- predict(preProcYeoCorr, StormdataTest)
#colnames(StormdataTrainYeoB)


### Arcsinh ----
StormdataTrainArcsinhPre <- StormdataTrain |>
  mutate(
    across(c(where(is.numeric), 
             -VMAX, #-HWRF,
             -StormElapsedTime, -LAT2, -LON2),
           ~predict(arcsinh_x(.x, standardize = FALSE), newdata = .x)
    ))

StormdataTestArcsinhPre <- StormdataTest |>
  mutate(
    across(c(where(is.numeric), 
             -VMAX, #-HWRF,
             -StormElapsedTime, -LAT2, -LON2),
           ~predict(arcsinh_x(.x, standardize = FALSE), newdata = .x)
    ))

preProcArcsinh <- preProcess(StormdataTrainArcsinhPre |> 
                               select(
                                 where(is.numeric),
                                 -VMAX, #-HWRF,
                                 -StormElapsedTime,
                                 -LAT2,
                                 -LON2
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
                                     -LAT2,
                                     -LON2
                                   ),
                                 method = c("scale", "center", "corr"))
#preProcArcsinhCorr
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
       #LAT + 
       #LON + 
       #basin +
       Land +  
       #MINSLP + 
       SHR_MAG +
       #STM_SPD + 
       #SST + 
       RHLO + 
       #CAPE1 + 
       CAPE3 + 
       #SHTFL2 + 
       #TCOND7002 + 
       #INST2 + 
       CP1 + 
       TCONDSYM2 +
       COUPLSYM3 + 
       HWFI +
       #VMAX_OP_T0 +  
       HWRF +  
       (1 | StormID)  # Random effect for storm-specific variation
     #sigma ~ HWRF #+ HWFI + STM_SPD
  ) + brmsfamily(family = "lognormal", link = "identity")

default_prior(formulaVMAX, data = StormdataTrainYeo)

# priorsVMAX <- c(
#   #prior(horseshoe(1), class = "b")
#   prior(normal(0, 5), class = "b"),
#   prior(inv_gamma(0.1, 0.1), class = "sigma")
#   #prior(inv_gamma(0.1, 0.1), class = "shape"),
#   #prior(inv_gamma(0.1, 0.1), class = "sd")
# )

## Fit brms ----
iters <- 4000
burn <- 2000
chains <- 4
sims <- (iters-burn)*chains

#system.time(
Fit <- brm(
  formulaVMAX,
  data = StormdataTrainYeo,
  prior = c(
    #prior(horseshoe(1), class = "b")
    prior(normal(0, 5), class = "b"),
    prior(inv_gamma(0.1, 0.1), class = "sigma"),
    #prior(inv_gamma(0.1, 0.1), class = "b", dpar = "sigma", lb = 0),
    #prior(inv_gamma(0.1, 0.1), class = "shape"),
    prior(inv_gamma(0.1, 0.1), class = "sd")
  ),
  save_pars = save_pars(all = TRUE), 
  chains = chains,
  iter = iters,
  cores = parallel::detectCores(),
  seed = 52,
  warmup = burn,
  #init = 0,
  normalize = TRUE,
  control = list(adapt_delta = 0.95),
  backend = "cmdstanr"
)
#)

#logNormalFit <- Fit
#logNormalRandFit <- Fit
#gammaFit <- Fit
#gammaRandFit <- Fit

#save(logNormalFit, file = "_data/logNormalFit.RData")
#save(logNormalRandFit, file = "_data/logNormalRandFit.RData")
#save(gammaFit, file = "_data/gammaFit.RData")
#save(gammaRandFit, file = "_data/gammaRandFit.RData")

## Diagnostics ----
fit <- 2
assign(paste0("Fit", fit), Fit)
logNormalFitFINAL <- Fit1
save(logNormalFitFINAL, file = "_data/logNormalFitFINAL.RData")

logNormalFit <- logNormalFitFINAL
plot(logNormalFit, ask = FALSE)
#prior_summary(logNormalFit)

print(logNormalFit, digits = 4)

# waicList <- list(
#   waic(logNormalFit),
#   waic(logNormalRandFit),
#   waic(gammaFit),
#   waic(gammaRandFit)
# )
#waic <- waic(logNormalFit)
#attributes(waic)$model_name <- paste0("logNormalFit", fit)
#waicList[[paste0("fit", 4)]] <- waic4

### Compare Candidate Models ----
# waicList
# 
waicListComp <- loo_compare(
  waicList
)
waicListComp <- waicListComp |>
  data.frame() |>
  rownames_to_column(var = "Model")

save(waicList, waicListComp, file = "_data/waicComps.RData")

pp_check(Fit1, ndraws = 100)
pp_check(Fit2, ndraws = 100)

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
    #Fit2,
    method = "posterior_predict",
    re_formula = NULL,
    robust = FALSE,
    probs = c(0.025, 0.975)) |>
  data.frame()
mean(abs(fitResiduals$Estimate))

predResiduals <- 
  residuals(
    logNormalFit,
    #Fit2,
    #newdata = StormdataTestArcsinh,
    newdata = StormdataTestYeo,
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
set.seed(52)
logNormalFitfinalFit <- posterior_predict(logNormalFit)
# logNormalFitfinalResiduals <- t(StormdataTrain3$VMAX - t(logNormalFitfinalFit))
# logNormalFitfinalResidualsMean <- colMeans(logNormalFitfinalResiduals)
logNormalFitfinalFitMean <- colMeans(logNormalFitfinalFit)
logNormalFitfinalFitMed <- apply(logNormalFitfinalFit, 2, function(x){quantile(x, 0.5)})
logNormalFitfinalFitLCB <- apply(logNormalFitfinalFit, 2, function(x){quantile(x, 0.025)})
logNormalFitfinalFitUCB <- apply(logNormalFitfinalFit, 2, function(x){quantile(x, 0.975)})

### Prediction ----
set.seed(52)
logNormalFitfinalPreds <- posterior_predict(logNormalFit, 
                                            newdata = StormdataTestYeo,
                                            #newdata = StormdataTestYeo,
                                            allow_new_levels = TRUE, 
                                            re_formula = NULL)
logNormalFitfinalPredsMean <- colMeans(logNormalFitfinalPreds)
logNormalFitfinalPredsMed <- apply(logNormalFitfinalPreds, 2, function(x){quantile(x, 0.5)})
logNormalFitfinalPredsLCB <- apply(logNormalFitfinalPreds, 2, function(x){quantile(x, 0.025)})
logNormalFitfinalPredsUCB <- apply(logNormalFitfinalPreds, 2, function(x){quantile(x, 0.975)})


## Plot Fit ----
### PPC ----
trainVMAX <- StormdataTrain$VMAX
testVMAX <- Actual_Yvec

ppc_dens_overlay(y = trainVMAX,
                 yrep = logNormalFitfinalFit[sample(1:sims, 100, replace = FALSE), ]) +
  labs(title = paste0("Fit", fit, " PPC"))

logNormalFitfinalppcFit <- ppc_dens_overlay(
  y = trainVMAX,
  yrep = logNormalFitfinalFit[sample(1:sims, 1000, replace = FALSE), ]
) +
  labs(title = "Simulated density of draws from the PPD vs Observed VMAX",
       subtitle = "n = 1000 draws",
       x = "VMAX",
       y = "Density") +
  scale_x_continuous(limits = c(0, 250), breaks = seq(0, 250, 25)) +
  theme_bw()
logNormalFitfinalppcFit

ppc_dens_overlay(y = testVMAX,
                 yrep = logNormalFitfinalPreds[sample(1:sims, 100, replace = FALSE), ]) +
  labs(title = paste0("Fit", fit, " PPD"))

#### Quantile 2.5 ----
logNormalFitfinalLCBsims <- apply(logNormalFitfinalFit, 
                                  MARGIN = 1,
                                  function(x){
                                    quantile(x, 0.025)
                                  })
logNormalFitfinalLCBpvalueVec <- logNormalFitfinalLCBsims < quantile(StormdataTrain$VMAX, 0.025)
logNormalFitfinalLCBpvalue <- sum(logNormalFitfinalLCBpvalueVec)
logNormalFitfinalLCBpvalue <- round(logNormalFitfinalLCBpvalue/sims, 3)
logNormalFitfinalLCBpvalue <- min(logNormalFitfinalLCBpvalue, 1 - logNormalFitfinalLCBpvalue)

logNormalFitfinal_ppcLCB <- 
  ppc_stat(StormdataTrain$VMAX,
           logNormalFitfinalFit,
           stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = paste0("2.5% Quantile (p-val = ", logNormalFitfinalLCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitfinal_ppcLCB

t <- ppc_stat(
  y = trainVMAX,
  yrep = logNormalFitfinalFit[sample(1:sims, 1000, replace = FALSE), ], 
  stat = function(y) quantile(y, 0.025)
)

#### Quantile 97.5 ----
logNormalFitfinalUCBsims <- apply(logNormalFitfinalFit, 
                                  MARGIN = 1,
                                  function(x){
                                    quantile(x, 0.975)
                                  })
logNormalFitfinalUCBpvalueVec <- logNormalFitfinalUCBsims < quantile(StormdataTrain$VMAX, 0.975)
logNormalFitfinalUCBpvalue <- as.numeric(sum(logNormalFitfinalUCBpvalueVec))
logNormalFitfinalUCBpvalue <- round(logNormalFitfinalUCBpvalue/sims, 3)
logNormalFitfinalUCBpvalue <- min(logNormalFitfinalUCBpvalue, 1 - logNormalFitfinalUCBpvalue)

logNormalFitfinal_ppcUCB <- 
  ppc_stat(StormdataTrain$VMAX,
           logNormalFitfinalFit,
           stat = function(y) quantile(y, 0.975), freq = FALSE) +
  labs(title = paste0("97.5% Quantile (p-val = ", logNormalFitfinalUCBpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitfinal_ppcUCB

#### Mean ----
logNormalFitfinalMEANsims <- apply(logNormalFitfinalFit, 
                                   MARGIN = 1,
                                   function(x){
                                     mean(x)
                                   })
logNormalFitfinalMEANpvalueVec <- logNormalFitfinalMEANsims < mean(StormdataTrain$VMAX)
logNormalFitfinalMEANpvalue <- sum(logNormalFitfinalMEANpvalueVec)
logNormalFitfinalMEANpvalue <- round(logNormalFitfinalMEANpvalue/sims, 3)
logNormalFitfinalMEANpvalue <- min(logNormalFitfinalMEANpvalue, 1 - logNormalFitfinalMEANpvalue)

logNormalFitfinal_ppcMEAN <- 
  ppc_stat(StormdataTrain$VMAX,
           logNormalFitfinalFit,
           stat = function(y) mean(y), freq = FALSE) +
  labs(title = paste0("Mean (p-val = ", logNormalFitfinalMEANpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitfinal_ppcMEAN
save(logNormalFitfinal_ppcMEAN, file = "_data/PPCplotMean.RData")

#### Med ----
logNormalFitfinalMEDsims <- apply(logNormalFitfinalFit, 
                                  MARGIN = 1,
                                  function(x){
                                    quantile(x, 0.5)
                                  })
logNormalFitfinalMEDpvalueVec <- logNormalFitfinalMEDsims < quantile(StormdataTrain$VMAX, 0.5)
logNormalFitfinalMEDpvalue <- sum(logNormalFitfinalMEDpvalueVec)
logNormalFitfinalMEDpvalue <- round(logNormalFitfinalMEDpvalue/sims, 3)
logNormalFitfinalMEDpvalue <- min(logNormalFitfinalMEDpvalue, 1 - logNormalFitfinalMEDpvalue)

logNormalFitfinal_ppcMED <- 
  ppc_stat(StormdataTrain$VMAX,
           logNormalFitfinalFit,
           stat = function(y) quantile(y, 0.5), freq = FALSE) +
  labs(title = paste0("Median (p-val = ", logNormalFitfinalMEDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitfinal_ppcMED

#### SD ----
logNormalFitfinalSDsims <- apply(logNormalFitfinalFit, 
                                 MARGIN = 1,
                                 function(x){
                                   sd(x)
                                 })
logNormalFitfinalSDpvalueVec <- logNormalFitfinalSDsims < sd(StormdataTrain$VMAX)
logNormalFitfinalSDpvalue <- sum(logNormalFitfinalSDpvalueVec)
logNormalFitfinalSDpvalue <- round(logNormalFitfinalSDpvalue/sims, 3)
logNormalFitfinalSDpvalue <- min(logNormalFitfinalSDpvalue, 1 - logNormalFitfinalSDpvalue)

logNormalFitfinal_ppcSD <- 
  ppc_stat(StormdataTrain$VMAX,
           logNormalFitfinalFit,
           stat = function(y) sd(y), freq = FALSE) +
  labs(title = paste0("SD (p-val = ", logNormalFitfinalSDpvalue, ")")) +
  theme_bw() +
  legend_none()
#logNormalFitfinal_ppcSD
save(logNormalFitfinal_ppcSD, file = "_data/PPCplotSD.RData")

#### Range ----
logNormalFitfinalRANGEsims <- apply(logNormalFitfinalFit, 
                                    MARGIN = 1,
                                    function(x){
                                      max(x)-min(x)
                                    })
logNormalFitfinalRANGEpvalueVec <- logNormalFitfinalRANGEsims < (max(StormdataTrain$VMAX)-min(StormdataTrain$VMAX))
logNormalFitfinalRANGEpvalue <- sum(logNormalFitfinalRANGEpvalueVec)
logNormalFitfinalRANGEpvalue <- round(logNormalFitfinalRANGEpvalue/sims, 3)
logNormalFitfinalRANGEpvalue <- min(logNormalFitfinalRANGEpvalue, 1 - logNormalFitfinalRANGEpvalue)

logNormalFitfinal_ppcRANGE <- 
  ppc_stat(StormdataTrain$VMAX,
           logNormalFitfinalFit,
           stat = function(y) max(y)-min(y), freq = FALSE) +
  labs(title = paste0("Range (p-val = ", logNormalFitfinalRANGEpvalue, ")")) +
  theme_bw() +
  legend_none()
save(logNormalFitfinal_ppcRANGE, file = "_data/PPCplotRange.RData")

### Combined Plot ----
logNormalFitfinal_ppcComb <- 
  logNormalFitfinalppcFit /
  (logNormalFitfinal_ppcLCB | logNormalFitfinal_ppcMED | logNormalFitfinal_ppcUCB) /
  (logNormalFitfinal_ppcRANGE | logNormalFitfinal_ppcMEAN | logNormalFitfinal_ppcSD)
logNormalFitfinal_ppcComb

logNormalFitfinal_ppcComb <- 
  logNormalFitfinalppcFit /
  (logNormalFitfinal_ppcMEAN | logNormalFitfinal_ppcSD | logNormalFitfinal_ppcRANGE)
logNormalFitfinal_ppcComb

save(logNormalFitfinal_ppcComb,
     file = "_data/PPCplot.RData")

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
logNormalFitkfoldgroup <- kfold(logNormalFit,
                                folds = kfoldID,
                                chains = 1,
                                save_fits = TRUE)
save(logNormalFitkfoldgroup,
     file = "_data/logNormalFitKfold.RData")
logNormalFitkfoldPreds <- kfold_predict(logNormalFitkfoldgroup)
logNormalFitkfoldPredsDat <- logNormalFitkfoldPreds$yrep
logNormalFitkfoldPredsMean <- colMeans(logNormalFitkfoldPredsDat)
logNormalFitkfoldPredsMed <- apply(logNormalFitkfoldPredsDat, 2, function(x){quantile(x, 0.5)})
logNormalFitkfoldPredsLCB <- apply(logNormalFitkfoldPredsDat, 2, function(x){quantile(x, 0.025)})
logNormalFitkfoldPredsUCB <- apply(logNormalFitkfoldPredsDat, 2, function(x){quantile(x, 0.975)})

logNormalFitkfoldMetrics <- tibble(
  Fit = paste0("logNormalFit", fit),
  MAE_HWRF = mean(abs(StormdataTrain$HWRF - StormdataTrain$VMAX)),
  MAE_kfold = mean(abs(logNormalFitkfoldPredsMean - logNormalFitkfoldPreds$y)),
  MAD_kfold = mean(abs(logNormalFitkfoldPredsMed - logNormalFitkfoldPreds$y)),
  COV_kfold = mean(logNormalFitkfoldPredsLCB < logNormalFitkfoldPreds$y & logNormalFitkfoldPreds$y < logNormalFitkfoldPredsUCB)
)
logNormalFitkfoldMetrics
save(logNormalFitkfoldMetrics,
     file = "_data/logNormalFitKfoldMetrics.RData")

## State-Space Plots ----  
### Training ----
logNormalFitFitDF <- bind_cols(
  StormdataTrain,
  LCB = logNormalFitfinalFitLCB,
  Mean = logNormalFitfinalFitMean,
  Med = logNormalFitfinalFitMed,
  UCB = logNormalFitfinalFitUCB
)

fillPPC <- "#d1e1ec"
colorPPC <- "#b3cde0"
fill2PPC <- "#011f4b"

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

fillPPC <- "#d1e1ec"
colorPPC <- "#b3cde0"
fill2PPC <- "#011f4b"
"lightskyblue2"

#### All ----
logNormalFitstormsPredplot <- ggplot(data = logNormalFitPredDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB, fill = "95% CI"), alpha = 0.5) +
  geom_line(aes(y = VMAX, color = "Observed"), linewidth = 1, alpha = 0.5) +
  #geom_line(aes(y = HWRF, color = "HWRF")) +
  geom_line(aes(y = Mean, color = "PPD Mean"), linewidth = 0.75) +
  scale_y_continuous(limits = c(0,250), breaks = seq(0,275,50)) +
  facet_wrap(vars(StormID))+#, ncol = 6)+
  labs(title = "PPD Mean vs Observed VMAX",
       subtitle = "95% Credible Interval about PPD Mean",
       x = "Storm Elapsed Time") +
  scale_color_manual(name = NULL, 
                     values = c(
                       "Observed" = "#011f4b",
                       "HWRF" = "coral",
                       "PPD Mean" = "springgreen3"
                     )
  ) +
  scale_fill_manual(name = NULL, 
                     values = c(
                       "95% CI" = "lightblue"
                     )
  ) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1), order = 2)
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )
logNormalFitstormsPredplot

#### Single ----
logNormalFitstormsPredplot2 <- 
  ggplot(data = logNormalFitPredDF |> 
           filter(StormID == "3612015"), 
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

#### Low, Mid, High ----
predStorms <- logNormalFitPredDF |> 
  group_by(StormID) |>
  summarise(
    MaxTime = max(StormElapsedTime)
  ) |>
  #arrange(desc(MaxTime))
  arrange(MaxTime)
print(predStorms, n = 30)
median(predStorms$MaxTime)
mean(predStorms$MaxTime)
midStorm <- which(predStorms$MaxTime == median(predStorms$MaxTime))[2]

predStormsPlot <- bind_rows(
  head(predStorms, 3),
  predStorms |> slice(midStorm + c(-1, 0, 1)),
  tail(predStorms, 3)
)

predStormsPlotID <- predStormsPlot |> 
  arrange(MaxTime) |>
  mutate(
    StormID = factor(as.character(StormID),
                     levels = as.character(StormID))
  ) |>
  mutate(
    StormLength = case_when(
      MaxTime < 50 ~ "Short",
      MaxTime > 200 ~ "Long",
      .default = "Medium"
      
    ),
    StormLength = factor(StormLength,
                         levels = c("Short", "Medium", "Long"))
  )
predStormsPlotID

logNormalFitPredDF3 <- left_join(
  predStormsPlotID,
  logNormalFitPredDF
)

logNormalFitstormsPredplot3 <- 
  ggplot(data = logNormalFitPredDF3, 
         aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = "lightblue", alpha = 0.7) +
  geom_line(aes(y = VMAX, color = "Observed")) +
  geom_line(aes(y = HWRF, color = "HWRF")) +
  geom_line(aes(y = Mean, color = "PPD Mean")) +
  scale_y_continuous(limits = c(0,250), breaks = seq(0,275,50)) +
  facet_wrap(vars(StormID), dir = "v")+#, ncol = 6)+
  labs(title = "logNormalFit PPD Mean vs Observed VMAX",
       subtitle = "95% Credible Interval about PPD Mean",
       x = "Storm Elapsed Time") +
  scale_color_manual(name = NULL, 
                     breaks = c(
                       "Observed",
                       "HWRF",
                       "PPD Mean"
                       #"HWRF" = "dodgerblue"
                     ), 
                     values = c(
                       "Observed" = "black",
                       "HWRF" = "dodgerblue",
                       "PPD Mean" = "red3"
                       #"HWRF" = "dodgerblue"
                     )
  ) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1))
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )
logNormalFitstormsPredplot3

#### Long 9 ----
predStormsPlot4 <- tail(predStorms, 9)

predStormsPlotID4 <- predStormsPlot4 |> 
  arrange(MaxTime) |>
  mutate(
    StormID = factor(as.character(StormID),
                     levels = as.character(StormID))
  )
predStormsPlotID4

logNormalFitPredDF4 <- left_join(
  predStormsPlotID4,
  logNormalFitPredDF
)

logNormalFitstormsPredplot4 <- 
  ggplot(data = logNormalFitPredDF4, 
         aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB, fill = "95% CI"), alpha = 0.5) +
  geom_line(aes(y = VMAX, color = "Observed"), linewidth = 1, alpha = .5) +
  geom_line(aes(y = HWRF, color = "HWRF"), linewidth = .75) +
  geom_line(aes(y = Mean, color = "PPD Mean"), linewidth = .75) +
  scale_y_continuous(limits = c(0,250), breaks = seq(0,275,50)) +
  facet_wrap(vars(StormID), dir = "h")+#, ncol = 6)+
  labs(title = "PPD Mean vs Observed VMAX for 9 Longest Storms",
       subtitle = "95% Credible Interval about PPD Mean",
       x = "Storm Elapsed Time") +
  scale_color_manual(name = NULL, 
                     breaks = c(
                       "Observed",
                       "HWRF",
                       "PPD Mean"
                       #"HWRF" = "dodgerblue"
                     ), 
                     values = c(
                       "Observed" = "#011f4b",
                       "HWRF" = "coral",
                       "PPD Mean" = "springgreen3"
                       #"HWRF" = "dodgerblue"
                     )
  )  +
  scale_fill_manual(name = NULL, 
                    values = c(
                      "95% CI" = "lightblue"
                    )
  ) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1), order = 2)
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )
logNormalFitstormsPredplot4


[1] "white"                "aliceblue"            "antiquewhite"         "antiquewhite1"       
[5] "antiquewhite2"        "antiquewhite3"        "antiquewhite4"        "aquamarine"          
[9] "aquamarine1"          "aquamarine2"          "aquamarine3"          "aquamarine4"         
[13] "azure"                "azure1"               "azure2"               "azure3"              
[17] "azure4"               "beige"                "bisque"               "bisque1"             
[21] "bisque2"              "bisque3"              "bisque4"              "black"               
[25] "blanchedalmond"       "blue"                 "blue1"                "blue2"               
[29] "blue3"                "blue4"                "blueviolet"           "brown"               
[33] "brown1"               "brown2"               "brown3"               "brown4"              
[37] "burlywood"            "burlywood1"           "burlywood2"           "burlywood3"          
[41] "burlywood4"           "cadetblue"            "cadetblue1"           "cadetblue2"          
[45] "cadetblue3"           "cadetblue4"           "chartreuse"           "chartreuse1"         
[49] "chartreuse2"          "chartreuse3"          "chartreuse4"          "chocolate"           
[53] "chocolate1"           "chocolate2"           "chocolate3"           "chocolate4"          
[57] "coral"                "coral1"               "coral2"               "coral3"              
[61] "coral4"               "cornflowerblue"       "cornsilk"             "cornsilk1"           
[65] "cornsilk2"            "cornsilk3"            "cornsilk4"            "cyan"                
[69] "cyan1"                "cyan2"                "cyan3"                "cyan4"               
[73] "darkblue"             "darkcyan"             "darkgoldenrod"        "darkgoldenrod1"      
[77] "darkgoldenrod2"       "darkgoldenrod3"       "darkgoldenrod4"       "darkgray"            
[81] "darkgreen"            "darkgrey"             "darkkhaki"            "darkmagenta"         
[85] "darkolivegreen"       "darkolivegreen1"      "darkolivegreen2"      "darkolivegreen3"     
[89] "darkolivegreen4"      "darkorange"           "darkorange1"          "darkorange2"         
[93] "darkorange3"          "darkorange4"          "darkorchid"           "darkorchid1"         
[97] "darkorchid2"          "darkorchid3"          "darkorchid4"          "darkred"             
[101] "darksalmon"           "darkseagreen"         "darkseagreen1"        "darkseagreen2"       
[105] "darkseagreen3"        "darkseagreen4"        "darkslateblue"        "darkslategray"       
[109] "darkslategray1"       "darkslategray2"       "darkslategray3"       "darkslategray4"      
[113] "darkslategrey"        "darkturquoise"        "darkviolet"           "deeppink"            
[117] "deeppink1"            "deeppink2"            "deeppink3"            "deeppink4"           
[121] "deepskyblue"          "deepskyblue1"         "deepskyblue2"         "deepskyblue3"        
[125] "deepskyblue4"         "dimgray"              "dimgrey"              "dodgerblue"          
[129] "dodgerblue1"          "dodgerblue2"          "dodgerblue3"          "dodgerblue4"         
[133] "firebrick"            "firebrick1"           "firebrick2"           "firebrick3"          
[137] "firebrick4"           "floralwhite"          "forestgreen"          "gainsboro"           
[141] "ghostwhite"           "gold"                 "gold1"                "gold2"               
[145] "gold3"                "gold4"                "goldenrod"            "goldenrod1"          
[149] "goldenrod2"           "goldenrod3"           "goldenrod4"           "gray"                
[153] "gray0"                "gray1"                "gray2"                "gray3"               
[157] "gray4"                "gray5"                "gray6"                "gray7"               
[161] "gray8"                "gray9"                "gray10"               "gray11"              
[165] "gray12"               "gray13"               "gray14"               "gray15"              
[169] "gray16"               "gray17"               "gray18"               "gray19"              
[173] "gray20"               "gray21"               "gray22"               "gray23"              
[177] "gray24"               "gray25"               "gray26"               "gray27"              
[181] "gray28"               "gray29"               "gray30"               "gray31"              
[185] "gray32"               "gray33"               "gray34"               "gray35"              
[189] "gray36"               "gray37"               "gray38"               "gray39"              
[193] "gray40"               "gray41"               "gray42"               "gray43"              
[197] "gray44"               "gray45"               "gray46"               "gray47"              
[201] "gray48"               "gray49"               "gray50"               "gray51"              
[205] "gray52"               "gray53"               "gray54"               "gray55"              
[209] "gray56"               "gray57"               "gray58"               "gray59"              
[213] "gray60"               "gray61"               "gray62"               "gray63"              
[217] "gray64"               "gray65"               "gray66"               "gray67"              
[221] "gray68"               "gray69"               "gray70"               "gray71"              
[225] "gray72"               "gray73"               "gray74"               "gray75"              
[229] "gray76"               "gray77"               "gray78"               "gray79"              
[233] "gray80"               "gray81"               "gray82"               "gray83"              
[237] "gray84"               "gray85"               "gray86"               "gray87"              
[241] "gray88"               "gray89"               "gray90"               "gray91"              
[245] "gray92"               "gray93"               "gray94"               "gray95"              
[249] "gray96"               "gray97"               "gray98"               "gray99"              
[253] "gray100"              "green"                "green1"               "green2"              
[257] "green3"               "green4"               "greenyellow"          "grey"                
[261] "grey0"                "grey1"                "grey2"                "grey3"               
[265] "grey4"                "grey5"                "grey6"                "grey7"               
[269] "grey8"                "grey9"                "grey10"               "grey11"              
[273] "grey12"               "grey13"               "grey14"               "grey15"              
[277] "grey16"               "grey17"               "grey18"               "grey19"              
[281] "grey20"               "grey21"               "grey22"               "grey23"              
[285] "grey24"               "grey25"               "grey26"               "grey27"              
[289] "grey28"               "grey29"               "grey30"               "grey31"              
[293] "grey32"               "grey33"               "grey34"               "grey35"              
[297] "grey36"               "grey37"               "grey38"               "grey39"              
[301] "grey40"               "grey41"               "grey42"               "grey43"              
[305] "grey44"               "grey45"               "grey46"               "grey47"              
[309] "grey48"               "grey49"               "grey50"               "grey51"              
[313] "grey52"               "grey53"               "grey54"               "grey55"              
[317] "grey56"               "grey57"               "grey58"               "grey59"              
[321] "grey60"               "grey61"               "grey62"               "grey63"              
[325] "grey64"               "grey65"               "grey66"               "grey67"              
[329] "grey68"               "grey69"               "grey70"               "grey71"              
[333] "grey72"               "grey73"               "grey74"               "grey75"              
[337] "grey76"               "grey77"               "grey78"               "grey79"              
[341] "grey80"               "grey81"               "grey82"               "grey83"              
[345] "grey84"               "grey85"               "grey86"               "grey87"              
[349] "grey88"               "grey89"               "grey90"               "grey91"              
[353] "grey92"               "grey93"               "grey94"               "grey95"              
[357] "grey96"               "grey97"               "grey98"               "grey99"              
[361] "grey100"              "honeydew"             "honeydew1"            "honeydew2"           
[365] "honeydew3"            "honeydew4"            "hotpink"              "hotpink1"            
[369] "hotpink2"             "hotpink3"             "hotpink4"             "indianred"           
[373] "indianred1"           "indianred2"           "indianred3"           "indianred4"          
[377] "ivory"                "ivory1"               "ivory2"               "ivory3"              
[381] "ivory4"               "khaki"                "khaki1"               "khaki2"              
[385] "khaki3"               "khaki4"               "lavender"             "lavenderblush"       
[389] "lavenderblush1"       "lavenderblush2"       "lavenderblush3"       "lavenderblush4"      
[393] "lawngreen"            "lemonchiffon"         "lemonchiffon1"        "lemonchiffon2"       
[397] "lemonchiffon3"        "lemonchiffon4"        "lightblue"            "lightblue1"          
[401] "lightblue2"           "lightblue3"           "lightblue4"           "lightcoral"          
[405] "lightcyan"            "lightcyan1"           "lightcyan2"           "lightcyan3"          
[409] "lightcyan4"           "lightgoldenrod"       "lightgoldenrod1"      "lightgoldenrod2"     
[413] "lightgoldenrod3"      "lightgoldenrod4"      "lightgoldenrodyellow" "lightgray"           
[417] "lightgreen"           "lightgrey"            "lightpink"            "lightpink1"          
[421] "lightpink2"           "lightpink3"           "lightpink4"           "lightsalmon"         
[425] "lightsalmon1"         "lightsalmon2"         "lightsalmon3"         "lightsalmon4"        
[429] "lightseagreen"        "lightskyblue"         "lightskyblue1"        "lightskyblue2"       
[433] "lightskyblue3"        "lightskyblue4"        "lightslateblue"       "lightslategray"      
[437] "lightslategrey"       "lightsteelblue"       "lightsteelblue1"      "lightsteelblue2"     
[441] "lightsteelblue3"      "lightsteelblue4"      "lightyellow"          "lightyellow1"        
[445] "lightyellow2"         "lightyellow3"         "lightyellow4"         "limegreen"           
[449] "linen"                "magenta"              "magenta1"             "magenta2"            
[453] "magenta3"             "magenta4"             "maroon"               "maroon1"             
[457] "maroon2"              "maroon3"              "maroon4"              "mediumaquamarine"    
[461] "mediumblue"           "mediumorchid"         "mediumorchid1"        "mediumorchid2"       
[465] "mediumorchid3"        "mediumorchid4"        "mediumpurple"         "mediumpurple1"       
[469] "mediumpurple2"        "mediumpurple3"        "mediumpurple4"        "mediumseagreen"      
[473] "mediumslateblue"      "mediumspringgreen"    "mediumturquoise"      "mediumvioletred"     
[477] "midnightblue"         "mintcream"            "mistyrose"            "mistyrose1"          
[481] "mistyrose2"           "mistyrose3"           "mistyrose4"           "moccasin"            
[485] "navajowhite"          "navajowhite1"         "navajowhite2"         "navajowhite3"        
[489] "navajowhite4"         "navy"                 "navyblue"             "oldlace"             
[493] "olivedrab"            "olivedrab1"           "olivedrab2"           "olivedrab3"          
[497] "olivedrab4"           "orange"               "orange1"              "orange2"             
[501] "orange3"              "orange4"              "orangered"            "orangered1"          
[505] "orangered2"           "orangered3"           "orangered4"           "orchid"              
[509] "orchid1"              "orchid2"              "orchid3"              "orchid4"             
[513] "palegoldenrod"        "palegreen"            "palegreen1"           "palegreen2"          
[517] "palegreen3"           "palegreen4"           "paleturquoise"        "paleturquoise1"      
[521] "paleturquoise2"       "paleturquoise3"       "paleturquoise4"       "palevioletred"       
[525] "palevioletred1"       "palevioletred2"       "palevioletred3"       "palevioletred4"      
[529] "papayawhip"           "peachpuff"            "peachpuff1"           "peachpuff2"          
[533] "peachpuff3"           "peachpuff4"           "peru"                 "pink"                
[537] "pink1"                "pink2"                "pink3"                "pink4"               
[541] "plum"                 "plum1"                "plum2"                "plum3"               
[545] "plum4"                "powderblue"           "purple"               "purple1"             
[549] "purple2"              "purple3"              "purple4"              "red"                 
[553] "red1"                 "red2"                 "red3"                 "red4"                
[557] "rosybrown"            "rosybrown1"           "rosybrown2"           "rosybrown3"          
[561] "rosybrown4"           "royalblue"            "royalblue1"           "royalblue2"          
[565] "royalblue3"           "royalblue4"           "saddlebrown"          "salmon"              
[569] "salmon1"              "salmon2"              "salmon3"              "salmon4"             
[573] "sandybrown"           "seagreen"             "seagreen1"            "seagreen2"           
[577] "seagreen3"            "seagreen4"            "seashell"             "seashell1"           
[581] "seashell2"            "seashell3"            "seashell4"            "sienna"              
[585] "sienna1"              "sienna2"              "sienna3"              "sienna4"             
[589] "skyblue"              "skyblue1"             "skyblue2"             "skyblue3"            
[593] "skyblue4"             "slateblue"            "slateblue1"           "slateblue2"          
[597] "slateblue3"           "slateblue4"           "slategray"            "slategray1"          
[601] "slategray2"           "slategray3"           "slategray4"           "slategrey"           
[605] "snow"                 "snow1"                "snow2"                "snow3"               
[609] "snow4"                "springgreen"          "springgreen1"         "springgreen2"        
[613] "springgreen3"         "springgreen4"         "steelblue"            "steelblue1"          
[617] "steelblue2"           "steelblue3"           "steelblue4"           "tan"                 
[621] "tan1"                 "tan2"                 "tan3"                 "tan4"                
[625] "thistle"              "thistle1"             "thistle2"             "thistle3"            
[629] "thistle4"             "tomato"               "tomato1"              "tomato2"             
[633] "tomato3"              "tomato4"              "turquoise"            "turquoise1"          
[637] "turquoise2"           "turquoise3"           "turquoise4"           "violet"              
[641] "violetred"            "violetred1"           "violetred2"           "violetred3"          
[645] "violetred4"           "wheat"                "wheat1"               "wheat2"              
[649] "wheat3"               "wheat4"               "whitesmoke"           "yellow"              
[653] "yellow1"              "yellow2"              "yellow3"              "yellow4"             
[657] "yellowgreen"  







