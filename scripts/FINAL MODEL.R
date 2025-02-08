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
library(rstanarm)
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
                             -StormElapsedTime
                             #-LAT,
                             #-LON
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
#preProcArcsinh
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
chains <- 2
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
  normalize = FALSE,
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
fit <- 1
assign(paste0("Fit", fit), Fit)
logNormalFitFINAL <- Fit
save(logNormalFitFINAL, file = "_data/logNormalFitFINAL.RData")

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
logNormalFitfinalFit <- posterior_predict(logNormalFit)
# logNormalFitfinalResiduals <- t(StormdataTrain3$VMAX - t(logNormalFitfinalFit))
# logNormalFitfinalResidualsMean <- colMeans(logNormalFitfinalResiduals)
logNormalFitfinalFitMean <- colMeans(logNormalFitfinalFit)
logNormalFitfinalFitMed <- apply(logNormalFitfinalFit, 2, function(x){quantile(x, 0.5)})
logNormalFitfinalFitLCB <- apply(logNormalFitfinalFit, 2, function(x){quantile(x, 0.025)})
logNormalFitfinalFitUCB <- apply(logNormalFitfinalFit, 2, function(x){quantile(x, 0.975)})

### Prediction ----
logNormalFitfinalPreds <- posterior_predict(logNormalFit, 
                                            newdata = StormdataTestYeo,
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
