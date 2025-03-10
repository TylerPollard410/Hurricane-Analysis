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


### Arcsinh ----
StormdataTrainArcsinhPre <- StormdataTrain |>
  mutate(
    across(c(where(is.numeric), 
             -VMAX, #-HWRF,
             -StormElapsedTime, 
             -LAT2, 
             -LON2),
           ~predict(arcsinh_x(.x, standardize = FALSE), newdata = .x)
    ))

StormdataTestArcsinhPre <- StormdataTest |>
  mutate(
    across(c(where(is.numeric), 
             -VMAX, #-HWRF,
             -StormElapsedTime, 
             -LAT2, 
             -LON2),
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

# Model Comparison ############################################################
## Shared Model Parameters ----
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

iters <- 4000
burn <- 2000
chains <- 4
sims <- (iters-burn)*chains

## Log-Normal ====
formula_logNormal <- 
  bf(VMAX ~ 
       LAT + 
       LON + 
       basin +
       Land +  
       MINSLP + 
       SHR_MAG +
       STM_SPD + 
       SST + 
       RHLO + 
       CAPE1 + 
       CAPE3 + 
       SHTFL2 + 
       #TCOND7002 + 
       #INST2 + 
       CP1 + 
       TCONDSYM2 +
       COUPLSYM3 + 
       HWFI +
       VMAX_OP_T0 +  
       HWRF #+  
     #(1 | StormID)  # Random effect for storm-specific variation
     #sigma ~ HWRF #+ HWFI + STM_SPD
  ) + brmsfamily(family = "lognormal", link = "identity")

default_prior(formula_logNormal, data = StormdataTrainArcsinh)

priors_logNormal <- c(
  #prior(horseshoe(1), class = "b")
  prior(normal(0, 5), class = "b"),
  prior(inv_gamma(0.1, 0.1), class = "sigma")
  #prior(inv_gamma(0.1, 0.1), class = "shape"),
  #prior(inv_gamma(0.1, 0.1), class = "sd")
)

system.time(
  logNormal_Fit <- brm(
    formula_logNormal,
    data = StormdataTrainArcsinh,
    prior = priors_logNormal,
    save_pars = save_pars(all = TRUE), 
    chains = chains,
    iter = iters,
    warmup = burn,
    cores = parallel::detectCores(),
    #init = 0,
    normalize = TRUE,
    control = list(adapt_delta = 0.95),
    backend = "rstan",
    seed = 52
  )
)

### Multicollinearity
check_collinearity(logNormal_Fit) # Drop TCOND7002, INST2

# Save Model
save(logNormal_Fit, file = "_data/logNormal_Fit.RData")


## Log-Normal Rand ====
formula_logNormal_Rand <- 
  bf(VMAX ~ 
       LAT + 
       LON + 
       basin +
       Land +  
       MINSLP + 
       SHR_MAG +
       STM_SPD + 
       SST + 
       RHLO + 
       CAPE1 + 
       CAPE3 + 
       SHTFL2 + 
       #TCOND7002 + 
       #INST2 + 
       CP1 + 
       TCONDSYM2 +
       COUPLSYM3 + 
       HWFI +
       VMAX_OP_T0 +  
       HWRF +  
       (1 | StormID)  # Random effect for storm-specific variation
     #sigma ~ HWRF #+ HWFI + STM_SPD
  ) + brmsfamily(family = "lognormal", link = "identity")

default_prior(formula_logNormal_Rand, data = StormdataTrainArcsinh)

priors_logNormal_Rand <- c(
  #prior(horseshoe(1), class = "b")
  prior(normal(0, 5), class = "b"),
  prior(inv_gamma(0.1, 0.1), class = "sigma"),
  #prior(inv_gamma(0.1, 0.1), class = "shape"),
  prior(inv_gamma(0.1, 0.1), class = "sd")
)

system.time(
  logNormal_Rand_Fit <- brm(
    formula_logNormal_Rand,
    data = StormdataTrainArcsinh,
    prior = priors_logNormal_Rand,
    save_pars = save_pars(all = TRUE), 
    chains = chains,
    iter = iters,
    warmup = burn,
    cores = parallel::detectCores(),
    #init = 0,
    normalize = TRUE,
    control = list(adapt_delta = 0.95),
    backend = "rstan",
    seed = 52
  )
)

### Multicollinearity
check_collinearity(logNormal_Rand_Fit) # Drop TCOND7002, INST2

# Save Model
save(logNormal_Rand_Fit, file = "_data/logNormal_Rand_Fit.RData")


## Gamma ====
formula_gamma <- 
  bf(VMAX ~ 
       LAT + 
       LON + 
       basin +
       Land +  
       MINSLP + 
       SHR_MAG +
       STM_SPD + 
       SST + 
       RHLO + 
       CAPE1 + 
       CAPE3 + 
       SHTFL2 + 
       #TCOND7002 + 
       #INST2 + 
       CP1 + 
       TCONDSYM2 +
       COUPLSYM3 + 
       HWFI +
       VMAX_OP_T0 +  
       HWRF #+  
     #(1 | StormID)  # Random effect for storm-specific variation
     #sigma ~ HWRF #+ HWFI + STM_SPD
  ) + brmsfamily(family = "Gamma", link = "log")

default_prior(formula_gamma, data = StormdataTrainArcsinh)

priors_gamma <- c(
  #prior(horseshoe(1), class = "b")
  prior(normal(0, 5), class = "b"),
  #prior(inv_gamma(0.1, 0.1), class = "sigma")
  prior(inv_gamma(0.1, 0.1), class = "shape")
  #prior(inv_gamma(0.1, 0.1), class = "sd")
)

system.time(
  gamma_Fit <- brm(
    formula_gamma,
    data = StormdataTrainArcsinh,
    prior = priors_gamma,
    save_pars = save_pars(all = TRUE), 
    chains = chains,
    iter = iters,
    warmup = burn,
    cores = parallel::detectCores(),
    #init = 0,
    normalize = TRUE,
    control = list(adapt_delta = 0.95),
    backend = "rstan",
    seed = 52
  )
)

### Multicollinearity
check_collinearity(gamma_Fit) # Drop TCOND7002, INST2

# Save Model
save(gamma_Fit, file = "_data/gamma_Fit.RData")


## Gamma Rand ====
formula_gamma_Rand <- 
  bf(VMAX ~ 
       LAT + 
       LON + 
       basin +
       Land +  
       MINSLP + 
       SHR_MAG +
       STM_SPD + 
       SST + 
       RHLO + 
       CAPE1 + 
       CAPE3 + 
       SHTFL2 + 
       #TCOND7002 + 
       #INST2 + 
       CP1 + 
       TCONDSYM2 +
       COUPLSYM3 + 
       HWFI +
       VMAX_OP_T0 +  
       HWRF +  
       (1 | StormID)  # Random effect for storm-specific variation
     #sigma ~ HWRF #+ HWFI + STM_SPD
  ) + brmsfamily(family = "Gamma", link = "log")

default_prior(formula_gamma_Rand, data = StormdataTrainArcsinh)

priors_gamma_Rand <- c(
  #prior(horseshoe(1), class = "b")
  prior(normal(0, 5), class = "b"),
  #prior(inv_gamma(0.1, 0.1), class = "sigma"),
  prior(inv_gamma(0.1, 0.1), class = "shape"),
  prior(inv_gamma(0.1, 0.1), class = "sd")
)

system.time(
  gamma_Rand_Fit <- brm(
    formula_gamma_Rand,
    data = StormdataTrainArcsinh,
    prior = priors_gamma_Rand,
    save_pars = save_pars(all = TRUE), 
    chains = chains,
    iter = iters,
    warmup = burn,
    cores = parallel::detectCores(),
    #init = 0,
    normalize = TRUE,
    control = list(adapt_delta = 0.95),
    backend = "rstan",
    seed = 52
  )
)

### Multicollinearity
check_collinearity(gamma_Rand_Fit) # Drop TCOND7002, INST2

# Save Model
save(gamma_Rand_Fit, file = "_data/gamma_Rand_Fit.RData")


## WAIC ------
waicList <- list(
  waic(logNormal_Fit),
  waic(logNormal_Rand_Fit),
  waic(gamma_Fit),
  waic(gamma_Rand_Fit)
)

### Compare Models ----
waicList

waicListComp <- loo_compare(
  waicList
)
waicListComp <- waicListComp |>
  data.frame() |>
  rownames_to_column(var = "Model") |>
  mutate(
    Model = case_when(
      Model == "logNormal_Fit" ~ "Log-normal (Fixed Effects)",
      Model == "logNormal_Rand_Fit" ~ "Log-normal (Random Effects)",
      Model == "gamma_Fit" ~ "Gamma (Fixed Effects)",
      Model == "gamma_Rand_Fit" ~ "Gamma (Random Effects)"
    )
  )

save(waicList, waicListComp, file = "_data/waicComps.RData")

### gt table ----
waicListComp |>
  select(
    Model, waic, p_waic, elpd_waic, elpd_diff, se_diff
  ) |>
  gt() |>
  # tab_header(
  #   title = "Table 1: Model Selection Criteria"
  # ) |>
  fmt_number(
    decimals = 2,
    use_seps = FALSE
  ) |>
  cols_align(
    align = "center",
    columns = -Model
  ) |>
  tab_footnote(
    footnote = "WAIC (waic): Widely Applicable Information Criterion, a model selection metric balancing fit and complexity. Lower values indicate better expected predictive accuracy.",
    locations = cells_column_labels(columns = waic)
  ) |>
  tab_footnote(
    footnote = "Effective Parameters (p_waic): An estimate of the number of effective parameters in the model; higher values indicate more flexibility.",
    locations = cells_column_labels(columns = p_waic)
  ) |>
  tab_footnote(
    footnote = "ELPD (elpd_waic): Expected log predictive density, quantifying out-of-sample predictive performance. Higher (less negative) values indicate better predictive accuracy.",
    locations = cells_column_labels(columns = elpd_waic)
  ) |>
  tab_footnote(
    footnote = "ELPD Difference (elpd_diff): The difference in elpd_waic relative to the best model (logNormal_Rand_Fit). The best model always has elpd_diff = 0.",
    locations = cells_column_labels(columns = elpd_diff)
  ) |>
  tab_footnote(
    footnote = "SE of Difference (se_diff): The standard error of elpd_diff, measuring uncertainty in the difference estimates. Large absolute elpd_diff values relative to se_diff indicate meaningful performance differences.",
    locations = cells_column_labels(columns = se_diff)
  ) |>
  tab_style(
    style = list(
      cell_borders(sides = "right")
    ),
    locations = list(
      cells_body(columns = Model),
      cells_column_labels(columns = Model)
    )
  ) |>
  tab_options(
    table.align = "center"
  ) |>
  tab_style(
    style = list(
      cell_fill(color = "#373737"),
      cell_text(color = "#fff")
    ),
    locations = list(
      cells_column_labels()
    )
  ) |>
  tab_options(
    column_labels.padding = "10px",
    table.background.color = "#f2f2f2"
  ) |>
  #opt_vertical_padding(scale = 0.25) |>
  tab_caption(caption = md("Table 1: Model Selection Criteria Using WAIC and ELPD")) |>
  as_raw_html()

## FINAL Model =============================================================
iters <- 4000
burn <- 2000
chains <- 4
sims <- (iters-burn)*chains

formula_Final <- 
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

default_prior(formula_Final, data = StormdataTrainArcsinh)

priors_Final <- c(
  #prior(horseshoe(1), class = "b")
  prior(normal(0, 5), class = "b"),
  prior(inv_gamma(0.1, 0.1), class = "sigma"),
  #prior(inv_gamma(0.1, 0.1), class = "shape"),
  prior(inv_gamma(0.1, 0.1), class = "sd")
)

system.time(
  final_Fit <- brm(
    formula_Final,
    data = StormdataTrainArcsinh,
    prior = priors_Final,
    save_pars = save_pars(all = TRUE), 
    chains = chains,
    iter = iters,
    warmup = burn,
    cores = parallel::detectCores(),
    #init = 0,
    normalize = TRUE,
    control = list(adapt_delta = 0.95),
    backend = "rstan",
    seed = 52
  )
)

plot(final_Fit)
print(final_Fit, digits = 4)
pp_check(final_Fit, ndraws = 100)

save(final_Fit, file = "_data/final_Fit.RData")

### Fixed Effects ----
fixedEff <- fixef(final_Fit)
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

### MAE ----
fitResiduals <- 
  residuals(
    final_Fit,
    #Fit2,
    method = "posterior_predict",
    re_formula = NULL,
    robust = FALSE,
    probs = c(0.025, 0.975)) |>
  data.frame()
mean(abs(fitResiduals$Estimate))


# Goodness of Fit ##########################################################
## Posteriors ----
### Training ----
set.seed(52)
posteriorFit <- posterior_predict(final_Fit)
posteriorFitMean <- colMeans(posteriorFit)
posteriorFitMed <- apply(posteriorFit, 2, function(x){quantile(x, 0.5)})
posteriorFitLCB <- apply(posteriorFit, 2, function(x){quantile(x, 0.025)})
posteriorFitUCB <- apply(posteriorFit, 2, function(x){quantile(x, 0.975)})

## PPC ----
trainVMAX <- StormdataTrain$VMAX
testVMAX <- Actual_Yvec

set.seed(52)
sampFitID <- sample(1:sims, 1000, replace = FALSE)
postSampFit <- posteriorFit[sampFitID, ]

fillPPC <- "#d1e1ec"
colorPPC <- "#b3cde0"
fill2PPC <- "#011f4b"

### Density -----
ppcDensPlot <- ppc_dens_overlay(
  y = trainVMAX,
  yrep = postSampFit
) +
  labs(
    # title = "Simulated density of draws from the PPD vs Observed VMAX",
    # subtitle = "n = 1000 draws",
    x = "VMAX",
    y = "Density"
  ) +
  scale_x_continuous(limits = c(0, 250), breaks = seq(0, 250, 25)) +
  theme_bw() +
  theme(
    legend.position = "none"
  )


#### Stats ----
# Make stat functions
meanFunc <- function(y){mean(y)}
sdFunc <- function(y){sd(y)}
rangeFunc <- function(y){max(y) - min(y)}

##### Mean ----
set.seed(52) # for reproducibility
meanY <- meanFunc(trainVMAX)

ppcMeanStat <- ppc_stat_data(
  y = trainVMAX,
  yrep = posteriorFit,
  group = NULL,
  stat = c("meanFunc")
) |>
  mutate(
    meanProbLow = value < meanY,
    meanProbHigh = value > meanY
  )

ppcMeanPlotGG <- ggplot() +
  geom_histogram(
    data = ppcMeanStat |> filter(variable != "y"),
    aes(x = value, color = "Posterior"),
    fill = fillPPC
  ) +
  geom_vline(
    data = ppcMeanStat |> filter(variable == "y"),
    aes(xintercept = value, color = "Observed"),
    linewidth = 1
  ) +
  scale_x_continuous(
    name = "VMAX"
  ) +
  scale_y_continuous(
    name = "Number of Posterior Draws",
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_color_manual(
    name = "Data",
    values = c(
      "Posterior" = colorPPC,
      "Observed" = "black"
    ),
    breaks = c("Posterior", "Observed")
  ) +
  labs(title = "Mean",
       subtitle = paste("p-value =", round(mean(ppcMeanStat$meanProbLow[-1]), 4))
  ) +
  theme_bw() +
  guides(color = guide_legend(byrow = TRUE)) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.key.spacing.y = unit(5, "points"),
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

##### SD ----
set.seed(52) # for reproducibility
sdY <- sdFunc(trainVMAX)

ppcSDStat <- ppc_stat_data(
  y = trainVMAX,
  yrep = posteriorFit,
  group = NULL,
  stat = c("sdFunc")
) |>
  mutate(
    sdProbLow = value < sdY,
    sdProbHigh = value > sdY
  )

ppcSDPlotGG <- ggplot() +
  geom_histogram(
    data = ppcSDStat |> filter(variable != "y"),
    aes(x = value, color = "Posterior"),
    fill = fillPPC
  ) +
  geom_vline(
    data = ppcSDStat |> filter(variable == "y"),
    aes(xintercept = value, color = "Observed"),
    linewidth = 1
  ) +
  scale_x_continuous(
    name = "VMAX"
  ) +
  scale_y_continuous(
    name = "Number of Posterior Draws",
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_color_manual(
    name = "Data",
    values = c(
      "Posterior" = colorPPC,
      "Observed" = "black"
    ),
    breaks = c("Posterior", "Observed")
  ) +
  labs(title = "SD",
       subtitle = paste("p-value =", round(mean(ppcSDStat$sdProbLow[-1]), 4))
  ) +
  theme_bw() +
  guides(color = guide_legend(byrow = TRUE)) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.key.spacing.y = unit(5, "points"),
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

##### Range ----
set.seed(52) # for reproducibility
rangeY <- rangeFunc(trainVMAX)

ppcRangeStat <- ppc_stat_data(
  y = trainVMAX,
  yrep = posteriorFit,
  group = NULL,
  stat = c("rangeFunc")
) |>
  mutate(
    rangeProbLow = value < rangeY,
    rangeProbHigh = value > rangeY
  )

ppcRangePlotGG <- ggplot() +
  geom_histogram(
    data = ppcRangeStat |> filter(variable != "y"),
    aes(x = value, color = "Posterior"),
    fill = fillPPC
  ) +
  geom_vline(
    data = ppcRangeStat |> filter(variable == "y"),
    aes(xintercept = value, color = "Observed"),
    linewidth = 1
  ) +
  scale_x_continuous(
    name = "VMAX"
  ) +
  scale_y_continuous(
    name = "Number of Posterior Draws",
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_color_manual(
    name = "Data",
    values = c(
      "Posterior" = colorPPC,
      "Observed" = "black"
    ),
    breaks = c("Posterior", "Observed")
  ) +
  labs(title = "Range",
       subtitle = paste("p-value =", round(mean(ppcRangeStat$rangeProbLow[-1]), 4))
  ) +
  theme_bw() +
  guides(color = guide_legend(byrow = TRUE)) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.key.spacing.y = unit(5, "points"),
    legend.position = "bottom",
    legend.direction = "horizontal"
  )


#### Combined plot ----
ppcCombPlot <- 
  ppcDensPlot /
  (ppcMeanPlotGG + ppcSDPlotGG + ppcRangePlotGG) +
  plot_layout(
    guides = "collect",
    axes = "collect_x"
  ) +
  plot_annotation(
    #title = "Posterior Predictive Checks for Distributional Statistics",
    #subtitle = paste("Bayesian predictive p-values for", sims, "Simulations"),
    theme = theme(
      legend.position = "bottom",
      legend.direction = "horizontal"
    )
  )
ppcCombPlot

# Variable Importance #########################################################
finalFitSum <- summary(final_Fit, digits = 4)
finalFitFixed <- finalFitSum$fixed
finalFitRand <- finalFitSum$random$StormID

## gt table ----
finalFitFixed |>
  select(
    Mean = Estimate, 
    SD = Est.Error, 
    `Q2.5` = `l-95% CI`,
    `Q97.5` = `u-95% CI`
  ) |>
  round(digits = 4) |>
  rownames_to_column(var = "Parameter") |>
  gt() |>
  # tab_header(
  #   title = "Table 2: Posterior summary for model parameters"
  # ) |>
  cols_align(
    align = "center",
    columns = c(everything(), -Parameter)
  ) |>
  tab_style(
    style = list(
      cell_borders(sides = "right")
    ),
    locations = list(
      cells_body(columns = Parameter),
      cells_column_labels(columns = Parameter)
    )
  ) |>
  tab_style(
    style = list(
      cell_fill(color = "#373737"),
      cell_text(color = "#fff")
    ),
    locations = list(
      cells_column_labels()
    )
  ) |>
  tab_options(
    column_labels.padding = "10px",
    table.background.color = "#f2f2f2"
  ) |>
  #opt_vertical_padding(scale = 0.25) |>
  tab_caption(caption = md("Table 2: Posterior summary for model parameters"))

# Prediction #################################################################
## Observed Data ====
### Cross Validation ----
set.seed(52)
kfoldID <- kfold_split_grouped(K = 5, StormdataTrain$StormID)
system.time(
  final_cv <- kfold(final_Fit,
                    folds = kfoldID,
                    chains = 1,
                    save_fits = TRUE,
                    seed = 52)
)

save(final_cv,
     file = "_data/final_cv.RData")

### Performance ----
trainHWRF <- StormdataTrain$HWRF

#### Model ----
final_Fit_predMetrics <- tibble(
  Method = "Model PPD",
  MAE_HWRF = mean(abs(trainHWRF - trainVMAX)),
  MAE_Model = mean(abs(posteriorFitMean - trainVMAX)),
  COV = mean(posteriorFitLCB < trainVMAX & trainVMAX < posteriorFitUCB)
)

#### CV ----
set.seed(52)
final_cv_Preds <- kfold_predict(final_cv)
final_cv_PredsDat <- final_cv_Preds$yrep
final_cv_PredsMean <- colMeans(final_cv_PredsDat)
final_cv_PredsMed <- apply(final_cv_PredsDat, 2, function(x){quantile(x, 0.5)})
final_cv_PredsLCB <- apply(final_cv_PredsDat, 2, function(x){quantile(x, 0.025)})
final_cv_PredsUCB <- apply(final_cv_PredsDat, 2, function(x){quantile(x, 0.975)})

final_cv_Metrics <- tibble(
  Fit = "logNormal_Rand",
  MAE_HWRF = mean(abs(trainHWRF - trainVMAX)),
  MAE_kfold = mean(abs(final_cv_PredsMean - final_cv_Preds$y)),
  MAD_kfold = mean(abs(final_cv_PredsMed - final_cv_Preds$y)),
  COV_kfold = mean(final_cv_PredsLCB < final_cv_Preds$y & final_cv_Preds$y < final_cv_PredsUCB)
)
final_cv_Metrics
save(final_cv_Metrics,
     file = "_data/final_cv_Metrics.RData")


### Compare ----
predMetricDF <- bind_rows(
  final_Fit_predMetrics,
  final_cv_Metrics |> 
    select(
      MAE_HWRF,
      MAE_Model = MAE_kfold,
      COV = COV_kfold
    ) |>
    round(digits = 3) |>
    mutate(
      Method = "5-fold CV",
      .before = 1
    )
)

#### gt table ----
predMetricDF |>
  gt() |>
  # tab_header(
  #   title = "Table 3: Prediction metrics"
  # ) |>
  fmt_number(
    decimals = 3
  ) |>
  cols_label(
    MAE_HWRF = "HWRF MAE",
    MAE_Model = "Model MAE"
  ) |>
  cols_align(
    align = "center",
    columns = -Method
  ) |>
  tab_style(
    style = list(
      cell_borders(sides = "right")
    ),
    locations = list(
      cells_body(columns = Method),
      cells_column_labels(columns = Method)
    )
  ) |>
  tab_style(
    style = list(
      cell_fill(color = "#373737"),
      cell_text(color = "#fff")
    ),
    locations = list(
      cells_column_labels()
    )
  ) |>
  tab_options(
    column_labels.padding = "10px",
    table.background.color = "#f2f2f2"
  ) |>
  #opt_vertical_padding(scale = 0.25) |>
  tab_caption(caption = md("Table 3: Prediction metrics")) 

### State-Space Plots ----  
observedPosteriorDF <- bind_cols(
  StormdataTrain,
  LCB = posteriorFitLCB,
  Mean = posteriorFitMean,
  Med = posteriorFitMed,
  UCB = posteriorFitUCB
)

fillPPD <- "lightblue"
colorObserved <- "#011f4b"
colorHWRF <- "coral"
colorPPD <- "springgreen3"
# fillPPD <- "#d1e1ec"
# colorPPD <- "#b3cde0"
# fill2PPD <- "#011f4b"

observedPPD_plot <- ggplot(data = observedPosteriorDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB), fill = fillPPD) +
  geom_line(aes(y = VMAX, color = "Observed")) +
  geom_line(aes(y = Mean, color = "PPD Mean")) +
  facet_wrap(vars(StormID))+
  scale_y_continuous(limits = c(0,275), breaks = seq(0,275,50)) +
  labs(title = "logNormalFit PPD Mean vs Observed VMAX",
       subtitle = "95% Credible Interval about PPD Mean") +
  scale_color_manual(name = NULL, values = c(colorObserved, colorPPD)) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1))
  ) +
  theme_bw()
observedPPD_plot



## Out-of-Sample =====
### Posterior ----
set.seed(52)
posteriorPred <- posterior_predict(final_Fit,
                                   newdata = StormdataTestArcsinh,
                                   allow_new_levels = TRUE, 
                                   re_formula = NULL)
posteriorPredMean <- colMeans(posteriorPred)
posteriorPredMed <- apply(posteriorPred, 2, function(x){quantile(x, 0.5)})
posteriorPredLCB <- apply(posteriorPred, 2, function(x){quantile(x, 0.025)})
posteriorPredUCB <- apply(posteriorPred, 2, function(x){quantile(x, 0.975)})

#### Prediction Metrics ----
trainVMAX <- StormdataTrain$VMAX
testVMAX <- Actual_Yvec

trainHWRF <- StormdataTrain$HWRF
testHWRF <- StormdataTest$HWRF

posterior_predMetrics <- tibble(
  Method = "Prediction",
  MAE_HWRF = mean(abs(testHWRF - testVMAX)),
  MAE_Model = mean(abs(posteriorPredMean - testVMAX)),
  MAD_Model = mean(abs(posteriorPredMed - testVMAX)),
  COV = mean(posteriorPredLCB < testVMAX & testVMAX < posteriorPredUCB)
)

### State-Space Plots ----  
oosPosteriorDF <- bind_cols(
  StormdataTest,
  LCB = posteriorPredLCB,
  Mean = posteriorPredMean,
  Med = posteriorPredMed,
  UCB = posteriorPredUCB
) |>
  mutate(
    VMAX = Actual_Yvec
  ) 

fillPPD <- "lightblue"
colorActual <- "#011f4b"
colorHWRF <- "coral"
colorPPD <- "springgreen3"

#### All ----
oosPPD_plot_All <- ggplot(data = oosPosteriorDF, aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB, fill = "95% CI"), alpha = 0.5) +
  geom_line(aes(y = VMAX, color = "Actual"), linewidth = 1, alpha = 0.5) +
  #geom_line(aes(y = HWRF, color = "HWRF")) +
  geom_line(aes(y = Mean, color = "PPD Mean"), linewidth = 0.75) +
  scale_y_continuous(limits = c(0,250), breaks = seq(0,275,50)) +
  facet_wrap(vars(StormID))+#, ncol = 6)+
  labs(title = "PPD Mean vs Observed VMAX",
       subtitle = "95% Credible Interval about PPD Mean",
       x = "Storm Elapsed Time") +
  scale_color_manual(name = NULL, 
                     breaks = c(
                       "Actual",
                       #"HWRF",
                       "PPD Mean"
                     ),
                     values = c(
                       "Actual" = colorActual,
                       #"HWRF" = colorHWRF,
                       "PPD Mean" = colorPPD
                     )
  ) +
  scale_fill_manual(name = NULL, 
                    values = c(
                      "95% CI" = fillPPD
                    )
  ) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1), order = 2)
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )
oosPPD_plot_All


#### Long 9 ----
predStorms <- oosPosteriorDF |> 
  group_by(StormID) |>
  summarise(
    MaxTime = max(StormElapsedTime)
  ) |>
  #arrange(desc(MaxTime))
  arrange(MaxTime)

predStormsLong <- tail(predStorms, 9)

predStormsLongID <- predStormsLong |> 
  arrange(MaxTime) |>
  mutate(
    StormID = factor(as.character(StormID),
                     levels = as.character(StormID))
  )
predStormsLongID

oosPosteriorDF_Long <- left_join(
  predStormsLongID,
  oosPosteriorDF
)

oosPPD_plot_Long <- 
  ggplot(data = oosPosteriorDF_Long, 
         aes(x = StormElapsedTime)) +
  geom_ribbon(aes(ymin = LCB, ymax = UCB, fill = "95% CI"), alpha = 0.5) +
  geom_line(aes(y = VMAX, color = "Actual"), linewidth = 1, alpha = .5) +
  geom_line(aes(y = HWRF, color = "HWRF"), linewidth = .75) +
  geom_line(aes(y = Mean, color = "PPD Mean"), linewidth = .75) +
  scale_y_continuous(limits = c(0,250), breaks = seq(0,275,50)) +
  facet_wrap(vars(StormID), dir = "h")+#, ncol = 6)+
  labs(title = "PPD Mean vs Observed VMAX for 9 Longest Storms",
       subtitle = "95% Credible Interval about PPD Mean",
       x = "Storm Elapsed Time") +
  scale_color_manual(name = NULL, 
                     breaks = c(
                       "Actual",
                       "HWRF",
                       "PPD Mean"
                       #"HWRF" = "dodgerblue"
                     ), 
                     values = c(
                       "Actual" = colorActual,
                       "HWRF" = colorHWRF,
                       "PPD Mean" = colorPPD
                     )
  )  +
  scale_fill_manual(name = NULL, 
                    values = c(
                      "95% CI" = fillPPD
                    )
  ) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1), order = 2)
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )
oosPPD_plot_Long











