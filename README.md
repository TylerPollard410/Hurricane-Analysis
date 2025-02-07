Hurricane Analysis
================
Tyler Pollard
2024-08-25

<!-- start custom head snippets, customize with your own _includes/head-custom.html file -->

<!-- Setup Google Analytics -->
<!-- {% include head-custom-google-analytics.html %} -->

<!-- You can set your favicon here -->
<!-- link rel="shortcut icon" type="image/x-icon" href="{{ '/favicon.ico' | relative_url }}" -->

<!-- Change content width onfull screen -->
<!-- <link rel="stylesheet" href="/Hurricane-Analysis/assets/css/custom.css"> -->


<!-- MathJax -->
<!-- inline config -->
<script>
  MathJax = {
    tex: {
      inlineMath: [['$', '$'], ['\\(', '\\)']],
      macros: {
      	RR: "{\\bf R}",
      	bold: ["{\\bf #1}", 1],
        indep: "{\\perp \\!\\!\\! \\perp}",
    	}
    },
    svg: {
    fontCache: 'global'
  	},
  };
</script>

<!-- load MathJax -->
<script type="text/javascript" id="MathJax-script" async
  src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js">
</script>

<!-- end custom head snippets -->

- [Introduction](#introduction)
- [Data Exploration (Shiny App)](#data-exploration-shiny-app)
- [Model Description](#model-description)
  - [Lognormal](#lognormal)
  - [Gamma](#gamma)
- [Model Comparisons](#model-comparisons)
  - [Model Selection Criteria](#model-selection-criteria)
- [Goodness of Fit](#goodness-of-fit)
- [Variable Importance](#variable-importance)
- [Prediction](#prediction)

# Introduction

The data for this analysis are from [A Feed Forward Neural Network Based
on Model Output Statistics for Short-Term Hurricane Intensity
Prediction](https://journals.ametsoc.org/view/journals/wefo/34/4/waf-d-18-0173_1.xml).
The data variables descriptions are available
[here](https://github.com/TylerPollard410/Hurricane-Analysis/blob/main/docs/data_description.pdf).
This paper uses deep learning to improve 24-hour ahead forecasts of
hurricane intensity (maximum wind velocity, VMAX). The main prediction
model is
[HWRF](https://www.aoml.noaa.gov/hurricane-weather-research-forecast-model/),
which is a mathematical model based on differential equations. In
addition to the forecast, HWRF has many other state variables such as
sea surface temperature, longitude, time of year, etc, that are usually
discarded. In this analysis, we will determine if including these state
variables can improve the HWRF prediction.

# Data Exploration (Shiny App)

To begin this analysis, I created a [Hurricane Analysis
app](https://tylerpollard410.shinyapps.io/Hurricane_EDA/) to better
understand and explore the data. The distribution of VMAX was plotted
using the “Histogram” plot type, which revealed positive right skewness.
The importance of spatial predictors like basin, Land, LAT, and LON can
be viewed using the “Map” plot type. VMAX was plotted against all
potential covariates and their transformed values to explore their
relationship with the response and identify the best transformations to
apply prior to model fitting.

# Model Description

The response variable, $Y_i$, is the observed VMAX of observation $i$
with StormID index $S_i$ of a given hurricane. The $p = 20$ covariates
$X_{ij}$ that were fit to the model were LAT, LON, basin, Land, MINSLP,
SHR_MAG, STM_SPD, SST, RHL0, CAPE1, CAPE3, SHTFL2, TCOND7002, INST2,
CP1, TCONDSYM2, COUPLSYM3, HWFI, VMAX_OP_T0, and HWRF. All covariates
were normalized using a Yeo-Johnson transformation with centering and
scaling. Due to the positive right skewness of VMAX, models were
compared for Log-normal and Gamma likelihoods both with and without
random effects.

## Lognormal

$$
\begin{aligned}
Y_{i} &\sim \text{LogNormal}(\mu_i, \sigma^2) \\
\mu_i &= \beta_{0} + \sum_{j=1}^{p}X_{ij}\beta_{j} + \theta_{S_i} \\
\theta_{S_i} &\sim \text{Normal}(0, \tau^2)
\end{aligned}
$$

## Gamma

$$
\begin{aligned}
Y_{i} &\sim \text{Gamma}(\alpha, \alpha/\mu_i) \\
log(\mu_i) &= \beta_{0} + \sum_{j=1}^{p}X_{ij}\beta_{j} + \theta_{S_i} \\
\theta_{S_i} &\sim \text{Normal}(0, \tau^2)
\end{aligned}
$$

where $\beta_{j}$ is the effect of covariate $j$ with weakly informative
priors $\beta_j \sim \text{Normal}(0,5)$. We set a non-informative prior
on $\alpha \sim \text{InvGamma}(0.1, 0.1)$. In the random effects
models, we set non-informative priors
$\tau \sim \text{InvGamma}(0.1, 0.1)$.

# Model Comparisons

All four models above were fit to the data and were first checked that
the models converged to an acceptable level. Once convergence was
confirmed, the models were compared based on WAIC to determine the model
with the best expected predictive accuracy. Based on Table 1, I chose to
go forward with the slopes as random effects model because it had the
smallest mean deviance which is a good sign that it will predict the
best, although it was not as simple as the first model with constant
slopes.

## Model Selection Criteria

# Goodness of Fit

To indicate the goodness of fit of the model, posterior predictive
checks were simulated as the data was being simulated. For every
iteration of the MCMC, a draw of size $n$ was made from the posterior
predictive distribution of $Y_i$ and the mean, standard deviation, and
max were calculated to also get a draw from their distribution to
compare directly back to the observed values of VMAX. The Bayesian
p-values for these posterior predictive checks and they are far from 0
or 1, so this is a good fitting model.

# Variable Importance

After fitting the model with random slopes with informative priors, the
importance of the covariates changed drastically in terms of what is
most important. The mean, standard deviation, and 95% credible set on
each parameter is displayed for both the posterior of model 1 and model
3 after pooling over the storm IDs. In the model that treated the slopes
as constant across storms, every parameter had a significant effect
other than MINSLP according the the 95% credible interval. However, in
model 3 when the slopes are treated as random effects and estimated from
the data, the only important variables are HWRF and HWFI. We were told
that HWRF would be a very important variable, but we may have another
one that can improve the predictions of VMAX across storms. Since VMAX
is log-normally distributed, both HWRF and HWFI increase VMAX about 8%.

# Prediction

The posterior predictive distributions were output for the $\beta_{js}$
and $\sigma_i$ parameters. Using these PPD’s to draw new values for the
new values for the covariates and StormIDs that were not present in the
training data set, new data was simulated. For each new observation, the
posterior predictive distribution was determined and posterior
predictive checks were conducted for the mean and 95% credible set from
$S = 20000$ iterations of MCMC sampling. The Bayesian p-values from the
posterior predictive checks were 0.35029940, 0.68263473, and 0.04041916
for the mean, lower bound, and upper bound, respectively when compared
to the observed values of the original data set. This is not a
one-to-one comparison, but gives a good idea of a decent fit and
prediction of the new observations.
