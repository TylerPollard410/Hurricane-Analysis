Hurricane Analysis
================
Tyler Pollard
2024-08-25

<!-- start custom head snippets, customize with your own _includes/head-custom.html file -->

<!-- Setup Google Analytics -->
<!-- {% include head-custom-google-analytics.html %} -->

<!-- You can set your favicon here -->
<!-- link rel="shortcut icon" type="image/x-icon" href="{{ '/favicon.ico' | relative_url }}" -->

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

- [Exploratory Data Analysis](#exploratory-data-analysis)
- [Model Description](#model-description)
- [Model Comparisons](#model-comparisons)
  - [Model Selection Criteria](#model-selection-criteria)
- [Goodness of Fit](#goodness-of-fit)
- [Variable Importance](#variable-importance)
- [Prediction](#prediction)

A Bayesian analysis of hurricane data from the National Hurricane
Center.

# Exploratory Data Analysis

To begin we start with using the [shiny
app](https://tylerpollard410.shinyapps.io/Hurricane_EDA/) to explore the
data.

# Model Description

The data for this analysis is related to 24-hour ahead forecasts of
hurricane intensity. The response variable, $Y_i$, is the observed
maximum wind velocity (VMAX) of a given hurricane StormID $S$ at
observation $i$. The objective of this Bayesian analysis is to predict
missing values of VMAX that are better predictions than the current
model which is captured by the HWRF covariate in the data. The $p = 11$
covariates $X_{ij}$ that were fit to the model were LAT, MINSLP,
SHR_MAG, STM_SPD, RHL0, CAPE3, INST2, TCONDSYM2, COUPLSYM3, HWFI, and
HWRF. For an observation in storm $S$, we assume the linear model

$$
\begin{aligned}
Y_i &= \beta_{0S} + \sum_{j=1}^{p}X_{ij}\beta_{jS} + \epsilon_{i} \\
\end{aligned}
$$

where $\beta_{jS}$ is the effect of covariate $j$ in Storm $S$. We
compare three models for the $\beta_{jS}$

1.) Constant slopes: $\beta_{jS} \equiv \beta_j$ for all Storms

2.) Varying slopes with uninformative priors:
$\beta_{jS} \sim Normal(0, 100)$

3.) Varying slopes with informative priors:
$\beta_{jS} \sim Normal(\mu_j, \sigma_j^2)$

Although the data for VMAX itself is not Gaussian, the errors are when
it is regressed onto the predictors. Also with flat priors and large $n$
the Bayesian central limit can be invoked to approximate the posterior
as Gaussian.

# Model Comparisons

All three models above were fit to the data and I first checked that the
models converged to an acceptable level. Once I confirmed that models
did in fact converge, as seen in Figure 1, I compared the models based
on DIC and WAIC to determine the best fitting model of the three. Based
on Table 1, I chose to go forward with the slopes as random effects
model because it had the smallest mean deviance which is a good sign
that it will predict the best, although it was not as simple as the
first model with constant slopes.

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
