#### ST540 Exam 2
## Tyler Pollard
## 4 April 2024

# Load Libraries ----
library(knitr)
library(data.table)
library(MASS)
library(rjags)
library(plyr)
library(stringr)
library(tidyverse)
library(tictoc)
library(caret)

library(DescTools)
library(bayesplot)
library(BayesFactor)
library(brms)
library(rstanarm)
library(tidybayes)

library(tidyverse)

# Read in data ----
## Clean data ----
Stormdata <- fread("E2_data.csv")
Stormdata$basin <- ifelse(Stormdata$basin == "atlantic", 1, 0)

# get unique values
lapply(Stormdata, function(x){length(unique(x))})

Stormdata2 <- Stormdata |>
  select(-lead_time)



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
    lead_time,
    VMAX,
    basin,
    LON,
    SST,
    CAPE1,
    SHTFL2,
    TCOND7002,
    CP1,
    VMAX_OP_T0
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

#### Linear Model ===================================================
##### Constant Slopes ########################################################################
model_string1 <- textConnection(
  "model{
  
  # Likelihood
  for(i in 1:n){
    Y[i] ~ dnorm(muY[i], taue)
    muY[i] <- inprod(X[i,], beta[])
  }
  
  # Priors
  for(j in 1:p){
    beta[j] ~ dnorm(0, 0.01)
  }
  taue ~ dgamma(0.1, 0.1)
  sigma <- sqrt(1/taue)
  
  # Posterior predictive checks
  for(i in 1:n){
    Y1[i] ~ dnorm(muY[i], taue)
  }
  DY1[1] <- min(Y1[])
  DY1[2] <- max(Y1[])
  DY1[3] <- max(Y1[]) - min(Y1[])
  DY1[4] <- mean(Y1[])
  DY1[5] <- sd(Y1[])
  
  # WAIC calculations
  for(i in 1:n){
    like[i] <- dnorm(Y[i], muY[i], taue)
  }
  }"
)

data1 <- list(Y = Y, X = Xscale, n = n, p = p)

###### Set seed for reproducibility ----
inits <- list(
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 52),
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 8)
)

model1 <- jags.model(model_string1, data = data1, inits = inits, 
                     n.chains = chains, n.adapt = 0, quiet = TRUE)

tic()
update(model1, n.iter = burn, progress.bar = "none")

params1 <- c("beta" 
             # "sigma", 
             # "Y1", 
             # "DY1"
)
samples1 <- coda.samples(model1, variable.names = params1, n.iter = iters, 
                         progress.bar = "none", quiet = TRUE)
toc()

# samp1Names <- colnames(samples1[[1]])
# #sampNames
# samples1df <- samples1[[1]]
# 
# ## Get column indexes 
# beta1Names <- which(str_detect(samp1Names, "beta"))
# sigma1Names <- which(str_detect(samp1Names, "sigma"))
# DY1Names <- which(str_detect(samp1Names, "DY1"))
# Y1Names <- which(str_detect(samp1Names, "Y1"))[-c(1:5)]
# 
# ## Subset samples
# beta1Samps <- samples1df[,beta1Names]
# sigma1Samps <- samples1df[,sigma1Names]
# DY1Samps <- samples1df[,DY1Names]
# Y1Samps <- samples1df[,Y1Names]
# 
# param1Names <- colnames(Xscale)
# DPrintnames <- c("Min of Y", "Max of Y", "Range of Y", "Mean of Y", "SD of Y")
# colnames(beta1Samps) <- param1Names
# param1Samps <- samples1df[,c(beta1Names, sigma1Names)]
# colnames(param1Samps) <- c(param1Names, "sigma")

###### Convergence Checks ----
ESS1 <- effectiveSize(samples1)
out1 <- summary(samples1)$quantiles
rownames(out1) <- colnames(Xscale)
which(sign(out1[,1]) != sign(out1[,5]))

# Compute DIC
dic1    <- dic.samples(model1,n.iter=iters,progress.bar="none")

# Compute WAIC
waic1   <- coda.samples(model1, 
                        variable.names=c("like"), 
                        n.iter=iters, progress.bar="none")
like1   <- waic1[[1]]
fbar1   <- colMeans(like1)
P1      <- sum(apply(log(like1),2,var))
WAIC1   <- -2*sum(log(fbar1))+2*P1

##### Fixed Effect Slopes########################################################################
model_string2 <- textConnection(
  "model{
  
  # Likelihood
  for(i in 1:n){
    Y[i] ~ dnorm(muY[i], taue)
    muY[i] <- inprod(X[i,], beta[StormID[i],])
  }
  
  # Slopes
  for(j in 1:p){
    for(i in 1:N){
      beta[i,j] ~ dnorm(0, 0.01)
    }
  }
  
  # Priors
  taue ~ dgamma(0.1, 0.1)
  sigma <- sqrt(1/taue)
  
  
  # Posterior predictive mean
  for(i in 1:n){
    Y2[i] ~ dnorm(muY[i], taue)
  }
  DY2[1] <- min(Y2[])
  DY2[2] <- max(Y2[])
  DY2[3] <- max(Y2[]) - min(Y2[])
  DY2[4] <- mean(Y2[])
  DY2[5] <- sd(Y2[])
  
   # WAIC calculations
   for(i in 1:n){
     like[i]    <- dnorm(Y[i],muY[i],taue)
   }
  }"
)

data2 <- list(Y = Y, X = Xscale, n = n, p = p,
              N = N, StormID = StormID)

###### Set seed for reproducibility ----
inits <- list(
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 52),
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 8)
)

model2 <- jags.model(model_string2, data = data2, inits = inits, 
                     n.chains = chains, n.adapt = 0, quiet = TRUE)

tic()
update(model2, n.iter = burn, progress.bar = "none")

params2 <- c("beta" 
             # "sigma", 
             # "Y2", "DY2"
)
samples2 <- coda.samples(model2, variable.names = params2, n.iter = iters, 
                         progress.bar = "none", quiet = TRUE)
toc()

# samp2Names <- colnames(samples2[[1]])
# #sampNames
# samples2df <- samples2[[1]]
# 
# ## Get column indexes 
# alpha2Names <- which(str_detect(samp2Names, "alpha"))
# beta2Names <- which(str_detect(samp2Names, "beta"))
# sigma2Names <- which(str_detect(samp2Names, "sigma"))
# DY2Names <- which(str_detect(samp2Names, "DY2"))
# Y2Names <- which(str_detect(samp2Names, "Y2"))[-c(1:5)]
# 
# ## Subset samples
# alpha2Samps <- samples2df[, alpha2Names]
# beta2Samps <- samples2df[,beta2Names]
# sigma2Samps <- samples2df[,sigma2Names]
# DY2Samps <- samples2df[,DY2Names]
# Y2Samps <- samples2df[,Y2Names]

# param2Names <- colnames(Xscale)
# DPrintnames <- c("Min of Y", "Max of Y", "Range of Y", "Mean of Y", "SD of Y")
# #Stormnames <- paste0("Storm", Storms$StormID)
# colnames(beta2Samps) <- param2Names
# param2Samps <- samples2df[,c(beta2Names, sigma2Names)]
# colnames(param2Samps) <- c("Intercept", param2Names, "sigma")

###### Convergence Checks ----
ESS2 <- effectiveSize(samples2)
# autocorr.diag(samples1, lag = 1)
# geweke.diag(samples1)
# gelman.diag(samples1)

#out2 <- summary(beta2Samps)$quantiles
sum2      <- summary(samples2)$stat
post_mn2 <- matrix(sum2[,1],N,p)
post_sd2 <- matrix(sum2[,2],N,p)

# Compute DIC
dic2   <- dic.samples(model2,n.iter=iters,progress.bar="none")

# Compute WAIC
waic2   <- coda.samples(model2, 
                        variable.names=c("like"), 
                        n.iter=iters, progress.bar="none")
like2   <- waic2[[1]]
fbar2   <- colMeans(like2)
P2      <- sum(apply(log(like2),2,var))
WAIC2   <- -2*sum(log(fbar2))+2*P2

##### Random Effect Slopes########################################################################
model_string3 <- textConnection(
  "model{
  
  # Likelihood
  for(i in 1:n){
    Y[i] ~ dnorm(muY[i], taue)
    muY[i] <- inprod(X[i,], beta[StormID[i],])
  }
  
  # Slopes
  for(j in 1:p){
    for(i in 1:N){
      beta[i,j] ~ dnorm(mu[j], taub[j])
    }
    mu[j] ~ dnorm(0, 0.01)
    taub[j] ~ dgamma(0.1, 0.1)
  }
  
  # Priors
  taue ~ dgamma(0.1, 0.1)
  sigma <- sqrt(1/taue)
  
  # Posterior predictive checks
  for(i in 1:n){
    Y3[i] ~ dnorm(muY[i], taue)
  }
  DY3[1] <- min(Y3[])
  DY3[2] <- max(Y3[])
  DY3[3] <- max(Y3[]) - min(Y3[])
  DY3[4] <- mean(Y3[])
  DY3[5] <- sd(Y3[])
  
   # WAIC calculations
   for(i in 1:n){
     like[i]    <- dnorm(Y[i],muY[i],taue)
   }
  }"
)

data3 <- list(Y = Y, X = Xscale, n = n, p = p,
              N = N, StormID = StormID)

###### Set seed for reproducibility ----
inits <- list(
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 52),
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 8)
)

model3 <- jags.model(model_string3, data = data3, inits = inits, 
                     n.chains = chains, n.adapt = 0, quiet = TRUE)

tic()
update(model3, n.iter = burn, progress.bar = "none")

params3 <- c("beta",
             "sigma", 
             "Y3", 
             "DY3"
)
samples3 <- coda.samples(model3, variable.names = params3, n.iter = iters, 
                         progress.bar = "none", quiet = TRUE)
toc()

# samp3Names <- colnames(samples3[[1]])
# #sampNames
# samples3df <- samples3[[1]]
# 
# ## Get column indexes 
# alpha3Names <- which(str_detect(samp3Names, "alpha"))
# beta3Names <- which(str_detect(samp3Names, "beta"))
# sigma3Names <- which(str_detect(samp3Names, "sigma"))
# DY3Names <- which(str_detect(samp3Names, "DY3"))
# Y3Names <- which(str_detect(samp3Names, "Y3"))[-c(1:5)]
# 
# ## Subset samples
# alpha3Samps <- samples3df[, alpha3Names]
# beta3Samps <- samples3df[,beta3Names]
# sigma3Samps <- samples3df[,sigma3Names]
# DY3Samps <- samples3df[,DY3Names]
# Y3Samps <- samples3df[,Y3Names]
# 
# param3Names <- colnames(Xscale)
# DPrintnames <- c("Min of Y", "Max of Y", "Range of Y", "Mean of Y", "SD of Y")
# #Stormnames <- paste0("Storm", Storms$StormID)
# colnames(beta3Samps) <- param3Names
# param3Samps <- samples3df[,c(beta3Names, sigma3Names)]
# colnames(param3Samps) <- c(param3Names, "sigma")

###### Convergence Checks ----
ESS3 <- effectiveSize(samples3)
# autocorr.diag(samples1, lag = 1)
# geweke.diag(samples1)
# gelman.diag(samples1)

#out3 <- summary(samples3)$quantiles
sum3      <- summary(samples3)$stat
post_mn3 <- matrix(sum3[,1],N,p)
post_sd3 <- matrix(sum3[,2],N,p)
quant <- summary(samples3)$quantiles

# Compute DIC
dic3   <- dic.samples(model3,n.iter=iters,progress.bar="none")

# Compute WAIC
waic3   <- coda.samples(model3, 
                        variable.names=c("like"), 
                        n.iter=iters, progress.bar="none")
like3   <- waic3[[1]]
fbar3   <- colMeans(like3)
P3      <- sum(apply(log(like3),2,var))
WAIC3   <- -2*sum(log(fbar3))+2*P3

dic1
WAIC1
P1

dic2
WAIC2
P2

dic3
WAIC3
P3


# 2. Model comparisons -############################################################
## Save model comparison data ----
fitsum_df <- data.frame(check.names = FALSE,
                        Model = c("Constant Slopes", "Effect Slopes", "Random Effect Slopes"),
                        "Mean Deviance" = c(13307, 12017,12069),
                        "DIC Penalty" = c(12.99, 472.2, 325.9),
                        "DIC Penalized Deviance" = c(13320, 12489, 12395),
                        "WAIC" = c(13322.25, 12502.82, 12399.66),
                        "WAIC Penalty" = c(15.34225, 396.272, 281.0766)
)

save(out1, file = "model1beta_estimates.RData")

save(ESS1, file = "ESS1.RData")
save(ESS2, file = "ESS2.RData")
save(ESS3, file = "ESS3.RData")

hist(ESS1, breaks = 30)
hist(ESS2, breaks = 30)
hist(ESS3, breaks = 30)

# final model means and sds
save(post_mn3, file = "posterior mean.RData")
save(post_sd3, file = "posterior sd.RData")

## Posterior predictive checks ----
### Log Data ----
#### Observed ----
D0_minZ <- min(Z)
D0_maxZ <- max(Z)
D0_rangeZ <- max(Z) - min(Z)
D0_meanZ <- mean(Z)
D0_sdZ <- sd(Z)
D0_Z <- cbind(
  D0_minZ,
  D0_maxZ,
  D0_rangeZ,
  D0_meanZ,
  D0_sdZ
)
colnames(D0_Z) <- DlogPrintnames
D0_Z

D0_minY <- min(Y)
D0_maxY <- max(Y)
D0_rangeY <- max(Y) - min(Y)
D0_meanY <- mean(Y)
D0_sdY <- sd(Y)
D0_Y <- cbind(
  D0_minY,
  D0_maxY,
  D0_rangeY,
  D0_meanY,
  D0_sdY
)
colnames(D0_Y) <- DPrintnames
D0_Y


#### Predicted ----
##### JAGS ----
D1_Z <- DZ1Samps
colnames(D1_Z) <- DlogPrintnames

D1_Y <- DY1Samps
colnames(D1_Y) <- DPrintnames

D2_Y <- DY2Samps
colnames(D2_Y) <- DPrintnames

D3_Y <- DY3Samps
colnames(D3_Y) <- DPrintnames

## Plot posterior predicitve checks ----
pval1 <- rep(0, 5)
names(pval1) <- DPrintnames
pval2 <- rep(0, 5)
names(pval2) <- DPrintnames
pval3 <- rep(0, 5)
names(pval3) <- DPrintnames

for(j in 1:5){
  plot(density(D3_Y[,j]), xlab = "D", ylab = "Posterior Probability",
       xlim = c(min(D1_Y[,j], D0_Y[j], D2_Y[,j], D3_Y[,j]), 
                max(D1_Y[,j], D0_Y[j], D2_Y[,j], D3_Y[,j])), 
       main = DPrintnames[j])
  lines(density(D2_Y[,j]), col = "red")
  lines(density(D1_Y[,j]), col = "blue")
  abline(v = D0_Y[j], col = "green", lwd = 2)
  legend("topleft", c("D3", "D2", "D1", "Observed"), 
         col = c("black", "red", "blue", "green"), lwd = 2)
  
  pval1[j] <- mean(D1_Y[,j] > D0_Y[j])
  pval2[j] <- mean(D2_Y[,j] > D0_Y[j])
  pval3[j] <- mean(D3_Y[,j] > D0_Y[j])
}
pval1
pval2
pval3


## ggplot ----
D0_YGG <- data.frame(D0_Y, check.names = FALSE)
D1_YGG <- data.frame(D1_Y, check.names = FALSE)

D0_YGG2 <- D0_YGG |>
  pivot_longer(everything(), cols_vary = "slowest", names_to = "Dvar", values_to = "D") |>
  add_column("Dataset" = "D0")
D1_YGG2 <- D1_YGG |>
  pivot_longer(everything(), cols_vary = "slowest", names_to = "Dvar", values_to = "D") |>
  add_column("Dataset" = "D1")
D_YGG <- rbind(D0_YGG2, D1_YGG2)

DallGG <- rbind(D_ZGG, D_YGG)
DallGG$Dvar <- factor(DallGG$Dvar,
                      levels = c(DlogPrintnames, DPrintnames))

ggplot() +
  geom_density(data = DallGG |> filter(Dataset %in% c("D1")),
               aes(x = D, color = Dataset)) +
  geom_vline(data = DallGG |> filter(Dataset == "D0"),
             aes(xintercept = D, color = Dataset),
             linewidth = 1) +
  #facet_wrap(vars(Dvar), nrow = 2, scales = "free") +
  facet_wrap(vars(Dvar), nrow = 5, scales = "free", dir = "v") +
  labs(x = "D",
       y = "Posterior Density",
       title = "PPD Check for no Error added") +
  theme_bw()


# 3. Goodness of fit ######################################################
##### Random Effect Slopes ----
Stormdata <- fread("E2_data.csv")
Stormdata$basin <- ifelse(Stormdata$basin == "atlantic", 1, 0)

## Create training and test data sets ----
training_data <- Stormdata |> filter(complete.cases(VMAX))
test_data <- Stormdata |> filter(!complete.cases(VMAX))

## Load data ----
Y <- training_data$VMAX
Z <- log(Y)
X <- training_data |>
  select(-c(
    StormID,
    Date,
    lead_time,
    VMAX,
    basin,
    LON,
    SST,
    CAPE1,
    SHTFL2,
    TCOND7002,
    CP1,
    VMAX_OP_T0
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


model_string3 <- textConnection(
  "model{
  
  # Likelihood
  for(i in 1:n){
    Y[i] ~ dnorm(muY[i], taue)
    muY[i] <- inprod(X[i,], beta[StormID[i],])
  }
  
  # Slopes
  for(j in 1:p){
    for(i in 1:N){
      beta[i,j] ~ dnorm(mu[j], taub[j])
    }
    mu[j] ~ dnorm(0, 0.01)
    taub[j] ~ dgamma(0.1, 0.1)
  }
  
  # Priors
  taue ~ dgamma(0.1, 0.1)
  sigma <- sqrt(1/taue)
  
  # Posterior predictive checks
  for(i in 1:n){
    Y3[i] ~ dnorm(muY[i], taue)
  }
  DY3[1] <- min(Y3[])
  DY3[2] <- max(Y3[])
  DY3[3] <- max(Y3[]) - min(Y3[])
  DY3[4] <- mean(Y3[])
  DY3[5] <- sd(Y3[])
  
   # WAIC calculations
   for(i in 1:n){
     like[i]    <- dnorm(Y[i],muY[i],taue)
   }
  }"
)

data3 <- list(Y = Y, X = Xscale, n = n, p = p,
              N = N, StormID = StormID)

###### Set seed for reproducibility ----
inits <- list(
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 52),
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 8)
)

model3 <- jags.model(model_string3, data = data3, inits = inits, 
                     n.chains = chains, n.adapt = 0, quiet = TRUE)

tic()
update(model3, n.iter = burn, progress.bar = "none")

params3 <- c("beta",
             "sigma", 
             "Y3", 
             "DY3"
)
samples3 <- coda.samples(model3, variable.names = params3, n.iter = iters, 
                         progress.bar = "none", quiet = TRUE)
toc()

samples3df1 <- samples3[[1]]
samples3df2 <- samples3[[2]]
samp3Names <- colnames(samples3[[1]])

## Get column indexes
beta3Names <- which(str_detect(samp3Names, "beta"))
sigma3Names <- which(str_detect(samp3Names, "sigma"))
DY3Names <- which(str_detect(samp3Names, "DY3"))
Y3Names <- which(str_detect(samp3Names, "Y3"))[-c(1:5)]

## Subset samples
beta3Samps1 <- samples3df1[,beta3Names]
sigma3Samps1 <- samples3df1[,sigma3Names]
DY3Samps1 <- samples3df1[,DY3Names]
Y3Samps1 <- samples3df1[,Y3Names]

beta3Samps2 <- samples3df2[,beta3Names]
sigma3Samps2 <- samples3df2[,sigma3Names]
DY3Samps2 <- samples3df2[,DY3Names]
Y3Samps2 <- samples3df2[,Y3Names]

param3Names <- colnames(Xscale)
DPrintnames <- c("Min of Y", "Max of Y", "Range of Y", "Mean of Y", "SD of Y")
param3Samps1 <- samples3df1[,c(beta3Names, sigma3Names)]
param3Samps2 <- samples3df2[,c(beta3Names, sigma3Names)]

###### Convergence Checks ----
ESS3 <- effectiveSize(samples3)
# autocorr.diag(samples1, lag = 1)
# geweke.diag(samples1)
# gelman.diag(samples1)

sum3      <- summary(samples3)$stat
post_mn3 <- matrix(sum3[,1],N,p)
post_sd3 <- matrix(sum3[,2],N,p)
quant3 <- summary(samples3)$quantiles

D0_minY <- min(Y)
D0_maxY <- max(Y)
D0_rangeY <- max(Y) - min(Y)
D0_meanY <- mean(Y)
D0_sdY <- sd(Y)
D0_Y <- cbind(
  D0_minY,
  D0_maxY,
  D0_rangeY,
  D0_meanY,
  D0_sdY
)
colnames(D0_Y) <- DPrintnames
D0_Y

save(D0_Y, file = "D0_Y.RData")


D1_Y <- DY3Samps1
colnames(D1_Y) <- DPrintnames
D2_Y <- DY3Samps2
colnames(D2_Y) <- DPrintnames
D3_Y <- rbind(D1_Y, D2_Y)

pval1 <- rep(0, 5)
names(pval1) <- DPrintnames
pval2 <- rep(0, 5)
names(pval2) <- DPrintnames
pval3 <- rep(0, 5)
names(pval3) <- DPrintnames

for(j in 1:5){
  plot(density(D3_Y[,j]), xlab = "D", ylab = "Posterior Probability",
       xlim = c(min(D1_Y[,j], D0_Y[j], D2_Y[,j], D3_Y[,j]), 
                max(D1_Y[,j], D0_Y[j], D2_Y[,j], D3_Y[,j])), 
       main = DPrintnames[j])
  lines(density(D2_Y[,j]), col = "red")
  lines(density(D1_Y[,j]), col = "blue")
  abline(v = D0_Y[j], col = "green", lwd = 2)
  legend("topleft", c("D3", "D2", "D1", "Observed"), 
         col = c("black", "red", "blue", "green"), lwd = 2)
  
  pval1[j] <- mean(D1_Y[,j] > D0_Y[j])
  pval2[j] <- mean(D2_Y[,j] > D0_Y[j])
  pval3[j] <- mean(D3_Y[,j] > D0_Y[j])
}
pval1
pval2
pval3
save(D0_Y, file = "D0_Y.RData")
save(D3_Y, file = "D3_Y.RData")

# 4. Variable importance #######################################################
beta3Sum <- sum3[beta3Names,]
beta3Quant <- quant3[beta3Names,]
post3Beta_mn3 <- matrix(beta3Sum[,1],N,p)
post3Beta_sd3 <- matrix(beta3Sum[,2],N,p)
post3Beta_q2.5 <- matrix(beta3Quant[,1],N,p)
post3Beta_q97.5 <- matrix(beta3Quant[,5],N,p)

post3BetaSumdf <- data.frame(row.names = param3Names,
                             "Mean" = colMeans(post3Beta_mn3),
                             "SD" = colMeans(post3Beta_sd3),
                             "Q2.5" = colMeans(post3Beta_q2.5),
                             "Q97.5" = colMeans(post3Beta_q97.5)
)
save(post3BetaSumdf, file = "Mod3 Post Beta Sum.RData")

load("samples1.RData")
beta1Samps <- samples1[[1]]
beta1SampsSum <- summary(beta1Samps)$stat
beta1SampsQuant <- summary(beta1Samps)$quantiles

post1BetaSumdf <- data.frame(row.names = param3Names,
                             "Mean" = beta1SampsSum[,1],
                             "SD" = beta1SampsSum[,2],
                             "Q2.5" = beta1SampsQuant[,1],
                             "Q97.5" = beta1SampsQuant[,5]
)
save(post1BetaSumdf, file = "Mod1 Post Beta Sum.RData")

sigma3Sum <- sum3[sigma3Names,]
sigma3Quant <- quant3[sigma3Names,]


# 5. Prediction #################################################################
Stormdata <- fread("E2_data.csv")
Stormdata$basin <- ifelse(Stormdata$basin == "atlantic", 1, 0)

### Create StormID indicator for random effects model ----
StormIDs <- Stormdata |> select(StormID)
Storms <- unique(StormIDs$StormID)
IDs <- 1:length(Storms)
Storms <- data.frame(
  StormID = Storms,
  ID = IDs
)
StormIDs <- left_join(StormIDs, Storms)
StormID <- StormIDs$ID

## Create training and test data sets ----
training_data <- Stormdata |> filter(complete.cases(VMAX))
test_data <- Stormdata |> filter(!complete.cases(VMAX))

## Load data ----
XT <- test_data |>
  select(-c(
    StormID,
    Date,
    lead_time,
    VMAX,
    basin,
    LON,
    SST,
    CAPE1,
    SHTFL2,
    TCOND7002,
    CP1,
    VMAX_OP_T0
  ))
XTscale <- scale(XT)
XTscale <- cbind("Intercept" = rep(1,nrow(XTscale)), XTscale)
HWRF <- test_data$HWRF

### Create StormID indicator for random effects model ----
StormIDs <- test_data |> select(StormID)
Storms <- unique(StormIDs$StormID)
IDs <- 1:length(Storms)
Storms <- data.frame(
  StormID = Storms,
  ID = IDs
)
StormIDs <- left_join(StormIDs, Storms)
StormID <- StormIDs$ID

n <- nrow(XTscale)
p <- ncol(Xscale)
N <- length(IDs)

burn <- 10000
iters <- 20000
chains <- 2
thin <- 5

set.seed(52)
betapreds <- matrix(data = 0, nrow = N, ncol = p)
sigmapreds <- rnorm(iters, sigma3Sum[1], sigma3Sum[2])
Ypreds <- matrix(data = NA, nrow = iters, ncol = n)
DYpreds <- matrix(data = NA, nrow = n, ncol = 4)

for(i in 1:N){
  betapreds[i,] <- rnorm(p, post3BetaSumdf[,1], post3BetaSumdf[,2])
}

for(i in 1:n){
  ID <- StormID[i]
  Ypred3 <- XTscale[i,]%*%betapreds[ID,] + rnorm(iters, 0, sigmapreds)
  Ypreds[,i] <- Ypred3
  DYpreds[i,1] <- mean(Ypred3)
  DYpreds[i,2] <- median(Ypred3)
  DYpreds[i,3] <- quantile(Ypred3, 0.025)
  DYpreds[i,4] <- quantile(Ypred3, 0.975)
}

DYpredNames <- c("Mean of YPred", "Median of YPred", "Q2.5 of YPred", "Q97.5 of YPred")
colnames(DYpreds) <- DYpredNames
D0Ypred <- c(mean(Y), median(Y), quantile(Y, 0.025), quantile(Y, 0.975))
names(D0Ypred) <- DYpredNames
pvalPred <- rep(0, 4)
names(pval1) <- colnames(DYpreds)

for(j in 1:4){
  plot(density(DYpreds[,j]), xlab = "D", ylab = "Posterior Probability",
       xlim = c(min(DYpreds[,j], D0Ypred[j]), 
                max(DYpreds[,j], D0Ypred[j])), 
       main = DYpredNames[j])
  abline(v = D0Ypred[j], col = "green", lwd = 2)
  legend("topleft", c("D3", "Observed"), 
         col = c("black", "green"), lwd = 2)
  
  pvalPred[j] <- mean(DYpreds[,j] > D0Ypred[j])
}
pvalPred

DYpredMean <- apply(Ypreds, 2, mean)
DYpredMedian <- apply(Ypreds, 2, median)
DYpredQ2.5 <- apply(Ypreds, 2, quantile, 0.025)
DYpredQ97.5 <- apply(Ypreds, 2, quantile, 0.975)

DYpreds2 <- data.frame(check.names = FALSE,
                       "Post Pred Mean" = DYpredMean,
                       "Post Pred Median" = DYpredMedian,
                       "Post Pred Q2.5" = DYpredQ2.5,
                       "Post Pred Q97.5" = DYpredQ97.5
)

DYpredNames <- colnames(DYpreds2)
D0Ypred <- c(mean(Y), median(Y), quantile(Y, 0.025), quantile(Y, 0.975))
names(D0Ypred) <- DYpredNames
pvalPred <- rep(0, 4)
names(pval1) <- colnames(DYpreds2)

for(j in 1:4){
  plot(density(DYpreds2[,j]), xlab = "D", ylab = "Posterior Probability",
       xlim = c(min(DYpreds2[,j], D0Ypred[j]), 
                max(DYpreds2[,j], D0Ypred[j])), 
       main = DYpredNames[j])
  abline(v = D0Ypred[j], col = "green", lwd = 2)
  legend("topleft", c("D3", "Observed"), 
         col = c("black", "green"), lwd = 2)
  
  pvalPred[j] <- mean(DYpreds[,j] > D0Ypred[j])
}
pvalPred

VMAXnew <- rep(NA, nrow(Stormdata))
VMAXnewL <- rep(NA, nrow(Stormdata))
VMAXnewU <- rep(NA, nrow(Stormdata))

outputDF <- data.frame(
  HWRF = Stormdata$HWRF,
  VMAX = VMAXnew,
  L = VMAXnewL,
  U = VMAXnewU
)

missingOG <- which(is.na(Stormdata$VMAX))
knownOG <- which(!is.na(Stormdata$VMAX))

outputDF$VMAX[knownOG] <- Stormdata$VMAX[knownOG]
outputDF$VMAX[missingOG] <- DYpreds2$`Post Pred Mean`
outputDF$L[missingOG] <- DYpreds2$`Post Pred Q2.5`
outputDF$U[missingOG] <- DYpreds2$`Post Pred Q97.5`

write.csv(outputDF, file = "PollardTyler.csv")

## VMAX - VMAX_pred --------
outputDF <- read.csv("PollardTyler.csv")
## Load data ----
Stormdata2 <- Stormdata
Stormdata2$VMAX <- outputDF$VMAX
Ynew <- Stormdata2$VMAX
Xnew <- Stormdata2 |>
  select(-c(
    StormID,
    Date,
    lead_time,
    VMAX,
    basin,
    LON,
    SST,
    CAPE1,
    SHTFL2,
    TCOND7002,
    CP1,
    VMAX_OP_T0
  ))
XNscale <- scale(Xnew)
XNscale <- cbind("Intercept" = rep(1,length(Ynew)), XNscale)

### Create StormID indicator for random effects model ----
StormIDs <- Stormdata2 |> select(StormID)
Storms <- unique(StormIDs$StormID)
IDs <- 1:length(Storms)
Storms <- data.frame(
  StormID = Storms,
  ID = IDs
)
StormIDs <- left_join(StormIDs, Storms)
StormID <- StormIDs$ID

n <- length(Ynew)
p <- ncol(XNscale)
N <- length(IDs)

burn <- 10000
iters <- 20000
chains <- 2
thin <- 5

model_string4 <- textConnection(
  "model{
  
  # Likelihood
  for(i in 1:n){
    Y[i] ~ dnorm(muY[i], taue)
    muY[i] <- inprod(X[i,], beta[StormID[i],])
  }
  
  # Slopes
  for(j in 1:p){
    for(i in 1:N){
      beta[i,j] ~ dnorm(mu[j], taub[j])
    }
    mu[j] ~ dnorm(0, 0.01)
    taub[j] ~ dgamma(0.1, 0.1)
  }
  
  # Priors
  taue ~ dgamma(0.1, 0.1)
  sigma <- sqrt(1/taue)
  
  # Posterior predictive checks
  for(i in 1:n){
    Y4[i] ~ dnorm(muY[i], taue)
  }
  DY4[1] <- min(Y4[])
  DY4[2] <- max(Y4[])
  DY4[3] <- max(Y4[]) - min(Y4[])
  DY4[4] <- mean(Y4[])
  DY4[5] <- sd(Y4[])
  
   # WAIC calculations
   for(i in 1:n){
     like[i]    <- dnorm(Y[i],muY[i],taue)
   }
  }"
)

inits <- list(
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 52),
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 8)
)

data4 <- list(Y = Ynew, X = XNscale, n = n, p = p, 
              N = N, StormID = StormID)

model4 <- jags.model(model_string4, data = data4, inits = inits, 
                     n.chains = chains, n.adapt = 0, quiet = TRUE)

tic()
update(model4, n.iter = burn, progress.bar = "none")

params4 <- c("beta",
             "sigma", 
             "Y4", 
             "DY4"
)
samples4 <- coda.samples(model4, variable.names = params4, n.iter = iters, 
                         progress.bar = "none", quiet = TRUE)
toc()

# Compute DIC
dic4   <- dic.samples(model4,n.iter=iters,progress.bar="none")

# Compute WAIC
waic4   <- coda.samples(model4, 
                        variable.names=c("like"), 
                        n.iter=iters, progress.bar="none")
like4   <- waic4[[1]]
fbar4   <- colMeans(like4)
P4      <- sum(apply(log(like4),2,var))
WAIC4   <- -2*sum(log(fbar4))+2*P4

dic4
WAIC4
P4

### Cross Validation ----
set.seed(52)
fold <- sample(1:5, 2373, replace = TRUE)

Y_mean   <- matrix(NA,2373)
Y_median <- matrix(NA,2373)
Y_low    <- matrix(NA,2373)
Y_high   <- matrix(NA,2373)

for(f in 1:5){
  Stormdataf <- Stormdata2[fold != f]
  Yf <- Stormdataf$VMAX
  Xf <- Stormdataf |>
    select(-c(
      StormID,
      Date,
      lead_time,
      VMAX,
      basin,
      LON,
      SST,
      CAPE1,
      SHTFL2,
      TCOND7002,
      CP1,
      VMAX_OP_T0
    ))
  Xfscale <- scale(Xf)
  Xfscale <- cbind("Intercept" = rep(1,length(Yf)), Xfscale)
  
  ### Create StormID indicator for random effects model ----
  StormIDs <- Stormdataf |> select(StormID)
  Storms <- unique(StormIDs$StormID)
  IDs <- 1:length(Storms)
  Storms <- data.frame(
    StormID = Storms,
    ID = IDs
  )
  StormIDs <- left_join(StormIDs, Storms)
  StormID <- StormIDs$ID
  
  n <- length(Yf)
  p <- ncol(Xfscale)
  N <- length(IDs)
  
  # Select training data with fold not equal to f
  dataf <- list(Y = Yf, 
                X = Xfscale, 
                n = n, 
                p = p, 
                N = N, 
                StormID = StormID
  )
  
  # Fit model 1   
  mf <- textConnection(
    "model{
  
  # Likelihood
  for(i in 1:n){
    Y[i] ~ dnorm(muY[i], taue)
    muY[i] <- inprod(X[i,], beta[StormID[i],])
  }
  
  # Slopes
  for(j in 1:p){
    for(i in 1:N){
      beta[i,j] ~ dnorm(mu[j], taub[j])
    }
    mu[j] ~ dnorm(0, 0.01)
    taub[j] ~ dgamma(0.1, 0.1)
  }
  
  # Priors
  taue ~ dgamma(0.1, 0.1)
  sigma <- sqrt(1/taue)
  
  # Posterior predictive checks
  for(i in 1:n){
    Y4[i] ~ dnorm(muY[i], taue)
  }
  
   # WAIC calculations
   for(i in 1:n){
     like[i]    <- dnorm(Y[i],muY[i],taue)
   }
  }"
  )
  
  modelf <- jags.model(mf,data = dataf, n.chains=1,quiet=TRUE)
  update(modelf, 10000, progress.bar="none")
  b1     <- coda.samples(modelf, 
                         variable.names=c("beta", "sigma"),
                         thin=5, 
                         n.iter=20000, 
                         progress.bar="none")[[1]]
  b1Sum <- summary(b1)$stat
  b1_beta <- b1Sum[-nrow(b1Sum),]
  b1_sigma <- b1Sum[nrow(b1Sum),]
  b1_BetaMn <- matrix(b1_beta[,1], ncol = p)
  b1_BetaMnSum <- colMeans(b1_BetaMn)
  
  # Make predictions
  for(i in 1:2373){if(fold[i]==f){
    Y_mod1        <- rnorm(nrow(b1), XNscale[i,]%*%b1_BetaMnSum, b1_sigma[1])
    Y_mean[i,1]   <- mean(Y_mod1)
    Y_median[i,1] <- median(Y_mod1)
    Y_low[i,1]    <- quantile(Y_mod1,0.025)
    Y_high[i,1]   <- quantile(Y_mod1,0.975)
    
    # ppd1 <- table(Y_mod1-0.1)
    # 
    # plot(ppd1,main=paste("Observation", i))
    # abline(v=Y[i],lwd=2,col=3) 
    # 
    # legend("topleft",c("PPD Model 1","PPD Model 2","Observed"),lwd=c(1,1,2),col=1:3,bty="n")
  }} 
}

BIAS  <- colMeans(Y_mean-Ynew)
MSE   <- colMeans((Y_mean-Ynew)^2)
MAD   <- colMeans(abs(Y_mean-Ynew))
COV   <- colMeans( (Y_low <= Ynew) & (Ynew <= Y_high))
WIDTH <- colMeans(Y_high-Y_low)

plot(Ynew,Y_mean[,1],pch=19,
     xlim=c(0,100),ylim=c(0,100),
     xlab="Observed",ylab="Predicted")
legend("topleft",c("Model 1"),pch=19,col=1:2,bty="n")

OUT   <- cbind(BIAS,MSE,MAD,COV,WIDTH)
kable(OUT)

fitP <- glm()

Yknown <- data.frame(
  Y = Y,
  Yk = Y_mean[knownOG,]
)
Ymissing <- data.frame(
  Ym = Y_mean[missingOG,]
)

ggplot() +
  geom_histogram(data = Yknown,
                 aes(x = Y, after_stat(density)), 
                 color = "black", fill = "lightblue") +
  geom_density(data = Yknown,
               aes(x = Y), 
               color = "blue", linewidth = 1)  +
  geom_density(data = Yknown,
               aes(x = Yk), 
               color = "green", linewidth = 1)  +
  geom_density(data = Ymissing,
               aes(x = Ym), 
               color = "red", linewidth = 1)  +
  theme_bw()







