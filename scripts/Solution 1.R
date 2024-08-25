#### Solution 1

# Porblem 1 ==========
data <- read.csv("E2_data.csv",header=TRUE,stringsAsFactors = TRUE)
basin <-as.numeric(data$basin)
Y<-data$VMAX
X<-scale(data[,5:22])
library(rjags)
## remove observations with missing Vmax
junk <- is.na(Y)
Yr <- Y[!junk]
Xr <- X[!junk,]
basin_r <- basin[!junk]
n <- length(Yr)
p <- ncol(Xr)
dat <- list(Xr=Xr,Yr=Yr,n=n,p=p,basin_r=basin_r)
model <- "model{
# Likelihood
for(i in 1:n){
Yr[i] ~ dnorm(alpha+inprod(Xr[i,],beta[])+theta[basin_r[i]],taue)
}
# Priors
for(j in 1:p){
beta[j] ~ dnorm(0,0.001)
}
for(j in 1:2){
theta[j]~dnorm(0,tau)
}
alpha ~ dnorm(0,0.01)
tau ~ dgamma(0.1,0.1)
taue ~ dgamma(0.1, 0.1)
}"
library(rjags)
iters <- 5000
chains <- 2
model <- jags.model(textConnection(model),n.chains=chains,data = dat,quiet=TRUE)
update(model, 10000, progress.bar="none")
samp <- coda.samples(model,
                     variable.names=c("beta","alpha","theta"),
                     n.iter=20000, progress.bar="none")
sum <- summary(samp)
gelman.diag(samp)



# Problem 2 ================
## Cross Validation for Model 1 ----
data <- read.csv("E2_data.csv",header=TRUE,stringsAsFactors = TRUE)
basin <-as.numeric(data$basin)
Y<-data$VMAX
X<-scale(data[,5:22])
library(rjags)
## remove observations with missing Vmax
junk <- is.na(Y)
Yr <- Y[!junk]
Xr <- X[!junk,]
## Prepare Data for JAGS
basin_r <- basin[!junk]
n <- length(Yr)
p <- ncol(Xr)
## Break Data into Groups
fold <- rep(1:5,341)
fold <- sample(fold)
length(fold)
## Vectors to store statistics
Y_mean <- c()
Y_median <- c()
Y_low <- c()
Y_high <- c()
error<- c()
abs_error <- c()
## Cross-validation
for(f in 1:5){
  # Select training data with fold not equal to f
  Yp<-Yr[fold==f]
  Xp<-Xr[fold==f,]
  basin_p<-basin_r[fold==f]
  data <- list(Yr=Yr[fold!=f],Xr=Xr[fold!=f,],n=sum(fold!=f),basin_r=basin_r[fold!=f],p=p)
  # Fit model 1
  m1 <- textConnection("model{
    # Likelihood
    for(i in 1:n){
      Yr[i] ~ dnorm(alpha+inprod(Xr[i,],beta[])+theta[basin_r[i]],taue)
    }
    
    # Priors
    for(j in 1:p){
      beta[j] ~ dnorm(0,taue)
    }
    for(j in 1:2){
      theta[j]~dnorm(0,tau)
    }
    alpha ~ dnorm(0,0.01)
    tau ~ dgamma(0.1,0.1)
    taue ~ dgamma(0.1, 0.1)
  }")
  model1 <- jags.model(m1,data = data, n.chains=1,quiet=TRUE)
  update(model1, 10000, progress.bar="none")
  samp <- coda.samples(model1,
                       variable.names=c("alpha","beta","theta","taue"),
                       thin=5, n.iter=20000, progress.bar="none")
  b1 <- samp[[1]]
  ## Make Predictions
  for (i in 1:length(Yr)){if(fold[i]==f){
    Y_mod<- rnorm(nrow(b1),b1[,1]+Xr[i,1]*b1[,2]+Xr[i,2]*b1[,3]+Xr[i,3]*b1[,4]+
                    Xr[i,4]*b1[,5]+
                    Xr[i,5]*b1[,6]+Xr[i,6]*b1[,7]+Xr[i,7]*b1[,8]+
                    Xr[i,8]*b1[,9]+Xr[i,9]*b1[,10]+
                    Xr[i,10]*b1[,11]+Xr[i,10]*b1[,11]+Xr[i,11]*b1[,12]+
                    Xr[i,12]*b1[,13]+Xr[i,13]*b1[,14]+Xr[i,14]*b1[,15]+
                    Xr[i,15]*b1[,16]+Xr[i,16]*b1[,17]+Xr[i,17]*b1[,18]+Xr[i,18]*b1[,19]+
                    b1[,20+basin_r[i]],1/b1[,20])
    Y_mean[i] <- mean(Y_mod)
    Y_low[i] <- quantile(Y_mod,0.025)
    Y_high[i] <- quantile(Y_mod,0.975)
    error[i]<-abs(Y_mean[i]-Yr[i])
  }}
  rm(model1)
}
y <- cbind(Yr,Yr) # Make data the same size/format as predictions
BIAS <- colMeans(Y_mean-y)
MSE <- colMeans((Y_mean-y)^2)
MAD <- colMeans(abs(Y_mean-y))
COV <- colMeans( (Y_low <= y) & (y <= Y_high))

## Cross Validation for Model 2 ----
data <- read.csv("E2_data.csv",header=TRUE,stringsAsFactors = TRUE)
basin <-as.numeric(data$basin)
Y<-data$VMAX
X<-scale(data[,5:22])
library(rjags)
## remove observations with missing Vmax
junk <- is.na(Y)
Yr <- Y[!junk]
Xr <- X[!junk,]
## Prepare Data for JAGS
n <- length(Yr)
p <- ncol(Xr)
## Split data into groups
fold <- rep(1:5,341)
fold <- sample(fold)
length(fold)
## Vectors to store statistics
Y_mean <- c()
Y_median <- c()
Y_low <- c()
Y_high <- c()
error<- c()
abs_error <- c()
for(f in 1:5){
  # Select training data with fold not equal to f
  Yp<-Yr[fold==f]
  Xp<-Xr[fold==f,]
  basin_p<-basin_r[fold==f]
  data <- list(Yr=Yr[fold!=f],Xr=Xr[fold!=f,],n=sum(fold!=f),p=p)
  # Fit model 1
  m1 <- textConnection("model{
    # Likelihood
    for(i in 1:n){
      Yr[i] ~ dnorm(alpha+inprod(Xr[i,],beta[]),taue)
    }
    # Priors
    for(j in 1:p){
      beta[j] ~ dnorm(0,0.001)
    }
    alpha ~ dnorm(0,0.01)
    taue ~ dgamma(0.1, 0.1)
  }")
  model1 <- jags.model(m1,data = data, n.chains=1,quiet=TRUE)
  update(model1, 10000, progress.bar="none")
  samp <- coda.samples(model1,
                       variable.names=c("alpha","beta","taue"),
                       thin=5, n.iter=20000, progress.bar="none")
  b1 <- samp[[1]]
  ## Make Predictions
  for (i in 1:length(Yr)){if(fold[i]==f){
    Y_mod<- rnorm(nrow(b1),b1[,1]+Xr[i,1]*b1[,2]+Xr[i,2]*b1[,3]+Xr[i,3]*b1[,4]+Xr[i,4]*b1[,5]+
                    Xr[i,5]*b1[,6]+
                    Xr[i,6]*b1[,7]+Xr[i,7]*b1[,8]+Xr[i,8]*b1[,9]+
                    Xr[i,9]*b1[,10]+Xr[i,10]*b1[,11]+Xr[i,10]*b1[,11]+
                    Xr[i,11]*b1[,12]+
                    Xr[i,12]*b1[,13]+Xr[i,13]*b1[,14]+Xr[i,14]*b1[,15]+
                    Xr[i,15]*b1[,16]+Xr[i,16]*b1[,17]+Xr[i,17]*b1[,18]+
                    Xr[i,18]*b1[,19],
                  1/b1[,20])
    Y_mean[i] <- mean(Y_mod)
    Y_low[i] <- quantile(Y_mod,0.025)
    Y_high[i] <- quantile(Y_mod,0.975)
    error[i]<-abs(Y_mean[i]-Yr[i])
  }}
  rm(model1)
}
y <- cbind(Yr,Yr) # Make data the same size/format as predictions
BIAS <- colMeans(Y_mean-y)
MSE <- colMeans((Y_mean-y)^2)
MAD <- colMeans(abs(Y_mean-y))
COV <- colMeans( (Y_low <= y) & (y <= Y_high))


# Problem 3 ===================
data <- read.csv("E2_data.csv",header=TRUE,stringsAsFactors = TRUE)
basin <-as.numeric(data$basin)
Y<-data$VMAX
X<-scale(data[,5:22])
library(rjags)
## remove observations with missing Vmax
junk <- is.na(Y)
Yr <- Y[!junk]
Xr <- X[!junk,]
basin_r <- basin[!junk]
n <- length(Yr)
p <- ncol(Xr)
dat <- list(Xr=Xr,Yr=Yr,n=n,p=p,basin_r=basin_r)
model <- "model{
# Likelihood
for(i in 1:n){
  Yr[i] ~ dnorm(alpha+inprod(Xr[i,],beta[])+theta[basin_r[i]],taue)
}

# Priors
for(j in 1:p){
  beta[j] ~ dnorm(0,0.001)
}
for(j in 1:2){
  theta[j]~dnorm(0,tau)
}
alpha ~ dnorm(0,0.01)
tau ~ dgamma(0.1,0.1)
taue ~ dgamma(0.1, 0.1)

# Posterior preditive checks
for(i in 1:n){
  Y2[i] ~ dnorm(alpha+inprod(Xr[i,],beta[])+theta[basin_r[i]],taue)
}
D[1] <- max(Y2[])-min(Y2[])
}"
library(rjags)
iters <- 5000
chains <- 2
model <- jags.model(textConnection(model),
                    n.chains=chains,data = dat,quiet=TRUE)
update(model, 10000, progress.bar="none")
samp <- coda.samples(model,
                     variable.names=c("D","beta"),
                     n.iter=20000, progress.bar="none")
sum <- summary(samp)
D1 <- samp[[1]]
# Compute the test stats for the data
D0 <- c(max(Yr)-min(Yr))
Dnames <- c("Range Y")
# Compute the test stats for the models
pval1 <- 0
names(pval1)<-Dnames
plot(density(D1[,1]), xlim=range(c(D0[1],D1[,1])),
     xlab="D",ylab="Posterior probability",
     main=Dnames[1])
abline(v=D0[1],col=2)
legend("topleft",c("Model","Data"),lty=1,col=1:2,bty="n")
pval1[1] <- mean(D1[,1]>D0[1])


# Problem 4 ====================
data <- read.csv("E2_data.csv",header=TRUE,stringsAsFactors = TRUE)
basin <-as.numeric(data$basin)
Y<-data$VMAX
X<-scale(data[,5:22])
library(rjags)
## remove observations with missing Vmax
junk <- is.na(Y)
Yr <- Y[!junk]
Xr <- X[!junk,]
basin_r <- basin[!junk]
n <- length(Yr)
p <- ncol(Xr)
dat <- list(Xr=Xr,Yr=Yr,n=n,p=p,basin_r=basin_r)
model <- "model{
# Likelihood
for(i in 1:n){
Yr[i] ~ dnorm(alpha+inprod(Xr[i,],beta[])+
theta[basin_r[i]],taue)
}
# Priors
for(j in 1:p){
beta[j] <- gamma[j]*delta[j]
gamma[j] ~ dbern(0.5)
delta[j] ~ dnorm(0,tau)
}
for(j in 1:2){
theta[j]~dnorm(0,tau)
}
alpha ~ dnorm(0,0.01)
tau ~ dgamma(0.1,0.1)
taue ~ dgamma(0.1, 0.1)
}"
library(rjags)
iters <- 5000
chains <- 2
model <- jags.model(textConnection(model),n.chains=chains,data = dat,quiet=TRUE)
update(model, 10000, progress.bar="none")
samp <- coda.samples(model,
                     variable.names=c("beta"),
                     n.iter=20000, progress.bar="none")
sum <- summary(samp)
beta <- NULL
for(l in 1:chains){
  beta <- rbind(beta,samp[[l]])
}
for(j in 1:p){
  hist(beta[,j],xlab=expression(beta[j]),ylab="Posterior density",
       breaks=100)
}
Inc_Prob <- apply(beta!=0,2,mean)
Q <- t(apply(beta,2,quantile,c(0.5,0.05,0.95)))
out <- cbind(Inc_Prob,Q)


# Problem 5 ====================
## (code for cross-validation is under Problem 3)
data <- read.csv("E2_data.csv",header=TRUE,stringsAsFactors = TRUE)
basin <-as.numeric(data$basin)
Y<-data$VMAX
X<-scale(data[,5:22])
nt<-length(Y)
library(rjags)
## remove observations with missing Vmax
junk <- is.na(Y)
Yr <- Y[!junk]
Xr <- X[!junk,]
basin_r <- basin[!junk]
Xp <-X[junk,]
Yp <-Y[junk]
basin_p<-basin[junk]
n <- length(Yr)
p <- ncol(Xr)
np <- length(Yp)
dat <- list(Xr=Xr,Yr=Yr,n=n,p=p,np=np,Xp=Xp,basin_r=basin_r,basin_p=basin_p)
model <- "model{
# Likelihood
for(i in 1:n){
Yr[i] ~ dnorm(alpha+inprod(Xr[i,],beta[])+
theta[basin_r[i]],taue)
}
# Priors
for(j in 1:p){
beta[j] ~ dnorm(0,0.001)
}
for(j in 1:2){
theta[j]~dnorm(0,tau)
}
alpha ~ dnorm(0,0.01)
tau ~ dgamma(0.1,0.1)
taue ~ dgamma(0.1, 0.1)
# Prediction
for(i in 1:np){
Yp[i] ~ dnorm(alpha+inprod(Xp[i,],beta[])+
theta[basin_p[i]],taue)
}
}"
library(rjags)
iters <- 5000
chains <- 2
model <- jags.model(textConnection(model),n.chains=chains,data = dat,quiet=TRUE)
update(model, 10000, progress.bar="none")
samp <- coda.samples(model,
                     variable.names=c("Yp","beta"),
                     n.iter=20000, progress.bar="none")
sum <- summary(samp)
D1 <- samp[[1]]
Yp.samps <- D1[,1:np]
Y.predict <-colMeans(Yp.samps)
low_v <- apply(Yp.samps,2,quantile,0.025)
high_v <- apply(Yp.samps,2,quantile,0.975)
L<-rep(NA,nt)
U<-rep(NA,nt)
j<-1
for(i in 1:nt){if (junk[i]){
  Y[i]<-Y.predict[j]
  L[i]<-low_v[j]
  U[i]<-high_v[j]
  j<-j+1
}
}
predict<-cbind(data$HWRF,Y,L,U)
write.csv(predict,"EbbyLouisa.csv")

load("PollardTyler.csv")
colnames(outputDF) <- c("HWRF", "VMAXtlyer", "Ltyler", "Utyler")
combPred_df <- cbind(predict, outputDF)
combPred_df <- combPred_df |>
  select(-1) |>
  select(HWRF, everything()) 

ActualY <- read.csv("Actual Y.csv")

TotalY <- cbind(combPred_df, ActualY)
TotalYpreds <- TotalY |> filter(!is.na(x))

MAD1 <- mean(abs(TotalYpreds$Y - TotalYpreds$x))
COV1 <- mean( (TotalYpreds$L <= TotalYpreds$x) & (TotalYpreds$x <= TotalYpreds$U))

MADtyler <- mean(abs(TotalYpreds$VMAXtlyer - TotalYpreds$x))
COVtyler <- mean( (TotalYpreds$Ltyler <= TotalYpreds$x) & (TotalYpreds$x <= TotalYpreds$Utyler))

MAD1
MADtyler

COV1
COVtyler

