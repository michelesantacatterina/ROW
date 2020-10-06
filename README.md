# ROW: an R package for obtaining Robust Orthoghonality Weights - a set of weights that optimally balance confounders for estimating the effect of binary and continuous treatments

This repository contains the R code to obtain ROW. Citation.

## How to install the package

1. Install devtools from CRAN

```
install.packages("devtools")
```
2. Load the devtools package

```
library(devtools)
```
3. Install ROW package

```
install_github("michelesantacatterina/ROW")
```

## What solver should you use?

ROW are obtained by solving a quadratic constrained optimization problem. Several solvers can be used to solve the optimization problem. In this package we provide the following options: Gurobi, quadprog, Dykstra and ipop. We suggest to use Gurobi. Please find detailed instruction on how to install the R interface of Gurobi [here](https://cran.r-project.org/web/packages/prioritizr/vignettes/gurobi_installation.html). Gurobi offers free academic licenses ([visit the solver website](https://www.gurobi.com/downloads/end-user-license-agreement-academic/)). The following Figure shows the computational time (computed using the system.time function in R) for solving the same optimization problem with n=50, 200, 1,000 and 5,000 across the four solvers. 

![](images/qp_solvers.png?raw=true)

## Examples

We provide some examples on how to estimate hazard ratios and ATE of binary and continuous treatments.

### Hazard Ratio - Binary Treatment

```
set.seed(4)

########################################################################
# Generate Covariates/Counfounders
########################################################################

n  <- 1000
X1 <- rnorm(n,.1,1)
X2 <- rnorm(n,.1,1)
X3 <- rlnorm(n,0,.5)
X4 <- 5*rbeta(n,3,1)

X5 <- sample( 1:4, n, replace=TRUE, prob=c(0.35, 0.25, 0.05, 0.35) )
X6 <- sample( 0:1, n, replace=TRUE, prob=c(0.75, 0.25) )

X <- cbind(X1,X2,X3,X4,X5,X6)


########################################################################
# Time-to-event data - Survival analysis
########################################################################

lambda  <- 0.01
beta_ou <- 1.4
beta_tr <- 1.5
d       <- 0.01
rho     <- 1
true_haz<- 0.2 #true log hazard ratio, exp(hr)= 1.22
rateC   <- 0.1 #rate of censoring

# --------------------------------------------------------------------------------
# ---------- Binary treatment

mbeta   <- mean(beta_tr*(X[,1] + X[,2] + X[,3] + X[,4] + X[,5]) + X[,6])
prt     <- 1./(1. + exp( mbeta - beta_tr * (X[,1] + X[,2] + X[,3] + X[,4] + X[,5]) + X[,6]))
Tr      <- rbinom(n,1,prt)
v       <- runif(n=n)
Tlat    <- (- log(v) / (lambda * exp(Tr*d + X[,2] + beta_ou*(X[,4] + X[,5] ) + X[,6] )))^(1 / rho)

# Censoring times
C       <- rexp(n=n, rate=rateC)

# Follow-up times and event indicators
time      <- pmin(Tlat, C)
status    <- as.numeric(Tlat <= C)
censoring <- 1-status


########## Compute ROW

row     <- ROW(intervention = Tr, confounders = X, delta = 0.0001, solver = 'gurobi')
summary(row$weights*n) # the weights sum up to 1


########## Checking covariate balance (absolute standardized mean difference)

library(cobalt)
lprow     <- love.plot(bal.tab(Tr~X,weights=row$weights,method="weighting",
                              s.d.denom="pooled"),colors = c("grey", "black"),var.order = "unadjusted",
                              shapes=c("circle","square"),stars = "raw", position=1, abs=TRUE)

lprow

abal_row  <- abs(bal.tab(Tr~X,weights=row$weights,method="weighting")$Balance$Diff.Adj)
c(min(abal_row), median(abal_row), max(abal_row))


########## Estimating hazard ratio
library(survey)
data          <- data.frame(Tr,X,time,status)
surveyDesign  <- svydesign(ids=1:n,weights=~row$weights,data=data)
modelw        <- svycoxph(Surv(time, status) ~ Tr,design=surveyDesign)

summary(modelw)

# Compare results with naive model
summary(coxph(Surv(time, status) ~ Tr, data=data))
```

### Hazard Ratio - Continuous Treatment

```
########################################################################
# Generate Covariates/Counfounders
########################################################################

n  <- 1000
X1 <- rnorm(n,.1,1)
X2 <- rnorm(n,.1,1)
X3 <- rlnorm(n,0,.5)
X4 <- 5*rbeta(n,3,1)

X5 <- sample( 1:4, n, replace=TRUE, prob=c(0.35, 0.25, 0.05, 0.35) )
X6 <- sample( 0:1, n, replace=TRUE, prob=c(0.75, 0.25) )

X <- cbind(X1,X2,X3,X4,X5,X6)


########################################################################
# Time-to-event data - Survival analysis
########################################################################

lambda  <- 0.01
beta_ou <- 1.4
beta_tr <- 1.5
d       <- 0.01
rho     <- 1
true_haz<- 0.2 #true log hazard ratio, exp(hr)= 1.22
rateC   <- 0.1 #rate of censoring
# --------------------------------------------------------------------------------
# ---------- Continuous treatment

temp  <- 0.1
mTr   <- mean((X[,1] + X[,2] + X[,3] + X[,4] + X[,5] + X[,6]) )
Tr    <- - mTr + (X[,1] + X[,2] + X[,3] + X[,4] + X[,5] + X[,6]) + rnorm(n)
v     <- runif(n=n)
Tlat  <- (- log(v) / (lambda * exp(Tr*d + X[,2] + beta_ou*(X[,4] + X[,5] ) + X[,6] )))^(1 / rho)

# Censoring times
C       <- rexp(n=n, rate=rateC)

# Follow-up times and event indicators
time      <- pmin(Tlat, C)
status    <- as.numeric(Tlat <= C)
censoring <- 1-status


########## Compute ROW

row     <- ROW(intervention = Tr, confounders = X, delta = 0.0001, solver = 'gurobi')
summary(row$weights*n) # the weights sum up to 1


########## Checking covariate balance (absolute treatment-covariate correlation)

library(cobalt)
lprow     <- love.plot(bal.tab(Tr~X,weights=row$weights,method="weighting",
                              s.d.denom="all"),colors = c("grey", "black"),var.order = "unadjusted",
                              shapes=c("circle","square"),stars = "raw", position=1, abs=TRUE)

lprow

abal_row  <- abs(bal.tab(Tr~X,weights=row$weights,method="weighting")$Balance$Corr.Adj)
c(min(abal_row), median(abal_row), max(abal_row))


########## Estimating hazard ratio

library(survey)
data          <- data.frame(Tr,X,time,status)
surveyDesign  <- svydesign(ids=1:n,weights=~row$weights,data=data)
modelw        <- svycoxph(Surv(time, status) ~ Tr,design=surveyDesign)

summary(modelw)

# Compare results with naive model
summary(coxph(Surv(time, status) ~ Tr, data=data))
```


### Average Treatment Effect - Binary Treatment

```
########################################################################
# Generate Covariates/Counfounders
########################################################################

n  <- 1000
X1 <- rnorm(n,.1,1)
X2 <- rnorm(n,.1,1)
X3 <- rlnorm(n,0,.5)
X4 <- 5*rbeta(n,3,1)

X5 <- sample( 1:4, n, replace=TRUE, prob=c(0.35, 0.25, 0.05, 0.35) )
X6 <- sample( 0:1, n, replace=TRUE, prob=c(0.75, 0.25) )

X <- cbind(X1,X2,X3,X4,X5,X6)

#########################################################################
# Continuous Outcome - Average Treatment Effect (ATE)
##########################################################################

beta_ou <- 1.4
beta_tr <- 1.5

# --------------------------------------------------------------------------------
# ---------- Binary treatment

mbeta   <- mean(beta_tr*(X[,1] + X[,2] + X[,3] + X[,4] + X[,5]) + X[,6])
prt     <- 1./(1. + exp( mbeta - beta_tr * (X[,1] + X[,2] + X[,3] + X[,4] + X[,5]) + X[,6]))
Tr      <- rbinom(n,1,prt)
t_effect<- 2 #true effect for binary treatment
Y0      <- 1  + (X[,1] + X[,2] + X[,3] + X[,4] + X[,5] + X[,6]) + rnorm(n)
Y1      <- Y0 + t_effect
Y       <- Y0*(1-Tr) + Y1*Tr

########## Compute ROW

row     <- ROW(intervention = Tr, confounders = X, delta = 0.0001, solver = 'gurobi')
summary(row$weights*n) # the weights sum up to 1


########## Checking covariate balance (absolute standardized mean difference)

library(cobalt)
lprow     <- love.plot(bal.tab(Tr~X,weights=row$weights,method="weighting",
                              s.d.denom="pooled"),colors = c("grey", "black"),var.order = "unadjusted",
                              shapes=c("circle","square"),stars = "raw", position=1, abs=TRUE)

lprow

abal_row  <- abs(bal.tab(Tr~X,weights=row$weights,method="weighting")$Balance$Diff.Adj)
c(min(abal_row), median(abal_row), max(abal_row))


########## Estimating ATE (true effect = 2)

library(sandwich)
data      <- data.frame(Y,Tr)
fit       <- lm(Y~Tr,weights=row$weights,data=data)
fit$coef[2]
sqrt(diag(sandwich(fit)))[2]

# Compare with Naive estimator
fitn      <- lm(Y~Tr,data=data)
fitn$coef[2]
sqrt(diag(sandwich(fitn)))[2]
```

### Effect of a continuous treatment

```
########################################################################
# Generate Covariates/Counfounders
########################################################################

n  <- 1000
X1 <- rnorm(n,.1,1)
X2 <- rnorm(n,.1,1)
X3 <- rlnorm(n,0,.5)
X4 <- 5*rbeta(n,3,1)

X5 <- sample( 1:4, n, replace=TRUE, prob=c(0.35, 0.25, 0.05, 0.35) )
X6 <- sample( 0:1, n, replace=TRUE, prob=c(0.75, 0.25) )

X <- cbind(X1,X2,X3,X4,X5,X6)

#########################################################################
# Continuous Outcome - Average Treatment Effect (ATE)
##########################################################################

beta_ou <- 1.4
beta_tr <- 1.5

# --------------------------------------------------------------------------------
# ---------- Continuous treatment

temp  <- 0.1
t_effect<- 1 #true effect for continuous treatment
mTr   <- mean((X[,1] + X[,2] + X[,3] + X[,4] + X[,5] + X[,6]) )
Tr    <- - mTr + (X[,1] + X[,2] + X[,3] + X[,4] + X[,5] + X[,6]) + rnorm(n)
Y     <- t_effect*Tr + (X[,1] + X[,2] + X[,3] + X[,4] + X[,5] + X[,6]) + rnorm(n)


########## Compute ROW

row     <- ROW(intervention = Tr, confounders = X, delta = 0.0001, solver = 'gurobi')
summary(row$weights*n) # the weights sum up to 1

########## Checking covariate balance (absolute treatment-covariate correlation)

library(cobalt)
lprow     <- love.plot(bal.tab(Tr~X,weights=row$weights,method="weighting",
                              s.d.denom="all"),colors = c("grey", "black"),var.order = "unadjusted",
                              shapes=c("circle","square"),stars = "raw", position=1, abs=TRUE)

lprow

abal_row  <- abs(bal.tab(Tr~X,weights=row$weights,method="weighting")$Balance$Corr.Adj)
c(min(abal_row), median(abal_row), max(abal_row))


########## Estimating effect of continuous treatment (true effect = 1)

library(sandwich)
data      <- data.frame(Y,Tr)
fit       <- lm(Y~Tr,weights=row$weights,data=data)
fit$coef[2]
sqrt(diag(sandwich(fit)))[2]

# Compare results with naive estimator
fitn      <- lm(Y~Tr,data=data)
fitn$coef[2]
sqrt(diag(sandwich(fitn)))[2]

```