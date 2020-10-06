#' Robust Orthogonality Weights for balancing covariates for estimating the effect of binary and continuous treatments
#'
#' \code{ROW} estimates weights with minimal variance, thus controlling for extreme weights, that optimally balance covariates by satisfing a constraints on the sample correlation between
#' covariates and treatment. \code{ROW} can be used to estimate the effects of binary and continuous treatments on continuous outcomes (ATE), binary outcomes (odds ratio) and
#' time-to-event outcomes (hazard ratio). \code{ROW} is based on a constrained quadratic optimization problem.
#' While also \code{quadprog} can be used to obtain the weights, we suggest to use \code{gurobi}.
#' To install the \code{R} interface of \code{Gurobi} please follow instructions
#' \href{https://cran.r-project.org/web/packages/prioritizr/vignettes/gurobi_installation.html}{here}. \code{gurobi} provides
#' \href{https://www.gurobi.com/downloads/end-user-license-agreement-academic/}{free academic licenses}
#' and \href{https://www.gurobi.com/downloads/eula-academic-online-course-license-request/}{free academic online course licenses}.
#'
#' @param intervention binary or continuous intervention/treatment/exposure.
#' @param confounders a matrix containing all confounders. Each column represents a confounder.
#' @param delta upper bound for the absolute value of the sample correlation between treatment and confounders. Default \code{delta=0.01}
#' @param tol lower bound for the weights, \i.e., weights must be positive. Default \code{tol=1e-08}
#' @param solver the solver used to obtain the weights. Choices are between \code{gurobi} (default), \code{quadprog}, \code{Dykstra}, and \code{ipop}.
#' @param Presolve presolve parameter for \code{gurobi} optimizer. More info \href{https://www.gurobi.com/documentation/8.1/refman/presolve.html}{here}. Default \code{Presolve=2} - aggressive.
#' @param OutputFlag enables or disables \code{gurobi} solver output. More info \href{https://www.gurobi.com/documentation/8.1/refman/outoutflag.html}{here}. Default \code{OutputFlag=0}.
#' @return A list the following.
#' \item{weights}{a vector of size n containing the set of Robust Ortoghonality Weights}
#' \item{res}{\code{gurobi} or \code{quadprog} type object}
#' @author Michele Santacatterina, \email{msanta@gwu.edu}
#' @references Santacatterina, M., Robust weights that optimally balance confounders for estimating the effect of binary and continuous treatments with time-to-event data, url: TBA
#' @examples
#' \dontrun{
#'set.seed(4)
#'
#'########################################################################
#'# Generate Covariates/Counfounders
#'########################################################################
#'
#'n  <- 1000
#'X1 <- rnorm(n,.1,1)
#'X2 <- rnorm(n,.1,1)
#'X3 <- rlnorm(n,0,.5)
#'X4 <- 5*rbeta(n,3,1)
#'
#'X5 <- sample( 1:4, n, replace=TRUE, prob=c(0.35, 0.25, 0.05, 0.35) )
#'X6 <- sample( 0:1, n, replace=TRUE, prob=c(0.75, 0.25) )
#'
#'X <- cbind(X1,X2,X3,X4,X5,X6)
#'
#'
#'########################################################################
#'# Time-to-event data - Survival analysis
#'########################################################################
#'
#'lambda  <- 0.01
#'beta_ou <- 1.4
#'beta_tr <- 1.5
#'d       <- 0.01
#'rho     <- 1
#'true_haz<- 0.2 #true log hazard ratio, exp(hr)= 1.22
#'rateC   <- 0.1 #rate of censoring
#'
#'# --------------------------------------------------------------------------------
#'# ---------- Binary treatment
#'
#'mbeta   <- mean(beta_tr*(X[,1] + X[,2] + X[,3] + X[,4] + X[,5]) + X[,6])
#'prt     <- 1./(1. + exp( mbeta - beta_tr * (X[,1] + X[,2] + X[,3] + X[,4] + X[,5]) + X[,6]))
#'Tr      <- rbinom(n,1,prt)
#'v       <- runif(n=n)
#'Tlat    <- (- log(v) / (lambda * exp(Tr*d + X[,2] + beta_ou*(X[,4] + X[,5] ) + X[,6] )))^(1 / rho)
#'
#'# Censoring times
#'C       <- rexp(n=n, rate=rateC)
#'
#'# Follow-up times and event indicators
#'time      <- pmin(Tlat, C)
#'status    <- as.numeric(Tlat <= C)
#'censoring <- 1-status
#'
#'
#'########## Compute ROW
#'
#'row     <- ROW(intervention = Tr, confounders = X, delta = 0.0001, solver = 'gurobi')
#'summary(row$weights*n) # the weights sum up to 1
#'
#'
#'########## Checking covariate balance (absolute standardized mean difference)
#'
#'library(cobalt)
#'lprow     <- love.plot(bal.tab(Tr~X,weights=row$weights,method="weighting",
#'                               s.d.denom="pooled"),colors = c("grey", "black"),var.order = "unadjusted",
#'                               shapes=c("circle","square"),stars = "raw", position=1, abs=TRUE)
#'
#'lprow
#'
#'abal_row  <- abs(bal.tab(Tr~X,weights=row$weights,method="weighting")$Balance$Diff.Adj)
#'c(min(abal_row), median(abal_row), max(abal_row))
#'
#'
#'########## Estimating hazard ratio
#'library(survey)
#'data          <- data.frame(Tr,X,time,status)
#'surveyDesign  <- svydesign(ids=1:n,weights=~row$weights,data=data)
#'modelw        <- svycoxph(Surv(time, status) ~ Tr,design=surveyDesign)
#'
#'summary(modelw)
#'
#'# Compare results with naive model
#'summary(coxph(Surv(time, status) ~ Tr, data=data))
#'
#'
#'
#'
#'# --------------------------------------------------------------------------------
#'# ---------- Continuous treatment
#'
#'temp  <- 0.1
#'mTr   <- mean((X[,1] + X[,2] + X[,3] + X[,4] + X[,5] + X[,6]) )
#'Tr    <- - mTr + (X[,1] + X[,2] + X[,3] + X[,4] + X[,5] + X[,6]) + rnorm(n)
#'v     <- runif(n=n)
#'Tlat  <- (- log(v) / (lambda * exp(Tr*d + X[,2] + beta_ou*(X[,4] + X[,5] ) + X[,6] )))^(1 / rho)
#'
#'# Censoring times
#'C       <- rexp(n=n, rate=rateC)
#'
#'# Follow-up times and event indicators
#'time      <- pmin(Tlat, C)
#'status    <- as.numeric(Tlat <= C)
#'censoring <- 1-status
#'
#'
#'########## Compute ROW
#'
#'row     <- ROW(intervention = Tr, confounders = X, delta = 0.0001, solver = 'gurobi')
#'summary(row$weights*n) # the weights sum up to 1
#'
#'
#'########## Checking covariate balance (absolute treatment-covariate correlation)
#'
#'library(cobalt)
#'lprow     <- love.plot(bal.tab(Tr~X,weights=row$weights,method="weighting",
#'                               s.d.denom="all"),colors = c("grey", "black"),var.order = "unadjusted",
#'                               shapes=c("circle","square"),stars = "raw", position=1, abs=TRUE)
#'
#'lprow
#'
#'abal_row  <- abs(bal.tab(Tr~X,weights=row$weights,method="weighting")$Balance$Corr.Adj)
#'c(min(abal_row), median(abal_row), max(abal_row))
#'
#'
#'########## Estimating hazard ratio
#'
#'library(survey)
#'data          <- data.frame(Tr,X,time,status)
#'surveyDesign  <- svydesign(ids=1:n,weights=~row$weights,data=data)
#'modelw        <- svycoxph(Surv(time, status) ~ Tr,design=surveyDesign)
#'
#'summary(modelw)
#'
#'# Compare results with naive model
#'summary(coxph(Surv(time, status) ~ Tr, data=data))
#'
#'##########################################################################
#'
#'
#'#########################################################################
#'# Continuous Outcome - Average Treatment Effect (ATE)
#'##########################################################################
#'
#'beta_ou <- 1.4
#'beta_tr <- 1.5
#'
#'# --------------------------------------------------------------------------------
#'# ---------- Binary treatment
#'
#'mbeta   <- mean(beta_tr*(X[,1] + X[,2] + X[,3] + X[,4] + X[,5]) + X[,6])
#'prt     <- 1./(1. + exp( mbeta - beta_tr * (X[,1] + X[,2] + X[,3] + X[,4] + X[,5]) + X[,6]))
#'Tr      <- rbinom(n,1,prt)
#'t_effect<- 2 #true effect for binary treatment
#'Y0      <- 1  + (X[,1] + X[,2] + X[,3] + X[,4] + X[,5] + X[,6]) + rnorm(n)
#'Y1      <- Y0 + t_effect
#'Y       <- Y0*(1-Tr) + Y1*Tr
#'
#'########## Compute ROW
#'
#'row     <- ROW(intervention = Tr, confounders = X, delta = 0.0001, solver = 'gurobi')
#'summary(row$weights*n) # the weights sum up to 1
#'
#'
#'########## Checking covariate balance (absolute standardized mean difference)
#'
#'library(cobalt)
#'lprow     <- love.plot(bal.tab(Tr~X,weights=row$weights,method="weighting",
#'                               s.d.denom="pooled"),colors = c("grey", "black"),var.order = "unadjusted",
#'                               shapes=c("circle","square"),stars = "raw", position=1, abs=TRUE)
#'
#'lprow
#'
#'abal_row  <- abs(bal.tab(Tr~X,weights=row$weights,method="weighting")$Balance$Diff.Adj)
#'c(min(abal_row), median(abal_row), max(abal_row))
#'
#'
#'########## Estimating ATE (true effect = 2)
#'
#'library(sandwich)
#'data      <- data.frame(Y,Tr)
#'fit       <- lm(Y~Tr,weights=row$weights,data=data)
#'fit$coef[2]
#'sqrt(diag(sandwich(fit)))[2]
#'
#'# Compare with Naive estimator
#'fitn      <- lm(Y~Tr,data=data)
#'fitn$coef[2]
#'sqrt(diag(sandwich(fitn)))[2]
#'
#'# --------------------------------------------------------------------------------
#'# ---------- Continuous treatment
#'
#'temp  <- 0.1
#'t_effect<- 1 #true effect for continuous treatment
#'mTr   <- mean((X[,1] + X[,2] + X[,3] + X[,4] + X[,5] + X[,6]) )
#'Tr    <- - mTr + (X[,1] + X[,2] + X[,3] + X[,4] + X[,5] + X[,6]) + rnorm(n)
#'Y     <- t_effect*Tr + (X[,1] + X[,2] + X[,3] + X[,4] + X[,5] + X[,6]) + rnorm(n)
#'
#'
#'########## Compute ROW
#'
#'row     <- ROW(intervention = Tr, confounders = X, delta = 0.0001, solver = 'gurobi')
#'summary(row$weights*n) # the weights sum up to 1
#'
#'########## Checking covariate balance (absolute treatment-covariate correlation)
#'
#'library(cobalt)
#'lprow     <- love.plot(bal.tab(Tr~X,weights=row$weights,method="weighting",
#'                               s.d.denom="all"),colors = c("grey", "black"),var.order = "unadjusted",
#'                               shapes=c("circle","square"),stars = "raw", position=1, abs=TRUE)
#'
#'lprow
#'
#'abal_row  <- abs(bal.tab(Tr~X,weights=row$weights,method="weighting")$Balance$Corr.Adj)
#'c(min(abal_row), median(abal_row), max(abal_row))
#'
#'
#'########## Estimating effect of continuous treatment (true effect = 1)
#'
#'library(sandwich)
#'data      <- data.frame(Y,Tr)
#'fit       <- lm(Y~Tr,weights=row$weights,data=data)
#'fit$coef[2]
#'sqrt(diag(sandwich(fit)))[2]
#'
#'# Compare results with naive estimator
#'fitn      <- lm(Y~Tr,data=data)
#'fitn$coef[2]
#'sqrt(diag(sandwich(fitn)))[2]
#'
#' }
#'
#devtools::load_all()

#library(gurobi)

ROW <- function(intervention,
                confounders,
                delta=0.01,
                tol=1e-08,
                solver='gurobi',
                Presolve=2,
                OutputFlag=0){

  Xs <- scale(cbind(confounders))
  trs <- as.numeric(scale(intervention))
  if(is.null(dim(confounders)[1])!=T){
    n <- dim(confounders)[1]
  }
  else{
    n <- length(confounders)
  }#end else


  reptimes          <- dim(Xs)[2]
  tol               <- tol

  #min(1/2 b^T D b -d^T b, A^T b >= b_0)
  if(solver == 'quadprog' ){
    Dmat  <- diag(n)
    Amat  <- t(rbind(matrix(c(rep(1,n),Xs*trs,-Xs*trs), nrow = (reptimes*2+1), byrow = T), diag(n) ))
    bvec  <- c(1,rep(delta,reptimes),rep(-delta,reptimes),rep(tol, n))
    dvec  <- rep(1/n,n)
    meq   <- 1
    res   <- quadprog::solve.QP(Dmat,dvec,Amat,bvec,meq=meq)

    weights <- res$solution
  }#end if solver == quadprog

  if(solver == 'Dykstra' ){
    Dmat  <- diag(n)
    Amat  <- t(rbind(matrix(c(rep(1,n),Xs*trs,-Xs*trs), nrow = (reptimes*2+1), byrow = T), diag(n) ))
    bvec  <- c(1,rep(delta,reptimes),rep(-delta,reptimes),rep(tol, n))
    dvec  <- rep(1/n,n)
    meq   <- 1
    res   <- Dykstra::dykstra(Dmat,dvec,Amat,bvec,meq=meq)

    weights <- res$solution
  }#end if solver == Dykstra

  if(solver == 'ipop' ){
    H <- diag(n)
    A <- matrix(c(rep(Xs*trs,1),rep(1,n)), nrow = (reptimes+1), byrow = T)
    b <- c(rep(-delta,reptimes),1)
    r <- c(rep(2*delta,reptimes),0)
    l <- matrix(rep(tol, n))
    u <- matrix(rep(n*100, n))
    c <- matrix(rep(1/n,n))
    res <- kernlab::ipop(c, H, A, b, l, u, r)

    weights <- primal(res)
  }#end if solver == 'ipop'


  if(solver == 'gurobi' ){
    model             <- list()
    model$A           <- matrix(c(rep(Xs*trs,2),rep(1,n)), nrow = (reptimes*2+1), byrow = T)
    model$rhs         <- c(rep(delta,reptimes),rep(-delta,reptimes),1)
    model$modelsense  <- "min"
    model$Q           <- diag(n)
    model$obj         <- rep(1/n,n)
    model$sense       <- c(rep("<=",reptimes),rep(">=",reptimes),"=")
    model$lb          <- rep(tol, n)
    model$vtypes      <- "C"
    params            <- list(Presolve = Presolve, OutputFlag = OutputFlag)
    res               <- gurobi::gurobi(model,params)

    weights <- res$x
  }#end if solver == gurobi



  return(list(weights = weights, res = res))

}
