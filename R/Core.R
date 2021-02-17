#' @name CoxICPen
#' @title Variable Selection for Cox's Model with Interval-Censored Data
#' @description Perform variable selection for Cox regression model with interval-censored data by using the methods proposed in Zhao et al. (2020) and Wu et al. (2020). Can deal with both low-dimensional and high-dimensional data.
#' @usage CoxICPen(LR = LR, x = x, BernD = 3, beta.initial = rep(0,ncol(x)),
#'          lamb = log(nrow(x))/2-2, pen = "BAR")
#' @param LR An n by 2 matrix that contains interval-censored failure times (L, R]. Please set time point R to "Inf" if a subject is right-censored.
#' @param x An n by p covariate matrix.
#' @param BernD The degree of Bernstein polynomials. Default is 3.
#' @param beta.initial The initial values for the regression coefficients in the Cox's model. Default value is 0.
#' @param lamb The value of the tuning parameter of the penalty term. Default is log(n)/2-2.
#' @param pen The penalty function. Choices include "RIDGE", "BAR", "LASSO", "SCAD", "MCP", "SICA", "SELO". Default is "BAR".
#' @return beta: Penalized estimates of the regression coefficients in the Cox's model.
#' @return phi: Estimates of the coefficients in Bernstein Polynomials.
#' @return logL: Log likelihood function based on current estimates.
#' @importFrom stats optim
#' @importFrom foreach foreach
#' @importFrom foreach %do%
#' @export
#' @examples
#'
#' # Generate an example data
#'
#' require(foreach)
#'
#' n <- 300  # Sample size
#' p <- 50   # Number of covariates
#'
#' bet0 <- c(1, -1, 1, -1, rep(0,p-4))  # True values of regression coefficients
#'
#' set.seed(1)
#' x.example <- matrix(rnorm(n*p,0,1),n,p)  # Generate covariates matrix
#'
#' T.example <- c()
#' for (i in 1:n){
#'   T.example[i] <- rexp(1,exp(x.example%*%bet0)[i])  # Generate true failure times
#' }
#'
#' timep <- seq(0,3,,10)
#' LR.example <- c()
#' for (i in 1:n){
#'   obsT <- timep*rbinom(10,1,0.5)
#'   if (max(obsT) < T.example[i]) {LR.example <- rbind(LR.example,c(max(obsT), Inf))} else {
#'     LR.example <- rbind(LR.example,c(max(obsT[obsT<T.example[i]]), min(obsT[obsT>=T.example[i]])))
#'   }
#' }  # Generate interval-censored failure times
#'
#'
#' # Fit Cox's model with penalized estimation
#'
#' model1 <- CoxICPen(LR = LR.example, x = x.example, lamb = 100, pen = "RIDGE")
#' beta.initial <- model1$beta
#'
#' model2 <- CoxICPen(LR = LR.example, x = x.example, beta.initial = beta.initial, pen = "BAR")
#' model2
#'
#' @references Zhao, H., Wu, Q., Li, G., Sun, J. (2020) Simultaneous Estimation and Variable Selection for Interval-Censored Data with Broken Adaptive Ridge Regression. Journal of the American Statistical Association. 115(529):204-216.
#' @references Wu, Q., Zhao, H., Zhu, L., Sun, J. (2020) Variable Selection for High-dimensional Partly Linear Additive Cox Model with Application to Alzheimer's disease. Statistics in Medicines.39(23):3120-3134.
#' @references Zhao, H., Wu, Q., Gilbert, P. B., Chen, Y. Q., Sun, J. (2020) A Regularized Estimation Approach for Case-cohort Periodic Follow-up Studies with An Application to HIV Vaccine Trials. Biometrical Journal. 62(5):1176-1191.

require(foreach)


# Core function

CoxICPen <- function(LR = LR, x = x, BernD = 3, beta.initial = rep(0,ncol(x)),
                     lamb = log(nrow(x))/2-2, pen = "BAR") {

  LR <- LR
  x <- x
  m <- BernD
  p <- ncol(x)
  n <- nrow(x)
  tau <- max(LR[LR[,2]!=Inf,])

  if (length(beta.initial)!=ncol(x)) {
    bet.all <- bet <- rep(0,ncol(x))
  } else {bet.all <- bet <- beta.initial}

  phi.all <- phi <- rep(0,m+1)

  ln_phi <- function(phi) {
    LN <- function(t) {
      u <- ifelse(max(LR[,2])==Inf,tau,max(LR[,2]))
      c <- min(LR[,1])
      tt <- ifelse(t>tau,tau,t)
      phi_s <- cumsum(exp(phi))
      B <- c()
      for (i in 1:(m+1)) B <- cbind(B, matrix(choose(m,(i-1))*((tt-c)/(u-c))^(i-1)*(1-(tt-c)/(u-c))^(m-i+1),ncol=1))
      return (B%*%matrix(phi_s,ncol=1))
    }
    expL <- ifelse(LR[,2]==0 , 1, exp(-LN(LR[,1])*exp(x%*%bet)))
    expR <- ifelse(LR[,2]==Inf,0, exp(-LN(LR[,2])*exp(x%*%bet)))
    return( sum( log(ifelse(expL-expR==0,1e-250,expL-expR)) ) )
  }

  LambNull <- function(t) {
    u <- ifelse(max(LR[,2])==Inf,tau,max(LR[,2]))
    c <- min(LR[,1])
    tt <- ifelse(t>tau,tau,t)
    phi_s <- cumsum(exp(phi))
    B <- c()
    for (i in 1:(m+1)) B <- cbind(B, matrix(choose(m,(i-1))*((tt-c)/(u-c))^(i-1)*(1-(tt-c)/(u-c))^(m-i+1),ncol=1))
    return (B%*%matrix(phi_s,ncol=1))
  }

  dln_phi <- function(phi) {
    LambL <- LambNull(LR[,1])*exp(x%*%bet)
    LambR <- LambNull(LR[,2])*exp(x%*%bet)
    SL <- exp(-LambL)
    SR <- ifelse(LR[,2]==Inf,0,exp(-LambR))
    tt <- ifelse(t>tau,tau,t)
    u <- ifelse(max(LR[,2])==Inf,tau,max(LR[,2]))
    c <- min(LR[,1])
    Bern <- function(t) {
      B <- c()
      for (i in 1:(m+1)) B <- cbind(B, matrix(choose(m,(i-1))*((tt-c)/(u-c))^(i-1)*(1-(tt-c)/(u-c))^(m-i+1),ncol=1))
      return(B)
    }
    dln <- numeric(m+1)
    for (i in 1:(m+1)) {
      dln[i] <- sum( exp(phi[i]+x%*%bet)*(SR*apply(as.matrix(Bern(LR[,2])[,i:(m+1)]),1,sum)-SL*apply(as.matrix(Bern(LR[,1])[,i:(m+1)]),1,sum))/
                       ifelse(abs(SL-SR) < 1e-150,1e-149,SL-SR) )
    }
    return(dln)
  }

  dg_beta <- function(bet, x, LR, dd) {
    LambL <- LambNull(LR[,1])*exp(x%*%bet)
    LambR <- LambNull(LR[,2])*exp(x%*%bet)
    SL <- ifelse(LR[,1]==0 ,1,exp(-LambL))
    SR <- ifelse(LR[,2]==Inf,0,exp(-LambR))

    result <- sum(x[,dd]*(SR*LambR-SL*LambL)/ifelse(abs(SL-SR) < 1e-150,1e-149,SL-SR))
    return(result)
  }

  ddg_beta <- function(bet, x, LR, dd) {
    LambL <- LambNull(LR[,1])*exp(x%*%bet)
    LambR <- LambNull(LR[,2])*exp(x%*%bet)
    SL <- ifelse(LR[,1]==0 ,1,exp(-LambL))
    SR <- ifelse(LR[,2]==Inf,0,exp(-LambR))

    result <- sum( x[,dd]^2*( (SR*LambR-SL*LambL)/ifelse(abs(SL-SR) < 1e-150,1e-149,SL-SR) - SR*SL*(LambL-LambR)^2/(ifelse(abs(SL-SR) < 1e-150,1e-149,SL-SR))^2 ) )
    return(result)
  }

  logL <- function(bet, x, LR) {
    expL <- ifelse(LR[,2]==0 ,1,exp(-LambNull(LR[,1])*exp(x%*%bet)))
    expR <- ifelse(LR[,2]==Inf,0,exp(-LambNull(LR[,2])*exp(x%*%bet)))
    return( sum( log(ifelse(expL-expR==0,1e-250,expL-expR)) ) )
  }

  dLASSO <- function(z) lamb
  dSCAD <- function(z,alph) (z<=lamb)+ifelse(alph*lamb-z>0,alph*lamb-z,0)/(lamb*(alph-1))*(z>lamb)
  dMCP <- function(z,tau2) (z<=tau2*lamb)*(lamb-z/tau2)
  dSICA <- function(z,tau3) lamb*tau3*(tau3+1)/(z+tau3)^2
  dSELO <- function(z,gam2) lamb*gam2/(log(2)*(2*z+gam2)*(z+gam2))
  eps <- 1e-200
  eps0 <- 1e-6

  act.set <- c(1:p)
  for (NN in 1:50) {
    phi <- optim(phi,ln_phi,dln_phi,control=c(fnscale=-1, maxit=20))$par
    phi.all <- rbind(phi.all,c(phi))

    if (pen=="SCAD") {
      foreach::foreach(dd=act.set, .combine='c') %do% {
        bet[dd] <- (bet[dd]*ddg_beta(bet, x, LR, dd)-dg_beta(bet, x, LR, dd))/(ddg_beta(bet, x, LR, dd)-n*dSCAD(abs(bet[dd]),3.7)/(abs(bet[dd])+eps))
      }
    }
    else if (pen=="LASSO") {
      foreach::foreach(dd=act.set, .combine='c') %do% {
        bet[dd] <- (bet[dd]*ddg_beta(bet, x, LR, dd)-dg_beta(bet, x, LR, dd))/(ddg_beta(bet, x, LR, dd)-n*dLASSO(abs(bet[dd]))/(abs(bet[dd])+eps))
      }
    }
    else if (pen=="RIDGE") {
      foreach::foreach(dd=act.set, .combine='c') %do% {
        bet[dd] <- bet[dd] - (dg_beta(bet, x, LR, dd)-2*lamb*bet[dd])/(ddg_beta(bet, x, LR, dd)-2*lamb)
      }
    }
    else if (pen=="BAR") {
      foreach::foreach(dd=act.set, .combine='c') %do% {
        bet[dd] <- bet[dd] - (dg_beta(bet, x, LR, dd)-2*lamb*bet[dd]/(bet[dd]^2+eps))/(ddg_beta(bet, x, LR, dd)-2*lamb/(bet[dd]^2+eps))
      }
    }
    else if (pen=="MCP") {
      foreach::foreach(dd=act.set, .combine='c') %do% {
        bet[dd] <- (bet[dd]*ddg_beta(bet, x, LR, dd)-dg_beta(bet, x, LR, dd))/(ddg_beta(bet, x, LR, dd)-n*dMCP(abs(bet[dd]),1.1)/(abs(bet[dd])+eps))
      }
    }
    else if (pen=="SICA") {
      foreach::foreach(dd=act.set, .combine='c') %do% {
        bet[dd] <- (bet[dd]*ddg_beta(bet, x, LR, dd)-dg_beta(bet, x, LR, dd))/(ddg_beta(bet, x, LR, dd)-n*dSICA(abs(bet[dd]),0.01)/(abs(bet[dd])+eps))
      }
    }
    else if (pen=="SELO") {
      foreach::foreach(dd=act.set, .combine='c') %do% {
        bet[dd] <- (bet[dd]*ddg_beta(bet, x, LR, dd)-dg_beta(bet, x, LR, dd))/(ddg_beta(bet, x, LR, dd)-n*dSELO(abs(bet[dd]),0.01)/(abs(bet[dd])+eps))
      }
    }

    act.set <- which(abs(bet) > eps0)
    bet[abs(bet) < eps0] <- 0
    bet.all <- rbind(bet.all,c(bet))
    if (mean(abs(bet.all[NN,]-bet)) < 1e-4 &&
        mean(abs(phi.all[NN,]-phi)) < 1e-2) break
  }

  time.points <- unique(sort(LR[LR[,2]!=Inf,]))
  Lamb0 <- LambNull(time.points)

  out <- NULL
  out$beta <- bet
  out$phi <- phi
  out$logL <- logL(bet, x, LR)
  out$Lamb0 <- data.frame(time.points, Lamb0)
  return(out)
}

globalVariables("dd", package="CoxICPen")
