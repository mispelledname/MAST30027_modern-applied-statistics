### Question 1

### Fit a binomial regression model to the orings data 
### Complementary log-log link
### Use R but not glm function

## a) Compute MLEs of the parameters in the model 

# load data
library(faraway)
data(orings)
str(orings)

# -- insert derivation of log-likelihood of binomial with clogloglink function --
# log likelihood function
# logL <- function(beta, orings) {
#   eta <- cbind(1, orings$temp) %*% beta
#   return (sum(orings$damage * log(1 - exp(-exp(eta))) - (6-orings$damage) * exp(eta)))
# }
logL <- function(beta, orings) {
  eta <- cbind(1, orings$temp) %*% beta
  return (sum(orings$damage * log(exp(exp(eta)) - 1) - 6*exp(eta)))
}
# determine estimates 
(betahat <- optim(c(10, -0.1), logL, orings=orings, control=list(fnscale=-1))$par)

## b) Compute 95% CI for the estimates of the parameters 

# derive standard error from first principles 
# -- insert derivation of Fisher's information --
# derive standard errors for betahat 
library(VGAM)
(phat <- clogloglink(betahat[1] + orings$temp * betahat[2], inverse=TRUE))
mult <- (6 * (1-phat) * log(1-phat) / phat) * log(1-phat)
I11 <- sum(mult)
I12 <- sum(orings$temp * mult)
I22 <- sum(orings$temp^2 * mult)
Iinv <- solve(matrix(c(I11, I12, I12, I22), 2, 2))
(se.betahat1 <- sqrt(Iinv[1,1]))
(se.betahat2 <- sqrt(Iinv[2,2]))

# compute CI for betahat
betahat[1] + c(-1,1) * qnorm(0.975) * se.betahat1
betahat[2] + c(-1,1) * qnorm(0.975) * se.betahat2

## c) Likelihood ratio test for the significance of the temperature coefficient

# -- insert set up for likelihood ratio test --
# compute the maximum log-likelihood for the full model
(mll.full <- logL(betahat, orings))

# compute the maximum log-likelihood for the reduced model
y <- orings$damage
n <- rep(6, length(y))
phatN <- sum(y) / sum(n)
(mll.red <- sum(y) * log(phatN) + sum(6-y) * log(1-phatN))

# LR test statistic
(LR <- -2 * (mll.red - mll.full))

# p-value 
pchisq(LR, df=1, lower=FALSE)
# -- insert conclusion -- 

## d) Compute an estimate of the probability of damage when
##    temperature is 31 F
##    (95% CI of estimate as well)

# estimate of probability
etahat <- betahat[1] + betahat[2] * 31
(p.31f <- clogloglink(etahat, inverse=TRUE))

# compute CI for eta
si2 <- matrix(c(1, 31), 1, 2) %*% Iinv %*% matrix(c(1, 31), 2, 1)
eta_l = etahat - qnorm(0.975) * sqrt(si2)
eta_r = etahat + qnorm(0.975) * sqrt(si2)
c(clogloglink(eta_l, inverse=TRUE), clogloglink(eta_r, inverse=TRUE))

## e) Make a plot comparing the fitted complementary model to the fitted logit
##    model. 

# logit model 
logitmodel <- glm(cbind(damage, 6-damage) ~ temp, family=stats::binomial, orings)

# plot fitted logit model
plot(damage/6 ~ temp, orings, xlim=c(25,85), ylim=c(0,1), xlab="Temperature", ylab="Prob of damage")
x <- seq(25,85,1)
ilogit <- function(x) exp(x) / (1+exp(x))
lines(x, ilogit(logitmodel$coefficients[1] + logitmodel$coefficients[2]*x), col="red")

# plot fitted cloglog model
icloglog <- function(x) 1 - exp(-exp(x))
lines(x, icloglog(betahat[1] + betahat[2] * x), col="blue")

### Question 2

# load pima set
library(faraway)
missing <- with(pima, missing <- glucose==0 | diastolic==0 | triceps==0 | bmi==0)
pima_subset = pima[!missing, c(6,9)]
str(pima_subset)

# fit a binomial regression model with logit link
pima_model <- glm(cbind(test, 1-test) ~ bmi, family=stats::binomial, pima)

## a) Amount of increase in the log(odds) when the bmi increases by 5. 


## b) 95% CI for the estimate. 





