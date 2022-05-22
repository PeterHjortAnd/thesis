library(tidyverse)
library(fitdistrplus)
library(ztpln)
library(poilog)
library(boot)
library(rsample)
library(actuar)
library(mixComp)

## 2015 analyse

data <- read.csv("https://raw.githubusercontent.com/mskoldSU/bachelor_bear/main/data/captures.csv")
head(data)
hundata2015 <- subset(data, sex == "Hona" & year == 2015)  ##hunbjørne
tail(hundata2015)

m <- n_distinct(hundata2015$id)   ##Antal obs individer
nid <- hundata2015 %>% count(id)
ndata2 <- nid$n  ##Data vi analyserer

##Simpel poisson fit
xbar <- mean(ndata2)
lambda_pois <- uniroot(f = function(lambda) lambda / (1-exp(-lambda)) - xbar, interval = c(0.01, 100))[["root"]]  ##Løser scoreligningen
(M_pois <- m/(1-dpois(0,lambda_pois)))
g <- nid %>% count(n) %>%
  ggplot(aes(x = n, y = nn)) + geom_col() + theme_bw() +
  labs(x = "Number of captures (poisson, 2015)", y = "Number of individuals") +
  geom_point(data = tibble(x = 1:40, y = dztpois(x,lambda_pois) * m), aes(x = x, y = y), color = "red", size = 2) +
  xlim(c(0, 20))
g + theme_grey(base_size = 18)

var_pois <- 1/(-exp(lambda_pois)/(exp(lambda_pois)-1)^2+exp(lambda_pois)/(lambda_pois*(exp(lambda_pois)-1)))  #Inverse fischer information
var_Mpois <- (-m*exp(-lambda_pois)/(1-exp(-lambda_pois))^2)^2*var_pois  ##Delta-metode
z025 <- 1.96
(konf_Mpois <- c(M_pois-z025*sqrt(var_Mpois), M_pois+z025*sqrt(var_Mpois)))  #95% konfidensinterval for M

## Pois-Gamma distribution
dpoisgamma <- function(x,r,p) gamma(x+r)/(factorial(x)*gamma(r))*((1-p)^r*p^x)/(1-(1-p)^r)
logL_gamma <- function(theta,data){
  -sum(log(gamma(data+theta[1]))-log(factorial(data))-log(gamma(theta[1]))+theta[1]*log(1-theta[2])+data*log(theta[2])-log(1-(1-theta[2])^theta[1]))
}
theta0 <- c(1,1/2)
poisgammafit <- optim(theta0,logL_gamma,hessian = TRUE, data = ndata2) ##Maksimerer loglikelihood
r_mle <- poisgammafit$par[[1]]
p_mle <- poisgammafit$par[[2]]
p0_gamma <- (1-p_mle)^r_mle   ###De ikke-observerede bjørne ved disse estimater
(M_gamma <- m/(1-p0_gamma)) ##Estimat af antal hunbjørne

g <- nid %>% count(n) %>%
  ggplot(aes(x = n, y = nn)) + geom_col() + theme_bw() +
  labs(x = "Number of captures (poisgamma, 2015)", y = "Number of individuals") +
  geom_point(data = tibble(x = 1:40, y = dpoisgamma(x,r_mle,p_mle) * m), aes(x = x, y = y), color = "red", size = 2) +
  xlim(c(0, 20)) ## Plot af fittet
g + theme_grey(base_size = 18)

(var_gamma <- solve(poisgammafit$hessian))  ##SE for fittet
jacobian_gamma <- matrix(0,nrow=1, ncol = 2)
jacobian_gamma[1,] <- c(-r_mle*(1-p_mle)^(r_mle-1),(1-p_mle)^r_mle*log(1-p_mle))
varp0_gamma <- jacobian_gamma%*%var_gamma%*%t(jacobian_gamma)
SE_Mgamma <- sqrt((m/(1-p0_gamma)^2)^2*varp0_gamma)
(konf_Mgamma <- c(M_gamma-z025*SE_Mgamma, M_gamma+z025*SE_Mgamma))  ## 95% konfidensinterval for M


## Pois-Lognormal distribution
logLlogn <- function(theta,data){
  -sum(dztpln(data,theta[1],theta[2],log = TRUE))
}
theta00 <- c(1/2,1)
fitlogn <- optim(theta00, logLlogn, hessian = TRUE, data = ndata2) ## Fit. Maksimerer loglike

(mu_mle <- fitlogn$par[[1]])
(sig_mle <- fitlogn$par[[2]])
p0_logn <- dpoilog(0,mu_mle,sig_mle)   ###De ikke-observerede bjørne ved disse estimater
(M_logn <- m/(1-p0_logn)) ## Antal hunbjørne i 2020

g <- nid %>% count(n) %>%
  ggplot(aes(x = n, y = nn)) + geom_col() + theme_bw() +
  labs(x = "Number of captures (poislognormal, 2015)", y = "Number of individuals") +
  geom_point(data = tibble(x = 1:40, y = dztpln(x,mu_mle,sig_mle) * m), aes(x = x, y = y), color = "red", size = 2) +
  xlim(c(0, 20))  ##Plot
g + theme_grey(base_size = 18)


##2-point discrete pois distribution
# initial values
pi1_2p <-0.5
pi2_2p <-0.5
lambda1_2p <- 0.5
lambda2_2p <-1
loglik_2p <- rep(NA, 1000)

dzeropois <- function(x,lambda) {
  lambda^x/((exp(lambda)-1)*factorial(x))
}
mysum <- function(x) {
  sum(x[is.finite(x)])
}
Emax<- function(dat,tau) {
  uniroot(f = function(lambda) mysum(tau*(dat/lambda -1/(1-exp(-lambda)))), interval = c(0.00001, 1000))[["root"]]
} ##Funktion til at maksimere i E-steppet

loglik_2p[1] <-0
loglik_2p[2] <- mysum(pi1_2p*(log(pi1_2p)+log(dzeropois(ndata2,lambda1_2p))))+mysum(pi2_2p*(log(pi2_2p)+log(dzeropois(ndata2,lambda2_2p))))
tau1<-0
tau2<-0
#k<-1
k<-2

# loop
while(abs(loglik_2p[k]-loglik_2p[k-1]) >= 0.00001) {
  # E step
  tau1<-pi1_2p*dzeropois(ndata2,lambda1_2p)/(pi1_2p*dzeropois(ndata2,lambda1_2p)+pi2_2p*dzeropois(ndata2,lambda2_2p))
  tau2<-pi2_2p*dzeropois(ndata2,lambda2_2p)/(pi1_2p*dzeropois(ndata2,lambda1_2p)+pi2_2p*dzeropois(ndata2,lambda2_2p))
  tau1[is.na(tau1)] <- 0.5
  tau2[is.na(tau2)] <- 0.5
  
  # M step
  pi1_2p <-mysum(tau1)/length(ndata2)
  pi2_2p <-mysum(tau2)/length(ndata2)
  
  lambda1_2p <- Emax(ndata2,tau1)
  lambda2_2p <-Emax(ndata2,tau2)
  
  loglik_2p[k+1]<-mysum(tau1*(log(pi1_2p)+log(dzeropois(ndata2,lambda1_2p))))+mysum(tau2*(log(pi2_2p)+log(dzeropois(ndata2,lambda2_2p))))
  k<-k+1
}
ddisc2mix <- function(x,lambda,w) w[1]*dztpois(x,lambda[1])+w[2]*dztpois(x,lambda[2])
lambdaest_2p <- c(lambda1_2p,lambda2_2p)
piest_2p <- c(pi1_2p,pi2_2p)
disc2pois_p0 <- pi1_2p*dpois(0,lambda1_2p)+pi2_2p*dpois(0,lambda2_2p)  ## Ikke obs bjørne
(disc2pois_M <- m/(1-disc2pois_p0)) ##Estimat af hunbjørne 2020

g <- nid %>% count(n) %>%
  ggplot(aes(x = n, y = nn)) + geom_col() + theme_bw() +
  labs(x = "Number of captures (2p-mixture, 2015)", y = "Number of individuals") +
  geom_point(data = tibble(x = 1:40, y = ddisc2mix(x,lambdaest_2p,piest_2p) * m), aes(x = x, y = y), color = "red", size = 2) +
  xlim(c(0, 20))
g + theme_grey(base_size = 18)


##3-point discrete pois distribution
# initial values
pi1<-0.5
pi2<-0.25
pi3 <- 0.25
lambda1<- 0.5
lambda2<-1
lambda3 <- 5
loglik<- rep(NA, 1000)

loglik[1]<-0
loglik[2]<- mysum(pi1*(log(pi1)+log(dzeropois(ndata2,lambda1))))+mysum(pi2*(log(pi2)+log(dzeropois(ndata2,lambda2))))+mysum(pi3*(log(pi3)+log(dzeropois(ndata2,lambda3))))
tau1<-0
tau2<-0
tau3<-0
#k<-1
k<-2

# loop
while(abs(loglik[k]-loglik[k-1]) >= 0.00001) {
  # E step
  tau1<-pi1*dzeropois(ndata2,lambda1)/(pi1*dzeropois(ndata2,lambda1)+pi2*dzeropois(ndata2,lambda2)+pi3*dzeropois(ndata2,lambda3))
  tau2<-pi2*dzeropois(ndata2,lambda2)/(pi1*dzeropois(ndata2,lambda1)+pi2*dzeropois(ndata2,lambda2)+pi3*dzeropois(ndata2,lambda3))
  tau3<- pi3*dzeropois(ndata2,lambda3)/(pi1*dzeropois(ndata2,lambda1)+pi2*dzeropois(ndata2,lambda2)+pi3*dzeropois(ndata2,lambda3))
  tau1[is.na(tau1)] <- 0.5
  tau2[is.na(tau2)] <- 0.5
  tau3[is.na(tau3)] <- 0.5
  
  # M step
  pi1<-mysum(tau1)/length(ndata2)
  pi2<-mysum(tau2)/length(ndata2)
  pi3 <- mysum(tau3)/length(ndata2)
  
  lambda1<- Emax(ndata2,tau1)
  lambda2<-Emax(ndata2,tau2)
  lambda3 <- Emax(ndata2,tau3)
  
  loglik[k+1]<-mysum(tau1*(log(pi1)+log(dzeropois(ndata2,lambda1))))+mysum(tau2*(log(pi2)+log(dzeropois(ndata2,lambda2))))+mysum(tau3*(log(pi3)+log(dzeropois(ndata2,lambda3))))
  k<-k+1
}

ddisc3mix <- function(x,lambda,w) w[1]*dztpois(x,lambda[1])+w[2]*dztpois(x,lambda[2])+w[3]*dztpois(x,lambda[3])
(lambdaest <- c(lambda1,lambda2,lambda3))  ## Estimater af lambda
(piest <- c(pi1,pi2,pi3))  ## Estimater af pi
disc3pois_p0 <- pi1*dpois(0,lambda1)+pi2*dpois(0,lambda2)+pi3*dpois(0,lambda3)  ## Ikke obs bjørne
(disc3pois_M <- m/(1-disc3pois_p0)) ##Estimat af hunbjørne 2020

g <- nid %>% count(n) %>%
  ggplot(aes(x = n, y = nn)) + geom_col() + theme_bw() +
  labs(x = "Number of captures (3p-mixture, 2015)", y = "Number of individuals") +
  geom_point(data = tibble(x = 1:40, y = ddisc3mix(x,lambdaest,piest) * m), aes(x = x, y = y), color = "red", size = 2) +
  xlim(c(0, 20))
g + theme_grey(base_size = 18)


## AIC
(value_pois <- sum(dztpois(ndata2,lambda_pois,log = TRUE)))
(value_gamma <- -poisgammafit$value)
(value_logn <- -fitlogn$value)
(value_2p<- sum(log(ddisc2mix(ndata2,lambdaest_2p,piest_2p))))
(value_3p <- sum(log(ddisc3mix(ndata2,lambdaest,piest))))


(AICpoi <- 2*1-2*value_pois)
(AICgamma <- 2*2-2*value_gamma)
(AIClogn <- 2*2-2*value_logn)
(AIC2p <- 2*4-2*value_2p)
(AIC3p <- 2*6-2*value_3p)
minAIC <- min(AICpoi,AICgamma,AIClogn,AIC2p,AIC3p)
dAICpoi <- AICpoi-minAIC
dAICgamma <- AICgamma-minAIC
dAIClogn <- AIClogn-minAIC
dAIC2p <- AIC2p-minAIC
dAIC3p <- AIC3p-minAIC


## Goodness of fit og p-værdier
hist <- hist(ndata2,breaks=1:100, right=FALSE, plot = FALSE) #Histogram

ford_pois <- c(dztpois(1:99,lambda_pois))
(gof_pois <- chisq.test(hist$counts, p=ford_pois, rescale.p=TRUE, simulate.p.value=TRUE))
ford_gamma <- c(dpoisgamma(1:99,r_mle,p_mle))
(gof_gamma <- chisq.test(hist$counts, p=ford_gamma, rescale.p=TRUE, simulate.p.value=TRUE))
ford_logn <- c(dztpln(1:99,mu_mle,sig_mle))
(gof_logn <- chisq.test(hist$counts, p=ford_logn, rescale.p=TRUE, simulate.p.value=TRUE))
ford_2p <- c(ddisc2mix(1:99,lambdaest_2p,piest_2p))
(gof_2p <- chisq.test(hist$counts, p=ford_2p, rescale.p=TRUE, simulate.p.value=TRUE))
ford_3p <- c(ddisc3mix(1:99,lambdaest,piest))
(gof_3p <- chisq.test(hist$counts, p=ford_3p, rescale.p=TRUE, simulate.p.value=TRUE))
gof_gamma$statistic[[1]]
gof_gamma$p.value

tab <- matrix(c(M_pois, M_gamma, M_logn, disc2pois_M, disc3pois_M, dAICpoi, dAICgamma, dAIClogn, dAIC2p, dAIC3p, gof_pois$statistic[[1]], gof_gamma$statistic[[1]],gof_logn$statistic[[1]],gof_2p$statistic[[1]],gof_3p$statistic[[1]], gof_pois$p.value, gof_gamma$p.value,gof_logn$p.value,gof_2p$p.value,gof_3p$p.value), ncol=4, byrow=FALSE)
colnames(tab) <- c('N','dAIC','X-squared','p-value')
rownames(tab) <- c('Pois', 'Poisgamma','Poislognormal','2-point mixture','3-point mixture')
tab <- as.table(tab)
tab

estimater <- matrix(rep(0,24),nrow = 6, ncol = 5)
estimater[,1] <- c(1,2,3,4,5,6)
estimater[,2] <- c(684-m, M_pois-m,M_gamma-m,M_logn-m,disc2pois_M-m,disc3pois_M-m)
estimater[,3] <- c(645-m,589-m,721-m[], 681-m[], 623-m[], 585-m[])
estimater[,4] <- c(723-m, 598-m[],955-m[],754-m[],656-m[],669-m[])
estimater[,5] <- c(rep(2015,6))
est15 <- estimater ##Estimater af ikke-observerede bjørne
esttyve <- as.data.frame(estimater)
colnames(esttyve) <- c(c('model', 'est','lower','upper'))
Models <- c('off.','poisson','poisgamma','poislogn','2-point','3-point')

g <- ggplot(esttyve, aes(x = model, y = est, color = Models)) +        # ggplot2 plot with confidence intervals
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  labs(x = "Models (2015)", y = "Estimates of non-observed bears (2015)") +
  scale_x_continuous(breaks = 1:6)
gplot2015 <- g + theme_bw(base_size = 18)
gplot2015