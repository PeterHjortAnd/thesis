library(boot)
library(rsample)
library(fitdistrplus)

##Run after one of the survey codes have been run
bt_samples <- bootstraps(nid,times = 2000)
analysis(bt_samples$splits[[1]]) %>% as_tibble()

#Poisson
##Delta metode
var_pois <- 1/(-exp(lambda_pois)/(exp(lambda_pois)-1)^2+exp(lambda_pois)/(lambda_pois*(exp(lambda_pois)-1)))*1/m  #Inverse fischer information
var_Mpois <- (-m*exp(-lambda_pois)/(1-exp(-lambda_pois))^2)^2*var_pois  ##Delta-metode
z025 <- 1.96
(konf_Mpois <- c(M_pois-z025*sqrt(var_Mpois), M_pois+z025*sqrt(var_Mpois)))

##Bootstrap metode
bt_sampois <- bootstraps(nid,times = 3000)
statistic_pois <- function(splits){
  x <- analysis(splits)
  ndata <- x$n
  xbar <- mean(ndata)
  lambda_pois <- uniroot(f = function(lambda) lambda / (1-exp(-lambda)) - xbar, interval = c(0.01, 100))[["root"]]  ##Løser scoreligningen
  return(m/(1-dpois(0,lambda_pois)))
}
bt_sampois$est <- map_dbl(bt_sampois$splits, statistic_pois)
ggplot(bt_sampois, aes(x = est)) + 
  geom_line(stat = "density", adjust = 1.25) + 
  xlab("Estimat af population hunbjørne (poisson)")
quantile(bt_sampois$est, probs = c(.025, .5, .975))
sd(bt_sampois$est)
c(2*M_pois-quantile(bt_sampois$est, 0.975)[[1]],2*M_pois-quantile(bt_sampois$est, 0.025)[[1]])

#Poisgamma
##Delta metode
(var_gamma <- solve(poisgammafit$hessian))  ##SE for fittet
jacobian_gamma <- matrix(0,nrow=1, ncol = 2)
jacobian_gamma[1,] <- c(-m*r_mle*(1-p_mle)^(r_mle-1)/(1-(1-p_mle)^r_mle)^2,m*(1-p_mle)^r_mle*log(1-p_mle)/(1-(1-p_mle)^r_mle)^2)
var_Mgamma <- jacobian_gamma%*%var_gamma%*%t(jacobian_gamma)
(konf_Mgamma <- c(M_gamma-z025*sqrt(var_Mgamma), M_gamma+z025*sqrt(var_Mgamma)))  ## 95% konfidensinterval for M


#Bootstrap metode
set.seed(3)
bt_samgamma <- bootstraps(nid,times = 3000)
statistic_gamma <- function(splits){
  x <- analysis(splits)
  ndata <- x$n
  fit <- mledist(ndata, "poisgamma",start=list(r=1,p=1/2))
  r <- fit$estimate[[1]]
  p <- fit$estimate[[2]]
  return(m/(1-(1-p)^r))
}
bt_samgamma$est <- map_dbl(bt_samgamma$splits, statistic_gamma)
ggplot(bt_samgamma, aes(x = est)) + 
  geom_line(stat = "density", adjust = 1.25) + 
  xlab("Estimat af population hunbjørne (poisgamma)")
quantile(bt_samgamma$est, probs = c(.025, .5, .975))
c(2*M_gamma-quantile(bt_samgamma$est, 0.975)[[1]],2*M_gamma-quantile(bt_samgamma$est, 0.025)[[1]])


#Poislognormal
set.seed(2)
bt_samlogn <- bootstraps(nid,times = 3000)
statistic_logn <- function(splits){
  x <- analysis(splits)
  ndata <- x$n
  fit <- optim(theta00, logLlogn, hessian = TRUE, data = ndata)
  return(m/(1-dpoilog(0,fit$par[[1]],fit$par[[2]])))
}
statistic_logn <- function(splits){
  x <- analysis(splits)
  ndata <- x$n
  fit <- ztplnMLE(ndata)
  mu_mle <- fit$mu
  sig_mle <- fit$sig
  p0 <- dpoilog(0,mu_mle,sig_mle)
  return(m/(1-p0))
}
bt_samlogn$est <- map_dbl(bt_samlogn$splits, statistic_logn)
ggplot(bt_samlogn, aes(x = est)) + 
  geom_line(stat = "density", adjust = 1.25) + 
  xlab("Estimat af population hunbjørne (poislognormal)")
quantile(bt_samlogn$est, probs = c(.025, .5, .975))
c(2*M_logn-quantile(bt_samlogn$est, 0.975)[[1]],2*M_logn-quantile(bt_samlogn$est, 0.025)[[1]])

#2-point mixture
set.seed(1)
bt_sam2p <- bootstraps(nid,times = 3000)
statistic_2p <- function(splits){
  x <- analysis(splits)
  ndata <- x$n
  pi1<-0.5
  pi2<-0.5
  lambda1<- 0.5
  lambda2<-1
  loglik<- rep(NA, 1000)
  loglik[1]<-0
  loglik[2]<- mysum(pi1*(log(pi1)+log(dzeropois(ndata,lambda1))))+mysum(pi2*(log(pi2)+log(dzeropois(ndata,lambda2))))
  tau1<-0
  tau2<-0
  #k<-1
  k<-2
  # loop
  while(abs(loglik[k]-loglik[k-1]) >= 0.0001) {
    # E step
    tau1<-pi1*dzeropois(ndata,lambda1)/(pi1*dzeropois(ndata,lambda1)+pi2*dzeropois(ndata,lambda2))
    tau2<-pi2*dzeropois(ndata,lambda2)/(pi1*dzeropois(ndata,lambda1)+pi2*dzeropois(ndata,lambda2))
    tau1[is.na(tau1)] <- 0.5
    tau2[is.na(tau2)] <- 0.5
    # M step
    pi1<-mysum(tau1)/length(ndata)
    pi2<-mysum(tau2)/length(ndata)
    lambda1<- Emax(ndata,tau1)
    lambda2<-Emax(ndata,tau2)
    #  loglik[k]<-sum(tau1*(log(pi1)+log(dnorm(x,mu1,sigma1))))+sum(tau2*(log(pi2)+log(dnorm(x,mu2,sigma2))))
    loglik[k+1]<-mysum(tau1*(log(pi1)+log(dzeropois(ndata,lambda1))))+mysum(tau2*(log(pi2)+log(dzeropois(ndata,lambda2))))
    k<-k+1
    if(k > 1000) {
      break
      pi1 <- -0.001
    }
  }
  ddisc2mix <- function(x,lambda,w) w[1]*dztpois(x,lambda[1])+w[2]*dztpois(x,lambda[2])
  lambdaest <- c(lambda1,lambda2)
  piest <- c(pi1,pi2)
  disc2pois_p0 <- pi1*dpois(0,lambda1)+pi2*dpois(0,lambda2)
  if(pi1<0){
    return(2000)
  } else{
    return(m/(1-disc2pois_p0))
  }
}

bt_sam2p$est <- map_dbl(bt_sam2p$splits, statistic_2p)
bt_sam2pv2 <- subset(bt_sam2p, est < 2000)
count(bt_sam2pv2)
ggplot(bt_sam2pv2, aes(x = est)) + 
  geom_line(stat = "density", adjust = 1.25) + 
  xlab("Estimat af population hunbjørne (2-p mixture)")
quantile(bt_sam2p$est, probs = c(.025, .5, .975))
c(2*disc2pois_M-quantile(bt_sam2p$est, 0.975)[[1]],2*disc2pois_M-quantile(bt_sam2p$est, 0.025)[[1]])

#3-point mixture
set.seed(0)
bt_sam3p <- bootstraps(nid,times = 3000)
statistic_3p <- function(splits){
  x <- analysis(splits)
  ndata <- x$n
  pi1<-0.5
  pi2<-0.25
  pi3 <- 0.25
  lambda1<- 0.5
  lambda2<-1
  lambda3 <- 5
  loglik<- rep(NA, 1000)
  loglik[1]<-0
  loglik[2]<- mysum(pi1*(log(pi1)+log(dzeropois(ndata,lambda1))))+mysum(pi2*(log(pi2)+log(dzeropois(ndata,lambda2))))+mysum(pi3*(log(pi3)+log(dzeropois(ndata,lambda3))))
  tau1<-0
  tau2<-0
  tau3<-0
  #k<-1
  k<-2
  
  # loop
  while(abs(loglik[k]-loglik[k-1]) >= 0.0001) {
    # E step
    tau1<-pi1*dzeropois(ndata,lambda1)/(pi1*dzeropois(ndata,lambda1)+pi2*dzeropois(ndata,lambda2)+pi3*dzeropois(ndata,lambda3))
    tau2<-pi2*dzeropois(ndata,lambda2)/(pi1*dzeropois(ndata,lambda1)+pi2*dzeropois(ndata,lambda2)+pi3*dzeropois(ndata,lambda3))
    tau3<- pi3*dzeropois(ndata,lambda3)/(pi1*dzeropois(ndata,lambda1)+pi2*dzeropois(ndata,lambda2)+pi3*dzeropois(ndata,lambda3))
    tau1[is.na(tau1)] <- 0.5
    tau2[is.na(tau2)] <- 0.5
    tau3[is.na(tau3)] <- 0.5
    # M step
    pi1<-mysum(tau1)/length(ndata)
    pi2<-mysum(tau2)/length(ndata)
    pi3 <- mysum(tau3)/length(ndata)
    lambda1<- Emax(ndata,tau1)
    lambda2<-Emax(ndata,tau2)
    lambda3 <- Emax(ndata,tau3)
    #  loglik[k]<-sum(tau1*(log(pi1)+log(dnorm(x,mu1,sigma1))))+sum(tau2*(log(pi2)+log(dnorm(x,mu2,sigma2))))
    loglik[k+1]<-mysum(tau1*(log(pi1)+log(dzeropois(ndata,lambda1))))+mysum(tau2*(log(pi2)+log(dzeropois(ndata,lambda2))))+mysum(tau3*(log(pi3)+log(dzeropois(ndata,lambda3))))
    k<- k+1
    if(k > 1000) {
      break
      pi1 <- -0.001
    }
  }
  ddisc3mix <- function(x,lambda,w) w[1]*dztpois(x,lambda[1])+w[2]*dztpois(x,lambda[2])+w[3]*dztpois(x,lambda[3])
  (lambdaest <- c(lambda1,lambda2,lambda3))  ## Estimater af lambda
  (piest <- c(pi1,pi2,pi3))  ## Estimater af pi
  disc3pois_p0 <- pi1*dpois(0,lambda1)+pi2*dpois(0,lambda2)+pi3*dpois(0,lambda3)
    if(pi1<0){
    return(2000)
    } else{
    return(m/(1-disc3pois_p0)) 
  }

}

bt_sam3p$est <- map_dbl(bt_sam3p$splits, statistic_3p)
bt_sam3pv2 <- subset(bt_sam3p, est < 2000)
count(bt_sam3pv2)
ggplot(bt_sam3pv2, aes(x = est)) + 
  geom_line(stat = "density", adjust = 1.25) + 
  xlab("Estimat af population hunbjørne (3-p mixture)")
quantile(bt_sam3p$est, probs = c(.025, .5, .975))
c(2*disc3pois_M-quantile(bt_sam3p$est, 0.975)[[1]],2*disc3pois_M-quantile(bt_sam3p$est, 0.025)[[1]])

