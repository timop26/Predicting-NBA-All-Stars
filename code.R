library(coda)
library(tcltk)
library(ggplot2)
setwd("~/Fall 2019/Stat 651/project")
all.star <- read.csv("all_star.csv")

vars <- c("PTS","TRB","AST","WS.48","prevAS","W.L.","tv")
means <- apply(all.star[,vars],2,mean)
sds <- apply(all.star[,vars],2,sd)
all.star[,c("PTS","TRB","AST","WS.48","W.L.","tv")] <-  scale(all.star[,c("PTS","TRB","AST","WS.48","W.L.","tv")])
#all.star$east <- scale(all.star$east,scale=1)
#all.star$prevAS <- scale(all.star$prevAS,scale=1)
subset <- all.star[all.star$Year!=2019,]
test <- all.star[all.star$Year==2019,]

# Bayesian Analysis
source("mcmc_sampler.R")

#########################################################
########  Getting Posterior Draws for all Data   ########
#########################################################
vars <- c("PTS","TRB","AST","WS.48","prevAS","W.L.","tv")
# X matrix with performance and non-performance variables
X <- cbind(1,as.matrix(all.star[,vars]))
# Getting independent draws
mdl1chain1 <- get_draws(0.1,1,1,0.5,10000,X,all.star)
mdl1chain2 <- get_draws(0.1,1,1,0.5,10000,X,all.star)
mdl1chain3 <- get_draws(0.1,1,1,0.5,10000,X,all.star)
mdl1chain4 <- get_draws(0.1,1,1,0.5,10000,X,all.star)
mdl1chain5 <- get_draws(0.1,1,1,0.5,10000,X,all.star)
# Removing burn-in
mdl1chain1 <- mdl1chain1[1001:10000,]
mdl1chain2 <- mdl1chain2[1001:10000,]
mdl1chain3 <- mdl1chain3[1001:10000,]
mdl1chain4 <- mdl1chain4[1001:10000,]
mdl1chain5 <- mdl1chain5[1001:10000,]
mdl1 <- rbind(mdl1chain1,mdl1chain2,mdl1chain3,mdl1chain4,mdl1chain5)
# Model diagnostics
# Trace plots for betas
plot.ts(mdl1[,1:8])
# Trace plots of mus
plot.ts(mdl1[,9:16])
# Trace plots of taus
plot.ts(mdl1[,17:24])
# Effective sample size
effectiveSize(mdl1)
# Rhat
calc_rhat <- function(chain1,chain2,chain3,chain4,chain5,parameter){
  n <- nrow(chain1)
  par_draws <- cbind(chain1[,parameter],chain2[,parameter],chain3[,parameter],chain4[,parameter],chain5[,parameter])
  # Chain averages
  chain.avgs <- apply(par_draws,2,mean)
  # Overall average
  overall.avg <- mean(chain.avgs)
  # Between variance
  between.var <- n*var(chain.avgs)
  # Within variances
  within.var <- apply(par_draws,2,var)
  W <- mean(within.var)
  varhat <- (n-1)/n*W+between.var/n
  return(sqrt(varhat/W))
}
# Rhat for each parameter
Rhats <- sapply(1:ncol(mdl1),function(x) calc_rhat(mdl1chain1,mdl1chain2,mdl1chain3,mdl1chain4,mdl1chain5,x))
# Number of accepted draws
accept.draws <- apply(sapply(2:nrow(mdl1),function(x) mdl1[x-1,]!=mdl1[x,]),1,sum)/nrow(mdl1)

########## Posterior Inference ###########
betaDraws <- mdl1[,1:8]
# Proportion of draws above 0 for each coefficient
apply(betaDraws,2,function(x) mean(x>0))
# Mean and CI for Bayesian model
means <- t(round(apply(mdl1,2,mean)[1:8],3))
cis <- round(apply(mdl1,2,function(x) quantile(x,c(0.025,0.975)))[,1:8],3)
# Frequentist estimate and CI
fit.all <- glm(All.Star~PTS+TRB+AST+WS.48+prevAS+W.L.+tv,data=all.star,family="binomial")
coef(fit.all)
confint(fit.all)

############ Posterior Predictive Test #############
# Drawn betas multiplied by observed X
XB <- X%*%t(betaDraws)
logitInv <- function(x) 1/(1+exp(-x))
probs <- logitInv(XB)
# Simulating the number of All-Stars
num.all.stars <- apply(probs,2,function(x) rbinom(nrow(probs),1,x))
# Average simulated All-Stars by year
all.stars.by.year <- sapply(2009:2019,function(x) apply(num.all.stars[all.star$Year==x,],2,sum))
apply(all.stars.by.year,2,mean)
# True number of All-Stars by year
true.AS <- sapply(2009:2019,function(x) sum(all.star$All.Star[all.star$Year==x]))
# Posterior predictive p-values
sapply(1:11,function(x) mean(all.stars.by.year[,x]>true.AS[x]))

############ Visualizing effects of non-performance variables #############
# Function to plot posterior predictive probabilities under various situations
density_plots <- function(post.pred,labels="",constant1,size,var,colors=c("#8980F5","#ff0000","#218928")){
  const <- ifelse(constant1==TRUE,"(PTS=20.9, TRB=6.5, AST=4.0, WS/48=0.1)","(PTS=20.9, TRB=8.9, AST=5.8, WS/48=0.2)")
  ggplot(post.pred,aes(x=prob,fill=f))+geom_density(adjust=1.5,alpha=0.5)+xlab("Probability")+
    ylab("Density")+labs(title=paste("Posterior Predictive Probability of Being Selected\n",const))+
    theme_light()+theme(plot.title=element_text(hjust=.5,size=size,face="bold"))+
    scale_fill_manual(values=colors,name=var,labels=labels)
}
# Effect of previous All-Star status for player 2 sds above average for PTS, 1 for AST, REB, and WS/48
post.pred1 <- logitInv(cbind(1,2,1,1,1,0,0,0)%*%t(betaDraws))
post.pred2 <- logitInv(cbind(1,2,1,1,1,1,0,0)%*%t(betaDraws))
post.pred <- data.frame(f=as.factor(rep(c("a","b"),each=nrow(mdl1))),
                        prob=c(post.pred1,post.pred2))
freq.pred <- predict(fit.all,newdata=data.frame(PTS=c(2,2),TRB=c(1,1),AST=c(1,1),WS.48=c(1,1),
                                                prevAS=c(0,1),W.L.=c(0,0),tv=c(0,0)),type="response")
plot11 <- density_plots(post.pred,c("No","Yes"),TRUE,11,"Previously \nSelected")+
  geom_vline(xintercept=freq.pred,color=c(I("blue"),I("red")),lty=c(2,3),lwd=I(1))
# Effect of previous All-Star status for player 2 sds above average for PTS AST, REB, and WS/48
post.pred1 <- logitInv(cbind(1,2,2,2,2,0,0,0)%*%t(betaDraws))
post.pred2 <- logitInv(cbind(1,2,2,2,2,1,0,0)%*%t(betaDraws))
post.pred <- data.frame(f=as.factor(rep(c("a","b"),each=nrow(mdl1))),
                        prob=c(post.pred1,post.pred2))
freq.pred <- predict(fit.all,newdata=data.frame(PTS=c(2,2),TRB=c(2,2),AST=c(2,2),WS.48=c(2,2),
                                                prevAS=c(0,1),W.L.=c(0,0),tv=c(0,0)),type="response")
plot12 <- density_plots(post.pred,c("No","Yes"),FALSE,11,"Previously \nSelected")+
  geom_vline(xintercept=freq.pred,color=c(I("blue"),I("red")),lty=c(2,3),lwd=I(1))
# Effect of win proprotion for player 2 sds above average for PTS, 1 for AST, REB, and WS/48
post.pred1 <- logitInv(cbind(1,2,1,1,1,0,-1,0)%*%t(betaDraws))
post.pred2 <- logitInv(cbind(1,2,1,1,1,0,0,0)%*%t(betaDraws))
post.pred3 <- logitInv(cbind(1,2,1,1,1,0,1,0)%*%t(betaDraws))
post.pred <- data.frame(f=as.factor(rep(c("a","b","c"),each=nrow(mdl1))),
                        prob=c(post.pred1,post.pred2,post.pred3))
freq.pred <- predict(fit.all,newdata=data.frame(PTS=c(2,2,2),TRB=c(1,1,1),AST=c(1,1,1),WS.48=c(1,1,1),
                                                prevAS=c(0,0,0),W.L.=c(-1,0,1),tv=c(0,0,0)),type="response")
plot21 <- density_plots(post.pred,c("35%","50%","65%"),TRUE,11,"Win %")+
  geom_vline(xintercept=freq.pred,color=c(I("blue"),I("red"),I("forestgreen")),lty=c(2,3,4),lwd=I(1))
# Effect of win proprotion for player 2 sds above average for PTS AST, REB, and WS/48
post.pred1 <- logitInv(cbind(1,2,2,2,2,0,-1,0)%*%t(betaDraws))
post.pred2 <- logitInv(cbind(1,2,2,2,2,0,0,0)%*%t(betaDraws))
post.pred3 <- logitInv(cbind(1,2,2,2,2,0,1,0)%*%t(betaDraws))
post.pred <- data.frame(f=as.factor(rep(c("a","b","c"),each=nrow(mdl1))),
                        prob=c(post.pred1,post.pred2,post.pred3))
freq.pred <- predict(fit.all,newdata=data.frame(PTS=c(2,2,2),TRB=c(2,2,2),AST=c(2,2,2),WS.48=c(2,2,2),
                                                prevAS=c(0,0,0),W.L.=c(-1,0,1),tv=c(0,0,0)),type="response")
plot22 <- density_plots(post.pred,c("35%","50%","65%"),FALSE,11,"Win %")+
  geom_vline(xintercept=freq.pred,color=c(I("blue"),I("red"),I("forestgreen")),lty=c(2,3,4),lwd=I(1))
# Effect of market size for player 2 sds above average for PTS, 1 for AST, REB, and WS/48
post.pred1 <- logitInv(cbind(1,2,1,1,1,0,0,-1)%*%t(betaDraws))
post.pred2 <- logitInv(cbind(1,2,1,1,1,0,0,0)%*%t(betaDraws))
post.pred3 <- logitInv(cbind(1,2,1,1,1,0,0,1)%*%t(betaDraws))
post.pred <- data.frame(f=as.factor(rep(c("a","b","c"),each=nrow(mdl1))),
                        prob=c(post.pred1,post.pred2,post.pred3))
freq.pred <- predict(fit.all,newdata=data.frame(PTS=c(2,2,2),TRB=c(1,1,1),AST=c(1,1,1),WS.48=c(1,1,1),
                                                prevAS=c(0,0,0),W.L.=c(0,0,0),tv=c(-1,0,1)),type="response")
plot31 <- density_plots(post.pred,c(0.71,1.62,3.95),TRUE,11,"Market Size\n  (Millions)")+
  geom_vline(xintercept=freq.pred,color=c(I("blue"),I("red"),I("forestgreen")),lty=c(2,3,4),lwd=I(1))
# Effect of market size for player 2 sds above average for PTS, AST, REB, and WS/48
post.pred1 <- logitInv(cbind(1,2,2,2,2,0,0,-1)%*%t(betaDraws))
post.pred2 <- logitInv(cbind(1,2,2,2,2,0,0,0)%*%t(betaDraws))
post.pred3 <- logitInv(cbind(1,2,2,2,2,0,0,1)%*%t(betaDraws))
post.pred <- data.frame(f=as.factor(rep(c("a","b","c"),each=nrow(mdl1))),
                        prob=c(post.pred1,post.pred2,post.pred3))
freq.pred <- predict(fit.all,newdata=data.frame(PTS=c(2,2,2),TRB=c(2,2,2),AST=c(2,2,2),WS.48=c(2,2,2),
                                                prevAS=c(0,0,0),W.L.=c(0,0,0),tv=c(-1,0,1)),type="response")
plot32 <- density_plots(post.pred,c(0.71,1.62,3.95),FALSE,11,"Market Size\n  (Millions)")+
  geom_vline(xintercept=freq.pred,color=c(I("blue"),I("red"),I("forestgreen")),lty=c(2,3,4),lwd=I(1))
# Plotting each together
plot_grid(plot11,plot21,plot31,plot12,plot22,plot32)

#####################################################
######## Entertaining other possible priors #########
#####################################################
# b = 64, c = 40, d = 20
vars <- c("PTS","TRB","AST","WS.48","prevAS","W.L.","tv")
data.list <- lapply(c(vars,"All.Star"),function(x) all.star[,x])
names(data.list) <- c(vars,"y")
data.list <- c(data.list,list(n=nrow(subset)))
betas <- paste0("b",0:7)
mus <- paste0("mu",0:7)
taus <- paste0("tau",0:7)
mdl <- "
model{
for(i in 1:n){
  y[i] ~ dbern(pi[i])
  logit(pi[i]) <- b0+b1*PTS[i]+b2*TRB[i]+b3*AST[i]+b4*WS.48[i]+b5*prevAS[i]+b6*W.L.[i]+b7*tv[i]
}
b0~dnorm(mu0,tau0)
b1~dnorm(mu1,tau1)
b2~dnorm(mu2,tau2)
b3~dnorm(mu3,tau3)
b4~dnorm(mu4,tau4)
b5~dnorm(mu5,tau5)
b6~dnorm(mu6,tau6)
b7~dnorm(mu7,tau7)
mu0~dnorm(0,64)
mu1~dnorm(0,64)
mu2~dnorm(0,64)
mu3~dnorm(0,64)
mu4~dnorm(0,64)
mu5~dnorm(0,64)
mu6~dnorm(0,64)
mu7~dnorm(0,64)
tau0~dgamma(40,20)
tau1~dgamma(40,20)
tau2~dgamma(40,20)
tau3~dgamma(40,20)
tau4~dgamma(40,20)
tau5~dgamma(40,20)
tau6~dgamma(40,20)
tau7~dgamma(40,20)
}
"
parms <- c(paste0("b",0:7),paste0("mu",0:7),paste0("tau",0:7))
niter <- 10000
nburn <- 1000
library(R2jags)
jags.sim <- jags(data=data.list,inits=NULL,parameters.to.save=parms,model.file=textConnection(mdl),
                 n.iter=niter,n.burnin=0,n.chains=5,n.thin=1)
sims <- as.mcmc(jags.sim)
chains <- as.matrix(sims)
chains <- chains[-c(1:nburn,(niter+1):(niter+nburn),(2*niter+1):(2*niter+nburn),(3*niter+1):(3*niter+nburn),(4*niter+1):(4*niter+nburn)),]
# Trace plots for each coefficient
plot.ts(chains[,1:8])
# Effective size
effectiveSize(chains)
# Rhat
gelman.diag(sims)

########### Posterior Inference #############
# Posterior inference on Betas
means <- apply(chains[,betas],2,mean)
cis <- apply(chains[,betas],2,function(x) quantile(x,c(0.025,0.975)))

########### Posterior Predictive Test ############
# Drawn betas multiplied by observed X
XB <- X%*%t(chains[,betas])
probs <- logitInv(XB)
# Simulating the number of All-Stars
num.all.stars <- apply(probs,2,function(x) rbinom(nrow(probs),1,x))
# Average simulated All-Stars by year
all.stars.by.year <- sapply(2009:2019,function(x) apply(num.all.stars[all.star$Year==x,],2,sum))
apply(all.stars.by.year,2,mean)
# True number of All-Stars by year
true.AS <- sapply(2009:2019,function(x) sum(all.star$All.Star[all.star$Year==x]))
# Posterior predictive p-values
sapply(1:length(true.AS),function(x) mean(all.stars.by.year[,x]>true.AS[x]))

#######################################################
####### Predicting 2019 All-Stars holding it out ######
#######################################################
########## Using both types of variables ###########
vars <- c("PTS","TRB","AST","WS.48","prevAS","W.L.","tv")
data.list <- lapply(c(vars,"All.Star"),function(x) subset[,x])
names(data.list) <- c(vars,"y")
data.list <- c(data.list,list(n=nrow(subset)))
betas <- paste0("b",0:length(vars))
mus <- paste0("mu",0:length(vars))
taus <- paste0("tau",0:length(vars))
mdl <- "
model{
for(i in 1:n){
  y[i] ~ dbern(pi[i])
  logit(pi[i]) <- b0+b1*PTS[i]+b2*TRB[i]+b3*AST[i]+b4*WS.48[i]+b5*prevAS[i]+b6*W.L.[i]+b7*tv[i]
}
b0~dnorm(mu0,tau0)
b1~dnorm(mu1,tau1)
b2~dnorm(mu2,tau2)
b3~dnorm(mu3,tau3)
b4~dnorm(mu4,tau4)
b5~dnorm(mu5,tau5)
b6~dnorm(mu6,tau6)
b7~dnorm(mu7,tau7)
mu0~dnorm(0,0.1)
mu1~dnorm(0,0.1)
mu2~dnorm(0,0.1)
mu3~dnorm(0,0.1)
mu4~dnorm(0,0.1)
mu5~dnorm(0,0.1)
mu6~dnorm(0,0.1)
mu7~dnorm(0,0.1)
tau0~dgamma(1,1)
tau1~dgamma(1,1)
tau2~dgamma(1,1)
tau3~dgamma(1,1)
tau4~dgamma(1,1)
tau5~dgamma(1,1)
tau6~dgamma(1,1)
tau7~dgamma(1,1)
}
"
parms <- c(betas,mus,taus)
niter <- 50000
jags.sim <- jags(data=data.list,inits=NULL,parameters.to.save=parms,
                  model.file=textConnection(mdl),n.iter=niter,n.burnin=0,n.chains=10,
                  n.thin=10)
sims <- as.mcmc(jags.sim)
chains <- as.matrix(sims)
niter <- 5000
nburn <- 500
chains <- chains[-c(1:nburn,(niter+1):(niter+nburn),(2*niter+1):(2*niter+nburn),(3*niter+1):(3*niter+nburn),
                    (4*niter+1):(4*niter+nburn),(5*niter+1):(5*niter+nburn),(6*niter+1):(6*niter+nburn),
                    (7*niter+1):(7*niter+nburn),(8*niter+1):(8*niter+nburn),(9*niter+1):(9*niter+nburn)),]
plot.ts(chains[,1:8])
# Effective size
effectiveSize(chains)
# Rhat
gelman.diag(sims)

########### Predicting 2019 All-Stars ##############
XB <- cbind(1,as.matrix(test[,vars]))%*%t(chains[,betas])
probs <- logitInv(XB)
probs.over.5 <- apply(probs,1,function(x) mean(x>0.50))
# Most likely to make the All-Star Game
results1 <- data.frame(1:403,name=test$Player[order(probs.over.5,decreasing=TRUE)],
                        prob=round(sort(probs.over.5,decreasing=TRUE),3),
                        AS=test$All.Star[order(probs.over.5,decreasing=TRUE)])
# Of the players predicted to be All-Stars, the number that were selected
results1[results1$prob>0.3,] # 23/28
# Of the true All-Stars, the number that were predicted as All-Stars
results1[results1$AS==1,] # 23/27

########## Using only performance variables ###########
vars <- c("PTS","TRB","AST","WS.48")
data.list <- lapply(c(vars,"All.Star"),function(x) subset[,x])
names(data.list) <- c(vars,"y")
data.list <- c(data.list,list(n=nrow(subset)))
betas <- paste0("b",0:length(vars))
mus <- paste0("mu",0:length(vars))
taus <- paste0("tau",0:length(vars))
mdl <- "
model{
for(i in 1:n){
  y[i] ~ dbern(pi[i])
  logit(pi[i]) <- b0+b1*PTS[i]+b2*TRB[i]+b3*AST[i]+b4*WS.48[i]
}
b0~dnorm(mu0,tau0)
b1~dnorm(mu1,tau1)
b2~dnorm(mu2,tau2)
b3~dnorm(mu3,tau3)
b4~dnorm(mu4,tau4)

mu0~dnorm(0,0.1)
mu1~dnorm(0,0.1)
mu2~dnorm(0,0.1)
mu3~dnorm(0,0.1)
mu4~dnorm(0,0.1)

tau0~dgamma(1,1)
tau1~dgamma(1,1)
tau2~dgamma(1,1)
tau3~dgamma(1,1)
tau4~dgamma(1,1)
}
"
parms <- c(betas,mus,taus)
niter <- 50000
jags.sim <- jags(data=data.list,inits=NULL,parameters.to.save=parms,
                  model.file=textConnection(mdl),n.iter=niter,n.burnin=0,n.chains=10,
                  n.thin=10)
sims <- as.mcmc(jags.sim)
chains2 <- as.matrix(sims)
niter <- 5000
nburn <- 500
chains2 <- chains2[-c(1:nburn,(niter+1):(niter+nburn),(2*niter+1):(2*niter+nburn),(3*niter+1):(3*niter+nburn),
                    (4*niter+1):(4*niter+nburn),(5*niter+1):(5*niter+nburn),(6*niter+1):(6*niter+nburn),
                    (7*niter+1):(7*niter+nburn),(8*niter+1):(8*niter+nburn),(9*niter+1):(9*niter+nburn)),]
# Trace plots for coefficients
plot.ts(chains2[,1:5])
# Effective size
effectiveSize(chains)
# Rhat

########### Predicting 2019 All-Stars ##############
XB2 <- cbind(1,as.matrix(test[,vars]))%*%t(chains2[,betas])
probs2 <- logitInv(XB2)
probs2.over.5 <- apply(probs2,1,function(x) mean(x>0.50))
# Most likely to make the All-Star Game
results2 <- data.frame(1:403,name=test$Player[order(probs2.over.5,decreasing=TRUE)],
                       prob=round(sort(probs2.over.5,decreasing=TRUE),3),
                       AS=test$All.Star[order(probs2.over.5,decreasing=TRUE)])
# Of the players predicted to be All-Stars, the number that were selected
results2[results2$prob>0.3,] # 20/30
# Of the true All-Stars, the number that were predicted as All-Stars
results2[results2$AS==1,] # 20/27

############# Frequentist predictions ###############
# Using both types of variables
fit.both <- glm(All.Star~PTS+TRB+AST+WS.48+prevAS+W.L.+tv,data=subset,family="binomial")
test$predictions <- predict.glm(fit.both,newdata=test,type="response")
test[test$predictions>0.5,c("Player",vars,"All.Star","predictions")] # 23/28
test[test$All.Star==1,c("Player",vars,"All.Star","predictions")] # 23/27
# Using only performance varaibles
fit.pp <- glm(All.Star~PTS+TRB+AST+WS.48,data=subset,family="binomial")
test$predictions <- predict.glm(fit.pp,newdata=test,type="response")
test[test$predictions>0.5,c("Player",vars,"All.Star","predictions")] # 20/30
test[test$All.Star==1,c("Player",vars,"All.Star","predictions")] # 20/27

