# Function to sample from full conditional of mu_i
complete.conditional.mu <- function(beta,tau,b) rnorm(1,beta/(1+b/tau),1/sqrt(tau+b)) 
# Function to sample from full conditional of tau_i
complete.conditional.tau <- function(beta,mu,c,d) rgamma(1,c+1/2,(beta-mu)^2/2+d)
# Log density of full conditional of beta_i
lcomplete.conditional.beta <- function(betas,tau,mu,i,X,data){
  sum(sapply(1:nrow(data),function(x) 
    X[x,i]*betas[i]*data$All.Star[x]-log(1+exp(X[x,]%*%t(t(betas))))))-tau/2*(betas[i]-mu)^2
}

# b is the prior precision of the mu_i's
# c is the prior shape parameter of the tau_i's
# d is the prior rate parameter of the tau_i's
# tune is the tuning sd of the proposal distribution
# B is the number of iterations
# X is the X matrix
# Data is the full data to be fit
get_draws <- function(b,c,d,tune,B,X,data){
  p <- ncol(X)
  betasInd <- 1:p
  muInd <- (p+1):(2*p)
  tauInd <- (2*p+1):(3*p)
  draws <- matrix(0,nrow=B,ncol=p*3)
  colnames(draws) <- c(paste("b",0:(p-1),sep=""),paste("mu",1:p,sep=""),paste("tau",1:p,sep=""))
  # Proposal distribution
  g <- function(x) rnorm(1,x,tune)
  # Initializing draws for mu
  draws[1,muInd] <- rnorm(p,0,1/sqrt(b)) 
  # Initializing draws from tau
  draws[1,tauInd] <- rgamma(p,c,d)
  # Using those draws to initialize draws from beta
  draws[1,betasInd] <- rnorm(p,draws[1,muInd],1/sqrt(draws[1,tauInd]))
  pb <- txtProgressBar(0,B,style=3)
  for(i in 2:B){
    setTxtProgressBar(pb,i)
    for(j in 1:p){
      # Proposal for the jth beta
      proposal <- g(draws[i-1,betasInd[j]])
      if(j==1){
        new.betas <- c(proposal,draws[i-1,2:p])
        old.betas <- draws[i-1,betasInd]
      } else if(j%in%2:(p-1)){
        new.betas <- c(draws[i,betasInd[1:(j-1)]],proposal,draws[i-1,betasInd[(j+1):p]])
        old.betas <- c(draws[i,betasInd[1:(j-1)]],draws[i-1,betasInd[j:p]])
      } else{
        new.betas <- c(draws[i,betasInd[1:(j-1)]],proposal)
        old.betas <- c(draws[i,betasInd[1:(j-1)]],draws[i-1,betasInd[j:p]])
      }
      # Metropolis ratio on log scale
      ratio <- lcomplete.conditional.beta(new.betas,draws[i-1,tauInd[j]],draws[i-1,muInd[j]],j,X,data)-
        lcomplete.conditional.beta(old.betas,draws[i-1,tauInd[j]],draws[i-1,muInd[j]],j,X,data)
      a <- min(ratio,0)
      if(log(runif(1))<a){
        draws[i,j] <- proposal
      } else{
        draws[i,j] <- draws[i-1,j]
      }
      # Update state of each mu_i
      draws[i,muInd[j]] <- complete.conditional.mu(draws[i,betasInd[j]],draws[i-1,tauInd[j]],b)
      # Update state of each tau_i
      draws[i,tauInd[j]] <- complete.conditional.tau(draws[i,betasInd[j]],draws[i,muInd[j]],c,d)
    }
  }
  return(draws)
}
