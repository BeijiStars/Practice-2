## initialize relative parameters
n <- 5500000 ## the initial susceptible population
lambda <- 0.4/n ## an overall viral infectivity parameter
ne <- 10 ## initially exposed individuals
pei <- 1/3 ## a daily probability that exposed persons enter the infectious state
pir <- 1/5 ## a daily probablity that infectious people leave this state
nd <- 150 ## number of days
S <- E <- I <- R <- inew <- rep(0,nd) ## set up storage for pop in each state
inew[1] <- NA ## we cannot caculate the new infections in day1
## beta[i] is a relative contact rate with others the ith person in the whole population has.
beta <- rlnorm(n,0,0.5); beta <- beta/mean(beta)
S[1] <- n-ne; E[1] <- ne

## 0 stands for S, 1 for E, 2 for I, 3 for R
x <- rep(0,n) ## initialize to susceptible state
x[sample((1:length(x)), 10)] <- 1 ## start the epidemic by setting ne randomly chosen people to the E stat                                                  

system.time(for (i in 2:nd) {
  ##on the second day
  u <- runif(n) ## uniform random deviates
  x[x==2&u<pir] <- 3 ## I -> R
  x[x==0&u<lambda*beta*sum(beta[x==2])] <- 1 ## S -> E, now we allow a person to enter I on the day entering E
  x[x==1&u<pei] <- 2 ## E -> I
  S[i] <- sum(x==0); E[i] <- sum(x==1)
  I[i] <- sum(x==2); R[i] <- sum(x==3)
  inew[i] <- I[i] - I[i-1] ## store the new infections of dayi in the whole pop
})




uu <- runif(length(xx))
xx <- c(0,2,2,0,0,2)
bb <- rlnorm(length(xx),0,0.5); bb <- bb/mean(bb)
xx[xx==0&uu<lambda*bb*sum(bb[x==2])] <- 1 ## S -> E
## xx = 0 is true in 1,4,5; bb[xx=2] is 
sum(bb[xx==2])
lambda*bb[1]*sum(bb[xx==2])