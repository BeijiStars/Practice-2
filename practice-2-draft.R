## initialize relative parameters
n <- 5500000 ## the initial susceptible population
lambda <- 0.4/n ## an overall viral infectivity parameter
ne <- 10 ## initially exposed individuals
pei <- 1/3 ## a daily probability that exposed persons enter the infectious state
pir <- 1/5 ## a daily probablity that infectious people leave this state
nd <- 100 ## number of days
S <- E <- I <- R <- I_lb <- I_r <- inew <- inew_lb <- inew_r<- rep(0,nd) ## set up storage for pop in each state
## beta[i] is a relative contact rate with others the ith person in the whole population has.
beta <- rlnorm(n,0,0.5); beta <- beta/mean(beta)
S[1] <- n-ne; E[1] <- ne

## 0 stands for S, 1 for E, 2 for I, 3 for R
x <- rep(0,n) ## initialize to susceptible state
x[sample((1:length(x)), 10)] <- 1 ## start the epidemic by setting ne randomly chosen people to the E stat                                                  
random_sample <- sample(1:n, size = 0.001*n) ## a random sample of 0:1% of the population

system.time(for (i in 2:nd) {
  ##on the second day
  u <- runif(n) ## uniform random deviates
  x[x==2&u<pir] <- 3 ## I -> R
  x[x==1&u<pei] <- 2 ## E -> I
  x[x==0&u<lambda*beta*sum(beta[x==2])] <- 1 ## S -> E
  I[i] <- sum(x==2); ## infections of the whole pop on this day
  I_lb[i] <- sum(x==2&beta<quantile(beta,0.1)) ## infections of the 10% of the pop with lowest beta
  I_r[i] <- sum(x[random_sample]==2) ## infections of a random sample of 0.1% of the population
  inew[i] <- I[i] - I[i-1] ## store the new infections of dayi in the whole pop
  inew_lb[i] <- I_lb[i] - I_lb[i-1] ## store the new infections of dayi in the 10% of the pop with lowest beta
  inew_r[i] <- I_r[i] -I_r[i-1] ## store the new infections of dayi of a random sample of 0.1% of the population
})

## if S->E is ahead of IR and EI, the epidemic develops and stops much more quickly,
## with peaks over 280000 appear between 30th and 40th days. Conversely, the peaks
## appear after the 40th days, and are about only 150000.
