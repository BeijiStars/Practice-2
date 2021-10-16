epidemic <- function(n=5500000,ne=10,lambda=0.4/n,pei=1/3,pir=1/5,nd=100) {
  ## Epidemic simulation model to compare the daily infection trajectories of the
  ## 10% of the population with the lowest individual transmission probabilities with
  ## those of the whole population.
  ## n = the initial susceptible population; ne = number of initially exposed individuals
  ## lambda = an overall viral infectivity parameter
  ## pei = daily probability that exposed persons enter the infectious state
  ## pir = daily probablity that infectious people leave this state; nd = number of days
  
  
  beta <- rlnorm(n,0,0.5); beta <- beta/mean(beta)## beta[i] is a relative contact rate with others the ith person in the whole population has.
  
  S <- E <- I <- R <- I_lb <- I_r <- R_lb <- R_r<- rep(0,nd) ## set up storage for pop in each state
  inew <- inew_lb <- inew_r <- rep(0,nd) ## set up storage for pop change in I and R
  x <- rep(0,n) ## initialize to susceptible state
  x[sample((1:length(x)), ne)] <- 1 ## start the epidemic by setting ne randomly chosen people to the E stat
  S[1] <- n-ne; E[1] <- ne
  random_sample <- sample(1:n, size = 0.001*n) ## a random sample of 0.1% of the population to be recorded
  
  for (i in 2:nd) { ## loop over days
    u <- runif(n) ## uniform random deviates
    x[x==2&u<pir] <- 3 ## I -> R
    x[x==0&u<lambda*beta*sum(beta[x==2])] <- 1 ## S -> E
    x[x==1&u<pei] <- 2 ## E -> I
    I[i] <- sum(x==2) ## infections of the whole pop on this day
    R[i] <- sum(x==3) ## recoveries of the whole pop on this day
    I_lb[i] <- sum(x==2&beta<quantile(beta,0.1)) ## infections of the 10% of the pop with lowest beta
    R_lb[i] <- sum(x==3&beta<quantile(beta,0.1)) ## recoveries of the 10% of the pop with lowest beta
    I_r[i] <- sum(x[random_sample]==2) ## infections of a random sample of 0.1% of the population
    R_r[i] <- sum(x[random_sample]==3) ## recoveries of a random sample of 0.1% of the population
    inew[i] <- I[i] - I[i-1] + (R[i]-R[i-1]) ## new infections of dayi in the whole pop
    inew_lb[i] <- I_lb[i] - I_lb[i-1] + (R_lb[i]-R_lb[i-1]) ## new infections of dayi in the 10% of the pop with lowest beta
    inew_r[i] <- I_r[i] -I_r[i-1]+ (R_r[i]-R_r[i-1]) ## new infections of dayi of a random sample of 0.1% of the population
  }
  list(inew=inew,inew_lb=inew_lb,inew_r=inew_r)
} ## Epidemic

epi <- epidemic()
epi$inew

