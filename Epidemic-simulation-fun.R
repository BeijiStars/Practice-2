epidemic <- function(n=5500000,ne=10,lambda=0.4/n,pei=1/3,pir=1/5,nd=100) {
  ## Epidemic simulation model to compare the daily infection trajectories of the
  ## 10% of the population with the lowest individual transmission probabilities with
  ## those of the whole population.
  ## n = the initial susceptible population; ne = number of initially exposed individuals
  ## lambda = an overall viral infectivity parameter
  ## pei = daily probability that exposed persons enter the infectious state
  ## pir = daily probablity that infectious people leave this state; nd = number of days
  
  
  beta <- rlnorm(n,0,0.5); beta <- beta/mean(beta)## beta[i] is a relative contact rate with others the ith person in the whole population has.
  
  E <- I <- R <- rep(0,nd) ## set up storage for pop in Exposed, Infectious and Recovered state
  E_lb <- E_r <- I_lb <- I_r <- R_lb <- R_r<- rep(0,nd) ## set up storage for pop with specific conditions
  inew <- inew_lb <- inew_r <- rep(0,nd) ## set up storage for infected and recovered pop change
  x <- rep(0,n) ## initialize to susceptible state
  x[sample((1:length(x)), ne)] <- 1 ## start the epidemic by setting ne randomly chosen people to the E stat
  E[1] <- ne
  random_sample <- sample(1:n, size = 0.001*n) ## a random sample of 0.1% of the population to be recorded
  ind_lb <- which(beta < quantile(beta, 0.1))
  for (i in 2:nd) { ## loop over days
    u <- runif(n) ## uniform random deviates
    x[x==2&u<pir] <- 3 ## I -> R
    x[x==0&u<lambda*beta*sum(beta[x==2|x==1])] <- 1 ## S -> E
    x[x==1&u<pei] <- 2 ## E -> I
    E[i] <- sum(x==1) ## exposed people of the whole pop on this day
    I[i] <- sum(x==2) ## infectious people of the whole pop on this day
    R[i] <- sum(x==3) ## recovered people of the whole pop on this day
    E_lb[i] <- sum(x[ind_lb]==1) ## exposed people of the 10% of the pop with lowest beta
    I_lb[i] <- sum(x[ind_lb]==2) ## infectious people of the 10% of the pop with lowest beta
    R_lb[i] <- sum(x[ind_lb]==3) ## recoveries of the 10% of the pop with lowest beta
    E_r[i] <- sum(x[random_sample]==1) ## exposed people of a random sample of 0.1% of the population
    I_r[i] <- sum(x[random_sample]==2) ## infectious people of a random sample of 0.1% of the population
    R_r[i] <- sum(x[random_sample]==3) ## recoveries of a random sample of 0.1% of the population
  }
  inew[-1] <- E[-1] - E[-nd] + I[-1] - I[-nd] + (R[-1]-R[-nd]) ## new infections of each day in the whole pop
  inew_lb[-1] <- E_lb[-1] - E_lb[-nd] + I_lb[-1] - I_lb[-nd] + (R_lb[-1]-R_lb[-nd]) ## new infections of each day in the 10% of the pop with lowest beta
  inew_r[-1] <- E_r[-1] - E_r[-nd] + I_r[-1] -I_r[-nd]+ (R_r[-1]-R_r[-nd]) ## new infections of each day of a random sample of 0.1% of the populatio
  list(inew=inew,inew_lb=inew_lb,inew_r=inew_r)
} ## Epidemic

epi <- epidemic()

## standardise the data
inew_s <- epi$inew/1e4
inew_lb_s <- epi$inew_lb/1e3
inew_r_s <- epi$inew_r/5e1

par(mfcol=c(1,1))
plot(inew_s,ylim=c(0,max(inew_s)),xlab="day",ylab="N", type = 'l', col=1) ## pop daily new infections (black)
lines(inew_lb_s,col=4) ## cautious 10% new daily infections (blue)
lines(inew_r_s,col=2) ## 0.1% random sample new daily infections
legend("topleft", legend = c("whole poplation per 10000","cautious 10% per 1000","0.1% random sample per 50"), 
lwd = 4, col = c(1,4,2), cex = 1.2, title = "New daily infections", bty = "n")
## label the day on which each trajectory peaks.
text(which.max((inew_s)),max(inew_s), labels = paste("peak on day", which.max((inew_s))), col=1, adj=0)
text(which.max((inew_lb_s)),max(inew_lb_s), labels = paste("peak on day", which.max((inew_lb_s))), col=4, adj=0)
text(which.max((inew_r_s)),max(inew_r_s), labels = paste("peak on day", which.max((inew_r_s))), col=2, adj=0)
