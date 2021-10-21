##By Group 02： Huanyu Su(s2171849), Mia Wu(s2025818), Peiran Zhang(s2184444)
##Github Repo URL: https://github.com/BeijiStars/Practice-2

## Overview
## ********
## The aim is to simulate an epidemic in a large population for 100 days, and try 
## to answer the following questions: what do the daily infection trajectories look like
## for the 10% of the population with the lowest individual transmission possibilities,
## and how do they compare to the infection trajectories for the whole population? To 
## address these questions, the number of new infections each day among the whole population,
## the 10% population with the lowest individual transmission possibilities and a radom
## sample of 0.1% of the population are recorded respectively. After that, a plot showing
## how the daily infections trajectories compare between the whole population and the two
## samples are produced. Finally, by running 10 replicate simulations, the variability in
## the results are visualized.
##
## Algorithm details
## *****************
## 1. The initial susceptible population of size 5.5 million is set, and 10 randomly chosen people 
##    are set to the exposed (E) state.
## 2. Once exposed, people have a daily daily probability 1/3 entering the infectious (I) state.
## 3. Each day an uninfected person, j, has a probability lambda*beta[j]*the sum of betas of all
##    currently of infectious people of being infected,and entering the E state, where lambda is
##    an overall viral infectivity parameter, vector beta contains relative contact rates for all
##    people, and the average value of betap[i] over the whole population is 1.
## 4. The infectious have a daily probability of 1/5 of leaving the infectious state (recovery or 
##    transition to serious disease). Over the course of simulation re-infection is not considered.
## 5. 9 vectors are created to record the number of people in the exposed (E), infectious (I) and 
##    recovered (R) states for the whole population, the 'cautious 10%' and the 0.1% of random sample
##    each day.
## 6. The new infections on day i (i>1) for three populations are then caculated by summing the differences
##    between E state and I state, and then plus the increase of number of individuals in R state.
## 7.


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
  ind_lb <- which(beta < quantile(beta, 0.1)) ## the 10% of the population with lowest betas
  for (i in 2:nd) { ## loop over days
    u <- runif(n) ## uniform random deviates
    x[x==2&u<pir] <- 3 ## I -> R
    x[x==0&u<lambda*beta*sum(beta[x==2])] <- 1 ## S -> E
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


## standardise the data
population <- 5.5e6
inew_s <- epi$inew/population
inew_lb_s <- epi$inew_lb/(population * 0.1)
inew_r_s <- epi$inew_r/(population * 0.001)

par(mfcol=c(1,1))
plot(inew_s,ylim=c(0,max(inew_s,inew_lb_s,inew_r_s)),xlab="day",ylab="N", type = 'l', col=1) ## pop daily new infections (black)
lines(inew_lb_s,col=4) ## cautious 10% new daily infections (blue)
lines(inew_r_s,col='brown') ## 0.1% random sample new daily infections
legend("topleft", legend = c("whole poplation","cautious 10%","0.1% random sample"), 
lwd = 4, col = c(1,4,'brown'), cex = 0.8, title = "New daily infections", bty = "n")
## label the day on which each trajectory peaks.
text(which.max((inew_s)),max(inew_s), labels = paste("Peak on day", which.max((inew_s))), col=1, adj=-0.25, cex = 0.8)
text(which.max((inew_lb_s)),max(inew_lb_s), labels = paste("Peak on day", which.max((inew_lb_s))), col=4, adj=-0.5, cex = 0.8)
text(which.max((inew_r_s)),max(inew_r_s), labels = paste("Peak on day", which.max((inew_r_s))), col='brown', adj=-0.05, cex = 0.8)


#10 replicate simulations and plot 
set.seed(10)
simulate_10 <- lapply(rep(5.5e6,10), epidemic)
par(mfcol=c(3,1))
plot(simulate_10[1][[1]]$inew,xlab="day",ylab="N", type = 'l', col=1)
for (i in 2:10){
  lines(simulate_10[i][[1]]$inew)
}
plot(simulate_10[1][[1]]$inew_lb,xlab="day",ylab="N", type = 'l', col=4)
for (i in 2:10){
  lines(simulate_10[i][[1]]$inew_lb, col=4)
}
plot(simulate_10[1][[1]]$inew_r,xlab="day",ylab="N", type = 'l', col='brown')
for (i in 2:10){
  lines(simulate_10[i][[1]]$inew_r, col='brown')
}