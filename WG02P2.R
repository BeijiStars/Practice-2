##By Group 02： Huanyu Su(s2171849), Mia Wu(s2025818), Peiran Zhang(s2184444)
##Github Repo URL: https://github.com/BeijiStars/Practice-2

## Overview
## ********
## To come up with an epidemic simulation（function, Plots and analysis) in large population in 100 days
## and try to answer the following questions: 
## 1  what do the daily infection trajectories look like for the 10% of the population 
##    with the lowest individual transmission possibilities;
## 2  and how do they compare to the infection trajectories for the whole population? 

## Steps
## ********
## 1  Write one function to simulate the number of new infections each day among the whole population;
##    and simulate the 10% population with the lowest individual transmission possibilities;
##    and simulate a radom sample of 0.1% of the population. 
## 2  Standardise and plot the daily infections trajectories and trajectory peaks compare 
##    among 'whole population', '10% cautious' and '0.1% random sample'.
## 3  Running 10 replicate simulations and plot results in three parts.
## 4  Analyze the implications of these results to the ZOE app data.

## Algorithm details
## *****************
## 1  The initial susceptible population of size 5.5 million is set, and 10 randomly chosen people 
##    are set to the exposed (E) state.
## 2  Once exposed, people have a daily probability 1/3 entering the infectious (I) state.
## 3  Each day an uninfected person, j, has a probability lambda*beta[j]*the sum of betas of all
##    currently of infectious people of being infected,and entering the E state, where lambda is
##    an overall viral infectivity parameter, vector beta contains relative contact rates for all
##    people, and the average value of betap[i] over the whole population is 1.
## 4  The infectious have a daily probability of 1/5 of leaving the infectious state (recovery or 
##    transition to serious disease). Over the course of simulation re-infection is not considered.
## 5  9 vectors are created to record the number of people in the exposed (E), infectious (I) and 
##    recovered (R) states for the whole population, the 'cautious 10%' and the 0.1% of random sample
##    each day.
## 6  The new infections on day i (i>1) for three populations are then caculated by summing the differences
##    between E state and I state, and then plus the increase of number of individuals in R state.
## 7  In order for the 3 curves to be displayed in one graph, they need to be divided 
##    by their own sample size in order to be normalised.
## 8  Simulation 10 was repeated, placing the 3 curves in 3 separate graphs to compare their fluctuations 
##    and thus explain the impact of using the ZOE app on the infection trajectory.

##  Step1 
##  write the function epidemic

epidemic <- function(n=5500000,ne=10,lambda=0.4/n,pei=1/3,pir=1/5,nd=100) {
  ## whole population (n = 5,500,000),
  ## randomly chosen people in E state (ne = 10),
  ## constant number (lambda = 0.4 / n),
  ## E <- I daily probability (pei = 1/3),
  ## I <- R or serious disease daily probabiliyu (pir = 1/5)
  ## simulate 100 days (nd = 100)
  beta <- rlnorm(n,0,0.5); beta <- beta/mean(beta)
  ## beta[i] is a relative contact rate with others the ith person in the whole population has
  E <- I <- R <- rep(0,nd) ## set up storage for pop in Exposed, Infectious and Recovered state
  E_lb <- E_r <- I_lb <- I_r <- R_lb <- R_r<- rep(0,nd) ## set up storage for pop with specific conditions
  inew <- inew_lb <- inew_r <- rep(0,nd) ## set up storage for infected and recovered pop change
  x <- rep(0,n) ## initialize to susceptible state
  
  ## start the epidemic by setting ne randomly chosen people to the E state
  x[sample((1:length(x)), ne)] <- 1 
  E[1] <- ne
  random_sample <- sample(1:n, size = 0.001*n) ## a random sample of 0.1% of the population to be recorded
  ind_lb <- which(beta < quantile(beta, 0.1)) ## the 10% of the population with lowest betas
  
  ## loop over days
  for (i in 2:nd) { 
    u <- runif(n) ## uniform random deviates
    x[x==2&u<pir] <- 3 ## I -> R
    x[x==0&u<lambda*beta*sum(beta[x==2])] <- 1 ## S -> E
    x[x==1&u<pei] <- 2 ## E -> I
    E[i] <- sum(x==1) ## exposed people of the whole pop on thid day
    I[i] <- sum(x==2) ## infectious people of the whole pop on this day
    R[i] <- sum(x==3) ## recoveries of the whole pop on this day
    E_lb[i] <- sum(x[ind_lb]==1) ## exposed people of the 10% of the pop with lowest beta
    I_lb[i] <- sum(x[ind_lb]==2) ## infectious people of the 10% of the pop with lowest beta
    R_lb[i] <- sum(x[ind_lb]==3) ## recoveries of the 10% of the pop with lowest beta
    E_r[i] <- sum(x[random_sample]==1) ## exposed people of a random sample of 0.1% of the population
    I_r[i] <- sum(x[random_sample]==2) ## infectious people of a random sample of 0.1% of the population
    R_r[i] <- sum(x[random_sample]==3) ## recoveries of a random sample of 0.1% of the population
  }
  
  inew[-1] <- E[-1] + I[-1] - (E[-nd] + I[-nd]) + (R[-1]-R[-nd]) 
  ## new infections of each day in the whole pop
  inew_lb[-1] <- E_lb[-1] + I_lb[-1] - (E_lb[-nd] + I_lb[-nd]) + (R_lb[-1]-R_lb[-nd])
  ## new infections of each day in the 10% of the pop with lowest beta
  inew_r[-1] <- E_r[-1] + I_r[-1] - (E_r[-nd] + I_r[-nd]) + (R_r[-1]-R_r[-nd]) 
  ## new infections of each day of a random sample of 0.1% of the populatio
  
  list(inew=inew,inew_lb=inew_lb,inew_r=inew_r)
} ## Epidemic

epi <- epidemic() ## simulate epidemic


##  Step2
##  Standardise and plot the daily infections trajectories

##  Standarlize new infections in three samples（whole pop ,10% cautious and 0.1% random sample)
## It means the number of infections per 100 people.
population <- 5.5e6
inew_s <- 100*epi$inew/population
inew_lb_s <- 100*epi$inew_lb/(population * 0.1)
inew_r_s <- 100*epi$inew_r/(population * 0.001)

##  plot the daily infections trajectories in three samples (one simulate)
par(mfcol=c(2,2),mar=c(2,2.5,4,2)) 
plot(inew_s,ylim=c(0,max(inew_s,inew_lb_s,inew_r_s)),xlim=c(0,120),xlab="Day",ylab="Daily new infections per 100", type = 'l', col=1) 
## pop daily new infections (black)
lines(inew_lb_s,col=4) ## cautious 10% new daily infections (blue)
lines(inew_r_s,col='brown') ## 0.1% random sample new daily infections
legend("topleft", legend = c("whole poplation","cautious 10%","0.1% random sample"), 
       lwd = 4, col = c(1,4,'brown'), cex = 0.4, bty = "n")
title("New daily infections", cex.main=0.9)

## label the day on which each trajectory peaks
text(which.max((inew_s)),max(inew_s), labels = paste("Peak on day", which.max((inew_s))), col=1, adj=-0.25, cex = 0.5)
text(which.max((inew_lb_s)),max(inew_lb_s), labels = paste("Peak on day", which.max((inew_lb_s))), col=4, adj=-0.25, cex = 0.5)
text(which.max((inew_r_s)),max(inew_r_s), labels = paste("Peak on day", which.max((inew_r_s))), col='brown', adj=-0.05, cex = 0.5)



##  Step3
##  Running 10 replicate simulations and plot results

simulate_10 <- lapply(rep(5.5e6,10), epidemic) ## running 10 replicate simulations and return a list

##  plot daily infections trajectories in three samples respectively
##  plot 10 waves of the whole pop
plot(simulate_10[[1]]$inew,xlab="day",ylab="Daily new infections", ylim=c(0, max(sapply(simulate_10,function(x) max(x$inew)))),
     type = 'l', col=1)
for (i in 2:10){
  lines(simulate_10[[i]]$inew)
}
title("10 times simulation of\n the whole poplation",cex.main=0.9)
##  plot 10 waves of the 10% cautious pop
plot(simulate_10[[1]]$inew_lb,xlab="day",ylab="daily new infections",ylim=c(0, max(sapply(simulate_10,function(x) max(x$inew_lb)))),
     type = 'l', col=4)
for (i in 2:10){
  lines(simulate_10[[i]]$inew_lb, col=4)
}
title(strwrap("10 times simulation of the 10% of the population with lowest beta",width = 30),cex.main=0.9)
##  plot 10 waves of a random sample of 0.1% whole pop
plot(simulate_10[[1]]$inew_r,xlab="day",ylab="daily new infections", ylim=c(0, max(sapply(simulate_10,function(x) max(x$inew_r)))), 
     type = 'l', col='brown')
for (i in 2:10){
  lines(simulate_10[[i]]$inew_r, col='brown')
}
title(strwrap("10 times simulation of a random sample of 0.1% of the population",width = 30),cex.main=0.9)


##  Step4
##  Analyze the implications of these results interpreting daily 
##  infection trajectories reconstructed using the ZOE app data.

##   The results show that：

#### 1 The daily infection peak reconstructed using the ZOE app data are likely to be later than others. 
##     The infection peak of 10% cautious population with low beta values is about 
##     2-5 days later than the peak of whole population and the peak of 0.1% random sample.

#### 2 The number of infections per 100 people in the population reconstructed from the ZOE app data was 
##     significantly smaller than in the total population. 
##     The infection curve for the 10% cautious population with low beta valiue is much lower than 
##     the infection curves for the whole population and the 0.1% random sample (especially near the infection peak).



