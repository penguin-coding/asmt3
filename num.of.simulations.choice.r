#### NEEDS TO BE RUN ON A CLEAN ENVIRONMENT AND ENVIRONMENT SHOULD BE CLEANED
#### BEFORE RUNNING FULL SIMULATION

### Note: running this code on any other machines will not produce the same
### predicted times, refer to the figure in the report for the specific
### graph this produced on my machine. 

x <- seq(3)  # Setting this to 4 takes 2 minutes on my machine, but 3 is enough
y <- rep(NA, length(x))

for (j in x){
  # Our x values are the number of different simulations we run for each case.
  # Our y values is the computation time required for the respective x.
  # We will attempt to draw a graph of simulations per case and computation 
  # time, so that we can select a compromise between high simulations and 
  # available time.
  
  pilot.sims <- x[j]
  y[j] <- system.time(source('driver.R'))[3]
}

A = lm(y~x-1)    # We fit a lm and check it
summary(A)       # produces a reasonable answer. 

plot(x,y)                    # We visualise the fitted line
abline(0,A$coefficients)

pred.time <- function(s){     # We write a very simple function
  secs <- A$coefficients*s    # which returns the expected
  mins <- secs/60             # calculation time for s simulations
  hours <- mins/60            # per case, based on our observed
  return(hours)               # values
}

plot.x = seq(1,2000)        # we plot expected calculation times
plot.y = pred.time(plot.x)  # so that we can choose how many to use

plot(plot.x, plot.y, xlab = 'number of simulations per case',
     ylab = 'Time for all simulations to run - hours',
     type='l',
     main = 'Predicted calculation times for 1.51GHz CPU')