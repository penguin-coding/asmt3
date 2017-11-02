x = seq(1,15)
y = c(4.39, 9.32, 13.77, 18.28, 23.48,   # Some observed times for
      31.85,32.38, 36.75, 40.69, 55,     # calculation times
      52.75, 58.27, 60.5, 62.58, 75.11)

A = lm(y~x-1)  # We fit a lm with no intercept, and check
summary(A)     # it produces a reasonable answer. 

plot(x,y)                    # We visualise the fitted line
abline(0,A$coefficients)

pred.time <- function(s){     # We write a very simple function
  secs <- 2*A$coefficients*s  # which returns the expected
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