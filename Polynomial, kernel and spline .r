
library(KernSmooth)
#---------------------------------------------
# My Local Polynomial Smoother Bootstrap Code 
#---------------------------------------------
# Fit a polynomial smoothers and calculates a symmetric nonparametric bootstrap confidence intervals.
# Arguments:
#    x, y       : data values
#    nreps      : number of bootstrap replicates

require(KernSmooth)
lpsboot <- function(x,y,nreps=1000, band=5, confidence = 0.95){
  # Put input data into a data frame, sorted by x, with no missing values.
  dat <- na.omit(data.frame(x=x,y=y))
  if(nrow(dat) == 0) {
    print("Error:l No data left after dropping NAs")
    print(dat)
    return(NULL)
  }
  ndx <- order(dat$x)
  dat$x <- dat$x[ndx]
  dat$y <- dat$y[ndx]
  # Fit curve to data
  len <- length(dat$x)
  f0 <- locpoly(x, y, kernel = "epanechnikov", bandwidth = band, gridsize = len)
  y.fit <- f0$y
  # Generate bootstrap replicates
  mat <- matrix(0,NROW(dat), nreps)
  for(i in seq(nreps)){
    ndx <- sample(len,replace=T)
    x.repl <- x[ndx]
    y.repl <- y[ndx]
    f <- locpoly(x.repl, y.repl, kernel = "epanechnikov", bandwidth =  5, gridsize = len) #change kernel to "normal"
    mat[, i] <- f$y
  }
  # calculating confidence intervals
  ci <- t(apply(mat, 1, quantile, probs = c((1-confidence)/2, (1+confidence)/2)))
  res <- cbind(as.data.frame(f0), ci)
  colnames(res) <- c('x','y', 'lwr.limit','upr.limit')
  res
}
#---------------------------------------------
# My Kernel Smoother Bootstrap Code 
#---------------------------------------------
# Fit a kernel smoothers and calculates a symmetric nonparametric bootstrap confidence intervals.
# Arguments:
#    x, y       : data values
#    nreps      : number of bootstrap replicates

ksboot <- function(x,y,nreps=1000, band=5, confidence = 0.95){
  # Put input data into a data frame, sorted by x, with no missing values.
  dat <- na.omit(data.frame(x=x,y=y))
  if(nrow(dat) == 0) {
    print("Error: No data left after dropping NAs")
    print(dat)
    return(NULL)
  }
  ndx <- order(dat$x)
  dat$x <- dat$x[ndx]
  dat$y <- dat$y[ndx]
  # Fit curve to data
  require(KernSmooth)
  len <- length(dat$x)
  f0 <- ksmooth(x, y, kernel = "normal", bandwidth = band, n.points = len)
  y.fit <- f0$y
  # Generate bootstrap replicates
  mat <- matrix(0,NROW(dat), nreps)
  for(i in seq(nreps)){
    ndx <- sample(len,replace=T)
    x.repl <- x[ndx]
    y.repl <- y[ndx]
    f <- ksmooth(x.repl, y.repl, kernel = "normal", bandwidth =  5, n.points = len)
    mat[, i] <- f$y
  }
  # calculating confidence intervals
  ci <- t(apply(mat, 1, quantile, probs = c((1-confidence)/2, (1+confidence)/2)))
  res <- cbind(as.data.frame(f0), ci)
  colnames(res) <- c('x','y', 'lwr.limit','upr.limit')
  res
}
#---------------------------------------------
# My Spline Smoother Bootstrap Code 
#---------------------------------------------
splineboot<-function(x,y,nreps=1000,confidence=0.95,spar=0.75){
  dat<-na.omit(data.frame(x=x,y=y))
  # original data 
  f0<- smooth.spline(x=dat$x, y = dat$y,spar=spar)
  len<-length(dat$x)
  # store replicates
  N <- 500
  mat <- matrix(0,N, nreps)
  new_x<-seq(min(x),max(x), length.out = N)
  y.fit<-predict(f0,new_x)
  # bootstraps
  for (i in 1:nreps){
    ndx<-sample(len,replace=T)
    x.repl<-x[ndx]
    y.repl<-y[ndx]
    f<-smooth.spline(x=x.repl, y = y.repl,spar=spar)
    mat[,i]<-predict(f,new_x)$y
  }
  ci <- t(apply(mat, 1, quantile, probs = c((1-.95)/2, (1+.95)/2)))
  res <- cbind(as.data.frame(y.fit), ci)
  colnames(res) <- c('x','y', 'lwr.limit','upr.limit')
  res
}

# Load your own data
data <- read.csv("/Users/shedrachezihe/Downloads/c1_bdhs.csv")
d1 <- na.omit(data)

# Replace 'x' and 'y' with your variables
x <- d1$AOR
y <- d1$RBMI

d2 <- data.frame(x, y)


# Step 2: Take a random sample of 100 rows
set.seed(8912)
d3 <- d2[sample(nrow(d2), 100), ]
d3

# local polynomial

lp_data <- with(d3, lpsboot(x, y, nreps = 5000))
with(d3, plot(x, y, las = 1,main = "Local Polynomial smoothing"))
with(lp_data, matpoints(x, lp_data[, - 1], type = 'l', col = c(1, 2, 2), lty = c(1, 2, 2)))


# kernel smoothing

ks_data <- with(d3, ksboot(x, y, nreps = 5000))
with(d3, plot(x, y, las = 1,main="Kernel Smoothing"))
with(ks_data, matpoints(x, ks_data[, - 1], type = 'l', col = c(1, 2, 2), lty = c(1, 2, 2)))

# spline smoothing

spline_data <- with(d3, splineboot(x, y, nreps = 5000))
with(d3, plot(x, y, las = 1, main="Spline smoothing"))
with(spline_data, matpoints(x, spline_data[, - 1], type = 'l', col = c(1, 2, 2), lty = c(1, 2, 2))) 

# Step 6: Plot the three smoothing curves

# Set up the 3 by 1 grid
par(mfrow=c(3, 1))

lp_data <- with(d3, lpsboot(x, y, nreps = 5000))
with(d3, plot(x, y, las = 1,main = "Local Polynomial smoothing"))
with(lp_data, matpoints(x, lp_data[, - 1], type = 'l', col = c(1, 2, 2), lty = c(1, 2, 2)))

ks_data <- with(d3, ksboot(x, y, nreps = 5000))
with(d3, plot(x, y, las = 1,main="Kernel Smoothing"))
with(ks_data, matpoints(x, ks_data[, - 1], type = 'l', col = c(1, 2, 2), lty = c(1, 2, 2)))

spline_data <- with(d3, splineboot(x, y, nreps = 5000))
with(d3, plot(x, y, las = 1, main="Spline smoothing"))
with(spline_data, matpoints(x, spline_data[, - 1], type = 'l', col = c(1, 2, 2), lty = c(1, 2, 2))) 

# Create data frames with the results
lp_values <- data.frame(x = lp_data$x, lp_smooth_values = lp_data$y)
ks_values <- data.frame(x = ks_data$x, ks_smooth_values = ks_data$y)
spline_values <- data.frame(x = spline_data$x, spline_smooth_values = spline_data$y)

# Combine the results into a single data frame
smoothed_values <- data.frame(x = d3$x, lp_values, ks_values, spline_values)

# Display the first 30 rows of the table
head(smoothed_values, 30)

