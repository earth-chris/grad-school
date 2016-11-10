# ESS-211 assignment 5, due 11-11-16
#  cba 11/2016

#############################
# below are the functions provided by the teaching group

# define the linear oscillator equation for determining bungee height
hmin <- function (x){
  return(x[1] - 2 * x[2] * 9.8 / (1.5 * x[3]))
}

#############################
# below are the exercises provided by the teaching group

###
# first, local SA

# define bungee names and parameters
par.names <- c("Height", "Mass", "# Strands")
par.mins <- c(40 , 67 , 20)
par.maxs <- c(60 , 74 , 40)

# define midpoints for parameters
par.mids <- 0.5 * (par.mins + par.maxs)

# set n parameters
npars <- length (par.mins)

# create array to hold results
difs <- numeric(length = npars)

# loop through each parameter set and
for (par in 1:npars){

  # create two input vectors to assign min/max values
  input1 <- par.mids ; input2 <- par.mids
  
  # assign min value to first vector
  input1[par] <- par.mins[par]
  
  # assign max value to second vector
  input2[par] <- par.maxs[par]
  
  # compute differences between model outputs using different min/max vals
  difs[par] = hmin(input2) - hmin(input1)
}

# plot the output
barplot(difs, xlab = "Parameter", ylab = "Max - Min", names.arg = par.names)

###
# next, one-at-a-time global SA

# define n parameters
k <- npars

# n times to compute effect for each parameter
r <- 10

# n possible levels for each parameter (should be an even #)
p <- 4

# the increment to adjust parameter values by
delta <- p / (2 * (p - 1))

# define function to sample from set {0, 1/(p -1), ... 1}
base.samp <- function (p , n){
  
  # limit to range of (0, 1-delta)
  x <- (p - 1 - p / 2) * runif(n)
  return(round(x) / (p -1))
}

# define range for parameters
par.mins <- c(40, 67, 20)
par.maxs <- c(60, 74, 40)

# create the morris function
morris <- function(k, r, p, delta, par.mins, par.maxs, fun){
  
  # set the range for parameters
  par.range <- par.maxs - par.mins
  
  # create an array to save the effects
  effects <- array(dim = c(k, r))
  
  # loop through each computed effect
  for (r.ind in 1:r){
    J = array (1 , dim = c ( k +1 , k ) ) #AN ARRAY OF ONES
    B = lower . tri ( J ) * 1 #A LOWER TRIANGULAR ARRAY OF ONES
    xstar = base . samp (p , k ) # BASE VECTOR
    D = array (0 , dim = c (k , k ) ) #AN ARRAY WITH EITHER 1 OR -1 IN DIAGONAL
    diag ( D ) = 1.0 -2 * ( runif ( k ) < .5)
    P = array (0 , dim = c (k , k ) ) #A RANDOM PERMUTATION ARRAY
    diag ( P ) = 1
    P = P [ , sample ( c (1: k ) ,k , replace = F ) ]
    Bstar = ( t ( xstar * t ( J ) ) + ( delta / 2) * ((2 * B - J ) % * % D + J ) ) % * % P
    y = numeric ( length =( k +1) ) # VECTOR TO HOLD MODEL OUTPUT
    for ( i in 1:( k +1) ) y [ i ] = fun ( par . mins + par . range * Bstar [i ,] )
    for ( i in 1: k ) {
      par . change = Bstar [ i +1 ,] - Bstar [i ,] # FIND WHICH PARAMETER CHANGED AND WHETHER UP OR DOWN
      i2 = which ( par . change ! = 0)
      effects [ i2 , r . ind ] = ( y [ i +1] - y [ i ]) * (1 - 2 * ( par . change [ i2 ] < 0) )
    }
  }
  return(effects) # RETURN THE k x r ARRAY OF COMPUTED EFFECTS
}
#NOW RUN THE MORRIS METHOD FOR hmin
effects = morris (k ,r ,p , delta , par . mins , par . maxs , hmin )
mu = apply ( effects ,1 , mean )
sigma = apply ( effects ,1 , sd )
mu . star = apply ( abs ( effects ) ,1 , mean )
plot ( mu . star , sigma , pch = c ( " H " ," M " ," s " ) )