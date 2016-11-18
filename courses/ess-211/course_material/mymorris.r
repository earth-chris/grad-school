# Assignment 3 - Example Code
#QUESTION 1
base.samp = function(p,n) {   #DEFINE FUNCTION TO SAMPLE FROM SET {0,1/(p-1),...1}
  x = (p-1-p/2) * runif(n)    #NEED TO LIMIT TO RANGE OF (0, 1-delta)
  round(x) / (p-1)
}

mymorris = function(k,r,p,delta,par.mins,par.maxs,fun) {  #NOTE THE LAST INPUT IS THE MODEL NAME FOR SENS ANALYSIS

  par.range = par.maxs - par.mins
  effects = array(dim = c(k,r))

  for (r.ind in 1:r) {
    J = array(1,dim=c(k+1,k)) #AN ARRAY OF ONES
    B = lower.tri(J) * 1    #A LOWER TRIANGULAR ARRAY OF ONES
    xstar = base.samp(p,k)    #BASE VECTOR
    D = array(0,dim=c(k,k))   #AN ARRAY WITH EITHER 1 OR -1 IN DIAGONAL
    diag(D) = 1.0 -2*(runif(k) < .5)
    P = array(0,dim=c(k,k))   #A RANDOM PERMUTATION ARRAY
    diag(P) = 1
    P = P[,sample(c(1:k),k,replace=F)]
    Bstar =  (t(xstar*t(J)) + (delta/2) * ((2*B - J)%*%D + J)) %*% P

    y = numeric(length=(k+1)) #VECTOR TO HOLD MODEL OUTPUT
    for (i in 1:(k+1)) y[i] = fun(par.mins + par.range * Bstar[i,] )
    for (i in 1:k) {
      par.change = Bstar[i+1,] - Bstar[i,]  #FIND WHICH PARAMETER CHANGED AND WHETHER UP OR DOWN
      i2 = which(par.change != 0)
      effects[i2,r.ind] = (y[i+1] - y[i]) * (1 - 2 * (par.change[i2] < 0))
    }
  }
  effects   #RETURN THE k x r ARRAY OF COMPUTED EFFECTS
}


