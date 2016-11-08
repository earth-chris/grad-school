vsa=function(fun,par.mins,par.maxs,nrun){
  par.range = par.maxs - par.mins
  k=length(par.range)
  M = array(runif(nrun*k),dim=c(nrun,k))      #CREATE U(0,1) n x k MATRIX
  for (i in 1:k) M[,i] = M[,i] * par.range[i] + par.mins[i] #TRANSFORM TO FIT RANGE OF EACH PAR

  #Now we are going to make a different input matrix, M2, the same way
  M2 = array(runif(nrun*k),dim=c(nrun,k))     #CREATE U(0,1) n x k MATRIX
  for (i in 1:k) M2[,i] = M2[,i] * par.range[i] + par.mins[i] #TRANSFORM TO FIT RANGE OF EACH PAR

  #Now we are going to make K different matrices, where in the jth matrix, all parameters are taken from M2 except parameter j is taken from M
  NJ.list = list()
  for (j in 1:k) {
    temp = M2
    temp[,j] = M[,j]
    NJ.list[[j]] = temp
  }

  #Now, to compute the first-order effects (S) of factor J, we need to compute the values UJ, as discussed in class
  S = numeric(length=k)
  ST = numeric(length=k)
  y1 = y2 = y3 = numeric(length = nrun)
  for (i in 1:nrun) {
      y1[i] = fun(M[i,])
      y2[i] = fun(M2[i,])
  }
 #Now to compute the expected value (EY) and total variance (VY) of Y, for which we will use M2
  EY=mean(y1)
 # VY=var(y3)
 #For computational reasons, we will estimate EY^2 using both M and M2
 # EY2=mean(y1*y2)
#  EY2=EY^2
  VY=sum(y1*y1)/(nrun-1)-EY^2

  for (j in 1:k) {
    NJ = NJ.list[[j]]
    for (i in 1:nrun) y3[i] = fun(NJ[i,])
    UJ = sum(y1*y3) / (nrun-1)   #here everything but factor j is resampled
    UJ2 = sum(y2*y3) / (nrun-1)  #here only factor j is resampled
    S[j] = (UJ - EY^2) / VY
    ST[j] = 1.0 - (UJ2 - EY^2) / VY
  }

  list(S,ST)
}