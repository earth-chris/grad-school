#first set directory and read in necessary functions and data
setwd('~/teaching/es211/2016/hw/hw3_calibrationAssignment_2015/')
source('RRmodels.R')
load('runoff.Rdata') # loads P, PET, and runoff

#first you might like to see how the bucket1 and bucket3 models differ
#run them both and plot together
partrue = c(160,.53,.26)
plot(bucket1(partrue,P,PET),col=1,type='l')
lines(bucket3(partrue,P,PET),col=2,type='l')

# define various loss functions
loss.rmse = function(y,yhat) sqrt(mean((y-yhat)^2,na.rm=T))
loss.abs = function(y,yhat) mean(abs(y-yhat),na.rm=T)

# define a function that calculates the objective value, given a loss function and a model
# make sure you understand what this beauty is doing!
minfunc = function(pars,loss,model,y,P,PET) loss(y, model(pars,P,PET))

#let's first try bucket1. 
#define "truth"
partrue = 160
ytrue = bucket1(partrue,P,PET)
par0 = 100
opt = optim(par0,minfunc,loss=loss.rmse,model=bucket1,y=ytrue,P=P,PET=PET,method='L-BFGS-B')
print(opt$par)

#try for different starting values
for (par0 in seq(0,200,by=50)){
  opt = optim(par0,minfunc,loss=loss.rmse,model=bucket1,y=ytrue,P=P,PET=PET,method='L-BFGS-B')
  print(opt$par)
}
#you should see that it always recovers true value, even for very different starting points

#now let's try for bucket3, which takes 3 parameters
partrue = c(120,.53,.26)
ytrue = bucket3(partrue,P,PET)

par0 = c(100,0,0)
opt = optim(par0,minfunc,loss=loss.rmse,model=bucket3,y=ytrue,P=P,PET=PET,method='L-BFGS-B')
print(rbind(true=partrue,est=opt$par))
#compare the "true" values with outputs from the calibrated model
matplot(cbind(ytrue, bucket3(opt$par,P,PET)), type='s', ylab="mm",xlab="months", main=paste("pars =",paste(round(opt$par,2),collapse=", "))); legend('topleft',leg=c('obs','model'),fill=1:2)

#you can see that the parameter values were quite
#different, but the runoff values were pretty close

#let's try some different starting points
par0=c(150,0,0)
opt = optim(par0,minfunc,loss=loss.rmse,model=bucket3,y=ytrue,P=P,PET=PET,method='L-BFGS-B')
print(rbind(true=partrue,est=opt$par))

#these end up much closer, which says we are somewhat sensitive to starting
#value. this happens because lots of parameter values give pretty good fits

#let's see if it helps at all if we restrict the search
lb = c(10,0,0) # lower parameter bounds
ub = c(200,.5,.5) # upper parameter bounds
par0 = c(100,0,0)
#no constraints
opt = optim(par0,minfunc,loss=loss.rmse,model=bucket3,y=ytrue,P=P,PET=PET,method='L-BFGS-B')
print(rbind(true=partrue,est=opt$par))
#with constraints
opt = optim(par0,minfunc,loss=loss.rmse,model=bucket3,y=ytrue,P=P,PET=PET,method='L-BFGS-B',lower=lb,upper=ub)
print(rbind(true=partrue,est=opt$par))
#does a bit better with constraints

#try different search methods
par0 = c(100,0,0)
#
for (method in c('L-BFGS-B','BFGS','Nelder-Mead')){
  opt = optim(par0,minfunc,loss=loss.rmse,model=bucket3,y=ytrue,P=P,PET=PET,method=method)
  print(paste('method:',method))
  print(rbind(true=partrue,est=opt$par))
}
#also matters a bit

# overall, this exercise should help you understand:
# 1) how to calibrate using optim()
# 2) how the results for some models will be robust (e.g. bucket1) but
# for other models (e.g. bucket3) the parameters may be hard to find
# because different parameter combinations can give very similar outputs
# 3) the way the search is done will also affect results, such as what
# constraints are used and the gradient search method

#now, on to the takehome portion. it involves data we generated
#with secret values of the parameters.

