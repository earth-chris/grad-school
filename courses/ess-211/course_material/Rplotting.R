## -------------- X-Y PLOTS --------------

# --- sample data ---
data(airquality)   # loads the airquality dataset into the workspace
attach(airquality) # extracts each variable as a vector



## -------------- PLOT FUNCTIONS --------------
# --- point plot ---
plot(Temp, Ozone)



# --- adding labels ---
plot(Temp, Ozone, main="Ozone vs. Temperature in NY, 1973", 
     xlab="Temperature (degrees F)", ylab="Ozone (ppb)")



# --- plotting one set of points ---
plot(Temp, main="Temp plotted against its indices (i.e. days)", 
     xlab="day")



# --- plotting lines ---
plot(Temp, type="l")  # "l" for lines
plot(Temp, type="s")  # "s" for stair steps



plot(Temp, Ozone, type="l")



# --- point types ---
x<-rep(seq(1,5),5) 
y<-sort(x,decreasing=TRUE) 
pch<-seq(1,25) 
plot(x,y,pch=pch,cex=2,xlim=c(1,5.4), axes=FALSE,xlab="R symbols",ylab="")
text(x+0.25,y,'pch')



# --- color and size ---
plot(Temp, Ozone, pch=17, col="blue", cex=1.5)
plot(Temp, Ozone, main="Default Sizes", pch=17, col='blue')
plot(Temp, Ozone, col='blue', main="Inflated Sizes", pch=17, cex=1.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, las=1)
plot(Temp, type='l', lwd=4,  main="Thicker line")



## -------------- ADDING POINTS OR LINES TO AN EXISTING PLOT --------------
plot(Ozone, type="l", xlab="", ylab="")
plot(Ozone, type="l", xlab="", ylab="")
points(Temp, col='blue')           # adds points
points(Temp, col='blue', type='l') # adds a line through these points
lines(Temp, col='blue')            # a shortcut for points(Temp, type='l'); does exactly the same thing



plot(Ozone, type="l", xlab="", ylab="")
points(Temp, col='blue')           # adds points
points(Temp, col='blue', type='l') # adds a line through these points
lines(Temp, col='blue')            # a shortcut for points(Temp, type='l'); does exactly the same thing
lines(Solar.R, col='red')



range(Temp)                        # gives the minimum and maximum values
range(Solar.R, na.rm=T)            # why do we need the na.rm=T argument here?
range(Ozone, na.rm=T)
# It looks like we need our plot window to range between 1 and 334 if we're to cover everything
# We'll also change the x-limit, just to show its effect



plot(Ozone, ylab="", type="l", xlim=c(1,170), ylim=c(1,334)) 
lines(Temp, col='blue')
lines(Solar.R, col='red')



## -------------- ADDING A LEGEND --------------
plot(Ozone, xlab="days", ylab="", type="l", xlim=c(1,170), ylim=c(1,334)) 
lines(Temp, col='blue')
lines(Solar.R, col='red')
legend('topright',legend=c('ozone (ppb)','temp (F)','rad (lang)'),col=c('black','blue','red'), lty='solid')



plot(Ozone, ylab="", type="l", xlim=c(1,170), ylim=c(1,334)) 
lines(Temp, col='blue')
lines(Solar.R, col='red')
# add a legend at x=50, y=250
# note: this isn't a brilliant place for it; just showing how it works
# bg="white" forces a white (as opposed to transparent) legend background
legend(50,250,legend=c('ozone (ppb)','temp (F)','rad (lang)'), col=c('black','blue','red'), lty='dashed', bg="white")



## -------------- MULTIPLE LINES ON ONE PLOT --------------
toplot = matrix(1:12, 4) + rnorm(12,sd=3)



toplot = t(matrix(1:12, 4)) + rnorm(12,sd=3)
matplot(toplot,  lty=1, lwd=2,type='l')     # Each column plotted as a separate line



## -------------- ADDING TEXT --------------
# Evaluate these lines one-by-one so you understand what each is doing!
value           = mean(Temp,na.rm=T)        # calculate the mean
rounded_value   = round(value, 1)           # round to 1 decimal place, for display purposes
text_to_display = paste("Mean Temp =", rounded_value)



# now plot, and add text_to_display at desired x,y coords
plot(Temp, type='l', ylab="degrees F", xlab="Days", main="Temperature")
text(25,95,text_to_display)



## -------------- VECTOR ARGUMENTS TO CEX, COL --------------
colors       = rep('blue',length(Temp))   # all elements are 'blue'
odds         = seq(1, length(Temp), by=2) # the odd indices: c(1,3,5, ...)
colors[odds] = 'red'
# check that "colors" has alternating "red" and "blue"
plot(Temp, col=colors, pch=19)



# help them with the more straightforward way
# e.g. colors[Month==5] = 'blue'; colors[Month==6] = 'red', etc.

n                = length(Temp)

colors           = rep(NA,n)
palette          = c('blue','red','green','brown','purple')
colors[Month==5] = palette[1]
colors[Month==6] = palette[2]
colors[Month==7] = palette[3]
colors[Month==8] = palette[4]
colors[Month==9] = palette[5]
## 
## sizes = rep(NA,n)
## sizes[Solar.R<100] = 1
## sizes[Solar.R>=200 & Solar.R<300] = 2
## sizes[Solar.R>300] = 3
## 
##plot(Temp,Ozone,col=colors,pch=19)
## legend(55,170, legend=unique(Month), fill=palette, title='Month')
## legend(63,170, pch=19, legend=c("<100","100 - 200",">300"), pt.cex=1:3, title="Solar.R")
## 
## 
## # then try it a little more programatically
## cols = rep(NA,n)
## months = unique(Month)
## for (i in 1:length(months)){
##   cols[Month==months[i]] = palette[i]
## }
##
# or better yet
## nmonths = length(unique(Month))
## library(RColorBrewer)
## palette = brewer.pal(nmonths, 'Set2')
## colors = as.character(cut(Month, nmonths, labels=palette))
## plot(Temp,Ozone,pch=19,col=colors,cex=sizes)
## legend(55,170, legend=unique(Month), fill=palette, title='Month')
## legend(63,170, pch=19, legend=c("<100","100 - 200",">300"), pt.cex=1:3, title="Solar.R")
## 
## 
## dates          = as.Date(paste("1973",Month,Day, sep="-"))
## dow            = format(dates, "%a")
## weekends       = which(dow %in% c('Sat','Sun'))
## cols           = rep('red',length(Temp))
## cols[weekends] = 'blue'
## library(scales)
## sizes          = scales::rescale(Solar.R,c(1,3))
## plot(Temp,Ozone, col=cols, cex=sizes)
## legend(60,170,legend=c('weekends','weekdays'),fill=c('blue','red'))
## leg.vals       = c(100,200,300)
## leg.sizes      = rescale(leg.vals, to=c(1,3), from=range(Solar.R,na.rm=T))
## legend(60,120,legend=leg.vals,pt.cex=leg.sizes, pch=19, title='Solar.R')
## 
  


## -------------- HISTOGRAM PLOTS --------------
hist(Temp, main="Default number of bins based on number of values")



hist(Temp, breaks=50, main="50 bins")



hist(Temp, breaks=c(0,60,85,100), main="Explicitly-set breaks")
# Be sure you understand exactly what this is doing!
# This is the sort of task you'll have to do a lot.



## -------------- BOXPLOTS --------------
boxplot(Ozone ~ Month, main="Boxplot of Ozone values",xlab='Month')



## -------------- BAR PLOTS --------------
# Bad way (this would be awful if we had, say, 100 elements to deal with instead of just 5)
monthly_averages    = c() # an empty vector
monthly_averages[1] = mean(Temp[Month==5])
monthly_averages[2] = mean(Temp[Month==6])
monthly_averages[3] = mean(Temp[Month==7])
monthly_averages[4] = mean(Temp[Month==8])
monthly_averages[5] = mean(Temp[Month==9])



# Better way
monthly_averages = c() # an empty vector
for (i in 5:9) {
  monthly_averages = c(monthly_averages, mean(Temp[Month==i]))
}



# Even better (why?)
months           = unique(Month)
monthly_averages = rep(NA, length(months))
for (i in 1:length(months)) {
  monthly_averages[i] = mean(Temp[Month==months[i]])
}



# Best way (stay tuned ...)
monthly_averages = sapply(unique(Month), function(x) mean(Temp[Month==x]))



# Plot monthly_averages, no matter how you obtained it
plot(monthly_averages, pch=19, xlab="Month", ylab="degrees T", main="Temperature")



barplot(monthly_averages, 
        names.arg = c('May','June','July','Aug','Sep'), 
        ylab="degrees F", main="Monthly Average Temperature")



a <- barplot(monthly_averages, names.arg = c('May','June','July','Aug','Sep'), ylab="degrees F", main="Monthly Average Temperature")
#now we can add some text easily
text(a, monthly_averages,round(monthly_averages,1),pos=1)



# --- Grouped bar plots ---


# Note a couple new tricks introduced here:
# 1) Periods are legal character in variable names.
# 2) You can assign two variables to equal the same value in one line.



months   = unique(Month)
temp.avg = ozone.avg = rep(NA, length(months)) 
for (i in 1:length(months)) {
  temp.avg[i] = mean(Temp[Month==months[i]])
  ozone.avg[i] = mean(Ozone[Month==months[i]], na.rm=T) 
  # Why did we need na.rm=T for Ozone, but not Temp?
}
avgs = cbind(temp.avg, ozone.avg) # bind temp.avg and ozone.avg as columns



# barplot groups on columns, so temp.avg and ozone.avg each get their own bar
barplot(avgs)
barplot(t(avgs)) # after transposing, each month is its own group



# beside=T unstacks the bars and places them side-by-side
# barplot() knows that if there are groups, we'll need a legend to specify which bar is which.
# So, it provide an automatic-legend generator, andall we have so supply is legend.text.
barplot(t(avgs), names.arg = c('May','June','July','Aug','Sep'), beside=T, legend.text=c('Temp','Ozone'))



## -------------- MULTI PANEL PLOTS --------------

# --- Regular grids ---
par(mfrow=c(2,2)) # make a 2x2 grid in which to place plots
hist(Ozone)
hist(Temp)
hist(Solar.R)
plot(Temp,Ozone)



# --- Irregular grids ---
layout(matrix(c(1,1,2,3), 2, 2, byrow=T))
plot(Temp,Ozone)
hist(Temp)
hist(Ozone)



# --- Two lines with different scales ---
plot(Ozone, ylab="", type="l", xlim=c(1,170)) 
par(new=T)
plot(Solar.R, col='red',type='l',axes=F,xlab='',ylab='', xlim=c(1,170))
axis(4,col.axis='red')
