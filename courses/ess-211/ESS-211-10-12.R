# an introduction to plotting

# load sample data
data(airquality)
attach(airquality)

# plot the data
plot(Temp, Ozone)

# add labels
plot(Temp, Ozone, main-"Ozone vs Temperature in NY, 1973", xlab="day")

# plot one set of points
plot(Temp, main="Temp plotted against its indices")

# plot with lines or with steps
plot(Temp, type="l")
plot(Temp, type="s")

# if you plot a scatter plot as a line plot, it will show up as spaghetti since it is trying to connect the dots
#  you would want to sort the data first
plot(Temp, Ozone, type="l")

# pch is the plot character - the default is a circle, but there are many different symbols to use (see ?pch)
plot(Temp,Ozone, pch=19)

# can add color, too
plot(Temp, Ozone, pch=19, col='magenta')

# and can change size, too
plot(Temp, Ozone, pch=19, col='magenta', cex=2)

# can change other params, too
plot(Temp, Ozone, pch=19, col='magenta', cex=1.5, cex.axis=1.5)

# and change y-labels to horizontal using las=1
plot(Temp, Ozone, pch=19, col='magenta', cex=1.5, cex.axis=1.5, las=1)

# and chanage line thickness
plot(Temp, type='l')
plot(Temp, type='l', lwd=3)

# and we can overlay plots
plot(Ozone, type='l', xlab='', ylab='')
points(Temp, col='blue') 

# to add another line, could do points, with type='l', or use lines()
lines(Temp, col='blue')

# can add solar radiation, but it'll go off plot
lines(Solar.R, col='red')

# want to find min/max for each object you're going to plot
range(Temp)
range(Solar.R, na.rm = True)
range(Ozone, na.rm = True)

# you'll want to do this programatically, but hardcode it here
#  you don't really want to do this, since these are all in different units
plot(Ozone, ylab='', type='l', xlim=c(1,170), ylim=c(1,334))
lines(Temp, col='blue')
lines(Solar.R, col='red')

# we can create legends, too
legend('topright', legend=c('ozone (ppb)', 'temp (F)', 'solar radiation (lang)'),
		col=c('black','blue','red'))

# you can specify the x and y location of the legend layout using a vector of x and y position
legend(50,250, legend=c('ozone (ppb)', 'temp (F)', 'solar radiation (lang)'),
		col=c('black','blue','red'))

# we can plot multiple items at once. first, create some matrices
toplot = matrix(1:12, 4) + rnorm(12, sd=3)
matplot(toplot, lty=1) # didnt catch full function from eva

# to add text
value = mean(Temp, na.rm=T) #calculate the mean
rounded_value = round(value,1) # round to one decimel place
text_to_display = paste("Mean Temp =", rounded_value)

# now plot and add text_to_display at desired x/y coords
plot(Temp, type='l', ylab='degrees F', xlab='days', main='Temperature')
text(25,95, text_to_display)

# we can plot unique colors per element in a vector
colors = rep('blue', length(Temp)) # make sure it is same length
odds = seq(1, length(Temp), by=2) # set an odd-numbered index to set every other color as red/blue
colors[odds] = 'red'
plot(Temp, col=colors, pch=19)

# work on an exercise
plot(Ozone, Temp, pch=19)

# color by month
# set a vector of colors for each month
cvec = c("dark red", "red", "burnt orange", "orange", "yellow", "green", "turquoise", "blue", "purple", "magenta", "burgundy")

# set a color vector for each element of the temps
colors = rep(NA, length(Temp))

# guess that was unnecessary, as only has data for may-sept
colors[Month==5] = cvec[1]
colors[Month==6] = cvec[2]
colors[Month==7] = cvec[3]
colors[Month==8] = cvec[4]
colors[Month==9] = cvec[5]

# plot it!
plot(Temp, col=colors, pch=19)

# histogram plots
hist(Temp, main='default # bins based on number of values')

# can set the number of bins
hist(Temp, breaks=50, main='50 bins')

# can set a vector for where you want breaks to occur
hist(Temp, breaks=c(0,60,85,100))

# lets do some box plots!
boxplot(Ozone ~ Month, main='Boxplot of Ozone values', xlab='Month')

# or bar plots!
# below is a bad way to do it
monthly_averages = c()
monthly_averages[1] = mean(Temp[Month==5])
monthly_averages[2] = mean(Temp[Month==6])
monthly_averages[3] = mean(Temp[Month==7])
monthly_averages[4] = mean(Temp[Month==8])
monthly_averages[5] = mean(Temp[Month==9])

plot(montly_averages, pch=10, xlab="month", ylab='degrees', main='temperature')

# or use a bar plot
barplot(monthly_averages, names.arg = c("May", "June", "July", "August", "September"),
		ylab = 'degrees F')

# set as a variable so we can display text easily
a <- barplot(monthly_averages, names.arg = c("May", "June", "July", "August", "September"),
		ylab = 'degrees F')

# pos tells you where the displayed text should go
text(a, monthly_averages, round(monthly_averages,1), pos=1)

# par sets the options to adjust plots (?par)

# you can creat multi panel grids with par
par(mfcol=c(2,2))
hist(Ozone)
hist(Temp)
hist(Solar.R)
hist(Temp, Ozone)

# and irregular grids
layout(matrix(c(1,1,2,3), 2, 2, byrow = True))
plot(Temp, Ozone)
hist(Temp)
hist(Ozone)

# or add two items with different scales
plot(Ozone)
par(new = True)
plot(Temp)

# and check out rbokeh for making nice plots, or ggplot
# http://hafen.github.io/rbokeh