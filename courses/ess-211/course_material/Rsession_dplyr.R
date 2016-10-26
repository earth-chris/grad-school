library(tidyr)
library(dplyr)
data(airquality) # load the dataset into the workspace
head(airquality) # look at its first few rows to understand its structure

## -------------- APPLYING A FUNCTION TO A SUBSET OF DATA --------------

# A loop as one solution
# Sloppy and slow, but easy to write (which can be a virtue!)
monthly.avgs   = c()                                  # initialize as an empty object
months         = unique(airquality$Month)
for (month in months) {
  month.data   = airquality[airquality$Month==month,] # split
  avgs         = apply(month.data,2,mean,na.rm=T)     # apply
  monthly.avgs = rbind(monthly.avgs, avgs)            # combine
}
monthly.avgs



# Takes a bit more time and care to write, but will run much faster
months       = unique(airquality$Month)
n_months     = length(months)
n_vars       = ncol(airquality)
monthly.avgs = matrix(nrow=n_months, ncol=n_vars)
t1=Sys.time()
for (i in 1:n_months) {
  month.data       = airquality[airquality$Month==months[i],] # split
  avgs             = apply(month.data,2,mean,na.rm=T)         # apply
  monthly.avgs[i,] = avgs                                     # combine
}
monthly.avgs



## -------------- SUMMARISE AND MUTATE FUNCTIONS --------------

# -------------- Basics --------------
# Toy dataset
air_day1 = airquality[airquality$Day==1,]
print(air_day1)



# Basic example of mutate
# Note how columns simply refered to by name 
# Note how mutate left the original columns intact and simply added the new ones
mutate(air_day1, Ozone_squared=Ozone^2, Temp_C=(Temp-32)*5/9)



# Basic example of summarise
# Note how columns simply refered to by name 
# Note how summarise only produced the newly requested columns
summarise(air_day1, Ozone_mean=mean(Ozone), Temp_sd=sd(Temp)) # can you remove the NA's?



# -------------- Using intermediate results --------------
# saturation pressure function from the first assignment
e0 = function(temp) 0.6108*exp(17.27*temp/(temp+237.3))
mutate(air_day1, Temp_C=(Temp-32)*5/9, sat_pressure=e0(Temp_C))



mutate(air_day1, sat_pressure=e0(Temp_C), Temp_C=(Temp-32)*5/9) # but order does matter!



# Traditional way
air_day1$Temp_C       = (air_day1$Temp-32)*5/9
air_day1$sat_pressure = e0(air_day1$Temp_C)



# -------------- Allowable size of columns --------------

mutate(air_day1, Temp_mean=mean(Temp)) # ok



mutate(air_day1, gt80=which(Temp>80))  # error




# -------------- Exercise --------------





tmp = mutate(air_day1, anoms=(Temp-mean(Temp)), z=anoms/sd(Temp), rmse=sqrt(mean(anoms^2)))
tmp
summarise(tmp, sig1=sum(abs(z)>1, na.rm=T))





# -------------- Transform all columns - summarise_each and mutate_each --------------

summarise(air_day1, Ozone   = mean(Ozone), 
                    Temp    = mean(Temp), 
                    Solar.R = mean(Solar.R),   
                    Wind    = mean(Wind))



# Summarise all columns
summarise_each(air_day1, "mean") # all columns



# Summarise on the Temp and Wind columns
summarise_each(air_day1, "mean", Temp, Wind) 



# Summarise all columns EXCEPT Month and Day
summarise_each(air_day1, "mean", -c(Month,Day))



# divide every column by 2
mutate_each(air_day1, funs(./2)) 



summarise_each(air_day1, funs(mean(., na.rm=T)), -c(Month,Day))



# Above line shorthand for 
summarise(air_day1, Ozone   = mean(Ozone, na.rm=T), 
                    Solar.R = mean(Solar.R, na.rm=T), 
                    Wind    = mean(Wind, na.rm=T), 
                    Temp    = mean(Temp, na.rm=T) )



## -------------- GROUPED OPERATIONS WITH GROUP_BY --------------

airquality %>% filter(Month==5) %>% head(3)
airquality %>% filter(Month==6) %>% head(3)
airquality %>% filter(Month==7) %>% head(3)



# group airquality according to unique values of its Month column
air_grouped = group_by(airquality, Month)



# summarise Ozone within each of these groups
summarise(air_grouped, Ozone_mean = mean(Ozone,na.rm=T)) # take care of NA's



## -------------- SUBSETTING ROWS AND COLUMNS WITH FILTER AND SELECT --------------

# Filter - select a subset of rows in a data frame
filter(airquality, Temp==67)           # keep rows where Temp equals 67



filter(airquality, Temp==67, Month==5) # keep rows where Temp==67 AND Month==5



# equivalent "manual" syntax
airquality[airquality$Temp==67,]
airquality[airquality$Temp==67 & airquality$Month==5,]



# Select - select a subset of columns in a data frame
select(air_day1, Ozone, Solar.R) # keep columns Ozone and Solar.R



select(air_day1, -Wind, -Temp)   # keep all columns except Wind and Temp



## -------------- CHAINING OPERATIONS TOGETHER WITH THE %>% OPERATOR --------------

# 1) Keep only rows where Month equals either 6, 7, or 8
# 2) group the result by month
# 3) take the by-month means of each columns
summarise_each(group_by(filter(airquality,Month %in% 6:8), Month), funs(mean(.,na.rm=T)))



df1 = filter(airquality, Month %in% 6:8)
df2 = group_by(df1, Month)
df3 = summarise_each(df2, funs(mean(.,na.rm=T)))
print(df3)



result = filter(airquality, Month %in% c(6,7,8)) %>% # keep months c(6,7,8), then ...
         group_by(Month) %>%                         # group by Month, then ...
         summarise_each(funs(mean(.,na.rm=T)))       # summarise each column
print(result)



## -------------- RESHAPING BETWEEN WIDE AND LONG DATAFRAMES WITH TIDYR --------------

avgs = airquality %>% group_by(Month) %>% summarise_each(funs(mean(.,na.rm=T)), -Day)
wide =  avgs %>% gather(variable, value, -Month) %>% spread(Month,value)
names(wide)[-1] = c('May','Jun','Jul','Aug','Sep')
# Example of data in wide format
wide



# Example of data in long format
wide %>% gather(Month,value, -variable) %>% spread(variable, value)



# -------------- Gather multiple columns into one --------------

avgs = airquality %>% group_by(Month) %>% summarise_each(funs(mean(.,na.rm=T)), -Day)
avgs



long = gather(avgs, variable, value, Ozone, Solar.R, Wind, Temp)
print(long)



# -------------- Spread one column into multiple columns --------------

wide = spread(long, Month, value)
print(wide)



# Rename month columns
names(wide)[-1] = c('May','Jun','Jul','Aug','Sep')


# Convert wide data back into original format
gather(wide, Month, value, -variable) %>% spread(variable,value)
