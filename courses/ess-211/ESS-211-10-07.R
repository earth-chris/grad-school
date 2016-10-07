# calculate root-mean square
rms <- function(x){
	return(sqrt(mean(x^2)))
}

# can use variables in defining functions

new_rms <- function(vec, func){
	return(sqrt(func(vec^2)))
}

# run rms and rss
root_sum_squares <- new_rms(rep(3,4), sum)
root_mean_squares <- new_rms(rep(3,4), mean)

# get to work on learning apply for functions

cities = matrix(sample(seq(80, 110), 15), nrow=3, ncol=5)
rownames(cities) = c("San Francisco", "London", "Marrakesh")
colnames(cities) + seq(2011, 2015)

cities

# get the mean for each city (using dimension 1, or by row)
apply(cities, 1, mean)

# get the max for each year (using dimension 2, or by column)
apply(cities, 2, max)

# you can deal with NAs. na.rm is an argument to mean(), not necessarily to apply() 
apply(cities, 1, mean, na.rm = TRUE)

# 
result = rep(NA, nrow(cities))
for (i in 1:nrow(cities)){
	x = cities[i,]
	result[i] = sqrt(mean(x^2))
}
result

# talk about lapply
# applies a function to every element of a 'list'
a_list = list(c(19,2738, 11395), c(2375, 1023), c(1234, 1243))

# returns a list of where the maximum values are
lapply(a_list, which.max)
a_list[lapply(a_list, which.max)]

# sapply returns a vecor, instead
maxs.vector = sapply(a_list, which.max)

a_list[maxs.vector]