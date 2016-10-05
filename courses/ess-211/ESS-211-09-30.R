# an introduction to vectors

# set up string vector
cities = c("San Francisco", "London", "Marrakesh")

# set up populations
population = c(837000, 8.31e6, 929000)

# set boolean vector for whether a city has been visited
have_visited = c(TRUE, TRUE, FALSE)

# assign a name to the boolean vector to reference later
names(have_visited) <- cities

# do some weird indexing to compare two vectors of (TRUE, TRUE) to (TRUE, FALSE) and should return (TRUE, FALSE) since the second elements do not match
c(3<4, TRUE) & c(TRUE, "one"==1)

# do some weird addition to demonstrate that you should always perform calculations on same sized arrays
#  it will wrap around once it gets to the end of the shorter vector, even if it is not an even multiple of the longer vector
c(1,2,3,4) + c(1,2)
c(1,2,3,4) + c(1,2,3)

# evaluate if certain elements are within another vector
#  will return FALSE, TRUE since 1 is not in the second vector, but 3 is.
c(1,3) %in% c(5,2,3)

# seq() is like range() in python
seq(1,10) # increments by 1
seq(1,10, by=2) # increments by 2
seq(1,10, length=5) # sets up even spacing for output length 5. can return float values
1:10 # shorthand for seq(1,10)

# rep() repeats the same value a bunch of times
#  recommended using rep() for vector multiplication instead of multiplying by a scalar for readability
rep(1,10) # repeat the value 1, 10 times
rep(1:4, times=3) # repeats the vector 1:4 three times (so 1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4)
rep(1:4, each=3) # repeats each element three times (so 1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4)

# analyze population for which cities have been visited of certain pop sizes
x <- population
which(population < 1e6)

visited_cities <- which(population[have_visited] < 1e6)
class(visited_cities)

# creates a random sampling of a normal distribution with mean 1000 and sd of 500
x <- rnorm(50, mean=1000, sd=500)