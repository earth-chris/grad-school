# session on for and while loops

# looping through a string vector
for (word in c("hello", "world")){
	print(word)
}

# loop through a sequence
some_vector = seq(1,10, by=3)
for (i in 1:length(some_vector)){
	print(paste("the value of i is: ", i))
	print(paste("the i'th element of the vector is: ", some_vector[i]))
}

# loop through a fibonacci sequence
#  create an na array first
fib_seq = rep(NA, 10)
fib_seq[c(1,2)] = 1 # add 1 to elements 1 and 2
for (i in 3:length(fib_seq)){
	fib_seq[i] = fib_seq[i-2] + fib_seq[i-1]
}

# while loops
#  be careful to not run an infinite while loop!
#  the below code will run, because it meets the condition, if you do x-1, it runs FOREVER (until int rollever :o)
x = 6
while(x < 10){
	print(x)
	x = x+1
}

# data types within a vector must be consistent
#  in this case, it converts boolean TRUE to 1
#  will always cast to the highest data type
c(10,TRUE)
c(10,TRUE, "hello")

# can explicitly cast data types using as.character, as.integer, etc.
as.character(c(10,TRUE))

# matrices 
matrix(1:6, nrow=2, ncol=3)

# r fills matrices in column-major order
matrix(1:6, nrow=3, ncol=2)
matrix(1:6, nrow=2, ncol=6)

# can do row-major stuff, too
mat1 = matrix(1:6, nrow=2)
mat2 = matrix(1:6, nrow=2, byrow=TRUE)
mat = mat1

# can transpose matrices using t() function
tmat = t(mat1)

# indexing
mat[2,3]
mat[2,c(1,3)]
mat[,c(2,3)]
mat[1]
mat[-1]

# data frames
# data frames can be considered a series of vectors
cities = c("San Francisco", "London", "Marrake")
have_visited = c(TRUE, TRUE, FALSE)
population = c(837000, 8310000, 929000)
df = data.frame(cities, population, have_visited)

df[,"cities"]
df$cities

# find cities that have been visited
df$cities[df$have_visited]

# but that shows up as a factor. can explicitlycast as a string 
as.character(df$cities[df$have_visited])

# default naming for columns is based on the variable name it was assigned
#  in this case the names are 'cities', 'population', and 'have_visited'
names(df)

# or explicitly name them
names(df) = c('city', 'pop', 'visited')
names(df)

# can use bind() and cbind() to add rows or columns to a data frame


# lists
my_list = list(a_matrix=matrix(1:4,nrow=2),
			char_vec = c("hello", "World"),
			df = data.fram(cbind(matrix(1:4, nrow=2), c("Hello", "world"))))
	
#