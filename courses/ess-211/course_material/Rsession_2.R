knitr::opts_chunk$set(comment=NA, error=T, warning=T) # for html
# knitr::opts_chunk$set(comment=NA, error=T, warning=T, tidy=T) # tidy=T sometimes helpful for pdf
cities = c("San Francisco","London","Marrakesh") # a character vector
population = c(8.37e5,8.31e6,9.29e5) # a numeric vector (using scientific notation)
have_visited = c(TRUE,TRUE,FALSE) # a logical vector
cities[2] # 2nd element of the "cities" vector
population[c(1,3)] # 1st and 3rd elements of the "population" vector
have_visited # before
have_visited[3] = TRUE # read: "assign TRUE to the third element of have_visited"
have_visited # after
have_visited
# The elements of have_visited that are TRUE are c(1,2), so indexing cities
# with have_visited is the same as indexing cities with c(1,2).
cities[c(1,2)]
cities[have_visited]
names(have_visited) # no names yet
names(have_visited) = cities
names(have_visited) # now the values of "cities" are indices of "have_visited"
have_visited["London"] # and we can use them just like integer indices
c(1,2,3)*2  # multiply each element of c(1,2,3) by 2
c(1,2,3)>=2 # evaluate whether each element of c(1,2,3) is greater than or equal to 2
c(2,3,4) + c(10,20,30) # means c(2+3, 3+3, 4+3)
c(2,3,4) < c(10,20,30) # means c(2<3, 3<3, 4<3)
c(3<4, TRUE) & c(TRUE, "one"==1) # means c(TRUE & TRUE, TRUE & FALSE)
# These ...
c(1,2,3,4) + c(1,2)
c(1,2,3,4) + c(1,2,3)

# are equivalent to these ...
c(1,2,3,4) + c(1,2,1,2) # c(1,2) gets recycled twice to have length 4
c(1,2,3,4) + c(1,2,3,1) # c(1,2,3) starts to be recycled, but stops at length 4
c(1,3) %in% c(5,3,2) 
# 1 is NOT in the vector c(5,3,2)
# 3 IS in the vector c(5,3,2)
1==5 | 1==3 | 1==2
3==5 | 3==2 | 3==2
seq(1,10)  # sequence from 1 to 10, in steps of 1 (the default)
seq(1,10, by=2)  # sequence from 1 to 10, in steps of 2
seq(1,10, length=5) # sequence of 5 evenly-spaced values between 1 and 10
1:10 # shorthand for seq(1,10); assumes integer sequences with steps of 1
rep(1,10)  # repeat the value 1, 10 times
rep(1:4, times=3)  # repeat the whole vector 1:4 three times
rep(1:4, each=3) # repeat each successive element of 1:4 three times
1:4-1 # First forms the vector 1:4, then subtracts 1 from every element
1:(4-1) # Maybe this is what you wanted instead?
x = population
which(population < 1e6) # which elements are TRUE (applies only to logical vectors)
length(x) # the length of the vector (i.e. how many elements it has)
class(x) # what variable type (e.g. "numeric"")
str(x) # prints helpful information, including the type, length, and first few values
any(have_visited) # are any elements TRUE?
all(have_visited) # are all elements TRUE?
sum(have_visited) # how many elements are TRUE? Know why this works!
