## -------------- BASIC CALCULATOR FUNCTIONALITY --------------
2+2       # Anything after a # is a "comment", and won't be executed.
3*8/2^2   # The order of opertions is exactly as you learned in algebra.
3*(8/2)^2 # Parentheses can be used to enforce a different order of operations.

## -------------- ASSIGNING VALUES TO VARIABLES --------------
x = 2+2   # evaluates 2+2, and stores ("assigns") the result to the variable "x".
x         # prints the value of x to the console
y = 2*x   # evaluates 2*x, i.e. 2*4, and stores the result in y. 
(y/x)^x

## -------------- CHARACTER (TEXT) AND LOGICAL (TRUE or FALSE) VARIABLES --------------
name = "Winston"
date = "06-11-2015"
declarationOfIndependence = "When in the course of human events ..."
x>2

## -------------- LOGICAL OPERATORS --------------
## & and # TRUE if both sides are TRUE
## | or # TRUE if either side is TRUE
## == equal to (not to be confused with =, the assignment operator!)
## != not equal to
## > greater than
## >= greater than or equal to
## < less than
## <= less than or equal to

## x = 2 # Assign the value 2 to the variable x. Does the order matter? I.e. would 2=x be ok?
## x == 3 # "Is the value of x equal to 3?" Does order matter here? Is 3==x ok?
## day = "Monday"
## day == "monday"
## x>0 & (day=="Saturday" | day=="Monday") # Is x positive and is it the weekend?
## x<3 != 3>x # Why does this give an error?
## (x<3) != (3>x) # Why does this fix it?
x = 2
!(x==2)           # x==2 is TRUE, and the negation of TRUE is FALSE
!!(5!=4 | !(x<7)) # Why? Start with the innermost expressions and work your way outward.
5!=4
!(x<7)
5!=4 | !(x<7)
!(5!=4 | !(x<7))
!!(5!=4 | !(x<7))

## -------------- IF-ELSE STATEMENTS -------------- 
## if (condition) {
##   code to executed if condition==TRUE
## }
## 
## # For multiple conditions, there's if-else
## if (condition1) {
##   code to execute if condition1==TRUE
## } else if (condition2) {
##   code to execute if condition1==FALSE & condition2==TRUE
## } else if (condition3) {
##   code to execute if condition1==FALSE & conditon2==FALSE & condition3==TRUE
## } .
##   .
##   .
## } else {
##   code to execute if all conditions are FALSE
## }
name = "Winston"
bday = "06-11"
date = "06-11"
if (date==bday) {
  print(paste("Happy Birthday,", name))
} else if (date=="12-25") {
  print(paste("Merry Christmas,", name))
} else {
  print(paste("Sorry,", name))
}

## -------------- SAVING AND RELOADING YOUR WORK -------------- 
## # display the path of the current working directory
## getwd()
## 
## # set the working directory
## setwd("path/to/desired/working/directory")
## 
## # list all variables in the workspace
## ls()
## 
## # list all files in the current working directory
## list.files()
## 
## # list all files in some other directory
## list.files("path/to/directory")
## 
## # save the variables var1, var2, etc. to an .Rdata file
## save(var1,var2, ..., file="vars_to_save.Rdata")
## 
## # save all variables in the workspace
## save.image(file="whole_workspace.Rdata")
## # Alternatively, a dialog box should pop up when you quit Rstudio,
## # asking whether you want to save your workspace.
## 
## # Load everything in this Rdata file into the workspace.
## # Note: this assumes "my_rdata_file.Rdata" is in your working directory.
## # Otherwise, you have to specify the full path.
## load("my_rdata_file.Rdata")