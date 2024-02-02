# these are equivalent:
# note that a "#" starts a comment line that will not be run as code
x = 1
x <- 1

# "x" can be re-assigned to a different variable
x <- 2
x

# here, we are reading in a data file called "bird_morphometrics.csv" 
# that is available to download from Canvas
# first, we need to set the "working directory", which is the folder on your computer
# where R will access files. Remember when we talked about folder structure in the 
# first week of class

setwd("/Users/ojohnson/Documents/Bioinformatics/Data")
getwd()
bm <- read.csv("bird_morphometrics.csv")

# once we read in the file and assig it to the variable "bm", we can look at the data
# this looks at the first few lines of the data file
head(bm) 

# next we can plot it using the 'plot' function
# we use the dollar sign to access column names in the data
# data of this type is called a 'data frame' (columns and rows of data)
plot(bm$Wing_Chord, bm$Tail_Length)

# we could plot this with log-transformed data
plot(log(bm$Wing_Chord), 
     log(bm$Tail_Length))

# R has some datasets that are already available within the language
# basically for users to play around with and test functions
# the most widely use one is 'mtcars', which is information from different makes
# and models of cars. Instead of reading this in from a file, we use the 'attach' function
attach(mtcars)
head(mtcars)

# The 'View' function opens the entire data frame in a new tab
# The 't' function "transforms" the data frame (swaps the rows and the columns)
View(t(mtcars))



