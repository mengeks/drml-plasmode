#---------------------------------------#
#                                       #
#   Data Wrangling in R                 #
#                                       #
#---------------------------------------#

#Data Science Workflow

# 1. Import
# 2. Clean (Tidy)
# 3. Understand:
#     - Transform
#     - Visualize
#     - Model
# 4. Communicate



# The tidyverse is an opinionated collection of R packages designed for data science. 
# All packages share an underlying design philosophy, grammar, and data structures. 
# install.packages("tidyverse")
# library(tidyverse)

# *** Core tidyverse packages ***
# readr - import data
# tidyr - for data "tidying"
# dplyr - data manipulation
# stringr - string manipulation
# forcats -factors handling
# ggplot2 - data visualization
# tibble - modern representation of a data frame
# purrr - functional programming
#
# *** Other related packages ***
# hms - time manipulation
# lubridate - dates/times manimpulation
# DBI - databases
# haven - data imported from SPSS, SAS and Stata 
# readxl - excel data import
# rvest - web scaping

#----------------------------------------------






# If the packages have not been installed, they can be installed using the following commands:
#install.packages("readxl")
#install.packages("dplyr")
#install.packages("stringr")

# Load packages into R environment
library(readxl)
library(dplyr)

# Reading data saved in Excel format often some challenges
visit <- read_excel("MedData.xlsx", sheet = "Visits")
visit

# We can see that the file contains some extra records at the top and at the bottom of the sheet.
# Use options to read_excel() function to skip these extra lines:
visit <- read_excel("MedData.xlsx", 
                    sheet = "Visits", 
                    skip = 3,  # how many lines to skip at the top of the sheet
                    n_max = 13)  # how many observations to read
visit


# Let's take a look at the type of each record to make sure that each column was read correctly
# This can be done with str() function.
# However dplyr package has also handy glimpse() function that has a nicely formatted output
str(visit)
glimpse(visit)

# Notice that DBP column is marked as "character" column even though it contains numeric values.
# This normally means that there are some missing values that were coded in non-standard way, so 
# R did not recognize them as missing values.
# In this case you should examine the column and try to identify the missing values:
min(visit$DBP)
#[1] "##N/A"

# Now we can go back and read file again adding "##N/A" as a Missing Data code:
visit <- read_excel("MedData.xlsx", 
                    sheet = "Visits", 
                    skip = 3,
                    n_max = 13,
                    na=c("","NA","##N/A"))
glimpse(visit)

# "pipe" symbol in R was introduced in magrittr package and then it was also implemented in dplyr package we will use today
# It allows to send the output of one function as an input objec to the next function:
visit %>% glimpse()
visit %>% summary()

# The following chain of functions:
sort(unique(visit$`Patient ID`), decreasing=TRUE)
# can be rewritten as
visit$`Patient ID` %>% 
  unique()  %>%            # find unique values
  sort(decreasing=TRUE)    # sort  in descending order

#-------------------------#
#    exercise
#-------------------------
# Read and explore "Patient Data" sheet from the same excel file
# pinfo <- read_excel(-----)
pinfo <- read_excel("MedData.xlsx", 
                    sheet = "Patient Data")

#  Use glimpse, head, summary, and other functions to explore the dataset
pinfo %>% glimpse()
pinfo %>% head()
pinfo %>% summary()



#-------------------------------------------------------------------------










#---------------------------------------------------------------
## A dataset is a collection of values, usually either 
# - numbers (if quantitative),  
# - strings (if qualitative),
# - logical (if binary).
# Every value belongs to a variable and an observation. 
# A variable contains all values that measure the same underlying attribute (like phone, age) across units. 
# An observation contains all values measured on the same unit (like a person, or a day, or an event) across attributes

## Clean or Tidy Data:
# - Each variable forms a column
# - Each observation forms a row
# Each type of observational unit forms a table

## Examples of messy datasets
# - Column headers are values, not variable names
# - Multiple variables are stored in one column
# - Variables stored in both rows and columns
# - Multiple types of observational units are stored in the same table
# - A single observational unit is stored in multiple tables.


#---------------------------------------------------------------------------

#
# Using dplyr package to 

# - Rename variables ( rename() )
# - Filter observations by their values ( filter() )
# - Reorder rows ( arrange() )
# - Select specific columns/variables ( select() )
# - Create new variables with functions of existing variables ( mutate() )
# - Summarise (summarise() )


# We can select columns in a Data Frame using "$" symbol, i.e:
visit$DBP

# It is not as easy to do so when the column name contains spaces or other special characters:
visit$`Patient ID`

# In some cases in makes sense to rename some columns to make it easy to work 
visit.clean <- visit %>%
  rename(id="Patient ID", 
         admission="Admission date",
         discharge="Discharge date",
         pulse="Heart Rate")

visit.clean %>% head()






#---------------------------------------------
# Working with strings
#---------------------------------------------
visit.clean$Allergies
# Converting strings to upper and lower case
toupper(visit.clean$Allergies)
# or using pipe operator
visit.clean$Allergies %>% toupper()

#Search character vector for a specific pattern:
grepl ("pain", visit.clean$Symptoms,  ignore.case=TRUE)

# Find all strings that start with "chest"
grepl ("^chest", visit.clean$Symptoms,  ignore.case=TRUE)

# Find all strings that end with "pain"
grepl ("pain$", visit.clean$Symptoms,  ignore.case=TRUE)

# Find all strings that contain either fever or pain (or both)
grepl ("fever|pain", visit.clean$Symptoms,  ignore.case=TRUE)


# grep( value=FALSE): returns a vector of indices of element where pattern is found
# grep( value=TRUE): returns a vector of element where pattern is found
grep ("fever", visit.clean$Symptoms,  ignore.case=TRUE, value=FALSE)
grep ("fever", visit.clean$Symptoms,  ignore.case=TRUE, value=TRUE)

# sub() function can be used to substitute the first occurance of a pattern with another string
visit.clean$Symptoms %>% head()
sub (",", ";", visit.clean$Symptoms)
# gsub() function can be used to substitute all occurances of a pattern with another string
gsub (",", ";", visit.clean$Symptoms)








#------------------------------------------
# Working with dates
#------------------------------------------
class(visit.clean$admission)

curr.time <- Sys.time()  # get current date and time
curr.date <- Sys.Date()  # get current date 
str(curr.time)           # view the structure of an object
class(curr.time)         # view type of an object

# some systems do not have timezone set up
Sys.timezone()

# convert character string to POSIXlt
class("2019-01-29 11:30:00")
t1=as.POSIXct("2019-01-29 11:30:00", "%Y-%m-%d %H:%M:%S", tz="EST")
OlsonNames()  # list of Time Zones
str(t1) 
class(t1)

# install.packages("lubridate")
library(lubridate)
mydates <- c("2/25/2018", "3/5/2017", "4-18-2018", "7.5.2017")
class(mydates)
newdates <- mdy(mydates)
class(newdates)
newdates

# There are a number of handy packages that have various functions to work with dates:
# lubridate
# chron

# Some more examples working with dates and times can be found:
# http://rcs.bu.edu/examples/r/timesSeries/dateTime.R

# Here we will calculate the length of stay of each patient in the hospital
visit.clean$discharge - visit.clean$admission
#The result variable has "difftime" class:
stay <- visit.clean$discharge - visit.clean$admission
class(stay)
#To convert it to a numeric value use as.numeric() function:
as.numeric(visit.clean$discharge - visit.clean$admission)


# ********************************************************************************








#----------------------
#Filtering rows
#----------------------
visit.clean %>% filter( !is.na(DBP) )

## Filter observations
# R logical operators
# ? Comparoson
# >  >=
# ==  !=
# <   <=
# is.na   !is.na
#  %in%


#************
# Exericise:
#************
# Select only those records for which pulse columns have values 100 and greater
# visit.clean %>% filter( --- )
visit.clean %>% filter( pulse >= 100)
visit.clean$pulse
# Select only those records for which DBP is less than 60 or SBP is greater than 120
# visit.clean %>% filter( --- | ---)
visit.clean %>% filter( DBP < 60 | SBP > 100)

# Select only those records for which Temperature is greater than 99 and Symptoms include "fever"
# visit.clean %>% filter( Temperature --- &  grepl("---", Symptoms, ignore.case=T)   )
feverfilter <- visit.clean %>% filter(Temperature>99 & grepl("fever", Symptoms, ignore.case = T))
feverfilter$Symptoms
feverfilter$Temperature

#***********









#----------------------
# Selecting Columns
#----------------------
visit.clean %>% select(id, Temperature:pulse)

# dplyr comes with a set of helper functions that can help you select groups of variables inside a select() call:
#    starts_with("XY"): every name that starts with "XY",
#    ends_with("XY"): every name that ends with "XY",
#    contains("XY"): every name that contains "XY",
#    matches("XY"): every name that matches "XY", where "XY" can be a regular expression,
#    num_range("x", 1:15): the variables named x01, x02, x03,..., x15,
#    one_of(XY): every name that appears in x, where XY is a character vector.

visit.clean %>% select( id, ends_with("BP") )








#---------------------------------------------------
# Modifying existing and/or creating  new variables
#---------------------------------------------------
visit.clean %>% mutate( Temperature = (Temperature-32)*5/9, stay = as.numeric(discharge - admission) )

#************
# Exericise:
#************

# Create a new column MAP which is equal to SBP/3 + 2*DBP/3
#visit.clean %>% mutate( --- ) 
visit.w.map <- visit.clean %>% mutate(MAP = SBP/3 + 2*DBP/3)
visit.w.map$MAP

# Let's put it all together:
# Use visit.clean dataframe as input and
# - select only those columns where Temperature is greater than 99F
# - select columns ID, DBP and SBP
# - calculate new variable MAP
#res <- visit.clean %>%
#            filter( --- ) %>%
#            select( --- ) %>%
#            mutate( --- ) 
# res

visit.clean %>% 
  filter(Temperature > 99)  %>% 
  select(id, contains("BP")) %>%
  mutate(MAP = SBP/3 + 2*DBP/3)







#---------------------------------------------------
#  Calculating summaries
#---------------------------------------------------
visit.clean %>% summarise( NumVisits = n(), max.t = max(Temperature), min.T = min(Temperature) )


# built-in functions often used within summarise()
# averages: mean(), median()
# spread: sd(), IQR(), mad()
# range: min(), max()
# count: n(), n_distinct()

#************
# Exericise:
#************

# Calculate the number of distinct patients
#visit.clean %>% summarise( N = --- ) 
visit.clean %>% summarise( N = n_distinct(id)) 








#---------------------------------------------------
#  Group by one or more variables
#---------------------------------------------------
visit.clean %>% group_by(id) %>% summarise( ave.pulse = mean(pulse) )


# Useful functions often used with group_by()
# first(), last(), nth()

# For each patient select the first record and find the lenght of stay
visit.clean %>% group_by(id) %>% summarise( first(discharge - admission) )
visit.clean %>% group_by(id) %>% summarise( sum(discharge - admission) )








#---------------------------------------------------
#  Sorting dataframe by one or more variables
#---------------------------------------------------
visit.clean %>% arrange(id, admission)

# sort in descending order
visit.clean %>% arrange(id, desc(admission))





#---------------------------------------------------
#  Joining 2 dataframes
#---------------------------------------------------


# First let's read the data for each patient:
pinfo <- read_excel("MedData.xlsx", sheet = "Patient Data")

#let's make sure the date is well formatted
pinfo %>%head()
pinfo %>%glimpse()
pinfo %>%summary()

# Now we want to have a single dataframe that contains patient infomration and patient visit informtaion
# There are a number of join* functions in dplyr package:
# inner_join
# left_join
# right_join
# full_join
# semi_join  - return all rows from x where there are matching values in y, keeping just columns from x.
# anti_join  - return all rows from x where there are not matching values in y, keeping just columns from x.

#Let's try full join of both dataframes we have
full_join(visit.clean, pinfo, by = c("id"="ID"))

# left join
result <- visit.clean  %>% left_join(pinfo, by = c("id"="ID"))
result

#******************
# Exercise
#******************

# Calculate the mean stay of the hospital for Male and Female patients:
# result %>%
#   filter( --- ) %>%            # select observations where Gender is not missing
#   mutate( stay = --- ) %>%     # create a new variable stay equal to the length of stay in the hospital
#   group_by ( --- ) %>%         # group by gender
#   summarize( ave.stay = --- )  # calculate mean length of stay
  
result %>%
  filter( !is.na(Gender) ) %>%            # select observations where Gender is not missing
  mutate( stay = as.numeric(discharge - admission)) %>%     # create a new variable stay equal to the length of stay in the hospital
  group_by ( Gender ) %>%         # group by gender
  summarize( ave.stay = mean(stay) )  # calculate mean length of stay








#************************************************************
## Converting Wide dataframes to long dataframes and back
#*************************************************************
library(tidyr)
city.temps <- data.frame( time = as.Date('2018-09-03') + 0:4,
                          Boston = rnorm( 5, 75, 3),
                          Denver = rnorm( 5, 85, 5),
                          SanFrancisco = rnorm(5, 80, 5),
                          Austin = rnorm( 5, 90, 5) )
city.temps

# Use gather to rearrange a table into a tidy data.frame:
# gather(data, key = "key", value = "value", ..., na.rm = FALSE,
# convert = FALSE, factor_key = FALSE)
city.temps2 <- gather( city.temps, 
                       key = "City",
                       value = "Temperature",
                       -time,  # collection of the columns,
                       factor_key = TRUE) # if key variable needs to be converted to a factor
city.temps2
# or
city.temps2 <- gather( city.temps, 
                       key = "City",
                       value = "Temperature",
                       Boston: Austin,  # collection of the columns,
                       factor_key = TRUE) # if key variable needs to be converted to a factor
#or
city.temps2 <- gather( city.temps, 
                       key = "City",
                       value = "Temperature",
                       c("Boston","Denver","SanFrancisco","Austin"),  # collection of the columns,
                       factor_key = TRUE) # if key variable needs to be converted to a factor


city.temps2
glimpse(city.temps2)


#  Sometimes it is useful to be able to
#perform the opposite operation: convert long format dataframe into wide representation
city.temps3 <- spread( city.temps2, City, Temperature)
city.temps3




