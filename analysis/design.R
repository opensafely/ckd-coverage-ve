####################################################### 
### analysis parameters
####################################################### 
# Set analysis intervals and last follow-up day
postvaxcuts <- 56*0:5
postvax_periods = c("1-56", "57-112", "113-168", "169-224", "225-280")
lastfupday <- max(postvaxcuts)

# Date parameters
library(lubridate)
start_date = as_date("2020-12-01") # index date
end_date = as_date("2022-02-16")

first_pfizer = as_date("2020-12-08")
first_az = as_date("2021-01-04")
first_moderna = as_date("2021-04-13")
