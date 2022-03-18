####################################################### 
### analysis parameters
####################################################### 
# Set analysis intervals and last follow-up day
postvaxcuts <- 56*0:5
postvax_periods = c("1-56", "57-112", "113-168", "169-224", "225-280")
lastfupday <- max(postvaxcuts)