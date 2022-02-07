######################################

# This script:
# - redefines vaccination dates according as follows:
# - (1) 4 sequential vaccination dates are assigned at intervals of 12-20 weeks
# - (2) dates are then assigned to products according to approximate product distribution from OpenSAFELY coverage reports

# simplifications for purposes of dummy data:
# - (1) all products used from 08/12/2020 onwards
# - (2) first and second dose are always of same product

######################################

# rcat function from https://github.com/wjchulme/dd4d/blob/main/R/rcat.R
rcat <- function(n, levels, p){
  sample(x=levels, size=n, replace=TRUE, prob=p)
}

## Import data
#data_processed <- read_rds(here::here("output", "data", "data_processed.rds"))
nsamples <- nrow(data_processed)

# set product-specific start dates
start_date <- date("2020-12-08")
end_date <- date("2021-12-31")

# set vaccine coverage for primary doses, third, and fourth dose
primary_coverage <- 0.95
third_coverage <- 0.6 # assignments only retained if primary doses given, so final prevalence ~0.6*0.95=0.57
fourth_coverage <- 0.2 # assignments only retained if primary doses given, so final prevalence ~0.57*0.2=0.11

# set distribution of vaccine type for primary doses, third, and fourth doses (https://reports.opensafely.org/reports/vaccine-coverage/)
primary_vax_type <- rcat(n=nsamples, c("pfizer","az","moderna",""), c(0.47*primary_coverage,0.5*primary_coverage,0.03*primary_coverage,1-primary_coverage))
third_vax_type <- rcat(n=nsamples, c("pfizer","az","moderna",""), c(0.76*third_coverage,0.01*third_coverage,0.23*third_coverage,1-third_coverage))
fourth_vax_type <- rcat(n=nsamples, c("pfizer","az","moderna",""), c(0.76*fourth_coverage,0.01*fourth_coverage,0.23*fourth_coverage,1-fourth_coverage))

# set distribution of dates for dose 1 of each vaccine - all administered within 12 weeks of product initiation for purposes of dummy data
data_processed$covid_vax_date_1 <- as_date(runif(nsamples, start_date, start_date+84)) 

# set distribution of dates for doses 2-4 of each vaccine - all administered 12-20 weeks after preceding dose for purposes of dummy data
for (i in 1:nsamples) { 
  data_processed$covid_vax_date_2[i] <- as_date(sample((data_processed$covid_vax_date_1[i]+84):(data_processed$covid_vax_date_1[i]+140),1))
  data_processed$covid_vax_date_3[i] <- as_date(sample((data_processed$covid_vax_date_2[i]+84):(data_processed$covid_vax_date_2[i]+140),1))
  data_processed$covid_vax_date_4[i] <- as_date(sample((data_processed$covid_vax_date_3[i]+84):(data_processed$covid_vax_date_3[i]+140),1))
}

# maximum theoretical dates by dose are
# dose 1: 2021-03-02 (2020-12-08 + 84)
# dose 2: 2021-07-20 (2021-03-02 + 140)
# dose 3: 2021-12-07 (2021-07-20  + 140)
# dose 4: 2022-04-26 (2021-12-07  + 140)

# set any dates above end date to NA
data_processed$covid_vax_date_4[data_processed$covid_vax_date_4>end_date] = NA

# remove vaccination dates for proportion with unassigned products
data_processed$covid_vax_date_1[primary_vax_type==""] = NA
data_processed$covid_vax_date_2[primary_vax_type==""] = NA
data_processed$covid_vax_date_3[third_vax_type=="" | primary_vax_type==""] = NA
data_processed$covid_vax_date_4[fourth_vax_type=="" | third_vax_type=="" | primary_vax_type==""] = NA

# reset all product-specific dates to NA, then select global vaccine date for corresponding doses according to product distribution defined above
# pfizer
data_processed$pfizer_date_1 = data_processed$pfizer_date_2 = data_processed$pfizer_date_3 = data_processed$pfizer_date_4 = as_date(NA)
data_processed$pfizer_date_1[primary_vax_type=="pfizer"] = data_processed$covid_vax_date_1[primary_vax_type=="pfizer"]
data_processed$pfizer_date_2[primary_vax_type=="pfizer"] = data_processed$covid_vax_date_2[primary_vax_type=="pfizer"]
data_processed$pfizer_date_3[third_vax_type=="pfizer"] = data_processed$covid_vax_date_3[third_vax_type=="pfizer"]
data_processed$pfizer_date_4[fourth_vax_type=="pfizer"] = data_processed$covid_vax_date_4[fourth_vax_type=="pfizer"]

# az
data_processed$az_date_1 = data_processed$az_date_2 = data_processed$az_date_3 = data_processed$az_date_4 = as_date(NA)
data_processed$az_date_1[primary_vax_type=="az"] = data_processed$covid_vax_date_1[primary_vax_type=="az"]
data_processed$az_date_2[primary_vax_type=="az"] = data_processed$covid_vax_date_2[primary_vax_type=="az"]
data_processed$az_date_3[third_vax_type=="az"] = data_processed$covid_vax_date_3[third_vax_type=="az"]
data_processed$az_date_4[fourth_vax_type=="az"] = data_processed$covid_vax_date_4[fourth_vax_type=="az"]

# moderna
data_processed$moderna_date_1 = data_processed$moderna_date_2 = data_processed$moderna_date_3 = data_processed$moderna_date_4 = as_date(NA)
data_processed$moderna_date_1[primary_vax_type=="moderna"] = data_processed$covid_vax_date_1[primary_vax_type=="moderna"]
data_processed$moderna_date_2[primary_vax_type=="moderna"] = data_processed$covid_vax_date_2[primary_vax_type=="moderna"]
data_processed$moderna_date_3[third_vax_type=="moderna"] = data_processed$covid_vax_date_3[third_vax_type=="moderna"]
data_processed$moderna_date_4[fourth_vax_type=="moderna"] = data_processed$covid_vax_date_4[fourth_vax_type=="moderna"]




