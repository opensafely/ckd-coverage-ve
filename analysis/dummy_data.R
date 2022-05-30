######################################

# This script:
# - redefines CKD subgroups to approximate expected distribution
# - redefines vaccination dates according as follows:
# - (1) 5 sequential vaccination dates are assigned at intervals of 10-20 weeks
# - (2) Dates are then assigned to products according to approximate product distribution from OpenSAFELY coverage reports

# Simplifications for purposes of dummy data:
# - (1) all products used from 08/12/2020 onwards
# - (2) first and second dose are always of same product

######################################

# rcat function from https://github.com/wjchulme/dd4d/blob/main/R/rcat.R
rcat <- function(n, levels, p){
  sample(x=levels, size=n, replace=TRUE, prob=p)
}

## Import data
nsamples <- nrow(data_processed)

# Set vaccination start date
start_date <- date("2020-12-08")
end_date <- date("2022-05-11")

# Set distribution of CKD groups
data_processed$ckd_6cat <- rcat(n=nsamples, c("No CKD", "CKD3a", "CKD3b", "CKD4", "CKD5",
                                              "RRT (dialysis)", "RRT (Tx)"), c(0.10,0.50,0.25,0.06,0.03,0.03,0.03))
# Recalculate CKD5 categories and RRT mismatch based on random assignments above
data_processed <- data_processed %>%
  mutate(
    ckd_5cat = ckd_6cat,
    ckd_5cat = ifelse(ckd_6cat == "CKD4" | ckd_6cat == "CKD5", "CKD4-5", ckd_5cat),
    rrt_mismatch = ifelse((ckd_5cat=="CKD3a" | ckd_5cat=="CKD3b" | ckd_5cat=="CKD4-5") & (dialysis==1 | kidney_transplant==1), 1, 0)
  )

# Increase proportion of cev_other to 10% to ensure effective running of later models
data_processed$cev_other <- as.numeric(rcat(n=nsamples, c("0", "1"), c(0.9,0.1)))

# Set cev_other to 0 if in dialysis or Tx group
data_processed$cev_other[data_processed$ckd_6cat %in% c("RRT (dialysis)", "RRT (Tx)")] = 0

# Set vaccine coverage for primary doses, third, and fourth dose
primary_coverage <- 0.95
third_coverage <- 0.90 # assignments only retained if primary doses given, so final prevalence ~0.95*0.90=0.855
fourth_coverage <- 0.40 # assignments only retained if primary doses given, so final prevalence ~0.855*0.40=0.342
fifth_coverage <- 0.05 # assignments only retained if primary doses given, so final prevalence ~0.342*0.05=0.0171

# Set distribution of vaccine type for primary doses, third, and fourth doses (https://reports.opensafely.org/reports/vaccine-coverage/)
primary_vax_type <- rcat(n=nsamples, c("pfizer","az","moderna",""), c(0.47*primary_coverage,0.5*primary_coverage,0.03*primary_coverage,1-primary_coverage))
third_vax_type <- rcat(n=nsamples, c("pfizer","az","moderna",""), c(0.76*third_coverage,0.01*third_coverage,0.23*third_coverage,1-third_coverage))
fourth_vax_type <- rcat(n=nsamples, c("pfizer","az","moderna",""), c(0.76*fourth_coverage,0.01*fourth_coverage,0.23*fourth_coverage,1-fourth_coverage))
fifth_vax_type <- rcat(n=nsamples, c("pfizer","az","moderna",""), c(0.76*fifth_coverage,0.01*fifth_coverage,0.23*fifth_coverage,1-fifth_coverage))

# Set distribution of dates for dose 1 of each vaccine - all administered within 10 weeks of product initiation for purposes of dummy data
data_processed$covid_vax_date_1 <- as_date(runif(nsamples, start_date, start_date+70)) 

# Set distribution of dates for doses 2-5 of each vaccine - all administered 10-20 weeks after preceding dose for purposes of dummy data
for (i in 1:nsamples) { 
  data_processed$covid_vax_date_2[i] <- as_date(sample((data_processed$covid_vax_date_1[i]+70):(data_processed$covid_vax_date_1[i]+140),1))
  data_processed$covid_vax_date_3[i] <- as_date(sample((data_processed$covid_vax_date_2[i]+70):(data_processed$covid_vax_date_2[i]+140),1))
  data_processed$covid_vax_date_4[i] <- as_date(sample((data_processed$covid_vax_date_3[i]+70):(data_processed$covid_vax_date_3[i]+140),1))
  data_processed$covid_vax_date_5[i] <- as_date(sample((data_processed$covid_vax_date_4[i]+70):(data_processed$covid_vax_date_4[i]+140),1))
}

# Middle dose date of theoretical distribution:
# Dose 1: 2021-01-12 #as_date("2020-12-08") + 35
# Dose 2: 2021-04-27 #as_date("2021-01-12") + 105
# Dose 3: 2021-08-10 #as_date("2021-04-27") + 105
# Dose 4: 2021-11-23 #as_date("2021-08-10") + 105
# Dose 5: 2022-03-08 #as_date("2021-11-23") + 105

# Set any dates above end date to NA
data_processed$covid_vax_date_1[data_processed$covid_vax_date_1>end_date] = NA
data_processed$covid_vax_date_2[data_processed$covid_vax_date_2>end_date] = NA
data_processed$covid_vax_date_3[data_processed$covid_vax_date_3>end_date] = NA
data_processed$covid_vax_date_4[data_processed$covid_vax_date_4>end_date] = NA
data_processed$covid_vax_date_5[data_processed$covid_vax_date_5>end_date] = NA

# Remove vaccination dates for proportion with unassigned products
data_processed$covid_vax_date_1[primary_vax_type==""] = NA
data_processed$covid_vax_date_2[primary_vax_type==""] = NA
data_processed$covid_vax_date_3[third_vax_type=="" | primary_vax_type==""] = NA
data_processed$covid_vax_date_4[fourth_vax_type=="" | third_vax_type=="" | primary_vax_type==""] = NA
data_processed$covid_vax_date_5[fifth_vax_type=="" | fourth_vax_type=="" | third_vax_type=="" | primary_vax_type==""] = NA

# Reset all product-specific dates to NA, then select global vaccine date for corresponding doses according to product distribution defined above
# pfizer
data_processed$pfizer_date_1 = data_processed$pfizer_date_2 = data_processed$pfizer_date_3 = data_processed$pfizer_date_4 = data_processed$pfizer_date_5 = as_date(NA)
data_processed$pfizer_date_1[primary_vax_type=="pfizer"] = data_processed$covid_vax_date_1[primary_vax_type=="pfizer"]
data_processed$pfizer_date_2[primary_vax_type=="pfizer"] = data_processed$covid_vax_date_2[primary_vax_type=="pfizer"]
data_processed$pfizer_date_3[third_vax_type=="pfizer"] = data_processed$covid_vax_date_3[third_vax_type=="pfizer"]
data_processed$pfizer_date_4[fourth_vax_type=="pfizer"] = data_processed$covid_vax_date_4[fourth_vax_type=="pfizer"]
data_processed$pfizer_date_5[fifth_vax_type=="pfizer"] = data_processed$covid_vax_date_5[fifth_vax_type=="pfizer"]

# az
data_processed$az_date_1 = data_processed$az_date_2 = data_processed$az_date_3 = data_processed$az_date_4 = data_processed$az_date_5 = as_date(NA)
data_processed$az_date_1[primary_vax_type=="az"] = data_processed$covid_vax_date_1[primary_vax_type=="az"]
data_processed$az_date_2[primary_vax_type=="az"] = data_processed$covid_vax_date_2[primary_vax_type=="az"]
data_processed$az_date_3[third_vax_type=="az"] = data_processed$covid_vax_date_3[third_vax_type=="az"]
data_processed$az_date_4[fourth_vax_type=="az"] = data_processed$covid_vax_date_4[fourth_vax_type=="az"]
data_processed$az_date_5[fifth_vax_type=="az"] = data_processed$covid_vax_date_5[fifth_vax_type=="az"]

# moderna
data_processed$moderna_date_1 = data_processed$moderna_date_2 = data_processed$moderna_date_3 = data_processed$moderna_date_4 = data_processed$moderna_date_5 = as_date(NA)
data_processed$moderna_date_1[primary_vax_type=="moderna"] = data_processed$covid_vax_date_1[primary_vax_type=="moderna"]
data_processed$moderna_date_2[primary_vax_type=="moderna"] = data_processed$covid_vax_date_2[primary_vax_type=="moderna"]
data_processed$moderna_date_3[third_vax_type=="moderna"] = data_processed$covid_vax_date_3[third_vax_type=="moderna"]
data_processed$moderna_date_4[fourth_vax_type=="moderna"] = data_processed$covid_vax_date_4[fourth_vax_type=="moderna"]
data_processed$moderna_date_5[fifth_vax_type=="moderna"] = data_processed$covid_vax_date_5[fifth_vax_type=="moderna"]
