# # # # # # # # # # # # # # # # # # # # #
# This script creates a table of incidence rates for diffferent outcomes, stratified by vaccine type and times since vaccination.
# # # # # # # # # # # # # # # # # # # # #

# Preliminaries ----

## Import libraries ----
library('tidyverse')
library('here')
library('glue')
library('survival')
library('gt')
library('gtsummary')
library('scales')
library('lubridate')

## Set whether or not to import matched data - ONLY UPDATE THIS LINE BETWEEN MATCHED/UNMATCHED OUTPUTS
matched=TRUE

## Define input and output file names
if (matched) {
  input_name = "data_cohort_VE_matched.rds"
  output_csv = "table_irr_matched_redacted.csv"
  output_rds = "table_irr_matched_redacted.rds"
} else {
  input_name = "data_cohort_VE.rds"
  output_csv = "table_irr_redacted.csv"
  output_rds = "table_irr_redacted.rds"
}

## Import custom user functions
source(here::here("analysis", "functions.R"))

## Set analysis intervals
postvaxcuts <- 56*0:5
lastfupday <- max(postvaxcuts)

## Import processed data ----
data_cohort <- read_rds(here::here("output", "data", input_name))

# Split into vaccine-specific dataframes
data_cohort_AZ <- subset(data_cohort, vax2_type=="az")
data_cohort_BNT <- subset(data_cohort, vax2_type=="pfizer")
#vec = data_cohort_AZ$tte_covid_hosp[data_cohort_AZ$ind_covid_hosp]
#vec[order(vec)]

# Function to calculate IRRs and return redacted table
redacted_irr_table = function(ind_endpoint, endpoint_date, tte_endpoint) {
  # add outcome date
  data_cohort_irr = data_cohort %>%
    mutate(
      outcome_date = get(endpoint_date), 
      valid_outcome = get(ind_endpoint),
  
      # select minimum date of censorship or outcome - used to exclude from analysis
      exclusion_date = pmin(censor_date, outcome_date, na.rm=TRUE),
      
      # time to censorship or outcome
      tte_exclusion = tte(vax2_date - 1, exclusion_date, exclusion_date),
      
      # person-days contributed to window 1
      persondays_window1 = ifelse(tte_exclusion>postvaxcuts[1] & tte_exclusion<=postvaxcuts[2], tte_exclusion, 56),
      
      # person-days contributed to window 2
      persondays_window2 = ifelse(tte_exclusion>postvaxcuts[2] & tte_exclusion<=postvaxcuts[3], tte_exclusion-postvaxcuts[2], 56),
      persondays_window2 = ifelse(tte_exclusion<=postvaxcuts[2], 0, persondays_window2),
      
      # person-days contributed to window 3
      persondays_window3 = ifelse(tte_exclusion>postvaxcuts[3] & tte_exclusion<=postvaxcuts[4], tte_exclusion-postvaxcuts[3], 56),
      persondays_window3 = ifelse(tte_exclusion<=postvaxcuts[3], 0, persondays_window3),
      
      # person-days contributed to window 4
      persondays_window4 = ifelse(tte_exclusion>postvaxcuts[4] & tte_exclusion<=postvaxcuts[5], tte_exclusion-postvaxcuts[4], 56),
      persondays_window4 = ifelse(tte_exclusion<=postvaxcuts[4], 0, persondays_window4),
      
      # person-days contributed to window 5
      persondays_window5 = ifelse(tte_exclusion>postvaxcuts[5] & tte_exclusion<=postvaxcuts[6], tte_exclusion-postvaxcuts[5], 56),
      persondays_window5 = ifelse(tte_exclusion<=postvaxcuts[5], 0, persondays_window5),
      
      # person-days contributed to window 56 (all)
      persondays_window6 = tte_exclusion
    )
  
  # split into vaccine-specific dataframes
  data_cohort_AZ <- subset(data_cohort_irr, vax2_type=="az")
  data_cohort_BNT <- subset(data_cohort_irr, vax2_type=="pfizer")
  
  # Create IRR dataframe to fill in
  postvax_irr <- data.frame(
    period = c("1-56", "57-112", "113-168", "169-224", "225-280", "All"),
    period_start = c(postvaxcuts[1:5]+1,1),
    period_end = c(postvaxcuts[2:6], 280)
  )
  
# Calculate vaccine-specific N, event rates, person-time, and rates
  for (i in 1:nrow(postvax_irr)) {
    postvax_irr$BNT_n[i] = sum(data_cohort_BNT$tte_exclusion>=postvax_irr$period_start[i])
    postvax_irr$BNT_events[i] = sum(data_cohort_BNT$valid_outcome==TRUE & data_cohort_BNT[,tte_endpoint]>=postvax_irr$period_start[i] & data_cohort_BNT[,tte_endpoint]<=postvax_irr$period_end[i], na.rm = T)
    postvax_irr$BNT_personyears[i] = round(sum(data_cohort_BNT[,paste0("persondays_window",i)])/365.25,2)
    postvax_irr$BNT_rate[i] = round(postvax_irr$BNT_events[i]/postvax_irr$BNT_personyears[i],2)

    postvax_irr$AZ_n[i] = sum(data_cohort_AZ$tte_exclusion>=postvax_irr$period_start[i])
    postvax_irr$AZ_events[i] = sum(data_cohort_AZ$valid_outcome==TRUE & data_cohort_AZ[,tte_endpoint]>=postvax_irr$period_start[i] & data_cohort_AZ[,tte_endpoint]<=postvax_irr$period_end[i], na.rm = T)
    postvax_irr$AZ_personyears[i] = round(sum(data_cohort_AZ[,paste0("persondays_window",i)])/365.25,2)
    postvax_irr$AZ_rate[i] = round(postvax_irr$AZ_events[i]/postvax_irr$AZ_personyears[i]*1000,2)
  } 
  
  # Calculate comparative IRRs
  postvax_irr <- postvax_irr %>%
    mutate(
      rrE = scales::label_number(accuracy=0.01, trim=FALSE)(AZ_rate/BNT_rate),
      rrCI = rrCI_exact(AZ_events, AZ_personyears, BNT_events, BNT_personyears, 0.01),
      rr = paste(rrE,rrCI)
    ) %>%
    select(-rrE, -rrCI)
  
  #Redact values <=5 and corresponding IRRs
  for (i in 1:nrow(postvax_irr)) {
     if (as.numeric(postvax_irr$AZ_n[i])>0 & as.numeric(postvax_irr$AZ_n[i])<=5) { postvax_irr[i,c("AZ_n", "AZ_rate")] = "[Redacted]" }
     if (as.numeric(postvax_irr$BNT_n[i])>0 & as.numeric(postvax_irr$BNT_n[i])<=5) { postvax_irr[i,c("BNT_n", "BNT_rate")] = "[Redacted]" }

     if (as.numeric(postvax_irr$AZ_events[i])>0 & as.numeric(postvax_irr$AZ_events[i])<=5) { postvax_irr[i,c("AZ_events", "AZ_rate")] = "[Redacted]" }
     if (as.numeric(postvax_irr$BNT_events[i])>0 & as.numeric(postvax_irr$BNT_events[i])<=5) { postvax_irr[i,c("BNT_events", "BNT_rate")] = "[Redacted]" }

     if (postvax_irr$AZ_rate[i]=="[Redacted]" | postvax_irr$BNT_rate[i]=="[Redacted]" ) { postvax_irr[i,c("rr")] = "--" }
   }
  
  # Add endpoint
  postvax_irr$outcome = tte_endpoint
  
  # Return redacted results
  postvax_irr
}

irr_collated = rbind(
  redacted_irr_table(ind_endpoint = "ind_covid_postest", endpoint_date = "postvax_positive_test_date", tte_endpoint = "tte_covid_postest"),
  redacted_irr_table(ind_endpoint = "ind_covid_emergency", endpoint_date = "postvax_covid_emergency_date", tte_endpoint = "tte_covid_emergency"),
  redacted_irr_table(ind_endpoint = "ind_covid_hosp", endpoint_date = "postvax_covid_hospitalisation_date", tte_endpoint = "tte_covid_hosp"),
  redacted_irr_table(ind_endpoint = "ind_covid_death", endpoint_date = "postvax_covid_death_date", tte_endpoint = "tte_covid_death")
)

irr_collated$outcome_clean = "Positive SARS-CoV-2 test"
irr_collated$outcome_clean[irr_collated$outcome=="tte_covid_emergency"] = "COVID-related A&E admission"
irr_collated$outcome_clean[irr_collated$outcome=="tte_covid_hosp"] = "COVID-related hospitalisation"
irr_collated$outcome_clean[irr_collated$outcome=="tte_covid_death"] = "COVID-related death"

# Save output
write_csv(irr_collated, here::here("output", "tables", output_csv))
write_rds(irr_collated, here::here("output", "tables", output_rds), compress="gz")
