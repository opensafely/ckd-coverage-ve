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

## Import custom user functions
source(here::here("analysis", "functions.R"))

## Set analysis intervals
postvaxcuts <- 56*0:5
lastfupday <- max(postvaxcuts)

## Import processed data ----
data_cohort <- read_rds(here("output", "data", "data_cohort_VE.rds"))

## Set analysis end date
data_cohort$end_date = as_date("2021-11-30")

## Create tte data ----
data_tte <- data_cohort %>%
  transmute(
    patient_id,
    vax2_date,
    vax2_type,
    postvax_positive_test_date,
    postvax_covid_emergency_date,
    postvax_covid_hospitalisation_date,
    postvax_covid_death_date,
    end_date,
    
    # set censor date
    censor_date = pmin(vax2_date + lastfupday, end_date, dereg_date, death_date, vax3_date, na.rm=TRUE),

    # time to last follow up day
    # tte function format: origin date, event date, censor date
    tte_enddate = tte(origin_date = vax2_date, # removed -1 here and below
                      event_date = end_date, 
                      censor_date = end_date),
    
    # time to last follow up day or death or deregistration
    tte_censor = tte(origin_date = vax2_date,
                     event_date = censor_date, 
                     censor_date = censor_date),
    
    # time to positive test
    tte_covid_postest = tte(origin_date = vax2_date, 
                            event_date = postvax_positive_test_date, 
                            censor_date = censor_date, na.censor=TRUE),
    ind_covid_postest = dplyr::if_else((postvax_positive_test_date>censor_date) | is.na(postvax_positive_test_date), FALSE, TRUE),
    
    # time to COVID-19 A&E attendance
    tte_covid_emergency = tte(origin_date = vax2_date, 
                              event_date = postvax_covid_emergency_date, 
                              censor_date  = censor_date, na.censor=TRUE),
    ind_covid_emergency = dplyr::if_else((postvax_covid_emergency_date>censor_date) | is.na(postvax_covid_emergency_date), FALSE, TRUE),
    
    # time to COVID-19 hospitalisation
    tte_covid_hosp = tte(origin_date = vax2_date, 
                         event_date = postvax_covid_hospitalisation_date, 
                         censor_date = censor_date, na.censor=TRUE),
    ind_covid_hosp = dplyr::if_else((postvax_covid_hospitalisation_date>censor_date) | is.na(postvax_covid_hospitalisation_date), FALSE, TRUE),
    
    # time to COVID-19 death
    tte_covid_death = tte(origin_date = vax2_date, 
                          event_date = postvax_covid_death_date, 
                          censor_date = censor_date, na.censor=TRUE),
    ind_covid_death = dplyr::if_else((postvax_covid_death_date>censor_date) | is.na(postvax_covid_death_date), FALSE, TRUE)
  ) %>%
   filter(
     # filter if censor event occurs before second dose
     tte_censor>0 | is.na(tte_censor)
   ) 

# Split into vaccine-specific dataframes
data_tte_AZ <- subset(data_tte, vax2_type=="az")
data_tte_BNT <- subset(data_tte, vax2_type=="pfizer")

# Function to calculate IRRs and return redacted table
redacted_irr_table = function(ind_endpoint, endpoint_date, tte_endpoint) {
  # add outcome date
  data_tte_irr = data_tte %>%
    rename(outcome_date = all_of(endpoint_date), 
           valid_outcome = all_of(ind_endpoint)) %>%
    mutate(
      # select minimum date of censorship or outcome - used to exclude from analysis
      exclusion_date = pmin(censor_date, outcome_date, na.rm=TRUE),
      
      # time to censorship or outcome
      tte_exclusion = tte(origin_date = vax2_date, event_date = exclusion_date, censor_date = exclusion_date),
      
      # person-days contributed to window 1
      persondays_window1 = ifelse(tte_exclusion<postvaxcuts[2], tte_exclusion, 56),
      
      # person-days contributed to window 2
      persondays_window2 = ifelse(tte_exclusion>=postvaxcuts[2] & tte_exclusion<postvaxcuts[3], tte_exclusion-postvaxcuts[2], 56),
      persondays_window2 = ifelse(tte_exclusion<postvaxcuts[2], 0, persondays_window2),
      
      # person-days contributed to window 3
      persondays_window3 = ifelse(tte_exclusion>=postvaxcuts[3] & tte_exclusion<postvaxcuts[4], tte_exclusion-postvaxcuts[3], 56),
      persondays_window3 = ifelse(tte_exclusion<postvaxcuts[3], 0, persondays_window3),
      
      # person-days contributed to window 4
      persondays_window4 = ifelse(tte_exclusion>=postvaxcuts[4] & tte_exclusion<postvaxcuts[5], tte_exclusion-postvaxcuts[4], 56),
      persondays_window4 = ifelse(tte_exclusion<postvaxcuts[4], 0, persondays_window4),
      
      # person-days contributed to window 5
      persondays_window5 = ifelse(tte_exclusion>=postvaxcuts[5] & tte_exclusion<postvaxcuts[6], tte_exclusion-postvaxcuts[5], 56),
      persondays_window5 = ifelse(tte_exclusion<postvaxcuts[5], 0, persondays_window5),
      
      # person-days contributed to window 56 (all)
      persondays_window6 = tte_exclusion
    )
  
  # split into vaccine-specific dataframes
  data_tte_AZ <- subset(data_tte_irr, vax2_type=="az")
  data_tte_BNT <- subset(data_tte_irr, vax2_type=="pfizer")
  
  # Create IRR dataframe to fill in
  postvax_irr <- data.frame(
    period = c("0-55", "56-111", "112-167", "168-223", "224-279", "All"),
    period_start = c(postvaxcuts[1:5],0),
    period_end = c(postvaxcuts[2:6], 280)
  )
  
# Calculate vaccine-specific N, event rates, person-time, and rates
  for (i in 1:nrow(postvax_irr)) {
    postvax_irr$AZ_n[i] = sum(data_tte_AZ$tte_exclusion>=postvax_irr$period_start[i])
    postvax_irr$AZ_events[i] = sum(data_tte_AZ$valid_outcome==TRUE & data_tte_AZ[,tte_endpoint]>=postvax_irr$period_start[i] & data_tte_AZ[,tte_endpoint]<postvax_irr$period_end[i], na.rm = T)
    postvax_irr$AZ_personyears[i] = sum(data_tte_AZ[,paste0("persondays_window",i)])/365.25
    postvax_irr$AZ_rate[i] = postvax_irr$AZ_events[i]/postvax_irr$AZ_personyears[i]
    
    postvax_irr$BNT_n[i] = sum(data_tte_BNT$tte_exclusion>=postvax_irr$period_start[i])
    postvax_irr$BNT_events[i] = sum(data_tte_BNT$valid_outcome==TRUE & data_tte_BNT[,tte_endpoint]>=postvax_irr$period_start[i] & data_tte_BNT[,tte_endpoint]<postvax_irr$period_end[i], na.rm = T)
    postvax_irr$BNT_personyears[i] = sum(data_tte_BNT[,paste0("persondays_window",i)])/365.25
    postvax_irr$BNT_rate[i] = postvax_irr$BNT_events[i]/postvax_irr$BNT_personyears[i]
  } 
  
  # Calculate comparative IRRs
  postvax_irr <- postvax_irr %>%
    mutate(
      rrE = scales::label_number(accuracy=0.01, trim=FALSE)(AZ_rate/BNT_rate),
      rrCI = rrCI_exact(AZ_events, AZ_personyears, BNT_events, BNT_personyears, 0.01)
    )
  
  #Redact values <=5 and corresponding IRRs
  for (i in 1:nrow(postvax_irr)) {
     if (as.numeric(postvax_irr$AZ_n[i])>0 & as.numeric(postvax_irr$AZ_n[i])<=5) { postvax_irr[i,c("AZ_n", "AZ_rate")] = "[Redacted]" }
     if (as.numeric(postvax_irr$BNT_n[i])>0 & as.numeric(postvax_irr$BNT_n[i])<=5) { postvax_irr[i,c("BNT_n", "BNT_rate")] = "[Redacted]" }
     
     if (as.numeric(postvax_irr$AZ_events[i])>0 & as.numeric(postvax_irr$AZ_events[i])<=5) { postvax_irr[i,c("AZ_events", "AZ_rate")] = "[Redacted]" }
     if (as.numeric(postvax_irr$BNT_events[i])>0 & as.numeric(postvax_irr$BNT_events[i])<=5) { postvax_irr[i,c("BNT_events", "BNT_rate")] = "[Redacted]" }
     
     if (postvax_irr$AZ_rate[i]=="[Redacted]" | postvax_irr$BNT_rate[i]=="[Redacted]" ) { postvax_irr[i,12:13] = "--" }
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
write_csv(irr_collated, here("output", "tables", "table_irr_redacted.csv"))

