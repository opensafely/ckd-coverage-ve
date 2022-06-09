# # # # # # # # # # # # # # # # # # # # #
# This script creates a table of incidence rates for diffferent outcomes, stratified by vaccine type and time since vaccination
# # # # # # # # # # # # # # # # # # # # #

## Import libraries ----
library('tidyverse')
library('here')
library('glue')
library('survival')
library('gt')
library('gtsummary')
library('scales')
library('lubridate')

## Set rounding (TRUE/FALSE) and threshold
round_logical = TRUE
round_threshold = 5
redaction_threshold = 10

## Import command-line arguments (specifying whether or not to run matched analysis)
args <- commandArgs(trailingOnly=TRUE)

## Set input and output pathways for matched/unmatched data - default is unmatched
if(length(args)==0){
  ## Default (unmatched) file names
  input_name = "data_cohort_VE.rds"
  output_csv = "table_irr_redacted.csv"
  output_rds = "table_irr_redacted.rds"
} else {
  if (args[[1]]=="unmatched") { 
    ## Unmatched file names    
    input_name = "data_cohort_VE.rds"
    output_csv = "table_irr_redacted.csv"
    output_rds = "table_irr_redacted.rds"
  } else if (args[[1]]=="matched") {
    ## Matched file names
    input_name = "data_cohort_VE_matched.rds"
    output_csv = "table_irr_matched_redacted.csv"
    output_rds = "table_irr_matched_redacted.rds"  
  } else {
    ## Print error if no argument specified
    print("No matching argument specified")
  }
}

## Import custom user functions
source(here::here("analysis", "functions.R"))

## Set analysis intervals
postvaxcuts <- 28*0:6
window_length = 28
lastfupday <- max(postvaxcuts)

## Import processed data
data_cohort <- read_rds(here::here("output", "data", input_name))

## Split into vaccine-specific data frames
data_cohort_AZ <- subset(data_cohort, vax2_type=="az")
data_cohort_BNT <- subset(data_cohort, vax2_type=="pfizer")

## Function to calculate IRRs and return redacted table
redacted_irr_table = function(ind_endpoint, endpoint_date, tte_endpoint) {
  ## Add outcome date
  data_cohort_irr = data_cohort %>%
    mutate(
      outcome_date = get(endpoint_date), 
      valid_outcome = get(ind_endpoint),
  
      ## Select minimum date of censorship or outcome - used to exclude from analysis
      exclusion_date = pmin(censor_date, outcome_date, na.rm=TRUE),
      
      ## Time to censorship or outcome
      tte_exclusion = tte(vax2_date - 1, exclusion_date, exclusion_date),
      
      ## Person-days contributed to window 1
      persondays_window1 = ifelse(tte_exclusion>postvaxcuts[1] & tte_exclusion<=postvaxcuts[2], tte_exclusion, window_length),
      
      ## Person-days contributed to window 2
      persondays_window2 = ifelse(tte_exclusion>postvaxcuts[2] & tte_exclusion<=postvaxcuts[3], tte_exclusion-postvaxcuts[2], window_length),
      persondays_window2 = ifelse(tte_exclusion<=postvaxcuts[2], 0, persondays_window2),
      
      ## Person-days contributed to window 3
      persondays_window3 = ifelse(tte_exclusion>postvaxcuts[3] & tte_exclusion<=postvaxcuts[4], tte_exclusion-postvaxcuts[3], window_length),
      persondays_window3 = ifelse(tte_exclusion<=postvaxcuts[3], 0, persondays_window3),
      
      ## Person-days contributed to window 4
      persondays_window4 = ifelse(tte_exclusion>postvaxcuts[4] & tte_exclusion<=postvaxcuts[5], tte_exclusion-postvaxcuts[4], window_length),
      persondays_window4 = ifelse(tte_exclusion<=postvaxcuts[4], 0, persondays_window4),
      
      ## Person-days contributed to window 5
      persondays_window5 = ifelse(tte_exclusion>postvaxcuts[5] & tte_exclusion<=postvaxcuts[6], tte_exclusion-postvaxcuts[5], window_length),
      persondays_window5 = ifelse(tte_exclusion<=postvaxcuts[5], 0, persondays_window5),
      
      ## Person-days contributed to window 6
      persondays_window6 = ifelse(tte_exclusion>postvaxcuts[6] & tte_exclusion<=postvaxcuts[7], tte_exclusion-postvaxcuts[6], 56),
      persondays_window6 = ifelse(tte_exclusion<=postvaxcuts[6], 0, persondays_window6),

      ## Person-days contributed to window 7 (all follow-up)
      persondays_window7 = tte_exclusion
    )
  
  ## Split into vaccine-specific data frames
  data_cohort_AZ <- subset(data_cohort_irr, vax2_type=="az")
  data_cohort_BNT <- subset(data_cohort_irr, vax2_type=="pfizer")
  
  ## Create IRR dataframe to fill in
  postvax_irr <- data.frame(
    period = c("1-28", "29-56", "57-84", "85-112", "113-140", "141-168", "1-168"),
    period_start = c(postvaxcuts[1:6]+1,1),
    period_end = c(postvaxcuts[2:7], 168)
  )
  
## Calculate vaccine-specific N, event rates, person-time, and rates
  for (i in 1:nrow(postvax_irr)) {
    postvax_irr$BNT_n[i] = sum(data_cohort_BNT$tte_exclusion>=postvax_irr$period_start[i])
    postvax_irr$BNT_events[i] = sum(data_cohort_BNT$valid_outcome==TRUE & data_cohort_BNT[,tte_endpoint]>=postvax_irr$period_start[i] & data_cohort_BNT[,tte_endpoint]<=postvax_irr$period_end[i], na.rm = T)
    postvax_irr$BNT_personyears[i] = round(sum(data_cohort_BNT[,paste0("persondays_window",i)])/365.25,2)
    postvax_irr$BNT_rate[i] = round(postvax_irr$BNT_events[i]/postvax_irr$BNT_personyears[i]*1000,2)

    postvax_irr$AZ_n[i] = sum(data_cohort_AZ$tte_exclusion>=postvax_irr$period_start[i])
    postvax_irr$AZ_events[i] = sum(data_cohort_AZ$valid_outcome==TRUE & data_cohort_AZ[,tte_endpoint]>=postvax_irr$period_start[i] & data_cohort_AZ[,tte_endpoint]<=postvax_irr$period_end[i], na.rm = T)
    postvax_irr$AZ_personyears[i] = round(sum(data_cohort_AZ[,paste0("persondays_window",i)])/365.25,2)
    postvax_irr$AZ_rate[i] = round(postvax_irr$AZ_events[i]/postvax_irr$AZ_personyears[i]*1000,2)
  } 
  
  if (round_logical==TRUE) {
    postvax_irr$BNT_n = plyr::round_any(postvax_irr$BNT_n,round_threshold)
    postvax_irr$BNT_events = plyr::round_any(postvax_irr$BNT_events,round_threshold)
    postvax_irr$BNT_rate[i] = round(postvax_irr$BNT_events[i]/postvax_irr$BNT_personyears[i]*1000,2)
    postvax_irr$AZ_n = plyr::round_any(postvax_irr$AZ_n,round_threshold)
    postvax_irr$AZ_events = plyr::round_any(postvax_irr$AZ_events,round_threshold)
    postvax_irr$AZ_rate[i] = round(postvax_irr$AZ_events[i]/postvax_irr$AZ_personyears[i]*1000,2)
  }
  
  ## Calculate comparative IRRs
  postvax_irr <- postvax_irr %>%
    mutate(
      rrE = scales::label_number(accuracy=0.01, trim=FALSE)(AZ_rate/BNT_rate),
      rrCI = rrCI_exact(AZ_events, AZ_personyears, BNT_events, BNT_personyears, 0.01),
      rr = paste(rrE,rrCI)
    ) %>%
    select(-rrE, -rrCI)
  
  ## Redact values and corresponding IRRs for non-zero counts under redaction threshold
  for (i in 1:nrow(postvax_irr)) {
     if (as.numeric(postvax_irr$AZ_n[i])>0 & as.numeric(postvax_irr$AZ_n[i])<=redaction_threshold) { postvax_irr[i,c("AZ_n", "AZ_rate")] = "[Redacted]" }
     if (as.numeric(postvax_irr$BNT_n[i])>0 & as.numeric(postvax_irr$BNT_n[i])<=redaction_threshold) { postvax_irr[i,c("BNT_n", "BNT_rate")] = "[Redacted]" }

     if (as.numeric(postvax_irr$AZ_events[i])>0 & as.numeric(postvax_irr$AZ_events[i])<=redaction_threshold) { postvax_irr[i,c("AZ_events", "AZ_rate")] = "[Redacted]" }
     if (as.numeric(postvax_irr$BNT_events[i])>0 & as.numeric(postvax_irr$BNT_events[i])<=redaction_threshold) { postvax_irr[i,c("BNT_events", "BNT_rate")] = "[Redacted]" }

     if (postvax_irr$AZ_rate[i]=="[Redacted]" | postvax_irr$BNT_rate[i]=="[Redacted]" ) { postvax_irr[i,c("rr")] = "--" }
     if (postvax_irr$AZ_events[i]==0 & postvax_irr$BNT_events[i]==0 ) { postvax_irr[i,c("rr")] = "--" }
  }
  
  ## Add endpoint
  postvax_irr$outcome = tte_endpoint
  
  ## Return redacted results
  postvax_irr
}

## Collate IRR table based on function above
irr_collated = rbind(
  redacted_irr_table(ind_endpoint = "ind_covid_postest", endpoint_date = "postvax_positive_test_date", tte_endpoint = "tte_covid_postest"),
  redacted_irr_table(ind_endpoint = "ind_covid_emergency", endpoint_date = "postvax_covid_emergency_date", tte_endpoint = "tte_covid_emergency"),
  redacted_irr_table(ind_endpoint = "ind_covid_hosp", endpoint_date = "postvax_covid_hospitalisation_date", tte_endpoint = "tte_covid_hosp"),
  redacted_irr_table(ind_endpoint = "ind_covid_death", endpoint_date = "postvax_covid_death_date", tte_endpoint = "tte_covid_death")
)

## Add clean names
irr_collated$outcome_clean = "Positive SARS-CoV-2 test"
irr_collated$outcome_clean[irr_collated$outcome=="tte_covid_emergency"] = "COVID-related A&E admission"
irr_collated$outcome_clean[irr_collated$outcome=="tte_covid_hosp"] = "COVID-related hospitalisation"
irr_collated$outcome_clean[irr_collated$outcome=="tte_covid_death"] = "COVID-related death"

## Remove event counts and person-years for summary metric so that redacted counts cannot be back-calculated
irr_collated$BNT_events[irr_collated$period=="1-168"] = "--"
irr_collated$BNT_personyears[irr_collated$period=="1-168"] = "--"
irr_collated$AZ_events[irr_collated$period=="1-168"] = "--"
irr_collated$AZ_personyears[irr_collated$period=="1-168"] = "--"

## Save output
write_csv(irr_collated, here::here("output", "tables", output_csv))
write_rds(irr_collated, here::here("output", "tables", output_rds), compress="gz")
