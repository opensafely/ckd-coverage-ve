# # # # # # # # # # # # # # # # # # # # #
# This script creates a table of incidence rates for diffferent outcomes, stratified by vaccine type and time since vaccination
# # # # # # # # # # # # # # # # # # # # #

## Import libraries 
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

if(length(args)==0){
  # default (unmatched) file names
  matching_status = "unmatched"
  subgroup = "all"
  vaccine = "primary"
} else {
  matching_status = args[[1]] # can be unmatched or matched
  subgroup = args[[2]] # can be all / CKD3 / CKD4-5 / RRT / Tx / dialysis
  vaccine = args[[3]] # can be primary or boost
}

## Import data
if (matching_status=="unmatched" & vaccine=="primary") { 
  data_cohort <- read_rds(here::here("output", "data", "data_cohort_VE_primary.rds"))
  
} else if (matching_status=="matched" & vaccine=="primary") { 
  data_cohort <- read_rds(here::here("output", "data", "data_cohort_VE_primary_matched.rds"))
  
} else if (matching_status=="unmatched" & vaccine=="boost") { 
  data_cohort <- read_rds(here::here("output", "data", "data_cohort_VE_boost.rds"))
  
} else { 
  data_cohort <- read_rds(here::here("output", "data", "data_cohort_VE_boost_matched.rds"))
}

# Modify dates to avoid timestamps causing inequalities for dates on the same day
data_cohort <- data_cohort %>%
  mutate(across(where(is.Date), 
                ~ floor_date(
                  as.Date(.x, format="%Y-%m-%d"),
                  unit = "days")))

## Specify output path names
output_csv = paste0("table_irr_",vaccine,"_redacted_",matching_status,"_",subgroup,".csv")
output_rds = paste0("table_irr_",vaccine,"_redacted_",matching_status,"_",subgroup,".rds")
fs::dir_create(here::here("output", "tables"))

## Select subset
if (subgroup=="all") {
  data_cohort = data_cohort
  ## Kidney disease subgroups
} else if (subgroup=="CKD3") {
  data_cohort = subset(data_cohort, ckd_3cat == "CKD3")
} else if (subgroup=="CKD4-5") {
  data_cohort = subset(data_cohort, ckd_3cat == "CKD4-5")
} else if (subgroup=="RRT") {
  data_cohort = subset(data_cohort, ckd_3cat == "RRT (any)")
} else if (subgroup=="Tx") {
  data_cohort = subset(data_cohort, ckd_5cat == "RRT (Tx)")
} else if (subgroup=="dialysis") {
  data_cohort = subset(data_cohort, ckd_5cat == "RRT (dialysis)")
} else {
  stop ("Arguments not specified correctly.")
}

## Import custom user functions
source(here::here("analysis", "functions.R"))

## Import outcome time periods
if (vaccine=="primary") {
  postvax_list <- read_rds(
    here::here("output", "lib", "postvax_list.rds")
  )
} else {
  postvax_list <- read_rds(
    here::here("output", "lib", "postboost_list.rds")
  )
}
postvaxcuts = postvax_list$postvaxcuts
postvax_periods = postvax_list$postvax_periods
period_length_1 = postvaxcuts[2]-postvaxcuts[1]
period_length = postvaxcuts[3]-postvaxcuts[2]
lastfupday = max(postvaxcuts)

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
      stop_date = pmin(censor_date, outcome_date, na.rm=TRUE),
    )
  
  ## Specify time to censorship or outcome
  if (vaccine=="primary") {
    data_cohort_irr = data_cohort_irr %>%
      mutate(
        tte_stop = tte(vax2_date - 1, stop_date, na.censor=TRUE),
      )
  } else {
    data_cohort_irr = data_cohort_irr %>%
      mutate(
        tte_stop = tte(vax3_date - 1, stop_date, na.censor=TRUE),
      )
  }
  
  ## Calculate person-days per window  
  data_cohort_irr = data_cohort_irr %>%
    mutate(
      ## Person-days contributed to window 1
      persondays_window1 = ifelse(tte_stop>postvaxcuts[1] & tte_stop<=postvaxcuts[2], tte_stop, period_length_1),
      
      ## Person-days contributed to window 2
      persondays_window2 = ifelse(tte_stop>postvaxcuts[2] & tte_stop<=postvaxcuts[3], tte_stop-postvaxcuts[2], period_length),
      persondays_window2 = ifelse(tte_stop<=postvaxcuts[2], 0, persondays_window2),
      
      ## Person-days contributed to window 3
      persondays_window3 = ifelse(tte_stop>postvaxcuts[3] & tte_stop<=postvaxcuts[4], tte_stop-postvaxcuts[3], period_length),
      persondays_window3 = ifelse(tte_stop<=postvaxcuts[3], 0, persondays_window3),

      ## Person-days contributed to window 4
      persondays_window4 = ifelse(tte_stop>postvaxcuts[4] & tte_stop<=postvaxcuts[5], tte_stop-postvaxcuts[4], period_length),
      persondays_window4 = ifelse(tte_stop<=postvaxcuts[4], 0, persondays_window4),

      ## Person-days contributed to window 5 (all follow-up)
      persondays_window5 = tte_stop
    )

  ## Split into vaccine-specific data frames
  data_cohort_AZ <- subset(data_cohort_irr, vax2_type=="az")
  data_cohort_BNT <- subset(data_cohort_irr, vax2_type=="pfizer")
  
  ## Create IRR dataframe to fill in
  postvax_irr <- data.frame(
      period = c(postvax_periods, paste0(postvaxcuts[1]+1,"-",lastfupday)),
      period_start = c(postvaxcuts[1:4]+1,postvaxcuts[1]+1),
      period_end = c(postvaxcuts[2:5], lastfupday)
  )

  
  ## Calculate vaccine-specific N, event rates, person-time, and rates
  for (i in 1:nrow(postvax_irr)) {
    postvax_irr$BNT_n[i] = sum(data_cohort_BNT$tte_stop>=postvax_irr$period_start[i])
    postvax_irr$BNT_events[i] = sum(data_cohort_BNT$valid_outcome==TRUE & data_cohort_BNT[,tte_endpoint]>=postvax_irr$period_start[i] & data_cohort_BNT[,tte_endpoint]<=postvax_irr$period_end[i], na.rm = T)
    postvax_irr$BNT_personyears[i] = round(sum(data_cohort_BNT[,paste0("persondays_window",i)])/365.25,2)
    postvax_irr$BNT_rate[i] = round(postvax_irr$BNT_events[i]/postvax_irr$BNT_personyears[i]*1000,2)
    
    postvax_irr$AZ_n[i] = sum(data_cohort_AZ$tte_stop>=postvax_irr$period_start[i])
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
  #redacted_irr_table(ind_endpoint = "ind_covid_emergency", endpoint_date = "postvax_covid_emergency_date", tte_endpoint = "tte_covid_emergency"),
  redacted_irr_table(ind_endpoint = "ind_covid_hosp", endpoint_date = "postvax_covid_hospitalisation_date", tte_endpoint = "tte_covid_hosp"),
  redacted_irr_table(ind_endpoint = "ind_covid_death", endpoint_date = "postvax_covid_death_date", tte_endpoint = "tte_covid_death"),
  redacted_irr_table(ind_endpoint = "ind_noncovid_death", endpoint_date = "noncoviddeath_date", tte_endpoint = "tte_noncovid_death")
)

## Add clean names
irr_collated$outcome_clean = "Positive SARS-CoV-2 test"
#irr_collated$outcome_clean[irr_collated$outcome=="tte_covid_emergency"] = "COVID-related A&E admission"
irr_collated$outcome_clean[irr_collated$outcome=="tte_covid_hosp"] = "COVID-19-related hospitalisation"
irr_collated$outcome_clean[irr_collated$outcome=="tte_covid_death"] = "COVID-19-related death"
irr_collated$outcome_clean[irr_collated$outcome=="tte_noncovid_death"] = "Non-COVID-19 death"

## Simplify output to report on whole study period in subgroup analyses
if (subgroup!="all") {
  irr_collated = subset(irr_collated, period=="1-182")
}

## Save output
write_csv(irr_collated, here::here("output", "tables", output_csv))
write_rds(irr_collated, here::here("output", "tables", output_rds), compress="gz")