######################################

# This script:
# - imports processed data
# - fit univariate and multivariable stratified cox model(s) using the coxph package
# - saves models

######################################


### Preliminaries ----

## Import libraries
library('here')
library('tidyr')
library('tidyverse')
library('lubridate')
library('survival')
library('gtsummary')
library('gt')
library('survminer')
library('glue')
library('fs')
library('splines')
sessionInfo()

## Import command-line arguments (specifying whether or not to run matched analysis)
args <- commandArgs(trailingOnly=TRUE)

## Set input and output pathways for matched/unmatched data - default is unmatched
if(length(args)==0){
  # default (unmatched) file names
  db = "VE"
  input_name = "data_cohort_VE.rds"
  irr_name = "table_irr_redacted.rds"
  outcome = "covid_postest"
} else {
  outcome = args[[2]]
  if (args[[1]]=="unmatched") { 
    # unmatched file names    
    db = "VE"
    input_name = "data_cohort_VE.rds"
    irr_name = "table_irr_redacted.rds"
  } else if (args[[1]]=="matched") {
    # matched file names
    db = "VE_matched"
    input_name = "data_cohort_VE_matched.rds"
    irr_name = "table_irr_matched_redacted.rds"
  } else if (args[[1]]=="calendar_time") {
    # file names for calendar time models
    db = "VE_calendar_time"
    input_name = "data_cohort_VE.rds"
    irr_name = "table_irr_redacted.rds" # calculated on person-time
  } else {
    # print error if no argument specified
    print("No matching argument specified")
  }
}

## Import data
data_cohort <- read_rds(here::here("output", "data", input_name))

## Import custom user functions and packages
source(here::here("analysis", "functions.R"))

## Create directory for full model outputs
dir.create(here::here("output", "model", db), showWarnings = FALSE, recursive=TRUE)

## Set analysis intervals and last follow-up day
# TBC
postvaxcuts <- 56*0:5
postvax_periods = c("1-56", "57-112", "113-168", "169-224", "225-280")
lastfupday <- max(postvaxcuts)

## create special log file ----
cat(glue("## script info for cox models ##"), "  \n", file = here::here("output", "model", db, glue("modelcox_log.txt")), append = FALSE)

## function to pass additional log text
logoutput <- function(...){
  cat(..., file = here::here("output", "model", db, glue("modelcox_log.txt")), sep = "\n  ", append = TRUE)
  cat("\n", file = here::here("output", "model", db, glue("modelcox_log.txt")), sep = "\n  ", append = TRUE)
}

### print dataset size ----
logoutput(
  glue("data_cohort data size = ", nrow(data_cohort)),
  glue("data_cohort memory usage = ", format(object.size(data_cohort), units="GB", standard="SI", digits=3L))
)

# create dataset containing one row per patient per post-vaccination period
postvax_time <- data_cohort %>%
  select(patient_id) %>%
  uncount(weights = length(postvaxcuts), .id="id_postvax") %>%
  mutate(
    fup_day = postvaxcuts[id_postvax],
    timesincevax_pw = timesince_cut(fup_day+1, postvaxcuts)
  ) %>%
  droplevels()# %>%
# 6 rows per patient (one per follow-up period)
# id_postvax = comparison period
# fup_day = final day of preceding follow-up period
# timesincevax_pw = days included in follow-up period
# select(patient_id, fup_day, timesincevax_pw)



####################################################### 
### formulae for unadjusted/adjusted models
#######################################################
if (db=="VE") {
  # cox models stratified by follow-up window
  formula0 <- Surv(tstart, tstop, ind_outcome) ~ vax2_az:strata(timesincevax_pw)
  formula1 <- formula0 %>% update(. ~ . + strata(region)*ns(vax2_day, 3))
  formula2 <- formula1 %>% update(. ~ . + poly(age, degree = 2, raw = TRUE) + ckd_5cat + immunosuppression + care_home + sex + imd + ethnicity + 
                                    rural_urban_group + prior_covid_cat + prevax_tests_cat + multimorb + sev_mental_ill)
  
  # cox models for full follow-up time
  formula0_full <- Surv(follow_up_time, ind_outcome) ~ vax2_az
  formula1_full <- formula0_full %>% update(. ~ . + strata(region)*ns(vax2_day, 3))
  formula2_full <- formula1_full %>% update(. ~ . + poly(age, degree = 2, raw = TRUE) + ckd_5cat + immunosuppression + care_home + sex + imd + ethnicity + 
                                              rural_urban_group + prior_covid_cat + prevax_tests_cat + multimorb + sev_mental_ill)
  
} else if (db=="VE_matched") {
  # cox models stratified by follow-up window - matched analysis
  formula0 <- Surv(tstart, tstop, ind_outcome) ~ vax2_az:strata(timesincevax_pw)
  formula1 <- formula0 %>% update(. ~ . + ns(vax2_day, 3)) # no longer need to adjust for region
  formula2 <- formula1 %>% update(. ~ . + sex + imd + ethnicity + rural_urban_group + prevax_tests_cat + multimorb + sev_mental_ill) # no longer need to adjust for age or prior COVID
  
  # cox models for full follow-up time
  formula0_full <- Surv(follow_up_time, ind_outcome) ~ vax2_az
  formula1_full <- formula0_full %>% update(. ~ . + ns(vax2_day, 3)) # no longer need to adjust for region
  formula2_full <- formula1_full %>% update(. ~ . + sex + imd + ethnicity + rural_urban_group + prevax_tests_cat + multimorb + sev_mental_ill) # no longer need to adjust for age or prior COVID
} else if (db=="VE_calendar_timw") {
  #TODO
  # cox models stratified by follow-up window
  formula0 <- Surv(tstart, tstop, ind_outcome) ~ vax2_az:strata(id_postvax)
  formula1 <- formula0 %>% update(. ~ . + strata(strata_var)) # strata_var = region * jcvi_group
  formula2 <- formula1 %>% update(. ~ . + 
                                    poly(age, degree = 2, raw = TRUE) +
                                    # poly(age_1, degree = 2, raw = TRUE) + # care home (no restriction on age)
                                    # poly(age_2, degree = 2, raw = TRUE) + # 80+
                                    # poly(age_3, degree = 1, raw = TRUE) + # 75-79
                                    # poly(age_4a, degree = 1, raw = TRUE) + # 70-74
                                    # poly(age_4b, degree = 2, raw = TRUE) + # 16-69
                                    # poly(age_5, degree = 1, raw = TRUE) + # 65-69
                                    # poly(age_6, degree = 2, raw = TRUE) + # 16-64
                                    ckd_5cat + immunosuppression + 
                                    # care_home + # remove carehome as only jcvi group 1? 
                                    sex + imd + ethnicity + rural_urban_group + prior_covid_cat + prevax_tests_cat + multimorb + sev_mental_ill)
}

####################################################### 
# prepare time to event data for the given outcome
####################################################### 

outcome_list = c("covid_postest", "covid_emergency", "covid_hosp", "covid_death")
clean_list = c("Positive SARS-CoV-2 test", "COVID-related A&E admission", "COVID-related hospitalisation", "COVID-related death")
date_list = c("postvax_positive_test_date", "postvax_covid_emergency_date", "postvax_covid_hospitalisation_date", "postvax_covid_death_date")
i = which(outcome_list == outcome)

# read in incidence rate ratio table
irr_table = read_rds(here::here("output", "tables", irr_name))


  
  selected_outcome = outcome_list[i]
  selected_outcome_clean = clean_list[i]
  irr_sub = subset(irr_table, outcome_clean==selected_outcome_clean)[,c("period", "BNT_n", "BNT_events", "AZ_n", "AZ_events")]
  
  data_tte <- data_cohort %>% 
    mutate(
      # select dates for outcome in question
      outcome_date = get(date_list[i]),
      
      # censor date already defined in data_selection_VE.R script 
      
      # calculate tte and ind for outcome in question
      tte_outcome = tte(vax2_date-1, outcome_date, censor_date, na.censor=TRUE),
      ind_outcome = get(paste0("ind_",selected_outcome)),
      tte_stop = pmin(tte_censor, tte_outcome, na.rm=TRUE),
      
      # calculate follow-up time (censor/event)
      follow_up_time = tte(vax2_date-1, get(date_list[i]), censor_date) 
    )
  
  data_cox <- tmerge(
    data1 = data_tte %>% select(-starts_with("ind_"), -ends_with("_date")),
    data2 = data_tte,
    id = patient_id,
    tstart = 0L,
    tstop = pmin(tte_censor, tte_outcome, na.rm=TRUE),
    ind_outcome = event(tte_outcome)
  ) %>%
    tmerge( # create treatment timescale variables
      data1 = .,
      data2 = postvax_time,
      id = patient_id,
      timesincevax_pw = tdc(fup_day, timesincevax_pw)
    )
  # set factor levels for postvaccination periods
  data_cox$timesincevax_pw = factor(data_cox$timesincevax_pw, levels = postvax_periods)
  
  # derive strata_var if using calendar timescale
  if (db == "VE_calendar_time") {
    data_cox <- data_cox %>% mutate(strata_var = as.character(glue("{jcvi_group}_{region}")))
  }
  
  ### print dataset size and save ----
  logoutput(
    glue(selected_outcome_clean, "\ndata_cox data size = ", nrow(data_cox)),
    glue("data_cox memory usage = ", format(object.size(data_cox), units="GB", standard="SI", digits=3L))
  )
  

