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
} else {
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
  droplevels() %>%
  # 6 rows per patient (one per follow-up period)
  # fup_day = final day of preceding follow-up period
  # timesincevax_pw = days included in follow-up period
  select(patient_id, fup_day, timesincevax_pw)



####################################################### 
### formulae for unadjusted/adjusted models
#######################################################
if (db=="VE") {
  # cox models stratified by follow-up window
  formula0 <- Surv(tstart, tstop, ind_outcome) ~ vax2_az:strata(timesincevax_pw)
  formula1 <- formula0 %>% update(. ~ . + strata(region)*ns(vax2_day, 3))
  formula2 <- formula1 %>% update(. ~ . + poly(age, degree = 2, raw = TRUE) + ckd_7cat + immunosuppression + care_home + sex + imd + ethnicity + 
                                    rural_urban_group + prior_covid_cat + prevax_tests_cat + multimorb + sev_mental_ill)
  
  # cox models for full follow-up time
  formula0_full <- Surv(follow_up_time, ind_outcome) ~ vax2_az
  formula1_full <- formula0_full %>% update(. ~ . + strata(region)*ns(vax2_day, 3))
  formula2_full <- formula1_full %>% update(. ~ . + poly(age, degree = 2, raw = TRUE) + ckd_7cat + immunosuppression + care_home + sex + imd + ethnicity + 
                                              rural_urban_group + prior_covid_cat + prevax_tests_cat + multimorb + sev_mental_ill)

} else {
  # cox models stratified by follow-up window - matched analysis
  formula0 <- Surv(tstart, tstop, ind_outcome) ~ vax2_az:strata(timesincevax_pw)
  formula1 <- formula0 %>% update(. ~ . + ns(vax2_day, 3)) # no longer need to adjust for region
  formula2 <- formula1 %>% update(. ~ . + sex + imd + ethnicity + rural_urban_group + prevax_tests_cat + multimorb + sev_mental_ill) # no longer need to adjust for age or prior COVID
  
  # cox models for full follow-up time
  formula0_full <- Surv(follow_up_time, ind_outcome) ~ vax2_az
  formula1_full <- formula0_full %>% update(. ~ . + ns(vax2_day, 3)) # no longer need to adjust for region
  formula2_full <- formula1_full %>% update(. ~ . + sex + imd + ethnicity + rural_urban_group + prevax_tests_cat + multimorb + sev_mental_ill) # no longer need to adjust for age or prior COVID
}



####################################################### 
### function to fit cox model for specified formula 
####################################################### 

#formula_cox = formula0
#number = 0

cox_model_VE <- function(number, formula_cox, stratified=TRUE) {
  if (number==0) { model_type = "unadjusted" } 
  if (number==1) { model_type = "region/date adjusted" } 
  if (number==2) { model_type = "fully adjusted" } 
  
  if (stratified) {
    # if stratified = FALSE, fit cox model stratified by follow-up window
    coxmod <- coxph(
      formula = formula_cox,
      data = data_cox, # input data
      robust = TRUE, # compute robust variance
      id = patient_id, # required since multiple rows per subject
      na.action = "na.fail",
      control = coxph.control(iter.max = 50)
    )
  } else {
    # if stratified = FALSE, fit model on full dataset
    coxmod <- coxph(
      formula_cox, 
      data = data_tte,
      control = coxph.control(iter.max = 50)
      )
  }
  # print warnings
  print(warnings())
  
  # print output status to log file
  logoutput(
    glue("model{number} data size = ", coxmod$n),
    glue("model{number} memory usage = ", format(object.size(coxmod), units="GB", standard="SI", digits=3L))
  )
  
  if (stratified) {
    logoutput(
      glue("convergence status: ", coxmod$info[["convergence"]])
    )
  }
  
  # return model summary
  tidy <- broom.helpers::tidy_plus_plus(coxmod, exponentiate = FALSE) %>%
    add_column(model_name = model_type, .before=1) %>%
    add_column(model = number, .before=1)
  
  # return brief model summary
  glance <- broom::glance(coxmod) %>%
    add_column(
      model_name = model_type,
      model = number,
      ram = format(object.size(coxmod), units="GB", standard="SI", digits=3L),
      .before = 1
    )
  
  if (stratified) {
    tidy$level = glance$level = "stratified"
    glance$convergence = coxmod$info[["convergence"]]
  } else {
    tidy$level = glance$level = "full"
    glance$convergence = NA
  }
  
  # model outputs
  coxmod$data <- NULL
  #write_rds(coxmod, here("output", "model","VE", glue("modelcox_model{number}_",selected_outcome,".rds")), compress="gz")
  lst(glance, tidy)
}




####################################################### 
# loop to fit cox model across outcome
####################################################### 

outcome_list = c("covid_postest", "covid_emergency", "covid_hosp", "covid_death")
clean_list = c("Positive SARS-CoV-2 test", "COVID-related A&E admission", "COVID-related hospitalisation", "COVID-related death")
date_list = c("postvax_positive_test_date", "postvax_covid_emergency_date", "postvax_covid_hospitalisation_date", "postvax_covid_death_date")

# read in incidence rate ratio table
irr_table = read_rds(here::here("output", "tables", irr_name))

for (i in 1:length(outcome_list)) {
  selected_outcome = outcome_list[i]
  selected_outcome_clean = clean_list[i]
  irr_sub = subset(irr_table, outcome_clean==selected_outcome_clean)[,c("period", "BNT_n", "BNT_events", "AZ_n", "AZ_events")]
  
  data_tte <- data_cohort %>% 
    mutate(
      # select dates for outcome in question
      outcome_date = get(date_list[i]),
      
      # censor date already defined in data_selection_VE.R script 
      
      # calculate tte and ind for outcome in question
      tte_outcome = tte(vax2_date-1, outcome_date, censor_date, na.censor=FALSE),
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
  
  ### print dataset size and save ----
  logoutput(
    glue(selected_outcome_clean, "\ndata_cox data size = ", nrow(data_cox)),
    glue("data_cox memory usage = ", format(object.size(data_cox), units="GB", standard="SI", digits=3L))
  )

  # run unadjusted and adjusted models, stratified and unstratified
  #assign("last.warning", NULL, envir = baseenv()) # clear warnings
  summary0 <- cox_model_VE(0, formula0)
  summary1 <- cox_model_VE(1, formula1)
  summary2 <- cox_model_VE(2, formula2)
  summary0_full <- cox_model_VE(0, formula0_full, stratified=FALSE)
  summary1_full <- cox_model_VE(1, formula1_full, stratified=FALSE)
  summary2_full <- cox_model_VE(2, formula2_full, stratified=FALSE)
  
  ## Combine and save model summary (brief) 
  model_glance <- data.frame(
    bind_rows(summary0$glance, summary0_full$glance, 
              summary1$glance, summary1_full$glance, 
              summary2$glance, summary2_full$glance) %>%
    mutate(outcome = selected_outcome, 
           outcome_clean = selected_outcome_clean)
  )
  
  ## Redact statistical outputs if <=10 events
  redaction_columns = c("nevent", "statistic.log", "p.value.log", "statistic.sc", "p.value.sc", "statistic.wald", "p.value.wald", 
                        "statistic.robust", "p.value.robust", "r.squared", "r.squared.max", "concordance", "std.error.concordance", "logLik", "AIC", "BIC")
  for (i in 1:nrow(model_glance)) {
    if (model_glance$nevent[i]>0 & model_glance$nevent[i]<=10) { model_glance[i,names(model_glance)%in%redaction_columns] = "[Redacted]" }
  }
  write_csv(model_glance, here::here("output", "model", db, glue(paste0("modelcox_glance_",selected_outcome,".csv"))))
  
  # combine and save model summary (full) 
  model_tidy <- bind_rows(summary0$tidy, summary0_full$tidy, 
                          summary1$tidy, summary1_full$tidy, 
                          summary2$tidy, summary2_full$tidy) %>%
    mutate(outcome = selected_outcome, 
           outcome_clean = selected_outcome_clean)
  write_csv(model_tidy, here::here("output", "model", db, glue(paste0("modelcox_tidy_full_",selected_outcome,".csv"))))
  
  # combine and save model summary (full - simplified)
  model_tidy_reduced <- data.frame(model_tidy) %>%
    filter(str_detect(term, fixed("timesincevax_pw")) | str_detect(term, fixed("vax2_az"))) %>%
    mutate(
      term=str_replace(term, pattern=fixed("vax2_az:strata(timesincevax_pw)"), ""),
      term=fct_inorder(term),
      term_left = as.numeric(str_extract(term, "^\\d+"))-1,
      term_right = as.numeric(str_extract(term, "\\d+$"))-1,
      term_right = if_else(is.na(term_right), lastfupday, term_right),
      term_midpoint = term_left + (term_right+1-term_left)/2
    )
  
  # add event counts from IRR table to unadjusted model
  model_tidy_reduced$BNT_n = irr_sub$BNT_n
  model_tidy_reduced$BNT_events = irr_sub$BNT_events
  model_tidy_reduced$AZ_n = irr_sub$AZ_n
  model_tidy_reduced$AZ_events = irr_sub$AZ_events
  redaction_columns = c("n_event", "exposure", "estimate", "std.error", "robust.se", "statistic", "p.value", "conf.low", "conf.high")
  for (i in 1:nrow(model_tidy_reduced)) {
    if (model_tidy_reduced$BNT_events[i]=="[Redacted]" | model_tidy_reduced$AZ_events[i]=="[Redacted]") { model_tidy_reduced[i,names(model_tidy_reduced)%in%redaction_columns] = "[Redacted]" }
    if (model_tidy_reduced$BNT_events[i]=="0" & model_tidy_reduced$AZ_events[i]=="0") { model_tidy_reduced[i,names(model_tidy_reduced)%in%redaction_columns] = "[No events]" }
  
  
  write_csv(model_tidy_reduced, here::here("output", "model", db, glue(paste0("modelcox_tidy_reduced_",selected_outcome,".csv"))))
  write_rds(model_tidy_reduced, here::here("output", "model", db, glue(paste0("modelcox_tidy_reduced_",selected_outcome,".rds"))), compress="gz")
}

