######################################

# This script:
# - imports processed data
# - fit univariate and multivariable stratified cox model(s) using the coxph package
# - saves models

######################################


### Preliminaries ----

## Import libraries
# library('here')
# library('tidyr')
library('tidyverse')
library('lubridate')
library('survival')
# library('gtsummary')
# library('gt')
# library('survminer')
library('glue')
# library('fs')
# library('splines')
sessionInfo()

## Import command-line arguments (specifying whether or not to run matched analysis)
args <- commandArgs(trailingOnly=TRUE)

## Set input and output pathways for matched/unmatched data - default is unmatched
if(length(args)==0){
  # default (unmatched) file names
  db = "main"
  input_name = "data_cohort_VE.rds"
  irr_name = "table_irr_redacted.rds"
  outcome = "covid_postest"
} else {
  outcome = args[[2]]
  if (args[[1]]=="unmatched") { 
    # unmatched file names    
    db = "main"
    input_name = "data_cohort_VE.rds"
    irr_name = "table_irr_redacted.rds"
  } else if (args[[1]]=="matched") {
    # matched file names
    db = "matched"
    input_name = "data_cohort_VE_matched.rds"
    irr_name = "table_irr_matched_redacted.rds"
  } else if (args[[1]]=="calendar_time") {
    # file names for calendar time models
    db = "calendartime"
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
### variables for models
#######################################################
surv_strata = c("tstart", "tstop", "ind_outcome") 
surv_full = c("follow_up_time", "ind_outcome")
expo = "vax2_az"
vars0_strata = "timesincevax_pw"
vars2 = c("age", "ckd_5cat", "immunosuppression", "care_home", "sex", "imd", 
         "ethnicity", "rural_urban_group", "prior_covid_cat", "prevax_tests_cat", 
         "multimorb", "sev_mental_ill")

####################################################### 
# prepare time to event data for the given outcome
####################################################### 

outcome_list = c("covid_postest", "covid_emergency", "covid_hosp", "covid_death")
clean_list = c("Positive SARS-CoV-2 test", "COVID-related A&E admission", "COVID-related hospitalisation", "COVID-related death")
date_list = c("postvax_positive_test_date", "postvax_covid_emergency_date", "postvax_covid_hospitalisation_date", "postvax_covid_death_date")

i = which(outcome_list == outcome)
selected_outcome = outcome_list[i]
selected_outcome_clean = clean_list[i]

# read in incidence rate ratio table
irr_table = read_rds(here::here("output", "tables", irr_name))
irr_sub = subset(irr_table, outcome_clean==selected_outcome_clean)[,c("period", "BNT_n", "BNT_events", "AZ_n", "AZ_events")]

# derive strata_var if using calendar timescale
if (db == "calendartime") {
  data_cohort <- data_cohort %>% mutate(strata_var = as.character(glue("{jcvi_group}_{region}")))
}

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
  ) %>%
  mutate(across(vars2[vars2 != "age"], as.factor))

data_cox_strata <- tmerge(
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
  ) %>%
  # set factor levels for postvaccination periods
  mutate(across(timesincevax_pw, factor, levels = postvax_periods)) %>%
  select(all_of(c(surv_strata, expo, vars0_strata, vars1, vars2))) 

data_cox_full <- data_tte %>%
  select(all_of(c(surv_full, expo, vars1, vars2))) 

rm(data_tte, data_cohort)

### check enough events to model ----
events_threshold <- 2

# data_cox_full
# check if >2 events for each level of expo in full data
check_events_full <- data_cox_full %>%
  group_by(across(all_of(expo))) %>%
  summarise(total_events = sum(ind_outcome), .groups = "keep") %>%
  ungroup() 

if (any(check_events_full$total_events <= events_threshold)) {
  # stop the script as too few events to model
  # save empty outputs
  write_rds(
    list(),
    here::here(glue("data_cox_empty_{db}_{outcome}.rds"))
  )
  write_rds(
    list(),
    here::here(glue("formulas_empty_{db}_{outcome}.rds"))
  )
  # wrap in try() to avoid the job from exiting with error code 1
  try(stop("Following condition not met: >2 events in each vaccine group."))
} 

# only keep postvaxcuts during which there are >2 events for each level of expo
data_cox_strata_keep <- data_cox_strata %>%
  group_by(across(all_of(c(vars0_strata, expo)))) %>%
  mutate(check_events_strata = sum(ind_outcome)) %>%
  ungroup() %>%
  group_by(across(all_of(vars0_strata))) %>%
  mutate(across(check_events_strata, min)) %>%
  ungroup() %>%
  filter(check_events_strata > events_threshold) %>%
  droplevels()

# strata controls whether to process data for stratified models
strata <- nrow(data_cox_strata_keep) > 0

if (db == "calendartime" & !strata) {
  try(stop("Not enough data for calendartime model."))
}

### merge or drop variable levels if <=2 events in either category of expo ----
merge_levels <- function(.data, var, strata=FALSE) {
  
  # if more than 2 levels, prepare levels to reorder factor according to event frequency
  old_levs <- levels(.data[[var]])
  if (length(old_levs) > 2) {
    new_levs <- .data %>%
      filter(ind_outcome) %>%
      group_by(across(all_of(var))) %>%
      count() %>%
      ungroup() %>%
      right_join(tibble(!! sym(var) := old_levs), by = var) %>%
      arrange(desc(n)) # NAs will be last (i.e. levels with no events)
    new_levs <- new_levs[[var]]
  } else {
    new_levs <- old_levs
  }
  
  # define group_vars
  if (strata) group_vars <- c(expo, vars0_strata) else group_vars <- expo
  
  data <- .data %>%
    # select only the necessary variables
    select(all_of(c(group_vars, var, "ind_outcome"))) %>%
    # reorder factor according to event frequency
    mutate(across(all_of(var), 
                  ~factor(as.character(.x), levels = new_levs)))
  
  event_count_fun <- function(var) {
    
    data %>%
      filter(ind_outcome) %>%
      group_by(across(all_of(c(group_vars, var)))) %>%
      count() %>%
      ungroup() %>%
      summarise(min_n=min(n)) %>%
      unlist() %>%
      unname()
    
  }
  
  var_levs <- levels(data[[var]])
  
  for (i in rev(seq_along(var_levs))) {
    
    if (event_count_fun(var) > events_threshold) {
      # no merging required if minimum count is above threshold
      # return original variable with original order of levels
      data <- data %>% 
        select(all_of(var)) %>%
        # reorder factor according to event frequency
        mutate(across(all_of(var), 
                      ~factor(as.character(.x), levels = new_levs)))
      return(data_out)
      break
    } else if (i > 2) {
      # merge labels i and i-1
      var_levs[(i-1)] <- str_c(var_levs[(i-1)], var_levs[i], sep = " / ")
      var_levs <- var_levs[-i]
      # merge levels in data  
      data <- data %>%
        mutate(across(var,
                      ~ factor(
                        if_else(
                          as.integer(.x) == i,
                          as.integer(i-1),
                          as.integer(.x)),
                        levels = 1:(i-1),
                        labels = var_levs)))
      
      
    } else {
      # if loop reaches i=2, break out of loop and return NULL
      return(NULL)
      break
    }
  }
  data %>% select(all_of(var))
}

if (db != "calendartime") {
  # merge levels in full dataset
  data_cox_full_merged <- data_cox_full %>%
    select(-all_of(vars2[vars2 != "age"])) %>%
    # join merged variables
    bind_cols(
      lapply(
        vars2[vars2 != "age"], 
        function(x) data_cox_full %>% merge_levels(var = x)
      )
    )
}

if (strata) {
  # merge levels in stratified dataset
  data_cox_strata_merged <- data_cox_strata_keep %>%
    select(-all_of(vars2[vars2 != "age"])) %>%
    # join merged variables
    bind_cols(
      lapply(
        vars2[vars2 != "age"], 
        function(x) data_cox_strata_keep %>% merge_levels(var = x, strata = TRUE)
      )
    )
}

### formulas ---
# I moved this here as some variables may have been dropped in the merging step
# formula1 same for stratified and full model, depends on db
if (db=="main") {
  formula1_update <- as.formula(". ~ . + strata(region)*ns(vax2_day, 3)")
} else if (db=="matched") {
  formula1_update <- as.formula(". ~ . + ns(vax2_day, 3)")
} else if (db=="calendartime") {
  formula1_update <- as.formula(". ~ . + strata(strata_var)")
}

## formulas for stratified model
if (strata) {
  
  # categorical vars for formula2
  vars2_formula_strata <- str_c(
    vars2[vars2 %in% names(data_cox_strata_merged) & vars2 != "age"], 
    collapse = " + "
  )
  if (vars2_formula_strata == "") {
    formula2_update_strata <- as.formula(glue(". ~ . + poly(age, degree = 2, raw = TRUE)"))
  } else {
    formula2_update_strata <- as.formula(glue(". ~ . + poly(age, degree = 2, raw = TRUE) + {vars2_formula_strata}"))
  }
  
  # final formulas
  formula0 <- as.formula(glue("Surv(tstart, tstop, ind_outcome) ~ {expo}:strata({vars0_strata})"))
  formula1 <- formula0 %>% update(formula1_update)
  formula2 <- formula1 %>% update(formula2_update_strata)
  
  formulas_strata <- list(formula0, formula1, formula2)
  
}

## formulas for full model
if (db != "calendartime") {
  
  # categorical vars for formula2
  vars2_formula_full <- str_c(
    vars2[vars2 %in% names(data_cox_full_merged) & vars2 != "age"], 
    collapse = " + "
  )
  
  if (vars2_formula_full == "") {
    formula2_update_full <- as.formula(glue(". ~ . + poly(age, degree = 2, raw = TRUE)"))
  } else {
    formula2_update_full <- as.formula(glue(". ~ . + poly(age, degree = 2, raw = TRUE) + {vars2_formula_full}"))
  }
  
  # final formulas
  formula0_full <- as.formula(glue("Surv(tstart, tstop, ind_outcome) ~ {expo}"))
  formula1_full <- formula0 %>% update(formula1_update)
  formula2_full <- formula1_full %>% update(formula2_update_full)
  
  formulas_full <- list(formula0, formula1, formula2)
  
}

print.f <- function(f) { 
  cat(paste(deparse(f, width.cutoff=getOption("width")), collapse="\n")) 
}


# print formulas to log file
if (db != "calendartime") {
  logoutput(
    "Formulas for the full model:\n",
    print.f(formulas_full[[1]]),
    print.f(formulas_full[[2]]),
    print.f(formulas_full[[3]])
  )
} else {
  logoutput(
    "No full model when using calendar timescale."
  )
}

if (strata) {
  logoutput(
    "Formulas for the full model:\n",
    print.f(formulas_strata[[1]]),
    print.f(formulas_strata[[2]]),
    print.f(formulas_strata[[3]])
  )
} else {
  logoutput(
    "Not enough events for the stratified model."
  )
}

### print dataset size and save outputs ----
logoutput(
  glue(selected_outcome_clean, "\ndata_cox data size = ", nrow(data_cox_strata_merged)),
  glue("data_cox memory usage = ", format(object.size(data_cox_strata_merged), units="GB", standard="SI", digits=3L))
)

# save data and formulas
# stratified
write_rds(
  data_cox_strata_merged,
  here::here(glue("data_cox_strata_{db}_{outcome}.rds")),
  compress = "gz"
)
write_rds(
  formulas_strata,
  here::here(glue("formulas_strata_{db}_{outcome}.rds"))
)

# full
if (db != "calendartime") {
  write_rds(
    data_cox_full_merged,
    here::here(glue("data_cox_full_{db}_{outcome}.rds")),
    compress = "gz"
  )
  write_rds(
    write_rds(
      formulas_full,
      here::here(glue("formulas_full_{db}_{outcome}.rds"))
    )
  )
}


