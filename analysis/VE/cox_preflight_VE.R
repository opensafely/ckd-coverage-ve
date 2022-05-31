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
# arg1 = db = matched / unmatched
# arg2 = timescale = persontime / calendartime
# arg3 = outcome = covid_postest / covid_emergency / covid_hosp / covid_death
if(length(args)==0){
  # default (unmatched) file names
  db = "unmatched"
  timescale = "calendartime"
  selected_outcome = "covid_postest"
} else {
  db = args[[1]]
  timescale = args[[2]]
  selected_outcome = args[[3]]
}

if (db == "unmatched") {
  input_name = "data_cohort_VE.rds"
} else {
  input_name = glue("data_cohort_VE_{db}.rds")
}

if (db == "matched" & timescale == "calendartime") stop("Do not fit matched calendartime model.")

## Import data
data_cohort <- read_rds(here::here("output", "data", input_name)) %>%
  # to avoid timestamps causing inqualities for dates on the same day
  mutate(across(where(is.Date), 
                ~ floor_date(
                  as.Date(.x, format="%Y-%m-%d"),
                  unit = "days")))

## Import custom user functions and packages
source(here::here("analysis", "functions.R"))

## Create directory for full model outputs
dir.create(here::here("output", "model"), showWarnings = FALSE, recursive=TRUE)

## Import outcomes
#TODO
# probs want to do this in an upstream file, but I've put here for now
outcomes_list <- list(
  short_name = c("covid_postest", "covid_emergency", "covid_hosp", "covid_death"),
  clean_name = c("Positive SARS-CoV-2 test", "COVID-related A&E admission", "COVID-related hospitalisation", "COVID-related death"),
  date_name = c("postvax_positive_test_date", "postvax_covid_emergency_date", "postvax_covid_hospitalisation_date", "postvax_covid_death_date")
)
dir.create(here::here("output", "lib"), showWarnings = FALSE, recursive=TRUE)
write_rds(
  outcomes_list,
  here::here("output", "lib", "outcomes.rds")
)
###
outcomes_list <- read_rds(
  here::here("output", "lib", "outcomes.rds")
)

outcome_index = which(outcomes_list$short_name == selected_outcome)
selected_outcome_clean = outcomes_list$clean_name[outcome_index]


## Set analysis intervals and last follow-up day
# TODO
# specify this upstream?
period_length <- 56
n_periods <- 6
postvaxcuts <- period_length*0:(n_periods - 1)
postvax_periods = c("1-56", "57-112", "113-168", "169-224", "225-280")
lastfupday <- max(postvaxcuts)
write_rds(
  list(
    postvaxcuts = postvaxcuts,
    postvax_periods = postvax_periods
  ),
  here::here("output", "lib", "postvax_list.rds")
)

## create special log file ----
cat(
  glue("## script info for cox preflight ##"), 
  "  \n", 
  file = here::here("output", "model", glue("log_cox_preflight_{db}_{timescale}_{selected_outcome}.txt")), 
  append = FALSE
  )

## function to pass additional log text
logoutput <- function(...){
  cat(..., file = here::here("output", "model", glue("log_cox_preflight_{db}_{timescale}_{selected_outcome}.txt")), sep = "\n  ", append = TRUE)
  cat("\n", file = here::here("output", "model", glue("log_cox_preflight_{db}_{timescale}_{selected_outcome}.txt")), sep = "\n  ", append = TRUE)
}

logoutput(
  "args:", 
  glue("db = {db}"),
  glue("timescale = {timescale}"),
  glue("outcome = {selected_outcome_clean}")
)


### print dataset size ----
logoutput(
  "data_cohort:",
  glue("data size = ", nrow(data_cohort)),
  glue("memory usage = ", format(object.size(data_cohort), units="GB", standard="SI", digits=3L))
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
vars0 = "timesincevax_pw"
vars1 = c("region", "vax2_day", "strata_var")
vars2_cont = "age"
vars2_cat = c("ckd_5cat", "immunosuppression", "care_home", "sex", "imd", 
         "ethnicity", "rural_urban_group", "prior_covid_cat", "prevax_tests_cat", 
         "multimorb", "sev_mental_ill")

####################################################### 
# prepare time to event data for the given outcome
####################################################### 
# TODO

# derive strata_var if using calendar timescale
if (timescale == "calendartime") {
  
  data_cox_strata <- postvax_time %>%
    left_join(
      data_cohort %>%
        rename(
          outcome_date = outcomes_list$date_name[outcome_index]
        ) %>%
        select(patient_id, censor_date, vax2_date, outcome_date),
      by = "patient_id"
    ) %>%
    mutate(
      # origin date is the day before the earliest vax2 date in dataset
      min_vax_2_date = min(vax2_date) - days(1),
      # individual start date for a given period is the day before the vax2 date
      start_date = vax2_date + (id_postvax - 1)*period_length - days(1),
      # end date of a period is min of:
      end_date = pmin(start_date + period_length, censor_date, outcome_date, na.rm=TRUE)
      ) %>%
    # remove samples for which end_date<=start_date:
    filter(start_date < end_date) %>%
    # outcome indicator
    mutate(ind_outcome = !is.na(outcome_date) & outcome_date == end_date) %>%
    # replace outcome_date and censor_date with missing if they don't occur 
    # between start_date and end_date
    mutate(across(c(outcome_date, censor_date),
                  ~if_else(
                    !is.na(.x) & start_date < .x & .x <= end_date,
                    .x,
                    as.Date(NA_character_)
                  ))) %>%
    # calendar time-scale: time since min_vax_2_date
    mutate(across(c(start_date, end_date),
                  ~ as.integer(.x - min_vax_2_date))) %>%
    select(patient_id, id_postvax, timesincevax_pw, ind_outcome,
           tstart = start_date, tstop = end_date) %>%
    left_join(
      data_cohort %>%
        mutate(strata_var = as.character(glue("{jcvi_group}_{region}"))) %>%
        select(patient_id, all_of(c(expo, vars1, vars2_cont, vars2_cat))),
      by = "patient_id"
    )
  
  data_cox_full <- data_cox_strata
  
} else {
  
  data_tte <- data_cohort %>% 
    mutate(
      # select dates for outcome in question
      outcome_date = get(outcomes_list$date_name[outcome_index]),
      
      # censor date already defined in data_selection_VE.R script 
      
      # calculate tte and ind for outcome in question
      tte_outcome = tte(vax2_date-1, outcome_date, censor_date, na.censor=TRUE),
      ind_outcome = get(paste0("ind_",selected_outcome)),
      tte_stop = pmin(tte_censor, tte_outcome, na.rm=TRUE),
      
      # calculate follow-up time (censor/event)
      follow_up_time = tte(vax2_date-1, get(outcomes_list$date_name[outcome_index]), censor_date) 
    ) %>%
    mutate(across(all_of(vars2_cat), as.factor))
  
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
    # ind_outcome as logical to match data_tte
    mutate(across(ind_outcome, as.logical)) %>%
    select(
      patient_id,
      any_of(c(surv_strata, expo, vars0, vars1, vars2_cont, vars2_cat))
    ) %>%
    as_tibble()
  
  data_cox_full <- data_tte %>%
    select(
      patient_id,
      any_of(c(surv_full, expo, vars1, vars2_cont, vars2_cat))
    ) 
  
  rm(data_tte, data_cohort)
  
}


### check enough events to model ----
events_threshold = 2

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
    tibble(),
    here::here("output", "model", db, glue("data_cox_full_{db}_{timescale}_{selected_outcome}.rds"))
  )
  write_rds(
    list(),
    here::here("output", "model", db, glue("formulas_full_{db}_{timescale}_{selected_outcome}.rds"))
  )
  error_message = "Too few outcome events to model."
  logoutput(error_message)
  # wrap in try() to avoid the job from exiting with error code 1
  try(stop(error_message))
} 

# only keep postvaxcuts during which there are >2 events for each level of expo
data_cox_strata_keep <- data_cox_strata %>%
  group_by(across(all_of(c(vars0, expo)))) %>%
  mutate(check_events_strata = sum(ind_outcome)) %>%
  ungroup() %>%
  group_by(across(all_of(vars0))) %>%
  mutate(across(check_events_strata, min)) %>%
  ungroup() %>%
  filter(check_events_strata > events_threshold) %>%
  select(-check_events_strata) %>%
  droplevels() 

# strata controls whether to process data for stratified models
strata = nrow(data_cox_strata_keep) > 0

if (!strata) {
  
  error_message = "Not enough events for stratified model."
  logoutput(error_message)
  # save empty outputs
  write_rds(
    data_cox_strata_merged,
    here::here("output", "model", db, glue("data_cox_strata_{db}_{timescale}_{selected_outcome}.rds")),
    compress = "gz"
  )
  write_rds(
    formulas_strata,
    here::here("output", "model", db, glue("formulas_strata_{db}_{timescale}_{selected_outcome}.rds"))
  )
}

if (timescale == "calendartime" & !strata) {
  error_message = "Cannot fit calendar time model."
  logoutput(error_message)
  try(stop(error_message))
}

### merge or drop variable levels if <=2 events in either category of expo ----
merge_levels <- function(.data, var, strata=FALSE) {
  
  # if var not already factor, convert to factor
  .data <- .data %>% mutate(across(all_of(var), as.factor))
  
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
  if (strata) group_vars <- c(expo, vars0) else group_vars <- expo
  
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
                      ~factor(as.character(.x), levels = var_levs)))
      break
    } 
    if (event_count_fun(var) <= events_threshold & i == 2) {
      # if loop reaches i=2, break out of loop and return NULL
      data <- NULL
      break
    }
    if (event_count_fun(var) <= events_threshold & i > 2) {
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
      
      
    } 
  }
  if (!is.null(data)) data <- data %>% select(all_of(var))
  return(data)
}

if (timescale == "persontime") {
  # merge levels in full dataset
  data_cox_full_merged <- data_cox_full %>%
    select(-all_of(vars2_cat)) %>%
    # join merged variables
    bind_cols(
      lapply(
        vars2_cat, 
        function(x) data_cox_full %>% merge_levels(var = x)
      )
    )
}

if (strata) {
  # merge levels in stratified dataset
  data_cox_strata_merged <- data_cox_strata_keep %>%
    select(-all_of(vars2_cat)) %>%
    # join merged variables
    bind_cols(
      lapply(
        vars2_cat, 
        function(x) data_cox_strata_keep %>% merge_levels(var = x, strata = TRUE)
      )
    )
}

### formulas ---
# I moved this here as some variables may have been dropped in the merging step
# formula1 same for stratified and full model, depends on db and timescale
if (db == "unmatched" & timescale == "persontime") {
  formula1_update <- as.formula(". ~ . + strata(region)*ns(vax2_day, 3)")
} 
if (db == "matched" & timescale == "persontime") {
  formula1_update <- as.formula(". ~ . + ns(vax2_day, 3)")
} 
if (timescale == "calendartime") {
  formula1_update <- as.formula(". ~ . + strata(strata_var)")
}

## formulas for stratified model
if (strata) {
  
  # categorical vars for formula2
  vars2_formula_strata <- str_c(
    vars2_cat[vars2_cat %in% names(data_cox_strata_merged)], 
    collapse = " + "
  )
  if (vars2_formula_strata == "") {
    formula2_update_strata <- as.formula(glue(". ~ . + poly(age, degree = 2, raw = TRUE)"))
  } else {
    formula2_update_strata <- as.formula(glue(". ~ . + poly(age, degree = 2, raw = TRUE) + {vars2_formula_strata}"))
  }
  
  # final formulas
  formula0_strata <- as.formula(glue("Surv(tstart, tstop, ind_outcome) ~ {expo}:strata({vars0})"))
  formula1_strata <- formula0_strata %>% update(formula1_update)
  formula2_strata <- formula1_strata %>% update(formula2_update_strata)
  
  formulas_strata <- list(formula0_strata, formula1_strata, formula2_strata)
  
}

## formulas for full model
if (timescale == "persontime") {
  
  # categorical vars for formula2
  vars2_formula_full <- str_c(
    vars2_cat[vars2_cat %in% names(data_cox_full_merged)], 
    collapse = " + "
  )
  
  if (vars2_formula_full == "") {
    formula2_update_full <- as.formula(glue(". ~ . + poly(age, degree = 2, raw = TRUE)"))
  } else {
    formula2_update_full <- as.formula(glue(". ~ . + poly(age, degree = 2, raw = TRUE) + {vars2_formula_full}"))
  }
  
  # final formulas
  formula0_full <- as.formula(glue("Surv(follow_up_time, ind_outcome) ~ {expo}"))
  formula1_full <- formula0_full %>% update(formula1_update)
  formula2_full <- formula1_full %>% update(formula2_update_full)
  
  formulas_full <- list(formula0_full, formula1_full, formula2_full)
  
}

print.f <- function(f) { 
  paste(deparse(f, width.cutoff=getOption("width")), collapse="\n")
}


# print formulas to log file
logoutput("Formulas for the full model:")
if (timescale == "persontime") {
  logoutput(
    "\nformula0----\n",
    print.f(formulas_full[[1]]),
    "\nformula1----\n",
    print.f(formulas_full[[2]]),
    "\nformula2----\n",
    print.f(formulas_full[[3]]),
    "\n------\n"
  )
} else {
  logoutput(
    "No full model when using calendar timescale."
  )
}

logoutput("Formulas for the stratified model:")
if (strata) {
  logoutput(
    "\nformula0----\n",
    print.f(formulas_strata[[1]]),
    "\nformula1----\n",
    print.f(formulas_strata[[2]]),
    "\nformula2----\n",
    print.f(formulas_strata[[3]]),
    "\n------\n"
  )
} else {
  logoutput(
    "Not enough events for the stratified model."
  )
}

### print dataset size and save outputs ----
if (timescale == "persontime") {
  logoutput(
    glue("data_cox_full:"),
    glue("data size = ", nrow(data_cox_full_merged)),
    glue("memory usage = ", format(object.size(data_cox_full_merged), units="GB", standard="SI", digits=3L))
  )
}
if (strata) {
  logoutput(
    glue("data_cox_strata:"),
    glue("data size = ", nrow(data_cox_strata_merged)),
    glue("memory usage = ", format(object.size(data_cox_strata_merged), units="GB", standard="SI", digits=3L))
  )
}

# save data and formulas
# stratified
write_rds(
  data_cox_strata_merged,
  here::here("output", "model", glue("data_cox_strata_{db}_{timescale}_{selected_outcome}.rds")),
  compress = "gz"
)
write_rds(
  formulas_strata,
  here::here("output", "model", glue("formulas_strata_{db}_{timescale}_{selected_outcome}.rds"))
)

# full
if (timescale == "persontime") {
  write_rds(
    data_cox_full_merged,
    here::here("output", "model", glue("data_cox_full_{db}_{timescale}_{selected_outcome}.rds")),
    compress = "gz"
  )
  write_rds(
    formulas_full,
    here::here("output", "model", glue("formulas_full_{db}_{timescale}_{selected_outcome}.rds"))
  )
} else {
  write_rds(
    tibble(),
    here::here("output", "model", glue("data_cox_full_{db}_{timescale}_{selected_outcome}.rds"))
  )
  write_rds(
    list(),
    here::here("output", "model", glue("formulas_full_{db}_{timescale}_{selected_outcome}.rds"))
  )
}


