######################################

# This script:
# - imports processed data
# - runs preflight checks to ensure at least 2 events per variable level per outcome period
# - merges levels if required
# - saves updated data and formulas

######################################

## Import libraries
library('tidyverse')
library('lubridate')
library('survival')
library('glue')
sessionInfo()

### Specify event threshold for preflight checks
events_threshold = 2

## Import command-line arguments (specifying whether or not to run matched analysis)
args <- commandArgs(trailingOnly=TRUE)

## Set input and output pathways for matched/unmatched data - default is unmatched
# arg1: db = matched / unmatched
# arg2: timescale = persontime / calendartime 
# arg3: outcome = covid_postest / covid_emergency / covid_hosp / covid_death / noncovid_death
# arg4: subset = all / CKD3 / CKD4-5 / RRT
# arg5: vaccine = primary / boost

if(length(args)==0){
  # default (unmatched) file names
  db = "unmatched"
  timescale = "calendartime"
  selected_outcome = "covid_postest"
  subgroup = "all"
  vaccine = "primary"
} else {
  db = args[[1]]
  timescale = args[[2]]
  selected_outcome = args[[3]]
  subgroup = args[[4]]
  vaccine = args[[5]]
}

## Specify input data
if (db == "unmatched" & vaccine == "primary") { 
  input_name = "data_cohort_VE.rds"
  selected_vax_date = "vax2_date"
  selected_vax_day = "vax2_day"
} else if (db == "matched" & vaccine == "primary") { 
  input_name = "data_cohort_VE_matched.rds"
  selected_vax_date = "vax2_date"
  selected_vax_day = "vax2_day"
} else if (db == "unmatched" & vaccine == "boost") { 
  input_name = "data_cohort_VE_boost.rds"
  selected_vax_date = "vax3_date"
  selected_vax_day = "vax3_day"
} else if (db == "matched" & vaccine == "boost") { 
  input_name = "data_cohort_VE_boost_matched.rds"
  selected_vax_date = "vax3_date"
  selected_vax_day = "vax3_day"
} else {
  stop ("Arguments not specified correctly.")
}

## Stop if argument combination does not match analysis plan
if (db == "matched" & timescale == "calendartime") stop ("Do not fit matched calendar time model.")
if (subgroup != "all" & timescale == "calendartime") stop ("Do not fit calendar time model to subgroups.")

## Import data and process dates
data_cohort <- read_rds(here::here("output", "data", input_name)) %>%
  ## Define vaccination date for start of follow-up
  mutate(vax_date = get(selected_vax_date), 
         vax_day = get(selected_vax_day)
  )

## Select subset
if (subgroup=="all") {
  data_cohort = data_cohort
} else if (subgroup=="CKD3") {
  data_cohort = subset(data_cohort, ckd_3cat == "CKD3")
} else if (subgroup=="CKD4-5") {
  data_cohort = subset(data_cohort, ckd_3cat == "CKD4-5")
} else if (subgroup=="RRT") {
  data_cohort = subset(data_cohort, ckd_3cat == "RRT (any)")
} else {
  stop ("Arguments not specified correctly.")
}

## Import custom user functions and packages
source(here::here("analysis", "functions.R"))

## Create directory for full model outputs
dir.create(here::here("output", "model", paste0("VE_",vaccine)), showWarnings = FALSE, recursive=TRUE)

## Import outcome time periods
if (vaccine=="primary") {
  postvax_list <- read_rds(
    here::here("output", "lib", "postvax_list.rds")
  )
} else if (vaccine=="boost") {
  postvax_list <- read_rds(
    here::here("output", "lib", "postboost_list.rds")
  )
} else {
  stop ("Arguments not specified correctly.")
}

## Import outcomes
if (vaccine=="primary") {
  outcomes_list <- read_rds(
    here::here("output", "lib", "outcomes.rds")
  )
} else if (vaccine=="boost") {
  outcomes_list <- read_rds(
    here::here("output", "lib", "outcomes_boost.rds")
  )
} else {
  stop ("Arguments not specified correctly.")
}
outcome_index = which(outcomes_list$short_name == selected_outcome)
selected_outcome_clean = outcomes_list$clean_name[outcome_index]
selected_outcome_date_name = outcomes_list$date_name[outcome_index]

## Create log file
cat(
  glue("## Script info for cox preflight ##"), 
  "  \n", 
  file = here::here("output", "model", paste0("VE_",vaccine), glue("log_cox_preflight_{db}_{timescale}_{selected_outcome}_{subgroup}.txt")), 
  append = FALSE
)

## Function to pass additional log text
logoutput <- function(...){
  cat(..., file = here::here("output", "model", paste0("VE_",vaccine), glue("log_cox_preflight_{db}_{timescale}_{selected_outcome}_{subgroup}.txt")), sep = "\n  ", append = TRUE)
  cat("\n", file = here::here("output", "model", paste0("VE_",vaccine), glue("log_cox_preflight_{db}_{timescale}_{selected_outcome}_{subgroup}.txt")), sep = "\n  ", append = TRUE)
}

## Pass parameters to log file
logoutput(
  "args:", 
  glue("db = {db}"),
  glue("timescale = {timescale}"),
  glue("outcome = {selected_outcome_clean}"),
  glue("subset = {subgroup}"),
  glue("vaccine = {vaccine}")
)

## Print dataset size
logoutput(
  "data_cohort:",
  glue("data size = ", nrow(data_cohort)),
  glue("memory usage = ", format(object.size(data_cohort), units="GB", standard="SI", digits=3L))
)

## Create dataset containing one row per patient per post-vaccination period
postvaxcuts = postvax_list$postvaxcuts
postvax_periods = postvax_list$postvax_periods
period_length = postvaxcuts[2]-postvaxcuts[1]
lastfupday = max(postvaxcuts)

postvax_time <- data_cohort %>%
  select(patient_id) %>%
  uncount(weights = length(postvaxcuts), .id="id_postvax") %>%
  mutate(
    fup_day = postvaxcuts[id_postvax],
    timesincevax_pw = timesince_cut(fup_day+1, postvaxcuts)
  ) %>%
  droplevels()
# 4 rows per patient (one per follow-up period)
# id_postvax = comparison period
# fup_day = final day of preceding follow-up period
# timesincevax_pw = days included in follow-up period
# select(patient_id, fup_day, timesincevax_pw)



####################################################### 
### Specify variables for models
#######################################################
surv = c("tstart", "tstop", "ind_outcome") 
expo = "vax2_az"
vars0 = "timesincevax_pw"
vars1 = c("vax_day", "jcvi_region")
vars2_cont = "age"
vars2_cat = c("sex", "imd", "ethnicity", "rural_urban_group", "ckd_3cat", "multimorb",
              "learning_disability", "sev_mental_ill", "any_immunosuppression", "prior_covid_cat", "prevax_tests_cat")

## Drop ckd_3cat category from subgroup analyses
if (subgroup!="all") {
  vars2_cat <- vars2_cat[!(vars2_cat %in% c("ckd_3cat"))]
}

## Matched models to be run without additional confounder adjustment

####################################################### 
### Prepare time to event data for the given outcome
####################################################### 

## Derive time to event data
data_tte <- data_cohort %>%
  transmute(
    
    patient_id,
    
    ## Select dates for outcome in question
    outcome_date = get(outcomes_list$date_name[outcome_index]),
    ind_outcome = get(paste0("ind_",selected_outcome)),
    
    ## Select vax date
    vax_date,
    
    ## Censor date already defined in data_selection_VE.R script
    censor_date,
    tte_censor, 
    
    ## Calculate tte and ind for outcome in question
    tte_outcome = tte(
      origin_date = vax_date-1, 
      event_date = outcome_date, 
      censor_date = censor_date,
      na.censor = TRUE
    )
  )

## Define full data
data_cox_full <- data_tte %>%
  mutate(
    tstart = 0,
    tstop = tte_outcome
  )

if (timescale == "calendartime") {
  
  ## Function for rescaling variables to calendar timescale
  rescale_calendartime <- function(.data) {
    .data %>%
      # time since the earliest vax_date in the dataset
      mutate(rescale = as.integer(vax_date - min(vax_date)),
             unscaled_tte_censor = tte_censor,
             unscaled_tte_outcome = tte_outcome) %>%
      mutate(across(starts_with(c("tte", "tstart", "tstop")),
                    ~ .x + rescale)) %>%
      select(-rescale)
  }
  
  ## Apply rescale_calendartime function
  data_cox_full <- data_cox_full %>% rescale_calendartime()
  
}

## Function for adding the relevant covariates to the tte data
add_covars <- function(.data) {
  .data %>%
    left_join(
      data_cohort %>%
        mutate(jcvi_region = as.character(glue("{jcvi_group}_{region}"))) %>%
        select(patient_id, all_of(c(expo, vars1, vars2_cont, vars2_cat))),
      by = "patient_id"
    ) %>%
    mutate(across(all_of(vars2_cat), as.factor))
}

data_cox_full <- data_cox_full %>% add_covars()

## Create stratified dataset if running model for primary schedule in whole population
if (subgroup=="all" & vaccine=="primary") strata <- TRUE else strata <- FALSE

if (strata) {
  data_cox_strata <- tmerge(
    data1 = data_tte %>% select(-ind_outcome, -censor_date, -outcome_date),
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
    ## Set factor levels for postvaccination periods
    mutate(across(timesincevax_pw, factor, levels = postvax_periods)) %>%
    ## ind_outcome as logical to match data_tte
    mutate(across(ind_outcome, as.logical)) %>%
    as_tibble()
  
  if (timescale == "calendartime") {
    data_cox_strata <- data_cox_strata %>% rescale_calendartime()
  }
  
  data_cox_strata <- data_cox_strata %>% add_covars()
}

## Clean-up
rm(data_tte, data_cohort)

################################################################################
## Check if >2 events for each level of expo in full data
check_events_full <- data_cox_full %>%
  group_by(across(all_of(expo))) %>%
  summarise(total_events = sum(ind_outcome), .groups = "keep") %>%
  ungroup() 
full = all(check_events_full$total_events > events_threshold) 

if (!full) {
  ## Save empty outputs if not enough events to model
  write_rds(
    tibble(),
    here::here("output", "model", paste0("VE_",vaccine), glue("data_cox_full_{db}_{timescale}_{selected_outcome}_{subgroup}.rds"))
  )
  write_rds(
    list(),
    here::here("output", "model", paste0("VE_",vaccine), glue("formulas_full_{db}_{timescale}_{selected_outcome}_{subgroup}.rds"))
  )
  ## Stop the script if too few events to model
  error_message = "Too few outcome events to model."
  logoutput(error_message)
  ## Wrap in try() to avoid the job from exiting with error code 1
  try(stop(error_message))
}

## Specify stratified dataset
if (strata) {
  ## Only keep postvaxcuts during which there are >2 events for each level of expo
  data_cox_strata_keep <- data_cox_strata %>%
    filter(!is.na(timesincevax_pw) & timesincevax_pw!="183+") %>% ## Added to cut periods of 1-14d and 183d+ (primary)
    group_by(across(all_of(c(vars0, expo)))) %>%
    mutate(check_events_strata = sum(ind_outcome)) %>%
    ungroup() %>%
    group_by(across(all_of(vars0))) %>%
    mutate(across(check_events_strata, min)) %>%
    ungroup() %>%
    filter(check_events_strata > events_threshold) %>%
    select(-check_events_strata) %>%
    droplevels() 
  
  ## Strata controls whether to process data for stratified models
  strata = nrow(data_cox_strata_keep) > 0
} 

if (!strata) {
  if (vaccine=="primary" & subgroup=="all") {
    error_message = "Not enough events for stratified model."
  } else {
    error_message = "Stratified model excluded."
  }
  logoutput(error_message)
  ## Save empty outputs
  write_rds(
    tibble(),
    here::here("output", "model", paste0("VE_",vaccine), glue("data_cox_strata_{db}_{timescale}_{selected_outcome}_{subgroup}.rds")),
    compress = "gz"
  )
  write_rds(
    tibble(),
    here::here("output", "model", paste0("VE_",vaccine), glue("formulas_strata_{db}_{timescale}_{selected_outcome}_{subgroup}.rds"))
  )
}

####################################################### 
### Pre-flight functions
####################################################### 

## Function to merge or drop variable levels if <=2 events in either category of expo 
merge_levels <- function(.data, var, strata=FALSE) {
  
  ## If var not already factor, convert to factor
  .data <- .data %>% mutate(across(all_of(var), as.factor))
  
  ## If more than 2 levels, prepare levels to reorder factor according to event frequency
  old_levs <- levels(.data[[var]])
  ## Which variables do we want to reorder by frequency of events?
  if (var %in% c("ethnicity", "rural_urban_group")) {
    ## if (length(old_levs) > 2) {
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
  
  ## Define group_vars
  if (strata) group_vars <- c(expo, vars0) else group_vars <- expo
  
  data <- .data %>%
    ## Select only the necessary variables
    select(all_of(c(group_vars, var, "ind_outcome"))) %>%
    ## Reorder factor according to event frequency
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
  ## Loop in reverse order, so last two levels first to merge
  for (i in rev(seq_along(var_levs))) {
    
    if (event_count_fun(var) > events_threshold) {
      ## No merging required if minimum count is above threshold
      ## Return original variable with original order of levels
      data <- data %>% 
        select(all_of(var)) %>%
        ## Reorder factor according to event frequency
        mutate(across(all_of(var), 
                      ~factor(as.character(.x), levels = var_levs)))
      break
    } 
    if (event_count_fun(var) <= events_threshold & i == 2) {
      ## If loop reaches i=2, break out of loop and return NULL
      data <- NULL
      break
    }
    if (event_count_fun(var) <= events_threshold & i > 2) {
      ## Merge labels i and i-1
      var_levs[(i-1)] <- str_c(var_levs[(i-1)], var_levs[i], sep = " / ")
      var_levs <- var_levs[-i]
      ## Merge levels in data  
      data <- data %>%
        mutate(across(all_of(var),
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

## Function to summaries outputs of function above
merge_summary <- function(data_in) {
  
  name_data <- deparse(substitute(data_in))
  
  vars2_keep <- vars2_cat[vars2_cat %in% names(data_in)]
  vars2_drop <- str_c(vars2_cat[!(vars2_cat %in% vars2_keep)], collapse = ", ")
  if (length(vars2_drop) == 0) vars2_drop <- "NA"
  out_drop <- glue("Dropped variables ----\n{vars2_drop}")
  out_keep <- "Kept variables ----\nNA"
  if (length(vars2_keep) > 0) {
    out_keep <- "Kept variables with levels (merged levels separated by \"/\") ---- "
    for (i in vars2_keep) {
      out_keep <- c(out_keep, str_c(glue("levels({i}): "), str_c(levels(data_in[[i]]), collapse = "; ")))
    }
  }
  out <- str_c(
    c(glue("Merge summary for \"{name_data}\" ----"), " ", out_drop, " ", out_keep), 
    collapse = "\n"
  )
  return(out)
}


####################################################### 
### Run pre-flight functions
####################################################### 

## Run pre-flight merging for unmatched data
if (db == "unmatched") {
  ### Full unmatched person-time models (all and subgroup)
  # Merge levels in full dataset
  data_cox_full_merged <- data_cox_full %>%
    select(-all_of(vars2_cat)) %>%
    # Join merged variables
    bind_cols(
      lapply(
        vars2_cat, 
        function(x) data_cox_full %>% merge_levels(var = x)
      )
    )
  logoutput(merge_summary(data_cox_full_merged))
  
  ## Stratified unmatched models, excluding subgroup analyses or booster analyses (where strata = FALSE)
  if (strata) {
    ## Merge levels in stratified dataset
    data_cox_strata_merged <- data_cox_strata_keep %>%
      select(-all_of(vars2_cat)) %>%
      ## Join merged variables
      bind_cols(
        lapply(
          vars2_cat, 
          function(x) data_cox_strata_keep %>% merge_levels(var = x, strata = TRUE)
        )
      )
    logoutput(merge_summary(data_cox_strata_merged))
  }
  
  ## Set formula updates for calendar time vs person time models
  if (timescale == "persontime") {
    formula1_update <- as.formula(". ~ . + strata(jcvi_region)*ns(vax_day, 3)")
  } 
  if (timescale == "calendartime") {
    formula1_update <- as.formula(". ~ . + strata(jcvi_region)")
  }
}

## Specify stratified matched models
if (strata) {
  
  formula0_strata <- as.formula(glue("Surv(tstart, tstop, ind_outcome) ~ {expo}:strata({vars0})"))
  
  if (db == "matched") {
    
    formulas_strata <- list(formula0_strata)
    
    ## Save formulas
    write_rds(
      formulas_strata,
      here::here("output", "model", paste0("VE_",vaccine), glue("formulas_strata_{db}_{timescale}_{selected_outcome}_{subgroup}.rds"))
    )
    ## Save data
    write_rds(
      data_cox_strata_keep, # Full stratified dataset rather than merged dataset
      here::here("output", "model", paste0("VE_",vaccine), glue("data_cox_strata_{db}_{timescale}_{selected_outcome}_{subgroup}.rds")),
      compress = "gz"
    )
    
  ## Specify stratified unmatched models
  } else if (db == "unmatched") {
    
    ## Categorical vars for formula2
    vars2_formula_strata <- str_c(
      vars2_cat[vars2_cat %in% names(data_cox_strata_merged)], 
      collapse = " + "
    )
    if (vars2_formula_strata == "") {
      formula2_update_strata <- as.formula(glue(". ~ . + poly(age, degree = 2, raw = TRUE)"))
    } else {
      formula2_update_strata <- as.formula(glue(". ~ . + poly(age, degree = 2, raw = TRUE) + {vars2_formula_strata}"))
    }
    
    ## Final formulas
    formula1_strata <- formula0_strata %>% update(formula1_update)
    formula2_strata <- formula1_strata %>% update(formula2_update_strata)
    formulas_strata <- list(formula0_strata, formula1_strata, formula2_strata)
    
    ## Save formulas
    write_rds(
      formulas_strata,
      here::here("output", "model", paste0("VE_",vaccine), glue("formulas_strata_{db}_{timescale}_{selected_outcome}_{subgroup}.rds"))
    )
    ## Save data
    write_rds(
      data_cox_strata_merged,
      here::here("output", "model", paste0("VE_",vaccine), glue("data_cox_strata_{db}_{timescale}_{selected_outcome}_{subgroup}.rds")),
      compress = "gz"
    )
  }
}

## Specify full matched model
formula0_full <- as.formula(glue("Surv(tstart, tstop, ind_outcome) ~ {expo}"))

if (db == "matched") {
  
  formulas_full <- list(formula0_full)
  
  ## Save formulas
  write_rds(
    formulas_full,
    here::here("output", "model", paste0("VE_",vaccine), glue("formulas_full_{db}_{timescale}_{selected_outcome}_{subgroup}.rds"))
  )
  ## Save data
  write_rds(
    data_cox_full, # Full dataset rather than merged dataset
    here::here("output", "model", paste0("VE_",vaccine), glue("data_cox_full_{db}_{timescale}_{selected_outcome}_{subgroup}.rds")),
    compress = "gz"
  )
  
}

## Specify full unmatched model
if (db == "unmatched") {
  
  ## Categorical vars for formula2
  vars2_formula_full <- str_c(
    vars2_cat[vars2_cat %in% names(data_cox_full_merged)], 
    collapse = " + "
  )
  
  if (vars2_formula_full == "") {
    formula2_update_full <- as.formula(glue(". ~ . + poly(age, degree = 2, raw = TRUE)"))
  } else {
    formula2_update_full <- as.formula(glue(". ~ . + poly(age, degree = 2, raw = TRUE) + {vars2_formula_full}"))
  }
  
  ## Final formulas
  formula1_full <- formula0_full %>% update(formula1_update)
  formula2_full <- formula1_full %>% update(formula2_update_full)
  formulas_full <- list(formula0_full, formula1_full, formula2_full)
  
  ## Save formulas
  write_rds(
    formulas_full,
    here::here("output", "model", paste0("VE_",vaccine), glue("formulas_full_{db}_{timescale}_{selected_outcome}_{subgroup}.rds"))
  )
  ## Save data
  write_rds(
    data_cox_full_merged,
    here::here("output", "model", paste0("VE_",vaccine), glue("data_cox_full_{db}_{timescale}_{selected_outcome}_{subgroup}.rds")),
    compress = "gz"
  )
}

## Function for printing functions
print.f <- function(f) { 
  paste(deparse(f, width.cutoff=getOption("width")), collapse="\n")
}

## Print formulas to log file
logoutput("Formulas for the full model:")

if (db == "matched") {
  logoutput(
    "\nformula0----\n",
    print.f(formulas_full[[1]])
  )
} else if (db == "unmatched") {
  logoutput(
    "\nformula0----\n",
    print.f(formulas_full[[1]]),
    "\nformula1----\n",
    print.f(formulas_full[[2]]),
    "\nformula2----\n",
    print.f(formulas_full[[3]]),
    "\n------\n"
  )
}

if (strata) {
  
  logoutput("Formulas for the stratified model:")
  
  if (db == "matched") {
    logoutput(
      "\nformula0----\n",
      print.f(formulas_strata[[1]]),
      "\n------\n"
    )
  }
  
  if (db == "unmatched") {
    logoutput(
      "\nformula0----\n",
      print.f(formulas_strata[[1]]),
      "\nformula1----\n",
      print.f(formulas_strata[[2]]),
      "\nformula2----\n",
      print.f(formulas_strata[[3]]),
      "\n------\n"
    )
    
  }
}

## Print dataset size and save outputs
if (db == "unmatched" & timescale == "persontime") {
  logoutput(
    glue("data_cox_full:"),
    glue("data size = ", nrow(data_cox_full_merged)),
    glue("memory usage = ", format(object.size(data_cox_full_merged), units="GB", standard="SI", digits=3L))
  )
}
if (db == "matched" & timescale == "persontime") {
  logoutput(
    glue("data_cox_full:"),
    glue("data size = ", nrow(data_cox_full)),
    glue("memory usage = ", format(object.size(data_cox_full), units="GB", standard="SI", digits=3L))
  )
}
if (db == "unmatched" & strata) {
  logoutput(
    glue("data_cox_strata:"),
    glue("data size = ", nrow(data_cox_strata_merged)),
    glue("memory usage = ", format(object.size(data_cox_strata_merged), units="GB", standard="SI", digits=3L))
  )
}
if (db == "matched" & strata) {
  logoutput(
    glue("data_cox_strata:"),
    glue("data size = ", nrow(data_cox_strata_keep)),
    glue("memory usage = ", format(object.size(data_cox_strata_keep), units="GB", standard="SI", digits=3L))
  )
}
