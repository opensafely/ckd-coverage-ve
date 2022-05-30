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
  outcome = "covid_postest"
} else {
  db = args[[1]]
  timescale = args[[2]]
  outcome = args[[3]]
}

## Import data
full = timescale != "calendartime"
if (full) {
  # get data for full model
  data_cox_full <- read_rds(here::here("output", "model", glue("data_cox_full_{db}_{timescale}_{outcome}.rds")))
}
# check dataset not empty
if (nrow(data_cox_full) == 0) {
  # log
  try(stop("Not enough events to fit model."))
}
# read data for stratified model
data_cox_strata <- read_rds(here::here("output", "model", glue("data_cox_strata_{db}_{timescale}_{outcome}.rds")))
strata = nrow(data_cox_strata) > 0
# check dataset not empty
if (!strata & timescale == "calendartime") {
  # log
  try(stop("Not enough events to fit model with calendar timescale."))
}

## Import formulas 
if (full) {
  formulas_full <- read_rds(here::here("output", "model", glue("formulas_full_{db}_{timescale}_{outcome}.rds")))
}
if (strata) {
  formulas_strata <- read_rds(here::here("output", "model", glue("formulas_strata_{db}_{timescale}_{outcome}.rds")))
}

# # read in incidence rate ratio table
# irr_name = glue("table_irr_{db}_{timescale}_redacted.rds")
irr_name = "table_irr_redacted.rds"
irr_table = read_rds(here::here("output", "tables", irr_name))
irr_sub = subset(irr_table, outcome_clean==selected_outcome_clean)[,c("period", "BNT_n", "BNT_events", "AZ_n", "AZ_events")]

## Import outcomes
outcomes_list <- read_rds(
  here::here("output", "lib", "outcomes.rds")
)
outcomes_index = which(outcomes_list$short_name == outcome)
selected_outcome = outcomes_list$short_name[outcomes_index]
selected_outcome_clean = outcomes_list$clean_name[outcomes_index]

## Import custom user functions and packages
# source(here::here("analysis", "functions.R"))

## Create directory for full model outputs
dir.create(here::here("output", "model"), showWarnings = FALSE, recursive=TRUE)

## create special log file ----
cat(glue("## script info for cox models ##"), "  \n", file = here::here("output", "model", glue("cox_model_log_{db}_{timescale}_{output}.txt")), append = FALSE)

## function to pass additional log text
logoutput <- function(...){
  cat(..., file = here::here("output", "model", glue("cox_model_log_{db}_{timescale}_{output}.txt")), sep = "\n  ", append = TRUE)
  cat("\n", file = here::here("output", "model", glue("cox_model_log_{db}_{timescale}_{output}.txt")), sep = "\n  ", append = TRUE)
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
  
  if (timescale == "calendartime" & !stratified) stop("Must set stratified=TRUE if using calendar timescale")
  
  if (stratified) {
    data_cox <- data_cox_strata
    # if stratified = FALSE, fit cox model stratified by follow-up window
    coxmod <- coxph(
      formula = formula_cox,
      data = data_cox,
      robust = TRUE, # compute robust variance
      id = patient_id, # required since multiple rows per subject
      na.action = "na.fail",
      control = coxph.control(iter.max = 50)
    )
  } else {
    data_cox <- data_cox_full
    # if stratified = FALSE, fit model on full dataset
    coxmod <- coxph(
      formula_cox, 
      data = data_cox,
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
    glance$convergence = NA # why?
  }
  
  # model outputs
  coxmod$data <- NULL
  #write_rds(coxmod, here("output", "model","VE", glue("modelcox_model{number}_",selected_outcome,".rds")), compress="gz")
  lst(glance, tidy)
}

#######################################################
# fit cox models
#######################################################
## Redact statistical outputs if <=10 events in glance
redact_glance <- function(.data) {
  redaction_columns = c(
    "nevent", "statistic.log", "p.value.log", "statistic.sc", "p.value.sc",
    "statistic.wald", "p.value.wald", "statistic.robust", "p.value.robust", 
    "r.squared", "r.squared.max", "concordance", "std.error.concordance", 
    "logLik", "AIC", "BIC")
  for (i in 1:nrow(.data)) {
    if (.data$nevent[i]>0 & .data$nevent[i]<=10) { .data[i,names(.data) %in% redaction_columns] = "[Redacted]" }
  }
}


if (full) {
  summary0_full <- cox_model_VE(0, formula0_full, stratified=FALSE)
  summary1_full <- cox_model_VE(1, formula1_full, stratified=FALSE)
  summary2_full <- cox_model_VE(2, formula2_full, stratified=FALSE)
  
  # glance
  model_glance_full <- bind_rows(
    summary0_full$glance, 
    summary1_full$glance, 
    summary2_full$glance
    ) %>%
    mutate(
      outcome = outcome, 
      outcome_clean = clean_list[outcome_list == outcome]
      ) %>%
    redact_glance()
  
  # save model summary
  write_csv(
    model_glance,
    here::here("output", "model", glue("modelcox_glance_full_{db}_{timescale}_{selected_outcome}.csv"))
    )
  
  
  # tidy
  model_tidy_full <- bind_rows(
    summary0_full$tidy, 
    summary1_full$tidy, 
    summary2_full$tidy
  ) %>%
    mutate(
      outcome = outcome, 
      outcome_clean = clean_list[outcome_list == outcome]
    )
  
  write_csv(
    model_tidy, 
    here::here("output", "model", glue("modelcox_tidy_full_{db}_{timescale}_{selected_outcome}.csv"))
    )
  
}
if (strata) {
  summary0_strata <- cox_model_VE(0, formula0_strata, stratified=TRUE)
  summary1_strata <- cox_model_VE(1, formula1_strata, stratified=TRUE)
  summary2_strata <- cox_model_VE(2, formula2_strata, stratified=TRUE)
  
  # glance
  model_glance_strata <- bind_rows(
    summary0_strata$glance, 
    summary1_strata$glance, 
    summary2_strata$glance
  ) %>%
    mutate(
      outcome = outcome, 
      outcome_clean = clean_list[outcome_list == outcome]
    ) %>%
    redact_glance()
  
  # save model summary
  write_csv(
    model_glance,
    here::here("output", "model", glue("modelcox_glance_strata_{db}_{timescale}_{selected_outcome}.csv"))
  )
  
  # tidy
  model_tidy_strata <- bind_rows(
    summary0_strata$tidy, 
    summary1_strata$tidy, 
    summary2_strata$tidy
  ) %>%
    mutate(
      outcome = outcome, 
      outcome_clean = clean_list[outcome_list == outcome]
    )
  
  write_csv(
    model_tidy, 
    here::here("output", "model", glue("modelcox_tidy_strata_{db}_{timescale}_{selected_outcome}.csv"))
  )
}

reduce_tidy <- function(.data) {
  
  # combine and save model summary (full - simplified)
  model_tidy_reduced <- .data %>%
    filter(
      str_detect(term, fixed("timesincevax_pw")) | str_detect(term, fixed("vax2_az"))
      ) %>%
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
  
  redaction_columns = c(
    "n_event", "exposure", "estimate", "std.error", "robust.se", "statistic", 
    "p.value", "conf.low", "conf.high")
  
  for (i in 1:nrow(model_tidy_reduced)) {
    if (
      model_tidy_reduced$BNT_events[i]=="[Redacted]" |
      model_tidy_reduced$AZ_events[i]=="[Redacted]"
    ) { 
      model_tidy_reduced[i,names(model_tidy_reduced) %in% redaction_columns] = "[Redacted]"
      }
    if (
      model_tidy_reduced$BNT_events[i]=="0" & 
      model_tidy_reduced$AZ_events[i]=="0"
    ) {
      model_tidy_reduced[i,names(model_tidy_reduced)%in%redaction_columns] = "[No events]" 
      }
  }
}

write_csv(model_tidy_reduced, here::here("output", "model", db, glue(paste0("modelcox_tidy_reduced_",selected_outcome,".csv"))))
write_rds(model_tidy_reduced, here::here("output", "model", db, glue(paste0("modelcox_tidy_reduced_",selected_outcome,".rds"))), compress="gz")
