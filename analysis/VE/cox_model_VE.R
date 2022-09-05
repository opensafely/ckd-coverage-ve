######################################

# This script:
# - imports processed data
# - fits univariate and multivariable stratified cox model(s) using the coxph package
# - saves models

######################################

## Import libraries
library('tidyverse')
library('lubridate')
library('survival')
library('glue')
library('splines')
sessionInfo()

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
  timescale = "persontime"
  selected_outcome = "covid_postest"
  subgroup = "all"
  vaccine = "boost"
} else {
  db = args[[1]]
  timescale = args[[2]]
  selected_outcome = args[[3]]
  subgroup = args[[4]]
  vaccine = args[[5]]
}

## Import data for full model
data_cox_full <- read_rds(here::here("output", "model", paste0("VE_",vaccine), glue("data_cox_full_{db}_{timescale}_{selected_outcome}_{subgroup}.rds")))
full = nrow(data_cox_full) > 0

## Check dataset not empty
if (!full) {
  try(stop("Not enough events to fit model."))
}

## Read data for stratified model
data_cox_strata <- read_rds(here::here("output", "model", paste0("VE_",vaccine), glue("data_cox_strata_{db}_{timescale}_{selected_outcome}_{subgroup}.rds")))
strata = nrow(data_cox_strata) > 0

## Check dataset not empty
if (!strata & subgroup=="all") {
  try(stop("Not enough events to fit stratified model."))
}

## Import formulas 
if (full) {
  formulas_full <- read_rds(
    here::here("output", "model", paste0("VE_",vaccine), glue("formulas_full_{db}_{timescale}_{selected_outcome}_{subgroup}.rds")))
}

if (strata) {
  formulas_strata <- read_rds(
    here::here("output", "model", paste0("VE_",vaccine), glue("formulas_strata_{db}_{timescale}_{selected_outcome}_{subgroup}.rds")))
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

## Import postvax timepoints
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
lastfupday = max(postvax_list$postvaxcuts)

## Read in incidence rate ratio table
irr_table = read_rds(here::here("output", "tables", glue("table_irr_{vaccine}_redacted_{db}_{subgroup}.rds"))) %>%
  mutate_all(as.character)

## Specify irr strata and full levels
irr_sub = subset(irr_table, outcome_clean==selected_outcome_clean)[,c("period", "BNT_n", "BNT_events", "AZ_n", "AZ_events")]

if (vaccine=="primary") {
  irr_sub_strata <- irr_sub %>% filter(period %in% postvax_list$postvax_periods)
  irr_sub_full <- irr_sub %>% filter(!(period %in% postvax_list$postvax_periods))
} else if (vaccine=="boost") {
  irr_sub_strata <- irr_sub %>% filter(period %in% postvax_list$postvax_periods)
  irr_sub_full <- irr_sub %>% filter(!(period %in% postvax_list$postvax_periods))
} else {
  stop ("Arguments not specified correctly.")
}

## Create directory for full model outputs
dir.create(here::here("output", "model", paste0("VE_",vaccine)), showWarnings = FALSE, recursive=TRUE)

## Create log file
cat(
  glue("## Script info for cox models ##"), 
  "  \n", 
  file = here::here("output", "model", paste0("VE_",vaccine), glue("log_cox_model_{db}_{timescale}_{selected_outcome}_{subgroup}.txt")), 
  append = FALSE
  )

## Function to pass additional log text
logoutput <- function(...){
  cat(..., file = here::here("output", "model", paste0("VE_",vaccine), glue("log_cox_model_{db}_{timescale}_{selected_outcome}_{subgroup}.txt")), sep = "\n  ", append = TRUE)
  cat("\n", file = here::here("output", "model", paste0("VE_",vaccine), glue("log_cox_model_{db}_{timescale}_{selected_outcome}_{subgroup}.txt")), sep = "\n  ", append = TRUE)
}

logoutput(
  "args:", 
  glue("db = {db}"),
  glue("timescale = {timescale}"),
  glue("outcome = {selected_outcome_clean}"),
  glue("subset = {subgroup}"),
  glue("vaccine = {vaccine}")
)

####################################################### 
### Function to fit cox model for specified formula 
####################################################### 
# number = 0; stratified = TRUE

cox_model_VE <- function(number, stratified=TRUE) {
  if (number==0) model_type = "unadjusted"
  if (number==1) model_type = "region/date adjusted"
  if (number==2) model_type = "fully adjusted"
  formula_index = number + 1
  if (stratified) {
    formula_cox <- formulas_strata[[formula_index]]
  } else {
    formula_cox <- formulas_full[[formula_index]]
  }
  
  if (stratified) {
    ## If stratified = TRUE, fit cox model stratified by follow-up window
    coxmod <- coxph(
      formula = formula_cox,
      data = data_cox_strata,
      robust = TRUE, # compute robust variance
      id = patient_id, # required since multiple rows per subject
      na.action = "na.fail",
      control = coxph.control(iter.max = 50)
    )
  } else {
    ## If stratified = FALSE, fit model on full dataset
    coxmod <- coxph(
      formula = formula_cox,
      data = data_cox_full,
      control = coxph.control(iter.max = 50)
      )
  }
  ## Print warnings
  print(warnings())
  
  ## Print output status to log file
  logoutput(
    glue("model{number} data size = ", coxmod$n),
    glue("model{number} memory usage = ", format(object.size(coxmod), units="GB", standard="SI", digits=3L)),
    glue("model{number} iterations required = ", coxmod$iter)
  )

  ## Return model summary
  tidy <- broom.helpers::tidy_plus_plus(coxmod, exponentiate = FALSE) %>%
    add_column(model_name = model_type, .before=1) %>%
    add_column(model = number, .before=1)
  
  ## Return brief model summary
  glance <- broom::glance(coxmod) %>%
    add_column(
      model_name = model_type,
      model = number,
      ram = format(object.size(coxmod), units="GB", standard="SI", digits=3L),
      .before = 1
    )
  glance$convergence_iter = coxmod$iter
  
  ## Add whether model stratified or full
  if (stratified) { 
    tidy$level = glance$level = "stratified"
  } else {
    tidy$level = glance$level = "full"
  }

  ## Model outputs
  coxmod$data <- NULL
  #write_rds(coxmod, here("output", "model","VE", glue("modelcox_model{number}_",selected_outcome,".rds")), compress="gz")
  lst(glance, tidy)
}

#######################################################
## Fit cox models
#######################################################

## Redact statistical outputs if <=10 events in glance
redact_glance <- function(.data) {
  redaction_columns = c(
    "nevent", "statistic.log", "p.value.log", "statistic.sc", "p.value.sc",
    "statistic.wald", "p.value.wald", "statistic.robust", "p.value.robust", 
    "r.squared", "r.squared.max", "concordance", "std.error.concordance", 
    "logLik", "AIC", "BIC")
  .data = .data %>%
    ## All redaction columns as character, otherwise error in following loop
    mutate(across(any_of(redaction_columns), as.character))
  for (i in 1:nrow(.data)) {
    if (as.numeric(.data$nevent[i])>0 & as.numeric(.data$nevent[i])<=10) { .data[i,names(.data) %in% redaction_columns] = "[Redacted]" }
  }
  return(.data)
}

## Run full models
if (full) {
  if (db == "matched") {
    summary0_full <- cox_model_VE(0, stratified=FALSE)
    
    model_glance_full <- summary0_full$glance %>%
      mutate(
        outcome = selected_outcome, 
        outcome_clean = selected_outcome_clean
      ) %>%
      redact_glance()
    
    ## Collate tidy
    model_tidy_full <- summary0_full$tidy %>%
      mutate(
        outcome = selected_outcome, 
        outcome_clean = selected_outcome_clean
      )
    
  } else {
    summary0_full <- cox_model_VE(0, stratified=FALSE)
    summary1_full <- cox_model_VE(1, stratified=FALSE)
    summary2_full <- cox_model_VE(2, stratified=FALSE)
    
    ## Collate glance
    model_glance_full <- bind_rows(
      summary0_full$glance, 
      summary1_full$glance, 
      summary2_full$glance
    ) %>%
      mutate(
        outcome = selected_outcome, 
        outcome_clean = selected_outcome_clean
      ) %>%
      redact_glance()
    
    ## Collate tidy
    model_tidy_full <- bind_rows(
      summary0_full$tidy, 
      summary1_full$tidy, 
      summary2_full$tidy
    ) %>%
      mutate(
        outcome = selected_outcome, 
        outcome_clean = selected_outcome_clean
      )
  }
  
} else {
  ## Empty outputs
  model_glance_full <- tibble()
  model_tidy_full <- tibble()
}
## Save model summary
write_csv(
  model_glance_full,
  here::here("output", "model", paste0("VE_",vaccine), glue("modelcox_glance_full_{db}_{timescale}_{selected_outcome}_{subgroup}.csv"))
)
write_csv(
  model_tidy_full, 
  here::here("output", "model", paste0("VE_",vaccine), glue("modelcox_tidy_full_{db}_{timescale}_{selected_outcome}_{subgroup}.csv"))
)

## Run sratified models
if (strata) {
  if (db == "matched") {
    summary0_strata <- cox_model_VE(0, stratified=TRUE)

    ## Collate glance
    model_glance_strata <- summary0_strata$glance %>%
      mutate(
        outcome = selected_outcome, 
        outcome_clean = selected_outcome_clean
      ) %>%
      redact_glance()
    
    ## Collate tidy
    model_tidy_strata <- summary0_strata$tidy %>%
      mutate(
        outcome = selected_outcome, 
        outcome_clean = selected_outcome_clean
      )
    
  } else {
    summary0_strata <- cox_model_VE(0, stratified=TRUE)
    summary1_strata <- cox_model_VE(1, stratified=TRUE)
    summary2_strata <- cox_model_VE(2, stratified=TRUE)
    
    ## Collate glance
    model_glance_strata <- bind_rows(
      summary0_strata$glance, 
      summary1_strata$glance, 
      summary2_strata$glance
    ) %>%
      mutate(
        outcome = selected_outcome, 
        outcome_clean = selected_outcome_clean
      ) %>%
      redact_glance()
    
    ## Collate tidy
    model_tidy_strata <- bind_rows(
      summary0_strata$tidy, 
      summary1_strata$tidy, 
      summary2_strata$tidy
    ) %>%
      mutate(
        outcome = selected_outcome, 
        outcome_clean = selected_outcome_clean
      ) 
  }
  
} else {
  ## Empty outputs
  model_glance_strata <- tibble()
  model_tidy_strata <- tibble()
}

## Save model summary
write_csv(
  model_glance_strata,
  here::here("output", "model", paste0("VE_",vaccine), glue("modelcox_glance_strata_{db}_{timescale}_{selected_outcome}_{subgroup}.csv"))
)
write_csv(
  model_tidy_strata, 
  here::here("output", "model", paste0("VE_",vaccine), glue("modelcox_tidy_strata_{db}_{timescale}_{selected_outcome}_{subgroup}.csv"))
)

reduce_tidy <- function(stratified=FALSE) {
  
  if (stratified) {
    .data <- model_tidy_strata
    .irr <- irr_sub_strata
  } else {
    .data <- model_tidy_full 
    .irr <- irr_sub_full
  }
  
  model_tidy_reduced <- .data %>%
    filter(str_detect(term, "^vax2_az")) %>%
    mutate(across(term, 
                  ~if_else(
                    str_detect(.x, "timesincevax_pw"),
                    ## For the stratified model extract the time period from term:
                    str_extract(.x, "\\d+-\\d+$"), 
                    ## For the full model extract the time period from irr_sub_full:
                    irr_sub_full$period 
                    ))) 
  
  n_models <- n_distinct(model_tidy_reduced$model)
  n_periods <- n_distinct(.irr$period)
  model_tidy_reduced <- model_tidy_reduced %>%
    full_join(.irr %>% 
                 mutate(model = n_models) %>%
                 uncount(model) %>%
                 mutate(
                   model = rep(unique(model_tidy_reduced$model), n_periods),
                   model_name =  rep(unique(model_tidy_reduced$model_name),n_periods),
                 ), 
               by = c("term" = "period", "model", "model_name")) %>%
    mutate(
      # term=fct_inorder(term), # moved this to line 364, otherwise converted back to factor when bind_rows
      term_left = as.numeric(str_extract(term, "^\\d+"))-1,
      term_right = as.numeric(str_extract(term, "\\d+$"))-1,
      term_right = if_else(is.na(term_right), lastfupday, term_right),
      term_midpoint = term_left + (term_right+1-term_left)/2
    )
  
  redaction_columns = c(
    "n_event", "exposure", "estimate", "std.error", "robust.se", "statistic", 
    "p.value", "conf.low", "conf.high")
  
  ## Add event counts from IRR table to unadjusted model
  model_tidy_reduced <- model_tidy_reduced %>%
    select(
      model, model_name, term, variable, label, outcome,
      any_of(redaction_columns), starts_with(c("term_", "BNT_", "AZ_"))
    ) %>%
    ## All redaction columns as character, otherwise error in following loop
    mutate(across(any_of(redaction_columns), as.character))
  
  for (i in 1:nrow(model_tidy_reduced)) {
    if (
      model_tidy_reduced$BNT_events[i]=="[Redacted]" | model_tidy_reduced$AZ_events[i]=="[Redacted]"
    ) { 
      model_tidy_reduced[i,names(model_tidy_reduced) %in% redaction_columns] = "[Redacted]"
      }
    if (
      model_tidy_reduced$BNT_events[i]=="0" & model_tidy_reduced$AZ_events[i]=="0"
    ) {
      model_tidy_reduced[i,names(model_tidy_reduced) %in% redaction_columns] = "[No events]" 
      }
  }
  
  return(model_tidy_reduced)
  
}

## Combine and save model summary (full - simplified)
model_tidy_final <- tibble()
if (strata) {
  model_tidy_final <- bind_rows(model_tidy_final, reduce_tidy(TRUE))
} 
if (full) {
  model_tidy_final <- bind_rows(model_tidy_final, reduce_tidy(FALSE))
}
model_tidy_final <- model_tidy_final %>%
  mutate(across(term, fct_inorder)) %>%
  arrange(model, term)

## Save outputs
write_csv(
  model_tidy_final, 
  here::here("output", "model", paste0("VE_",vaccine), glue("modelcox_tidy_reduced_{db}_{timescale}_{selected_outcome}_{subgroup}.csv"))
  )
write_rds(
  model_tidy_final,
  here::here("output", "model", paste0("VE_",vaccine), glue("modelcox_tidy_reduced_{db}_{timescale}_{selected_outcome}_{subgroup}.rds")), 
  compress="gz"
  )
