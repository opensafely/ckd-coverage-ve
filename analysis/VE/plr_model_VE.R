######################################
# This script:
# imports processed data and restricts it to patients in "cohort"
# fits pooled logistic regression models, with different adjustment sets
######################################

### Preliminaries ----

## Import libraries
library('tidyverse')
library('here')
library('glue')
library('survival')
library('splines')
library('parglm')
library('sandwich')
library('lubridate')
sessionInfo()

## Specify number of non-outcomes to sample (5000 if running locally, 50,000 if on server)
if(Sys.getenv("OPENSAFELY_BACKEND") %in% c("", "expectations")){
  samplesize_nonoutcomes_n <- 5000
} else {
  samplesize_nonoutcomes_n <- 50000
}

## Import command-line arguments (specifying whether or not to run matched analysis)
args <- commandArgs(trailingOnly=TRUE)

## Set input and output pathways for matched/unmatched data - default is unmatched
if(length(args)==0){
  # default (unmatched) file names
  db = "VE_plr"
  input_name = "data_cohort_VE.rds"
  irr_name = "table_irr_redacted.rds"
} else {
  if (args[[1]]=="unmatched") { 
    # unmatched file names    
    db = "VE_plr"
    input_name = "data_cohort_VE.rds"
    irr_name = "table_irr_redacted.rds"
  } else if (args[[1]]=="matched") {
    # matched file names
    db = "VE_plr_matched"
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

## Function from HCW comparative effectiveness repo for retrieving time-dependent hazard ratio ----
get_HR <- function(.data, model, vcov, term_index){
  ## function to get AZ/BNT hazard ratio spline over time since first dose
  tstop <- .data$tstop
  
  tt <- terms(model) # this helpfully grabs the correct spline basis from the model, rather than recalculating based on `.data`
  Terms <- delete.response(tt)
  mat0 <- model.matrix(Terms, data=mutate(.data, vax2_az=0))
  mat1 <- model.matrix(Terms, data=mutate(.data, vax2_az=1))
  
  tstop_distinct <- unique(tstop)
  
  # take a single row for each follow up time, and only select columns relevant for time-dependent vaccine effect
  partialX0 <- mat0[match(tstop_distinct, tstop), term_index]
  partialX1 <- mat1[match(tstop_distinct, tstop), term_index]
  
  # calculate difference between vaccine types on linear scale
  diffX <- partialX1-partialX0
  
  # get vaccine relevant  estimates
  partialbeta <- coef(model)[term_index]
  partialV <- vcov[term_index, term_index]
  
  tibble(
    tstop = tstop_distinct,
    loghr = unname((diffX %*% partialbeta)[,1]),
    se = sqrt( rowSums(diffX * (diffX %*% partialV)) ),
    loghr.ul = loghr + se * qnorm(0.025),
    loghr.ll = loghr - se * qnorm(0.025)
  )
}

## Create directory for full model outputs
dir.create(here::here("output", "model", db), showWarnings = FALSE, recursive=TRUE)

## Set analysis intervals and last follow-up day
postvaxcuts <- 56*0:5
postvax_periods = c("1-56", "57-112", "113-168", "169-224", "225-280")
lastfupday <- max(postvaxcuts)

## create special log file ----
cat(glue("## script info for plr models ##"), "  \n", file = here::here("output", "model", db, glue("modelplr_log.txt")), append = FALSE)

## function to pass additional log text
logoutput <- function(...){
  cat(..., file = here::here("output", "model", db, glue("modelplr_log.txt")), sep = "\n  ", append = TRUE)
  cat("\n", file = here::here("output", "model", db, glue("modelplr_log.txt")), sep = "\n  ", append = TRUE)
}

### print dataset size ----
logoutput(
  glue("data_cohort data size = ", nrow(data_cohort)),
  glue("data_cohort memory usage = ", format(object.size(data_cohort), units="GB", standard="SI", digits=3L))
)


####################################################### 
### formulae for unadjusted/adjusted models
####################################################### 
# plr models stratified by follow-up window
formula0 <- outcome_event ~ vax2_az + timesincevax_pw + vax2_az:timesincevax_pw
formula1 <- formula0 %>% update(. ~ . + region*ns(vax2_day, 3))
formula2 <- formula1 %>% update(. ~ . + poly(age, degree = 2, raw = TRUE) + ckd_7cat + immunosuppression + care_home + sex + imd + ethnicity + rural_urban_group + prior_covid_cat + prevax_tests_cat + multimorb + sev_mental_ill)

## optimisation options ----
parglmparams <- parglm.control(
  method = "LINPACK",
  nthreads = 8,
  maxit = 40 # default = 25
)


####################################################### 
### function to fit cox model for specified formula 
####################################################### 

#formula_plr = formula0
#number = 0

plr_summary_VE <- function(plrmod, number, outcome) {
  if (number==0) { model_type = "unadjusted" } 
  if (number==1) { model_type = "region/date adjusted" } 
  if (number==2) { model_type = "fully adjusted" } 
  
  # print warnings
  print(warnings())
  
  # print output status to log file
  logoutput(
    glue("model{number} data size = ", plrmod$n),
    glue("model{number} memory usage = ", format(object.size(plrmod), units="GB", standard="SI", digits=3L)),
    glue("convergence status: ", plrmod$converged)
  )
  
  # if (stratified) {
  #   logoutput(
  #     glue("convergence status: ", coxmod$info[["convergence"]])
  #   )
  # }
  
  # return model summary
  tidy <- broom.helpers::tidy_plus_plus(
    plrmod, 
    exponentiate = FALSE,
    tidy_fun = tidy_plr,
    cluster = data_plr$patient_id) %>%
    add_column(model_name = model_type, .before=1) %>%
    add_column(model = number, .before=1)
  
  # return brief model summary
  glance <- glance_plr(plrmod) %>%
    add_column(
      model_name = model_type,
      model = number,
      convergence = plrmod$converged,
      ram = format(object.size(plrmod), units="GB", standard="SI", digits=3L),
      .before = 1
    )

  # if (stratified) {
  #   tidy$level = glance$level = "stratified"
  #   glance$convergence = coxmod$info[["convergence"]]
  # } else {
  #   tidy$level = glance$level = "full"
  #   glance$convergence = NA
  # }
  
  vcov <- vcovCL(plrmod, cluster = data_plr$patient_id, type = "HC0")
  write_rds(vcov, here("output", "model", db, glue("modelplr_vcov_{outcome}_{number}.rds")), compress="gz")
  
  plrmod$data <- NULL
  write_rds(plrmod, here("output", "model", db, glue("modelplr_model_{outcome}_{number}.rds")), compress="gz")
  
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

#for (i in 1:length(outcome_list)) {
for (i in 1) {
  selected_outcome = outcome_list[i]
  selected_outcome_clean = clean_list[i]
  irr_sub = subset(irr_table, outcome_clean==selected_outcome_clean)[,c("period", "BNT_n", "BNT_events", "AZ_n", "AZ_events")]
  
  ## Calculate time to event variables
  data_tte <- data_cohort %>%
    mutate(
      # select dates for outcome in question
      outcome_date = get(date_list[i]),
      
      # censor date already defined in data_selection_VE.R script 
      
      # calculate tte and ind for outcome in question
      tte_outcome = tte(vax2_date-1, outcome_date, censor_date, na.censor=TRUE),
      ind_outcome = get(paste0("ind_",selected_outcome)),
      tte_stop = pmin(tte_censor, tte_outcome, na.rm=TRUE),
      
      sample_outcome = sample_nonoutcomes_n(!is.na(tte_outcome), patient_id, samplesize_nonoutcomes_n),
      sample_weights = sample_weights(!is.na(tte_outcome), sample_outcome)
    ) %>%
    filter(
      # select all patients who experienced the outcome, and a sample of those who don't
      sample_outcome==1L
    )

  alltimes <- expand(data_tte, patient_id, times=as.integer(full_seq(c(1, tte_stop),1)))
  
  # one row per day
  data_plr <-
    tmerge(
      data1 = data_tte %>% select(everything()),
      data2 = data_tte,
      id = patient_id,
      
      outcome_status = tdc(tte_outcome),
      censor_status= tdc(tte_censor),
      
      outcome_event = event(tte_outcome),
      censor_event = event(tte_censor),
      
      tstart = 0L,
      tstop = tte_stop
      
    ) %>%
    tmerge(
      data2 = alltimes,
      id = patient_id,
      alltimes = event(times, times)
    ) %>%
    mutate(
      timesincevax_pw = droplevels(timesince_cut(tstop, postvaxcuts)),
      tstart_calendar = tstart + vax2_day - 1,
      tstop_calendar = tstop + vax2_day - 1
    )
  
  ### print dataset size and save ----
  logoutput(
    glue(selected_outcome_clean,"\ndata_pt data size = ", nrow(data_plr)),
    glue("data_pt memory usage = ", format(object.size(data_plr), units="GB", standard="SI", digits=3L))
  )
  
  ## Fit models with outcome grouped into 8-week windows
  ## Model fitting and processing not included in same function since sandwich::vcovCL unable to handle formulae 
  ## Model 0: vaccination + timescale only, unadjusted
  plrmod0 <- parglm(
    formula = formula0,
    data = data_plr,
    weights = sample_weights,
    family = binomial,
    control = parglmparams,
    na.action = "na.fail",
    model = FALSE
  )
  summary0 <- plr_summary_VE(plrmod0, 0, selected_outcome)

  ## Model 1: adjusting for region/date of vaccination
  plrmod1 <- parglm(
    formula = formula1,
    data = data_plr,
    weights = sample_weights,
    family = binomial,
    control = parglmparams,
    na.action = "na.fail",
    model = FALSE
  )
  summary1 <- plr_summary_VE(plrmod1, 1, selected_outcome)
  
  ## Model 2: fully adjusted model including demography and comorbidities
  plrmod2 <- parglm(
    formula = formula2,
    data = data_plr,
    weights = sample_weights,
    family = binomial,
    control = parglmparams,
    na.action = "na.fail",
    model = FALSE
  )
  summary2 <- plr_summary_VE(plrmod2, 2, selected_outcome)

  ## Combine and save model summary (brief) 
  model_glance <- data.frame(
    bind_rows(summary0$glance,  
              summary1$glance,
              summary2$glance) %>%
      mutate(outcome = selected_outcome, 
             outcome_clean = selected_outcome_clean)
  )
  
  ## Redact statistical outputs if <=5 events
  # redaction_columns = c("nevent", "statistic.log", "p.value.log", "statistic.sc", "p.value.sc", "statistic.wald", "p.value.wald", 
  #                       "statistic.robust", "p.value.robust", "r.squared", "r.squared.max", "concordance", "std.error.concordance", "logLik", "AIC", "BIC")
  # for (i in 1:nrow(model_glance)) {
  #   if (model_glance$nevent[i]>0 & model_glance$nevent[i]<=5) { model_glance[i,names(model_glance)%in%redaction_columns] = "[Redacted]" }
  # }
  write_csv(model_glance, here::here("output", "model", db, glue(paste0("modelplr_glance_",selected_outcome,".csv"))))

  ## Combine and save model summary (full) 
  model_tidy <- bind_rows(summary0$tidy,
                          summary1$tidy,
                          summary2$tidy) %>%
    mutate(outcome = selected_outcome, 
           outcome_clean = selected_outcome_clean)
  write_csv(model_tidy, here::here("output", "model", db, glue(paste0("modelplr_tidy_full_",selected_outcome,".csv"))))

  ## Import plrmod outputs for selected outcome
  plrmod0 <- read_rds(here("output", "model", db, glue("modelplr_model_{selected_outcome}_0.rds")))
  plrmod1 <- read_rds(here("output", "model", db, glue("modelplr_model_{selected_outcome}_1.rds")))
  plrmod2 <- read_rds(here("output", "model", db, glue("modelplr_model_{selected_outcome}_2.rds")))

  ## Import vcov outputs for selected outcome
  vcov0 <- read_rds(here("output", "model", db, glue("modelplr_vcov_{selected_outcome}_0.rds")))
  vcov1 <- read_rds(here("output", "model", db, glue("modelplr_vcov_{selected_outcome}_1.rds")))
  vcov2 <- read_rds(here("output", "model", db, glue("modelplr_vcov_{selected_outcome}_2.rds")))

  ## define term indices
  term_index0 <- str_detect(names(coef(plrmod0)), fixed("timesincevax")) | str_detect(names(coef(plrmod0)), fixed("vax2_az"))
  term_index1 <- str_detect(names(coef(plrmod1)), fixed("timesincevax")) | str_detect(names(coef(plrmod1)), fixed("vax2_az"))
  term_index2 <- str_detect(names(coef(plrmod2)), fixed("timesincevax")) | str_detect(names(coef(plrmod2)), fixed("vax2_az"))
  
  ## extract HR estimates
  effectsplr <-
    bind_rows(
      get_HR(data_plr, plrmod0, vcov0, term_index=term_index0) %>% mutate(model=0, model_type="unadjusted"),
      get_HR(data_plr, plrmod1, vcov1, term_index=term_index1) %>% mutate(model=1, model_type="region/date adjusted"),
      get_HR(data_plr, plrmod2, vcov2, term_index=term_index2) %>% mutate(model=2, model_type="fully adjusted")
    ) %>%
    select(model, everything()) %>%
    group_by(model) %>%
    mutate(
      lag_tstop=lag(tstop, 1, 0),
      hr=exp(loghr),
      hr.ll=exp(loghr.ll),
      hr.ul=exp(loghr.ul)
    ) %>%
    ungroup() %>%
    left_join(
      model_tidy %>% group_by(model_name, model) %>% summarise() %>% ungroup(), by="model"
    )

  ### extract estimates at midpoints of cuts
  cuts <-
    tibble(
      left = postvaxcuts[-length(postvaxcuts)],
      right = postvaxcuts[-1],
      midpoint = left + ((right - left ) / 2),
      midpoint_rounded = as.integer(midpoint)
    )

  effectsplr_pw <- effectsplr %>%
    filter(tstop %in% cuts$midpoint_rounded) %>%
    left_join(cuts, by=c("tstop"="midpoint_rounded"))

  write_rds(effectsplr_pw, path = here("output", "model", db, glue("modelplr_VE_estimates_{selected_outcome}.rds")))
  write_csv(effectsplr_pw, path = here("output", "model", db, glue("modelplr_VE_estimates_{selected_outcome}.csv")))
}
