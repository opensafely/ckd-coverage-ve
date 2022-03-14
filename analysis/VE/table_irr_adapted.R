# # # # # # # # # # # # # # # # # # # # #
# This script creates a table of incidence rates for diffferent outcomes, stratified by vaccine type and times since vaccination.
# # # # # # # # # # # # # # # # # # # # #

# Preliminaries ----

## Import libraries ----
library('tidyverse')
library('here')
library('glue')
library('survival')
library('gt')
library('gtsummary')
library('scales')
library('lubridate')

## Import custom user functions from lib

# source(here::here("analysis", "lib", "utility_functions.R"))
# source(here::here("analysis", "lib", "redaction_functions.R"))
# 

## Import custom user functions
source(here::here("analysis", "functions.R"))

# args <- commandArgs(trailingOnly=TRUE)
# if(length(args)==0){
#   # use for interactive testing
#   removeobs <- FALSE
# } else {
#   removeobs <- TRUE
# }
# 
# 
# 
# ## import metadata ----
# var_labels <- read_rds(here("output", "data", "metadata_labels.rds"))
# 
# list_formula <- read_rds(here("output", "data", "metadata_formulas.rds"))
# list2env(list_formula, globalenv())
postvaxcuts_8week <- 56*0:6
lastfupday <- max(postvaxcuts_8week)

#metadata_outcomes <- read_rds(here("output", "data", "metadata_outcomes.rds"))

## create output directory ----
# fs::dir_create(here("output", "descriptive", "tables"))

## Import processed data ----
data_cohort <- read_rds(here("output", "data", "data_cohort_VE.rds"))
#data_comp <- read_rds(here("output", "data", "data_cohort_comp.rds"))


# vaccine initiation dates
data_cohort$end_date = as_date("2021-11-30")

# create pt data ----

data_tte <- data_cohort %>%
  transmute(
    patient_id,
    vax2_type,
    postvax_positive_test_date,
    postvax_covid_emergency_date,
    postvax_covid_hospitalisation_date,
    postvax_covid_death_date,
    end_date,
    
    censor_date = pmin(vax2_date - 1 + lastfupday, end_date, vax3_date, dereg_date, death_date, na.rm=TRUE),
    
    # time to last follow up day
    # tte function format: origin date, event date, censor date
    tte_enddate = tte(vax2_date-1, end_date, end_date),
    
    # time to last follow up day or death or deregistration
    tte_censor = tte(vax2_date-1, censor_date, censor_date),
    
    # time to positive test
    tte_covid_postest = tte(vax2_date-1, postvax_positive_test_date, censor_date, na.censor=TRUE),
    ind_covid_postest = censor_indicator(postvax_positive_test_date, censor_date),
    
    # time to COVID-19 A&E attendance
    tte_covid_emergency = tte(vax2_date-1, postvax_covid_emergency_date, censor_date, na.censor=TRUE),
    ind_covid_emergency = censor_indicator(postvax_covid_emergency_date, censor_date),
    
    # time to COVID-19 hospitalisation
    tte_covid_hosp = tte(vax2_date-1, postvax_covid_hospitalisation_date, censor_date, na.censor=TRUE),
    ind_covid_hosp = censor_indicator(postvax_covid_hospitalisation_date, censor_date),
    
    # time to COVID-19 death
    tte_covid_death = tte(vax2_date-1, postvax_covid_death_date, censor_date, na.censor=TRUE),
    ind_covid_death = censor_indicator(postvax_covid_death_date, censor_date)
  ) %>%
   filter(
     # TDOD remove once study def rerun with new dereg date
     tte_censor>0 | is.na(tte_censor)
   )

# one row per patient per post-vaccination week
postvax_time <- data_tte %>%
  select(patient_id, tte_censor) %>%
  mutate(
    fup_day = list(postvaxcuts_8week),
    timesincevax = map(fup_day, ~droplevels(timesince_cut_end(.x+1, postvaxcuts_8week)))
  ) %>%
  unnest(c(fup_day, timesincevax))

# create dataset that splits follow-up time by
# time since vaccination (using postvaxcuts cutoffs)
data_cox_split <- tmerge(
  data1 = data_tte %>% select(-starts_with("ind_"), -ends_with("_date")),
  data2 = data_tte,
  id = patient_id,
  tstart = 0L,
  tstop = tte_censor,
  covid_postest = event(tte_covid_postest),
  covid_emergency = event(tte_covid_emergency),
  covid_hosp = event(tte_covid_hosp),
  covid_death = event(tte_covid_death),

  status_covid_postest = tdc(tte_covid_postest),
  status_covid_emergency = tdc(tte_covid_emergency),
  status_covid_hosp = tdc(tte_covid_hosp),
  status_covid_death = tdc(tte_covid_death)
) %>%
  tmerge( # create treatment timescale variables
    data1 = .,
    data2 = postvax_time,
    id = patient_id,
    timesincevax = tdc(fup_day, timesincevax)
  ) %>%
  mutate(
    pt = tstop - tstart
  )

## create person-time table ----

format_ratio = function(numer,denom, width=7){
  paste0(
    replace_na(scales::comma_format(accuracy=1)(numer), "--"),
    " /",
    str_pad(replace_na(scales::comma_format(accuracy=1)(denom),"--"), width=width, pad=" ")
  )
}



pt_summary <- function(data, event){

  #data=data_cox_split
  #event = "test"
  unredacted <- data %>%
    mutate(
      timesincevax,
      status_event = .[[paste0("status_", event)]],
      ind_event = .[[event]],
      event = event
    ) %>%
    group_by(vax2_type, event, timesincevax) %>%
    summarise(
      yearsatrisk=sum(pt*(1-status_event))/365.25,
      n=sum(ind_event),
      rate=n/yearsatrisk
    ) %>%
    ungroup()

  unredacted_all <- data %>%
    mutate(
      status_event = .[[paste0("status_", event)]],
      ind_event = .[[event]],
      event = event
    ) %>%
    group_by(vax2_type, event) %>%
    summarise(
      yearsatrisk=sum(pt*(1-status_event))/365.25,
      n=sum(ind_event),
      rate=n/yearsatrisk
    ) %>%
    ungroup()

  unredacted_add_all <-
    bind_rows(
      unredacted,
      unredacted_all
    ) %>%
    mutate(
      timesincevax=forcats::fct_explicit_na(timesincevax, na_level="All")
    )


  unredacted_wide <-
    unredacted_add_all %>%
    pivot_wider(
      id_cols =c(event, timesincevax),
      names_from = vax2_type,
      values_from = c(yearsatrisk, n, rate),
      names_glue = "{vax2_type}_{.value}"
    ) %>%
    select(
      event, timesincevax, starts_with("pfizer"), starts_with("az")
    ) %>%
    mutate(
      rr = az_rate / pfizer_rate,
      rrE = scales::label_number(accuracy=0.01, trim=FALSE)(rr),
      rrCI = rrCI_exact(az_n, az_yearsatrisk, pfizer_n, pfizer_yearsatrisk, 0.01)
    )

  redacted <- unredacted_wide %>%
    mutate(
      pfizer_rate = redactor2(pfizer_n, 5, pfizer_rate),
      az_rate = redactor2(az_n, 5, az_rate),

      rr = redactor2(pmin(az_n, pfizer_n), 5, rr),
      rrE = redactor2(pmin(az_n, pfizer_n), 5, rrE),
      rrCI = redactor2(pmin(az_n, pfizer_n), 5, rrCI),

      pfizer_n = redactor2(pfizer_n, 5),
      az_n = redactor2(az_n, 5),

      pfizer_q = format_ratio(pfizer_n, pfizer_yearsatrisk),
      az_q = format_ratio(az_n, az_yearsatrisk)
    )
}

data_summary <- bind_rows(
  temp1 <- pt_summary(data_cox_split, "covid_postest"),
  temp2 <- pt_summary(data_cox_split, "covid_emergency"),
  temp3 <- pt_summary(data_cox_split, "covid_hosp"),
  temp4 <- pt_summary(data_cox_split, "covid_death")
)

write_csv(data_summary, here("output", "tables", "table_irr.csv"))

# tab_summary <- data_summary %>%
#   select(-ends_with("_n"), -ends_with("_yearsatrisk"), -rrE) %>%
#   gt() %>%
#   cols_label(
#     event = "Outcome",
#     timesincevax = "Time since first dose",
# 
#     pfizer_q = "Events / person-years",
#     az_q   = "Events / person-years",
# 
#     pfizer_rate = "Incidence",
#     az_rate = "Incidence",
#     rr = "Incidence rate ratio",
# 
#     rrCI = "95% CI"
#   ) %>%
#   tab_spanner(
#     label = "BNT162b2",
#     columns = starts_with("pfizer")
#   ) %>%
#   tab_spanner(
#     label = "ChAdOx1",
#     columns = starts_with("az")
#   ) %>%
#   fmt_number(
#     columns = ends_with(c("rr", "_rate")),
#     decimals = 2
#   ) %>%
#   fmt_missing(
#     everything(),
#     missing_text="--"
#   ) %>%
#   cols_align(
#     align = "right",
#     columns = everything()
#   ) %>%
#   cols_align(
#     align = "left",
#     columns = "timesincevax"
#   )
# 
# write_rds(tab_summary, here("output", "descriptive", "tables", "table_irr.rds"))
# gtsave(tab_summary, here("output", "descriptive", "tables", "table_irr.html"))


tab_summary_simple <-
  data_summary %>%
  filter(
    event %in% c(
      "covid_postest",
      "covid_emergency",
      "covid_hosp",
      "covid_death"
    )
  ) %>%
  transmute(
    event, timesincevax,
    pfizer_q,
    pfizer_rate_fmt = if_else(pfizer_rate<0.001, "<0.001", number_format(0.001)(pfizer_rate)),
    pfizer_rate_fmt = if_else(pfizer_rate==0, "0", pfizer_rate_fmt),
    az_q,
    az_rate_fmt = if_else(az_rate<0.001, "<0.001", number_format(0.001)(az_rate)),
    az_rate_fmt = if_else(az_rate==0, "0", az_rate_fmt)
  ) %>%
  gt(
    groupname_col = "event"
  ) %>%
  cols_label(
    event = "Outcome",
    timesincevax = "Time since first dose",

    pfizer_q = "BNT162b2\nEvents / person-years",
    pfizer_rate_fmt = "BNT162b2\nIncidence",

    az_q   = "ChAdOx1\nEvents / person-years",
    az_rate_fmt = "ChAdOx1\nIncidence"
  ) %>%
  fmt_missing(
    everything(),
    missing_text="--"
  ) %>%
  cols_align(
    align = "right",
    columns = everything()
  ) %>%
  cols_align(
    align = "left",
    columns = "timesincevax"
  )

# write_rds(tab_summary_simple, here("output", "tables", "table_irr_simple.rds"))
gtsave(tab_summary_simple, here("output", "tables", "table_irr_simple.html"))

