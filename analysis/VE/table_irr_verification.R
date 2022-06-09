# # # # # # # # # # # # # # # # # # # # #
# This script creates a table of incidence rates for different outcomes, stratified by vaccine type and times since vaccination.
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

## Import custom user functions
source(here::here("analysis", "functions.R"))

## Specify follow-up windows
postvaxcuts_8week <- c(28*0:6)
lastfupday <- max(postvaxcuts_8week)

## Create output directory
fs::dir_create(here("output", "descriptive", "tables"))

## Import processed data ----
data_cohort <- read_rds(here("output", "data", "data_cohort_VE.rds"))

## Create pt data
data_tte <- data_cohort %>%
  transmute(
    patient_id,
    vax2_type,
    postvax_positive_test_date,
    postvax_covid_emergency_date,
    postvax_covid_hospitalisation_date,
    postvax_covid_death_date,
    end_date,
    censor_date,
    tte_censor,
    tte_covid_postest, 
    ind_covid_postest, 
    tte_covid_emergency, 
    ind_covid_emergency, 
    tte_covid_hosp, 
    ind_covid_hosp,
    tte_covid_death,
    ind_covid_death
    ) 

## One row per patient per post-vaccination week
postvax_time <- data_tte %>%
  select(patient_id, tte_censor) %>%
  mutate(
    fup_day = list(postvaxcuts_8week),
    timesincevax = map(fup_day, ~droplevels(timesince_cut_end(.x+1, postvaxcuts_8week)))
  ) %>%
  unnest(c(fup_day, timesincevax))

## Create dataset that splits follow-up time by time since vaccination (using postvaxcuts cutoffs)
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
  tmerge( ## Create treatment timescale variables
    data1 = .,
    data2 = postvax_time,
    id = patient_id,
    timesincevax = tdc(fup_day, timesincevax)
  ) %>%
  mutate(
    pt = tstop - tstart
  )

## Create person-time table 
format_ratio = function(numer,denom, width=7){
  paste0(
    replace_na(scales::comma_format(accuracy=1)(numer), "--"),
    " /",
    str_pad(replace_na(scales::comma_format(accuracy=1)(denom),"--"), width=width, pad=" ")
  )
}

## Create function to summarise IRRs
pt_summary <- function(data, event){

  unredacted <- data_cox_split %>%
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

## Collate IRRs by outcome
data_summary <- bind_rows(
  temp1 <- pt_summary(data_cox_split, "covid_postest"),
  temp2 <- pt_summary(data_cox_split, "covid_emergency"),
  temp3 <- pt_summary(data_cox_split, "covid_hosp"),
  temp4 <- pt_summary(data_cox_split, "covid_death")
)

## Write csv
write_csv(data_summary, here("output", "tables", "table_irr_redacted_verification.csv"))