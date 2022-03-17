
# # # # # # # # # # # # # # # # # # # # #
# This script:
# imports processed data
# filters out people who are excluded from the main analysis
# outputs inclusion/exclusions flowchart data
# # # # # # # # # # # # # # # # # # # # #

# Preliminaries ----

## Import libraries ----
library('tidyverse')
library('lubridate')
library('here')
library('glue')

## Import custom user functions
source(here::here("analysis", "functions.R"))

## Create output directory
fs::dir_create(here::here("output", "tables"))

## Import processed data ----
data_processed <- read_rds(here::here("output", "data", "data_processed.rds"))

## vaccine initiation dates
first_pfizer = as_date("2020-12-08")
first_az = as_date("2021-01-04")
first_moderna = as_date("2021-04-13")

## Set analysis end date
data_processed$end_date = as_date("2021-11-30")

## Set analysis intervals
postvaxcuts <- 56*0:5
lastfupday <- max(postvaxcuts)

## Create cohort data with tte calculations ----
data_processed <- data_processed %>%
  mutate(
    # set censor date
    censor_date = pmin(vax2_date - 1 + lastfupday, end_date, dereg_date, death_date, vax3_date, na.rm=TRUE),
    
    # time to last follow up day
    # tte function format: origin date, event date, censor date
    tte_enddate = tte(origin_date = vax2_date - 1, # day of vaccination = day 1 of tte
                      event_date = end_date, 
                      censor_date = end_date),
    
    # time to last follow up day or death or deregistration
    tte_censor = tte(vax2_date - 1, censor_date, censor_date),
    
    # time to positive test
    tte_covid_postest = tte(vax2_date - 1, postvax_positive_test_date, censor_date, na.censor=TRUE),
    ind_covid_postest = dplyr::if_else((postvax_positive_test_date>censor_date) | is.na(postvax_positive_test_date), FALSE, TRUE),
    
    # time to COVID-19 A&E attendance
    tte_covid_emergency = tte(vax2_date - 1, postvax_covid_emergency_date, censor_date, na.censor=TRUE),
    ind_covid_emergency = dplyr::if_else((postvax_covid_emergency_date>censor_date) | is.na(postvax_covid_emergency_date), FALSE, TRUE),
    
    # time to COVID-19 hospitalisation
    tte_covid_hosp = tte(vax2_date - 1, postvax_covid_hospitalisation_date, censor_date, na.censor=TRUE),
    ind_covid_hosp = dplyr::if_else((postvax_covid_hospitalisation_date>censor_date) | is.na(postvax_covid_hospitalisation_date), FALSE, TRUE),
    
    # time to COVID-19 death
    tte_covid_death = tte(vax2_date - 1, postvax_covid_death_date, censor_date, na.censor=TRUE),
    ind_covid_death = dplyr::if_else((postvax_covid_death_date>censor_date) | is.na(postvax_covid_death_date), FALSE, TRUE)
  )

# Define selection criteria ----
data_criteria <- data_processed %>%
  transmute(
    patient_id,
    
    # Age + CKD
    has_age = !is.na(age) & age >=16 & age<120,
    has_ckd_any = ckd_inclusion_any==1,
    has_ckd_strict = ckd_inclusion_strict==1,
    
    # Demography
    has_sex = !is.na(sex),
    has_imd = !is.na(imd),
    has_ethnicity = !is.na(ethnicity),
    has_region = !is.na(region),
    
    # Vaccine profile
    vax_pfi_az = !is.na(vax12_type) & (vax12_type=="az-az" | vax12_type=="pfizer-pfizer"),
    vax_date = !is.na(vax1_date) & vax1_date>first_az,
    vax_interval = !is.na(tbv1_2) & tbv1_2>=(8*7) & tbv1_2<(16*7),
    
    # Postvax events
    positive_test_date_check = is.na(postvax_positive_test_date) | postvax_positive_test_date>=vax2_date,
    emergency_date_check = is.na(postvax_covid_emergency_date) | postvax_covid_emergency_date>=vax2_date,
    hospitalisation_date_check = is.na(postvax_covid_hospitalisation_date) | postvax_covid_hospitalisation_date>=vax2_date,
    death_date_check = is.na(postvax_covid_death_date) | postvax_covid_death_date>=vax2_date,
    
    # Population exclusions
    isnot_hscworker = !hscworker,
    isnot_carehomeresident = !care_home,
    isnot_endoflife = !endoflife,
    isnot_housebound = !housebound,
    
    # not censored pre dose 2
    isnot_censored_early = tte_censor>0 | is.na(tte_censor),
    
    include = (
      has_age & has_ckd_any & has_ckd_strict & 
      has_sex & has_imd & has_ethnicity & has_region &
      vax_pfi_az & vax_date & vax_interval &
      positive_test_date_check & emergency_date_check & hospitalisation_date_check & death_date_check &
      isnot_hscworker & isnot_carehomeresident & isnot_endoflife & isnot_housebound &
      isnot_censored_early
     )
  )

## Create cohort data including patients fulfilling selection criteria
data_cohort <- data_criteria %>%
  filter(include) %>%
  select(patient_id) %>%
  left_join(data_processed, by="patient_id") %>%
  select(-c(ckd_inclusion_any, ckd_inclusion_strict)) %>%
  droplevels()

## Save data
write_rds(data_cohort, here::here("output", "data", "data_cohort_VE.rds"), compress="gz")
write_csv(data_cohort, here::here("output", "data", "data_cohort_VE.csv"))

## Create and save flow chart
data_flowchart <- data_criteria %>%
  transmute(
    c0 = has_age & has_ckd_any,
    c1 = c0 & has_ckd_strict,
    c2 = c1 & (has_sex & has_imd & has_ethnicity & has_region),
    c3 = c2 & (vax_pfi_az),
    c4 = c3 & (vax_date),
    c5 = c4 & (vax_interval),
    c6 = c5 & (positive_test_date_check & emergency_date_check & hospitalisation_date_check & death_date_check),
    c7 = c6 & (isnot_hscworker & isnot_carehomeresident & isnot_endoflife & isnot_housebound),
    c8 = c7 & isnot_censored_early
  ) %>%
  summarise(
    across(.fns=sum)
  ) %>%
  pivot_longer(
    cols=everything(),
    names_to="criteria",
    values_to="n"
  ) %>%
  mutate(
    n_exclude = lag(n) - n,
    pct_exclude = n_exclude/lag(n),
    pct_all = n / first(n),
    pct_step = n / lag(n),
    crit = str_extract(criteria, "^c\\d+"),
    criteria = fct_case_when(
      crit == "c0" ~ "Aged 16+ with eGFR<60 in 2 years before 01 Dec 2020 or dialysis code, kidney transplant code, CKD diagnostic code, or CKD3-5 code", # paste0("Aged 18+\n with 2 doses on or before ", format(study_dates$lastvax2_date, "%d %b %Y")),
      crit == "c1" ~ "  with eGFR<60 or dialysis/kidney transplant code",
      crit == "c2" ~ "  with no missing demographic information",
      crit == "c3" ~ "  received 2 x ChAdOx1-S or 2 x BNT162b2",
      crit == "c4" ~ "  received first dose after 04 January 2021",
      crit == "c5" ~ "  dose interval of 8-16 weeks",
      crit == "c6" ~ "  post-vaccination outcomes recorded after second dose",
      crit == "c7" ~ "  not healthcare worker, care home resident, receiving end-of-life care, or housebound",
      crit == "c8" ~ "  not censored before second dose",
      TRUE ~ NA_character_
    )
  )
write_csv(data_flowchart, here::here("output", "tables", "flowchart_VE.csv"))
