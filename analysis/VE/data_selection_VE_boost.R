
# # # # # # # # # # # # # # # # # # # # #
# This script:
# imports processed data
# filters out people who are excluded from the main analysis
# defines primary (unmatched) and secondary (matched) study populations
# outputs inclusion/exclusions flowchart data
# # # # # # # # # # # # # # # # # # # # #

## Import libraries ----
library('tidyverse')
library('lubridate')
library('here')
library('glue')
library('MatchIt')

## Import custom user functions
source(here::here("analysis", "functions.R"))

## Create output directory
fs::dir_create(here::here("output", "tables"))

## Import processed data
data_processed <- read_rds(here::here("output", "data", "data_processed_VE_boost.rds"))  %>%
  # to avoid timestamps causing inequalities for dates on the same day
  mutate(across(where(is.Date), 
                ~ floor_date(
                  as.Date(.x, format="%Y-%m-%d"),
                  unit = "days"))) 

## Vaccine initiation dates
first_dose_min = as_date("2021-01-04")
third_dose_min = as_date("2021-09-01")

## Set analysis end date
# 2 months after launch of spring 2022 booster campaign 
# https://www.england.nhs.uk/2022/03/nhs-covid-19-vaccine-programme-delivers-first-spring-boosters/
data_processed$end_date = as_date("2022-05-21") 
data_processed$end_date_testing = as_date("2022-03-31") 

## Set and store outcomes list
outcomes_list <- list(
  short_name = c("covid_postest", "covid_emergency", "covid_hosp", "covid_death", "noncovid_death"),
clean_name = c("Positive SARS-CoV-2 test", "COVID-19-related A&E admission", "COVID-19-related hospitalisation", "COVID-19-related death", "Non-COVID-19 death"),
  short_name = c("SARS-CoV-2+", "COVID-19 A&E", "COVID-19 hosp.", "COVID-19 death", "Non-COVID-19 death"),
  date_name = c("postvax_positive_test_date", "postvax_covid_emergency_date", "postvax_covid_hospitalisation_date", "postvax_covid_death_date", "noncoviddeath_date")
)
dir.create(here::here("output", "lib"), showWarnings = FALSE, recursive=TRUE)
write_rds(
  outcomes_list,
  here::here("output", "lib", "outcomes_boost.rds")
)

## Set and store analysis intervals and last follow-up day
period_length <- 56 # match primary VE analysis
n_periods <- 3
postvaxcuts <- c(0,period_length*0:(n_periods)+14)
postvax_periods = paste0(postvaxcuts[1:((length(postvaxcuts)-1))]+1,"-",postvaxcuts[2:length(postvaxcuts)])
postvax_periods_weeks = paste0((postvaxcuts/7)[1:((length(postvaxcuts)-1))]+1,"-",(postvaxcuts/7)[2:length(postvaxcuts)])
lastfupday <- max(postvaxcuts)
write_rds(
  list(
    postvaxcuts = postvaxcuts,
    postvax_periods = postvax_periods,
    postvax_periods_weeks = postvax_periods_weeks
  ),
  here::here("output", "lib", "postboost_list.rds")
)

## Create cohort data with tte calculations
data_processed <- data_processed %>%
  mutate(
    # set censor date (last follow-up day, end date, deregistration, death, 4th dose)
    censor_date = pmin(vax3_date - 1 + lastfupday, end_date, dereg_date, death_date, vax4_date, na.rm=TRUE),
    tte_censor = tte(vax3_date - 1, censor_date, censor_date),
    ind_censor = dplyr::if_else((censor_date>censor_date) | is.na(censor_date), FALSE, TRUE),
    
    # set censor date for testing (last follow-up day, end date, deregistration, death, 4th dose)
    censor_date_testing = pmin(vax3_date - 1 + lastfupday, end_date_testing, dereg_date, death_date, vax4_date, na.rm=TRUE),
    tte_censor_testing = tte(vax3_date - 1, censor_date_testing, censor_date),
    ind_censor_testing = dplyr::if_else((censor_date_testing>censor_date_testing) | is.na(censor_date_testing), FALSE, TRUE),
    
    # time to positive test
    tte_covid_postest = tte(vax3_date - 1, postvax_positive_test_date, censor_date_testing, na.censor=TRUE),
    tte_covid_postest_or_censor = tte(vax3_date - 1, postvax_positive_test_date, censor_date_testing, na.censor=FALSE),
    ind_covid_postest = dplyr::if_else((postvax_positive_test_date>censor_date_testing) | is.na(postvax_positive_test_date), FALSE, TRUE),
    
    # time to COVID-19 A&E attendance
    tte_covid_emergency = tte(vax3_date - 1, postvax_covid_emergency_date, censor_date, na.censor=TRUE),
    tte_covid_emergency_or_censor = tte(vax3_date - 1, postvax_covid_emergency_date, censor_date, na.censor=FALSE),
    ind_covid_emergency = dplyr::if_else((postvax_covid_emergency_date>censor_date) | is.na(postvax_covid_emergency_date), FALSE, TRUE),
    
    # time to COVID-19 hospitalisation
    tte_covid_hosp = tte(vax3_date - 1, postvax_covid_hospitalisation_date, censor_date, na.censor=TRUE),
    tte_covid_hosp_or_censor = tte(vax3_date - 1, postvax_covid_hospitalisation_date, censor_date, na.censor=FALSE),
    ind_covid_hosp = dplyr::if_else((postvax_covid_hospitalisation_date>censor_date) | is.na(postvax_covid_hospitalisation_date), FALSE, TRUE),
    
    # time to COVID-19 death
    tte_covid_death = tte(vax3_date - 1, postvax_covid_death_date, censor_date, na.censor=TRUE),
    tte_covid_death_or_censor = tte(vax3_date - 1, postvax_covid_death_date, censor_date, na.censor=FALSE),
    ind_covid_death = dplyr::if_else((postvax_covid_death_date>censor_date) | is.na(postvax_covid_death_date), FALSE, TRUE),
    
    # time to non-COVID-19 death
    tte_noncovid_death = tte(vax3_date - 1, noncoviddeath_date, censor_date, na.censor=TRUE),
    tte_noncovid_death_or_censor = tte(vax3_date - 1, noncoviddeath_date, censor_date, na.censor=FALSE),
    ind_noncovid_death = dplyr::if_else((noncoviddeath_date>censor_date) | is.na(noncoviddeath_date), FALSE, TRUE),

    # time dose 3 to cut-off
    tte_dose3_to_cutoff = as.numeric(date(end_date)-date(vax3_date-1)),
    tte_dose3_to_cutoff_testing = as.numeric(date(end_date_testing)-date(vax3_date-1))
  )

###################################
### Unmatched data selection
###################################

# Define selection criteria ----
data_criteria <- data_processed %>%
  transmute(
    patient_id,
    
    # Made it into into study population with valid age
    study_definition = TRUE,
    has_age = !is.na(age) & age >=16 & age<120,
    
    # CKD inclusion criteria 
    has_valid_creatinine_or_ukrr = (creatinine_date_issue==0 & creatinine_operator_issue==0) | ukrr_2020_group=="Tx" | ukrr_2020_group=="Dialysis",
    has_ckd_egfr_ukrr = ckd_inclusion_egfr_ukrr==1,
    has_no_rrt_mismatch = rrt_mismatch==0,
    
    # Demography
    has_sex = !is.na(sex),
    has_imd = !is.na(imd),
    has_ethnicity = !is.na(ethnicity),
    has_region = !is.na(region),
    
    # Vaccine profile
    vax_homol_heterol_pfi = (!is.na(vax123_type)) & (vax123_type=="az-az-pfizer" | vax123_type=="pfizer-pfizer-pfizer"),
    vax_date_valid = (!is.na(vax1_date)) & vax1_date>=first_dose_min & (!is.na(vax3_date)) & vax3_date>=third_dose_min,
    vax_interval_valid_1_2 = (!is.na(tbv1_2)) & tbv1_2>=(8*7) & tbv1_2<=(14*7),
    vax_interval_valid_2_3 = (!is.na(tbv2_3)) & tbv2_3>=(12*7),
    
    # Population exclusions
    isnot_hscworker = !hscworker,
    isnot_carehomeresident = !care_home,
    isnot_endoflife = !endoflife,
    isnot_housebound = !housebound,
    isnot_JCVI2 = jcvi_group != "2 (80+ or health/social care worker)",
    #isnot_inhospital = !inhospital,
    
    # Postvax events
    positive_test_date_check = is.na(postvax_positive_test_date) | postvax_positive_test_date>=vax3_date,
    hospitalisation_date_check = is.na(postvax_covid_hospitalisation_date) | postvax_covid_hospitalisation_date>=vax3_date,
    death_date_check = is.na(postvax_covid_death_date) | postvax_covid_death_date>=vax3_date,
    noncoviddeath_date_check = is.na(noncoviddeath_date) | noncoviddeath_date>=vax3_date,
    
    # Not censored pre dose 3
    isnot_censored_early = tte_censor>0 | is.na(tte_censor),
    
    # No COVID in window spanning dose 1 to dose 3
    nointervax_covid = intervax_covid_cat==0,
    
    # Primary outcome study population
    include = (
      has_age &
      has_valid_creatinine_or_ukrr & has_ckd_egfr_ukrr & has_no_rrt_mismatch &
      has_sex & has_imd & has_ethnicity & has_region &
      vax_homol_heterol_pfi & vax_date_valid & vax_interval_valid_1_2 & vax_interval_valid_2_3 &
      isnot_hscworker & isnot_carehomeresident & isnot_endoflife & isnot_housebound & isnot_JCVI2 & #isnot_inhospital &
      positive_test_date_check & hospitalisation_date_check & death_date_check & noncoviddeath_date_check & isnot_censored_early &
      nointervax_covid
     )
  )

## Create cohort data including patients fulfilling selection criteria
data_cohort <- data_criteria %>%
  filter(include) %>%
  select(patient_id) %>%
  left_join(data_processed, by="patient_id") %>%
  select(-starts_with("ckd_inclusion_")) %>%
  droplevels() %>%
  # Additional vaccine/time covariates
  mutate(
    vax1_day = as.integer(floor((vax1_date - first_dose_min))+1), # day 1 is the day first dose 1 given
    vax2_day = as.integer(floor((vax2_date - first_dose_min))+1), # day 1 is the day first dose 1 given
    vax3_day = as.integer(floor((vax3_date - first_dose_min))+1), # day 1 is the day first dose 1 given
    vax1_week = as.integer(floor((vax1_date - first_dose_min)/7)+1), # week 1 is days 1-7
    vax2_week = as.integer(floor((vax2_date - first_dose_min)/7)+1), # week 1 is days 1-7
    vax3_week = as.integer(floor((vax3_date - first_dose_min)/7)+1), # week 1 is days 1-7
    week_region = paste0(vax3_week, "_", region),
    vax2_az = (vax2_type=="az")*1 # since comparison group is related to primary schedule (homologous vs heterologous), vax2_az can still be used
  )

## Save data
write_rds(data_cohort, here::here("output", "data", "data_cohort_VE_boost.rds"), compress="gz")
write_csv(data_cohort, here::here("output", "data", "data_cohort_VE_boost.csv"))

## Create and save flow chart
data_flowchart <- data_criteria %>%
  transmute(
    c0 = study_definition & has_age,
    c1 = c0 & has_valid_creatinine_or_ukrr,
    c2 = c1 & has_ckd_egfr_ukrr,
    c3 = c2 & has_no_rrt_mismatch,
    c4 = c3 & (has_sex & has_imd & has_ethnicity & has_region),
    c5 = c4 & (vax_homol_heterol_pfi),
    c6 = c5 & (vax_date_valid),
    c7 = c6 & (vax_interval_valid_1_2 & vax_interval_valid_2_3),
    c8 = c7 & (isnot_hscworker & isnot_carehomeresident & isnot_endoflife & isnot_housebound & isnot_JCVI2),
    c9 = c8 &  (positive_test_date_check & hospitalisation_date_check & death_date_check & noncoviddeath_date_check & isnot_censored_early),
    c10 = c9 & nointervax_covid,
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
      crit == "c0" ~ "Aged >=16 years on 31st March 2021 with any serum creatinine measurement in 2 years preceding definition date or in UK Renal Registry on 31st December 2020",
      crit == "c1" ~ "Valid record for most recent creatinine measurement (with associated date and no linked operators) or in UK Renal Registry on 31st December 2020", 
      crit == "c2" ~ "eGFR <60 ml/min/1.73 m2 based on most recent creatinine measurement or in UK Renal Registry on 31st December 2020", 
      crit == "c3" ~ "No RRT status mismatch (primary care code indicating dialysis or kidney transplant but not in UK Renal Registry on 31st December 2020)", 
      crit == "c4" ~ "No missing demographic information (sex, region, index of multiple deprivation, or ethnicity)",
      crit == "c5" ~ "Received AZ-AZ-BNT/BNT-BNT-BNT",
      crit == "c6" ~ "Received first dose on or after 4th January 2021 and third dose on or after 1st September 2021",
      crit == "c7" ~ "Dose 1-2 interval of 8-14 weeks and dose 2-3 interval of >=12 weeks",
      crit == "c8" ~ "Not healthcare worker, care home resident, receiving end-of-life care, housebound, or in JCVI priority group 2",
      crit == "c9" ~ "No outcome or censoring events recorded before start of follow-up",
      crit == "c10" ~ "No documented SARS-CoV-2 infection between doses 1 and 3",
      TRUE ~ NA_character_
    )
  )
write_csv(data_flowchart, here::here("output", "tables", "flowchart_VE_boost.csv"))





###################################
### Matched data selection
###################################

## If running locally, inflate unmatched cohort and use inflated data for both matched/unmatched analyses
if(Sys.getenv("OPENSAFELY_BACKEND") %in% c("", "expectations")) {
  ## Remove 'No CKD' cat, which may be retained in dummy data but would be excluded if run on server
  data_cohort = subset(data_cohort, ckd_6cat!="No CKD") %>% droplevels()
  
  ## Increase size of data cohort by 10-fold to assist with later model fitting
  data_cohort = rbind(data_cohort, data_cohort, data_cohort, data_cohort, data_cohort, data_cohort,
                       data_cohort, data_cohort, data_cohort, data_cohort, data_cohort, data_cohort)
  
  ## Assign unique patient IDs to inflated cohort
  data_cohort$patient_id = sample.int(100000, nrow(data_cohort), replace = FALSE)
  
  ## Override data
  write_rds(data_cohort, here::here("output", "data", "data_cohort_VE_boost.rds"), compress="gz")
  write_csv(data_cohort, here::here("output", "data", "data_cohort_VE_boost.csv"))
  
  ## Save inflated cohort into unmatched data 
  write_rds(data_cohort, here::here("output", "data", "data_cohort_VE_boost_matched.rds"), compress="gz")
  write_csv(data_cohort, here::here("output", "data", "data_cohort_VE_boost_matched.csv"))
  
  ## Create dummy matched selection flowchart
  data_flowchart = data_flowchart[1:2,]
  data_flowchart$criteria = c("Unmatched VE cohort", "Matched VE cohort")
  data_flowchart$n = max(data_flowchart$n)
  write_csv(data_flowchart, here::here("output", "tables", "flowchart_VE_boost_matched.csv"))
  
} else {
  
  ## Specify exact matching variables
  exact_variables <- c(
    "region",
    "imd",
    "sex",
    "ckd_5cat",
    "cev",
    "prior_covid_cat",
    "any_immunosuppression",    
    NULL
  )
  
  ## Specify caliper variables
  caliper_variables <- c(
    age = 3,
    vax3_day = 3,
    vax2_day = 14,
    NULL
  )
  
  ## Define full list of matching variables
  matching_variables <- c(exact_variables, names(caliper_variables))
  
  ## Run matching algorithm
  safely_matchit <- purrr::safely(matchit)
  matching <-
    safely_matchit(
      formula = vax2_az ~ 1,
      data = data_cohort,
      method = "nearest", distance = "glm", # two options redundant since we are using exact + caliper matching
      replace = FALSE,
      estimand = "ATT",
      exact = exact_variables,
      caliper = caliper_variables, std.caliper=FALSE,
      m.order = "data", # data is sorted on (effectively random) patient ID
      ratio = 1L
    )[[1]]
  
  ## Print summary of matching process
  summary(matching)
  
  ## Pick out matched candidates
  data_matched <-
    as.data.frame(matching$X) %>%
    add_column(
      match_id = matching$subclass,
      treated = matching$treat,
      patient_id = data_cohort$patient_id,
      weight = matching$weights
    ) %>%
    filter(!is.na(match_id)) %>% # remove unmatched people. equivalent to weight != 0
    arrange(match_id, desc(treated)) %>%
    left_join(
      data_cohort %>% select(-matching_variables),
      by = "patient_id"
    ) 
  
  ## Save data
  write_rds(data_matched, here::here("output", "data", "data_cohort_VE_boost_matched.rds"), compress="gz")
  write_csv(data_matched, here::here("output", "data", "data_cohort_VE_boost_matched.csv"))
  
  ## Define selection criteria
  data_criteria <- data_cohort %>%
    transmute(
      patient_id,
      unmatched = TRUE,
      has_match = patient_id %in% data_matched$patient_id,
      include = (has_match)
    )
  
  ## Create and save flow chart
  data_flowchart <- data_criteria %>%
    transmute(
      c0 = unmatched,
      c1 = c0 & has_match
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
        crit == "c0" ~ "Unmatched VE cohort",
        crit == "c1" ~ "Matched VE cohort",
        TRUE ~ NA_character_
      )
    )
  write_csv(data_flowchart, here::here("output", "tables", "flowchart_VE_boost_matched.csv"))
}
