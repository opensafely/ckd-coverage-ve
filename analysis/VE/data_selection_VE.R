
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

# vaccine initiation dates
first_pfizer = as_date("2020-12-08")
first_az = as_date("2021-01-04")
first_moderna = as_date("2021-04-13")

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
    
    # Population exclusions
    isnot_hscworker = !hscworker,
    isnot_carehomeresident = !care_home,
    isnot_endoflife = !endoflife,
    isnot_housebound = !housebound,
    
    include = (
      has_age & has_ckd_any & has_ckd_strict & 
      has_sex & has_imd & has_ethnicity & has_region &
      vax_pfi_az & vax_date & vax_interval &
      isnot_hscworker & isnot_carehomeresident & isnot_endoflife & isnot_housebound
     )
  )

data_cohort <- data_criteria %>%
  filter(include) %>%
  select(patient_id) %>%
  left_join(data_processed, by="patient_id") %>%
  select(-c(ckd_inclusion_any, ckd_inclusion_strict)) %>%
  droplevels()

write_rds(data_cohort, here::here("output", "data", "data_cohort_VE.rds"), compress="gz")
write_csv(data_cohort, here::here("output", "data", "data_cohort_VE.csv"))

data_flowchart <- data_criteria %>%
  transmute(
    c0 = has_age & has_ckd_any,
    c1 = c0 & has_ckd_strict,
    c2 = c1 & (has_sex & has_imd & has_ethnicity & has_region),
    c3 = c2 & (vax_pfi_az),
    c4 = c3 & (vax_date),
    c5 = c4 & (vax_interval),
    c6 = c5 & (isnot_hscworker & isnot_carehomeresident & isnot_endoflife & isnot_housebound)
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
      crit == "c6" ~ "  not healthcare worker, care home resident, receiving end-of-life care, or housebound",
      TRUE ~ NA_character_
    )
  )
write_csv(data_flowchart, here::here("output", "tables", "flowchart_VE.csv"))
