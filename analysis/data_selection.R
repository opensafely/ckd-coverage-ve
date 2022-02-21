
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
    has_ckd_any = ckd_inclusion_any==1,
    has_ckd_strict = ckd_inclusion_strict==1,
    
    has_age = !is.na(age) & age >=16 & age<120,
    has_sex = !is.na(sex),
    has_imd = !is.na(imd),
    has_ethnicity = !is.na(ethnicity),
    has_region = !is.na(region),
    #isnot_hscworker = !hscworker,
    #isnot_carehomeresident = !care_home_combined,
    #isnot_endoflife = !endoflife,
    #isnot_housebound = !housebound,
    has_max_4_vax = n_vax <= 4, # maximum of 4 recorded doses
    has_vaxgap12 = tbv1_2 >= 14 | is.na(vax2_date), # at least 14 days between dose 1 and dose 2 if dose 2 given
    has_vaxgap23 = tbv2_3 >= 14 | is.na(vax3_date), # at least 14 days between dose 2 and dose 3 if dose 3 given
    has_vaxgap34 = tbv3_4 >= 14 | is.na(vax4_date), # at least 14 days between dose 3 and dose 4 f dose 3 given
    #has_knownvax1 = vax1_type %in% c("pfizer", "az", "moderna"),
    #has_knownvax2 = vax2_type %in% c("pfizer", "az", "moderna"),
    #has_knownvax3 = vax3_type %in% c("pfizer", "az", "moderna"),
    #has_knownvax4 = vax4_type %in% c("pfizer", "az", "moderna"),
    
    include = (
      #jcvi_group_6orhigher & # temporary until more data available
      has_ckd_any & has_age & 
      has_ckd_strict & 
      has_sex & has_imd & has_ethnicity & has_region &
      has_max_4_vax & 
      has_vaxgap12 & has_vaxgap23 & has_vaxgap34
     )
  )

data_cohort <- data_criteria %>%
  filter(include) %>%
  select(patient_id) %>%
  left_join(data_processed, by="patient_id") %>%
  select(-c(ckd_inclusion_any, ckd_inclusion_strict)) %>%
  droplevels()

write_rds(data_cohort, here::here("output", "data", "data_cohort.rds"), compress="gz")
write_csv(data_cohort, here::here("output", "data", "data_cohort.csv"))

data_flowchart <- data_criteria %>%
  transmute(
    c0 = has_age & has_ckd_any,
    c1 = c0 & has_ckd_strict,
    c2 = c1 & (has_sex & has_imd & has_ethnicity & has_region),
    c3 = c2 & (has_max_4_vax),
    c4 = c3 & (has_vaxgap12 & has_vaxgap23 & has_vaxgap34)
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
      crit == "c3" ~ "  with maximum of 4 doses recorded",
      crit == "c4" ~ "  with no vaccines administered at an interval of <14 days",
      TRUE ~ NA_character_
    )
  )
write_csv(data_flowchart, here::here("output", "data", "flowchart.csv"))
