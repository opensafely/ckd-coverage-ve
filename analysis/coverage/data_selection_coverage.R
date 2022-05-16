
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

# Define selection criteria ----
data_criteria <- data_processed %>%
  transmute(
    patient_id,
    
    # Made it into into study population
    study_definition = TRUE,
    
    # Age
    has_age = !is.na(age) & age >=16 & age<120,
    has_creatinine_age = !is.na(age_creatinine) | dialysis==1 | kidney_transplant==1 | ukrr_index_group=="Tx" | ukrr_index_group=="Dialysis",
    
    # Sequential CKD
    has_ckd_egfr_ukrr_D_T_3to5_diagnostic = ckd_inclusion_egfr_ukrr_D_T_3to5_diagnostic==1,
    has_ckd_egfr_ukrr_D_T_3to5 = ckd_inclusion_egfr_ukrr_D_T_3to5==1,
    has_ckd_egfr_ukrr_D_T = ckd_inclusion_egfr_ukrr_D_T==1,
    has_ckd_egfr_ukrr = ckd_inclusion_egfr_ukrr==1,
    
    # Demography
    has_sex = !is.na(sex),
    has_imd = !is.na(imd),
    has_ethnicity = !is.na(ethnicity),
    has_region = !is.na(region),

    # Vaccine profile
    has_max_5_vax = n_vax <= 5, # maximum of 5 recorded doses
    valid_vaxgap12 = tbv1_2 >= 14 | is.na(vax2_date), # at least 14 days between dose 1 and dose 2 if dose 2 given
    valid_vaxgap23 = tbv2_3 >= 14 | is.na(vax3_date), # at least 14 days between dose 2 and dose 3 if dose 3 given
    valid_vaxgap34 = tbv3_4 >= 14 | is.na(vax4_date), # at least 14 days between dose 3 and dose 4 f dose 3 given
    
    # Other
    alive_throughout = is.na(death_date),
    registered_throughout = is.na(dereg_date),
    
    include = (
      has_age & has_creatinine_age &
      has_ckd_egfr_ukrr_D_T_3to5_diagnostic & has_ckd_egfr_ukrr_D_T_3to5 & has_ckd_egfr_ukrr_D_T & has_ckd_egfr_ukrr &
      has_sex & has_imd & has_ethnicity & has_region &
      has_max_5_vax & valid_vaxgap12 & valid_vaxgap23 & valid_vaxgap34
    ),

    include_logistic = (
      include & alive_throughout & registered_throughout
    )
  )

## Define data cohort for primary analysis
data_cohort <- data_criteria %>%
  filter(include) %>%
  select(patient_id) %>%
  left_join(data_processed, by="patient_id") %>%
  select(-starts_with("ckd_inclusion_")) %>%
  droplevels()

write_rds(data_cohort, here::here("output", "data", "data_cohort_coverage.rds"), compress="gz")
write_csv(data_cohort, here::here("output", "data", "data_cohort_coverage.csv"))

## Define data cohort for logistic regression (sensitivity analyses)
data_cohort_logistic <- data_criteria %>%
  filter(include_logistic) %>%
  select(patient_id) %>%
  left_join(data_processed, by="patient_id") %>%
  select(-starts_with("ckd_inclusion_")) %>%
  droplevels()

write_rds(data_cohort_logistic, here::here("output", "data", "data_cohort_coverage_logistic.rds"), compress="gz")

data_flowchart <- data_criteria %>%
  transmute(
    c0 = study_definition,
    c1 = c0 & has_age,
    c2 = c1 & has_creatinine_age,
    c3 = c2 & has_ckd_egfr_ukrr_D_T_3to5_diagnostic,
    c4 = c3 & has_ckd_egfr_ukrr_D_T_3to5,
    c5 = c4 & has_ckd_egfr_ukrr_D_T,
    c6 = c5 & has_ckd_egfr_ukrr,
    c7 = c6 & (has_sex & has_imd & has_ethnicity & has_region),
    c8 = c7 & (has_max_5_vax),
    c9 = c8 & (valid_vaxgap12 & valid_vaxgap23 & valid_vaxgap34),
    # Logistic regression analysis cohort
    c10 = c9 & alive_throughout,
    c11 = c10 & registered_throughout
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
      crit == "c0" ~ "Aged 16+ with serum creatinine record in 2y before 01 Dec 2020, or UKRR 2019/2020, dialysis code, kidney transplant code, CKD diagnostic code, or CKD3-5 code",
      crit == "c1" ~ "  with age", 
      crit == "c2" ~ "  with creatinine-associated age or UKRR at index", 
      crit == "c3" ~ "  with eGFR<60 or UKRR at index or any CKD-related code (diagnostic/CDK3-5/dialysis/kidney transplant)", 
      crit == "c4" ~ "  with eGFR<60 or UKRR at index or CDK3-5/dialysis/kidney transplant code", 
      crit == "c5" ~ "  with eGFR<60 or UKRR at index or dialysis/kidney transplant code", 
      crit == "c6" ~ "  with eGFR<60 or UKRR at index", 
      crit == "c7" ~ "  with no missing demographic information",
      crit == "c8" ~ "  with maximum of 5 doses recorded up to 20 April 2022",
      crit == "c9" ~ "  with no vaccines administered at an interval of <14 days (primary analysis subset)",
      crit == "c10" ~ "  alive throughout follow-up period (logistic regression sensitivity subset only)",
      crit == "c11" ~ "  registered throughout follow-up period (logistic regression sensitivity subset only)",
      TRUE ~ NA_character_
    )
  )
write_csv(data_flowchart, here::here("output", "tables", "flowchart_coverage.csv"))
