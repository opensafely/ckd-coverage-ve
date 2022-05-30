
# # # # # # # # # # # # # # # # # # # # #
# This script:
# imports processed data
# filters out people who are excluded from the main analysis
# outputs inclusion/exclusions flowchart data
# # # # # # # # # # # # # # # # # # # # #

## Import libraries ----
library('tidyverse')
library('lubridate')
library('here')
library('glue')

## Import custom user functions
source(here::here("analysis", "functions.R"))

## Create output directory
fs::dir_create(here::here("output", "tables"))

## Import processed data
data_processed <- read_rds(here::here("output", "data", "data_processed.rds")) 

## Define selection criteria
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
    has_max_5_vax = n_vax <= 5, # maximum of 5 recorded doses
    valid_vaxgap12 = tbv1_2 >= 14 | is.na(vax2_date), # at least 14 days between dose 1 and dose 2 if dose 2 given
    valid_vaxgap23 = tbv2_3 >= 14 | is.na(vax3_date), # at least 14 days between dose 2 and dose 3 if dose 3 given
    valid_vaxgap34 = tbv3_4 >= 14 | is.na(vax4_date), # at least 14 days between dose 3 and dose 4 if dose 4 given
    
    # Dose 4 denominator
    eligible_dose4 = (age>=75 | care_home==1 | haem_cancer==1 | immunosuppression==1 | ckd_5cat=="RRT (Tx)" | other_transplant==1), 
    
    # Other
    alive_throughout = is.na(death_date),
    registered_throughout = is.na(dereg_date),
    
    # Primary outcome study population
    include = (
      has_age & 
      has_valid_creatinine_or_ukrr & has_ckd_egfr_ukrr & has_no_rrt_mismatch &
      has_sex & has_imd & has_ethnicity & has_region &
      has_max_5_vax & valid_vaxgap12 & valid_vaxgap23 & valid_vaxgap34
    ),

    # Logistic regression population (sensitivity analysis)
    include_logistic = (
      include & alive_throughout & registered_throughout
    ),
    
    # Secondary outcome population (dose 4)
    include_dose4 = (
      include & eligible_dose4
    ),
  )

## Define primary outcome population
data_cohort <- data_criteria %>%
  filter(include) %>%
  select(patient_id) %>%
  left_join(data_processed, by="patient_id") %>%
  select(-starts_with("ckd_inclusion_")) %>%
  droplevels()

## Save dataset
write_rds(data_cohort, here::here("output", "data", "data_cohort_coverage.rds"), compress="gz")
write_csv(data_cohort, here::here("output", "data", "data_cohort_coverage.csv"))

## Define logistic regression population (sensitivity analysis)
data_cohort_logistic <- data_criteria %>%
  filter(include_logistic) %>%
  select(patient_id) %>%
  left_join(data_processed, by="patient_id") %>%
  select(-starts_with("ckd_inclusion_")) %>%
  droplevels()

## Save dataset
write_rds(data_cohort_logistic, here::here("output", "data", "data_cohort_coverage_logistic.rds"), compress="gz")

## Define secondary outcome population (dose 4)
data_cohort_dose4 <- data_criteria %>%
  filter(include_dose4) %>%
  select(patient_id) %>%
  left_join(data_processed, by="patient_id") %>%
  select(-starts_with("ckd_inclusion_")) %>%
  droplevels()

## Save dataset
write_rds(data_cohort_dose4, here::here("output", "data", "data_cohort_coverage_dose4.rds"), compress="gz")

## Create flowchart of inclusions/exclusions
data_flowchart <- data_criteria %>%
  transmute(
    # Primary outcome study population
    c0 = study_definition & has_age,
    c1 = c0 & has_valid_creatinine_or_ukrr,
    c2 = c1 & has_ckd_egfr_ukrr,
    c3 = c2 & has_no_rrt_mismatch,
    c4 = c3 & (has_sex & has_imd & has_ethnicity & has_region),
    c5 = c4 & (has_max_5_vax),
    c6 = c5 & (valid_vaxgap12 & valid_vaxgap23 & valid_vaxgap34),
    # Logistic regression population (sensitivity analysis)
    c7 = c6 & alive_throughout,
    c8 = c7 & registered_throughout,
    # Secondary outcome population (dose 4)
    c9 = c6 & eligible_dose4
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
    crit = str_extract(criteria, "^c\\d+"),
    criteria = fct_case_when(
      crit == "c0" ~ "Aged 16+ with serum creatinine record in 2y before index (01-Dec-2020) or in UKRR 2020 population",
      crit == "c1" ~ "  with valid creatinine record (with associated age and no linked operators) or in UKRR population", 
      crit == "c2" ~ "  with eGFR<60 or UKRR 2020 population", 
      crit == "c3" ~ "  with no RRT mismatch (primary care dialysis/Tx code but absent from UKRR 2020 population)", 
      crit == "c4" ~ "  with no missing demographic information",
      crit == "c5" ~ "  with maximum of 5 doses recorded up to 20 April 2022",
      crit == "c6" ~ "  with no vaccines administered at an interval of <14 days (primary analysis subset)",
      crit == "c7" ~ "  alive throughout follow-up period (logistic regression sensitivity subset only)",
      crit == "c8" ~ "  registered throughout follow-up period (logistic regression sensitivity subset only)",
      crit == "c9" ~ "  eligible for dose 4 based on immunosuppression, transplant, haematologic cancer, care home residence, or age >=75y (dose 4 analysis subset)",
      TRUE ~ NA_character_
    )
  )
write_csv(data_flowchart, here::here("output", "tables", "flowchart_coverage.csv"))
