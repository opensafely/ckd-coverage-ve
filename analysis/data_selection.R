
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
    has_ckd = !is.na(ckd_inclusion) & ckd_inclusion==1,
    has_age = !is.na(age) & age >=16 & age<120,
    has_sex = !is.na(sex),
    has_imd = !is.na(imd),
    has_ethnicity = !is.na(ethnicity),
    has_region = !is.na(region),
    #isnot_hscworker = !hscworker,
    #isnot_carehomeresident = !care_home_combined,
    #isnot_endoflife = !endoflife,
    #isnot_housebound = !housebound,
    # vax_afterfirstvaxdate = case_when(
    #   (vax1_type=="pfizer") & (vax1_date >= first_pfizer) ~ TRUE,
    #   (vax1_type=="az") & (vax1_date >= first_az) ~ TRUE,
    #   (vax1_type=="moderna") & (vax1_date >= first_moderna) ~ TRUE,
    #   (vax2_type=="pfizer") & (vax2_date >= first_pfizer) ~ TRUE,
    #   (vax2_type=="az") & (vax2_date >= first_az) ~ TRUE,
    #   (vax2_type=="moderna") & (vax2_date >= first_moderna) ~ TRUE,
    #   (vax3_type=="pfizer") & (vax3_date >= first_pfizer) ~ TRUE,
    #   (vax3_type=="az") & (vax3_date >= first_az) ~ TRUE,
    #   (vax3_type=="moderna") & (vax3_date >= first_moderna) ~ TRUE,
    #   is.na(vax1_type) ~ TRUE,
    #   TRUE ~ FALSE
    # ),

    #vax2_beforelastvaxdate = !is.na(vax2_date) & (vax2_date <= study_dates$lastvax2_date),
    #vax3_afterstudystartdate = (vax3_date >= study_dates$studystart_date) | is.na(vax3_date),
    #vax3_beforelastvaxdate = (vax3_date <= study_dates$lastvax3_date) & !is.na(vax3_date),
    #vax12_homologous = vax1_type==vax2_type,
    has_vaxgap12 = tbv1_2 >= 14 | is.na(vax2_date), # at least 14 days between dose 1 and dose 2 if dose 2 given
    has_vaxgap23 = tbv2_3 >= 14 | is.na(vax3_date), # at least 14 days between dose 2 and dose 3 if dose 3 given
    has_vaxgap34 = tbv3_4 >= 14 | is.na(vax4_date), # at least 14 days between dose 3 and dose 4 f dose 3 given
    #has_knownvax1 = vax1_type %in% c("pfizer", "az", "moderna"),
    #has_knownvax2 = vax2_type %in% c("pfizer", "az", "moderna"),
    #has_knownvax3 = vax3_type %in% c("pfizer", "az", "moderna"),
    #has_knownvax4 = vax4_type %in% c("pfizer", "az", "moderna"),
    
    #jcvi_group_6orhigher = jcvi_group %in% as.character(1:6),

    include = (
      #jcvi_group_6orhigher & # temporary until more data available
      has_ckd & has_age & 
      has_sex & has_imd & has_ethnicity & has_region &
      #vax_afterfirstvaxdate &
      has_vaxgap12 & has_vaxgap23 & has_vaxgap34 #& 
      #has_knownvax1 & has_knownvax2 & has_knownvax3 & has_knownvax4 & #vax12_homologous &
      #isnot_hscworker &
      #isnot_carehomeresident & isnot_endoflife &
      #isnot_housebound
    )
  )

data_cohort <- data_criteria %>%
  filter(include) %>%
  select(patient_id) %>%
  left_join(data_processed, by="patient_id") %>%
  select(-ckd_inclusion) %>%
  droplevels()

write_rds(data_cohort, here("output", "data", "data_cohort.rds"), compress="gz")
write_csv(data_cohort, here::here("output", "data", "data_cohort.csv"))

data_flowchart <- data_criteria %>%
  transmute(
    c0 = has_age & has_ckd,
    c1 = c0 & (has_sex & has_imd & has_ethnicity & has_region),
    c2 = c1 & (has_vaxgap12 & has_vaxgap23 & has_vaxgap34)
    #c0 = vax1_afterfirstvaxdate & vax2_beforelastvaxdate & vax3_afterstudystartdate & jcvi_group_6orhigher,
    #c1_1yearfup = c0_all & (has_follow_up_previous_year),
    #c1 = c0 & (has_age & has_sex & has_imd & has_ethnicity & has_region),
    #c2 = c1 & (has_vaxgap12 & has_vaxgap23 & has_knownvax1 & has_knownvax2 & vax12_homologous),
    #c3 = c2 & (isnot_hscworker ),
    #c4 = c3 & (isnot_carehomeresident & isnot_endoflife & isnot_housebound),
    #c5 = c4 & vax3_beforelastvaxdate & has_expectedvax3type
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
      crit == "c0" ~ "Aged 16+ with eGFR<60 in 2 years before first dose or any prior dialysis/end-stage renal disease flag", # paste0("Aged 18+\n with 2 doses on or before ", format(study_dates$lastvax2_date, "%d %b %Y")),
      crit == "c1" ~ "  with no missing demographic information",
      crit == "c2" ~ "  with no vaccines administered at an interval of <14 days",
      TRUE ~ NA_character_
    )
  )
write_csv(data_flowchart, here::here("output", "data", "flowchart.csv"))
