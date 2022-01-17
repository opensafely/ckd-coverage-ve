######################################

# This script:
# - imports data extracted by the cohort extractor
# - combines ethnicity columns
# - standardises some variables (eg convert to factor) and derives some new ones
# - saves processed dataset

######################################


# Preliminaries ----

## Import libraries
library('tidyverse')
library('lubridate')

## Custom functions
fct_case_when <- function(...) {
  # uses dplyr::case_when but converts the output to a factor,
  # with factors ordered as they appear in the case_when's  ... argument
  args <- as.list(match.call())
  levels <- sapply(args[-1], function(f) f[[3]])  # extract RHS of formula
  levels <- levels[!is.na(levels)]
  factor(dplyr::case_when(...), levels=levels)
}

## Output processed data to rds
dir.create(here::here("output", "data"), showWarnings = FALSE, recursive=TRUE)


# Process data ----

## Print variable names
read_csv(here::here("output", "data", "input.csv"),
         n_max = 0,
         col_types = cols()) %>%
  names() %>%
  print()

## Read in data (don't rely on defaults)
data_extract0 <- read_csv(
  here::here("output", "data", "input.csv"),
  col_types = cols_only(
    
    # Identifier
    patient_id = col_integer(),
    
    # Vaccination dates
    covid_vacc_date_1 = col_date(format="%Y-%m-%d"),
    covid_vacc_date_2 = col_date(format="%Y-%m-%d"),
    covid_vacc_date_3 = col_date(format="%Y-%m-%d"),
    covid_vacc_date_4 = col_date(format="%Y-%m-%d"),
    pfizer_date_1 = col_date(format="%Y-%m-%d"),
    pfizer_date_2 = col_date(format="%Y-%m-%d"),
    pfizer_date_3 = col_date(format="%Y-%m-%d"),
    pfizer_date_4 = col_date(format="%Y-%m-%d"),
    astrazeneca_date_1 = col_date(format="%Y-%m-%d"),
    astrazeneca_date_2 = col_date(format="%Y-%m-%d"),
    astrazeneca_date_3 = col_date(format="%Y-%m-%d"),
    astrazeneca_date_4 = col_date(format="%Y-%m-%d"),
    moderna_date_1 = col_date(format="%Y-%m-%d"),
    moderna_date_2 = col_date(format="%Y-%m-%d"),
    moderna_date_3 = col_date(format="%Y-%m-%d"),
    moderna_date_4 = col_date(format="%Y-%m-%d"),
    
    # CKD groups
    chronic_kidney_disease_diagnostic = col_date(format="%Y-%m-%d"),
    chronic_kidney_disease_all_stages = col_date(format="%Y-%m-%d"),
    chronic_kidney_disease_all_stages_3_5 = col_date(format="%Y-%m-%d"),
    end_stage_renal = col_logical(), 
    
    # Priority groups
    care_home = col_logical(),
    shielded = col_logical(),
    hscworker = col_logical(),
    
    # Clinical/demographic variables
    age = col_integer(),
    sex = col_character(),
    bmi = col_character(),
    smoking_status =  col_character(),
    ethnicity_6 = col_character(),
    ethnicity_6_sus = col_character(),
    imd = col_character(),
    region = col_character(),
    asthma = col_logical(),
    asplenia = col_logical(),
    bp_sys = col_double(),
    bp_dias = col_double(),
    cancer = col_logical(),
    haem_cancer = col_logical(),
    chd = col_logical(),
    chronic_neuro_dis_inc_sig_learn_dis = col_logical(),
    chronic_resp_dis = col_logical(),
    cld = col_logical(),
    diabetes = col_logical(),
    immunosuppression_diagnosis_date = col_date(format="%Y-%m-%d"),
    immunosuppression_medication_date = col_date(format="%Y-%m-%d"),
    learning_disability = col_logical(),
    sev_mental_ill = col_date(format="%Y-%m-%d"),
    organ_transplant = col_logical(),
    
    # Other
    prior_positive_test_date = col_date(format="%Y-%m-%d"),
    prior_primary_care_covid_case_date = col_date(format="%Y-%m-%d"),
    prior_covidadmitted_date = col_date(format="%Y-%m-%d"),
    tests_conducted_any = col_double(),
    tests_conducted_positive = col_double()
    
  ),
  na = character() # more stable to convert to missing later
)

## Parse NAs
data_extract <- data_extract0 %>%
  mutate(across(
    .cols = where(is.character),
    .fns = ~na_if(.x, "")
  )) %>%
  mutate(across(
    .cols = c(where(is.numeric), -ends_with("_id")), #convert numeric+integer but not id variables
    .fns = ~na_if(.x, 0)
  )) %>%
  arrange(patient_id) 

## Format columns (i.e, set factor levels)
data_processed <- data_extract %>%
  mutate(
    # Care home (65+)
    care_home_65plus = ifelse(care_home == 1 & age >=65, 1, 0),
    
    # Shielding
    shielded = ifelse(shielded == 1 & (age >=16 & age < 70), 1, 0),
    
    # Age
    ageband = cut(
      age,
      breaks = c(16, 70, 80, Inf),
      labels = c("16-69", "70-79", "80+"),
      right = FALSE
    ),
    
    ageband = ifelse(care_home == 1, NA, ageband),
    
    ageband2 = cut(
      age,
      breaks = c(16, 80, 85, 90, 95, Inf),
      labels = c("16-79", "80-84", "85-89", "90-94", "95+"),
      right = FALSE
    ),
    
    # Sex
    sex = fct_case_when(
      sex == "F" ~ "Female",
      sex == "M" ~ "Male",
      #sex == "I" ~ "Inter-sex",
      #sex == "U" ~ "Unknown",
      TRUE ~ NA_character_
    ),
    
    # BMI
    bmi = ifelse(bmi == "Not obese","Not obese", "Obese"),
    
    # Smoking status
    smoking_status = ifelse(smoking_status == "E" | smoking_status == "S","S&E", smoking_status),
    smoking_status = ifelse(smoking_status == "S&E", "S&E", "N&M"),
    
    # Ethnicity
    ethnicity_filled = ifelse(is.na(ethnicity_6), ethnicity_6_sus, ethnicity_6),
    ethnicity = ifelse(is.na(ethnicity_filled), 6, ethnicity_filled),
    
    ethnicity = fct_case_when(
      ethnicity == "1" ~ "White",
      ethnicity == "2" ~ "Mixed",
      ethnicity == "3" ~ "Asian or Asian British",
      ethnicity == "4" ~ "Black or Black British",
      ethnicity == "5" ~ "Other ethnic groups",
      ethnicity == "6" ~ "Unknown",
      #TRUE ~ "Unknown"
      TRUE ~ NA_character_),
    
    # IMD
    imd = na_if(imd, "0"),
    imd = fct_case_when(
      imd == 1 ~ "1 most deprived",
      imd == 2 ~ "2",
      imd == 3 ~ "3",
      imd == 4 ~ "4",
      imd == 5 ~ "5 least deprived",
      #TRUE ~ "Unknown",
      TRUE ~ NA_character_
    ),
    
    # Region
    region = fct_case_when(
      region == "London" ~ "London",
      region == "East" ~ "East of England",
      region == "East Midlands" ~ "East Midlands",
      region == "North East" ~ "North East",
      region == "North West" ~ "North West",
      region == "South East" ~ "South East",
      region == "South West" ~ "South West",
      region == "West Midlands" ~ "West Midlands",
      region == "Yorkshire and The Humber" ~ "Yorkshire and the Humber",
      #TRUE ~ "Unknown",
      TRUE ~ NA_character_),
    
    ## Blood pressure
    bpcat = ifelse(bp_sys < 120 &  bp_dias < 80, 1, NA),
    bpcat = ifelse(bp_sys >= 120 & bp_sys < 130 & bp_dias < 80, 2, bpcat),
    bpcat = ifelse(bp_sys >= 130 &  bp_dias >= 90, 3, bpcat),
    bpcat = ifelse(is.na(bpcat), 4, bpcat),
    
    bpcat = fct_case_when(
      bpcat == 1 ~ "Normal",
      bpcat == 2 ~ "Elevated",
      bpcat == 3 ~ "High",
      bpcat == 4 ~ "Unknown",
      #TRUE ~ "Unknown",
      TRUE ~ NA_character_
    ),
    
    # CKD
    chronic_kidney_disease = case_when(
      !is.na(chronic_kidney_disease_diagnostic) ~ TRUE,
      is.na(chronic_kidney_disease_all_stages) ~ FALSE,
      !is.na(chronic_kidney_disease_all_stages_3_5) & (chronic_kidney_disease_all_stages_3_5 >= chronic_kidney_disease_all_stages) ~ TRUE,
      TRUE ~ FALSE
    ),
    
    # Immunosuppression
    immunosuppression_diagnosis_date = !is.na(immunosuppression_diagnosis_date),
    immunosuppression_medication_date = !is.na(immunosuppression_medication_date),
    immunosuppression = immunosuppression_diagnosis_date | immunosuppression_medication_date,
    
    # Mental illness
    sev_mental_ill = !is.na(sev_mental_ill),
    
    # Time between vaccinations
    tbv1_2 = as.numeric(covid_vacc_date_2 - covid_vacc_date_1),
    tbv2_3 = as.numeric(covid_vacc_date_3 - covid_vacc_date_2),
    tbv3_4 = as.numeric(covid_vacc_date_4 - covid_vacc_date_3),
    
    # Prior covid
    prior_covid_date = pmin(prior_positive_test_date, 
                           prior_primary_care_covid_case_date, 
                           prior_covidadmitted_date,
                           na.rm=TRUE), 
    
    prior_covid_cat = ifelse(prior_covid_date < covid_vacc_date_3 & prior_covid_date >= covid_vacc_date_2, 1, NA),
    prior_covid_cat = ifelse(prior_covid_date < covid_vacc_date_2 & prior_covid_date >= covid_vacc_date_1, 2, NA),
    prior_covid_cat = ifelse(prior_covid_date < covid_vacc_date_1, 3, prior_covid_cat),
    
    prior_covid_cat = fct_case_when(
      prior_covid_cat == 1 ~ "Between second and third dose",
      prior_covid_cat == 2 ~ "Between first and second dose",
      prior_covid_cat == 3 ~ "Prior to first dose",
      #TRUE ~ "Unknown",
      TRUE ~ NA_character_
    )
    
  ) %>%
  droplevels() %>%
  mutate(
    # converts TRUE/FALSE to 1/0
    across(
      where(is.logical),
      ~.x*1L
    )
  ) %>%
  filter(tbv1_2>0 | is.na(tbv1_2), 
         tbv2_3>0 | is.na(tbv2_3),
         tbv3_4>0 | is.na(tbv3_4),
         age >= 18,
         age < 120,
         !is.na(sex),
         chronic_kidney_disease == 1)

# Save dataset as .rds files ----
write_rds(data_processed, here::here("output", "data", "data_processed.rds"), compress = "gz")
write_csv(data_processed, here::here("output", "data", "data_processed.csv"))








