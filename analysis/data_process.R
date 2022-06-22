######################################

# This script:
# - imports data extracted by the cohort extractor
# - sets variable types (e.g. factors, dates, characters) and derives new variables (eGFR, CKD subgroup, etc)
# - applies additional dummy data processing steps via dummy_data.R script if being run locally (skipped if being run on OpenSAFELY server)
# - saves processed data as data_processed.csv abd data_processed.rds

######################################

## Packages
library('tidyverse')
library('lubridate')
library('here')

## Import custom user functions
source(here::here("analysis", "functions.R"))

## Create directory for processed data
dir.create(here::here("output", "data"), showWarnings = FALSE, recursive=TRUE)

## Print session info to metadata log file
sessionInfo()

## Print variable names
read_csv(here::here("output", "data", "input.csv"),
         n_max = 0,
         col_types = cols()) %>%
  names() %>%
  print()

## Read in data and set variable types
data_extract <- read_csv(
  here::here("output", "data", "input.csv"),
  col_types = cols_only(
    
    # Identifier
    patient_id = col_integer(),
    
    # Vaccination dates
    covid_vax_date_1 = col_date(format="%Y-%m-%d"),
    covid_vax_date_2 = col_date(format="%Y-%m-%d"),
    covid_vax_date_3 = col_date(format="%Y-%m-%d"),
    covid_vax_date_4 = col_date(format="%Y-%m-%d"),
    covid_vax_date_5 = col_date(format="%Y-%m-%d"),
    pfizer_date_1 = col_date(format="%Y-%m-%d"),
    pfizer_date_2 = col_date(format="%Y-%m-%d"),
    pfizer_date_3 = col_date(format="%Y-%m-%d"),
    pfizer_date_4 = col_date(format="%Y-%m-%d"),
    pfizer_date_5 = col_date(format="%Y-%m-%d"),
    az_date_1 = col_date(format="%Y-%m-%d"),
    az_date_2 = col_date(format="%Y-%m-%d"),
    az_date_3 = col_date(format="%Y-%m-%d"),
    az_date_4 = col_date(format="%Y-%m-%d"),
    az_date_5 = col_date(format="%Y-%m-%d"),
    moderna_date_1 = col_date(format="%Y-%m-%d"),
    moderna_date_2 = col_date(format="%Y-%m-%d"),
    moderna_date_3 = col_date(format="%Y-%m-%d"),
    moderna_date_4 = col_date(format="%Y-%m-%d"),
    moderna_date_5 = col_date(format="%Y-%m-%d"),
    
    # Dates for clinical variables
    creatinine_date = col_date(format="%Y-%m-%d"),
    death_date = col_date(format="%Y-%m-%d"),
    dereg_date = col_date(format="%Y-%m-%d"),
    immunosuppression_diagnosis_date = col_date(format="%Y-%m-%d"),
    immunosuppression_medication_date = col_date(format="%Y-%m-%d"),
    
    # Dates for Covid-related variables
    prior_positive_test_date = col_date(format="%Y-%m-%d"),
    prior_primary_care_covid_case_date = col_date(format="%Y-%m-%d"),
    prior_covid_emergency_date = col_date(format="%Y-%m-%d"),
    prior_covid_hospitalisation_date = col_date(format="%Y-%m-%d"),
    prior_positive_test_date_boost = col_date(format="%Y-%m-%d"),
    prior_primary_care_covid_case_date_boost = col_date(format="%Y-%m-%d"),
    prior_covid_emergency_date_boost = col_date(format="%Y-%m-%d"),
    prior_covid_hospitalisation_date_boost = col_date(format="%Y-%m-%d"),
    prevax_positive_test_date	= col_date(format="%Y-%m-%d"),
    prevax_primary_care_covid_case_date	= col_date(format="%Y-%m-%d"),
    prevax_covid_emergency_date	= col_date(format="%Y-%m-%d"),
    prevax_covid_hospitalisation_date	= col_date(format="%Y-%m-%d"),
    postvax_positive_test_date = col_date(format="%Y-%m-%d"),
    postvax_covid_emergency_date = col_date(format="%Y-%m-%d"),
    postvax_covid_hospitalisation_date = col_date(format="%Y-%m-%d"),
    postvax_covid_death_date = col_date(format="%Y-%m-%d"),
    postvax_any_test_date	= col_date(format="%Y-%m-%d"),
    postvax_any_emergency_date = col_date(format="%Y-%m-%d"),
    postvax_any_hospitalisation_date = col_date(format="%Y-%m-%d"),
    postvax_any_death_date = col_date(format="%Y-%m-%d"),
    
    # CKD groups
    creatinine = col_double(), 
    dialysis = col_logical(), 
    kidney_transplant = col_logical(), 
    chronic_kidney_disease_diagnostic = col_logical(), 
    chronic_kidney_disease_stages_3_5 = col_logical(), 
    ukrr_2019 = col_logical(), 
    ukrr_2019_mod = col_character(), 
    ukrr_2020 = col_logical(), 
    ukrr_2020_mod = col_character(), 
    creatinine_operator = col_character(), 
    age_creatinine = col_integer(),
    
    # Priority groups
    care_home_type =  col_character(),
    care_home_tpp = col_logical(),
    care_home_code = col_logical(),
    cev = col_logical(),
    hscworker = col_logical(),
    endoflife = col_logical(),
    housebound = col_logical(),
    
    # Clinical/demographic variables
    age = col_integer(),    
    age_index = col_integer(),
    age_august2021 = col_integer(),
    sex = col_character(),
    bmi = col_character(),
    ethnicity_6 = col_character(),
    ethnicity_6_sus = col_character(),
    imd = col_character(),
    stp = col_character(),
    region = col_character(),
    rural_urban = col_integer(),
    sev_obesity = col_logical(),
    asthma = col_logical(),
    asplenia = col_logical(),
    cancer = col_logical(),
    haem_cancer = col_logical(),
    chd = col_logical(),
    chronic_neuro_dis_inc_sig_learn_dis = col_logical(),
    chronic_resp_dis = col_logical(),
    cld = col_logical(),
    diabetes = col_logical(),
    learning_disability = col_logical(),
    sev_mental_ill = col_logical(),
    organ_transplant = col_logical(),
    non_kidney_transplant = col_logical(),
    
    # Other
    prevax_tests_conducted_any = col_double()
  ),
  na = character() # more stable to convert to missing later
)

## Parse NAs
data_extract <- data_extract %>%
  mutate(across(
    .cols = where(is.character),
    .fns = ~na_if(.x, "")
  )) %>%
  # Convert numerics and integers but not id variables to NAs if 0
  mutate(across(
    .cols = c(where(is.numeric), -ends_with("_id")), 
    .fns = ~na_if(.x, 0)
  )) %>%
  # Converts TRUE/FALSE to 1/0
  mutate(across(
      where(is.logical),
      ~.x*1L 
    )) %>%
  arrange(patient_id) 


### eGFR calculations: adapted from COVID-19-vaccine-breakthrough Feb 2022 branch
# https://github.com/opensafely/COVID-19-vaccine-breakthrough/blob/updates-feb/analysis/data_process.R

data_extract <- data_extract %>%
  mutate(
    ## Define variables needed for calculation
    creatinine = replace(creatinine, creatinine <20 | creatinine >3000, NA), # Set implausible creatinine values to missing
    SCR_adj = creatinine/88.4 # Divide by 88.4 (to convert umol/l to mg/dl)
  ) %>%
  rowwise() %>%
  mutate(
    min = case_when(sex == "M" ~ min(SCR_adj/0.9, 1, na.rm = F)^-0.411, 
                    sex == "F" ~ min(SCR_adj/0.7, 1, na.rm = F)^-0.329),
    max = case_when(sex == "M" ~ max(SCR_adj/0.9, 1, na.rm = F)^-1.209, 
                    sex == "F" ~ max(SCR_adj/0.7, 1, na.rm = F)^-1.209)) %>%
  ungroup() %>%
  mutate(
    egfr = (min*max*141)*(0.993^age_creatinine),
    egfr = case_when(sex == "F" ~ egfr*1.018, TRUE ~ egfr),
    
    # For purpose of inclusion/exclusion criteria, replace egfr NAs with arbitrary high value
    egfr = replace_na(egfr, 300),
    
    egfr_cat5 = cut(
      egfr,
      breaks = c(0, 15, 30, 45, 60, 5000),
      labels = c("Stage 5", "Stage 4", "Stage 3b", "Stage 3a", "No CKD"),
      right = FALSE
    ),
    # Replace NAs with 'No CKD'
    egfr_cat5 = replace_na(egfr_cat5, "No CKD"), 
  ) %>%
  # Drop extra variables
  select(-c(min, max, SCR_adj))

## Define UKRR groups
data_extract <- data_extract %>%
  mutate(
    # Add UKRR modality at index date - either 2020 modality, or 2019 modality if died between index date and end of 2020
    ukrr_index_mod = ifelse(((!is.na(death_date)) & death_date>=as_date("2020-12-01") & death_date<=as_date("2020-12-31")) |
                             ((!is.na(dereg_date)) & dereg_date>=as_date("2020-12-01") & dereg_date<=as_date("2020-12-31")), ukrr_2019_mod, ukrr_2020_mod),

    # 2019
    ukrr_2019_group = "None",
    ukrr_2019_group = ifelse((!is.na(ukrr_2019_mod)) & ukrr_2019_mod == "Tx", "Tx", ukrr_2019_group),
    ukrr_2019_group = ifelse((!is.na(ukrr_2019_mod)) & (ukrr_2019_mod == "ICHD" | ukrr_2019_mod == "HHD" | ukrr_2019_mod == "HD" | ukrr_2019_mod == "PD"), "Dialysis", ukrr_2019_group),
    
    # 2020
    ukrr_2020_group = "None",
    ukrr_2020_group = ifelse((!is.na(ukrr_2020_mod)) & ukrr_2020_mod == "Tx", "Tx", ukrr_2020_group),
    ukrr_2020_group = ifelse((!is.na(ukrr_2020_mod)) & (ukrr_2020_mod == "ICHD" | ukrr_2020_mod == "HHD" | ukrr_2020_mod == "HD" | ukrr_2020_mod == "PD"), "Dialysis", ukrr_2020_group),
    
    # Index (01-Dec-2020)
    ukrr_index_group = "None",
    ukrr_index_group = ifelse((!is.na(ukrr_index_mod)) & ukrr_index_mod == "Tx", "Tx", ukrr_index_group),
    ukrr_index_group = ifelse((!is.na(ukrr_index_mod)) & (ukrr_index_mod == "ICHD" | ukrr_index_mod == "HHD" | ukrr_index_mod == "HD" | ukrr_index_mod == "PD"), "Dialysis", ukrr_index_group),

    # Set modalities as 'None' instead of NA
    ukrr_2019_mod = ifelse(is.na(ukrr_2019_mod), "None", ukrr_2019_mod), 
    ukrr_2020_mod = ifelse(is.na(ukrr_2020_mod), "None", ukrr_2020_mod), 
    ukrr_index_mod = ifelse(is.na(ukrr_index_mod), "None", ukrr_index_mod), 

    # Flag issues with ambiguous creatinine entries - either no creatinine-associated age or creatinine-linked operator (impacting interpretation of numeric values)
    creatinine_date_issue = ifelse((!is.na(creatinine)) & is.na(age_creatinine),1,0),
    creatinine_operator_issue = ifelse((!is.na(creatinine_operator)) & creatinine_operator %in% c("~", ">=", ">", "<", "<="),1,0)
  )

## Print key processing metrics to log file
print("Number of creatinine records with no associated date")
print(sum(data_extract$creatinine_date_issue))
print("Number of creatinine records with operator")
print(sum(data_extract$creatinine_operator_issue))

print("Cross-tabulate 2020 modalities")
print(table(data_extract$ukrr_2020_mod))
print("Cross-tabulate index modalities")
print(table(data_extract$ukrr_index_mod))
print("Cross-tabulate 2019 vs 2020 modalities")
print(table(data_extract$ukrr_2019_mod, data_extract$ukrr_2020_mod))
print("Cross-tabulate index vs 2020 modalities")
print(table(data_extract$ukrr_index_mod, data_extract$ukrr_2020_mod))

print("Cross-tabulate 2020 groups")
print(table(data_extract$ukrr_2020_group))
print("Cross-tabulate index groups")
print(table(data_extract$ukrr_index_group))
print("Cross-tabulate 2019 vs 2020 groups")
print(table(data_extract$ukrr_2019_group, data_extract$ukrr_2020_group))
print("Cross-tabulate index vs 2020 groups")
print(table(data_extract$ukrr_index_group, data_extract$ukrr_2020_group))

print("Cross-tabulate UKRR at index (dialysis/Tx) vs primary care dialysis/Tx flag")
print(table((data_extract$ukrr_index_group=="Dialysis" | data_extract$ukrr_index_group=="Tx"), (data_extract$dialysis==1 | data_extract$kidney_transplant)))

print("Sum UKRR 2019 population with dereg/death between 01-Dec-2020 and 31-Dec-2020")
print(sum(data_extract$ukrr_2019==1 &
  (((!is.na(data_extract$death_date)) & data_extract$death_date>=as_date("2020-12-01") & data_extract$death_date<=as_date("2020-12-31")) |
            ((!is.na(data_extract$dereg_date)) & data_extract$dereg_date>=as_date("2020-12-01") & data_extract$dereg_date<=as_date("2020-12-31")))))

## Add derived variables
data_processed <- data_extract %>%
  mutate(
    # CKD inclusion criteria
    ckd_inclusion_egfr_ukrr_D_T_3to5_diagnostic = ifelse(egfr < 60 | ukrr_2020_group=="Tx" | ukrr_2020_group=="Dialysis" | dialysis==1 | kidney_transplant==1 | chronic_kidney_disease_stages_3_5==1 | chronic_kidney_disease_diagnostic==1, 1, 0),
    ckd_inclusion_egfr_ukrr_D_T_3to5 = ifelse(egfr < 60 | ukrr_2020_group=="Tx" | ukrr_2020_group=="Dialysis" | dialysis==1 | kidney_transplant==1 | chronic_kidney_disease_stages_3_5==1, 1, 0),
    ckd_inclusion_egfr_ukrr_D_T = ifelse(egfr < 60 | ukrr_2020_group=="Tx" | ukrr_2020_group=="Dialysis" | dialysis==1 | kidney_transplant==1, 1, 0),
    ckd_inclusion_egfr_ukrr = ifelse(egfr < 60 | ukrr_2020_group=="Tx" | ukrr_2020_group=="Dialysis", 1, 0),

    # CKD 6-levels
    ckd_6cat = "No CKD",
    ckd_6cat = ifelse(egfr_cat5 == "Stage 3a" & ukrr_2020_group=="None", "CKD3a", ckd_6cat),
    ckd_6cat = ifelse(egfr_cat5 == "Stage 3b" & ukrr_2020_group=="None", "CKD3b", ckd_6cat),
    ckd_6cat = ifelse(egfr_cat5 == "Stage 4" & ukrr_2020_group=="None", "CKD4", ckd_6cat),
    ckd_6cat = ifelse(egfr_cat5 == "Stage 5" & ukrr_2020_group=="None", "CKD5", ckd_6cat),
    ckd_6cat = ifelse(ukrr_2020_group == "Dialysis", "RRT (dialysis)", ckd_6cat),
    ckd_6cat = ifelse(ukrr_2020_group == "Tx", "RRT (Tx)", ckd_6cat),

    # CKD 5-levels (merging 4/5)
    ckd_5cat = ckd_6cat,
    ckd_5cat = ifelse(ckd_6cat == "CKD4" | ckd_6cat == "CKD5", "CKD4-5", ckd_5cat),
    
    # Flag individuals with mismatch between UKRR and primary care data
    rrt_mismatch = ifelse((ckd_5cat=="CKD3a" | ckd_5cat=="CKD3b" | ckd_5cat=="CKD4-5") & (dialysis==1 | kidney_transplant==1), 1, 0),
    
    # Age
    ageband = cut(
      age,
      breaks = c(16, 70, 80, Inf),
      labels = c("16-69", "70-79", "80+"),
      right = FALSE
    ),
    
    ageband2 = cut(
      age,
      breaks = c(16, 65, 70, 75, 80, Inf),
      labels = c("16-64", "65-69", "70-74", "75-79", "80+"),
      right = FALSE
    ),
    
    # Any care home flag
    care_home = ifelse(care_home_tpp==1 | care_home_code==1, 1, 0), 
    
    # Original JCVI priority groups: https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/1007737/Greenbook_chapter_14a_30July2021.pdf#page=15
    jcvi_group = fct_case_when(
      care_home==1 & age>=65 ~ "1 (65+ care home resident)",
      care_home==0 & (age>=80 | hscworker==1) ~ "2 (80+ or health/social care worker)",
      care_home==0 & (age>=75 & age<80) ~ "3 (75+)",
      care_home==0 & ((age>=70 & age<75) | (cev==1 & age>=16 & age<70)) ~ "4 (70+ or clinically extremely vulnerable)",
      care_home==0 & cev==0 & (age>=65 & age<70) ~ "5 (65+)",
      TRUE ~ "6 (16-65 and clinically vulnerable)" # Since CKD patients would be classified as clinically vulnerable, 6 is maximum JCVI group for this population
    ),
    
    # Updated JCVI priority groups: https://www.england.nhs.uk/coronavirus/wp-content/uploads/sites/52/2021/07/C1327-covid-19-vaccination-autumn-winter-phase-3-planning.pdf
    jcvi_group_update = fct_case_when(
      care_home==1 | hscworker==1  ~ "1 (care home resident or health/social care worker)",
      care_home==0 & age_august2021>=80 ~ "2 (80+)",
      care_home==0 & age_august2021>=75 ~ "3 (75+)",
      care_home==0 & (age_august2021>=70 | (cev & (age_august2021>=16))) ~ "4 (70+ or clinically extremely vulnerable)",
      care_home==0 & age_august2021>=65 ~ "5 (65+)",
      TRUE ~ "6 (16-65 and clinically vulnerable)" # Since CKD patients would be classified as clinically vulnerable, 6 is maximum JCVI group for this population
    ),
    
    # Sex
    sex = fct_case_when(
      sex == "F" ~ "Female",
      sex == "M" ~ "Male",
      TRUE ~ NA_character_
    ),
    
    # Obesity
    obesity = ifelse(bmi == "Not obese",0,1),
    mod_sev_obesity = ifelse(bmi == "Obese II (35-39.9)" | bmi == "Obese III (40+)", 1, 0),
    
    # Ethnicity
    ethnicity_filled = ifelse(is.na(ethnicity_6), ethnicity_6_sus, ethnicity_6),
    ethnicity = ifelse(is.na(ethnicity_filled), 6, ethnicity_filled),
    ethnicity = fct_case_when(
      ethnicity == "1" ~ "White",
      ethnicity == "4" ~ "Black",
      ethnicity == "3" ~ "South Asian",
      ethnicity == "2" ~ "Mixed",
      ethnicity == "5" ~ "Other",
      TRUE ~ NA_character_),
    
    # IMD
    imd = na_if(imd, "0"),
    imd = fct_case_when(
      imd == 1 ~ "1 most deprived",
      imd == 2 ~ "2",
      imd == 3 ~ "3",
      imd == 4 ~ "4",
      imd == 5 ~ "5 least deprived",
      TRUE ~ NA_character_
    ),
    
    # Region
    region = fct_collapse(
      region,
      `East of England` = "East",
      `London` = "London",
      `Midlands` = c("West Midlands", "East Midlands"),
      `North East and Yorkshire` = c("Yorkshire and The Humber", "North East"),
      `North West` = "North West",
      `South East` = "South East",
      `South West` = "South West"
    ),
    
    # Rurality
    rural_urban_group = fct_case_when(
      rural_urban %in% c(3,4) ~ "Urban city or town",
      rural_urban %in% c(1,2) ~ "Urban conurbation",
      rural_urban %in% c(5,6,7,8) ~ "Rural town or village",
      TRUE ~ "Unknown"
    ),
    
    # Immunosuppression
    immunosuppression = ifelse((!is.na(immunosuppression_diagnosis_date)) | (!is.na(immunosuppression_medication_date)), 1, 0),
    
    # Other transplant excluding kidney transplants
    other_transplant = ifelse((organ_transplant==1 | non_kidney_transplant==1) & kidney_transplant==0 & ukrr_2020_group!="Tx", 1, 0), 
    
    # Any respiratory disease
    any_resp_dis= ifelse(chronic_resp_dis==1 | asthma==1, 1, 0), 
    
    # Any cancer
    any_cancer = ifelse(cancer==1 | haem_cancer==1, 1, 0), 

    # CEV other
    rrt_2020 = ifelse(ukrr_2020_group=="Tx" | ukrr_2020_group=="Dialysis", 1, 0),
    any_comorb = pmax(rrt_2020, 
                      immunosuppression, sev_obesity, diabetes, any_resp_dis,
                      chd, cld, asplenia, any_cancer, other_transplant, chronic_neuro_dis_inc_sig_learn_dis, sev_mental_ill),
    cev_other = ifelse(cev==1 & any_comorb==0, 1, 0),
    
    # Multiple comorbidities (non-CKD-related) - 0, 1, or 2+
    multimorb =
      (sev_obesity) +
      (chd) +
      (diabetes) +
      (cld) +
      (any_resp_dis) +
      (chronic_neuro_dis_inc_sig_learn_dis),
    multimorb = cut(multimorb, breaks = c(0, 1, 2, Inf), labels=c("0", "1", "2+"), right=FALSE),
    
    # Any immunosuppression (transplant, cancer, haematologic cancer, asplenia)
    any_immunosuppression = ifelse(ukrr_2020_group=="Tx" | other_transplant==1 | immunosuppression==1 | haem_cancer==1 | asplenia==1, 1, 0), 
    
    # Prior COVID - index
    prior_covid_cat = as.numeric(!is.na(pmin(prior_positive_test_date, prior_primary_care_covid_case_date, prior_covid_emergency_date, prior_covid_hospitalisation_date, na.rm=TRUE))),
    
    # Prior COVID - boost index
    prior_covid_cat_boost = as.numeric(!is.na(pmin(prior_positive_test_date_boost, prior_primary_care_covid_case_date_boost, prior_covid_emergency_date_boost, prior_covid_hospitalisation_date_boost, na.rm=TRUE))),

    # COVID in window spanning 90 days pre dose 1
    prevax_covid_cat = as.numeric(!is.na(pmin(prevax_positive_test_date, prevax_primary_care_covid_case_date, prevax_covid_emergency_date, prevax_covid_hospitalisation_date, na.rm=TRUE))),

    # Number of tests in pre-vaccination period
    prevax_tests_conducted_any = ifelse(is.na(prevax_tests_conducted_any), 0, prevax_tests_conducted_any),
    prevax_tests_cat = cut(prevax_tests_conducted_any, breaks=c(0, 1, 2, 3, Inf), labels=c("0", "1", "2", "3+"), right=FALSE),
    
    # Non-COVID death date
    noncoviddeath_date = if_else(!is.na(death_date) & is.na(postvax_covid_death_date), death_date, as.Date(NA_character_)),
  ) %>%
  # Drop superseded variables
  select(-c(care_home_code, care_home_tpp, ethnicity_6, ethnicity_filled, 
            immunosuppression_diagnosis_date, immunosuppression_medication_date)) %>%
  droplevels() 

## Print number of RRT mismatches to log file
print("Number of individuals with CDK3-5 based on eGFR<60 and dialysis/Tx primary care flag but not in UKRR at index")
print(sum(data_processed$rrt_mismatch))

## Summarise CKD subgroup
print("Cross-tabulate CKD subgroup")
print(table(data_processed$ckd_5cat))


## Apply dummy data script if not running in the server
if(Sys.getenv("OPENSAFELY_BACKEND") %in% c("", "expectations")){
  source(here::here("analysis", "dummy_data.R"))
  # Add binary flag to signal whether dummy data processing steps have been applied
  data_processed$mock_data_flag=1
} else {
  data_processed$mock_data_flag=0
}

## Linearise the vaccination dates of each individual, then determine the product and index for each dose
data_vax <- local({
  
  data_vax_pfizer <- data_processed %>%
    select(patient_id, matches("pfizer\\_date\\_\\d+")) %>%
    pivot_longer(
      cols = -patient_id,
      values_to = "date",
      values_drop_na = TRUE
    ) %>%
    arrange(patient_id, date) %>%
    mutate(type = "pfizer") %>%
    select(-name)

  data_vax_az <- data_processed %>%
    select(patient_id, matches("az\\_date\\_\\d+")) %>%
    pivot_longer(
      cols = -patient_id,
      values_to = "date",
      values_drop_na = TRUE
    ) %>%
    arrange(patient_id, date) %>%
    mutate(type = "az") %>%
    select(-name)
  
  data_vax_moderna <- data_processed %>%
    select(patient_id, matches("moderna\\_date\\_\\d+")) %>%
    pivot_longer(
      cols = -patient_id,
      values_to = "date",
      values_drop_na = TRUE
    ) %>%
    arrange(patient_id, date) %>%
    mutate(type = "moderna") %>%
    select(-name)
  
  data_vax <- rbind(data_vax_pfizer, data_vax_az, data_vax_moderna) %>%
    arrange(patient_id, date) %>%
    group_by(patient_id) %>%
    mutate(
      vax_index=row_number(),
      n_vax=max(row_number())
    ) %>%
    ungroup()
  
  data_vax
  
})

## Pivot to wide-form table that lists dates and products for up to 5 sequential doses
data_vax_wide = data_vax %>%
  pivot_wider(
    id_cols= patient_id,
    names_from = c("vax_index"),
    values_from = c("date", "type"),
    names_glue = "covid_vax_derived_{vax_index}_{.value}"
  )

## Pick out vaccine count for each individual, then merge
data_vax_count = data_vax[!duplicated(data_vax$patient_id),c("patient_id", "n_vax")]
data_vax_wide = data_vax_wide %>% left_join(data_vax_count, by ="patient_id")

## Merge with full data-set
data_processed_updated <- data_processed %>%
  left_join(data_vax_wide, by ="patient_id") %>%
  mutate(
    vax1_type = covid_vax_derived_1_type,
    vax2_type = covid_vax_derived_2_type,
    vax3_type = covid_vax_derived_3_type,
    vax4_type = covid_vax_derived_4_type,
    vax5_type = covid_vax_derived_5_type,
    
    # Derive vaccine combinations for doses 2 and 3
    vax12_type = paste0(vax1_type, "-", vax2_type),
    vax123_type = paste0(vax12_type, "-", vax3_type),
    
    # Set combinations to NA if latter dose not administered
    vax12_type = ifelse(is.na(vax2_type), NA_character_, vax12_type),
    vax123_type = ifelse(is.na(vax3_type), NA_character_, vax123_type),
    
    # Set n_vax to 0 if none recorded
    n_vax = ifelse(is.na(n_vax), 0, n_vax),
    
    # Add new descriptor variables with formal product names
    vax1_type_descr = fct_case_when(
      vax1_type == "pfizer" ~ "BNT162b2",
      vax1_type == "az" ~ "ChAdOx1",
      vax1_type == "moderna" ~ "Moderna",
      TRUE ~ NA_character_
    ),
    vax2_type_descr = fct_case_when(
      vax2_type == "pfizer" ~ "BNT162b2",
      vax2_type == "az" ~ "ChAdOx1",
      vax2_type == "moderna" ~ "Moderna",
      TRUE ~ NA_character_
    ),
    vax3_type_descr = fct_case_when(
      vax3_type == "pfizer" ~ "BNT162b2",
      vax3_type == "az" ~ "ChAdOx1",
      vax3_type == "moderna" ~ "Moderna",
      TRUE ~ NA_character_
    ),
    vax4_type_descr = fct_case_when(
      vax4_type == "pfizer" ~ "BNT162b2",
      vax4_type == "az" ~ "ChAdOx1",
      vax4_type == "moderna" ~ "Moderna",
      TRUE ~ NA_character_
    ),
    vax5_type_descr = fct_case_when(
      vax5_type == "pfizer" ~ "BNT162b2",
      vax5_type == "az" ~ "ChAdOx1",
      vax5_type == "moderna" ~ "Moderna",
      TRUE ~ NA_character_
    ),
    
    # Add new date variables
    vax1_date = covid_vax_derived_1_date,
    vax2_date = covid_vax_derived_2_date,
    vax3_date = covid_vax_derived_3_date,
    vax4_date = covid_vax_derived_4_date, 
    vax5_date = covid_vax_derived_5_date, 
    
    # Calculate time between vaccinations
    tbv1_2 = as.numeric(vax2_date - vax1_date),
    tbv2_3 = as.numeric(vax3_date - vax2_date),
    tbv3_4 = as.numeric(vax4_date - vax3_date),
    tbv4_5 = as.numeric(vax5_date - vax4_date) 
  ) %>%
  # Remove unnecessary variables including overall and product-specific dates (now replaced by vax{N}_date)
  select(
    -starts_with(c("covid_vax_", "pfizer_date_", "az_date_", "moderna_date_"))
  )

## Set factor levels for derived variables
data_processed_updated <- data_processed_updated %>%
  mutate(
    ukrr_2019_mod = factor(ukrr_2019_mod, levels = c("None", "HD", "ICHD", "HHD", "PD", "Tx")),
    ukrr_2020_mod = factor(ukrr_2020_mod, levels = c("None", "HD", "ICHD", "HHD", "PD", "Tx")),
    ukrr_index_mod = factor(ukrr_index_mod, levels = c("None", "HD", "ICHD", "HHD", "PD", "Tx")),
    ukrr_2019_group = factor(ukrr_2019_group, levels = c("None", "Dialysis", "Tx")),
    ukrr_2020_group = factor(ukrr_2020_group, levels = c("None", "Dialysis", "Tx")),
    ukrr_index_group = factor(ukrr_index_group, levels = c("None", "Dialysis", "Tx")),
    ckd_6cat = factor(ckd_6cat, levels = c("No CKD", "CKD3a", "CKD3b", "CKD4", "CKD5", "RRT (dialysis)", "RRT (Tx)")),
    ckd_5cat = factor(ckd_5cat, levels = c("No CKD", "CKD3a", "CKD3b", "CKD4-5", "RRT (dialysis)", "RRT (Tx)")),
  )

## Save dataset
write_rds(data_processed_updated, here::here("output", "data", "data_processed.rds"), compress = "gz")
write_csv(data_processed_updated, here::here("output", "data", "data_processed.csv"))
