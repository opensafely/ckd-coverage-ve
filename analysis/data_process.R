######################################

# This script:
# - imports data extracted by the cohort extractor
# - combines ethnicity columns
# - standardises some variables (eg convert to factor) and derives some new ones
# - applies additional dummy data processing steps via dummy_data.R script if being run locally (skipped if being run on OpenSAFELY server)
# - saves processed dataset

######################################


# Preliminaries ----

## packages
library('tidyverse')
library('lubridate')
library('here')

## Import custom user functions and packages
source(here::here("analysis", "functions.R"))

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
    pfizer_date_1 = col_date(format="%Y-%m-%d"),
    pfizer_date_2 = col_date(format="%Y-%m-%d"),
    pfizer_date_3 = col_date(format="%Y-%m-%d"),
    pfizer_date_4 = col_date(format="%Y-%m-%d"),
    az_date_1 = col_date(format="%Y-%m-%d"),
    az_date_2 = col_date(format="%Y-%m-%d"),
    az_date_3 = col_date(format="%Y-%m-%d"),
    az_date_4 = col_date(format="%Y-%m-%d"),
    moderna_date_1 = col_date(format="%Y-%m-%d"),
    moderna_date_2 = col_date(format="%Y-%m-%d"),
    moderna_date_3 = col_date(format="%Y-%m-%d"),
    moderna_date_4 = col_date(format="%Y-%m-%d"),
    
    # Dates for clinical variables
    creatinine_date = col_date(format="%Y-%m-%d"),
    death_date = col_date(format="%Y-%m-%d"),
    dereg_date = col_date(format="%Y-%m-%d"),
    bp_sys_date_measured = col_date(format="%Y-%m-%d"),
    bp_dias_date_measured = col_date(format="%Y-%m-%d"),
    immunosuppression_diagnosis_date = col_date(format="%Y-%m-%d"),
    immunosuppression_medication_date = col_date(format="%Y-%m-%d"),
    
    # Dates for Covid-related variables
    prior_positive_test_date = col_date(format="%Y-%m-%d"),
    prior_primary_care_covid_case_date = col_date(format="%Y-%m-%d"),
    prior_covid_hospitalisation_date = col_date(format="%Y-%m-%d"),
    prevax_positive_test_date	= col_date(format="%Y-%m-%d"),
    prevax_primary_care_covid_case_date	= col_date(format="%Y-%m-%d"),
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
    sex = col_character(),
    bmi = col_character(),
    smoking_status =  col_character(),
    ethnicity_6 = col_character(),
    ethnicity_6_sus = col_character(),
    imd = col_character(),
    region = col_character(),
    rural_urban = col_integer(),
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
  # convert numerics and integers but not id variables to NAs if 0
  mutate(across(
    .cols = c(where(is.numeric), -ends_with("_id")), 
    .fns = ~na_if(.x, 0)
  )) %>%
  # converts TRUE/FALSE to 1/0
  mutate(across(
      where(is.logical),
      ~.x*1L 
    )) %>%
  arrange(patient_id) 


### eGFR calculations: adapted from cr_create_analysis_dataset.do in risk factor analysis repo
# https://github.com/opensafely/risk-factors-research/

data_extract <- data_extract %>%
  mutate(
    ## define variables needed for calculation
    creatinine = replace(creatinine, creatinine <20 | creatinine >3000, NA), # set implausible creatinine values to missing
    SCR_adj = creatinine/88.4 # divide by 88.4 (to convert umol/l to mg/dl)
  ) %>%
  rowwise() %>%
  mutate(
    min = case_when(sex == "M" ~ min((SCR_adj/0.9), 1, na.rm = F)^-0.411, 
                    sex == "F" ~ min(SCR_adj/0.7, 1, na.rm = F)^-0.329),
    max = case_when(sex == "M" ~ max(SCR_adj/0.9, 1, na.rm = F)^-1.209, 
                    sex == "F" ~ max(SCR_adj/0.7, 1, na.rm = F)^-1.209)) %>%
  ungroup() %>%
  mutate(
    egfr = (min*max*141)*(0.993^age),
    egfr = case_when(sex == "F" ~ egfr*1.018, TRUE ~ egfr),
    
    egfr_cat5 = cut(
      egfr,
      breaks = c(0, 15, 30, 45, 60, 5000),
      labels = c("Stage 5", "Stage 4", "Stage 3b", "Stage 3a", "No CKD"),
      right = FALSE
    ),
    
    egfr_cat3 = cut(
      egfr,
      breaks = c(0, 30, 60, 5000),
      labels = c("Stage 4/5 eGFR<30", "Stage 3a/3b eGFR 30-60", "No CKD"),
      right = FALSE
    ) 
  ) %>%
  # Drop extra variables
  select(-c(min, max, SCR_adj))

# Replace NAs with 'No CKD' in new categories
data_extract$egfr_cat5[is.na(data_extract$egfr_cat5)] = "No CKD"
data_extract$egfr_cat3[is.na(data_extract$egfr_cat3)] = "No CKD"

# Set NA egfr readings to arbitrary value of 300 to enable inclusion criteria to be applied
data_extract$egfr[is.na(data_extract$egfr)] = 300

## Format columns (i.e, set factor levels)
data_processed <- data_extract %>%
  mutate(
    # CKD inclusion criteria
    ckd_inclusion_any = ifelse(egfr < 60 | dialysis==1 | kidney_transplant==1 | chronic_kidney_disease_diagnostic==1 | chronic_kidney_disease_stages_3_5==1, 1, 0),
    ckd_inclusion_strict = ifelse(egfr < 60 | dialysis==1 | kidney_transplant==1, 1, 0),
    
    # CKD 7-levels
    ckd_7cat = ifelse(egfr_cat5 == "Stage 3a" & dialysis == 0 & kidney_transplant== 0, "CKD3a (D-T-)", NA),
    ckd_7cat = ifelse(egfr_cat5 == "Stage 3b" & dialysis == 0 & kidney_transplant== 0, "CKD3b (D-T-)", ckd_7cat),
    ckd_7cat = ifelse(egfr_cat5 == "Stage 4" & dialysis == 0 & kidney_transplant== 0, "CKD4 (D-T-)", ckd_7cat),
    ckd_7cat = ifelse(egfr_cat5 == "Stage 5" & dialysis == 0 & kidney_transplant== 0, "CKD5 (D-T-)", ckd_7cat),
    ckd_7cat = ifelse(dialysis == 1 & kidney_transplant == 0, "CKD (D+T-)", ckd_7cat),
    ckd_7cat = ifelse(dialysis == 0 & kidney_transplant == 1, "CKD (D-T+)", ckd_7cat),
    ckd_7cat = ifelse(dialysis == 1 & kidney_transplant == 1, "CKD (D+T+)", ckd_7cat),
    ckd_7cat = ifelse(is.na(ckd_7cat), "Unknown", ckd_7cat),
    
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
    
    # any care home flag
    care_home = ifelse(care_home_tpp==1 | care_home_code==1, 1, 0), 
    
    # https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/1007737/Greenbook_chapter_14a_30July2021.pdf#page=12
    jcvi_group = fct_case_when(
      care_home==1 ~ "1 (care home residents/staff)",
      care_home==0 & (age>=80 | hscworker==1) ~ "2 (80+ or health/social care worker)",
      care_home==0 & (age>=75 & age<80) ~ "3 (75+)",
      care_home==0 & ((age>=70 & age<75) | (cev==1 & age>=16 & age<70)) ~ "4 (70+ or clinically extremely vulnerable)",
      care_home==0 & cev==0 & (age>=65 & age<70) ~ "5 (65+)",
      TRUE ~ "6 (16-65 and clinically vulnerable)" # Since CKD patients would be classified as clinically vulnerable, 6 is maximum JCVI group for this population
    ),
    
    # Sex
    sex = fct_case_when(
      sex == "F" ~ "Female",
      sex == "M" ~ "Male",
      TRUE ~ NA_character_
    ),
    
    # obesity
    obesity = ifelse(bmi == "Not obese",0,1),
    mod_sev_obesity = ifelse(bmi == "Obese II (35-39.9)" | bmi == "Obese III (40+)", 1, 0),
    
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
      TRUE ~ NA_character_),
    
    # Rurality
    rural_urban_group = fct_case_when(
      rural_urban %in% c(1,2) ~ "Urban conurbation",
      rural_urban %in% c(3,4) ~ "Urban city or town",
      rural_urban %in% c(5,6,7,8) ~ "Rural town or village",
      TRUE ~ "Unknown"
    ),
    
    ## Blood pressure
    bpcat = ifelse(bp_sys < 120 & bp_dias < 80, 1, NA),
    bpcat = ifelse(bp_sys >= 120 & bp_sys < 130 & bp_dias < 80, 2, bpcat),
    bpcat = ifelse(bp_sys >= 130 & bp_dias >= 90, 3, bpcat),
    bpcat = ifelse(is.na(bpcat), 4, bpcat),
    
    bpcat = fct_case_when(
      bpcat == 1 ~ "Normal",
      bpcat == 2 ~ "Elevated",
      bpcat == 3 ~ "High",
      bpcat == 4 ~ "Unknown",
      TRUE ~ NA_character_
    ),
    
    # Immunosuppression
    immunosuppression = ifelse(!is.na(immunosuppression_diagnosis_date) | !is.na(immunosuppression_medication_date), 1, 0),
    
    # Multiple morbidities (non-CKD) - 0, 1, or 2+
    multimorb =
      (mod_sev_obesity) +
      (diabetes) +
      (chronic_resp_dis | asthma) +
      (chd) +
      (cld) +
      (asplenia),
    multimorb = cut(multimorb, breaks = c(0, 1, 2, Inf), labels=c("0", "1", "2+"), right=FALSE),
    
    # Any respiratory disease
    any_resp_dis= ifelse(chronic_resp_dis==1 | asthma==1, 1, 0), 
    
    # Any cancer
    any_cancer = ifelse(cancer==1 | haem_cancer==1, 1, 0), 

    # CEV other
    any_comorb = pmax(dialysis, kidney_transplant, immunosuppression, mod_sev_obesity, diabetes, chronic_resp_dis,
                      chd, cld, asplenia, cancer, haem_cancer, non_kidney_transplant, chronic_neuro_dis_inc_sig_learn_dis, sev_mental_ill),
    cev_other = ifelse(cev==1 & any_comorb==0, 1, 0),
    
    # Prior covid
    prior_covid_cat = !is.na(pmin(prior_positive_test_date, prior_primary_care_covid_case_date, prior_covid_hospitalisation_date, na.rm=TRUE)),
    
    # Covid in window spanning 90 days pre dose 1 up to dose 2
    prevax_covid_cat = !is.na(pmin(prevax_positive_test_date, prevax_primary_care_covid_case_date, prevax_covid_hospitalisation_date, na.rm=TRUE)),

    # Number of tests in pre-vaccination period
    prevax_tests_conducted_any = ifelse(is.na(prevax_tests_conducted_any), 0, prevax_tests_conducted_any),
    prevax_tests_cat = cut(prevax_tests_conducted_any, breaks=c(0, 1, 2, 3, Inf), labels=c("0", "1", "2", "3+"), right=FALSE)
  ) %>%
  # Drop superseded variables
  select(-c(care_home_code, care_home_tpp, ethnicity_6, ethnicity_filled, 
            immunosuppression_diagnosis_date, immunosuppression_medication_date)) %>%
  droplevels() 

# apply dummy data script if not running in the server
#Sys.getenv()
#Sys.getenv("USER")
#Sys.getenv("OPENSAFELY_BACKEND")

if(Sys.getenv("OPENSAFELY_BACKEND") %in% c("", "expectations")){
  source(here::here("analysis", "dummy_data.R"))
  # add binary flag to signal whether dummy data processing steps have been applied
  data_processed$mock_data_flag=1
} else {
  data_processed$mock_data_flag=0
}

# linearise the vaccination dates of each individual, then determine the product and index for each dose
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

# pivot to wide-form table that lists dates and products for up to 4 sequential doses
data_vax_wide = data_vax %>%
  pivot_wider(
    id_cols= patient_id,
    names_from = c("vax_index"),
    values_from = c("date", "type"),
    names_glue = "covid_vax_{vax_index}_{.value}"
  )

# pick out vaccine count for each individual, then merge
data_vax_count = data_vax[!duplicated(data_vax$patient_id),c("patient_id", "n_vax")]
data_vax_wide = data_vax_wide %>% left_join(data_vax_count, by ="patient_id")

# merge with full data-set
data_processed_updated <- data_processed %>%
  left_join(data_vax_wide, by ="patient_id") %>%
  mutate(
    vax1_type = covid_vax_1_type,
    vax2_type = covid_vax_2_type,
    vax3_type = covid_vax_3_type,
    vax4_type = covid_vax_4_type,
    
    vax12_type = paste0(vax1_type, "-", vax2_type),
    vax123_type = paste0(vax12_type, "-", vax3_type),

    # add new descriptor variables with formal product names
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
    
    vax1_date = covid_vax_1_date,
    vax2_date = covid_vax_2_date,
    vax3_date = covid_vax_3_date,
    vax4_date = covid_vax_4_date, 
    
    # Calculate time between vaccinations
    tbv1_2 = as.numeric(vax2_date - vax1_date),
    tbv2_3 = as.numeric(vax3_date - vax2_date),
    tbv3_4 = as.numeric(vax4_date - vax3_date) 

  ) %>%
  # Remove unneccessary variables including overall and product-specific dates (now replaced by vax{N}_date)
  select(
    -starts_with(c("covid_vax_", "pfizer_date_", "az_date_", "moderna_date_"))
  )

# remove type combinations if latter vaccine not administered
data_processed_updated$vax12_type[is.na(data_processed_updated$vax2_type)] = NA_character_
data_processed_updated$vax123_type[is.na(data_processed_updated$vax3_type)] = NA_character_

# set n_vax to 0 if none recorded
data_processed_updated$n_vax[is.na(data_processed_updated$n_vax)]=0

# Save dataset as .rds files ----
write_rds(data_processed_updated, here::here("output", "data", "data_processed.rds"), compress = "gz")
write_csv(data_processed_updated, here::here("output", "data", "data_processed.csv"))
