######################################
# This script provides the formal specification of the study data that will be extracted from the OpenSAFELY database
######################################

# IMPORT STATEMENTS ----

## Import code building blocks from cohort extractor package
from cohortextractor import (
    StudyDefinition,
    patients,
    codelist_from_csv,
    codelist,
    filter_codes_by_category,
    combine_codelists
)

## Import codelists from codelist.py (which pulls them from the codelist folder)
from codelists import *

### Set initial date parameters
start_date = "2020-12-01"
end_date = "2022-04-20" # date of most recent coverage report update (https://reports.opensafely.org/reports/vaccine-coverage/)
# if end_date updated, also update return_expectations in defined variables for consistency

study = StudyDefinition(
    # Configure the expectations framework
    default_expectations={
        "date": {"earliest": "1970-01-01", "latest": end_date},
        "rate": "uniform",
        "incidence": 0.2,
    },

    # Set index date to start date
    index_date = start_date,

    # Define study population
    population=patients.satisfying(
        """
        registered
        AND
        NOT has_died
        AND
        (age >= 16 AND age < 120)
        AND
        (creatinine>0 OR dialysis OR kidney_transplant OR chronic_kidney_disease_diagnostic OR chronic_kidney_disease_stages_3_5 OR
        ukrr_2019 OR ukrr_2020)
        """,
        # registered before vaccine campaign commenced
        registered = patients.registered_as_of(
            "index_date - 3 months",
        ),
        has_died = patients.died_from_any_cause(
            on_or_before = "index_date - 1 day",
            returning = "binary_flag",
        ),
    ),

###############################################################################
# COVID VACCINATION - ANY TYPE
###############################################################################

    # date of first COVID vaccination - source nhs-covid-vaccination-coverage
    covid_vax_date_1=patients.with_tpp_vaccination_record(
        target_disease_matches="SARS-2 CORONAVIRUS",
        between=["2020-12-01",end_date],  # any dose recorded after 01/12/2020
        find_first_match_in_period=True,
        returning="date",
        date_format="YYYY-MM-DD"
    ),
    
    # date of second COVID vaccination - source nhs-covid-vaccination-coverage
    covid_vax_date_2=patients.with_tpp_vaccination_record(
        target_disease_matches="SARS-2 CORONAVIRUS",
        between=["covid_vax_date_1 + 1 day",end_date], # from day after previous dose
        find_first_match_in_period=True,
        returning="date",
        date_format="YYYY-MM-DD"
    ),
    
    # date of third COVID vaccination (primary or booster) - modified from nhs-covid-vaccination-coverage
    # 01 Sep 2021: 3rd dose (primary) at interval of >=8w recommended for immunosuppressed
    # 14 Sep 2021: 3rd dose (booster) reommended for groups 1-9 at >=6m
    # 15 Nov 2021: 3rd dose (booster) recommended for 40–49y at >=6m
    # 29 Nov 2021: 3rd dose (booster) recommended for 18–39y at >=3m
    covid_vax_date_3=patients.with_tpp_vaccination_record(
        target_disease_matches="SARS-2 CORONAVIRUS",
        between=["covid_vax_date_2 + 1 day",end_date], # from day after previous dose
        find_first_match_in_period=True,
        returning="date",
        date_format="YYYY-MM-DD"
    ),
    
    # date of fourth COVID vaccination (booster) - newly added
    # booster for immunosuppressed individuals who recieved a 3-dose primary series
    # recommended at interval of 3 months
    # to capture above, include fourth doses from 77 days (84 - 7) to allow for early receipt within 12th week
    covid_vax_date_4=patients.with_tpp_vaccination_record(
        target_disease_matches="SARS-2 CORONAVIRUS",
        between=["covid_vax_date_3 + 1 day",end_date], # from day after previous dose
        find_first_match_in_period=True,
        returning="date",
        date_format="YYYY-MM-DD"
    ),
    
    # date of fifth COVID vaccination (extended primary series + two boosters
    covid_vax_date_5=patients.with_tpp_vaccination_record(
        target_disease_matches="SARS-2 CORONAVIRUS",
        between=["covid_vax_date_4 + 1 day",end_date], # from day after previous dose
        find_first_match_in_period=True,
        returning="date",
        date_format="YYYY-MM-DD"
    ),


###############################################################################
# COVID VACCINATION - pfizer/biontech
###############################################################################

    # date of first COVID vaccination - pfizer
    pfizer_date_1=patients.with_tpp_vaccination_record(
        product_name_matches="COVID-19 mRNA Vaccine Comirnaty 30micrograms/0.3ml dose conc for susp for inj MDV (Pfizer)",
        between=["2020-12-01",end_date], # any dose recorded after 01/12/2020
        find_first_match_in_period=True,
        returning="date",
        date_format="YYYY-MM-DD"
    ),
    
    # date of second COVID vaccination - pfizer
    pfizer_date_2=patients.with_tpp_vaccination_record(
        product_name_matches="COVID-19 mRNA Vaccine Comirnaty 30micrograms/0.3ml dose conc for susp for inj MDV (Pfizer)",
        between=["pfizer_date_1 + 1 day",end_date],  # from day after previous dose
        find_first_match_in_period=True,
        returning="date",
        date_format="YYYY-MM-DD"
    ),
    
    # date of third COVID vaccination - pfizer
    pfizer_date_3=patients.with_tpp_vaccination_record(
        product_name_matches="COVID-19 mRNA Vaccine Comirnaty 30micrograms/0.3ml dose conc for susp for inj MDV (Pfizer)",
        between=["pfizer_date_2 + 1 day",end_date],   # from day after previous dose
        find_first_match_in_period=True,
        returning="date",
        date_format="YYYY-MM-DD"
    ),
    
    # date of fourth COVID vaccination - pfizer
    pfizer_date_4=patients.with_tpp_vaccination_record(
        product_name_matches="COVID-19 mRNA Vaccine Comirnaty 30micrograms/0.3ml dose conc for susp for inj MDV (Pfizer)",
        between=["pfizer_date_3 + 1 day",end_date], # from day after previous dose
        find_first_match_in_period=True,
        returning="date",
        date_format="YYYY-MM-DD"
    ),
    
    # date of fifth COVID vaccination - pfizer
    pfizer_date_5=patients.with_tpp_vaccination_record(
        product_name_matches="COVID-19 mRNA Vaccine Comirnaty 30micrograms/0.3ml dose conc for susp for inj MDV (Pfizer)",
        between=["pfizer_date_4 + 1 day",end_date], # from day after previous dose
        find_first_match_in_period=True,
        returning="date",
        date_format="YYYY-MM-DD"
    ),
    
###############################################################################
# COVID VACCINATION - astrazeneca
###############################################################################

    # date of first COVID vaccination - astrazeneca
    az_date_1=patients.with_tpp_vaccination_record(
        product_name_matches="COVID-19 Vac AstraZeneca (ChAdOx1 S recomb) 5x10000000000 viral particles/0.5ml dose sol for inj MDV",
        between=["2020-12-01",end_date], # any dose recorded after 01/12/2020
        find_first_match_in_period=True,
        returning="date",
        date_format="YYYY-MM-DD"
    ),
    
    # date of second COVID vaccination - astrazeneca
    az_date_2=patients.with_tpp_vaccination_record(
        product_name_matches="COVID-19 Vac AstraZeneca (ChAdOx1 S recomb) 5x10000000000 viral particles/0.5ml dose sol for inj MDV",
        between=["az_date_1 + 1 day",end_date], # from day after previous dose
        find_first_match_in_period=True,
        returning="date",
        date_format="YYYY-MM-DD"
    ),
    
    # date of third COVID vaccination - astrazeneca
    az_date_3=patients.with_tpp_vaccination_record(
        product_name_matches="COVID-19 Vac AstraZeneca (ChAdOx1 S recomb) 5x10000000000 viral particles/0.5ml dose sol for inj MDV",
        between=["az_date_2 + 1 day",end_date], # from day after previous dose
        find_first_match_in_period=True,
        returning="date",
        date_format="YYYY-MM-DD"
    ),
    
    # date of fourth COVID vaccination - astrazeneca
    az_date_4=patients.with_tpp_vaccination_record(
        product_name_matches="COVID-19 Vac AstraZeneca (ChAdOx1 S recomb) 5x10000000000 viral particles/0.5ml dose sol for inj MDV",
        between=["az_date_3 + 1 day",end_date], # from day after previous dose
        find_first_match_in_period=True,
        returning="date",
        date_format="YYYY-MM-DD"
    ),
    
    # date of fourth COVID vaccination - astrazeneca
    az_date_5=patients.with_tpp_vaccination_record(
        product_name_matches="COVID-19 Vac AstraZeneca (ChAdOx1 S recomb) 5x10000000000 viral particles/0.5ml dose sol for inj MDV",
        between=["az_date_4 + 1 day",end_date], # from day after previous dose
        find_first_match_in_period=True,
        returning="date",
        date_format="YYYY-MM-DD"
    ),
    
###############################################################################
# COVID VACCINATION - moderna
###############################################################################

    # date of first COVID vaccination - moderna
    moderna_date_1=patients.with_tpp_vaccination_record(
        product_name_matches="COVID-19 mRNA Vaccine Spikevax (nucleoside modified) 0.1mg/0.5mL dose disp for inj MDV (Moderna)",
        between=["2020-12-01",end_date], # any dose recorded after 01/12/2020
        find_first_match_in_period=True,
        returning="date",
        date_format="YYYY-MM-DD"
    ),
    
    # date of second COVID vaccination - moderna
    moderna_date_2=patients.with_tpp_vaccination_record(
        product_name_matches="COVID-19 mRNA Vaccine Spikevax (nucleoside modified) 0.1mg/0.5mL dose disp for inj MDV (Moderna)",
        between=["moderna_date_1 + 1 day",end_date], # from day after previous dose
        find_first_match_in_period=True,
        returning="date",
        date_format="YYYY-MM-DD"
    ),
    
    # date of third COVID vaccination - moderna
    moderna_date_3=patients.with_tpp_vaccination_record(
        product_name_matches="COVID-19 mRNA Vaccine Spikevax (nucleoside modified) 0.1mg/0.5mL dose disp for inj MDV (Moderna)",
        between=["moderna_date_2 + 1 day",end_date], # from day after previous dose
        find_first_match_in_period=True,
        returning="date",
        date_format="YYYY-MM-DD"
    ),
    
    # date of fourth COVID vaccination - moderna
    moderna_date_4=patients.with_tpp_vaccination_record(
        product_name_matches="COVID-19 mRNA Vaccine Spikevax (nucleoside modified) 0.1mg/0.5mL dose disp for inj MDV (Moderna)",
        between=["moderna_date_3 + 1 day",end_date], # from day after previous dose
        find_first_match_in_period=True,
        returning="date",
        date_format="YYYY-MM-DD"
    ),
    
    # date of fifth COVID vaccination - moderna
    moderna_date_5=patients.with_tpp_vaccination_record(
        product_name_matches="COVID-19 mRNA Vaccine Spikevax (nucleoside modified) 0.1mg/0.5mL dose disp for inj MDV (Moderna)",
        between=["moderna_date_4 + 1 day",end_date], # from day after previous dose
        find_first_match_in_period=True,
        returning="date",
        date_format="YYYY-MM-DD"
    ),

###############################################################################
# CKD DEFINITIONS - adapted from https://github.com/opensafely/risk-factors-research
###############################################################################

  ## Creatinine level for eGFR calculation
  # https://github.com/ebmdatalab/tpp-sql-notebook/issues/17
  creatinine = patients.with_these_clinical_events(
    creatinine_codes,
    find_last_match_in_period=True,
    between=["index_date - 2 years","index_date - 1 day"],
    returning="numeric_value",
    #include_date_of_match=True,
    #date_format = "YYYY-MM-DD", #formerly include_month=True,
    return_expectations={
        "float": {"distribution": "normal", "mean": 60.0, "stddev": 15},
        "incidence": 0.95,
    },
  ),
  
  ## Determine date of creatinine test
  creatinine_date = patients.with_these_clinical_events(
    creatinine_codes,
    find_last_match_in_period=True,
    between=["index_date - 2 years","index_date - 1 day"],
    returning="date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "index_date - 2 years", "latest" : "index_date - 1 day"},
      "rate": "uniform",
      "incidence": 0.95,
    },
  ),
  
  ## Extract any operators associated with creatinine readings
  creatinine_operator = patients.comparator_from(
    "creatinine",
     return_expectations={
       "rate": "universal",
        "category": {
         "ratios": {  # ~, =, >= , > , < , <=
            None: 0.10,
            "~": 0.05,
            "=": 0.65,
            ">=": 0.05,
            ">": 0.05,
            "<": 0.05,
            "<=": 0.05,
            }
        },
        "incidence": 0.80,
      },
    ),
    
  ## Age at creatinine test
  age_creatinine = patients.age_as_of(
    "creatinine_date",
       return_expectations = {
      "rate": "universal",
      "int": {"distribution": "population_ages"},
    },
  ),
  
  ## CKD - dialysis
  dialysis = patients.with_these_clinical_events(
    dialysis_codes,
    returning = "binary_flag",
    find_last_match_in_period = True,
    on_or_before = "index_date - 1 day"
  ),
  
  ## CKD - kidney transplant
  kidney_transplant = patients.with_these_clinical_events(
    kidney_transplant_codes,
    returning = "binary_flag",
    find_last_match_in_period = True,
    on_or_before = "index_date - 1 day"
  ),
  
  ## CKD - diagnostic codes
  chronic_kidney_disease_diagnostic = patients.with_these_clinical_events(
    chronic_kidney_disease_diagnostic_codes,
    returning = "binary_flag",
    find_last_match_in_period = True,
    on_or_before = "index_date - 1 day",
  ),
  
  ## CKD codes - all stages
  # chronic_kidney_disease_all_stages = patients.with_these_clinical_events(
  #   chronic_kidney_disease_all_stages_codes,
  #   returning = "date",
  #   find_last_match_in_period = True,
  #   on_or_before = "index_date - 1 day",
  #   date_format = "YYYY-MM-DD",
  # ),
  
  ## CKD codes - stages 3 - 5
  chronic_kidney_disease_stages_3_5 = patients.with_these_clinical_events(
    chronic_kidney_disease_stages_3_5_codes,
    returning = "binary_flag",
    find_last_match_in_period = True,
    on_or_before = "index_date - 1 day",
  ),

###############################################################################
# UKRR VARIABLES - adapted from https://github.com/opensafely/renal-short-data-report
###############################################################################

  ## Present in UKRR cohort at end of 2019
  ukrr_2019 = patients.with_record_in_ukrr(
    from_dataset="2019_prevalence",
    returning="binary_flag",
    return_expectations={
      "incidence": 0.25
      },
    ),
  
  ## Modality at end of 2019
  ukrr_2019_mod = patients.with_record_in_ukrr(
    from_dataset="2019_prevalence",
    returning="treatment_modality_prevalence",
    return_expectations={
      "category": {"ratios": {"ICHD": 0.2, "HHD": 0.1, "HD": 0.1, "PD": 0.1, "Tx": 0.5}},
      "incidence": 0.25,
      },
    ),
    
  ## Present in UKRR cohort at end of 2020
  ukrr_2020 = patients.with_record_in_ukrr(
    from_dataset="2020_prevalence",
    returning="binary_flag",
    return_expectations={
            "incidence": 0.25
        },
    ),
    
  ## Modality at end of 2020
  ukrr_2020_mod = patients.with_record_in_ukrr(
    from_dataset="2020_prevalence",
    returning="treatment_modality_prevalence",
    return_expectations={
      "category": {"ratios": {"ICHD": 0.2, "HHD": 0.1, "HD": 0.1, "PD": 0.1, "Tx": 0.5}},
      "incidence": 0.25,
      },
    ),

###############################################################################
# CENSORING VARIABLES
###############################################################################

  ## Death
  death_date = patients.died_from_any_cause(
    between=["index_date",end_date],
    returning = "date_of_death",
    date_format = "YYYY-MM-DD",
  ),
  
  ## De-registration
  dereg_date = patients.date_deregistered_from_all_supported_practices(
    between=["index_date",end_date],
    date_format = "YYYY-MM-DD",
  ),
  
###############################################################################
# PRIORITY GROUPS
###############################################################################

  ## Care home status
  care_home_type=patients.care_home_status_as_of(
      "index_date - 1 day",
      categorised_as={
          "Carehome": """
            IsPotentialCareHome
            AND LocationDoesNotRequireNursing='Y'
            AND LocationRequiresNursing='N'
          """,
          "Nursinghome": """
            IsPotentialCareHome
            AND LocationDoesNotRequireNursing='N'
            AND LocationRequiresNursing='Y'
          """,
          "Mixed": "IsPotentialCareHome",
          "": "DEFAULT",  # use empty string
      },
      return_expectations={
          "category": {"ratios": {"Carehome": 0.05, "Nursinghome": 0.05, "Mixed": 0.05, "": 0.85, }, },
          "incidence": 1,
      },
  ),
  
  ## Simple care home flag
  care_home_tpp=patients.satisfying(
      """care_home_type""",
      return_expectations={"incidence": 0.01},
  ),
  
  ## Care home
  care_home_code = patients.with_these_clinical_events(
    carehome_primis_codes,
    on_or_before = "index_date - 1 day",
    returning="binary_flag",
  ),
  
  ## PRIMIS overall flag for clinically extremely vulnerable group
  cev = patients.satisfying(
    """
      severely_clinically_vulnerable
      AND
      NOT less_vulnerable
      """,
    return_expectations = {
      "incidence": 0.3,
    },
    
    ### SHIELDED GROUP - first flag all patients with "high risk" codes
    severely_clinically_vulnerable = patients.with_these_clinical_events(
      high_risk_codes,
      returning="binary_flag",
      on_or_before = "index_date - 1 day",
      find_last_match_in_period = True,
      return_expectations = {"incidence": 0.02,},
    ),
    
    # find date at which the high risk code was added
    date_severely_clinically_vulnerable = patients.date_of(
      "severely_clinically_vulnerable",
      date_format = "YYYY-MM-DD",
    ),
    
    ### NOT SHIELDED GROUP (medium and low risk) - only flag if later than 'shielded'
    less_vulnerable = patients.with_these_clinical_events(
      not_high_risk_codes,
      between=["date_severely_clinically_vulnerable + 1 day", "index_date - 1 day",],
      return_expectations = {"incidence": 0.01,},
    ),
  ),
  
  ## HCW
  hscworker = patients.with_healthcare_worker_flag_on_covid_vaccine_record(returning = "binary_flag"),
  
  ## End of life care
  endoflife = patients.satisfying(
    """
    midazolam OR
    endoflife_coding
    """,
  
    midazolam = patients.with_these_medications(
      midazolam,
      returning="binary_flag",
      on_or_before = "index_date - 1 day",
    ),
    
    endoflife_coding = patients.with_these_clinical_events(
      eol,
      returning="binary_flag",
      on_or_before = "index_date - 1 day",
      find_last_match_in_period = True,
    ),
        
  ),
  
  ## Housebound
  housebound = patients.satisfying(
    """housebound_date
    AND NOT no_longer_housebound
    AND NOT moved_into_care_home
    """,
        
    housebound_date=patients.with_these_clinical_events(
      housebound,
      on_or_before="index_date - 1 day",
      find_last_match_in_period = True,
      returning = "date",
      date_format = "YYYY-MM-DD",
    ),
    
    no_longer_housebound=patients.with_these_clinical_events(
      no_longer_housebound,
      between = ["housebound_date","index_date - 1 day"]
    ),
    
    moved_into_care_home=patients.with_these_clinical_events(
      carehome_primis_codes,
      between = ["housebound_date","index_date - 1 day"]
    ),
  ),
  
###############################################################################
# CLINICAL/DEMOGRAPHIC COVARIATES
###############################################################################
  
  ## Age for JCVI group definition
  age = patients.age_as_of(
    "2021-03-31", # Age defined on 31 March 2021 as per coverage paper (https://bjgp.org/content/72/714/e51)
     return_expectations = {
      "rate": "universal",
      "int": {"distribution": "population_ages"},
    },
  ),
  
  ## Age at index
  age_index = patients.age_as_of(
    "index_date - 1 day",
    return_expectations = {
      "rate": "universal",
      "int": {"distribution": "population_ages"},
    },
  ),
  
  ## Sex
  sex = patients.sex(
    return_expectations = {
      "rate": "universal",
      "category": {"ratios": {"M": 0.49, "F": 0.51}},
    }
  ),
  
  ## BMI
  bmi = patients.categorised_as(
    {
      "Not obese": "DEFAULT",
      "Obese I (30-34.9)": """ bmi_value >= 30 AND bmi_value < 35""",
      "Obese II (35-39.9)": """ bmi_value >= 35 AND bmi_value < 40""",
      "Obese III (40+)": """ bmi_value >= 40 AND bmi_value < 100""",
      # set maximum to avoid any impossibly extreme values being classified as obese
    },
    bmi_value = patients.most_recent_bmi(
      on_or_after = "index_date - 5 years",
      minimum_age_at_measurement = 16
    ),
    return_expectations = {
      "rate": "universal",
      "category": {
        "ratios": {
          "Not obese": 0.7,
          "Obese I (30-34.9)": 0.1,
          "Obese II (35-39.9)": 0.1,
          "Obese III (40+)": 0.1,
        }
      },
    },
  ),
  
  ## Smoking
  smoking_status = patients.categorised_as(
    {
      "S": "most_recent_smoking_code = 'S'",
      "E": """
                 most_recent_smoking_code = 'E' OR (
                 most_recent_smoking_code = 'N' AND ever_smoked
                 )
            """,
      "N": "most_recent_smoking_code = 'N' AND NOT ever_smoked",
      "M": "DEFAULT",
    },
    return_expectations = {
      "category": {"ratios": {"S": 0.6, "E": 0.1, "N": 0.2, "M": 0.1}}
    },
    
    most_recent_smoking_code = patients.with_these_clinical_events(
      clear_smoking_codes,
      find_last_match_in_period = True,
      on_or_before = "index_date - 1 day",
      returning="category",
    ),
    
    ever_smoked=patients.with_these_clinical_events(
      filter_codes_by_category(clear_smoking_codes, include=["S", "E"]),
      on_or_before = "index_date - 1 day",
    ),
  ),
  
  ## Ethnicity
  ethnicity_6 = patients.with_these_clinical_events(
    ethnicity_6_codes,
    returning = "category",
    find_last_match_in_period = True,
    include_date_of_match = False,
    return_expectations = {
      "category": {"ratios": {"1": 0.2, "2": 0.2, "3": 0.2, "4": 0.2, "5": 0.2}},
      "incidence": 0.75,
    },
  ),
  
  ethnicity_6_sus = patients.with_ethnicity_from_sus(
    returning = "group_6",
    use_most_frequent_code = True,
    return_expectations = {
      "category": {"ratios": {"1": 0.2, "2": 0.2, "3": 0.2, "4": 0.2, "5": 0.2}},
      "incidence": 0.8,
    },
  ),
  
  ## Index of multiple deprivation
  imd = patients.categorised_as(
    {"0": "DEFAULT",
      "1": """index_of_multiple_deprivation >=1 AND index_of_multiple_deprivation < 32844*1/5""",
      "2": """index_of_multiple_deprivation >= 32844*1/5 AND index_of_multiple_deprivation < 32844*2/5""",
      "3": """index_of_multiple_deprivation >= 32844*2/5 AND index_of_multiple_deprivation < 32844*3/5""",
      "4": """index_of_multiple_deprivation >= 32844*3/5 AND index_of_multiple_deprivation < 32844*4/5""",
      "5": """index_of_multiple_deprivation >= 32844*4/5 """,
    },
    index_of_multiple_deprivation = patients.address_as_of(
      "index_date - 1 day",
      returning = "index_of_multiple_deprivation",
      round_to_nearest = 100,
    ),
    return_expectations = {
      "rate": "universal",
      "category": {
        "ratios": {
          "0": 0.01,
          "1": 0.20,
          "2": 0.20,
          "3": 0.20,
          "4": 0.20,
          "5": 0.19,
        }},
    },
  ),
  
  ## Region - NHS England 9 regions
  region = patients.registered_practice_as_of(
    "index_date - 1 day",
    returning = "nuts1_region_name",
    return_expectations = {
      "rate": "universal",
      "category": {
        "ratios": {
          "North East": 0.1,
          "North West": 0.1,
          "Yorkshire and The Humber": 0.1,
          "East Midlands": 0.1,
          "West Midlands": 0.1,
          "East": 0.1,
          "London": 0.2,
          "South West": 0.1,
          "South East": 0.1,},},
    },
  ),
  
  ## Rurality
  rural_urban = patients.address_as_of(
    "index_date - 1 days",
    returning="rural_urban_classification",
    return_expectations={
      "rate": "universal",
      "category": {"ratios": {1: 0.125, 2: 0.125, 3: 0.125, 4: 0.125, 5: 0.125, 6: 0.125, 7: 0.125, 8: 0.125}},
    },
  ),
  
###############################################################################
# COMORBIDITIES
###############################################################################

  ## Asthma
  asthma = patients.with_these_clinical_events(
    asthma_codes,
    returning = "binary_flag",
    find_last_match_in_period = True,
    on_or_before = "index_date - 1 day",
  ),
  
  ## Asplenia or Dysfunction of the Spleen codes
  asplenia = patients.with_these_clinical_events(
    spln_codes,
    returning = "binary_flag",
    find_last_match_in_period = True,
    on_or_before = "index_date - 1 day",
  ),
  
  ## Blood pressure
  bp_sys = patients.mean_recorded_value(
    systolic_blood_pressure_codes,
    on_most_recent_day_of_measurement = True,
    on_or_before = "index_date - 1 day",
    include_measurement_date = True,
    date_format = "YYYY-MM-DD", #include_month=True,
    return_expectations = {
      "incidence": 0.6,
      "float": {"distribution": "normal", "mean": 80, "stddev": 10},
    },
  ),
  
  bp_dias = patients.mean_recorded_value(
    diastolic_blood_pressure_codes,
    on_most_recent_day_of_measurement = True,
    on_or_before = "index_date - 1 day",
    include_measurement_date = True,
    date_format = "YYYY-MM-DD", #include_month=True,
    return_expectations ={
      "incidence": 0.6,
      "float": {"distribution": "normal", "mean": 120, "stddev": 10},
    },
  ),
  
  ## Cancer (non-haematological)
  cancer = patients.with_these_clinical_events(
    combine_codelists(
      lung_cancer_codes,
      other_cancer_codes
    ),
    returning = "binary_flag",
    find_last_match_in_period = True,
    on_or_before = "index_date - 1 day",
  ),
  
  ## Cancer (haematological)
  haem_cancer = patients.with_these_clinical_events(
    haem_cancer_codes,
    returning = "binary_flag",
    find_last_match_in_period = True,
    on_or_before = "index_date - 1 day",
  ),
  
  ## Chronic heart disease codes
  chd = patients.with_these_clinical_events(
    chd_codes,
    returning = "binary_flag",
    find_last_match_in_period = True,
    on_or_before = "index_date - 1 day",
  ),
  
  ## Chronic neurological disease (including Significant Learning Disorder)
  chronic_neuro_dis_inc_sig_learn_dis = patients.with_these_clinical_events(
    cnd_inc_sig_learn_dis_codes,
    returning = "binary_flag",
    find_last_match_in_period = True,
    on_or_before = "index_date - 1 day",
  ),
  
  ## Chronic respiratory disease
  chronic_resp_dis = patients.with_these_clinical_events(
    crs_codes,
    returning = "binary_flag",
    find_last_match_in_period = True,
    on_or_before = "index_date - 1 day",
  ),
  
  ## Chronic Liver disease codes
  cld = patients.with_these_clinical_events(
    cld_codes,
    returning = "binary_flag",
    find_last_match_in_period = True,
    on_or_before = "index_date - 1 day",
    date_format = "YYYY-MM-DD",
  ),
  
  ## Diabetes diagnosis codes
  diabetes = patients.with_these_clinical_events(
    diab_codes,
    returning = "binary_flag",
    find_last_match_in_period = True,
    on_or_before = "index_date - 1 day",
  ),
  
  ## Immunosuppression diagnosis
  immunosuppression_diagnosis_date = patients.with_these_clinical_events(
    immunosuppression_diagnosis_codes,
    returning = "date",
    find_last_match_in_period = True,
    on_or_before = "index_date - 1 day",
    date_format = "YYYY-MM-DD",
  ),
  
  ## Immunosuppression medication
  immunosuppression_medication_date = patients.with_these_medications(
    immunosuppression_medication_codes,
    returning = "date",
    find_last_match_in_period = True,
    on_or_before = "index_date - 1 day",
    date_format = "YYYY-MM-DD",
  ),
  
  ## Learning disabilities
  learning_disability = patients.with_these_clinical_events(
    learning_disability_codes,
    returning = "binary_flag",
    find_last_match_in_period = True,
    on_or_before = "index_date - 1 day"
  ),
  
  ### Severe mental illness
  sev_mental_ill = patients.with_these_clinical_events(
    sev_mental_ill_codes,
    returning = "binary_flag",
    find_last_match_in_period = True,
    on_or_before = "index_date - 1 day",
    date_format = "YYYY-MM-DD",
  ),
  
  ## Organ transplant
  organ_transplant = patients.with_these_clinical_events(
    organ_transplant_codes,
    returning = "binary_flag",
    find_last_match_in_period = True,
    on_or_before = "index_date - 1 day"
  ),
  
  ## Non-kidney transplant
  non_kidney_transplant = patients.with_these_clinical_events(
    non_kidney_transplant_codes,
    returning = "binary_flag",
    find_last_match_in_period = True,
    on_or_before = "index_date - 1 day"
  ),
  
  ################################################
  ############ pre-index COVID events ############
  ################################################
  
  ## Positive test prior to study period
  prior_positive_test_date = patients.with_test_result_in_sgss(
    pathogen = "SARS-CoV-2",
    test_result = "positive",
    returning = "date",
    date_format = "YYYY-MM-DD",
    on_or_before = "index_date - 1 day",
    find_first_match_in_period = True,
    restrict_to_earliest_specimen_date = False,
    return_expectations = {
      "date": {"earliest": "2020-02-01", "latest": "index_date - 1 day"}, # need both earliest/latest to obtain expected incidence
      "rate": "uniform",
      "incidence": 0.02,
    },
  ),
  
  ## Positive case identification prior to study period
  prior_primary_care_covid_case_date = patients.with_these_clinical_events(
    combine_codelists(
      covid_primary_care_code,
      covid_primary_care_positive_test,
      covid_primary_care_sequalae,
    ),
    returning = "date",
    date_format = "YYYY-MM-DD",
    on_or_before = "index_date - 1 day",
    find_first_match_in_period=True,
    return_expectations = {
      "date": {"earliest": "2020-02-01", "latest": "index_date - 1 day"}, # need both earliest/latest to obtain expected incidence
      "rate": "uniform",
      "incidence": 0.02,
    },
  ),

  ## Positive covid admission prior to study period
  prior_covid_hospitalisation_date = patients.admitted_to_hospital(
    returning = "date_admitted",
    with_these_diagnoses = covid_icd10,
    on_or_before = "index_date - 1 day",
    date_format = "YYYY-MM-DD",
    find_first_match_in_period = True,
    return_expectations = {
      "date": {"earliest": "2020-02-01", "latest": "index_date - 1 day"}, # need both earliest/latest to obtain expected incidence
      "rate": "uniform",
      "incidence": 0.005,
    },
  ),
  
  ################################################
  ############ pre-boost COVID events ############
  ################################################
  
  ## Positive test prior to study period
  prior_positive_test_date_boost = patients.with_test_result_in_sgss(
    pathogen = "SARS-CoV-2",
    test_result = "positive",
    returning = "date",
    date_format = "YYYY-MM-DD",
    on_or_before = "2021-08-31", # day before JCVI extended series recommendations
    find_first_match_in_period = True,
    restrict_to_earliest_specimen_date = False,
    return_expectations = {
      "date": {"earliest": "2020-02-01", "latest": "2021-08-31"}, # need both earliest/latest to obtain expected incidence
      "rate": "uniform",
      "incidence": 0.02,
    },
  ),
  
  ## Positive case identification prior to study period
  prior_primary_care_covid_case_date_boost = patients.with_these_clinical_events(
    combine_codelists(
      covid_primary_care_code,
      covid_primary_care_positive_test,
      covid_primary_care_sequalae,
    ),
    returning = "date",
    date_format = "YYYY-MM-DD",
    on_or_before = "2021-08-31", # day before JCVI extended series recommendations
    find_first_match_in_period=True,
    return_expectations = {
      "date": {"earliest": "2020-02-01", "latest": "2021-08-31"}, # need both earliest/latest to obtain expected incidence
      "rate": "uniform",
      "incidence": 0.02,
    },
  ),

  ## Positive covid admission prior to study period
  prior_covid_hospitalisation_date_boost = patients.admitted_to_hospital(
    returning = "date_admitted",
    with_these_diagnoses = covid_icd10,
    on_or_before = "2021-08-31", # day before JCVI extended series recommendations
    date_format = "YYYY-MM-DD",
    find_first_match_in_period = True,
    return_expectations = {
      "date": {"earliest": "2020-02-01", "latest": "2021-08-31"}, # need both earliest/latest to obtain expected incidence
      "rate": "uniform",
      "incidence": 0.005,
    },
  ),

###############################################################################
# ADDITIONAL VARIABLES FOR COMPARATIVE VE
###############################################################################

  ################################################
  ############ pre-vaccine events ################
  ################################################
  
  ## Positive test in pre-vaccination period
  prevax_positive_test_date = patients.with_test_result_in_sgss(
    pathogen = "SARS-CoV-2",
    test_result = "positive",
    returning = "date",
    date_format = "YYYY-MM-DD",
    between=["covid_vax_date_1 - 90 days","covid_vax_date_2 - 1 day"],
    find_first_match_in_period = True,
    restrict_to_earliest_specimen_date = False,
    return_expectations = {
      "date": {"earliest": "2020-09-01", "latest": "2021-04-01"}, # need both earliest/latest to obtain expected incidence
      "rate": "uniform",
      "incidence": 0.02,
    },
  ),
  
  ## Positive case identification in pre-vaccination period
  prevax_primary_care_covid_case_date = patients.with_these_clinical_events(
    combine_codelists(
      covid_primary_care_code,
      covid_primary_care_positive_test,
      covid_primary_care_sequalae,
    ),
    returning = "date",
    date_format = "YYYY-MM-DD",
    between=["covid_vax_date_1 - 90 days","covid_vax_date_2 - 1 day"],
    find_first_match_in_period=True,
    return_expectations = {
      "date": {"earliest": "2020-09-01", "latest": "2021-04-01"}, # need both earliest/latest to obtain expected incidence
      "rate": "uniform",
      "incidence": 0.02,
    },
  ),
  
  ## Positive covid admission in pre-vaccination period
  prevax_covid_hospitalisation_date = patients.admitted_to_hospital(
    returning = "date_admitted",
    with_these_diagnoses = covid_icd10,
    between=["covid_vax_date_1 - 90 days","covid_vax_date_2 - 1 day"],
    date_format = "YYYY-MM-DD",
    find_first_match_in_period = True,
    return_expectations = {
      "date": {"earliest": "2020-09-01", "latest": "2021-04-01"}, # need both earliest/latest to obtain expected incidence
      "rate": "uniform",
      "incidence": 0.005,
    },
  ),
  
  ## Count of tests (any) in pre-vaccination period
  prevax_tests_conducted_any = patients.with_test_result_in_sgss(
    pathogen = "SARS-CoV-2",
    test_result = "any",
    returning = "number_of_matches_in_period",
    between=["covid_vax_date_1 - 90 days","covid_vax_date_2 - 1 day"],
    restrict_to_earliest_specimen_date = False,
    return_expectations={
      "int": {"distribution": "normal", "mean": 4, "stddev": 1},
      "incidence": 0.05,
    },
  ),

  ################################################
  ############ events during study period ########
  ################################################
  
  ## Covid-related positive test after second dose
  postvax_positive_test_date = patients.with_test_result_in_sgss(
    pathogen = "SARS-CoV-2",
    test_result = "positive",
    between=["covid_vax_date_2",end_date],
    find_first_match_in_period = True,
    restrict_to_earliest_specimen_date = False,
    returning = "date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-02-01", "latest" : "2022-04-20"},
      "rate": "uniform",
      "incidence": 0.15,
    },
  ),
  
  ## Covid-related emergency admission after second dose
  postvax_covid_emergency_date = patients.attended_emergency_care(
    returning = "date_arrived",
    with_these_diagnoses = covid_emergency,
    between=["covid_vax_date_2",end_date],
    date_format = "YYYY-MM-DD",
    find_first_match_in_period = True,
    return_expectations = {
      "date": {"earliest": "2021-02-01", "latest" : "2022-04-20"},
      "rate": "uniform",
      "incidence": 0.05,
    },
  ),
    
  ## Covid-related unplanned admission after second dose
  postvax_covid_hospitalisation_date = patients.admitted_to_hospital(
    returning = "date_admitted",
    with_these_diagnoses = covid_icd10,
    with_admission_method = ["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
    between=["covid_vax_date_2",end_date],
    date_format = "YYYY-MM-DD",
    find_first_match_in_period = True,
    return_expectations = {
      "date": {"earliest": "2021-02-01", "latest" : "2022-04-20"},
      "rate": "uniform",
      "incidence": 0.05,
    },
  ),
  
  # Covid-related death after second dose
  postvax_covid_death_date = patients.with_these_codes_on_death_certificate(
    covid_icd10,
    returning = "date_of_death",
    between=["covid_vax_date_2",end_date],
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-02-01", "latest" : "2022-04-20"},
      "rate": "uniform",
      "incidence": 0.02
    },
  ),
  
  ################################################
  ############ Alternative endpoints #############
  ################################################

  ## Any test after second dose
  postvax_any_test_date = patients.with_test_result_in_sgss(
    pathogen = "SARS-CoV-2",
    test_result = "any",
    between=["covid_vax_date_2",end_date],
    find_first_match_in_period = True,
    restrict_to_earliest_specimen_date = False,
    returning = "date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-02-01", "latest" : "2022-04-20"},
      "rate": "uniform",
      "incidence": 0.5,
    },
  ),
  
  ## Any emergency admission after second dose
  postvax_any_emergency_date = patients.attended_emergency_care(
    returning = "date_arrived",
    between=["covid_vax_date_2",end_date],
    date_format = "YYYY-MM-DD",
    find_first_match_in_period = True,
    return_expectations = {
      "date": {"earliest": "2021-02-01", "latest" : "2022-04-20"},
      "rate": "uniform",
      "incidence": 0.05,
    },
  ),
  
  ## Any unplanned admission after second dose
  postvax_any_hospitalisation_date = patients.admitted_to_hospital(
    returning = "date_admitted",
    with_admission_method = ["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
    between=["covid_vax_date_2",end_date],
    date_format = "YYYY-MM-DD",
    find_first_match_in_period = True,
    return_expectations = {
      "date": {"earliest": "2021-02-01", "latest" : "2022-04-20"},
      "rate": "uniform",
      "incidence": 0.05,
    },
  ),
  
  # Any death after second dose
  postvax_any_death_date = patients.died_from_any_cause(
    returning="date_of_death",
    between=["covid_vax_date_2",end_date],
    date_format="YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-02-01", "latest" : "2022-04-20"},
      "rate": "uniform",
      "incidence": 0.02
    },
  ),
  
)
