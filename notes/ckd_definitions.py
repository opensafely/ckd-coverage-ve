### 1 ### nhs-covid-vaccination-coverage

    # Chronic kidney disease flag
    # as per official COVID-19 vaccine reporting specification
    # IF CKD_COV_DAT <> NULL (diagnoses) | Select | Next
    # IF CKD15_DAT = NULL  (No stages)   | Reject | Next
    # IF CKD35_DAT>=CKD15_DAT            | Select | Reject
    # (i.e. any diagnostic code, or most recent stage recorded >=3)
    ckd = patients.satisfying(
        """ckd_cov_dat 
            OR
            ckd35_dat
        """,
        # Chronic kidney disease diagnostic codes
        ckd_cov_dat=patients.with_these_clinical_events(
            ckd_cov,
            returning="date",
            find_first_match_in_period=True,
            on_or_before="index_date",
            date_format="YYYY-MM-DD",
        ),
        # Chronic kidney disease codes - all stages
        ckd15_dat=patients.with_these_clinical_events(
            ckd15,
            returning="date",
            find_last_match_in_period=True,
            on_or_before="index_date",
            date_format="YYYY-MM-DD",
        ),
        # Chronic kidney disease codes-stages 3 - 5
        # only on or after latest CKD1-5 code
        ckd35_dat=patients.with_these_clinical_events(
            ckd35,
            returning="date",
            find_last_match_in_period=True,
            between=["ckd15_dat", "index_date"],
            date_format="YYYY-MM-DD",
        ),
    )
    
    

### 2 ###
  ## Chronic kidney disease diagnostic
  chronic_kidney_disease_diagnostic = patients.with_these_clinical_events(
    chronic_kidney_disease_diagnostic_codes,
    returning = "date",
    find_first_match_in_period = True,
    on_or_before = "covid_vax_2_date",
    date_format = "YYYY-MM-DD",
  ),
  
  ## Chronic kidney disease codes - all stages
  chronic_kidney_disease_all_stages = patients.with_these_clinical_events(
    chronic_kidney_disease_all_stages_codes,
    returning = "date",
    find_last_match_in_period = True,
    on_or_before = "covid_vax_2_date",
    date_format = "YYYY-MM-DD",
  ),
  
  ## Chronic kidney disease codes-stages 3 - 5
  chronic_kidney_disease_all_stages_3_5 = patients.with_these_clinical_events(
    chronic_kidney_disease_all_stages_3_5_codes,
    returning = "date",
    find_last_match_in_period = True,
    on_or_before = "covid_vax_2_date",
    date_format = "YYYY-MM-DD",
  ),
  
  ## Chronic kidney disease - end-stage renal disease
  end_stage_renal = patients.with_these_clinical_events(
    ckd_codes,
    returning = "binary_flag",
    find_last_match_in_period = True,
    on_or_before = "covid_vax_2_date"
  ),


### Corresponding codelists - 1 - ###
ckd_cov = "codelists/primis-covid19-vacc-uptake-ckd_cov.csv", system="snomed", column="code"
ckd15 = "codelists/primis-covid19-vacc-uptake-ckd15.csv", system="snomed", column="code"
ckd35 = "codelists/primis-covid19-vacc-uptake-ckd35.csv", system="snomed", column="code"


### Corresponding codelists - 2 - ###
chronic_kidney_disease_diagnostic_codes = "codelists/primis-covid19-vacc-uptake-ckd_cov.csv", system = "snomed", column = "code"
chronic_kidney_disease_all_stages_codes = "codelists/primis-covid19-vacc-uptake-ckd15.csv", system="snomed", column="code"
chronic_kidney_disease_all_stages_3_5_codes = "codelists/primis-covid19-vacc-uptake-ckd35.csv", system="snomed", column="code"
ckd_codes = "codelists/opensafely-chronic-kidney-disease.csv", system = "ctv3", column = "CTV3ID"
