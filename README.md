# OpenSAFELY: COVID-19 vaccine coverage and effectiveness in kidney disease patients

This is the code and configuration for OpenSAFELY.

* A preprint of the coverage analyses is available [here](https://www.medrxiv.org/content/10.1101/2022.06.14.22276391v1)
* The protocol draft is [here](https://docs.google.com/document/d/1w48W-bCMfn0RdkfxlU6fbRkPv3MnIYd3/edit#heading=h.gjdgxs)
* If you are interested in how we defined our variables, take a look at the [study definition](analysis/study_definition.py); this is written in `python`, but non-programmers should be able to get a relatively good idea of what is going on here
* If you are interested in how we defined our code lists, look in the [codelists folder](./codelists/).
* Developers and epidemiologists interested in the framework should review [the OpenSAFELY documentation](https://docs.opensafely.org)

### Summary of analysis code for factors associated with vaccine uptake
(1) **analysis/data_process.R**
* Reads in *inputdata.csv* (derived from *analysis/study_definition.py*)
* Sets variable classes
* Calculates eGFR from most recent creatinine reading
* Defines kidney disease subgroups based on eGFR and UK Renal Registry status (with the latter used for dialysis and kidney transplant recipients)
* Calculates sequential vaccination dates and assigns vaccine type to each date
* Outputs *data_processed.rds* and *data_processed.csv*

(2) **analysis/coverage/data_selection.R**
* Reads in *data_processed.rds*
* Applies following selection criteria to define final study population for coverage analyses:
(i) Aged 16+ with eGFR<60 in 2 years before 01 Dec 2020 or in UK Renal Registry on 31 Dec 2020;
(ii) No ambiguities in kidney disease status (primary care code indicating dialysis or kidney transplant but not in UK Renal Registry)
(iii) No missing demographic information (sex, IMD, ethnicity, region);
(iv) Maximum of 5 doses recorded and no vaccines administered at an interval of <14 days between doses 1 and 4.
* Also selects nested cohorts for logistic regression sensitivity analysis (uncensored throughout follow-up) and dose 4 analysis (at least one indicator for a third primary or second booster dose)
* Outputs *data_cohort.rds*, *data_cohort.csv*, *data_cohort_coverage_logistic.rds*, *data_cohort_coverage_dose4.rds*, and *flowchart.csv*

(3) **analysis/coverage/table_1_coverage.R**
* Generates table 1 of demographic/clinical covariates in study population (full and stratified by kidney disease subgroup)
* Outputs *table1_coverage_redacted_by_CKD.html* and *table1_coverage_redacted_by_CKD.csv* (full study population)
* Outputs *table1_coverage_redacted_by_CKD_dose4.html* and *table1_coverage_redacted_by_CKD_dose4.csv* (dose 4 subset)

(4) **analysis/coverage/cox_model.R**
* Reads in *data_cohort.rds*
* Calculates Cox proportional hazard models for time to completion of 3- or 4-dose primary vaccine series, censoring at death, de-registration, or analysis cut-off
* 25 variables of interest are defined in the vector *var_list*
* Minimally, partially, and fully adjusted models are calculated for each covariate, stratified by region
* Outputs *data_cox_coverage_{dose}_{subset}.rds*, *mod_strat_coxph_redacted_{dose}_{subset}.rds*, and *mod_strat_coxph_redacted_{dose}_{subset}.csv*

(5) **analysis/coverage/logistic_model.R** (replicates *cox_model.R* but with logistic regression)
* Reads in *data_cohort_coverage_logistic.rds*
* Calculates logistic regression models exploring factors associated with odds of receiving 3-dose primary series by analysis cut-off
* Outputs *data_lr.rds* , *mod_strat_logistic_redacted.rds*, and *mod_strat_logistic_redacted.csv*

(6) **analysis/coverage/vaccine_coverage.Rmd**
* Visualises the outputs of **table_1_coverage.R**
* Calculates and plots dose 1, 2, 3, 4, and 5 cumulative coverage over time (Kaplan-Meier estimates with step counts delayed until accrual of 10 events)
* Tabulates frequency tables for products, product combinations, and dosing intervals at doses 1–4 (rounding to nearest 5 with small number redaction)
* Tabulates and plots the outputs of **cox_model.R** (with rounding to nearest 5 and small number redaction)
* Tabulates and plots the outputs of **logistic_model.R** (with rounding to nearest 5 and small number redaction)

(6) **analysis/coverage/cox_model_check_coverage.Rmd**
* Run separately for each outcome (dose 3/dose 4) and each analysis subset
* Visualises survival curves and Schoenfeld residual plots for outputs of **cox_model.R**, enabling proportional hazards assumption to be checked

### Summary of analysis code for estimating comparative effectiveness of different vaccine regimens
Under development

# About the OpenSAFELY framework
The OpenSAFELY framework is a Trusted Research Environment (TRE) for electronic health records research in the NHS, with a focus on public accountability and research quality.

Read more at [OpenSAFELY.org](https://opensafely.org).

# Licences
As standard, research projects have a MIT license. 
