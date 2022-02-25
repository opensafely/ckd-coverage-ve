# OpenSAFELY: COVID-19 vaccine coverage and effectiveness in chronic kidney disease patients

This is the code and configuration for OpenSAFELY.

* The protocol draft is [here](https://docs.google.com/document/d/1w48W-bCMfn0RdkfxlU6fbRkPv3MnIYd3/edit#heading=h.gjdgxs)
* If you are interested in how we defined our variables, take a look at the [study definition](analysis/study_definition.py); this is written in `python`, but non-programmers should be able to understand what is going on there
* If you are interested in how we defined our code lists, look in the [codelists folder](./codelists/).
* Developers and epidemiologists interested in the framework should review [the OpenSAFELY documentation](https://docs.opensafely.org)

### Summary of analysis code
(1) **data_process.R**
* Reads in *inputdata.csv*
* Sets variable classes
* Calculates eGFR from most recent creatinine reading
* Calculates sequential vaccination dates and assigns vaccine type to each date
* Outputs *data_processed.rds* and *data_processed.csv*

(2) **data_selection.R**
* Reads in *data_processed.rds*
* Applies following selection criteria to define final study population for coverage analyses:
(i) Aged 16+ with eGFR<60 in 2 years before 01 Dec 2020 or has any previous dialysis code/kidney transplant code/ CKD diagnostic code/CKD3-5 code;
(ii) eGFR<60 or dialysis/kidney transplant code",
(iii) No missing demographic information (sex, IMD, ethnicity, region);
(iv) Maximum of 4 doses recorded;
(v) No vaccines administered at an interval of <14 days",
* Outputs *data_cohort.rds*, *data_cohort.csv*, and *flowchart.csv*

(3) **table_1.R**
* Generates table 1 of demographic/clinical covariates in study population
* Outputs *table1_redacted.html* and *table1_redacted.csv*

(4) **cox_model.R**
* Reads in *data_cohort.rds*
* Calculates Cox proportional hazard models for time to completion of 2-dose primary vaccine series, censoring at death, de-registration, or 01 July 2021
* 21 variables of interest are defined in the vector *var_list*
* Univariate and multivariate (fully adjusted) models are calculated for each covariate, stratified by region
* Outputs *data_cox.rds* (individuals with complete data for all covariates), *mod_strat_coxph_redacted.rds*, and *mod_strat_coxph_redacted.csv*

(5) **logistic_model.R** (replicates *cox_model.R* but with logistic regression)
* Reads in *data_cox.rds*
* Calculates univariate and multivariate (fully adjusted) logistic regression models exploring factors associated with odds of receiving 2-dose primary series by 01 July 2021
* Outputs *mod_strat_logistic_redacted.rds*, and *mod_strat_logistic_redacted.csv*

(6) **vaccine_coverage.Rmd**
* Visualises the outputs of **table_1.R**
* Calculates and plots dose 1, 2, 3, and 4 coverage over time (rounded to nearest 10)
* Calculates and plots product-specific (ChAdOx1-S, Moderna, Pfizer, or Heterologous) 2-dose coverage over time (rounded to nearest 10)
* Tabulates frequency tables for products, product combinations, and dosing intervals at doses 1–4 (rounding to nearest 5 with small number redaction)
* Tabulates and plots the outputs of **cox_model.R** (with rounding to nearest 5 and small number redaction)
* Tabulates and plots the outputs of **logistic_model.R** (with rounding to nearest 5 and small number redaction)

# About the OpenSAFELY framework

The OpenSAFELY framework is a Trusted Research Environment (TRE) for electronic
health records research in the NHS, with a focus on public accountability and
research quality.

Read more at [OpenSAFELY.org](https://opensafely.org).

# Licences
As standard, research projects have a MIT license. 
