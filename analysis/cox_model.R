######################################

# This script:
# - imports processed data
# - fit univariate and multivariable stratified cox model(s) using the coxph package
# - saves models

######################################


### Preliminaries ----

## Import libraries
library('here')
library('readr')
library('tidyr')
library('tidyverse')
library('lubridate')
library('survival')
library('gtsummary')
library('gt')
library('survminer')
library('glue')
library('fs')


## Import data
data_cohort <- read_rds(here::here("output", "data", "data_cohort_coverage.rds"))

## Import custom user functions and packages
source(here::here("analysis", "functions.R"))

## Amend data to cox ph format
data_cox <- data_cohort %>%
  mutate(
    # Start date
    start_date = as.Date("2020-12-08", format = "%Y-%m-%d"),
    
    # End date
    end_date = as.Date("2021-07-01", format = "%Y-%m-%d"),
    
    # Censoring
    censor_date = pmin(as.Date("2021-07-01", format = "%Y-%m-%d"), 
                       na.rm=TRUE),
    
    # Follow-up time
    follow_up_time = tte(start_date, vax2_date, censor_date),
    
    # COVID vaccination: 1 of vaccinated before end date, 0 otherwise
    covid_vax = as.integer(ifelse(is.na(vax2_date) | vax2_date>end_date, 0, 1)),
  )

## Create output directory
dir.create(here::here("output", "model"), showWarnings = FALSE, recursive=TRUE)

## Set variable list
var_list <- c("ageband2", "sex", "ethnicity", "imd", "rural_urban_group",
              "chronic_kidney_disease_diagnostic", "dialysis", "kidney_transplant", "chronic_kidney_disease_stages_3_5", "prior_covid_cat",
              "cev", "care_home", "hscworker", "endoflife", "housebound", 
              "immunosuppression", "chronic_resp_dis", "diabetes", "cld",
              "chd", "asplenia", "cancer", "obesity", "chronic_neuro_dis_inc_sig_learn_dis", "sev_mental_ill") 

## Pick out variables, recode list as factors, and pikc complete cases
data_cox <- data_cox %>% select(all_of(var_list), follow_up_time, covid_vax, region)
data_cox[,var_list] <- lapply(data_cox[,var_list], factor)
data_cox <- data_cox[complete.cases(data_cox),]


## Save cox model input data
write_rds(data_cox, here::here("output", "data", "data_cox.rds"), compress="gz")
write_csv(data_cox, here::here("output", "data", "data_cox.csv"))



### UNIVARIATE + MULTIVARIATE COXPH MODELS ----

## Fit stratified univariate models in loop
for (i in 1:length(var_list)) {
  uni_model = coxph(as.formula(paste0("Surv(follow_up_time, covid_vax) ~",var_list[i],"+ strata(region)")), data = data_cox)
  summary = extract_model(uni_model)
  if (i == 1) { cox_uni = summary } else { cox_uni = rbind(cox_uni,summary) }
}
names(cox_uni)[2:6] = paste0(names(cox_uni)[2:6],"_uni")

## Fit stratified multivariate model - adjusted for all covariates
adj_model <- coxph(as.formula(paste0("Surv(follow_up_time, covid_vax) ~", paste(var_list, collapse="+"),"+ strata(region)")),
                             data = data_cox)
cox_adj = extract_model(adj_model)
  




### PROCESS MODEL OUTPUTS ----

## Create summary table for plot
tbl_summary <- tbl_regression(
  x = adj_model,
  pvalue_fun = ~style_pvalue(.x, digits=3),
  tidy_fun = tidy_coxph,
  exponentiate= TRUE,
  label = list(ageband2 = "Age", sex = "Sex", ethnicity = "Ethnicity", 
               imd = "IMD", rural_urban_group = "Setting")
  ) %>%
  as_gt() %>%
  .$`_data` %>%
  filter(
    !is.na(term)
  ) %>%
  mutate(
    var_label = if_else(row_type=="label", "", var_label),
    label = if_else(reference_row %in% TRUE, paste0(label, " (ref)"),label),
    estimate = if_else(reference_row %in% TRUE, 1, estimate),
    variable = fct_inorder(variable),
    conf.low = round(conf.low, 2),
    conf.high = round(conf.high,2),
    p.value = round(p.value,4)
  ) %>%
  mutate(
    estimate = round(estimate,2)
    )
#write_csv(tbl_summary, here::here("output", "data", "tbl_summary.csv"))


## Remove outputs for binary (yes/no) variables where not observed
tbl_summary = subset(tbl_summary, var_nlevels>2 | variable=="sex" | !is.na(p.value)) %>% 
  # merge with outputs of univariate analysis for comparison
  left_join(., cox_uni, by="term")
tbl_summary$N_event = sum(data_cox$covid_vax)

## Relabel variables for plotting
tbl_summary$label[tbl_summary$var_label=="chronic_kidney_disease_diagnostic"] = "CKD diagnostic code"
tbl_summary$label[tbl_summary$var_label=="dialysis"] = "Dialysis"
tbl_summary$label[tbl_summary$var_label=="kidney_transplant"] = "Kidney transplant"
tbl_summary$label[tbl_summary$var_label=="chronic_kidney_disease_stages_3_5"] = "CKD stage 3-5 code"
tbl_summary$label[tbl_summary$var_label=="prior_covid_cat"] = "Prior COVID"
tbl_summary$label[tbl_summary$var_label=="cev"] = "Clinically extremely vulnerable"
tbl_summary$label[tbl_summary$var_label=="care_home"] = "Care home resident"
tbl_summary$label[tbl_summary$var_label=="hscworker"] = "Health/social care worker"
tbl_summary$label[tbl_summary$var_label=="endoflife"] = "End of life care"
tbl_summary$label[tbl_summary$var_label=="housebound"] = "Housebound"
tbl_summary$label[tbl_summary$var_label=="immunosuppression"] = "Immunosuppression"
tbl_summary$label[tbl_summary$var_label=="chronic_resp_dis"] = "Chronic respiratory disease"
tbl_summary$label[tbl_summary$var_label=="diabetes"] = "Diabetes"
tbl_summary$label[tbl_summary$var_label=="cld"] = "Chronic liver disease"
tbl_summary$label[tbl_summary$var_label=="chd"] = "Chronic heart disease"
tbl_summary$label[tbl_summary$var_label=="asplenia"] = "Asplenia"
tbl_summary$label[tbl_summary$var_label=="cancer"] = "Cancer"
tbl_summary$label[tbl_summary$var_label=="obesity"] = "Obesity"
tbl_summary$label[tbl_summary$var_label=="chronic_neuro_dis_inc_sig_learn_dis"] = "Chronic neurological disease (including learning disability)"
tbl_summary$label[tbl_summary$var_label=="sev_mental_ill"] = "Severe mental illness"

## Group variables for plotting
# var_label
tbl_summary$var_label[tbl_summary$label %in% c("CKD diagnostic code", "Dialysis", "Kidney transplant", "CKD stage 3-5 code")] = "CKD"
tbl_summary$var_label[tbl_summary$label %in% c("Prior COVID", "Clinically extremely vulnerable", "Care home resident", "Health/social care worker", 
                                               "End of life care", "Housebound", "Immunosuppression", "Chronic respiratory disease",
                                               "Diabetes", "Chronic liver disease", "Chronic heart disease", "Asplenia", "Cancer", "Obesity", 
                                               "Chronic neurological disease (including learning disability)", "Severe mental illness")] = "Other"
# var_group
tbl_summary$var_group = "Demography"
tbl_summary$var_group[tbl_summary$var_label %in% c("CKD")] = "Clinical (CKD-related)"
tbl_summary$var_group[tbl_summary$var_label %in% c("Other")] = "Clinical (other)"

# order factor levels for labels, var_label, and var_group
tbl_summary$label = factor(tbl_summary$label, levels = rev(tbl_summary$label))
tbl_summary$var_label = factor(tbl_summary$var_label, levels = unique(tbl_summary$var_label))
tbl_summary$var_group = factor(tbl_summary$var_group, levels = unique(tbl_summary$var_group))

## Round counts to nearest 5 to avoid disclosure issues
tbl_summary$N = plyr::round_any(tbl_summary$N,5)
tbl_summary$N_event = plyr::round_any(tbl_summary$N_event,5)
tbl_summary$n_obs = plyr::round_any(tbl_summary$n_obs,5)
tbl_summary$n_event = plyr::round_any(tbl_summary$n_event,5)
tbl_summary$N_uni = plyr::round_any(tbl_summary$N_uni,5)

## Add non-event count (overall and variable-specific)
tbl_summary$N_nonevent = tbl_summary$N - tbl_summary$N_event
tbl_summary$n_nonevent = tbl_summary$n_obs - tbl_summary$n_event

## Pick out key variables for outputs
tbl_reduced <- data.frame(tbl_summary) %>%
  select(variable, var_label, var_group, label, # columns 1:4
         N, N_event, N_nonevent, n_obs, n_event, n_nonevent, N_uni, # columns 5:11
         estimate, std.error, conf.low, conf.high, p.value, # columns 12:16
         estimate_uni, conf.low_uni, conf.high_uni, p.value_uni) # columns 17:20 

## Redact all model outputs if less than or equal to 10 events/non-events in group
for (i in 1:nrow(tbl_reduced)) {
  if (min(as.numeric(tbl_reduced[i,5:11]), na.rm=T)<=10) { tbl_reduced[i,5:20] = "[Redacted]" }
}

## Save outputs
write_rds(tbl_reduced, here::here("output", "model", "mod_strat_coxph_redacted.rds"), compress="gz")
write_csv(tbl_reduced, here::here("output", "model", "mod_strat_coxph_redacted.csv"))
