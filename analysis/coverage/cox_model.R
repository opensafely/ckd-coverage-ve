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
    censor_date = pmin(death_date, 
                       dereg_date, 
                       as.Date("2021-07-01", format = "%Y-%m-%d"), 
                       na.rm=TRUE),
    
    # Follow-up time
    follow_up_time = tte(start_date, vax2_date, censor_date),
    
    # COVID vaccination: 1 of vaccinated before end date, 0 otherwise
    covid_vax = dplyr::if_else((vax2_date>censor_date) | is.na(vax2_date), 0, 1),
    
    # Uncensored pre cut-off
    uncensored_at_cut_off = if_else( (is.na(death_date) & is.na(dereg_date)) | 
                                (!is.na(death_date) & death_date>end_date) | 
                                (!is.na(dereg_date) & dereg_date>end_date), 1, 0),
    
    # Unvaccinated and uncensored at cut-off
    unvax_uncensored_at_cut_off = if_else(covid_vax==0 & uncensored_at_cut_off==1, 1, 0)
  )

## Create output directory
dir.create(here::here("output", "model"), showWarnings = FALSE, recursive=TRUE)

## Set variable list
var_list <- c("ageband2", "care_home", "hscworker", "housebound", "endoflife", "rural_urban_group",
              "sex", "ethnicity", "imd", "ckd_7cat", 
              "prior_covid_cat", "immunosuppression", "mod_sev_obesity", "diabetes", "any_resp_dis", "chd", "cld", "asplenia", "cancer",
              "haem_cancer", "non_kidney_transplant", "chronic_neuro_dis_inc_sig_learn_dis","sev_mental_ill",
              "cev_other") 
adj_list <- c("ageband2", "care_home", "hscworker", "housebound", "endoflife", "rural_urban_group",
              "sex", "ethnicity", "imd", "prior_covid_cat")

## Pick out variables and recode list as factors
data_cox <- data_cox %>% select(all_of(var_list), follow_up_time, covid_vax, region, uncensored_at_cut_off, unvax_uncensored_at_cut_off)
data_cox[,var_list] <- lapply(data_cox[,var_list], factor)

## Set factor levels for CKD subgroup
data_cox$ckd_7cat <- factor(data_cox$ckd_7cat, levels = c("CKD3a (D-T-)", "CKD3b (D-T-)", "CKD4 (D-T-)", "CKD5 (D-T-)",
                                                          "CKD (D-T+)", "CKD (D+T-)", "CKD (D+T+)"))

## Check all cases complete
if (all(complete.cases(data_cox))==FALSE) stop('incomplete data for one or more patients in model') 

## Save cox model input data
write_rds(data_cox, here::here("output", "data", "data_cox.rds"), compress="gz")
write_csv(data_cox, here::here("output", "data", "data_cox.csv"))





### UNIVARIATE + MULTIVARIATE COXPH MODELS ----

## Fit stratified univariate models in loop
for (i in 1:length(var_list)) {
  cox_uni = coxph(as.formula(paste0("Surv(follow_up_time, covid_vax) ~",var_list[i],"+ strata(region)")), data = data_cox)
  summary = extract_model(cox_uni)
  if (i == 1) { cox_uni_collated = summary } else { cox_uni_collated = rbind(cox_uni_collated,summary) }
}
names(cox_uni_collated)[2:6] = paste0(names(cox_uni_collated)[2:6],"_uni")



## Fit stratified partially adjusted models in loop
for (i in 1:length(var_list)) {
  ## Select current variable from list, then select final list of model covariates
  var = var_list[i]
  if (var %in% adj_list) { final_list = adj_list } else { final_list = c(var, adj_list) }
  
  ## Fit model and extract output
  cox_partial = coxph(as.formula(paste0("Surv(follow_up_time, covid_vax) ~",paste(final_list, collapse="+"),"+ strata(region)")), 
                    data = data_cox)
  summary = extract_model(cox_partial) 
  
  ## Pick out adjusted outputs for term of interest
  summary_var = summary[grepl(var, summary$term),]
  
  ## Collate with previous model outputs
  if (i == 1) { cox_partial_collated = summary_var } else { cox_partial_collated = rbind(cox_partial_collated,summary_var) }
}
names(cox_partial_collated)[2:6] = paste0(names(cox_partial_collated)[2:6],"_partial")



## Fit stratified multivariate model - adjusted for all covariates
cox_full <- coxph(as.formula(paste0("Surv(follow_up_time, covid_vax) ~", paste(var_list, collapse="+"),"+ strata(region)")),
                             data = data_cox)
#cox_adj_full = extract_model(cox_full)
  




### PROCESS MODEL OUTPUTS ----

## Create summary table
tbl_full <- tbl_regression(
  x = cox_full,
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
    conf.low = round(conf.low,2),
    conf.high = round(conf.high,2),
    p.value = round(p.value,5)
  ) %>%
  mutate(
    estimate = round(estimate,2)
    )
#write_csv(tbl_full, here::here("output", "data", "tbl_full"))

## Add variable to summary table that tabulates number uncensored at analysis end date and number uncensored/unvaccinated
tbl_full$n_uncensored_at_cut_off = NA
tbl_full$n_unvax_at_cut_off = NA

## Loop over variable list
for (i in 1:length(var_list)) {
  # i = current variable; j = current row number in summary table
  if (i == 1) { j = 1 }
  # Tabulate variable vs uncensored at cut-off, then incorporate into summary table
  var_uncensored = table(data.frame(data_cox)[,var_list[i]], data_cox$uncensored_at_cut_off)[,2]
  tbl_full$n_uncensored_at_cut_off[j:(j+length(var_uncensored)-1)] = var_uncensored

  # Tabulate variable vs uncensored at cut-off, then incorporate into summary table
  var_unvax = table(data.frame(data_cox)[,var_list[i]], data_cox$unvax_uncensored_at_cut_off)[,2]
  tbl_full$n_unvax_at_cut_off[j:(j+length(var_unvax)-1)] = var_unvax
  
  # Update row number for next step in loop
  j = j + length(var_uncensored)
}

## Check unadjusted and adjusted models have outputs for same variables, then merge
if (all(cox_uni_collated$term %in% tbl_full$term)==FALSE) stop('univariate/multivariate model outputs have non-matching variables') 
if (all(cox_partial_collated$term %in% tbl_full$term)==FALSE) stop('partially/fully adjusted model outputs have non-matching variables') 

tbl_summary <- tbl_full %>% 
  # Merge with outputs of univariate analysis for comparison
  left_join(., cox_uni_collated, by="term") %>% 
  # Merge with outputs of partially adjusted analysis for comparison
  left_join(., cox_partial_collated, by="term")

## Remove outputs for binary (yes/no) variables where not observed
tbl_summary = subset(tbl_summary, var_nlevels>2 | variable=="sex" | !is.na(p.value)) 
tbl_summary$N_event = sum(data_cox$covid_vax)

## Add reference HRs to univariate and partially adjusted models
tbl_summary$estimate_uni[tbl_summary$reference_row==TRUE] = 1
tbl_summary$estimate_partial[tbl_summary$reference_row==TRUE] = 1

## Relabel variables for plotting
tbl_summary$label[tbl_summary$var_label=="care_home"] = "Care home resident"
tbl_summary$label[tbl_summary$var_label=="hscworker"] = "Health/social care worker"
tbl_summary$label[tbl_summary$var_label=="housebound"] = "Housebound"
tbl_summary$label[tbl_summary$var_label=="endoflife"] = "End of life care"
tbl_summary$label[tbl_summary$var_label=="prior_covid_cat"] = "Prior COVID"
tbl_summary$label[tbl_summary$var_label=="immunosuppression"] = "Immunosuppression"
tbl_summary$label[tbl_summary$var_label=="mod_sev_obesity"] = "Moderate/severe obesity"
tbl_summary$label[tbl_summary$var_label=="diabetes"] = "Diabetes"
tbl_summary$label[tbl_summary$var_label=="any_resp_dis"] = "Chronic respiratory disease (inc. asthma)"
tbl_summary$label[tbl_summary$var_label=="chd"] = "Chronic heart disease"
tbl_summary$label[tbl_summary$var_label=="cld"] = "Chronic liver disease"
tbl_summary$label[tbl_summary$var_label=="asplenia"] = "Asplenia"
tbl_summary$label[tbl_summary$var_label=="cancer"] = "Cancer (non-haematologic)"
tbl_summary$label[tbl_summary$var_label=="haem_cancer"] = "Haematologic cancer"
tbl_summary$label[tbl_summary$var_label=="non_kidney_transplant"] = "Organ transplant (non-kidney)"
tbl_summary$label[tbl_summary$var_label=="chronic_neuro_dis_inc_sig_learn_dis"] = "Chronic neurological disease (inc. learning disability)"
tbl_summary$label[tbl_summary$var_label=="sev_mental_ill"] = "Severe mental illness"
tbl_summary$label[tbl_summary$var_label=="cev_other"] = "Clinically extremely vulnerable (other)"

## Group variables for plotting
# var_label
tbl_summary$var_label[tbl_summary$var_label=="ckd_7cat"] = "CKD subgroup"
tbl_summary$var_label[tbl_summary$label %in% c("Care home resident", "Health/social care worker", "Housebound", "End of life care")] = "Risk group (occupation/access)"
tbl_summary$var_label[tbl_summary$label %in% c("Prior COVID", "Immunosuppression", "Moderate/severe obesity", "Diabetes", "Chronic respiratory disease (inc. asthma)",
                                               "Chronic heart disease", "Chronic liver disease","Asplenia", "Cancer (non-haematologic)", "Haematologic cancer", "Obesity", 
                                               "Organ transplant (non-kidney)", "Chronic neurological disease (inc. learning disability)", "Severe mental illness", 
                                               "Clinically extremely vulnerable (other)")] = "Risk group (clinical, non-CKD)"

# var_group
tbl_summary$var_group = "Demography"
tbl_summary$var_group[tbl_summary$var_label %in% c("CKD subgroup")] = "Risk group (CKD-related)"
tbl_summary$var_group[tbl_summary$var_label %in% c("Risk group (clinical, non-CKD)", "Prior exposure status")] = "Risk group (clinical, non-CKD)"

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
         N, N_event, N_nonevent, n_obs, n_event, n_nonevent, n_uncensored_at_cut_off, n_unvax_at_cut_off, N_uni, N_partial, # columns 5:14
         estimate, conf.low, conf.high, p.value, # columns 15:18
         estimate_uni, conf.low_uni, conf.high_uni, p.value_uni, # columns 19:22
         estimate_partial, conf.low_partial, conf.high_partial, p.value_partial) # columns 23:26 

## Redact all model outputs if less than or equal to 10 events/non-events in group
for (i in 1:nrow(tbl_reduced)) {
  if (min(as.numeric(tbl_reduced[i,5:14]), na.rm=T)<=10) { tbl_reduced[i,5:26] = "[Redacted]" }
}

## Save outputs
write_rds(tbl_reduced, here::here("output", "model", "mod_strat_coxph_redacted.rds"), compress="gz")
write_csv(tbl_reduced, here::here("output", "model", "mod_strat_coxph_redacted.csv"))
