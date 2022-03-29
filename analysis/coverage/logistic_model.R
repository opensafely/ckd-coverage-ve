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
data <- read_rds(here::here("output", "data", "data_cohort_coverage_logistic.rds")) %>%
  mutate(
    end_date = as.Date("2021-07-01", format = "%Y-%m-%d"),
    # COVID vaccination: 1 of vaccinated before end date, 0 otherwise
    covid_vax = dplyr::if_else(!is.na(vax2_date) & vax2_date<=end_date, 1, 0)
  )
 
## Import custom user functions and packages
source(here::here("analysis", "functions.R"))

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
data <- data %>% select(all_of(var_list), covid_vax, region)
data[,var_list] <- lapply(data[,var_list], factor)

## Set factor levels for CKD subgroup
data$ckd_7cat <- factor(data$ckd_7cat, levels = c("CKD3a (D-T-)", "CKD3b (D-T-)", "CKD4 (D-T-)", "CKD5 (D-T-)",
                                                          "CKD (D-T+)", "CKD (D+T-)", "CKD (D+T+)"))

## Check all cases complete
if (all(complete.cases(data))==FALSE) stop('incomplete data for one or more patients in model') 

## Save cox model input data
write_rds(data, here::here("output", "data", "data_lr.rds"), compress="gz")





### UNIVARIATE + MULTIVARIATE LOGISTIC REGRESSION MODELS ----

## Fit univariate models in loop
for (i in 1:length(var_list)) {
  lr_uni = glm(as.formula(paste0("covid_vax ~",var_list[i])), data = data, family = binomial)
  summary = extract_model_logistic(lr_uni)
  
  ## Determine factor levels for variable of interest, provided at least 1 observation
  var_levels = table(data[,var_list[i]])
  summary$term = paste0(var_list[i], names(var_levels)[var_levels>0])
  if (i == 1) { lr_uni_collated = summary } else { lr_uni_collated = rbind(lr_uni_collated,summary) }
}
names(lr_uni_collated)[2:6] = paste0(names(lr_uni_collated)[2:6],"_uni")



## Fit stratified partially adjusted models in loop
for (i in 1:length(var_list)) {
  ## Select current variable from list, then select final list of model covariates
  var = var_list[i]
  if (var %in% adj_list) { final_list = adj_list } else { final_list = c(var, adj_list) }
  
  ## Fit model and extract output
  lr_partial = glm(as.formula(paste0("covid_vax ~", paste(final_list, collapse="+"))),
                          data = data, family = binomial)
  summary = extract_model_logistic(lr_partial)
  
  ## Pick out adjusted outputs for term of interest and global intercept
  summary_var = summary[grepl("Intercept", summary$outcome) | grepl(var, summary$outcome),]
  
  ## Determine factor levels for variable of interest, provided at least 1 observation
  var_levels = table(data[,var_list[i]])
  summary_var$term = paste0(var_list[i], names(var_levels)[var_levels>0])
  
  ## Collate with previous model outputs
  if (i == 1) { lr_partial_collated = summary_var } else { lr_partial_collated = rbind(lr_partial_collated,summary_var) }
}
names(lr_partial_collated)[2:6] = paste0(names(lr_partial_collated)[2:6],"_partial")



## Fit multivariate model - adjusted for all covariates
lr_full <- glm(as.formula(paste0("covid_vax ~", paste(var_list, collapse="+"))),
                             data = data, family = binomial)





### PROCESS MODEL OUTPUTS ----

## Create summary table 
tbl_full <- tbl_regression(
  x = lr_full,
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

## Check unadjusted and adjusted models have outputs for same variables, then merge
if (!all(tbl_full$term %in% lr_uni_collated$term)) stop('univariate/multivariate model outputs have non-matching variables') 
if (!all(tbl_full$term %in% lr_partial_collated$term)) stop('partially/fully adjusted model outputs have non-matching variables') 

tbl_summary <- tbl_full %>% 
    # Merge with outputs of univariate analysis for comparison
    left_join(., lr_uni_collated %>% select(-outcome), by="term") %>% 
    # Merge with outputs of partially adjusted analysis for comparison
    left_join(., lr_partial_collated %>% select(-outcome), by="term")

## Remove outputs for binary (yes/no) variables where not observed
tbl_summary = subset(tbl_summary, var_nlevels>2 | variable=="sex" | !is.na(p.value)) 
tbl_summary$N_event = sum(data$covid_vax)

## Add reference ORs to univariate and partially adjusted models
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
         N, N_event, N_nonevent, n_obs, n_event, n_nonevent, N_uni, N_partial, # columns 5:12
         estimate, conf.low, conf.high, p.value, # columns 13:16
         estimate_uni, conf.low_uni, conf.high_uni, p.value_uni, # columns 17:20
         estimate_partial, conf.low_partial, conf.high_partial, p.value_partial) # columns 21:24 

## Redact all model outputs if less than or equal to 10 events/non-events in group
for (i in 1:nrow(tbl_reduced)) {
  if (min(as.numeric(tbl_reduced[i,5:12]), na.rm=T)<=10) { tbl_reduced[i,5:24] = "[Redacted]" }
}

## Save outputs
write_rds(tbl_reduced, here::here("output", "model", "mod_strat_logistic_redacted.rds"), compress="gz")
write_csv(tbl_reduced, here::here("output", "model", "mod_strat_logistic_redacted.csv"))
