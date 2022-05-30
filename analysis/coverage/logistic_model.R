######################################

# This script:
# - imports processed data
# - fit multivariable logistic regression model on uncensored data
# - saves model outputs

######################################

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

## Set rounding and redaction thresholds
rounding_threshold = 5
redaction_threshold = 10

## Import data
data <- read_rds(here::here("output", "data", "data_cohort_coverage_logistic.rds")) %>%
  mutate(
    end_date = as.Date("2022-05-11", format = "%Y-%m-%d"),
    covid_vax = ifelse((!is.na(vax3_date)) & vax3_date<=end_date, 1, 0)
  )
 
## Import custom user functions and packages
source(here::here("analysis", "functions.R"))

## Create output directory
dir.create(here::here("output", "model"), showWarnings = FALSE, recursive=TRUE)

## Set variable list
var_list <- c("ageband2", "sex", "ethnicity", "imd",  "rural_urban_group",
              "ckd_5cat", "chronic_kidney_disease_stages_3_5",
              "care_home", "hscworker", "housebound", "endoflife",
              "prior_covid_cat", "immunosuppression", "mod_sev_obesity", "diabetes", "any_resp_dis", "chd", "cld", "asplenia", "cancer",
              "haem_cancer", "other_transplant", "chronic_neuro_dis_inc_sig_learn_dis","sev_mental_ill", "cev_other",
              "region") 

## Set subset of variables to be used in minimally adjusted models (same as Cox models but with region as covariate)
min_adj_list <- c("ageband2", "care_home", "hscworker", "region")

## Set subset of variables to be used in partially adjusted models (same as Cox models but with region as covariate)
adj_list <- c("ageband2", "care_home", "hscworker", "housebound", "endoflife", "rural_urban_group",
              "sex", "ethnicity", "imd", "prior_covid_cat", "haem_cancer", "immunosuppression", "region")

## Pick out variables and recode list as factors
data <- data %>% select(all_of(var_list), covid_vax)
data[,var_list] <- lapply(data[,var_list], factor)

## Check all cases complete
if (all(complete.cases(data))==FALSE) stop('Incomplete data for one or more patients in model') 

## Save cox model input data
write_rds(data, here::here("output", "data", "data_lr.rds"), compress="gz")





## Fit minimally adjusted models in loop
for (i in 1:length(var_list)) {
  ## Select current variable from list, then select final list of model covariates
  var = var_list[i]
  if (var %in% min_adj_list) { final_list = min_adj_list } else { final_list = c(var, min_adj_list) }
  
  ## Fit model and extract output
  lr_minimal = glm(as.formula(paste0("covid_vax ~", paste(final_list, collapse="+"))),
                   data = data, family = binomial)
  summary = extract_model_logistic(lr_minimal)
  
  ## Pick out adjusted outputs for term of interest and global intercept
  summary_var = summary[grepl("Intercept", summary$outcome) | grepl(var, summary$outcome),]
  
  ## Determine factor levels for variable of interest, provided at least 1 observation
  var_levels = table(data[,var_list[i]])
  summary_var$term = paste0(var_list[i], names(var_levels)[var_levels>0])
  
  ## Collate with previous model outputs
  if (i == 1) { lr_minimal_collated = summary_var } else { lr_minimal_collated = rbind(lr_minimal_collated,summary_var) }
}
names(lr_minimal_collated)[2:6] = paste0(names(lr_minimal_collated)[2:6],"_minimal")



## Fit partially adjusted models in loop
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
  
  ## Pick out adjusted outputs for term of interest
  if (var!="cancer") { 
    summary_var = summary[grepl("Intercept", summary$outcome) | grepl(var, summary$outcome),]
  } else { 
    summary_var = summary[grepl("Intercept", summary$outcome) | summary$outcome=="cancer1",]
  }
  
  ## Determine factor levels for variable of interest, provided at least 1 observation
  var_levels = table(data[,var_list[i]])
  summary_var$term = paste0(var_list[i], names(var_levels)[var_levels>0])
  
  ## Collate with previous model outputs
  if (i == 1) { lr_partial_collated = summary_var } else { lr_partial_collated = rbind(lr_partial_collated,summary_var) }
}
names(lr_partial_collated)[2:6] = paste0(names(lr_partial_collated)[2:6],"_partial")





## Fit fully adjusted model
lr_full <- glm(as.formula(paste0("covid_vax ~", paste(var_list, collapse="+"))),
                             data = data, family = binomial)





### PROCESS MODEL OUTPUTS

## Create summary table 
tbl_full <- tbl_regression(
  x = lr_full,
  exponentiate= TRUE,
  label = list(ageband2 = "Age", sex = "Sex", ethnicity = "Ethnicity", 
               imd = "IMD", rural_urban_group = "Setting", region = "Region")
) %>%
  as_gt() %>%
  .$`_data` %>%
  filter(
    !is.na(term)
  ) %>%
  mutate(
    var_label = if_else(row_type=="label", "", var_label),
    label_clean = if_else(reference_row %in% TRUE, paste0(label, " (ref)"),label),
    estimate_full = if_else(reference_row %in% TRUE, 1, estimate),
    variable = fct_inorder(variable),
    conf.low_full = round(conf.low,2),
    conf.high_full = round(conf.high,2),
    p.value_full = round(p.value,5)
  ) %>%
  mutate(
    estimate_full = round(estimate_full,2)
    )

## Check adjusted models have outputs for same variables, then merge
if (!all(tbl_full$term %in% lr_minimal_collated$term)) stop('Minimally/fully adjusted model outputs have non-matching variables') 
if (!all(tbl_full$term %in% lr_partial_collated$term)) stop('Partially/fully adjusted model outputs have non-matching variables') 

tbl_summary <- tbl_full %>% 
  left_join(., lr_minimal_collated %>% select(-outcome), by="term") %>% 
  left_join(., lr_partial_collated %>% select(-outcome), by="term")

## Remove outputs for binary (yes/no) variables where not observed
tbl_summary = subset(tbl_summary, var_nlevels>2 | variable=="sex" | (!is.na(p.value_full))) 
tbl_summary$N_event = sum(data$covid_vax)

## Add reference ORs to minimally and partially adjusted models
tbl_summary$estimate_minimal[tbl_summary$reference_row==TRUE] = 1
tbl_summary$estimate_partial[tbl_summary$reference_row==TRUE] = 1

## Relabel variables for plotting
tbl_summary$label_clean[tbl_summary$var_label=="care_home"] = "Care home resident"
tbl_summary$label_clean[tbl_summary$var_label=="hscworker"] = "Health/social care worker"
tbl_summary$label_clean[tbl_summary$var_label=="housebound"] = "Housebound"
tbl_summary$label_clean[tbl_summary$var_label=="endoflife"] = "End of life care"
tbl_summary$label_clean[tbl_summary$var_label=="prior_covid_cat"] = "Prior SARS-CoV-2"
tbl_summary$label_clean[tbl_summary$var_label=="immunosuppression"] = "Immunosuppression"
tbl_summary$label_clean[tbl_summary$var_label=="mod_sev_obesity"] = "Moderate/severe obesity"
tbl_summary$label_clean[tbl_summary$var_label=="diabetes"] = "Diabetes"
tbl_summary$label_clean[tbl_summary$var_label=="any_resp_dis"] = "Chronic respiratory disease (inc. asthma)"
tbl_summary$label_clean[tbl_summary$var_label=="chd"] = "Chronic heart disease"
tbl_summary$label_clean[tbl_summary$var_label=="cld"] = "Chronic liver disease"
tbl_summary$label_clean[tbl_summary$var_label=="asplenia"] = "Asplenia"
tbl_summary$label_clean[tbl_summary$var_label=="cancer"] = "Cancer (non-haematologic)"
tbl_summary$label_clean[tbl_summary$var_label=="haem_cancer"] = "Haematologic cancer"
tbl_summary$label_clean[tbl_summary$var_label=="other_transplant"] = "Organ transplant (non-kidney)"
tbl_summary$label_clean[tbl_summary$var_label=="chronic_neuro_dis_inc_sig_learn_dis"] = "Chronic neurological disease (inc. learning disability)"
tbl_summary$label_clean[tbl_summary$var_label=="sev_mental_ill"] = "Severe mental illness"
tbl_summary$label_clean[tbl_summary$var_label=="cev_other"] = "Clinically extremely vulnerable (other)"
tbl_summary$label_clean[tbl_summary$var_label=="chronic_kidney_disease_stages_3_5"] = "CKD3-5 code"

## Group variables for plotting
# var_label
tbl_summary$var_label[tbl_summary$var_label=="ckd_5cat"] = "Kidney disease subgroup"
tbl_summary$var_label[tbl_summary$label_clean %in% c("CKD3-5 primary care code")] = "Primary care coding of kidney disease"
tbl_summary$var_label[tbl_summary$label_clean %in% c("Care home resident", "Health/social care worker", "Housebound", "End of life care")] = "Risk group (occupation/access)"
tbl_summary$var_label[tbl_summary$label_clean %in% c("Prior SARS-CoV-2", "Immunosuppression", "Moderate/severe obesity", "Diabetes", "Chronic respiratory disease (inc. asthma)",
                                                     "Chronic heart disease", "Chronic liver disease","Asplenia", "Cancer (non-haematologic)", "Haematologic cancer", "Obesity", 
                                                     "Organ transplant (non-kidney)", "Chronic neurological disease (inc. learning disability)", "Severe mental illness", 
                                                     "Clinically extremely vulnerable (other)")] = "Risk group (clinical)"
tbl_summary$group = tbl_summary$var_label

# plot_group
tbl_summary$plot_group = "Demography"
tbl_summary$plot_group[tbl_summary$group %in% c("Kidney disease subgroup", "Primary care coding of kidney disease")] = "Risk group (CKD/RRT-related)"
tbl_summary$plot_group[tbl_summary$group %in% c("Risk group (occupation/access)")] = "Risk group (occupation/access)"
tbl_summary$plot_group[tbl_summary$group %in% c("Risk group (clinical)")] = "Risk group (clinical)"

## Order factor levels for labels, var_label, and var_group
tbl_summary$label_clean = factor(tbl_summary$label_clean, levels = rev(tbl_summary$label_clean))
tbl_summary$group = factor(tbl_summary$group, levels = unique(tbl_summary$group))
tbl_summary$plot_group = factor(tbl_summary$plot_group, levels = unique(tbl_summary$plot_group))

## Round counts to rounding threshold to avoid disclosure issues and add non-event counts
tbl_summary$N = plyr::round_any(tbl_summary$N,rounding_threshold)
tbl_summary$N_event = plyr::round_any(tbl_summary$N_event,rounding_threshold)
tbl_summary$n_obs = plyr::round_any(tbl_summary$n_obs,rounding_threshold)
tbl_summary$n_event = plyr::round_any(tbl_summary$n_event,rounding_threshold)

## Add non-event count (overall and variable-specific)
tbl_summary$N_nonevent = tbl_summary$N - tbl_summary$N_event
tbl_summary$n_nonevent = tbl_summary$n_obs - tbl_summary$n_event

## Pick out key variables for outputs
tbl_reduced <- data.frame(tbl_summary) %>%
  select(variable, group, plot_group, label_clean, # columns 1:4
         N, N_event, n_obs, # columns 5:7
         n_event, n_nonevent, # columns 8:9 - references for redaction
         estimate_minimal, conf.low_minimal, conf.high_minimal, # columns 10:12 
         estimate_partial, conf.low_partial, conf.high_partial, # columns 13:15
         estimate_full, conf.low_full, conf.high_full) # columns 16:18 

## Redact all model outputs if events/non-events <= redaction threshold
for (i in 1:nrow(tbl_reduced)) {
  if (
    (as.numeric(tbl_reduced[i,"n_event"])>0 & as.numeric(tbl_reduced[i,"n_event"])<=redaction_threshold) |
    (as.numeric(tbl_reduced[i,"n_nonevent"])>0 & as.numeric(tbl_reduced[i,"n_nonevent"])<=redaction_threshold)
    ) {  
      tbl_reduced[i,5:18] = "[Redacted]" 
    }
}

## Save outputs
write_rds(tbl_reduced, here::here("output", "model", "mod_strat_logistic_redacted.rds"), compress="gz")
write_csv(tbl_reduced, here::here("output", "model", "mod_strat_logistic_redacted.csv"))
