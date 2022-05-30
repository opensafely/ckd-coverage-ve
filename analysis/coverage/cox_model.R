######################################

# This script:
# - imports processed data
# - fits multivariable stratified cox model(s) using the coxph package
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

## Import command-line arguments
args <- commandArgs(trailingOnly=TRUE)

## Set input and output pathways - default is dose 3 uptake by 20th April 2022
if(length(args)==0) {
  outcome = "vax3_date"
  outcome_label = "dose3"
  cutoff = as.Date("2022-05-11", format = "%Y-%m-%d")
} else {
  if (args[[1]]=="dose3") {
    outcome = "vax3_date"
    outcome_label = "dose3"
    cutoff = as.Date("2022-05-11", format = "%Y-%m-%d")
  } else if (args[[1]]=="dose4") {
    outcome = "vax4_date"
    outcome_label = "dose4"
    cutoff = as.Date("2022-05-11", format = "%Y-%m-%d")
  } else {
    # Print error if no argument specified
    stop("No outcome specified")
  }
}

## Import data
if (outcome_label == "dose3") {
  data_cohort <- read_rds(here::here("output", "data", "data_cohort_coverage.rds"))
} else {
  data_cohort <- read_rds(here::here("output", "data", "data_cohort_coverage_dose4.rds"))
}

## Import custom user functions and packages
source(here::here("analysis", "functions.R"))

## Amend data to cox ph format
data_cox <- data_cohort %>%
  mutate(
    # Start date
    start_date = as.Date("2020-12-08", format = "%Y-%m-%d"),
    
    # End date
    end_date = cutoff,
    
    # Determine date for selected outcome
    selected_outcome_date = get(outcome),
    
    # Censoring
    censor_date = pmin(death_date, 
                       dereg_date, 
                       end_date, 
                       na.rm=TRUE),
    
    # Follow-up time
    follow_up_time = tte(start_date, selected_outcome_date, censor_date),
    
    # COVID vaccination: 1 if vaccinated before end date, 0 otherwise
    covid_vax = dplyr::if_else((selected_outcome_date>censor_date) | is.na(selected_outcome_date), 0, 1),
    
    # Uncensored pre cut-off
    uncensored_at_cut_off = ifelse(
                              (is.na(death_date) & is.na(dereg_date)) | 
                              ((!is.na(death_date)) & death_date>end_date) | 
                              ((!is.na(dereg_date)) & dereg_date>end_date), 1, 0),
    
    # Vaccinated and uncensored at cut-off
    vax_uncensored_at_cut_off = ifelse(covid_vax==1 & uncensored_at_cut_off==1, 1, 0)
  )

## Create output directory
dir.create(here::here("output", "model"), showWarnings = FALSE, recursive=TRUE)

## Set variable list
var_list <- c("ageband2", "sex", "ethnicity", "imd",  "rural_urban_group",
              "ckd_5cat", "chronic_kidney_disease_stages_3_5",
              "care_home", "hscworker", "housebound", "endoflife",
              "prior_covid_cat", "immunosuppression", "mod_sev_obesity", "diabetes", "any_resp_dis", "chd", "cld", "asplenia", "cancer",
              "haem_cancer", "other_transplant", "chronic_neuro_dis_inc_sig_learn_dis","sev_mental_ill", "cev_other") 

## Set subset of variables to be used in minimally adjusted models
min_adj_list <- c("ageband2", "care_home", "hscworker")

## Set subset of variables to be used in partially adjusted models
adj_list <- c("ageband2", "care_home", "hscworker", "housebound", "endoflife", "rural_urban_group",
              "sex", "ethnicity", "imd", "prior_covid_cat", "haem_cancer", "immunosuppression")

## Pick out variables and recode list as factors
data_cox <- data_cox %>% select(all_of(var_list), follow_up_time, covid_vax, region, uncensored_at_cut_off, vax_uncensored_at_cut_off)
data_cox[,var_list] <- lapply(data_cox[,var_list], factor)

## Check all cases complete
if (all(complete.cases(data_cox))==FALSE) stop('Incomplete data for one or more patients in model') 

## Save cox model input data
write_rds(data_cox, here::here("output", "data", paste0("data_cox_coverage_",outcome_label,".rds")), compress="gz")





### MULTIVARIATE COX MODELS

## Define strata for full and stratified analyses
strata = c("full", "CKD3a", "CKD3b", "CKD4-5", "RRT (dialysis)", "RRT (Tx)")

## Run model in loop over each kidney disease stratum
for (s in 1:length(strata)) {
  ckd_group = strata[s]
  
  ## Select data subset
  if (ckd_group=="full") { 
    data_subset = data_cox 
    var_list_subset = var_list
    } else { 
    data_subset = subset(data_cox, ckd_5cat==ckd_group) 
    var_list_subset = var_list[var_list!="ckd_5cat"]
    }
  
  ## Exclude cev_other from RRT groups given that all will be 0 by definition
  if (ckd_group %in% c("RRT (dialysis)", "RRT (Tx)")) {
    var_list_subset = var_list_subset[var_list_subset!="cev_other"]
  }
  
  
  
  
    
  ## Fit stratified minimally adjusted models in loop
  for (i in 1:length(var_list_subset)) {
    # Select current variable from list, then select final list of model covariates
    var = var_list_subset[i]
    if (var %in% min_adj_list) { final_list = min_adj_list } else { final_list = c(var, min_adj_list) }
    
    # Fit model and extract output
    cox_minimal = coxph(as.formula(paste0("Surv(follow_up_time, covid_vax) ~",paste(final_list, collapse="+"),"+ strata(region)")), 
                        data = data_subset)
    summary = extract_model(cox_minimal) 
    
    # Pick out adjusted outputs for term of interest
    summary_var = summary[grepl(var, summary$term),]
    
    # Collate with previous model outputs
    if (i == 1) { cox_minimal_collated = summary_var } else { cox_minimal_collated = rbind(cox_minimal_collated,summary_var) }
  }
  names(cox_minimal_collated)[2:6] = paste0(names(cox_minimal_collated)[2:6],"_minimal")
  
  
  
  
  
  ## Fit stratified partially adjusted models in loop
  for (i in 1:length(var_list_subset)) {
    # Select current variable from list, then select final list of model covariates
    var = var_list_subset[i]
    if (var %in% adj_list) { final_list = adj_list } else { final_list = c(var, adj_list) }
    
    # Fit model and extract output
    cox_partial = coxph(as.formula(paste0("Surv(follow_up_time, covid_vax) ~",paste(final_list, collapse="+"),"+ strata(region)")), 
                        data = data_subset)
    summary = extract_model(cox_partial) 
    
    # Pick out adjusted outputs for term of interest
    if (var!="cancer") { summary_var = summary[grepl(var, summary$term),] } else { summary_var = summary[summary$term=="cancer1",] }
    
    # Collate with previous model outputs
    if (i == 1) { cox_partial_collated = summary_var } else { cox_partial_collated = rbind(cox_partial_collated,summary_var) }
  }
  names(cox_partial_collated)[2:6] = paste0(names(cox_partial_collated)[2:6],"_partial")
  
  
  
  
  
  ## Fit stratified fully adjusted model
  cox_full <- coxph(as.formula(paste0("Surv(follow_up_time, covid_vax) ~", paste(var_list_subset, collapse="+"),"+ strata(region)")),
                    data = data_subset)
  
  
  
  
  
  ### PROCESS MODEL OUTPUTS
  
  ## Create summary table derived from fully adjusted model
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

  ## Add variable to summary table that tabulates number uncensored at analysis end date and number vaccinated within uncensored population
  ## KM estimates will be primary measure of coverage, but these can be used as sense check to ensure KM approximates coverage at analysis cut-off
  tbl_full$n_uncensored_at_cut_off = NA
  tbl_full$n_vax_at_cut_off = NA
  
  ## Loop over variable list
  for (i in 1:length(var_list_subset)) {
    # i = current variable; j = current row number in summary table
    if (i == 1) { j = 1 }
    # Tabulate variable vs uncensored at cut-off, then incorporate into summary table
    var_uncensored = table(data.frame(data_subset)[,var_list_subset[i]], data_subset$uncensored_at_cut_off)[,2]
    tbl_full$n_uncensored_at_cut_off[j:(j+length(var_uncensored)-1)] = var_uncensored
    
    # Tabulate variable vs uncensored at cut-off, then incorporate into summary table
    var_vax = table(data.frame(data_subset)[,var_list_subset[i]], data_subset$vax_uncensored_at_cut_off)[,2]
    tbl_full$n_vax_at_cut_off[j:(j+length(var_vax)-1)] = var_vax
    
    # Update row number for next step in loop
    j = j + length(var_uncensored)
  }
  
  ## Add variable to summary table that tabulates KM estimates of coverage by group
  tbl_full$km_coverage = NA
  tbl_full$km_cml_event_floor = NA
  
  ## Convert data subset to dataframe to enable variable name processing within loop
  data_subset = data.frame(data_subset)
  
  ## Calculate KM coverage estimates via loop over variable list
  for (i in 1:length(var_list_subset)) {
    
    var_selected = var_list_subset[i]
    data_subset$var_selected = data_subset[,var_selected]
    levels = names(table(data_subset$var_selected))
    
    ## Calculate KM estimates within each unique level of covariate of interest
    for (j in 1:length(levels)) {
      data_level = subset(data_subset, var_selected==levels[j])
      surv = survfit(as.formula("Surv(follow_up_time, covid_vax) ~ 1"), data = data_level) %>% 
        broom::tidy() %>% 
        filter(estimate > 0) %>%
        mutate(
          N = plyr::round_any(max(n.risk, na.rm=TRUE),10),
          cml.event = cumsum(replace_na(n.event, 0)),
          cml.censor = cumsum(replace_na(n.censor, 0)),
          cml.event.floor = floor_any(cml.event, 10),
          cml.censor.floor = floor_any(cml.censor, 10),
          n.event.floor = diff(c(0,cml.event.floor)),
          n.censor.floor = diff(c(0,cml.censor.floor)),
          n.risk.floor = N - lag(cml.event.floor + cml.censor.floor,1,0),
          surv.floor = cumprod(1 - n.event.floor / n.risk.floor),
          cum.in.floor = 1 - surv.floor
        )
      # Merge KM estimates and N events (floor of 10) with main table
      tbl_full$km_coverage[tbl_full$variable==var_list_subset[i] & tbl_full$label==levels[j]] = round(max(surv$cum.in.floor)*100,1)
      tbl_full$km_cml_event_floor[tbl_full$variable==var_list_subset[i] & tbl_full$label==levels[j]] = max(surv$cml.event.floor)
    }
  }
    
  ## Check adjusted models have outputs for same variables, then merge
  if (all(cox_minimal_collated$term %in% tbl_full$term)==FALSE) stop('Minimally/fully adjusted model outputs have non-matching variables') 
  if (all(cox_partial_collated$term %in% tbl_full$term)==FALSE) stop('Partially/fully adjusted model outputs have non-matching variables') 
  
  tbl_summary <- tbl_full %>% 
    left_join(., cox_minimal_collated, by="term") %>% 
    left_join(., cox_partial_collated, by="term")
  
  ## Remove outputs for binary (yes/no) variables where not observed
  tbl_summary = subset(tbl_summary, var_nlevels>2 | variable=="sex" | (!is.na(p.value_full))) 
  tbl_summary$N_event = sum(data_subset$covid_vax)
  
  ## Add reference HRs to univariate and partially adjusted models
  tbl_summary$estimate_minimal[tbl_summary$reference_row==TRUE] = 1
  tbl_summary$estimate_partial[tbl_summary$reference_row==TRUE] = 1
  
  ## Relabel model output variables
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
  # var_label/group
  tbl_summary$var_label[tbl_summary$var_label=="ckd_5cat"] = "Kidney disease subgroup"
  tbl_summary$var_label[tbl_summary$label_clean %in% c("CKD3-5 code")] = "Primary care coding of kidney disease"
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

  ## Order factor levels for labels, group, and plot_group
  tbl_summary$label_clean = factor(tbl_summary$label_clean, levels = rev(tbl_summary$label_clean))
  tbl_summary$group = factor(tbl_summary$group, levels = unique(tbl_summary$group))
  tbl_summary$plot_group = factor(tbl_summary$plot_group, levels = unique(tbl_summary$plot_group))
  
  ## Round counts to rounding threshold to avoid disclosure issues and add non-event counts
  tbl_summary$N = plyr::round_any(tbl_summary$N, rounding_threshold)
  tbl_summary$N_event = plyr::round_any(tbl_summary$N_event, rounding_threshold)
  tbl_summary$n_obs = plyr::round_any(tbl_summary$n_obs, rounding_threshold)
  tbl_summary$n_event = plyr::round_any(tbl_summary$n_event, rounding_threshold)
  tbl_summary$n_uncensored_at_cut_off = plyr::round_any(tbl_summary$n_uncensored_at_cut_off, rounding_threshold)
  tbl_summary$n_vax_at_cut_off = plyr::round_any(tbl_summary$n_vax_at_cut_off, rounding_threshold)
  tbl_summary$perc_vax_cut_off = round(tbl_summary$n_vax_at_cut_off/tbl_summary$n_uncensored_at_cut_off*100,1)
  
  ## Add non-event count (overall and variable-specific)
  tbl_summary$n_nonevent = tbl_summary$n_obs - tbl_summary$n_event
  tbl_summary$n_unvax_uncensored = tbl_summary$n_uncensored_at_cut_off - tbl_summary$n_vax_at_cut_off
  
  ## Pick out key variables for outputs
  tbl_reduced <- data.frame(tbl_summary) %>%
    select(variable, group, plot_group, label_clean, # columns 1:4
           N, N_event, n_obs, # columns 5:7
           n_event, n_unvax_uncensored, # columns 8:9 - references for redaction
           n_uncensored_at_cut_off, n_vax_at_cut_off, perc_vax_cut_off, km_cml_event_floor, km_coverage, # columns 10:14
           estimate_minimal, conf.low_minimal, conf.high_minimal, p.value_minimal, # columns 15:18
           estimate_partial, conf.low_partial, conf.high_partial, p.value_partial, # columns 19:22
           estimate_full, conf.low_full, conf.high_full, p.value_full) # columns 23:26 
  
  ## Redact all model outputs if events/non-events <= redaction threshold (where non-events calculated among uncensored population)
  for (i in 1:nrow(tbl_reduced)) {
    if (
      (as.numeric(tbl_reduced[i,"n_event"])>0 & as.numeric(tbl_reduced[i,"n_event"])<=redaction_threshold) |
      (as.numeric(tbl_reduced[i,"n_unvax_uncensored"])>0 & as.numeric(tbl_reduced[i,"n_unvax_uncensored"])<=redaction_threshold)
      ) {  
      tbl_reduced[i,5:26] = "[Redacted]" 
      }
  }
  
  ## Save outputs
  write_rds(tbl_reduced, here::here("output", "model", paste0("mod_strat_coxph_redacted_",outcome_label,"_",ckd_group,".rds")), compress="gz")
  write_csv(tbl_reduced, here::here("output", "model", paste0("mod_strat_coxph_redacted_",outcome_label,"_",ckd_group,".csv")))
}
