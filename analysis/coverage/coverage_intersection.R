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
  cutoff = as.Date("2022-08-31", format = "%Y-%m-%d")
} else {
  if (args[[1]]=="dose3") {
    outcome = "vax3_date"
    outcome_label = "dose3"
    cutoff = as.Date("2022-08-31", format = "%Y-%m-%d")
  } else if (args[[1]]=="dose4") {
    outcome = "vax4_date"
    outcome_label = "dose4"
    cutoff = as.Date("2022-08-31", format = "%Y-%m-%d")
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
    covid_vax = dplyr::if_else((selected_outcome_date>censor_date) | is.na(selected_outcome_date), 0, 1)
  )


## Important covariate levels from table 1 output
tab <- read_rds(here::here("output", "tables", "table1_coverage_redacted_by_CKD.rds"))[,1:2] %>%
  filter(Group != "N" & Group != "Ethnicity" & Group != "IMD" & Group != "Region" & Group != "JCVI group" &
           Variable != "Dialysis code" & Variable != "Kidney transplant code")

## Specify variable levels
tab$var = NA
tab$var[tab$Group=="Age"] = "ageband2"
tab$var[tab$Group=="Sex"] = "sex"
tab$var[tab$Group=="Setting"] = "rural_urban_group"
tab$var[tab$Group=="Kidney disease subgroup"] = "ckd_5cat"
tab$var[tab$Variable=="CKD3-5 code"] = "chronic_kidney_disease_stages_3_5"
tab$var[tab$Variable=="Care home resident"] = "care_home"
tab$var[tab$Variable=="Health/social care worker"] = "hscworker"
tab$var[tab$Variable=="Housebound"] = "housebound"
tab$var[tab$Variable=="End of life care"] = "endoflife"
tab$var[tab$Variable=="Prior SARS-CoV-2"] = "prior_covid_cat"
tab$var[tab$Variable=="Immunosuppression"] = "immunosuppression"
tab$var[tab$Variable=="Moderate/severe obesity"] = "mod_sev_obesity"
tab$var[tab$Variable=="Diabetes"] = "diabetes"
tab$var[tab$Variable=="Chronic respiratory disease (inc. asthma)"] = "any_resp_dis"
tab$var[tab$Variable=="Chronic heart disease"] = "chd"
tab$var[tab$Variable=="Chronic liver disease"] = "cld"
tab$var[tab$Variable=="Asplenia"] = "asplenia"
tab$var[tab$Variable=="Cancer (non-haematologic)"] = "cancer"
tab$var[tab$Variable=="Haematologic cancer"] = "haem_cancer"
tab$var[tab$Variable=="Organ transplant (non-kidney)"] = "other_transplant"
tab$var[tab$Variable=="Chronic neurological disease (inc. learning disability)"] = "chronic_neuro_dis_inc_sig_learn_dis"
tab$var[tab$Variable=="Severe mental illness"] = "sev_mental_ill"
tab$var[tab$Variable=="Clinically extremely vulnerable (other)"] = "cev_other"

## Set variable list
var_list <- unique(tab$var)

## Pick out variables and recode list as factors
data_cox <- data_cox %>% select(all_of(var_list), ethnicity, imd, follow_up_time, covid_vax, region)
data_cox[,var_list] <- lapply(data_cox[,var_list], factor)

## Check all cases complete
if (all(complete.cases(data_cox))==FALSE) stop('Incomplete data for one or more patients in model') 

## Outer loop for IMD/ethnicity subset 
for (e in 1:10) {
  tab_group = tab
  tab_group$km_coverage_ul = tab_group$km_coverage_ll = tab_group$km_coverage = tab_group$km_cml_event = tab_group$N = NA
  
  if (e==1) { imd_subset = subset(data_cox, ethnicity=="White" & imd=="1 (most deprived)") }
  if (e==2) { imd_subset = subset(data_cox, ethnicity=="White" & imd=="2") }
  if (e==3) { imd_subset = subset(data_cox, ethnicity=="White" & imd=="3") }
  if (e==4) { imd_subset = subset(data_cox, ethnicity=="White" & imd=="4") }
  if (e==5) { imd_subset = subset(data_cox, ethnicity=="White" & imd=="5 (least deprived)") }
  if (e==6) { imd_subset = subset(data_cox, ethnicity!="White" & imd=="1 (most deprived)") }
  if (e==7) { imd_subset = subset(data_cox, ethnicity!="White" & imd=="2") }
  if (e==8) { imd_subset = subset(data_cox, ethnicity!="White" & imd=="3") }
  if (e==9) { imd_subset = subset(data_cox, ethnicity!="White" & imd=="4") }
  if (e==10) { imd_subset = subset(data_cox, ethnicity!="White" & imd=="5 (least deprived)") }

    ## Calculate KM coverage estimates via loop over variable list
    for (v in 1:length(var_list)) {
      
      var_subset = data.frame(imd_subset)
      var_selected = var_list[v]
      var_subset$var_selected = var_subset[,var_selected]
      levels = names(table(var_subset$var_selected))    
    
      ## Calculate KM estimates within each unique level of covariate of interest
      for (j in 1:length(levels)) {
        data_level = subset(var_subset, var_selected==levels[j])
        surv = survfit(as.formula("Surv(follow_up_time, covid_vax) ~ 1"), data = data_level) %>% 
          broom::tidy() %>% 
          filter(estimate > 0) %>%
          mutate(
            N = plyr::round_any(max(n.risk, na.rm=TRUE), rounding_threshold),
            cml.event = cumsum(replace_na(n.event, 0)),
            cml.censor = cumsum(replace_na(n.censor, 0)),
            cml.event = floor_any(cml.event, rounding_threshold),
            cml.censor = floor_any(cml.censor, rounding_threshold),
            n.event = diff(c(0,cml.event)),
            n.censor = diff(c(0,cml.censor)),
            n.risk = N - lag(cml.event + cml.censor,1,0),
            summand = n.event / ((n.risk - n.event) * n.risk),
            surv = cumprod(1 - n.event / n.risk),
            surv.se = surv * sqrt(cumsum(replace_na(summand, 0))),
            llsurv = log(-log(surv)),
            llsurv.se = sqrt((1 / log(surv)^2) * cumsum(summand)),
            surv.ll = exp(-exp(llsurv + qnorm(0.025)*llsurv.se)),
            surv.ul = exp(-exp(llsurv + qnorm(0.975)*llsurv.se)),
            cum.in = 1 - surv,
            cum.in.ul = 1 - surv.ul,
            cum.in.ll = 1 - surv.ll
          )
        
        if (levels[j]!="0" & levels[j]!="1") {
          # Merge KM estimates and N events (floor of 10) with main table
          tab_group$N[tab_group$var==var_list[v] & tab_group$Variable==levels[j]] = surv$N[1]
          tab_group$km_cml_event[tab_group$var==var_list[v] & tab_group$Variable==levels[j]] = max(surv$cml.event)
          tab_group$km_coverage[tab_group$var==var_list[v] & tab_group$Variable==levels[j]] = round(max(surv$cum.in)*100,1)
          tab_group$km_coverage_ll[tab_group$var==var_list[v] & tab_group$Variable==levels[j]] = round(max(surv$cum.in.ll, na.rm=TRUE)*100,1)
          tab_group$km_coverage_ul[tab_group$var==var_list[v] & tab_group$Variable==levels[j]] = round(max(surv$cum.in.ul, na.rm=TRUE)*100,1)
        } else if (levels[j]=="1") {
          tab_group$N[tab_group$var==var_list[v]] = surv$N[1]
          tab_group$km_cml_event[tab_group$var==var_list[v]] = max(surv$cml.event)
          tab_group$km_coverage[tab_group$var==var_list[v]] = round(max(surv$cum.in)*100,1)
          tab_group$km_coverage_ll[tab_group$var==var_list[v]] = round(max(surv$cum.in.ll, na.rm=TRUE)*100,1)
          tab_group$km_coverage_ul[tab_group$var==var_list[v]] = round(max(surv$cum.in.ul, na.rm=TRUE)*100,1)
        }
      } # closes j loop 
    } # closes v loop
  
  tab_group$N[tab_group$N<=redaction_threshold] = "[Redacted]"
  tab_group$km_cml_event[tab_group$km_cml_event<=redaction_threshold] = "[Redacted]"
  tab_group$km_coverage[tab_group$km_cml_event<=redaction_threshold | tab_group$N<=redaction_threshold] = "[Redacted]"
  tab_group$km_coverage_ll[tab_group$km_cml_event<=redaction_threshold | tab_group$N<=redaction_threshold] = "[Redacted]"
  tab_group$km_coverage_ul[tab_group$km_cml_event<=redaction_threshold | tab_group$N<=redaction_threshold] = "[Redacted]"
  
  ## Update column names
  tab_group$imd = imd_subset$imd[1]
  if (imd_subset$ethnicity[1]=="White") {  tab_group$ethnicity = "White" } else {  tab_group$ethnicity = "Non-white" }
  
  if (e == 1) { collated_tab = tab_group } else { collated_tab = rbind(collated_tab, tab_group) }
  
} # closes e loop

collated_tab$imd_alt = as.character(collated_tab$imd)
collated_tab$imd_alt[collated_tab$imd_alt=="1 (most deprived)"] = "1"
collated_tab$imd_alt[collated_tab$imd_alt=="5 (least deprived)"] = "5"

## Save outputs
write_csv(collated_tab, here::here("output", "tables", paste0("intersection_coverage_",outcome_label,".csv")))
