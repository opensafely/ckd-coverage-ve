######################################

# This script 
# - produces a table summarising selected clinical and demographic groups in study cohort (overall and stratified by kidney disease subgroup)
# - saves table as html

######################################

## Import libraries
library('tidyverse')
library('here')
library('glue')
library('gt')
library('gtsummary')
library('plyr')
library('reshape2')

## Import command-line arguments
args <- commandArgs(trailingOnly=TRUE)

## Set input and output pathways for matched/unmatched data - default is unmatched
if(length(args)==0) {
  outcome_label = "dose3"
} else {
  if (args[[1]]=="dose3") {
    outcome_label = "dose3"
  } else if (args[[1]]=="dose4") {
    outcome_label = "dose4"
  } else {
    # Print error if no argument specified
    stop("No outcome specified")
  }
}

## Set rounding and redaction thresholds
rounding_threshold = 5
redaction_threshold = 10

## Create output directory
fs::dir_create(here::here("output", "tables"))

## Import data
if (outcome_label=="dose3") {
  data_cohort <- read_rds(here::here("output", "data", "data_cohort_coverage.rds"))
} else {
  data_cohort <- read_rds(here::here("output", "data", "data_cohort_coverage_dose4.rds"))
}

## Format data
data_cohort <- data_cohort %>%
  mutate(
    N = 1,
    allpop = "All"
  ) 

## Define variables of interest
counts <- data_cohort %>% 
  select(
         N,
         allpop,
         
         ## Demographics
         ageband2,
         sex,
         ethnicity,
         imd,
         rural_urban_group,
         
         # CKD groups and coding
         ckd_5cat,
         chronic_kidney_disease_stages_3_5,
         dialysis,
         kidney_transplant,
         
         # Risk group (occupation/access)
         care_home,
         hscworker,
         housebound,
         endoflife,
         
         ## Risk group (clinical)
         prior_covid_cat,
         immunosuppression, 
         mod_sev_obesity,
         diabetes,
         any_resp_dis,
         chd, 
         cld,
         asplenia,
         cancer,
         haem_cancer,
         other_transplant,
         chronic_neuro_dis_inc_sig_learn_dis,
         sev_mental_ill,
         cev_other,

         ## Other descriptors of interest
         region,
         jcvi_group
         ) 

## Function to clean table names
clean_table_names = function(input_table) {
  # Relabel variables for plotting
  input_table$Variable[input_table$Variable=="care_home"] = "Care home resident"
  input_table$Variable[input_table$Variable=="hscworker"] = "Health/social care worker"
  input_table$Variable[input_table$Variable=="housebound"] = "Housebound"
  input_table$Variable[input_table$Variable=="endoflife"] = "End of life care"
  input_table$Variable[input_table$Variable=="prior_covid_cat"] = "Prior SARS-CoV-2"
  input_table$Variable[input_table$Variable=="immunosuppression"] = "Immunosuppression"
  input_table$Variable[input_table$Variable=="mod_sev_obesity"] = "Moderate/severe obesity"
  input_table$Variable[input_table$Variable=="diabetes"] = "Diabetes"
  input_table$Variable[input_table$Variable=="any_resp_dis"] = "Chronic respiratory disease (inc. asthma)"
  input_table$Variable[input_table$Variable=="chd"] = "Chronic heart disease"
  input_table$Variable[input_table$Variable=="cld"] = "Chronic liver disease"
  input_table$Variable[input_table$Variable=="asplenia"] = "Asplenia"
  input_table$Variable[input_table$Variable=="cancer"] = "Cancer (non-haematologic)"
  input_table$Variable[input_table$Variable=="haem_cancer"] = "Haematologic cancer"
  input_table$Variable[input_table$Variable=="other_transplant"] = "Organ transplant (non-kidney)"
  input_table$Variable[input_table$Variable=="chronic_neuro_dis_inc_sig_learn_dis"] = "Chronic neurological disease (inc. learning disability)"
  input_table$Variable[input_table$Variable=="sev_mental_ill"] = "Severe mental illness"
  input_table$Variable[input_table$Variable=="cev_other"] = "Clinically extremely vulnerable (other)"
  input_table$Variable[input_table$Variable=="chronic_kidney_disease_stages_3_5"] = "CKD3-5 code"
  input_table$Variable[input_table$Variable=="dialysis"] = "Dialysis code"
  input_table$Variable[input_table$Variable=="kidney_transplant"] = "Kidney transplant code"
  
  # Relabel groups for plotting
  input_table$Group[input_table$Group=="ageband2"] = "Age"
  input_table$Group[input_table$Group=="sex"] = "Sex"
  input_table$Group[input_table$Group=="ethnicity"] = "Ethnicity"
  input_table$Group[input_table$Group=="imd"] = "IMD"
  input_table$Group[input_table$Group=="region"] = "Region"
  input_table$Group[input_table$Group=="jcvi_group"] = "JCVI group"
  input_table$Group[input_table$Group=="rural_urban_group"] = "Setting"
  input_table$Group[input_table$Group=="ckd_5cat"] = "Kidney disease subgroup"
  input_table$Group[(input_table$Variable %in% c("CKD3-5 code", "Dialysis code", "Kidney transplant code"))] = "Primary care coding of kidney disease"
  input_table$Group[(input_table$Variable %in% c("Care home resident", "Health/social care worker", "Housebound", "End of life care"))] = "Risk group (occupation/access)"
  input_table$Group[input_table$Variable %in% c("Prior SARS-CoV-2", "Immunosuppression", "Moderate/severe obesity", "Diabetes", "Chronic respiratory disease (inc. asthma)",
                                              "Chronic heart disease", "Chronic liver disease","Asplenia", "Cancer (non-haematologic)", "Haematologic cancer", "Obesity", 
                                              "Organ transplant (non-kidney)", "Chronic neurological disease (inc. learning disability)", "Severe mental illness", 
                                              "Clinically extremely vulnerable (other)")] = "Risk group (clinical)"
  input_table$Group[input_table$Variable=="N"] = "N"
  return(input_table)
}

## Generate full and stratified table
ckd_levels = c("All", "CKD3a", "CKD3b",  "CKD4-5", "RRT (dialysis)", "RRT (Tx)")

## Generate table - full and stratified populations
for (i in 1:length(ckd_levels)) {
  
  if (i == 1) { 
    data_subset = counts
    counts_summary = data_subset %>% 
      tbl_summary(by = allpop)
    counts_summary$inputs$data <- NULL
  } else { 
    data_subset = subset(counts, ckd_5cat==ckd_levels[i]) 
    counts_summary = data_subset %>% 
      select(-ckd_5cat) %>% 
      tbl_summary(by = allpop)
    counts_summary$inputs$data <- NULL
  }

  table1 <- counts_summary$table_body %>%
    select(group = variable, variable = label, count = stat_1) %>%
    separate(count, c("count","perc"), sep = "([(])") %>%
    mutate(count = gsub(" ", "", count)) %>%
    mutate(count = as.numeric(gsub(",", "", count))) %>%
    filter(!(is.na(count))) %>%
    select(-perc)
  table1$percent = round(table1$count/nrow(data_cohort)*100,1)
  colnames(table1) = c("Group", "Variable", "Count", "Percent")
  
  ## Clean names
  table1_clean = clean_table_names(table1)
  
  ## Calculate rounded total
  rounded_n = plyr::round_any(nrow(data_subset), rounding_threshold)
  
  ## Round individual values to rounding threshold
  table1_redacted <- table1_clean %>%
    mutate(Count = plyr::round_any(Count, rounding_threshold))
  table1_redacted$Percent = round(table1_redacted$Count/rounded_n*100,1)
  table1_redacted$Non_Count = rounded_n - table1_redacted$Count
  
  ## Redact any rows with rounded cell counts or non-counts <= redaction threshold 
  table1_redacted$Summary = paste0(prettyNum(table1_redacted$Count, big.mark=",")," (",format(table1_redacted$Percent,nsmall=1),"%)")
  table1_redacted$Summary = gsub(" ", "", table1_redacted$Summary, fixed = TRUE) # Remove spaces generated by decimal formatting
  table1_redacted$Summary = gsub("(", " (", table1_redacted$Summary, fixed = TRUE) # Add first space before (
  table1_redacted$Summary[(table1_redacted$Count>0 & table1_redacted$Count<=redaction_threshold) | (table1_redacted$Non_Count>0 & table1_redacted$Non_Count<=redaction_threshold)] = "[Redacted]"
  table1_redacted$Summary[table1_redacted$Variable=="N"] = prettyNum(table1_redacted$Count[table1_redacted$Variable=="N"], big.mark=",")
  table1_redacted <- table1_redacted %>% select(-Non_Count, -Count, -Percent)
  names(table1_redacted)[3] = ckd_levels[i]
  
  if (i==1) { 
    collated_table = table1_redacted 
    } else { 
    collated_table = collated_table %>% left_join(table1_redacted[,2:3], by ="Variable") 
    collated_table[,i+2][is.na(collated_table[,i+2])] = "--"
    }
}

# Save as html/csv
if (outcome_label=="dose3") {
  gt::gtsave(gt(collated_table), here::here("output","tables", "table1_coverage_redacted_by_CKD.html"))
  write_rds(collated_table, here::here("output", "tables", "table1_coverage_redacted_by_CKD.rds"), compress = "gz")
} else {
  gt::gtsave(gt(collated_table), here::here("output","tables", "table1_coverage_redacted_by_CKD_dose4.html"))
  write_rds(collated_table, here::here("output", "tables", "table1_coverage_redacted_by_CKD_dose4.rds"), compress = "gz")
}