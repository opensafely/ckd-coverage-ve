######################################

# This script 
# - produces a table with the number of CKD patients in study population in relation to 
#   selected clinical and demographic groups
# - saves table as html

######################################

# Preliminaries ----

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
    # print error if no argument specified
    stop("No outcome specified")
  }
}

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
    allpop = "All",
    # Calculate time between vaccinations 1-2
    time_between_vaccinations1_2 = as.character(cut(tbv1_2,
                                                    breaks = c(0, 42, 56, 70, 84, 98, Inf),
                                                    labels = c("6 weeks or less", "6-8 weeks", "8-10 weeks", "10-12 weeks", "12-14 weeks", "14+ weeks"),
                                                    right = FALSE)),
    # Calculate time between vaccinations 2-3
    time_between_vaccinations2_3 = as.character(cut(tbv2_3,
                                                    breaks = c(0, 84, 168, 252, Inf),
                                                    labels = c("12 weeks or less", "12-24 weeks", "24-36 week", "36+ weeks"),
                                                    right = FALSE)),
    smoking_status = ifelse(is.na(smoking_status), "N&M", smoking_status),
    bpcat = ifelse(bpcat=="High" | bpcat=="Elevated", 1, 0)
  ) %>%
  mutate(
    smoking_status = ifelse(smoking_status=="S&E", 1, 0),
    time_between_vaccinations1_2 = ifelse(is.na(vax2_date), "Not applicable (no dose given)", time_between_vaccinations1_2),
    time_between_vaccinations2_3 = ifelse(is.na(vax3_date), "Not applicable (no dose given)", time_between_vaccinations2_3)
  )

counts <- data_cohort %>% 
  select(
         N,
         allpop,
         ## Factors associated with attitudes/confidence
         ageband2,
         care_home,
         hscworker,
         housebound,
         endoflife,
         rural_urban_group,
         
         ## Factors associated with attitudes/confidence
         sex,
         ethnicity,
         imd,
         prior_covid_cat,
         
         ## Clinical risk group (CKD-related)
         ckd_5cat,
         ckd_7cat,
         #removed: dialysis, kidney_transplant, chronic_kidney_disease_stages_3_5
         
         ## Clinical risk group (non-CKD-related)
         immunosuppression, 
         mod_sev_obesity,
         diabetes,
         any_resp_dis,
         chd, 
         cld,
         asplenia,
         cancer,
         haem_cancer,
         non_kidney_transplant,
         chronic_neuro_dis_inc_sig_learn_dis,
         sev_mental_ill,
         cev_other,
         #removed: smoking_status, asthma, bpcat
         
         ## Other descriptors of interest
         region,
         jcvi_group,
         chronic_kidney_disease_stages_3_5,
         time_between_vaccinations1_2,
         time_between_vaccinations2_3
         #removed: any_ckd_flag
         ) 

## Function to clean table names
clean_table_names = function(input_table) {
  # Relabel variables for plotting
  input_table$Variable[input_table$Variable=="care_home"] = "Care home resident"
  input_table$Variable[input_table$Variable=="hscworker"] = "Health/social care worker"
  input_table$Variable[input_table$Variable=="cev"] = "Clinically extremely vulnerable"
  input_table$Variable[input_table$Variable=="housebound"] = "Housebound"
  input_table$Variable[input_table$Variable=="endoflife"] = "End of life care"
  input_table$Variable[input_table$Variable=="prior_covid_cat"] = "Prior COVID"
  input_table$Variable[input_table$Variable=="immunosuppression"] = "Immunosuppression"
  input_table$Variable[input_table$Variable=="mod_sev_obesity"] = "Moderate/severe obesity"
  input_table$Variable[input_table$Variable=="diabetes"] = "Diabetes"
  input_table$Variable[input_table$Variable=="any_resp_dis"] = "Chronic respiratory disease (inc. asthma)"
  input_table$Variable[input_table$Variable=="chd"] = "Chronic heart disease"
  input_table$Variable[input_table$Variable=="cld"] = "Chronic liver disease"
  input_table$Variable[input_table$Variable=="asplenia"] = "Asplenia"
  input_table$Variable[input_table$Variable=="cancer"] = "Cancer (non-haematologic)"
  input_table$Variable[input_table$Variable=="haem_cancer"] = "Haematologic cancer"
  input_table$Variable[input_table$Variable=="non_kidney_transplant"] = "Organ transplant (non-kidney)"
  input_table$Variable[input_table$Variable=="chronic_neuro_dis_inc_sig_learn_dis"] = "Chronic neurological disease (inc. learning disability)"
  input_table$Variable[input_table$Variable=="sev_mental_ill"] = "Severe mental illness"
  input_table$Variable[input_table$Variable=="cev_other"] = "Clinically extremely vulnerable (other)"
  input_table$Variable[input_table$Variable=="chronic_kidney_disease_stages_3_5"] = "CKD3-5 diagnostic code"
  
  # Relabel groups for plotting
  # Demography
  input_table$Group[input_table$Group=="ageband2"] = "Age"
  input_table$Group[input_table$Group=="sex"] = "Sex"
  input_table$Group[input_table$Group=="ethnicity"] = "Ethnicity"
  input_table$Group[input_table$Group=="imd"] = "IMD"
  input_table$Group[input_table$Group=="region"] = "Region"
  input_table$Group[input_table$Group=="jcvi_group"] = "JCVI group"
  input_table$Group[input_table$Group=="rural_urban_group"] = "Setting"
  input_table$Group[input_table$Group=="ckd_7cat"] = "CKD subgroup"
  input_table$Group[input_table$Group=="time_between_vaccinations1_2"] = "Time between doses 1 and 2"
  input_table$Group[input_table$Group=="time_between_vaccinations2_3"] = "Time between doses 2 and 3"
  
  # Other
  input_table$Group[!(input_table$Group %in% c("Age", "Sex", "Ethnicity", "IMD", "Region", "JCVI group", "Setting", "CKD subgroup", 
                                     "Time between doses 1 and 2", "Time between doses 2 and 3"))] = "Other"
  input_table$Group[input_table$Variable=="N"] = "N"
  return(input_table)
}

## Generate full table
counts_summary = counts %>% 
  select(-ckd_5cat) %>%
  tbl_summary(by = allpop)
counts_summary$inputs$data <- NULL

table1 <- counts_summary$table_body %>%
  select(group = variable, variable = label, count = stat_1) %>%
  separate(count, c("count","perc"), sep = "([(])") %>%
  mutate(count = gsub(" ", "", count)) %>%
  mutate(count = as.numeric(gsub(",", "", count))) %>%
  filter(!(is.na(count))) %>%
  select(-perc)
table1$percent = round(table1$count/nrow(data_cohort)*100,1)
colnames(table1) = c("Group", "Variable", "Count", "Percent")

table1_clean = clean_table_names(table1)

# Redaction ----
rounded_n = plyr::round_any(nrow(data_cohort),5)

## Round to nearest 5
table1_redacted <- table1_clean %>%
  mutate(Count = plyr::round_any(Count, 5))
table1_redacted$Percent = round(table1_redacted$Count/rounded_n*100,1)
table1_redacted$Non_Count = rounded_n - table1_redacted$Count

## Redact any rows with rounded cell counts <=10 or within 10 of total population size
table1_redacted$Count[(table1_redacted$Count>0 & table1_redacted$Count<=10) | (table1_redacted$Non_Count>0 & table1_redacted$Non_Count<=10)] = "[Redacted]"
table1_redacted$Percent[(table1_redacted$Count>0 & table1_redacted$Count<=10) | (table1_redacted$Non_Count>0 & table1_redacted$Non_Count<=10)] = "[Redacted]"
table1_redacted <- table1_redacted %>% select(-Non_Count)

# Save as html ----
if (outcome_label=="dose3") {
  gt::gtsave(gt(table1_redacted), here::here("output","tables", "table1_coverage_redacted.html"))
  write_rds(table1_redacted, here::here("output", "tables", "table1_coverage_redacted.rds"), compress = "gz")
} else {
  gt::gtsave(gt(table1_redacted), here::here("output","tables", "table1_coverage_redacted_dose4.html"))
  write_rds(table1_redacted, here::here("output", "tables", "table1_coverage_redacted_dose4.rds"), compress = "gz")
}


## Set CKD levels for stratified table
ckd_levels = c( "CKD3a (D-T-)", "CKD3b (D-T-)",  "CKD4-5 (D-T-)", "CKD (T+)", "CKD (D+T-)")

## Generate CKD-statified table
for (i in 1:length(ckd_levels)) {
  data_subset = subset(counts, ckd_5cat==ckd_levels[i])
  counts_summary = data_subset %>% 
    select(-ckd_5cat, -ckd_7cat) %>%
    tbl_summary(by = allpop)
  counts_summary$inputs$data <- NULL

  table1 <- counts_summary$table_body %>%
    select(group = variable, variable = label, count = stat_1) %>%
    separate(count, c("count","perc"), sep = "([(])") %>%
    mutate(count = gsub(" ", "", count)) %>%
    mutate(count = as.numeric(gsub(",", "", count))) %>%
    filter(!(is.na(count))) %>%
    select(-perc)
  table1$percent = round(table1$count/nrow(data_cohort)*100,1)
  colnames(table1) = c("Group", "Variable", "Count", "Percent")
  
  # Clean names
  table1_clean = clean_table_names(table1)
  
  # Redaction ----
  rounded_n = plyr::round_any(nrow(data_subset),5)
  
  ## Round to nearest 5
  table1_redacted <- table1_clean %>%
    mutate(Count = plyr::round_any(Count, 5))
  table1_redacted$Percent = round(table1_redacted$Count/rounded_n*100,1)
  table1_redacted$Non_Count = rounded_n - table1_redacted$Count
  
  ## Redact any rows with rounded cell counts <=10 or within 10 of total population size
  table1_redacted$Summary =  paste0(prettyNum(table1_redacted$Count, big.mark=",")," (",table1_redacted$Percent,"%)")
  table1_redacted$Summary[(table1_redacted$Count>0 & table1_redacted$Count<=10) | (table1_redacted$Non_Count>0 & table1_redacted$Non_Count<=10)] = "[Redacted]"
  table1_redacted$Summary[table1_redacted$Variable=="N"] = prettyNum(table1_redacted$Count[table1_redacted$Variable=="N"], big.mark=",")
  
  table1_redacted <- table1_redacted %>%
    select(-Non_Count, -Count, -Percent)
  names(table1_redacted)[3] = ckd_levels[i]
  
  if (i==1) { collated_table = table1_redacted } else { collated_table[,2+i] = table1_redacted[,3] }
}

# Save as html ----
if (outcome_label=="dose3") {
  gt::gtsave(gt(collated_table), here::here("output","tables", "table1_coverage_redacted_by_CKD.html"))
  write_rds(collated_table, here::here("output", "tables", "table1_coverage_redacted_by_CKD.rds"), compress = "gz")
} else {
  gt::gtsave(gt(collated_table), here::here("output","tables", "table1_coverage_redacted_dose4_by_CKD.html"))
  write_rds(collated_table, here::here("output", "tables", "table1_coverage_redacted_dose4_by_CKD.rds"), compress = "gz")
}
