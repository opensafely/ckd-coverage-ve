######################################

# This script 
# - produces a table summarising selected clinical and demographic groups in study cohort, stratified by primary vaccine product
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

## Import command-line arguments (specifying whether or not to run matched analysis)
args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # default (unmatched) file names
  matching_status = "unmatched"
  subgroup = "all"
} else {
  matching_status = args[[1]] # can be unmatched or matched
  subgroup = args[[2]] # can be all, CKD, dialysis, or transplant
}

## Import data
if (matching_status=="unmatched") { 
  data_cohort <- read_rds(here::here("output", "data", "data_cohort_VE.rds"))
  } else { 
  data_cohort <- read_rds(here::here("output", "data", "data_cohort_VE_matched.rds"))
  }

## Specify output path names
output_html = paste0("table1_VE_redacted_",matching_status,"_",subgroup,".html")
output_rds = paste0("table1_VE_redacted_",matching_status,"_",subgroup,".rds")

## Set rounding and redaction thresholds
rounding_threshold = 5
redaction_threshold = 10

## Create output directory
fs::dir_create(here::here("output", "tables"))

## Select subset
if (subgroup=="all") {
  data_cohort = data_cohort
} else if (subgroup=="CKD3") {
  data_cohort = subset(data_cohort, ckd_3cat == "CKD3")
} else if (subgroup=="CKD4-5") {
  data_cohort = subset(data_cohort, ckd_3cat == "CKD4-5")
} else if (subgroup=="RRT") {
  data_cohort = subset(data_cohort, ckd_3cat == "RRT (any)")
}

## Format data
data_cohort <- data_cohort %>%
  mutate(
    N = 1,
  ) 

## Baseline variables
counts0 <- data_cohort %>%
  select(vax12_type,
         N,
         
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
         
         ## Risk group (clinical)
         prior_covid_cat,
         immunosuppression, 
         sev_obesity,
         diabetes,
         any_resp_dis,
         chd, 
         cld,
         asplenia,
         haem_cancer,
         other_transplant,
         chronic_neuro_dis_inc_sig_learn_dis,
         learning_disability,
         sev_mental_ill,
         cev,
         
         ## Summary comorbidity metrics
         any_immunosuppression,
         multimorb,
         
         ## Other descriptors of interest
         region,
         jcvi_group
         ) %>%
  tbl_summary(by = vax12_type) 
counts0$inputs$data <- NULL

## Create table 1
table1 <- counts0$table_body %>%
  select(group = variable, variable = label, count1 = stat_1, count2 = stat_2) %>%
  separate(count1, c("count1","perc1"), sep = "([(])") %>%
  mutate(count1 = gsub(" ", "", count1)) %>%
  mutate(count1 = as.numeric(gsub(",", "", count1))) %>%
  filter(!(is.na(count1))) %>%
  separate(count2, c("count2","perc2"), sep = "([(])") %>%
  mutate(count2 = gsub(" ", "", count2)) %>%
  mutate(count2 = as.numeric(gsub(",", "", count2))) %>%
  filter(!(is.na(count2))) %>%
  select(-perc1, -perc2)
colnames(table1) = c("Group", "Variable", "Count_az", "Count_pfizer")

# Relabel variables for plotting
table1$Variable[table1$Variable=="prior_covid_cat"] = "Prior SARS-CoV-2"
table1$Variable[table1$Variable=="immunosuppression"] = "Immunosuppression"
table1$Variable[table1$Variable=="sev_obesity"] = "Severe obesity"
table1$Variable[table1$Variable=="diabetes"] = "Diabetes"
table1$Variable[table1$Variable=="any_resp_dis"] = "Chronic respiratory disease (inc. asthma)"
table1$Variable[table1$Variable=="chd"] = "Chronic heart disease"
table1$Variable[table1$Variable=="cld"] = "Chronic liver disease"
table1$Variable[table1$Variable=="asplenia"] = "Asplenia"
table1$Variable[table1$Variable=="haem_cancer"] = "Haematologic cancer"
table1$Variable[table1$Variable=="other_transplant"] = "Organ transplant (non-kidney)"
table1$Variable[table1$Variable=="chronic_neuro_dis_inc_sig_learn_dis"] = "Chronic neurological disease"
table1$Variable[table1$Variable=="learning_disability"] = "Learning disability"
table1$Variable[table1$Variable=="sev_mental_ill"] = "Severe mental illness"
table1$Variable[table1$Variable=="cev"] = "Clinically extremely vulnerable"
table1$Variable[table1$Variable=="chronic_kidney_disease_stages_3_5"] = "CKD3-5 code"
table1$Variable[table1$Variable=="dialysis"] = "Dialysis code"
table1$Variable[table1$Variable=="kidney_transplant"] = "Kidney transplant code"
table1$Variable[table1$Variable=="any_immunosuppression"] = "Any immunosuppression"

# Relabel groups for plotting
table1$Group[table1$Group=="ageband2"] = "Age"
table1$Group[table1$Group=="sex"] = "Sex"
table1$Group[table1$Group=="ethnicity"] = "Ethnicity"
table1$Group[table1$Group=="imd"] = "IMD"
table1$Group[table1$Group=="region"] = "Region"
table1$Group[table1$Group=="jcvi_group"] = "JCVI group"
table1$Group[table1$Group=="rural_urban_group"] = "Setting"
table1$Group[table1$Group=="ckd_5cat"] = "Kidney disease subgroup"
table1$Group[table1$Group=="multimorb"] = "Comorbidity count"
table1$Group[(table1$Variable %in% c("CKD3-5 code", "Dialysis code", "Kidney transplant code"))] = "Primary care coding of kidney disease"
table1$Group[(table1$Variable %in% c("Care home resident", "Health/social care worker", "Housebound", "End of life care"))] = "Risk group (occupation/access)"
table1$Group[table1$Variable %in% c("Prior SARS-CoV-2", "Immunosuppression", "Any immunosuppression", "Severe obesity", "Diabetes", "Chronic respiratory disease (inc. asthma)",
                                              "Chronic heart disease", "Chronic liver disease","Asplenia", "Haematologic cancer", "Obesity", 
                                              "Organ transplant (non-kidney)", "Chronic neurological disease", "Learning disability", "Severe mental illness", 
                                              "Clinically extremely vulnerable")] = "Risk group (clinical)"
table1$Group[table1$Variable=="N"] = "N"

## Calculate rounded total
rounded_n_az = plyr::round_any(sum(data_cohort$vax12_type=="az-az"), rounding_threshold)
rounded_n_pfizer = plyr::round_any(sum(data_cohort$vax12_type=="pfizer-pfizer"), rounding_threshold)

## Round individual values to rounding threshold
table1_redacted <- table1 %>%
  mutate(Count_az = plyr::round_any(Count_az, 5)) %>%
  mutate(Count_pfizer = plyr::round_any(Count_pfizer, 5))

## Calculate rounded percentages
table1_redacted$Percent_az = round(table1_redacted$Count_az/rounded_n_az*100,1)
table1_redacted$Percent_pfizer = round(table1_redacted$Count_pfizer/rounded_n_pfizer*100,1)

## Calculate rounded non-counts
table1_redacted$Non_Count_az = rounded_n_az - table1_redacted$Count_az
table1_redacted$Non_Count_pfizer = rounded_n_pfizer - table1_redacted$Count_pfizer

## Redact any rows with rounded cell counts or non-counts <= redaction threshold 
table1_redacted$Summary_az = paste0(prettyNum(table1_redacted$Count_az, big.mark=",")," (",format(table1_redacted$Percent_az,nsmall=1),"%)")
table1_redacted$Summary_az = gsub(" ", "", table1_redacted$Summary_az, fixed = TRUE) # Remove spaces generated by decimal formatting
table1_redacted$Summary_az = gsub("(", " (", table1_redacted$Summary_az, fixed = TRUE) # Add first space before (
table1_redacted$Summary_az[(table1_redacted$Count_az>0 & table1_redacted$Count_az<=redaction_threshold) | (table1_redacted$Non_Count_az>0 & table1_redacted$Non_Count_az<=redaction_threshold)] = "[Redacted]"
table1_redacted$Summary_az[table1_redacted$Variable=="N"] = prettyNum(table1_redacted$Count_az[table1_redacted$Variable=="N"], big.mark=",")

table1_redacted$Summary_pfizer = paste0(prettyNum(table1_redacted$Count_pfizer, big.mark=",")," (",format(table1_redacted$Percent_pfizer,nsmall=1),"%)")
table1_redacted$Summary_pfizer = gsub(" ", "", table1_redacted$Summary_pfizer, fixed = TRUE) # Remove spaces generated by decimal formatting
table1_redacted$Summary_pfizer = gsub("(", " (", table1_redacted$Summary_pfizer, fixed = TRUE) # Add first space before (
table1_redacted$Summary_pfizer[(table1_redacted$Count_pfizer>0 & table1_redacted$Count_pfizer<=redaction_threshold) | (table1_redacted$Non_Count_pfizer>0 & table1_redacted$Non_Count_pfizer<=redaction_threshold)] = "[Redacted]"
table1_redacted$Summary_pfizer[table1_redacted$Variable=="N"] = prettyNum(table1_redacted$Count_pfizer[table1_redacted$Variable=="N"], big.mark=",")

## Drop unnecessary columns
table1_redacted <- table1_redacted %>% select(-Non_Count_az, -Count_az, -Percent_az, -Non_Count_pfizer, -Count_pfizer, -Percent_pfizer)

## Relabel vaccine groups
names(table1_redacted)[3:4] <- c("ChAdOx1-S", "BNT162b2")

## Save as html/rds
gt::gtsave(gt(table1_redacted), here::here("output","tables", output_html))
write_rds(table1_redacted, here::here("output", "tables", output_rds), compress = "gz")
