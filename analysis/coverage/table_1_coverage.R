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
  outcome_label = "dose2"
} else {
  if (args[[1]]=="dose2") {
    outcome_label = "dose2"
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
if (outcome_label=="dose2") {
  data_cohort <- read_rds(here::here("output", "data", "data_cohort_coverage.rds"))
} else {
  data_cohort <- read_rds(here::here("output", "data", "data_cohort_coverage_dose4.rds"))
}

## Format data
data_cohort <- data_cohort %>%
  mutate(
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
    bpcat = ifelse(bpcat=="High" | bpcat=="Elevated", 1, 0),
  ) %>%
  mutate(
    smoking_status = ifelse(smoking_status=="S&E", 1, 0),
    time_between_vaccinations1_2 = ifelse(is.na(vax2_date), "Not applicable (no dose given)", time_between_vaccinations1_2),
    time_between_vaccinations2_3 = ifelse(is.na(vax3_date), "Not applicable (no dose given)", time_between_vaccinations2_3)
  )

counts0 <- data_cohort %>% 
  select(
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
         ckd_7cat,
         #removed: dialysis, kidney_transplant, chronic_kidney_disease_stages_3_5,
         
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
         #removed: smoking_status, asthma, bpcat,
         
         ## Other descriptors of interest
         region,
         jcvi_group,
         any_ckd_flag,
         time_between_vaccinations1_2,
         time_between_vaccinations2_3
         ) %>%
  tbl_summary()
counts0$inputs$data <- NULL

table1 <- counts0$table_body %>%
  select(group = variable, variable = label, count = stat_0) %>%
  separate(count, c("count","perc"), sep = "([(])") %>%
  mutate(count = gsub(" ", "", count)) %>%
  mutate(count = as.numeric(gsub(",", "", count))) %>%
  filter(!(is.na(count))) %>%
  select(-perc)

table1$percent = round(table1$count/nrow(data_cohort)*100,1)
colnames(table1) = c("Group", "Variable", "Count", "Percent")

# Relabel variables for plotting
table1$Variable[table1$Variable=="care_home"] = "Care home resident"
table1$Variable[table1$Variable=="hscworker"] = "Health/social care worker"
table1$Variable[table1$Variable=="cev"] = "Clinically extremely vulnerable"
table1$Variable[table1$Variable=="housebound"] = "Housebound"
table1$Variable[table1$Variable=="endoflife"] = "End of life care"
table1$Variable[table1$Variable=="prior_covid_cat"] = "Prior COVID"
table1$Variable[table1$Variable=="immunosuppression"] = "Immunosuppression"
table1$Variable[table1$Variable=="mod_sev_obesity"] = "Moderate/severe obesity"
table1$Variable[table1$Variable=="diabetes"] = "Diabetes"
table1$Variable[table1$Variable=="any_resp_dis"] = "Chronic respiratory disease (inc. asthma)"
table1$Variable[table1$Variable=="chd"] = "Chronic heart disease"
table1$Variable[table1$Variable=="cld"] = "Chronic liver disease"
table1$Variable[table1$Variable=="asplenia"] = "Asplenia"
table1$Variable[table1$Variable=="cancer"] = "Cancer (non-haematologic)"
table1$Variable[table1$Variable=="haem_cancer"] = "Haematologic cancer"
table1$Variable[table1$Variable=="non_kidney_transplant"] = "Organ transplant (non-kidney)"
table1$Variable[table1$Variable=="chronic_neuro_dis_inc_sig_learn_dis"] = "Chronic neurological disease (inc. learning disability)"
table1$Variable[table1$Variable=="sev_mental_ill"] = "Severe mental illness"
table1$Variable[table1$Variable=="cev_other"] = "Clinically extremely vulnerable (other)"
table1$Variable[table1$Variable=="any_ckd_flag"] = "CKD diagnostic code"

# Relabel groups for plotting
# Demography
table1$Group[table1$Group=="ageband2"] = "Age"
table1$Group[table1$Group=="sex"] = "Sex"
table1$Group[table1$Group=="ethnicity"] = "Ethnicity"
table1$Group[table1$Group=="imd"] = "IMD"
table1$Group[table1$Group=="region"] = "Region"
table1$Group[table1$Group=="jcvi_group"] = "JCVI group"
table1$Group[table1$Group=="rural_urban_group"] = "Setting"
table1$Group[table1$Group=="ckd_7cat"] = "CKD subgroup"
table1$Group[table1$Group=="time_between_vaccinations1_2"] = "Time between doses 1 and 2"
table1$Group[table1$Group=="time_between_vaccinations2_3"] = "Time between doses 2 and 3"

# Other
table1$Group[!(table1$Group %in% c("Age", "Sex", "Ethnicity", "IMD", "Region", "JCVI group", "Setting", "CKD subgroup", 
                                   "Time between doses 1 and 2", "Time between doses 2 and 3"))] = "Other"

# Redaction ----
rounded_n = plyr::round_any(nrow(data_cohort),5)

## Round to nearest 5
table1_redacted <- table1 %>%
  mutate(Count = plyr::round_any(Count, 5))
table1_redacted$Percent = round(table1_redacted$Count/rounded_n*100,1)
table1_redacted$Non_Count = rounded_n - table1_redacted$Count

## Redact any rows with rounded cell counts <=10 or within 10 of total population size
table1_redacted$Count[table1_redacted$Count<=10 | table1_redacted$Non_Count<=10] = "[Redacted]"
table1_redacted$Percent[table1_redacted$Count<=10 | table1_redacted$Non_Count<=10] = "[Redacted]"
table1_redacted <- table1_redacted %>% select(-Non_Count)

# Save as html ----
if (outcome_label=="dose2") {
  gt::gtsave(gt(table1_redacted), here::here("output","tables", "table1_coverage_redacted.html"))
  write_rds(table1_redacted, here::here("output", "tables", "table1_coverage_redacted.rds"), compress = "gz")
} else {
  gt::gtsave(gt(table1_redacted), here::here("output","tables", "table1_coverage_redacted_dose4.html"))
  write_rds(table1_redacted, here::here("output", "tables", "table1_coverage_redacted_dose4.rds"), compress = "gz")
}
