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
library('reshape2')

## Import custom user functions
#source(here("analysis", "functions.R"))

## Create output directory
fs::dir_create(here::here("output", "tables"))

## Import data
data_processed <- read_rds(here::here("output", "data", "data_cohort.rds"))

## Format data
data_processed <- data_processed %>%
  mutate(group = ifelse(care_home_65plus == 1, 1, NA),
         group = ifelse(is.na(group) & ageband == 3, 2, group),
         group = ifelse(is.na(group) & hscworker == 1, 3, group),
         group = ifelse(is.na(group) & ageband == 2, 4, group),
         group = ifelse(is.na(group) & shielded == 1, 5, group),
         group = ifelse(is.na(group) & age >=50 & age <70, 6, group),
         group = ifelse(is.na(group), 7, group),
         group = factor(group),
         ageband3 = cut(
           age,
           breaks = c(16, 50, 60, 70, 80, Inf),
           labels = c("16-50", "50-59", "60-69", "70-79", "80+"),
           right = FALSE)) %>%
  group_by(patient_id) %>%
  ungroup()

## Counts
counts0 <- data_processed %>%
  mutate(time_between_vaccinations1_2 = cut(tbv1_2,
                                         breaks = c(0, 42, 84, Inf),
                                         labels = c("6 weeks or less", "6-12 weeks", "12 weeks or more"),
                                         right = FALSE),
         time_between_vaccinations2_3 = cut(tbv2_3,
                                            breaks = c(0, 84, 168, Inf),
                                            labels = c("12 weeks or less", "12-24 weeks", "24 weeks or more"),
                                            right = FALSE),
         
         smoking_status = ifelse(is.na(smoking_status), "N&M", smoking_status)) %>%
  select(ageband3, 
         sex,
         bmi,
         smoking_status,
         ethnicity,
         imd,
         region,
         asthma,
         asplenia,
         bpcat,
         cancer,
         diabetes,
         chd,
         haem_cancer,
         immunosuppression,
         ckd_diagnostic_code,
         ckd_code_35,
         learning_disability,
         cld,
         chronic_neuro_dis_inc_sig_learn_dis,
         chronic_resp_dis,
         dialysis, 
         sev_mental_ill, 
         organ_transplant,
         time_between_vaccinations1_2,
         time_between_vaccinations2_3,
         prior_covid_cat) %>%
  tbl_summary()

counts0$inputs$data <- NULL

table1 <- counts0$table_body %>%
  select(group = variable, variable = label, count = stat_0) %>%
  separate(count, c("count","perc"), sep = "([(])") %>%
  mutate(count = gsub(" ", "", count),
         count = as.numeric(gsub(",", "", count))) %>%
  filter(!(is.na(count))) %>%
  select(-perc) %>%
  filter(!(group == "prior_covid_cat" & variable == "Unknown"))


table1$percent = round(table1$count/nrow(data_processed)*100,1)
colnames(table1) = c("Group", "Variable", "Count", "Percent")

# Redaction ----

## Redact values < 8
threshold = 8

table1_redacted <- table1 %>%
  mutate(Count = ifelse(Count < threshold, NA, as.numeric(Count)),
         Percent = ifelse(is.na(Count), NA, Percent))
         
## Round to nearest 5
table1_redacted <- table1_redacted %>%
  mutate(Count = plyr::round_any(Count, 5))
table1_redacted$Percent = round(table1_redacted$Count/nrow(data_processed)*100,1)

# Save as html ----
gt::gtsave(gt(table1_redacted), here::here("output","tables", "table1_redacted.html"))

