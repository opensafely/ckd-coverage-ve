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

## Import custom user functions
#source(here("analysis", "functions.R"))

## Create output directory
fs::dir_create(here::here("output", "tables"))

## Import data
data_cohort <- read_rds(here::here("output", "data", "data_cohort_coverage.rds"))

## Format data
data_cohort <- data_cohort %>%
  mutate(time_between_vaccinations1_2 = as.character(cut(tbv1_2,
                                         breaks = c(0, 42, 70, 98, Inf),
                                         labels = c("6 weeks or less", "6-10 weeks", "10-14 weeks", "14 weeks or more"),
                                         right = FALSE)),
         time_between_vaccinations2_3 = as.character(cut(tbv2_3,
                                            breaks = c(0, 84, 168, Inf),
                                            labels = c("12 weeks or less", "12-24 weeks", "24 weeks or more"),
                                            right = FALSE)),
         smoking_status = ifelse(is.na(smoking_status), "N&M", smoking_status),
         bpcat = ifelse(bpcat=="High" | bpcat=="Elevated", 1, 0),
  ) %>%
  mutate(smoking_status = ifelse(smoking_status=="S&E", 1, 0),
        time_between_vaccinations1_2 = ifelse(is.na(vax2_date), "Not applicable (no dose given)", time_between_vaccinations1_2),
        time_between_vaccinations2_3 = ifelse(is.na(vax3_date), "Not applicable (no dose given)", time_between_vaccinations2_3)
  )

counts0 <- data_cohort %>% 
  select(ageband2, 
         sex,
         ethnicity,
         imd,
         region,
         jcvi_group,
         rural_urban_group,
         chronic_kidney_disease_diagnostic,
         dialysis, 
         kidney_transplant, 
         chronic_kidney_disease_stages_3_5,
         cev,
         care_home,
         hscworker,
         endoflife,
         housebound,
         smoking_status,
         asthma,
         bpcat,
         immunosuppression, 
         chronic_resp_dis, 
         diabetes, 
         cld, 
         chd, 
         asplenia, 
         cancer, 
         haem_cancer,
         obesity, 
         chronic_neuro_dis_inc_sig_learn_dis, 
         sev_mental_ill, 
         non_kidney_transplant,
         multimorb,
         prior_covid_cat,
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

## Relabel variables for plotting
table1$Variable[table1$Variable=="cev"] = "Clinically extremely vulnerable"
table1$Variable[table1$Variable=="care_home"] = "Care home resident"
table1$Variable[table1$Variable=="hscworker"] = "Health and social care worker"
table1$Variable[table1$Variable=="endoflife"] = "End of life care"
table1$Variable[table1$Variable=="housebound"] = "Housebound"
table1$Variable[table1$Variable=="chronic_kidney_disease_diagnostic"] = "CKD diagnostic code"
table1$Variable[table1$Variable=="dialysis"] = "Dialysis"
table1$Variable[table1$Variable=="kidney_transplant"] = "Kidney transplant"
table1$Variable[table1$Variable=="chronic_kidney_disease_stages_3_5"] = "CKD stage 3-5 code"
table1$Variable[table1$Variable=="smoking_status"] = "Current or former smoker"
table1$Variable[table1$Variable=="asthma"] = "Asthma"
table1$Variable[table1$Variable=="bpcat"] = "High or elevated blood pressure"
table1$Variable[table1$Variable=="immunosuppression"] = "Immunosuppression"
table1$Variable[table1$Variable=="chronic_resp_dis"] = "Chronic respiratory disease"
table1$Variable[table1$Variable=="diabetes"] = "Diabetes"
table1$Variable[table1$Variable=="cld"] = "Chronic liver disease"
table1$Variable[table1$Variable=="chd"] = "Chronic heart disease"
table1$Variable[table1$Variable=="asplenia"] = "Asplenia"
table1$Variable[table1$Variable=="cancer"] = "Cancer"
table1$Variable[table1$Variable=="haem_cancer"] = "Haematologic cancer"
table1$Variable[table1$Variable=="obesity"] = "Obesity"
table1$Variable[table1$Variable=="chronic_neuro_dis_inc_sig_learn_dis"] = "Chronic neurological disease (including learning disability)"
table1$Variable[table1$Variable=="sev_mental_ill"] = "Severe mental illness"
table1$Variable[table1$Variable=="non_kidney_transplant"] = "Organ transplant (non-kidney)"
table1$Variable[table1$Variable=="multimorb"] = "Comorbidity count (non-CKD)"
table1$Variable[table1$Variable=="prior_covid_cat"] = "Prior COVID"

# Relabel groups for plotting
# Demography
table1$Group[table1$Group=="ageband2"] = "Age"
table1$Group[table1$Group=="sex"] = "Sex"
table1$Group[table1$Group=="ethnicity"] = "Ethnicity"
table1$Group[table1$Group=="imd"] = "IMD"
table1$Group[table1$Group=="region"] = "Region"
table1$Group[table1$Group=="jcvi_group"] = "JCVI group"
table1$Group[table1$Group=="rural_urban_group"] = "Setting"

# Other
table1$Group[table1$Variable %in% c("CKD diagnostic code", "Dialysis", "Kidney transplant", "CKD stage 3-5 code")] = "Clinical (CKD-related)"
table1$Group[table1$Variable %in% c("Clinically extremely vulnerable", "Care home resident", "Healthcare worker", "End of life care", "Housebound", 
                                    "Current or former smoker", "Asthma", "High or elevated blood pressure", "Shielding", "Immunosuppression", 
                                    "Chronic respiratory disease", "Diabetes", "Chronic liver disease", "Chronic heart disease", "Asplenia", "Cancer",
                                    "Haematologic cancer", "Obesity", "Chronic neurological disease (including learning disability)", "Severe mental illness", 
                                    "Organ transplant (any)", "Organ transplant (non-kidney)", "Comorbidity count (non-CKD)", "Prior COVID")] = "Other"
table1$Group[table1$Group=="time_between_vaccinations1_2"] = "Time between doses 1 and 2"
table1$Group[table1$Group=="time_between_vaccinations2_3"] = "Time between doses 2 and 3"

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
gt::gtsave(gt(table1_redacted), here::here("output","tables", "table1_coverage_redacted.html"))
write_rds(table1_redacted, here::here("output", "tables", "table1_coverage_redacted.rds"), compress = "gz")

