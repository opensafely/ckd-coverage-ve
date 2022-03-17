######################################

# This script 
# - produces a table with the number of CKD patients in study population in relation to 
#   selected clinical and demographic groups, stratified by primary vaccine series
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

## Create output directory
fs::dir_create(here::here("output", "tables"))

## import metadata ----
#var_labels <- read_rds(here::here("analysis", "variable-labels.rds"))

## Import data
data_cohort <- read_rds(here::here("output", "data", "data_cohort_VE.rds"))

## Format data
data_cohort <- data_cohort %>%
  mutate(N=1,
         time_between_vaccinations1_2 = as.character(cut(tbv1_2,
                                                         breaks = c(0, 42, 70, 98, Inf),
                                                         labels = c("6 weeks or less", "6-10 weeks", "10-14 weeks", "14 weeks or more"),
                                                         right = FALSE)),
         smoking_status = ifelse(is.na(smoking_status), "N&M", smoking_status),
         bpcat = ifelse(bpcat=="High" | bpcat=="Elevated", 1, 0)
  ) %>%
  mutate(smoking_status = ifelse(smoking_status=="S&E", 1, 0),
         time_between_vaccinations1_2 = ifelse(is.na(vax2_date), "Not applicable (no dose given)", time_between_vaccinations1_2)
  )

## baseline variables
counts0 <- data_cohort %>%
  select(vax12_type,
         N,
         ageband2, 
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
         prior_covid_cat,
         time_between_vaccinations1_2 
         ) %>%
  tbl_summary(by = vax12_type) 
counts0$inputs$data <- NULL

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

table1$percent1 = round(table1$count1/sum(data_cohort$vax12_type=="az-az")*100,1)
table1$percent2 = round(table1$count2/sum(data_cohort$vax12_type=="pfizer-pfizer")*100,1)
colnames(table1) = c("Group", "Variable", "Count_az", "Count_pfizer", "Percent_az", "Percent_pfizer")

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
                                    "Organ transplant (any)", "Organ transplant (non-kidney)", "Prior COVID", "Health and social care worker")] = "Other"
table1$Group[table1$Group=="time_between_vaccinations1_2"] = "Time between doses 1 and 2"

# Redaction ----
rounded_n_az = plyr::round_any(sum(data_cohort$vax12_type=="az-az"),5)
rounded_n_pfizer = plyr::round_any(sum(data_cohort$vax12_type=="pfizer-pfizer"),5)

## Round to nearest 5
table1_redacted <- table1 %>%
  mutate(Count_az = plyr::round_any(Count_az, 5)) %>%
  mutate(Count_pfizer = plyr::round_any(Count_pfizer, 5))

## Calculate rounded percentages
table1_redacted$Percent_az = round(table1_redacted$Count_az/rounded_n_az*100,1)
table1_redacted$Percent_pfizer = round(table1_redacted$Count_pfizer/rounded_n_pfizer*100,1)

## Calculate rounded non-counts
table1_redacted$Non_Count_az = rounded_n_az - table1_redacted$Count_az
table1_redacted$Non_Count_pfizer = rounded_n_pfizer - table1_redacted$Count_pfizer

## Redact any rows with rounded cell counts <=10 or within 10 of total population size
table1_redacted$Count_az[table1_redacted$Count_az<=10 | (table1_redacted$Non_Count_az<=10 & table1_redacted$Non_Count_az>0)] = "[Redacted]"
table1_redacted$Count_pfizer[table1_redacted$Count_pfizer<=10 | (table1_redacted$Non_Count_pfizer<=10 & table1_redacted$Non_Count_pfizer>0)] = "[Redacted]"
table1_redacted <- table1_redacted %>% select(-Non_Count_az, -Non_Count_pfizer)

for (i in 1:nrow(table1_redacted)) {
  if (table1_redacted$Count_az[i]!="[Redacted]" & table1_redacted$Group[i]!="N") { 
    table1_redacted$Count_az[i] = paste0(table1_redacted$Count_az[i]," (",table1_redacted$Percent_az[i],"%)") 
  }
  if (table1_redacted$Count_pfizer[i]!="[Redacted]" & table1_redacted$Group[i]!="N") { 
    table1_redacted$Count_pfizer[i] = paste0(table1_redacted$Count_pfizer[i]," (",table1_redacted$Percent_pfizer[i],"%)") 
  }
}
names(table1_redacted)[3:4] <- c("ChAdOx1-S", "BNT162b2")
table1_redacted <- table1_redacted %>% select(-Percent_az, -Percent_pfizer)

# Save as html ----
gt::gtsave(gt(table1_redacted), here::here("output","tables", "table1_VE_redacted.html"))
write_rds(table1_redacted, here::here("output", "tables", "table1_VE_redacted.rds"), compress = "gz")
