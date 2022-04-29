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

## Import command-line arguments (specifying whether or not to run matched analysis)
args <- commandArgs(trailingOnly=TRUE)

## Set input and output pathways for matched/unmatched data - default is unmatched
if(length(args)==0){
  # default (unmatched) file names
  input_name = "data_cohort_VE.rds"
  output_html = "table1_VE_redacted.html"
  output_rds = "table1_VE_redacted.rds"
} else {
  if (args[[1]]=="unmatched") { 
    # unmatched file names
    input_name = "data_cohort_VE.rds"
    output_html = "table1_VE_redacted.html"
    output_rds = "table1_VE_redacted.rds"
  } else if (args[[1]]=="matched") {
    input_name = "data_cohort_VE_matched.rds"
    output_html = "table1_VE_matched_redacted.html"
    output_rds = "table1_VE_matched_redacted.rds"
  } else {
    # print error if no argument specified
    print("No matching argument specified")
  }
}

## Create output directory
fs::dir_create(here::here("output", "tables"))

## Import data
data_cohort <- read_rds(here::here("output", "data", input_name))

## Format data
data_cohort <- data_cohort %>%
  mutate(
    N = 1,
    #allpop = "All",
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

## baseline variables
counts0 <- data_cohort %>%
  select(vax12_type,
         N,
         #allpop,
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
         #removed: dialysis, kidney_transplant
         
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
table1$Variable[table1$Variable=="chronic_kidney_disease_stages_3_5"] = "CKD3-5 diagnostic code"

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
table1$Group[table1$Variable=="N"] = "N"

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
gt::gtsave(gt(table1_redacted), here::here("output","tables", output_html))
write_rds(table1_redacted, here::here("output", "tables", output_rds), compress = "gz")
