# # # # # # # # # # # # # # # # # # # # #
# This script:
# imports processed data from unmatched VE analysis
# matches AZ and BNT recipients by pre-specified exact/caliper variables
# outputs inclusion/exclusions flowchart data
# # # # # # # # # # # # # # # # # # # # #

# Preliminaries ----

## Import libraries ----
library('tidyverse')
library('lubridate')
library('here')
library('glue')
library('MatchIt')

## Import custom user functions
source(here::here("analysis", "functions.R"))

## Import processed data ----
data_processed <- read_rds(here::here("output", "data", "data_cohort_VE.rds")) 

## Specify exact matching variables
exact_variables <- c(
  "jcvi_group",
  "region",
  "kidney_transplant",
  "dialysis",
  "egfr_cat3",
  "prior_covid_cat",
  NULL #"cev",
)

## Specify caliper variables
caliper_variables <- c(
  age = 5,
  vax2_day = 7,
  tbv1_2 = 7,
  NULL
)

## Define full list of matching variabbles
matching_variables <- c(exact_variables, names(caliper_variables))

## Run matching algorithm
safely_matchit <- purrr::safely(matchit)
matching <-
  safely_matchit(
    formula = vax2_az ~ 1,
    data = data_processed,
    method = "nearest", distance = "glm", # two options redundant since we are using exact + caliper matching
    replace = FALSE,
    estimand = "ATT",
    exact = exact_variables,
    caliper = caliper_variables, std.caliper=FALSE,
    m.order = "data", # data is sorted on (effectively random) patient ID
    ratio = 1L
  )[[1]]

## Print summary of matching process
summary(matching)

## Pick out matched candidates
data_matched <-
  as.data.frame(matching$X) %>%
  add_column(
    match_id = matching$subclass,
    treated = matching$treat,
    patient_id = data_processed$patient_id,
    weight = matching$weights
  ) %>%
  filter(!is.na(match_id)) %>% # remove unmatched people. equivalent to weight != 0
  arrange(match_id, desc(treated)) %>%
  left_join(
    data_processed %>% select(-matching_variables),
    by = "patient_id"
  ) 

## Save data
write_rds(data_matched, here::here("output", "data", "data_cohort_VE_matched.rds"), compress="gz")
write_csv(data_matched, here::here("output", "data", "data_cohort_VE_matched.csv"))

# Define selection criteria ----
data_criteria <- data_processed %>%
  transmute(
    patient_id,
    unmatched = TRUE,
    has_match = patient_id %in% data_matched$patient_id,
    include = (has_match)
  )

## Create and save flow chart
data_flowchart <- data_criteria %>%
  transmute(
    c0 = unmatched,
    c1 = c0 & has_match
  ) %>%
  summarise(
    across(.fns=sum)
  ) %>%
  pivot_longer(
    cols=everything(),
    names_to="criteria",
    values_to="n"
  ) %>%
  mutate(
    n_exclude = lag(n) - n,
    pct_exclude = n_exclude/lag(n),
    pct_all = n / first(n),
    pct_step = n / lag(n),
    crit = str_extract(criteria, "^c\\d+"),
    criteria = fct_case_when(
      crit == "c0" ~ "Unmatched VE cohort",
      crit == "c1" ~ "Matched VE cohort",
      TRUE ~ NA_character_
    )
  )
write_csv(data_flowchart, here::here("output", "tables", "flowchart_VE_matched.csv"))
