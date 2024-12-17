# External Helper Function
process_data <- function(file_path) {
  library(tidyverse)
  library(readxl)
  
  # Load Main Data
  data <- read_excel(file_path, sheet = "merge_data_for_model_input") %>%
    select(1:31)
  
  # Clean Column Names
  data <- data %>%
    janitor::clean_names() # Clean column names to ensure compatibility
  
  # Add Additional Sheets and Columns
  u5_data <- read_excel(file_path, sheet = "U5_mortality_rate") %>%
    janitor::clean_names()
  fu_data <- read_excel(file_path, sheet = "study_end_point_year") %>%
    janitor::clean_names()
  ab_data <- read_excel(file_path, sheet = "antibiotic_from_model_estimate") %>%
    janitor::clean_names()
  opv_data <- read_excel(file_path, sheet = "opv_coverage") %>%
    janitor::clean_names()
  
  data <- data %>%
    mutate(
      u5_mortality_rate = as.numeric(u5_data$u5_mortality_rate),
      fu = fu_data$end_point_in_years,
      fu_cat = case_when(fu <= 2 ~ 1, TRUE ~ 2),
      ab = as.numeric(ab_data$mean),
      opv = as.numeric(opv_data$opv_coverage)
    )
  
  # Remove Unwanted Studies
  data <- data %>%
    mutate(study_code = case_when(
      study_code %in% c(50, 22, 14) & num_of_doses %in% c(2, 3) ~ -99,
      TRUE ~ study_code
    )) %>%
    filter(study_code != -99)
  
  # Adjust StudyCode for Specific Conditions
  data <- data %>%
    mutate(study_code = case_when(
      study_code == 4 & country_code == "GHA" ~ "4A",
      study_code == 4 ~ "4B",
      study_code == 9 & ve_severe_rv == 53.6 ~ "9A",
      study_code == 9 ~ "9B",
      study_code == 78 & efficacy_1_or_effectiveness_0 == 1 ~ "78A",
      study_code == 78 ~ "78B",
      study_code == 81 & country_code == "BGD" ~ "81A",
      study_code == 81 ~ "81B",
      TRUE ~ as.character(study_code)
    ))
  
  # Correct WHO Regions
  data <- data %>%
    mutate(who_regions = case_when(
      country_code %in% c("TWN", "AUS") ~ "WPR",
      TRUE ~ who_regions
    ))
  
  # Adjust Weight Values
  data <- data %>%
    mutate(weight = case_when(
      study_code %in% c("4A", "4B", "4C", "81A", "81B", "61", "64") & weight > 0 ~ 1,
      study_code == 75 & weight == 1 ~ 1 / 6,
      TRUE ~ weight
    ))
  
  # Create Dataset 1: Vaccine Efficacy/Effectiveness
  data1 <- data %>%
    mutate(ve_severe_rv = case_when(
      study_code == 80 ~ 100 * (1 - ((0.5 / 7700) / (23.5 / 5831))),
      study_code == 40 ~ 100 * (1 - ((0 + 0.5) / 380) / ((10 + 0.5) / 381)),
      TRUE ~ ve_severe_rv
    )) %>%
    transmute(
      study_code,
      ve_severe_rv = ve_severe_rv / 100,
      ve_severe_rv_lower95ci = ve_severe_rv_lower95ci / 100,
      ve_severe_rv_upper95ci = ve_severe_rv_upper95ci / 100,
      ve_severe_rv_cidifference = ve_severe_rv_upper95ci - ve_severe_rv_lower95ci,
      logrr_severe_rv = log(1 - ve_severe_rv),
      efficacy1_effectiveness0 = efficacy_1_or_effectiveness_0,
      vax_type = vaccine_type,
      num_of_doses
    )
  
  # Add Weights for Each Country
  country_weights <- data %>%
    select(study_code, alpha_country_code, weight) %>%
    pivot_wider(names_from = alpha_country_code, values_from = weight, values_fill = 0)
  
  data1 <- data1 %>% left_join(country_weights, by = "study_code")
  
  # Check for Errors
  data1 <- data1 %>% mutate(check = rowSums(select(., starts_with("a"))))
  
  # Create Dataset 2: Country-Level Variables
  data2 <- data %>%
    select(
      alpha_country_code, who_regions, avg_diarrhea, avg_popdensity, avg_gdppercapita,
      avg_percentpop_poverty, fu_cat, u5_mortality_rate, opv,
      avg_water_atleastbasic, avg_sanitation_atleastbasic, ab
    ) %>%
    rename(countrycode = alpha_country_code) %>%
    group_by(countrycode) %>%
    summarize(across(everything(), mean, na.rm = TRUE))
  
  list(data1 = data1, data2 = data2)
}

