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
  colnames(data)[27]<-"avg_gdppercapita"
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
  data1 <-unique(data %>%
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
    )) 
  # Add Weights for Each Country
  countrycode <- data.frame(countrycode=unique(data$alpha_country_code))
  for(i in 1:nrow(countrycode)){
    countrycode.i <- countrycode[i,1]
    country_weights <- data %>%
      dplyr::select(., c(study_code, alpha_country_code, weight)) %>%
      filter(alpha_country_code == countrycode.i) %>%
      rename(!!sprintf(countrycode.i):=weight) %>%
      dplyr::select(.,-alpha_country_code)
    data1 <- left_join(data1, country_weights,by="study_code") %>%
      replace(is.na(.),0)
  }
  data1 <- data1 %>%
    distinct()
  # Check for Errors
  #data1 <- data1 %>% mutate(check = rowSums(select(., where(is.numeric)), na.rm = TRUE))
  
  # Create Dataset 2: Country-Level Variables
  data2 <- data %>%
    group_by(alpha_country_code,who_regions) %>%
    summarize(
      u5_mortality_rate_new = mean(u5_mortality_rate, na.rm = TRUE),
      ab = mean(ab, na.rm = TRUE),
      fu_cat = min(fu_cat, na.rm = TRUE),
      avg_diarrhea = mean(avg_diarrhea, na.rm = TRUE),
      avg_popdensity = mean(avg_popdensity, na.rm = TRUE),
      avg_gdppercapita = mean(avg_gdppercapita, na.rm = TRUE),
      avg_percentpop_poverty = mean(avg_percentpop_poverty, na.rm = TRUE),
      avg_water_atleastbasic = mean(avg_water_atleastbasic, na.rm = TRUE),
      avg_sanitation_atleastbasic = mean(avg_sanitation_atleastbasic, na.rm = TRUE),
      opv = mean(opv, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(countrycode = alpha_country_code)
  
  list(data1 = data1, data2 = data2)
}

prepare_global_prediction_data <- function(data_pred, data2, save_dir = "./Results") {
  # Step 1: Replace invalid values with NA and filter rows
  data_pred <- data_pred %>%
    mutate(across(everything(), ~ replace(., . %in% c("-", "#N/A"), NA))) %>%
    dplyr::filter(!is.na(WHO_region)) # Keep rows with WHO region
  
  # Convert specific columns to numeric
  data_pred <- data_pred %>%
    mutate(across(c(LogPopDens_2020, poverty, diarrhea, avg_opv, 
                    Water_atleastbasic, Sanitation_atleastbasic, 
                    U5_mortality, Antibiotic), as.numeric))
  
  # Log transformations
  data2$avg_popdensity <- log(data2$avg_popdensity)
  data2$avg_gdppercapita <- log(data2$avg_gdppercapita)
  
  # Match and update FU_cat
  # Initialize `fu_cat` column
  if (!"fu_cat" %in% colnames(data_pred)) {
    data_pred$fu_cat <- NA
  }
  matched_indices <- match(data_pred$CountryCode, data2$countrycode)
  data_pred$fu_cat[!is.na(matched_indices)] <- data2$fu_cat[matched_indices[!is.na(matched_indices)]]
  data_pred$fu_cat[is.na(matched_indices)] <- 1
  
  # Step 2: Impute missing values by WHO region
  data_pred <- data_pred %>%
    group_by(WHO_region) %>%
    mutate(
      poverty_imp = ifelse(is.na(poverty), mean(poverty, na.rm = TRUE), poverty),
      GDP_percap_2021_imp = ifelse(is.na(GDP_percap_2021), mean(GDP_percap_2021, na.rm = TRUE), GDP_percap_2021),
      LogPopDens_2020_imp = ifelse(is.na(LogPopDens_2020), mean(LogPopDens_2020, na.rm = TRUE), LogPopDens_2020),
      diarrhea_imp = ifelse(is.na(diarrhea), mean(diarrhea, na.rm = TRUE), diarrhea),
      avg_opv_imp = ifelse(is.na(avg_opv), mean(avg_opv, na.rm = TRUE), avg_opv),
      Water_atleastbasic_imp = ifelse(is.na(Water_atleastbasic), mean(Water_atleastbasic, na.rm = TRUE), Water_atleastbasic),
      Sanitation_atleastbasic_imp = ifelse(is.na(Sanitation_atleastbasic), mean(Sanitation_atleastbasic, na.rm = TRUE), Sanitation_atleastbasic),
      Antibiotic = ifelse(is.na(Antibiotic), mean(Antibiotic, na.rm = TRUE), Antibiotic)
    ) %>%
    ungroup()
  
  data_pred$logGDP_percap_2021_imp <- log(data_pred$GDP_percap_2021_imp)
  
  # Step 3: Merge and clean data
  data2$countrycode[data2$countrycode == "Bel"] <- "BEL"
  colnames(data2)[1] <- "CountryCode"
  data2 <- data2 %>% distinct(CountryCode, .keep_all = TRUE)
  
  delete <- rep(NA, times = nrow(data_pred))
  for (j in 1:nrow(data_pred)) {
    delete[j] <- as.numeric(sum(data2$CountryCode == data_pred$CountryCode[j]) > 0)
  }
  data_pred <- data_pred[delete == 0, ]
  
  # Select and reorder columns
  data_pred <- data_pred[, c(2, 8, 12, 13, 15, 19, 18, 23, 16, 21, 22, 9)]
  colnames(data2) <- colnames(data_pred)
  
  # Combine datasets
  data_pred <- rbind(data2, data_pred)
  
  # Save datasets
  saveRDS(data_pred, file.path(save_dir, "data_pred.RDS"))
  saveRDS(data2, file.path(save_dir, "data2_CV.RDS"))
  
  # Step 4: Compute SD for predictors
  SD_pred <- apply(data_pred[, c(3, 4, 6, 7, 8, 9, 10, 11, 12)], 2, sd)
  
  return(list(data_pred = data_pred, data2 = data2, SD_pred = SD_pred))
}
prep_data_pred <- function(data1, data_pred, interactions = FALSE) {
  #######################################################################
  # Step 1: Outcomes
  #######################################################################
  theta_hat <- data1$logrr_severe_rv
  
  # Calculate variance and its inverse using the delta method
  se_ve <- (data1$ve_severe_rv_upper95ci - data1$ve_severe_rv_lower95ci) / 3.92
  var_log_rr <- (se_ve / (1.00 - data1$ve_severe_rv))^2
  sigma2_hat_inv <- 1.00 / var_log_rr
  
  # Define z based on efficacy/effectiveness
  z <- 1 - data1$efficacy1_effectiveness0
  n <- length(theta_hat)
  
  #######################################################################
  # Step 2: Weights
  #######################################################################
  # Extract weight columns and add placeholder columns (manual adjustment required if structure changes)
  w <- data1[, c(10:ncol(data1))]
  w <- cbind(w, matrix(0, nrow = n, ncol = 139))
  
  #######################################################################
  # Step 3: Predictors
  #######################################################################
  n_c <- nrow(data_pred)
  
  # Map WHO regions to unique numeric IDs
  super_region <- rep(NA, times = n_c)
  unique_regions <- unique(data_pred$WHO_region)
  for (j in seq_along(unique_regions)) {
    super_region[data_pred$WHO_region == unique_regions[j]] <- j
  }
  n_s <- max(super_region, na.rm = TRUE)
  
  # Define predictors
  x1 <- data_pred$diarrhea_imp
  x2 <- data_pred$LogPopDens_2020_imp
  x3 <- data_pred$logGDP_percap_2021_imp
  x4 <- data_pred$Water_atleastbasic_imp
  x5 <- data_pred$Sanitation_atleastbasic_imp
  x6 <- data_pred$poverty_imp
  x7 <- data_pred$avg_opv
  x8 <- model.matrix(~0 + as.factor(data_pred$fu_cat))[, 2] # Remove one column
  x9 <- data_pred$U5_mortality
  x10 <- data_pred$Antibiotic
  
  # Handle interactions if specified
  if (interactions) {
    # Example interaction terms (commented out to save space, add as needed)
    # x11 <- x1 * x1
    # x12 <- x1 * x2
    # Add additional interactions as needed
    
    preds <- cbind(
      scale(x1), scale(x2), scale(x3), scale(x4), scale(x5), 
      scale(x6), scale(x7), scale(x8) # Add interaction terms here
    )
  } else {
    preds <- cbind(
      scale(x1), scale(x2), scale(x3), scale(x4), scale(x5), 
      scale(x6), scale(x7), x8, scale(x9), scale(x10)
    )
  }
  
  #######################################################################
  # Step 4: Principal Component Analysis (Optional)
  #######################################################################
  # Uncomment if PCA is needed
  # pca <- princomp(preds, scores = TRUE)
  # plot(pca, type = "l")
  # summary(pca)
  # preds <- cbind(scale(pca$scores[, 1:3])) # Use top 3 principal components
  
  #######################################################################
  # Step 5: Output the results as a list
  #######################################################################
  result <- list(
    "theta_hat"      = theta_hat,
    "sigma2_hat_inv" = sigma2_hat_inv,
    "z"              = z,
    "n"              = n,
    "n_s"            = n_s,
    "n_c"            = nrow(preds),
    "w"              = w,
    "super_region"   = super_region,
    "x"              = preds,
    "p_x"            = ncol(preds)
  )
  
  return(result)
}
call_jags_full_model <- function(n, n_c, n_s, theta_hat, sigma2_hat_inv, x, p_x, z, w, super_region,n.chains,n.iter,thin) {
  
  ##############
  # Seed for Reproducibility
  ##############
  set.seed(8984)
  
  ##################
  # Statistical Model Definition
  ##################
  model_string <- "
  model {
    # Study-Level Data: Log(Relative Risk)/Log(Odds Ratio) Estimates
    for (i in 1:n) {
      theta_hat[i] ~ dnorm(theta[i], sigma2_hat_inv[i])
      theta[i] ~ dnorm(mu_theta[i], sigma2_theta_inv[z[i] + 1])
      mu_theta[i] <- (w[i,] %*% eta0) * (1 - z[i]) + (w[i,] %*% eta1) * z[i]
    }
    
    # Country-Level Processes
    for (j in 1:n_c) {
      # Efficacy
      eta0[j] ~ dnorm(mu_eta0[j], sigma2_eta0_inv)
      mu_eta0[j] <- mu + x[j,] %*% beta + phi[super_region[j]]
      
      # Effectiveness
      eta1[j] ~ dnorm(mu_eta1[j], sigma2_eta1_inv)
      mu_eta1[j] <- gamma0 + gamma1 * eta0[j]
    }
    
    # Super Region Random Effects
    for (j in 1:n_s) {
      phi[j] ~ dnorm(0.0, sigma2_phi_inv)
    }
    
    # Priors
    for (j in 1:2) {
      sigma2_theta_inv[j] ~ dgamma(0.01, 0.01)
    }
    mu ~ dnorm(0.0, 0.0001)
    for (j in 1:p_x) {
      beta[j] ~ dnorm(0.0, sigma2_beta_inv)
    }
    sigma2_beta_inv ~ dgamma(0.01, 0.01)
    sigma2_phi_inv ~ dgamma(0.01, 0.01)
    gamma0 ~ dnorm(0.0, 0.0001)
    gamma1 ~ dnorm(0.0, 0.0001)
    sigma2_eta0_inv ~ dgamma(0.01, 0.01)
    sigma2_eta1_inv ~ dgamma(0.01, 0.01)
  }
  "
  
  ######################################################################
  # Model Fitting
  ######################################################################
  # Initialize the JAGS model
  model_jags <- jags.model(
    textConnection(model_string),
    data = list(
      'n' = n,
      'n_c' = n_c,
      'n_s' = n_s,
      'theta_hat' = theta_hat,
      'sigma2_hat_inv' = sigma2_hat_inv,
      'x' = x,
      'p_x' = p_x,
      'z' = z,
      'w' = w,
      'super_region' = super_region
    ),
    n.chains = n.chains
  )
  
  # Burn-in period
  update(model_jags, n.iter = n.iter)
  
  # Sample posterior distributions
  posterior_samples <- coda.samples(
    model_jags,
    variable.names = c(
      "theta", "eta0", "eta1", "sigma2_theta_inv", "mu", 
      "beta", "sigma2_beta_inv", "phi", "sigma2_phi_inv", 
      "gamma0", "gamma1", "sigma2_eta0_inv", "sigma2_eta1_inv"
    ),
    thin = thin,
    n.iter = n.iter
  )
  
  # Uncomment below if DIC is needed
  # dic <- dic.samples(model_jags, thin = 10, n.iter = 100000)
  # result <- list("posterior_samples" = posterior_samples, "DIC" = dic)
  
  return(posterior_samples)
}


