prep_data_CV <- function(data1, data2, interactions = FALSE) {
  #######################################################################
  # Step 1: Outcomes
  #######################################################################
  theta_hat <- data1$logrr_severe_rv
  
  # Calculate variance and its inverse using the Delta method
  se_ve <- (data1$ve_severe_rv_upper95ci - data1$ve_severe_rv_lower95ci) / 3.92
  var_log_rr <- (se_ve / (1.00 - data1$ve_severe_rv))^2
  sigma2_hat_inv <- 1.00 / var_log_rr
  
  # Define binary indicator (efficacy = 0, effectiveness = 1)
  z <- 1 - data1$efficacy1_effectiveness0
  n <- length(theta_hat)
  
  #######################################################################
  # Step 2: Study-Specific Weights
  #######################################################################
  # Extract weights for studies (must check column names manually)
  w <- data1[, 10:ncol(data1)]
  n_c <- ncol(w)  # Number of countries
  
  #######################################################################
  # Step 3: Predictors
  #######################################################################
  # Map WHO regions to numeric IDs
  unique_regions <- unique(data2$WHO_region)
  super_region <- sapply(data2$WHO_region, function(region) which(unique_regions == region))
  n_s <- max(super_region)
  
  # Define predictors
  predictors <- list(
    x1 = data2$diarrhea_imp,
    x2 = data2$LogPopDens_2020_imp,
    x3 = data2$logGDP_percap_2021_imp,
    x4 = data2$Water_atleastbasic_imp,
    x5 = data2$Sanitation_atleastbasic_imp,
    x6 = data2$poverty_imp,
    x7 = data2$avg_opv,
    x8 = model.matrix(~0 + as.factor(data2$fu_cat))[, 2],  # Remove one column
    x9 = data2$U5_mortality,
    x10 = data2$Antibiotic
  )
  
  # If interactions are TRUE, define interaction terms
  if (interactions) {
    stop("Interaction terms generation is not implemented in this simplified version.")
  }
  
  # Scale predictors and combine into a matrix
  predictors_scaled <- lapply(predictors, scale)
  x <- do.call(cbind, predictors_scaled)
  p_x <- ncol(x)
  
  #######################################################################
  # Step 4: Return Results
  #######################################################################
  result <- list(
    "theta_hat"      = theta_hat,
    "sigma2_hat_inv" = sigma2_hat_inv,
    "z"              = z,
    "n"              = n,
    "n_s"            = n_s,
    "n_c"            = n_c,
    "w"              = w,
    "super_region"   = super_region,
    "x"              = x,
    "p_x"            = p_x
  )
  
  return(result)
}

build_model_string_loco <- function() {
  return("
    model {
      for (i in 1:n) {
        theta_hat[i] ~ dnorm(theta[i], sigma2_hat_inv[i])
        theta[i] ~ dnorm(mu_theta[i], sigma2_theta_inv[z[i] + 1])
        mu_theta[i] <- (w[i, ] %*% eta0) * (1 - z[i]) + (w[i, ] %*% eta1) * z[i]
      }
      for (j in 1:n_c) {
        eta0[j] ~ dnorm(mu_eta0[j], sigma2_eta0_inv)
        mu_eta0[j] <- mu + phi[super_region[j]]
        eta1[j] ~ dnorm(mu_eta1[j], sigma2_eta1_inv)
        mu_eta1[j] <- gamma0 + gamma1 * eta0[j]
      }
      for (j in 1:n_s) {
        phi[j] ~ dnorm(0.0, sigma2_phi_inv)
      }
      mu ~ dnorm(0.0, 0.0001)
      gamma0 ~ dnorm(0.0, 0.0001)
      gamma1 ~ dnorm(0.0, 0.0001)
      for(j in 1:2){
   sigma2_theta_inv[j] ~ dgamma(0.01, 
                                0.01)
      }
      sigma2_phi_inv ~ dgamma(0.01, 0.01)
      sigma2_eta0_inv ~ dgamma(0.01, 0.01)
      sigma2_eta1_inv ~ dgamma(0.01, 0.01)
    }
  ")
}

build_model_string_full_model <- function() {
  # Define the model string for the full model (similar structure to above)
  return(build_model_string_loco())
}

run_leave_one_country_out_validation <- function(model_input, model_string, results_path) {
  cv_results <- foreach(cv = 1:model_input$n, .packages = c("rjags", "coda")) %dopar% {
    # Create country-specific theta_hat with one country excluded
    theta_hat_cv <- model_input$theta_hat
    excluded_indices <- which(model_input$w[cv, ] > 0)
    theta_hat_cv[excluded_indices] <- NA
    
    # Run JAGS model
    model_jags <- jags.model(
      textConnection(model_string),
      data = list(
        n = model_input$n,
        n_c = model_input$n_c,
        n_s = model_input$n_s,
        theta_hat = theta_hat_cv,
        sigma2_hat_inv = model_input$sigma2_hat_inv,
        z = model_input$z,
        w = model_input$w,
        super_region = model_input$super_region
      ),
      n.chains = 1
    )
    update(model_jags, n.iter = 5000)  # Burn-in
    posterior_samples <- coda.samples(
      model_jags,
      variable.names = c("theta"),
      n.iter = 10000,
      thin = 10
    )
    return(posterior_samples)
  }
  saveRDS(cv_results, results_path)
  return(cv_results)
}

run_80_20_validation <- function(model_input, model_string, results_path, k_folds) {
  folds <- sample(rep(1:k_folds, length.out = model_input$n))
  cv_results <- list()
  for (fold in 1:k_folds) {
    train_indices <- which(folds != fold)
    test_indices <- which(folds == fold)
    
    theta_hat_train <- model_input$theta_hat
    theta_hat_train[test_indices] <- NA  # Exclude test data
    
    # Run JAGS model
    model_jags <- jags.model(
      textConnection(model_string),
      data = list(
        n = model_input$n,
        n_c = model_input$n_c,
        n_s = model_input$n_s,
        theta_hat = theta_hat_train,
        sigma2_hat_inv = model_input$sigma2_hat_inv,
        z = model_input$z,
        w = model_input$w,
        super_region = model_input$super_region
      ),
      n.chains = 1
    )
    update(model_jags, n.iter = 5000)
    posterior_samples <- coda.samples(
      model_jags,
      variable.names = c("theta"),
      n.iter = 10000,
      thin = 10
    )
    cv_results[[fold]] <- posterior_samples
  }
  saveRDS(cv_results, results_path)
  return(cv_results)
}

post_process_cv_results <- function(cv_results_path, theta_hat, output_fig_path, plot_title) {
  cv_results <- readRDS(cv_results_path)
  
  # Extract posterior medians and credible intervals
  RR_post_median <- sapply(cv_results, function(res) median(exp(res[[1]])))
  RR_ci <- t(sapply(cv_results, function(res) hdi(exp(res[[1]]), credMass = 0.95)))
  
  # Calculate performance metrics
  corr <- cor(exp(theta_hat), RR_post_median)
  rmse <- sqrt(mean((exp(theta_hat) - RR_post_median)^2))
  coverage <- mean((RR_ci[, 1] <= exp(theta_hat)) & (RR_ci[, 2] >= exp(theta_hat)))
  
  # Plot results
  df <- data.frame(VE_obs = 1 - exp(theta_hat), VE_est = 1 - RR_post_median)
  ggplot(df, aes(x = VE_obs, y = VE_est)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    labs(title = plot_title, x = "Observed VE", y = "Predicted VE (Median)") +
    ggsave(output_fig_path)
}
