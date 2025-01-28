# Main Script

# Load necessary library
library(rjags)
library(readxl)

# Define file paths
data_pred_path <- "./Data/data_pred_U5MR_antibiotic.xlsx"
processed_data_path <- "./Data/updated_data_with_doi_and_antibiotic_for_model_opv_cov.xlsx"
results_dir <- "./Results"
date_suffix <- format(Sys.Date(), "%Y%m%d") # Automatically generate today's date for file naming

# Load data
data_pred <- read_excel(data_pred_path, sheet = "data_globalpredictors")
processed_data <- process_data(processed_data_path)
data1 <- processed_data$data1
data2 <- processed_data$data2

# Export processed data
write.csv(data1, file.path("./Data", paste0("RV_VE_data1_", date_suffix, ".csv")), row.names = FALSE)
write.csv(data2, file.path("./Data", paste0("RV_VE_data2_", date_suffix, ".csv")), row.names = FALSE)

# Prepare the global prediction data
result <- prepare_global_prediction_data(data_pred, data2, save_dir = results_dir)
data_pred <- result$data_pred
data2 <- result$data2
SD_pred <- result$SD_pred

#######################################################################
# Manipulate data to obtain input for model
#######################################################################
# Prepare input data for JAGS model
# NOTE: Ensure column names for weights in `prep_data_pred` are checked and correct.
model_input <- prep_data_pred(data1, data_pred)

#######################################################################
# Full model call using JAGS
#######################################################################
# Call the full model function with the prepared model input
posterior_samples <- call_jags_full_model(
  n = model_input$n,
  n_c = model_input$n_c,
  n_s = model_input$n_s,
  theta_hat = model_input$theta_hat,
  sigma2_hat_inv = model_input$sigma2_hat_inv,
  x = model_input$x,
  p_x = model_input$p_x,
  z = model_input$z,
  w = model_input$w,
  super_region = model_input$super_region,
  n.chains = 2,    # Number of MCMC chains
  n.iter = 5000,   # Number of iterations (adjust based on your requirements)
  thin = 10        # Thinning interval
)

# Save posterior samples to file
saveRDS(posterior_samples, file.path(results_dir, "posterior_samples_full_model_pred_DIC.RDS"))

# Optional: Combine all posterior samples into a single object for further analysis
# posterior_samples_all <- do.call(rbind, posterior_samples)

#######################################################################
# End of script
#######################################################################
