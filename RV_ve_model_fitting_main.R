# Main Script

# Load necessary library
library(rjags)
library(readxl)
library(matrixStats)
library(ggplot2)
library(HDInterval)
# Load helper functions
source("./R/helper_functions.R")

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

# Prepare the global prediction data
result <- prepare_global_prediction_data(data_pred, data2, save_dir = results_dir)
data_pred <- result$data_pred
data2 <- result$data2
write.csv(data2, file.path("./Data", paste0("RV_VE_data2_", date_suffix, ".csv")), row.names = FALSE)
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
  n.chains = 2,    # Number of MCMC chains (>=2)
  n.iter = 500,   # Number of iterations (adjust based on your requirements)
  thin = 10        # Thinning interval
)

# Save posterior samples to file
saveRDS(posterior_samples, file.path(results_dir, "posterior_samples_full_model_pred_DIC.RDS"))

# Optional: Combine all posterior samples into a single object for further analysis
# posterior_samples_all <- do.call(rbind, posterior_samples)

#######################################################################
# Post-Processing of Results
#######################################################################
# Load required results and posterior samples
data_pred <- readRDS("./Results/data_pred.RDS")
posterior_samples <- readRDS("./Results/posterior_samples_full_model_pred_DIC.RDS")

# Summarize eta estimates
eta_summary <- eta_summary(data_pred, posterior_samples)

# Save eta summaries to file
saveRDS(eta_summary, "./Results/VEs_est_full_model_pred.RDS")

#######################################################################
# Combine Posterior Samples from Chains
#######################################################################
chains <- length(posterior_samples)
final <- posterior_samples[[1]]

for (j in 2:chains) {
  final <- rbind(final, posterior_samples[[1]][[j]])
}

#######################################################################
# Summarize Posterior Distributions of Betas
#######################################################################

# Extract beta coefficients and compute posterior means and credible intervals
beta <- final[, substring(colnames(final), 1, 4) == "beta"]
post_means <- apply(exp(beta), 2, mean)
ci <- t(hdi(exp(beta), credMass = 0.95))
ci <- matrix(ci, ncol = 2)

# Combine posterior means and credible intervals into a summary dataframe
combined <- cbind.data.frame(post_means, ci)
names(combined) <- c("mean", "lcl", "ucl")
rownames(combined) <- c(
  "diarrhea", "popdensity", "GDPpercapita", "water", "sanit",
  "percentpop_poverty", "opv", "FU_cat", "U5_mortality_rate_new", "Ab"
)

# Save betas summary to file
write.csv(combined, "./Results/betas_summary_full_model.csv")

# Visualize posterior distributions of betas
beta_long <- as.data.frame(beta) %>% gather()
ggplot(beta_long, aes(exp(value))) + 
  geom_histogram(bins = 50) + 
  facet_wrap(~key, scales = "free_x")

# Quantile summary of betas
betas_summary <- colQuantiles(exp(beta), probs = c(0.0275, 0.5, 0.975), na.rm = FALSE)
rownames(betas_summary) <- c(
  "diarrhea", "popdensity", "GDPpercapita", "water", "sanit",
  "percentpop_poverty", "opv", "FU_cat", "U5_mortality_rate_new", "Ab"
)

#######################################################################
# Summarize Posterior Distributions of Gamma and Eta
#######################################################################

# Extract gamma and eta values
gamma1 <- final[, substring(colnames(final), 1, 6) == "gamma1"]
gamma0 <- final[, substring(colnames(final), 1, 6) == "gamma0"]
eta0 <- final[, substring(colnames(final), 1, 4) == "eta0"]

# Summarize gamma1
post_means_gamma1 <- mean(exp(gamma1))
ci_gamma1 <- t(hdi(exp(gamma1), credMass = 0.95))
ci_gamma1 <- matrix(ci_gamma1, ncol = 2)
combined_gamma1 <- cbind.data.frame(post_means_gamma1, ci_gamma1)
names(combined_gamma1) <- c("mean", "lcl", "ucl")

# Summarize gamma0
post_means_gamma0 <- mean(exp(gamma0))
ci_gamma0 <- t(hdi(exp(gamma0), credMass = 0.95))
ci_gamma0 <- matrix(ci_gamma0, ncol = 2)
combined_gamma0 <- cbind.data.frame(post_means_gamma0, ci_gamma0)
names(combined_gamma0) <- c("mean", "lcl", "ucl")

# Compute fitted eta1
fitted_eta1 <- 1 - exp(post_means_gamma0 + post_means_gamma1 * eta0)
fitted_eta1 <- colMedians(fitted_eta1)

# Print summaries for gamma0 and gamma1
print(combined_gamma0)
print(combined_gamma1)

