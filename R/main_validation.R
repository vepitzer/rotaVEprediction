# Main Script for Leave-One-Country-Out and 80-20 Validation

# Load libraries
library(rjags)
library(coda)
library(doParallel)
library(ggplot2)
library(gridExtra)
library(readr)
source("./R/helper_functions_validation.R")  # Load supporting functions

# Set file paths
results_dir <- "./Results"
figures_dir <- "./Figures"

# Load data and prepare inputs
data1<-read_csv("./Data/RV_VE_data1_20250128.csv")
data2<-read_csv("Data/RV_VE_data2_20250128.csv")
model_input <- prep_data_CV(data1,data2)
n <- model_input$n
n_c <- model_input$n_c
n_s <- model_input$n_s
theta_hat <- model_input$theta_hat
sigma2_hat_inv <- model_input$sigma2_hat_inv
z <- model_input$z
w <- model_input$w
super_region <- model_input$super_region

# Set seed for reproducibility
set.seed(3268)

# Run Leave-One-Country-Out validation
cv_results_loco <- run_leave_one_country_out_validation(
  model_input = model_input,
  model_string = build_model_string_loco(),  # Function to create JAGS model
  results_path = file.path(results_dir, "cv_results_LOCO.RDS")
)

# Run 80-20 cross-validation
cv_results_8020 <- run_80_20_validation(
  model_input = model_input,
  model_string = build_model_string_full_model(),  # Function to create JAGS model
  results_path = file.path(results_dir, "cv_results_8020.RDS"),
  k_folds = 5  # Specify the number of folds for 80-20 CV
)

# Post-process results for LOCO
post_process_cv_results(
  cv_results_path = file.path(results_dir, "cv_results_LOCO.RDS"),
  theta_hat = theta_hat,
  output_fig_path = file.path(figures_dir, "Full_Figure_LOCOV.pdf"),
  plot_title = "Leave-One-Country-Out Validation"
)

# Post-process results for 80-20 CV
post_process_cv_results(
  cv_results_path = file.path(results_dir, "cv_results_8020.RDS"),
  theta_hat = theta_hat,
  output_fig_path = file.path(figures_dir, "Full_Figure_8020.pdf"),
  plot_title = "80-20 Cross-Validation"
)
