library(dplyr)
library(knitr)
library(MASS)
library(ggplot2)

# Load the data
nvda_data <- read.csv("z5417441_Vishnu_Yannam_NVDA.csv")
amd_data <- read.csv("z5417441_Vishnu_Yannam_AMD.csv")

# Select relevant columns and calculate daily log returns
nvda_data_clean <- nvda_data %>% select(Date, Adj.Close) %>%
  mutate(LogReturn = c(NA, diff(log(Adj.Close)))) %>%
  filter(!is.na(LogReturn)) %>%
  mutate(Loss = 1 - exp(LogReturn))

amd_data_clean <- amd_data %>% select(Date, Adj.Close) %>%
  mutate(LogReturn = c(NA, diff(log(Adj.Close)))) %>%
  filter(!is.na(LogReturn)) %>%
  mutate(Loss = 1 - exp(LogReturn))

# Define the alpha levels
alpha_levels <- c(0.95, 0.99, 0.995)

# Function to calculate VaR and ES
VaR_ES <- function(losses, alpha) {
  VaR <- quantile(losses, alpha)
  ES <- mean(losses[losses > VaR])
  return(c(VaR = VaR, ES = ES))
}

# Calculate VaR and ES for NVDA and AMD
nvda_VaR_ES <- sapply(alpha_levels, function(alpha) VaR_ES(nvda_data_clean$Loss, alpha))
amd_VaR_ES <- sapply(alpha_levels, function(alpha) VaR_ES(amd_data_clean$Loss, alpha))

# Empirical Portfolio VaR and ES
portfolio_loss <- 0.5 * (nvda_data_clean$Loss + amd_data_clean$Loss)
portfolio_VaR_ES <- sapply(alpha_levels, function(alpha) VaR_ES(portfolio_loss, alpha))

#----------------------------TASK 2-------------------------------------------#
n_samples <- 5000

# Function to simulate samples from the independence copula
simulate_independence_copula <- function(n) {
  u_samples <- matrix(runif(n * 2), ncol = 2)
  return(u_samples)
}

# Simulate samples from the independence copula
set.seed(5417441)
u_independence <- simulate_independence_copula(n_samples)

# Convert to original scale using empirical distribution functions
x_tilde_nvda_independence <- quantile(nvda_data_clean$LogReturn, u_independence[, 1])
x_tilde_amd_independence <- quantile(amd_data_clean$LogReturn, u_independence[, 2])

# Calculate potential losses
l_tilde_nvda_independence <- 1 - exp(x_tilde_nvda_independence)
l_tilde_amd_independence <- 1 - exp(x_tilde_amd_independence)

# Portfolio loss
portfolio_loss_independence <- 0.5 * (l_tilde_nvda_independence + l_tilde_amd_independence)

# Calculate VaR and ES for the portfolio using independence copula
portfolio_VaR_ES_independence <- sapply(alpha_levels, function(alpha) VaR_ES(portfolio_loss_independence, alpha))

#----------------------------TASK 3-------------------------------------------#
# Calculate the probability of concordance and discordance
calculate_concordance_discordance <- function(x, y) {
  n <- length(x)
  concordant_pairs <- 0
  discordant_pairs <- 0
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if ((x[i] - x[j]) * (y[i] - y[j]) > 0) {
        concordant_pairs <- concordant_pairs + 1
      } else if ((x[i] - x[j]) * (y[i] - y[j]) < 0) {
        discordant_pairs <- discordant_pairs + 1
      }
    }
  }
  
  total_pairs <- 0.5 * n * (n - 1)
  prob_concordance <- concordant_pairs / total_pairs
  prob_discordance <- discordant_pairs / total_pairs
  
  return(list(prob_concordance = prob_concordance, prob_discordance = prob_discordance))
}

# Calculate probabilities of concordance and discordance
concordance_discordance_probs <- calculate_concordance_discordance(nvda_data_clean$Loss, amd_data_clean$Loss)

# Extract the probabilities
prob_concordance <- concordance_discordance_probs$prob_concordance
prob_discordance <- concordance_discordance_probs$prob_discordance

# Calculate Kendall's Tau using the probabilities
kendall_tau <- prob_concordance - prob_discordance

# Function to calculate Spearman's Rho using empirical distribution function
compute_spearman_rho <- function(x, y) {
  rank_x <- rank(x)
  rank_y <- rank(y)
  rho <- cor(rank_x, rank_y)
  return(rho)
}

# Calculate Spearman's Rho using the computed losses
spearman_rho <- compute_spearman_rho(nvda_data_clean$Loss, amd_data_clean$Loss)

#----------------------------TASK 4-------------------------------------------#
# Define the sample size
n_samples <- 5000

# Calculate the sample correlation coefficient (rho) for the returns
rho <- cor(nvda_data_clean$LogReturn, amd_data_clean$LogReturn)

# Function to apply the Gaussian copula model
gaussian_copula <- function(n_samples, rho, x1, x2) {
  # Simulate multivariate normal random deviates
  set.seed(123)
  y <- mvrnorm(n = n_samples, mu = c(0, 0), Sigma = matrix(c(1, rho, rho, 1), ncol = 2))
  
  # Simulate uniform samples from normal deviates
  u <- pnorm(y, mean = 0, sd = 1)
  
  # Transform the uniform samples to empirical marginals
  x_tilde_nvda <- quantile(nvda_data_clean$LogReturn, u[, 1])
  x_tilde_amd <- quantile(amd_data_clean$LogReturn, u[, 2])
  
  # Calculate the potential loss
  l_tilde_nvda <- 1 - exp(x_tilde_nvda)
  l_tilde_amd <- 1 - exp(x_tilde_amd)
  
  # Calculate the loss for the portfolio
  portfolio_loss_tilde <- 0.5 * (l_tilde_nvda + l_tilde_amd)
  
  # Calculate VaR and ES for the portfolio using Gaussian copula
  VaR_ES_gaussian <- sapply(alpha_levels, function(alpha) VaR_ES(portfolio_loss_tilde, alpha))
  
  # Return the results
  list(VaR = VaR_ES_gaussian[1, ], ES = VaR_ES_gaussian[2, ], l_tilde_nvda = l_tilde_nvda, l_tilde_amd = l_tilde_amd)
}


# Apply the Gaussian copula model
results <- gaussian_copula(n_samples, rho, nvda_data_clean$LogReturn, amd_data_clean$LogReturn)
l_tilde_nvda <- results$l_tilde_nvda
l_tilde_amd <- results$l_tilde_amd

#----------------------------TASK 5-------------------------------------------#

# Function to apply the t-copula model
t_copula <- function(n_samples, rho, df, x1, x2) {
  # Simulate multivariate t distribution random deviates
  set.seed(123)
  y <- mvrnorm(n = n_samples, mu = c(0, 0), Sigma = matrix(c(1, rho, rho, 1), ncol = 2))
  chi_samples <- rchisq(n_samples, df)
  t_samples <- y / sqrt(chi_samples / df)
  
  # Simulate uniform samples from t distribution deviates
  u <- pt(t_samples, df)
  
  # Transform the uniform samples to empirical marginals
  x_tilde_nvda <- quantile(nvda_data_clean$LogReturn, u[, 1])
  x_tilde_amd <- quantile(amd_data_clean$LogReturn, u[, 2])
  
  # Calculate the potential loss
  l_tilde_nvda <- 1 - exp(x_tilde_nvda)
  l_tilde_amd <- 1 - exp(x_tilde_amd)
  
  # Calculate the loss for the portfolio
  portfolio_loss_tilde <- 0.5 * (l_tilde_nvda + l_tilde_amd)
  
  # Calculate VaR and ES for the portfolio using t-copula
  VaR_ES_t_copula <- sapply(alpha_levels, function(alpha) VaR_ES(portfolio_loss_tilde, alpha))
  
  # Return the results
  list(VaR = VaR_ES_t_copula[1, ], ES = VaR_ES_t_copula[2, ], l_tilde_nvda = l_tilde_nvda, l_tilde_amd = l_tilde_amd)}


# Apply the t-copula model with different degrees of freedom
t3_results <- t_copula(n_samples, rho, 3, nvda_data_clean$LogReturn, amd_data_clean$LogReturn)
t10_results <- t_copula(n_samples, rho, 10, nvda_data_clean$LogReturn, amd_data_clean$LogReturn)
t10000_results <- t_copula(n_samples, rho, 10000, nvda_data_clean$LogReturn, amd_data_clean$LogReturn)

l_tilde_nvda_t3 <- t3_results$l_tilde_nvda
l_tilde_amd_t3 <- t3_results$l_tilde_amd
l_tilde_nvda_t10 <- t10_results$l_tilde_nvda
l_tilde_amd_t10 <- t10_results$l_tilde_amd
l_tilde_nvda_t10000 <- t10000_results$l_tilde_nvda
l_tilde_amd_t10000 <- t10000_results$l_tilde_amd

#----------------------------TASK 6-------------------------------------------#

# Function to simulate comonotonicity copula samples
simulate_comonotonic_copula <- function(n) {
  u_samples <- matrix(runif(n), ncol = 1)
  u_samples <- cbind(u_samples, u_samples)
  return(u_samples)
}

# Function to simulate countermonotonicity copula samples
simulate_countermonotonic_copula <- function(n) {
  u_samples <- matrix(runif(n), ncol = 1)
  u_samples <- cbind(u_samples, 1 - u_samples)
  return(u_samples)
}

# Simulate comonotonicity copula samples
set.seed(5417441)
u_comonotonic <- simulate_comonotonic_copula(n_samples)

# Simulate countermonotonicity copula samples
u_countermonotonic <- simulate_countermonotonic_copula(n_samples)

# Convert to original scale using empirical distribution functions
x_tilde_nvda_comonotonic <- quantile(nvda_data_clean$LogReturn, u_comonotonic[, 1])
x_tilde_amd_comonotonic <- quantile(amd_data_clean$LogReturn, u_comonotonic[, 2])

x_tilde_nvda_countermonotonic <- quantile(nvda_data_clean$LogReturn, u_countermonotonic[, 1])
x_tilde_amd_countermonotonic <- quantile(amd_data_clean$LogReturn, u_countermonotonic[, 2])

# Calculate potential losses
l_tilde_nvda_comonotonic <- 1 - exp(x_tilde_nvda_comonotonic)
l_tilde_amd_comonotonic <- 1 - exp(x_tilde_amd_comonotonic)

l_tilde_nvda_countermonotonic <- 1 - exp(x_tilde_nvda_countermonotonic)
l_tilde_amd_countermonotonic <- 1 - exp(x_tilde_amd_countermonotonic)

# Portfolio loss
portfolio_loss_comonotonic <- 0.5 * (l_tilde_nvda_comonotonic + l_tilde_amd_comonotonic)
portfolio_loss_countermonotonic <- 0.5 * (l_tilde_nvda_countermonotonic + l_tilde_amd_countermonotonic)

# Calculate VaR and ES for comonotonicity and countermonotonicity copulas
VaR_ES_comonotonic <- sapply(alpha_levels, function(alpha) VaR_ES(portfolio_loss_comonotonic, alpha))
VaR_ES_countermonotonic <- sapply(alpha_levels, function(alpha) VaR_ES(portfolio_loss_countermonotonic, alpha))

#----------------------------TASK 7-------------------------------------------#

# Define the function to estimate the conditional probability P(F(X) > alpha | F(Y) > alpha)
estimate_probability <- function(losses_x, losses_y, alpha) {
  exceedances_x <- losses_x > quantile(losses_x, alpha)
  exceedances_y <- losses_y > quantile(losses_y, alpha)
  prob <- mean(exceedances_x & exceedances_y) / mean(exceedances_y)
  return(prob)
}

# Calculate the conditional probabilities for each copula model
prob_empirical <- sapply(alpha_levels, function(alpha) estimate_probability(nvda_data_clean$Loss, amd_data_clean$Loss, alpha))
prob_independence <- sapply(alpha_levels, function(alpha) estimate_probability(l_tilde_nvda_independence, l_tilde_amd_independence, alpha))
prob_gaussian <- sapply(alpha_levels, function(alpha) estimate_probability(l_tilde_nvda, l_tilde_amd, alpha))
prob_t3 <- sapply(alpha_levels, function(alpha) estimate_probability(l_tilde_nvda_t3, l_tilde_amd_t3, alpha))
prob_t10 <- sapply(alpha_levels, function(alpha) estimate_probability(l_tilde_nvda_t10, l_tilde_amd_t10, alpha))
prob_t10000 <- sapply(alpha_levels, function(alpha) estimate_probability(l_tilde_nvda_t10000, l_tilde_amd_t10000, alpha))
prob_comonotonic <- sapply(alpha_levels, function(alpha) estimate_probability(l_tilde_nvda_comonotonic, l_tilde_amd_comonotonic, alpha))
prob_countermonotonic <- sapply(alpha_levels, function(alpha) estimate_probability(l_tilde_nvda_countermonotonic, l_tilde_amd_countermonotonic, alpha))


# Create a data frame to summarise the probabilities
prob_summary <- data.frame(
  Alpha = rep(alpha_levels, 8),
  Model = rep(c("Empirical" , "Independence", "Gaussian", "t-copula df=3", "t-copula df=10", "t-copula df=10000", "Comonotonic", "Countermonotonic"), each = length(alpha_levels)),
  Probability = c(prob_empirical, prob_independence, prob_gaussian, prob_t3, prob_t10, prob_t10000, prob_comonotonic, prob_countermonotonic)
)

# Print the probability summary
print("Summary of probabilities P(F(X) > alpha | F(Y) > alpha) for different copulas:")
print(prob_summary)

# Plot the probabilities
ggplot(prob_summary, aes(x = Alpha, y = Probability, color = Model)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Empirical" = "green", 
                                "Independence" = "cyan", 
                                "Gaussian" = "blue", 
                                "t-copula df=3" = "magenta", 
                                "t-copula df=10" = "purple", 
                                "t-copula df=10000" = "brown",
                                "Comonotonic" = "red", 
                                "Countermonotonic" = "orange")) +
  labs(title = "Conditional Probabilities for Different Copulas",
       x = "Alpha Level",
       y = "Conditional Probability",
       color = "Copula Model") +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))

# Function to calculate mean squared error for VaR and ES
calculate_MSE <- function(estimated_VaR_ES, empirical_VaR_ES) {
  MSE_VaR <- mean((estimated_VaR_ES["VaR", ] - empirical_VaR_ES["VaR", ])^2)
  MSE_ES <- mean((estimated_VaR_ES["ES", ] - empirical_VaR_ES["ES", ])^2)
  return(c(MSE_VaR = MSE_VaR, MSE_ES = MSE_ES))
}

# Extract empirical VaR and ES as named vectors
empirical_nvda_VaR_ES <- matrix(c(empirical_VaR_ES(nvda_data_clean$Loss, alpha_levels), empirical_VaR_ES(amd_data_clean$Loss, alpha_levels)), nrow = 2, byrow = TRUE)
rownames(empirical_nvda_VaR_ES) <- c("VaR", "ES")

# Extract estimated VaR and ES as named vectors for each model
estimated_independence_VaR_ES <- matrix(portfolio_VaR_ES_independence, nrow = 2, byrow = TRUE)
rownames(estimated_independence_VaR_ES) <- c("VaR", "ES")

estimated_gaussian_VaR_ES <- matrix(results$VaR, nrow = 2, byrow = TRUE)
rownames(estimated_gaussian_VaR_ES) <- c("VaR", "ES")

estimated_t3_VaR_ES <- matrix(t3_results$VaR, nrow = 2, byrow = TRUE)
rownames(estimated_t3_VaR_ES) <- c("VaR", "ES")

estimated_t10_VaR_ES <- matrix(t10_results$VaR, nrow = 2, byrow = TRUE)
rownames(estimated_t10_VaR_ES) <- c("VaR", "ES")

estimated_t10000_VaR_ES <- matrix(t10000_results$VaR, nrow = 2, byrow = TRUE)
rownames(estimated_t10000_VaR_ES) <- c("VaR", "ES")

estimated_comonotonic_VaR_ES <- matrix(VaR_ES_comonotonic, nrow = 2, byrow = TRUE)
rownames(estimated_comonotonic_VaR_ES) <- c("VaR", "ES")

estimated_countermonotonic_VaR_ES <- matrix(VaR_ES_countermonotonic, nrow = 2, byrow = TRUE)
rownames(estimated_countermonotonic_VaR_ES) <- c("VaR", "ES")

# Calculate MSE for each model
MSE_independence <- calculate_MSE(estimated_independence_VaR_ES, empirical_nvda_VaR_ES)
MSE_gaussian <- calculate_MSE(estimated_gaussian_VaR_ES, empirical_nvda_VaR_ES)
MSE_t3 <- calculate_MSE(estimated_t3_VaR_ES, empirical_nvda_VaR_ES)
MSE_t10 <- calculate_MSE(estimated_t10_VaR_ES, empirical_nvda_VaR_ES)
MSE_t10000 <- calculate_MSE(estimated_t10000_VaR_ES, empirical_nvda_VaR_ES)
MSE_comonotonic <- calculate_MSE(estimated_comonotonic_VaR_ES, empirical_nvda_VaR_ES)
MSE_countermonotonic <- calculate_MSE(estimated_countermonotonic_VaR_ES, empirical_nvda_VaR_ES)

# Summarise MSE results
MSE_summary <- data.frame(
  Model = c("Independence", "Gaussian", "t-copula df=3", "t-copula df=10", "t-copula df=10000", "Comonotonic", "Countermonotonic"),
  MSE_VaR = c(MSE_independence[1], MSE_gaussian[1], MSE_t3[1], MSE_t10[1], MSE_t10000[1], MSE_comonotonic[1], MSE_countermonotonic[1]),
  MSE_ES = c(MSE_independence[2], MSE_gaussian[2], MSE_t3[2], MSE_t10[2], MSE_t10000[2], MSE_comonotonic[2], MSE_countermonotonic[2])
)

print("Mean Squared Error for VaR and ES for different models:")
print(MSE_summary)

# Plot MSE results
ggplot(MSE_summary, aes(x = Model)) +
  geom_bar(aes(y = MSE_VaR), stat = "identity", position = "dodge", fill = "blue") +
  geom_bar(aes(y = MSE_ES), stat = "identity", position = "dodge", fill = "red") +
  labs(title = "Mean Squared Error for VaR and ES", y = "MSE", x = "Model") +
  theme_bw()
