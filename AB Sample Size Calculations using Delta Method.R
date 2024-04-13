# Load necessary libraries
# library(truncnorm) # For truncated normal distributions
# library(dplyr)     # For data manipulation
# library(pwr)       # For power of test calculations



# Data generation process -----------------------------------------------------

set.seed(12)

# Define parameters
N <- 5000           # Total number of users
lambda <- 5         # Average number of sessions per user
a_mu <- 0.6         # Mean of the truncated normal distribution for conversion rate
b_mu <- 0.65
ate <- b_mu - a_mu
range <- 0.35       # Range for truncation to keep probabilities between 0 and 1
alpha <- 0.05       # Significance level for hypothesis testing
beta <- 0.2         # Power of the test

# Initialize vectors for session and conversion data
sess <- rpois(n = N, lambda = lambda) # sessions per user
conv <- numeric(N)                    # conversions per user
pi_values <- numeric(N)               # user-specific conversion probabilities

# Simulate user behavior for session conversions
for (i in 1:N) {
  pi <- rtruncnorm(n = 1, a = max(0, a_mu - range), b = min(1, a_mu + range), mean = a_mu, sd = range/2)
  pi_values[i] <- pi
  conv[i] <- sum(rbinom(n = sess[i], size = 1, prob = pi))
}




# Estimating 'h' and 'k' ------------------------------------------------------


# Compute 'h' and 'k' for determining sample sizes
h <- 1 / mean(sess)^2 * (var(conv) - 2 * mean(conv) / mean(sess) * cov(sess, conv) + mean(conv)^2 / mean(sess)^2 * var(sess))
k <- 2 * h * (qnorm(1 - alpha / 2) + qnorm(1 - beta))^2 / (ate^2)
k
max_sess <- k * mean(sess)




# Sample size calculation for independent binary data -------------------------

p_pool <- (a_mu + b_mu) / 2   #  pooled proportion

ate <- abs(a_mu - b_mu)       # absolute treatment effect

n <- 2 * p_pool * (1 - p_pool) * (qnorm(1 - alpha / 2) + qnorm(1 - beta))^2 / ate^2
n



# Power_corr function  -------------------------
power_corr <- function(N, lambda, a_mu, b_mu, range, num_simulations, alpha) {
  t_stats <- numeric(num_simulations)
  
  for (sim in 1:num_simulations) {
    a_sess <- rpois(n = N, lambda = lambda)
    b_sess <- rpois(n = N, lambda = lambda)
    a_conv <- numeric(N)
    b_conv <- numeric(N)
    
    for (i in 1:N) {
      a_pi <- rtruncnorm(n = 1, a = max(0, a_mu - range), b = min(1, a_mu + range), mean = a_mu, sd = range/2)
      b_pi <- rtruncnorm(n = 1, a = max(0, b_mu - range), b = min(1, b_mu + range), mean = b_mu, sd = range/2)
      a_conv[i] <- sum(rbinom(n = a_sess[i], size = 1, prob = a_pi))
      b_conv[i] <- sum(rbinom(n = b_sess[i], size = 1, prob = b_pi))
    }
    
    
    X_mean <- (sum(a_conv)/N) / (sum(a_sess)/N)
    Y_mean <- (sum(b_conv)/N) / (sum(b_sess)/N)
    a_delta <- 1 / (mean(a_sess)^2 * N) * (var(a_conv) - 2 * mean(a_conv)/mean(a_sess) * cov(a_sess, a_conv) + mean(a_conv)^2/mean(a_sess)^2 * var(a_sess))
    b_delta <- 1 / (mean(b_sess)^2 * N) * (var(b_conv) - 2 * mean(b_conv)/mean(b_sess) * cov(b_sess, b_conv) + mean(b_conv)^2/mean(b_sess)^2 * var(b_sess))
    # Calculate t_stat for this simulation
    t_stat <- (Y_mean - X_mean) / sqrt(a_delta + b_delta)
    t_stats[sim] <- t_stat
  }
  
  power <- sum(abs(t_stats) > qnorm(1 - alpha / 2))
  power_corr <-  power / num_simulations
  
  return(power_corr)
}




# Power_independent i function  -------------------------

power_corr_a <- function(N, lambda, a_mu, b_mu, range, num_simulations, alpha, target_sessions) {
  t_stats <- numeric(num_simulations)
  
  for (sim in 1:num_simulations) {
    # Initialize variables
    a_sess <- rpois(n = N, lambda = lambda)
    b_sess <- rpois(n = N, lambda = lambda)
    a_conv <- numeric(N)
    b_conv <- numeric(N)
    
    # Initialize cumulative session count and user index
    cumulative_sessions <- 0
    last_included_user <- N
    
    # Simulate for each user
    for (i in 1:N) {
      # Check if target sessions reached and break loop if so
      if (cumulative_sessions >= target_sessions) {
        last_included_user <- i - 1
        break
      }
      
      # Update cumulative sessions
      cumulative_sessions <- cumulative_sessions + a_sess[i] + b_sess[i]
      # Simulate conversions
      a_pi <- rtruncnorm(n = 1, a = max(0, a_mu - range), b = min(1, a_mu + range), mean = a_mu, sd = range/2)
      b_pi <- rtruncnorm(n = 1, a = max(0, b_mu - range), b = min(1, b_mu + range), mean = b_mu, sd = range/2)
      a_conv[i] <- sum(rbinom(n = a_sess[i], size = 1, prob = a_pi))
      b_conv[i] <- sum(rbinom(n = b_sess[i], size = 1, prob = b_pi))
    }
    
    # Limit data to users up to last_included_user
    a_sess <- a_sess[1:last_included_user]
    b_sess <- b_sess[1:last_included_user]
    a_conv <- a_conv[1:last_included_user]
    b_conv <- b_conv[1:last_included_user]
    
    X_mean <- (sum(a_conv)/N) / (sum(a_sess)/N)
    Y_mean <- (sum(b_conv)/N) / (sum(b_sess)/N)
    a_delta <- 1 / (mean(a_sess)^2 * N) * (var(a_conv) - 2 * mean(a_conv)/mean(a_sess) * cov(a_sess, a_conv) + mean(a_conv)^2/mean(a_sess)^2 * var(a_sess))
    b_delta <- 1 / (mean(b_sess)^2 * N) * (var(b_conv) - 2 * mean(b_conv)/mean(b_sess) * cov(b_sess, b_conv) + mean(b_conv)^2/mean(b_sess)^2 * var(b_sess))
    
    # Calculate t_stat for this simulation
    t_stat <- (Y_mean - X_mean) / sqrt(a_delta + b_delta)
    t_stats[sim] <- t_stat
  }
  
  power <- sum(abs(t_stats) > qnorm(1 - alpha / 2))
  power_corr_a <-  power / num_simulations
  
  return(power_corr_a)
}



# Power_independent ii function  -------------------------

power_corr_b <- function(N, lambda, a_mu, b_mu, range, num_simulations, alpha, max_sessions) {
  
  t_stats <- numeric(num_simulations)
  # Maximum number of sessions per group, from initial k calculations
  max_sess <- max_sessions
  
  for (sim in 1:num_simulations) {
    a_sess <- rpois(n = N, lambda = lambda)
    b_sess <- rpois(n = N, lambda = lambda)
    
    # Proportionally reduce sessions for each group if they exceed the maximum
    if (sum(a_sess) > max_sess) {
      reduction_factor_a <- max_sess/ sum(a_sess)
      a_sess <- round(a_sess * reduction_factor_a)
    }
    if (sum(b_sess) > max_sess) {
      reduction_factor_b <- max_sess / sum(b_sess)
      b_sess <- round(b_sess * reduction_factor_b)
    }
    
    a_conv <- numeric(N)
    b_conv <- numeric(N)
    
    for (i in 1:N) {
      a_pi <- rtruncnorm(n = 1, a = max(0, a_mu - range), b = min(1, a_mu + range), mean = a_mu, sd = range/2)
      b_pi <- rtruncnorm(n = 1, a = max(0, b_mu - range), b = min(1, b_mu + range), mean = b_mu, sd = range/2)
      
      a_conv[i] <- sum(rbinom(n = a_sess[i], size = 1, prob = a_pi))
      b_conv[i] <- sum(rbinom(n = b_sess[i], size = 1, prob = b_pi))
    }
    
    X_mean <- (sum(a_conv)/N) / (sum(a_sess)/N)
    Y_mean <- (sum(b_conv)/N) / (sum(b_sess)/N)
    
    a_delta <- 1 / (mean(a_sess)^2 * N) * (var(a_conv) - 2 * mean(a_conv)/mean(a_sess) * cov(a_sess, a_conv) + mean(a_conv)^2/mean(a_sess)^2 * var(a_sess))
    b_delta <- 1 / (mean(b_sess)^2 * N) * (var(b_conv) - 2 * mean(b_conv)/mean(b_sess) * cov(b_sess, b_conv) + mean(b_conv)^2/mean(b_sess)^2 * var(b_sess))
    
    # Calculate t_stat for this simulation
    t_stat <- (Y_mean - X_mean) / sqrt(a_delta + b_delta)
    t_stats[sim] <- t_stat
  }
  
  power_corr_b <- sum(abs(t_stats) > qnorm(1 - alpha / 2)) / num_simulations
  
  return(power_corr_b)
}




# Parameters  -------------------------
# Example usage
N <- k
lambda <- 5
a_mu <- 0.6
b_mu <- 0.65
range <-  0.35
alpha <-  0.05
beta <- 0.2
num_simulations <- 100
target_sessions <- n*2
max_sessions <- max_sess

power_corr <- power_corr(N, lambda, a_mu, b_mu, range, num_simulations, alpha)

power_corr_a <- power_corr_a(N, lambda, a_mu, b_mu, range, num_simulations, alpha, target_sessions)

power_corr_b <- power_corr_b(N, lambda, a_mu, b_mu, range, num_simulations, alpha, max_sessions)