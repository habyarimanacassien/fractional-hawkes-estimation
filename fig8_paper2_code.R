rm(list=ls())  # Clear all objects from the current R environment

# Useful libraries
library(rstudioapi)  # Library to interact with the RStudio API
cur_dir = dirname(getSourceEditorContext()$path)  # Get the current script directory
setwd(cur_dir)  # Set working directory to the current script location
library(MittagLeffleR)  # Mittag-Leffler random variable generation package
library(pracma)  # Useful functions for numerical computation
library("rlist")  # Work with lists in R
require(parallel)  # Parallel processing library
library(ggplot2)  # Data visualization package
library(LaplacesDemon)  # Bayesian inference using MCMC

################################################################################
#                SUMMARY STATISTICS & DISTANCE FUNCTIONS TO USE                #
################################################################################

# 1. Function to compute_diff_log_nbr_distance (LN)
compute_diff_log_nbr_distance <- function(distr1, distr2){
  # It computes and returns the log number difference distance between two distributions
  n2 <- length(distr2)  # Get the length of the second distribution
  n1 <- length(distr1)  # Get the length of the first distribution
  distance = abs(log(n2)-log(n1))  # Compute the absolute difference of log lengths
  return (distance)  # Return the computed distance
  # End function
}

# 2. Function to compute_wasserstein_distance (WS)
compute_wasserstein_distance <- function(distr1, distr2, Time){
  # It computes and returns the Wasserstein distance between two distributions
  n = min(length(distr1),length(distr2))  # Get the smaller length of the two distributions
  m = max(length(distr1),length(distr2))  # Get the larger length of the two distributions
  
  # Handle the case when distr1 is longer than distr2
  if (length(distr1) > length(distr2)) {
    if (n == 0) {
      distance = (m - n) * Time - sum(na.omit(distr1[n+1:m]))  # Calculate distance for excess points
    } else {
      distance = sum(abs(distr1[1:n] - distr2[1:n])) + 
        (m - n) * Time - sum(na.omit(distr1[n+1:m]))  # Calculate distance and adjust for excess points
    }
  # Handle the case when distr2 is longer than distr1
  } else if (length(distr1) < length(distr2)) {
    if (n == 0) {
      distance = (m - n) * Time - sum(na.omit(distr2[n+1:m]))  # Calculate distance for excess points
    } else {
      distance = sum(abs(distr1[1:n] - distr2[1:n])) + 
        (m - n) * Time - sum(na.omit(distr2[n+1:m]))  # Calculate distance and adjust for excess points
    }
  } else {
    distance = sum(abs(distr1 - distr2))  # Calculate distance if both distributions have the same length
  }
  return (distance)  # Return the computed distance
  # End function
}

# 3. Function to compute_meandiff_grt_q90_distance (MD90)
compute_meandiff_grt_q90_distance <- function(distr1, distr2){
  # It computes and returns the difference between the means of event time differences greater than q90 for two distributions
  q1 <- unname(quantile(diff(distr1), probs = .90))  # Compute 90th percentile for the first distribution
  s1 <- mean(diff(distr1)[diff(distr1) > q1])  # Compute mean of values greater than q90 for the first distribution
  if (length(distr2) > 2){
    q2 <- unname(quantile(diff(distr2), probs = .90))  # Compute 90th percentile for the second distribution
    s2 <- mean(diff(distr2)[diff(distr2) > q2])  # Compute mean of values greater than q90 for the second distribution
    distance <- abs(s2 - s1)  # Compute the absolute difference between the two means
  } else {
    distance <- Inf  # If distr2 has fewer than 3 elements, return infinite distance
  }
  return (distance)  # Return the computed distance
  # End function
}

# 4. Function to compute_meandiff_less_q50_distance (MD50)
compute_meandiff_less_q50_distance <- function(distr1, distr2){
  # It computes and returns the difference between the means of event time differences less than q50 for two distributions
  q1 <- unname(quantile(diff(distr1), probs = .50))  # Compute 50th percentile for the first distribution
  s1 <- mean(diff(distr1)[diff(distr1) < q1])  # Compute mean of values less than q50 for the first distribution
  
  if (length(distr2) > 2){
    q2 <- unname(quantile(diff(distr2), probs = .50))  # Compute 50th percentile for the second distribution
    s2 <- mean(diff(distr2)[diff(distr2) < q2])  # Compute mean of values less than q50 for the second distribution
    distance <- abs(s2 - s1)  # Compute the absolute difference between the two means
  } else {
    distance <- Inf  # If distr2 has fewer than 3 elements, return infinite distance
  }
  return (distance)  # Return the computed distance
  # End function
}

# 5. Function to compute digamma distance (DIGAD)
digamma_distance <- function(distr1, distr2){
  # It computes and returns the difference in digamma sums for two distributions
  if (length(distr2) > 1){
    s1 = sum(digamma(diff(distr1)), na.rm = TRUE)  # Compute the digamma sum for the first distribution
    s2 = sum(digamma(diff(distr2)), na.rm = TRUE)  # Compute the digamma sum for the second distribution
    if (is.finite(s1) & is.finite(s2)){
      distance = abs(s1 - s2)  # Compute the absolute difference between the two sums
    } else {
      distance = Inf  # If the digamma sum is not finite, return infinite distance
    }
  } else {
    distance = Inf  # If distr2 has fewer than 2 elements, return infinite distance
  }
  return (distance)  # Return the computed distance
  # End function
}

################################################################################
# function to convert time (in seconds) into days, hours, minutes, and seconds #
################################################################################

time_converter <- function(start_time, end_time){
  # Convert the time difference into seconds
  time <- as.numeric((end_time - start_time), units = "secs")
  # Convert seconds into days, hours, minutes, and seconds
  days <- floor(time / (24 * 3600))  # Calculate number of days
  hours <- floor((time %% (24 * 3600)) / 3600)  # Calculate number of hours
  minutes <- floor((time %% 3600) / 60)  # Calculate number of minutes
  seconds <- round(time %% 60, 2)  # Calculate number of seconds
  
  # Create a list to store the time breakdown
  time_converted <- list(
    # Store the time breakdown in a string format
    time = paste0(days, "d", hours, "h", minutes, "m", seconds, "s")
  )
  return(time_converted)  # Return the converted time
}

################################################################################
#             FRACTIONAL HAWKES PROCESS (FHP) PARAMETER ESTIMATION             #
################################################################################

# Function to simulate observation points for FHP
simulate_points <- function(lambda, alpha, beta, Time){
  # It simulates and returns fractional Hawkes process number of observations and time realization points
  # Hawkes process time realization
  library(MittagLeffleR)
  library(LaplacesDemon)
  library(hawkesbow)  # Load the required packages for simulation
  
  x <- hawkes(Time, fun = lambda, repr = alpha, family = function(n) {rml(n, beta)})  # Simulate Hawkes process
  SimPts <- c(0, x$p)  # Combine time points with 0 as the initial point
  
  return (SimPts)  # Return the simulated points
  # End function
}

################################################################################
#                      FUNCTION TO COMPUTE THE LIKELIHOOD                      #
################################################################################

# Likelihood function for FRACTIONAL HAWKES PROCESS
compute_logLikelihood <- function(distr, lambda, alpha, beta, Time){
  library(MittagLeffleR)  # Load the required library
  k = length(distr)  # Get the number of observation points
  ll1 = 0  # Initialize the log-likelihood
  
  # Loop through the points to compute the log-likelihood
  for (i in 1:k){
    si = 0
    if (i > 1){
      # Compute the sum of past events contributing to the likelihood
      for(j in 1:(i-1)){
        si = si + ((distr[i] - distr[j])^(beta - 1)) * mlf(-(distr[i] - distr[j])^beta, beta, beta, 1)
      }
    }
    # Update the log-likelihood for the current point
    ll1 = ll1 + log(lambda + alpha * si) - alpha * (1 - mlf(-(distr[k] - distr[i])^beta, beta, 1, 1))
  }
  
  # Adjust the log-likelihood by subtracting the base rate contribution
  ll <- ll1 - lambda * Time
  # Return the computed log-likelihood
  return (ll)
  # End function
}

################################################################################
#                              1. FHP2PP_LAMBDA                                #
################################################################################

# Process identification
proc <- "FHP2PP"  # Process name
par <- "Lambda"  # Parameter to estimate
dataset <- "SimulatedData"  # Data being used

# Cores, Number of iterations, and max of simulated data sets setting
number_of_cores <- 24  # Number of cores or clusters for parallel computation
Niter <- 100000  # Number of iterations for the simulation
Max_sim <- 5  # Max number of different simulated datasets for each process

# Inputs
Timev <- c(100, 1000)  # Time intervals for the process
truelambda <- 0.5  # True value of the parameter to estimate
alpha <- 0.0001  # Small assumed alpha value
beta <- 0.5  # Assumed beta value

# List to store results of execution time for each simulation
result_list <- list()

# Loop for each simulation (Max_sim times)
for (sim in 1:Max_sim){
  start_sim <- Sys.time()  # Start timing the simulation
  
  # Save figures in a PDF file
  pdf(file = paste0("Fig", proc, "", par, "", dataset, "_1", "_N", Niter, "_i", sim, ".pdf"), 
      width = 12.5, height = 10.5)  # Create a PDF file for output figures
  par(mar = c(4, 4, 3, 1))  # Adjust margins for the plot
  par(mfrow = c(4,5), oma = c(3,3,3,3))  # Set up a multi-plot layout
  
  # Loop for each time interval
  for (iT in 1:length(Timev)){
    Time <- Timev[iT]  # Get the current time interval
    
    # Simulate observation points with true parameter value
    SimPoints <- simulate_points(lambda = truelambda, alpha = alpha, beta = beta, Time = Time)
    
    ################################################################################
    #                           ABC PARAMETER ESTIMATION                           #
    ################################################################################
    
    # Start timing the ABC estimation
    ABC_start_time <- Sys.time()
    
    # Generate a uniformly distributed prior for the parameter sample
    ParameterSample_values <- runif(Niter, min = 0, max = 1)
    
    # Simulate data with uniform prior for lambda parameter
    cl <- makeCluster(number_of_cores)  # Open the cluster for parallel computation
    SimPointsth_list <- parSapply(cl , ParameterSample_values, simulate_points, alpha = alpha, beta = beta, Time = Time)
    stopCluster(cl)  # Stop the cluster after the computation
    
    # End timing the ABC estimation
    ABC_end_time <- Sys.time()
    
    # Compute the ABC execution time
    time_ABC <- time_converter(ABC_start_time, ABC_end_time)
    
    ################################################################################
    #                      COMPUTE DISTANCE & SUMMARY STATISTICS                   #
    ################################################################################
    
    cl <- makeCluster(number_of_cores)  # Open the cluster for parallel computation
    dist1 <- "LN"  # Log number distance (LN)
    out_distance1 <- parLapply(cl, SimPointsth_list, compute_diff_log_nbr_distance, distr1 = SimPoints)
    stopCluster(cl)  # Stop the cluster after the computation
    
    cl <- makeCluster(number_of_cores)
    dist2 <- "WS"  # Wasserstein distance (WS)
    out_distance2 <- parLapply(cl, SimPointsth_list, compute_wasserstein_distance, distr1 = SimPoints, Time = Time)
    stopCluster(cl)
    
    cl <- makeCluster(number_of_cores)
    dist3 <- "MDG90"  # Mean difference greater than q90 (MD90)
    out_distance3 <- parLapply(cl, SimPointsth_list, compute_meandiff_grt_q90_distance, distr1 = SimPoints)
    stopCluster(cl)
    
    cl <- makeCluster(number_of_cores)
    dist4 <- "MDL50"  # Mean difference less than q50 (MD50)
    out_distance4 <- parLapply(cl, SimPointsth_list, compute_meandiff_less_q50_distance, distr1 = SimPoints)
    stopCluster(cl)
    
    cl <- makeCluster(number_of_cores)
    dist5 <- "DIGAD"  # Digamma distance (DIGAD)
    out_distance5 <- parLapply(cl, SimPointsth_list, digamma_distance, distr1 = SimPoints)
    stopCluster(cl)
    
    ################################################################################
    #                          MLE PARAMETER ESTIMATION                            #
    ################################################################################
    
    # Start timing the MLE estimation
    MLE_start_time <- Sys.time()
    
    # Set the bin size and interval for likelihood computation
    dx <- 1/200  # Bin size
    LL <- 1/200  # Lower limit of the interval
    UL <- 199/200  # Upper limit of the interval
    breaks_interval <- seq(LL - dx, UL + dx, by = dx)  # Create the breaks interval
    mids_interval <- seq(LL - dx/2, UL + dx/2, by = dx)  # Create the mids interval
    
    # Compute the log-likelihood for each mid value
    cl <- makeCluster(number_of_cores)  # Open the cluster for parallel computation
    logLikelihood <- parSapply(cl, mids_interval, compute_logLikelihood, distr = SimPoints, beta = beta, alpha = alpha, Time = Time)
    stopCluster(cl)  # Stop the cluster after the computation
    
    # Compute the likelihood by exponentiating the log-likelihood
    likelihood <- exp(logLikelihood - mean(logLikelihood))  # Remove impossibilities
    
    # End timing the MLE estimation
    MLE_end_time <- Sys.time()
    
    # Compute the MLE execution time
    time_MLE <- time_converter(MLE_start_time, MLE_end_time)
    
    # Store the results for ABC and MLE execution times
    result_list[[paste0("Time_", Time, "Sim_", sim)]] <- list(
      Time = unlist(c("ABC" = time_ABC, "MLE" = time_MLE))
    )
    
    ############################################################################
    # Selection rate
    rate_v <- c(0.05, 0.1)  # Selection rates
    
    for (j in 1:length(rate_v)){
      
      # DATAFRAMES FOR DISTANCES & SUMMARY STATISTICS, AND THEIR NAMES 
      data_frame1 <- list(out_distance1, out_distance2, out_distance3, out_distance4, out_distance5)
      data_frame2 <- list(dist1, dist2, dist3, dist4, dist5)
      
      for (i in 1:length(data_frame2)){ # Iterate through distances
        out_distance <- as.numeric(data_frame1[[i]])
        dist <- as.character(data_frame2[[i]])
        
        rate <- rate_v[j]
        d <- quantile(out_distance, probs = rate, na.rm = TRUE)  # Compute the threshold distance
        
        # Select parameter values corresponding to distances less than threshold distance d
        abc_posterior <- ParameterSample_values[which(out_distance < d)]
        
        ########################################################################
        ## Kullback–Leibler divergence
        ########################################################################
        # With classical histogram
        h <- hist(abc_posterior, breaks = breaks_interval, plot = FALSE)  # Create a histogram for the posterior
        density_ABC <- h$density
        
        posterior <- likelihood  # Likelihood for uniform prior in [0, 1]
        exact_posterior <- posterior / sum(posterior) / (dx)  # Normalization
        density_LIK <- exact_posterior
        
        # Penalize zero values in the density
        eps <- 1e-10
        for (v in 1:length(density_ABC)){
          if (density_ABC[v] == 0){
            density_ABC[v] = eps  # Replace zero with small epsilon value
          }
          if (density_LIK[v] == 0){
            density_LIK[v] = eps  # Replace zero with small epsilon value
          }
        }
        
        # Compute Kullback-Leibler divergence using histogram densities
        KLD_value <- KLD(density_LIK, density_ABC)  # Compute KLD(Likelihood, ABC)
        
        # Prepare text for plot
        dist_used = paste0("Dist : ", dist)
        threshold_dist = bquote("Threshold "*epsilon == .(round(d, 3)))
        sel_rate = paste0("Sel Rate = ", round(rate * 100, 1), "%")
        dkl_LIK_ABC <- paste0("KLD(LIK, ABC) = ", round(KLD_value[[4]], 3))
        dkl_ABC_LIK <- paste0("KLD(ABC, LIK) = ", round(KLD_value[[5]], 3))
        average <- paste0("Post Mean = ", round(mean(abc_posterior), 3))
        stand_dev <- paste0("Post STD = ", round(std(abc_posterior), 3))
        
        min_h1 <- min(abc_posterior)
        max_h1 <- max(abc_posterior)
        dx_h1 <- 0.001
        h1 <- hist(abc_posterior, breaks = seq(min_h1 - dx_h1, max_h1 + dx_h1, dx_h1), plot = FALSE) 
        ymax <- max(max(h1$density), max(exact_posterior))
        
        # Create plot
        fig <- function(){
          h2 <- hist(abc_posterior, probability = TRUE, ylim = c(0, ceiling(1.35 * ymax)), 
                     col = "darkolivegreen", border = "#333333", xlab = "", ylab = "", main = "", xlim = c(0, 1))
          lines(h$mids, exact_posterior, col = 'red', type = "l", lwd = 2)  # Add posterior line
          segments(x0 = truelambda, y0 = 0, x1 = truelambda, y1 = ymax, col = 'blue', lwd = 2)  # True lambda
          text(-0.015, ymax * 1.25, dkl_ABC_LIK , pos = 4, col = "#330000")  # Add KLD text
          text(-0.015, ymax * 1.15, average , pos = 4, col = "#330000")  # Add average text
          text(-0.015, ymax * 1.05, stand_dev , pos = 4, col = "#330000")  # Add standard deviation text
        }
        
        fig()  # Plot the figure
        mtext(paste0(dist), side = 3, line = 0, outer = FALSE, las = 1, col = "darkblue", lwd = 0.5)  # Add distance name
      }
      mtext(expression(Lambda[0]), side = 1, line = 0, outer = TRUE, las = 0)  # Label for lambda
      mtext("Density", side = 2, line = 0, outer = TRUE, las = 0)  # Label for density
      mtext(paste0("\n T = ", Time, "\n Sel Rate = ", round(rate * 100, 1), "%"), side = 3, line = 2, outer = FALSE, las = 1, col = "darkred")  # Add selection rate info
    }
    
  }
  dev.off()  # Close the PDF
  
  ############################
  end_sim <- Sys.time()  # End the simulation timing
  
  # Compute the simulation execution time
  time_sim <- time_converter(start_sim, end_sim)
  
  # Print the execution summary for each simulation
  cat(paste0("Iteration ", sim, " of ", proc, "_", par, " is completed.", " (", format(start_sim, "%Y-%m-%d %H:%M:%S"), 
             " - ", format(end_sim, "%Y-%m-%d %H:%M:%S"), "). Duration: ", time_sim, "\n"))
}

# Save the results for all simulations
saveRDS(result_list, file = paste0("Execution_Time_", proc, "", par, "", dataset, "_N", Niter, ".rds"))  # Save the result list to an RDS file
