PARAMETER ESTIMATION FOR THE FRACTIONAL HAWKES PROCESS

This repository contains the code used in the paper "Parameter Estimation for Fractional Hawkes Process". The paper presents methods for simulating, analyzing, and estimating parameters for fractional Hawkes processes using various statistical and computational techniques.

In this study, we focus on parameter estimation for the Fractional Hawkes Process (FHP) using Approximate Bayesian Computation (ABC) and Maximum Likelihood Estimation (MLE) techniques . The codebase provided here includes the following key components:

1. Simulation of Fractional Hawkes Process data: Generating time realizations for the Hawkes process with Mittag-Leffler distributed waiting times using the package 'hawkesbow'.
2. Distance Functions: Various summary statistics are used to calculate distances between simulated and observed data points for the ABC method, including:
        - Log-number distance (LN)
        - Wasserstein distance (WS)
        - Mean difference greater than q90 (MD90)
        - Mean difference less than q50 (MD50)
        - Digamma distance (DIGAD)
3. Parameter Estimation Methods: Implementing ABC and MLE for estimating key parameters (lambda, alpha, beta).
4. Results: Output plots and results from the simulations and parameter estimation are saved.

Each file has R code to produce a single figure with posterior distribution estimates and performance metrics (e.g., Kullback-Leibler divergence) that have been reported in the paper. It is recommended to run the code using a high performance computer to reduce the computation time.

Inside the R script, you can modify the value of the key parameters like lambda, alpha, beta, and Timev to test different setups or run more simulations (Max_sim).

After running the script, the estimated posterior distributions for the parameters (lambda, alpha, and beta) will be generated using both ABC and MLE.
Kullback-Leibler Divergence: The comparison between the posterior distribution estimated via ABC and MLE will be summarized through KLD values. The execution time for each method will be generated to compare the ABC and MLE computation efficiency.
