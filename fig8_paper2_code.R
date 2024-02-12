
#set directory as folder containing current file
#install.packages("rstudioapi")
library(rstudioapi)
cur_dir = dirname(getSourceEditorContext()$path)
setwd(cur_dir)

number_of_cores <- 2  # number of cores or clusters
Niter <- 100000  # number of iterations
Max_iter <- 1  # max number of different simulated data sets for each process

# Useful libraries
#install.packages("rlist")
library(MittagLeffleR)
library(pracma)
library("rlist")
# load parallel package
require(parallel)
library(ggplot2)
library(LaplacesDemon)



#############################################################################
#############################################################################

### SUMMARY STATISTICS & DISTANCE FUNCTIONS TO USE

# 1. Function to compute_diff_log_nbr_distance (LN)
compute_diff_log_nbr_distance <- function(distr1, distr2){
  # It compute and return the wasserstein distance between two distributions
  n2 <- length(distr2)
  n1 <- length(distr1)
  distance = abs(log(n2)-log(n1))
  return (distance)
  #End function
}


# 2. Function to compute_wasserstein_distance (WS)
compute_wasserstein_distance <- function(distr1, distr2, Time){
  # It compute and return the wasserstein distance between two distributions
  
  n = min(length(distr1),length(distr2))
  m = max(length(distr1),length(distr2))
  
  if (length(distr1) > length(distr2)) {
    if (n==0){
      distance = (m-n)*Time - sum(na.omit(distr1[n+1:m]))
    } else {
      distance = sum(abs(distr1[1:n]-distr2[1:n])) + (m-n)*Time - sum(na.omit(distr1[n+1:m]))
    }
  } else if (length(distr1) < length(distr2)) {
    if (n==0){
      distance = (m-n)*Time - sum(na.omit(distr2[n+1:m]))
    } else {
      distance = sum(abs(distr1[1:n]-distr2[1:n])) + (m-n)*Time - sum(na.omit(distr2[n+1:m]))
    }
  } else {
    distance = sum(abs(distr1-distr2))
  }
  return (distance)
  #End function
}


# 3. Function to compute_meandiff_grt_q90_distance (MD90)
compute_meandiff_grt_q90_distance <- function(distr1, distr2){
  # It compute and return the diff between mean of event time diff
  # that great than their q90 for two distributions 
  q1 <- unname(quantile(diff(distr1), probs = .90))
  s1 <- mean(diff(distr1)[diff(distr1) > q1])
  if (length(distr2) > 2){
    q2 <- unname(quantile(diff(distr2), probs = .90))
    s2 <- mean(diff(distr2)[diff(distr2) > q2])
    distance <- abs(s2 - s1)
  } else {distance <- Inf}
  return (distance)
  #End function
}


# 4. Function to compute_meandiff_less_q50_distance (MD50)
compute_meandiff_less_q50_distance <- function(distr1, distr2){
  # It compute and return the diff between mean of event time diff
  # that great than their q90 for two distributions 
  q1 <- unname(quantile(diff(distr1), probs = .50))
  s1 <- mean(diff(distr1)[diff(distr1) < q1])
  
  if (length(distr2) > 2){
    q2 <- unname(quantile(diff(distr2), probs = .50))
    s2 <- mean(diff(distr2)[diff(distr2) < q2])
    distance <- abs(s2 - s1)
  } else {distance <- Inf}
  return (distance)
  #End function
}


# 5. Function to digamma_distance (DIGAD)
digamma_distance <- function(distr1, distr2){
  if (length(distr2) > 1){
    s1 = sum(digamma(diff(distr1)), na.rm = TRUE)
    s2 = sum(digamma(diff(distr2)), na.rm = TRUE)
    if(is.finite(s1) & is.finite(s2)){
      distance = abs(s1-s2)
    }else{
      distance = Inf
    }
  } else {
    distance = Inf
  }
  return (distance)
  #End function
}


#############################################################################
#############################################################################

#FRACTIONAL HAWKES PROCESS (FHP) PARAMETER ESTIMATION

#Function to simulate observation points for FHP
simulate_points <- function(lambda, alpha, beta, Time){
  #It simulates and returns fractional Hawkes process number of observations
  #and time realization points
  # Hawkes process time realization
  library(MittagLeffleR)
  library(LaplacesDemon)
  library(hawkesbow)
  
  x <- hawkes(Time,fun=lambda,repr=alpha,family=function(n) {rml(n,beta)})
  SimPts <- c(0, x$p)
  
  return (SimPts)
  #end function
}


#############################################################################
##### 1.  FHP2PP_LAMBDA

proc <- "FHP2PP"
par <- "Lambda"
dataset <- "SimulatedData"


### PARAMETER ESTIMATION

#Inputs
Timev <- c(100, 1000) # Final time interval
truelambda <- 0.5  # true value of the parameter to estimate  
alpha <- 0.0001      ## ASSUMED SMALL, CLOSER TO ZERO (0)
beta <- 0.5 

for (iT in 1:length(Timev)){  #iT in 1:2
  Time <- Timev[iT]
  
  for (iter in 1:Max_iter){ # iter in 1:5
    # Observed data (SimPoints) with a true parameter value
    SimPoints <- simulate_points(lambda=truelambda, alpha=alpha, beta=beta, Time=Time) 
    
    # We generate a uniformly distributed prior
    ParameterSample_values <- runif(Niter,min=0,max=1)    # We generate a uniformly distributed prior
    
    #Simulated data (SimPointsth) with uniform prior for LAMBDA parameter
    cl <- makeCluster(number_of_cores)  # 'number_of_cores' is the number of clusters (cores)
    SimPointsth_list <- parSapply(cl , ParameterSample_values, simulate_points, alpha=alpha, beta=beta, Time=Time)
    stopCluster(cl)
    
    
    ### Run ABC algorithm
    ### 1. compute_diff_log_nbr_distance (LN)
    cl <- makeCluster(number_of_cores)  # 'number_of_cores' is the number of clusters (cores)
    out_distance1 <- parSapply(cl , SimPointsth_list, compute_diff_log_nbr_distance, distr1 = SimPoints)
    stopCluster(cl)  
    ### 2. compute_wasserstein_distance (WS)
    cl <- makeCluster(number_of_cores)  
    out_distance2 <- parSapply(cl , SimPointsth_list, compute_wasserstein_distance, distr1 = SimPoints, Time = Time)
    stopCluster(cl)  
    ### 3. compute_meandiff_grt_q90_distance (MDG90)
    cl <- makeCluster(number_of_cores)
    out_distance3 <- parSapply(cl , SimPointsth_list, compute_meandiff_grt_q90_distance, distr1 = SimPoints)
    stopCluster(cl)
    ### 4. compute_meandiff_less_q50_distance (MDL50)
    cl <- makeCluster(number_of_cores)
    out_distance4 <- parSapply(cl , SimPointsth_list, compute_meandiff_less_q50_distance, distr1 = SimPoints)
    stopCluster(cl)
    ## 5. digamma_distance (DIGAD)
    cl <- makeCluster(number_of_cores)
    out_distance5 <- parSapply(cl , SimPointsth_list, digamma_distance, distr1 = SimPoints)
    stopCluster(cl)
    
    # Create dataframes and Save into files
    df1 <- data.frame(out_distance1, out_distance2, out_distance3, out_distance4, out_distance5, ParameterSample_values)
    save(df1, file=paste0("data",proc,"_",par,"_",dataset,"_N",Niter,"_T",Time,"_D1_i",iter,".rdata"))
    df2 <- data.frame(SimPoints)
    save(df2, file=paste0("data",proc,"_",par,"_",dataset,"_N",Niter,"_T",Time,"_D2_i",iter,".rdata"))
    #    save(SimPointsth_list, file=paste0("data",proc,"_",par,"_",dataset,"_SimPointsth_list_N",Niter,"_T",Time,"_i",iter,".rdata"))    
  }
}


############################################################################
### PLOT ###
############################################################################

rate_v <- c(0.05, 0.1)


#Function to compute the likelihood function  for FRACTIONAL HAWKES PROCESS
compute_logLikelihood <- function(distr, lambda, alpha, beta, Time){
  library(MittagLeffleR)
  k = length(distr)
  ll1 = 0
  for (i in 1:k){
    si = 0
    if (i > 1){
      for(j in 1:(i-1)){
        si = si + ((distr[i]-distr[j])^(beta-1))*mlf(-(distr[i]-distr[j])^beta,beta,beta,1)
      }}
    ll1 = ll1 + log(lambda + alpha*si) - alpha*(1 - mlf(-(distr[k]-distr[i])^beta,beta,1,1))
  }
  ll <- ll1 - lambda*Time
  #  likelihood <- exp(ll)
  return (ll)
  #end function
}


for (iter in 1:Max_iter){
  
  # Save FIGURES in a pdf or png file
  pdf(file = paste0("Fig",proc,"_",par,"_",dataset,"_1","_N",Niter,"_i",iter,".pdf"), width=12.5, height=10.5)
  par(mar = c(4, 4, 3, 1))        # Reduce space around plots
  par(mfrow=c(4,5),oma=c(3,3,3,3))
  
  for (iT in 1:length(Timev)){
    Time <- Timev[iT]
    
    # LOAD THE DATA SAVED
    load(file=paste0("data",proc,"_",par,"_",dataset,"_N",Niter,"_T",Time,"_D1_i",iter,".rdata"))
    load(file=paste0("data",proc,"_",par,"_",dataset,"_N",Niter,"_T",Time,"_D2_i",iter,".rdata"))
    # Extract data
    out_distance1 <- df1$out_distance1
    out_distance2 <- df1$out_distance2
    out_distance3 <- df1$out_distance3
    out_distance4 <- df1$out_distance4
    out_distance5 <- df1$out_distance5
    ParameterSample_values <- df1$ParameterSample_values
    SimPoints <- df2$SimPoints
    
    
    ### DISTANCE NAMES 
    dist1 <- "LN"
    dist2 <- "WS"
    dist3 <- "MDG90"
    dist4 <- "MDL50"
    dist5 <- "DIGAD"
    
    data_frame1 <- data.frame(out_distance1, out_distance2, out_distance3, out_distance4, out_distance5)
    data_frame2 <- data.frame(dist1, dist2, dist3, dist4, dist5)
    
    ## parameter values interval: [0, 1]
    dx <- 1/60 #0.025/3*2 #bin size
    LL <- 1/60 #0.025/3*2 # lower limit of the breaks interval
    UL <- 59/60 # upper limit of the breaks interval
    breaks_interval <- seq(LL-dx,UL+dx,by=dx)
    mids_interval <- seq(LL-dx/2, UL+dx/2, by=dx)
    cl <- makeCluster(number_of_cores)  # 'number_of_cores' is the number of clusters (cores)
    logLikelihood <- parSapply(cl, mids_interval, compute_logLikelihood, distr = SimPoints, beta = beta, alpha = alpha, Time = Time)
    stopCluster(cl)
    
    likelihood <- exp(logLikelihood - mean(logLikelihood))  # translation of the logLikelihood to remove impossibilities
    
    
    for (j in 1:length(rate_v)){
      
      
      for (i in 1:length(data_frame1)){ # 
        out_distance <- as.numeric(data_frame1[[i]])
        dist <- as.character(data_frame2[[i]])
        
        
        rate <- rate_v[j]
        d <- quantile(out_distance, probs = rate, na.rm = TRUE)
        
        # Select parameter values corresponding to distances 
        # less than threshold distance d
        abc_posterior <- ParameterSample_values[which(out_distance < d)]
        
        ## kullback–Leibler divergence
        
        ### With Classical HISTOGRAM
        h <- hist(abc_posterior,breaks = breaks_interval, plot = FALSE) # No need 
        ## of setting 'probability=TRUE' when used 'plot=FALSE'
        density_ABC <- h$density
        
        
        posterior <- likelihood # For uniform prior in [0,1]
        exact_posterior <- posterior
        exact_posterior <- exact_posterior/sum(exact_posterior)/(dx) #Normalization
        density_LIK <- exact_posterior
        
        ### Penalise zero (0) values of the density
        eps <- 1e-10
        for (v in 1:length(density_ABC)){
          if (density_ABC[v]==0){
            density_ABC[v]=eps
          }
          if (density_LIK[v]==0){
            density_LIK[v]=eps
          }
        }
        
        ## We can find the KLD using hist densities
        KLD_value <- KLD(density_LIK,density_ABC)  ## KLD(Likelihood, ABC)
        
        dist_used = paste0("Dist : ",dist)
        threshold_dist = bquote("Threshold "*epsilon == .(round(d,3)))
        sel_rate = paste0("Sel Rate = ",round(rate*100,1),"%")
        dkl_LIK_ABC <- paste0("KLD(LIK,ABC) = ",round(KLD_value[[4]], 3))
        dkl_ABC_LIK <- paste0("KLD(ABC,LIK) = ",round(KLD_value[[5]], 3))
        average <- paste0("Post Mean = ",round(mean(abc_posterior),3))
        stand_dev <- paste0("Post STD = ",round(std(abc_posterior),3))
        
        min_h1 <- min(abc_posterior)
        max_h1 <- max(abc_posterior)
        dx_h1 <- 0.001
        h1 <- hist(abc_posterior, breaks = seq(min_h1-dx_h1, max_h1+dx_h1, dx_h1), plot = FALSE) 
        ymax <- max(max(h1$density), max(exact_posterior))
        
        #create plot
        fig <- function(){
          h2 <- hist(abc_posterior, probability=TRUE, ylim=c(0,ceiling(1.35*ymax)), 
                     col="darkolivegreen", border="#333333", xlab="", ylab="", main="", xlim=c(0,1))
          #          col="darkolivegreen", border="#333333", xlab="", ylab="", main="", xlim=c(0,1))
          lines(h$mids,exact_posterior,col='red',type="l", lwd=2)
          # abline(v=truelambda, col="blue", lty=1, lwd=2)
          segments(x0=truelambda,y0=0,x1=truelambda,y1=ymax,col='blue',lwd=2)
          text(-0.015, ymax*1.25, dkl_ABC_LIK , pos = 4, col = "#330000")
          text(-0.015, ymax*1.15, average , pos = 4, col = "#330000")
          text(-0.015, ymax*1.05, stand_dev , pos = 4, col = "#330000")
          ###return(hist_gr) 
          
        }
        
        fig()
        mtext(paste0(dist),side=3,line=0,outer=FALSE,las=1, col="darkblue", lwd = 0.5)        
      }
      mtext(expression(Lambda[0]),side=1,line=0,outer=TRUE,las=0)
      mtext("Density",side=2,line=0,outer=TRUE,las=0)
      mtext(paste0("\n T = ",Time,"\n Sel Rate = ",round(rate*100,1),"%"),side=3,line=2,outer=FALSE,las=1, col="darkred")
    }
    
  }
  dev.off()
}
