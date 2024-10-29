# Required libraries
library(plotly)

library(GA)    

library(MittagLeffleR)
library(pracma)

#library(webshot)
#webshot::install_phantomjs()

Time <- c(60,120,180)
seedv <- c(3, 739, 123)

for(s in 2:2){
  Ti <- Time[s]
  set.seed(seedv[s])
  epsilon <- 1e-10
  lambda1 <- lambda <- 0.9
  alpha1 <- alpha <- 0.5
  beta <- 0.5
  gamma <- 0.5     
  
  SimPoints <- c(0)
  t <- 0
  n <- 0
  while (t<Ti){
    M <- lambda + sum(alpha*mlf(-gamma*(t+epsilon-SimPoints)^beta,beta,beta,1)*gamma*(t+epsilon-SimPoints)^(beta-1))
    dt <- rexp(1,M)
    t <- t+dt
    U <- runif(1)
    if (U<(((lambda+sum(alpha*mlf(-gamma*(t-SimPoints)^beta,beta,beta,1)*gamma*(t-SimPoints)^(beta-1))))/M)) 
    {
      n <- n+1
      SimPoints <- c(SimPoints,t)		
    }
  }
  #SimPoints
  
  
  ### Likelihood function
  lik_fun <- function(gamma, beta, lambda=lambda1, alpha=alpha1, t = SimPoints) {
    # beta in (0,1) , alpha in (0,1),  lambda>0, gamma>0
    library(MittagLeffleR)
    library(pracma)
    k = length(t)
    logL1 = rep(0,length(beta))
    for (m in 1:length(beta)){
      for (i in 1:k){
        si = 0
        if (i == 1){
          si = 0
        } else {
          for(j in 1:(i-1)){
            #        for (m in 1:length(beta)){
            si = si + ((t[i]-t[j])^(beta[m]-1))*mlf(-gamma[m]*(t[i]-t[j])^beta[m],beta[m],beta[m],1)
            #        }
          }}
        
        logL1[m] = logL1[m] + log(lambda + alpha*gamma[m]*si) - alpha*(1 - mlf(-gamma[m]*(t[k]-t[i])^beta[m],beta[m],1,1))
      }
    }
    
    logL <- logL1 - lambda*Ti
    #  logL
    L <- exp(logL)
    return(L)
  }
  
  # Variables
  gamma <- seq(0.01, 1, by=.01)   
  beta <- seq(0.01, 1, by=.01)   
  
  
  # Generate data for the plot
  
  
  likelihood <- t(outer(gamma, beta, lik_fun))
  
  # Save as PDF file the plot
  
  # Create the plot
  
  #x11(width=20.5,height=13.5)
  font_x = list(
    family = "Courier New",
    size = 30,
    color = "RebeccaPurple")
  
  fig <- plot_ly(x=gamma,y=beta,z=likelihood,type = "contour",colorbar = list(title = "Likelihood", tickfont = list(size = 30))) |>
    layout(xaxis=list(title="Gamma",range = c(0,1), tickfont = list(size = 30)),font=font_x,yaxis=list(title=list(text="Beta"),range = c(0,1), tickfont = list(size = 30)))
  fig
  plotly::export(p = fig, file = paste0("Contour_plot_T",Ti,"_i",d,".png")) 
}     
