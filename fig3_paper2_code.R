#set directory as folder containing current file
#install.packages("rstudioapi")
library(rstudioapi)
cur_dir = dirname(getSourceEditorContext()$path)
setwd(cur_dir)

library(plot3D)
library(plotly)

#library(GA)    

library(MittagLeffleR)
library(pracma)

#library(webshot)
#webshot::install_phantomjs()

Time <- c(10,20,30)
seed <- 321

#windows(width = 23.5, height = 13.5)
par(mar = c(4, 4, 1, 2))        # Reduce space around plots
par(mfrow = c(1, 3),oma=c(.5,0,0,0))

for(s in 1:3){
  Ti <- Time[s]
  set.seed(seed)
  epsilon <- 1e-10
  lambda1 <- lambda <- 0.5
  alpha <- 0.5
  beta <- 0.5
  gamma1 <- gamma <- 0.9
  
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
  lik_fun <- function(alpha, beta, lambda=lambda1, gamma=gamma1, t = SimPoints) {
    # beta in (0,1) , alpha in (0,1),  lambda>0, gamma>0
    library(MittagLeffleR)
    library(pracma)
    k = length(t)
    logL1 = rep(0,length(beta))
    for (n in 1:length(beta)){
      for (i in 1:k){
        si = 0
        if (i == 1){
          si = 0
        } else {
          for(j in 1:(i-1)){
            #        for (m in 1:length(beta)){
            si = si + ((t[i]-t[j])^(beta[n]-1))*mlf(-gamma*(t[i]-t[j])^beta[n],beta[n],beta[n],1)
            #        }
          }}
        
        logL1[n] = logL1[n] + log(lambda + alpha[n]*gamma*si) - alpha[n]*(1 - mlf(-gamma*(t[k]-t[i])^beta[n],beta[n],1,1))
      }
    }
    
    logL <- logL1 - lambda*Ti
    #  logL
    L <- exp(logL)
    return(L)
  }
  
  # Variables
  alpha <- seq(0, .9, by=.1)   
  beta <- seq(0, 1, by=.1)   
  
  
  # Generate data for the plot
  
  
  likelihood <- t(outer(alpha, beta, lik_fun))
  
  # Save as PDF file the plot
  
  # Create the plot
  
  #x11(width=20.5,height=13.5)
  font_x = list(
    family = "Courier New",
    size = 30,
    color = "RebeccaPurple")
  
  fig <- plot_ly(x=alpha,y=beta,z=likelihood,type = "contour",colorbar = list(title = "Likelihood", tickfont = list(size = 30))) |>
    layout(xaxis=list(title="Alpha",range = c(0,.9), tickfont = list(size = 30)),font=font_x,yaxis=list(title=list(text="Beta"),range = c(0,1), tickfont = list(size = 30)))
  fig
  plotly::export(p = fig, file = paste0("Contour_plot_T",Ti,"_i",d,".png")) 
}

