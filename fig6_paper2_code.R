#set directory as folder containing current file
#install.packages("rstudioapi")
library(rstudioapi)
cur_dir = dirname(getSourceEditorContext()$path)
setwd(cur_dir)

library(plotly)

library(GA)

library(MittagLeffleR)
library(pracma)

#library(webshot)
#webshot::install_phantomjs()

Time <- c(60,120,180)
seedv <- c(112233, 112233, 122) 

for(s in 1:3){
  Ti <- Time[s]
  set.seed(seedv[s])
  epsilon <- 1e-10
  lambda <- 0.5
  alpha1 <- alpha <- 0.5
  beta1 <- beta <- 0.8
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
  lik_fun <- function(lambda, gamma, beta=beta1, alpha=alpha1, t = SimPoints) {
    # beta in (0,1) , alpha in (0,1),  lambda>0, gamma>0
    library(MittagLeffleR)
    library(pracma)
    k = length(t)
    
    logL1 = 0
    for (i in 1:k){
      si = 0
      if (i == 1){
        si = 0
      } else {
        for(j in 1:(i-1)){
          si = si + ((t[i]-t[j])^(beta-1))*mlf(-gamma*(t[i]-t[j])^beta,beta,beta,1)
        }}
      logL1 = logL1 + log(lambda + alpha*gamma*si) - alpha*(1 - mlf(-gamma*(t[k]-t[i])^beta,beta,1,1))
    }
    
    logL <- logL1 - lambda*Ti
    #  logL
    L <- exp(logL)
    return(L)
  }
  
  # Variables
  lambda <- seq(0.025, 1, by=.025)   
  gamma <- seq(0, 1, by=.025)   
  
  
  # Generate data for the plot
  
  likelihood <- t(outer(lambda, gamma, lik_fun))
  
  # Save as PDF file the plot
  
  # Create the plot
  font_x = list(
    family = "Courier New",
    size = 30,
    color = "RebeccaPurple")
  fig <- plot_ly(x=lambda,y=gamma,z=likelihood,type = "contour",colorbar = list(title = "Likelihood")) |>
    layout(xaxis=list(title="Lambda",range = c(0,1)),font=font_x,yaxis=list(title=list(text="Gamma"),range = c(0,1)))
  fig
  plotly::export(p = fig, file = paste0("Contour_plot_T",Ti,"_i",d,".png")) 
}
