#set directory as folder containing current file
#install.packages("rstudioapi")
library(rstudioapi)
cur_dir = dirname(getSourceEditorContext()$path)
setwd(cur_dir)

library(GA) 
library(plot3D)
 
library(MittagLeffleR)
library(pracma)



Ti <- 120
seed <- 1234
set.seed(seed)
epsilon <- 1e-10
lambda1 <- lambda <- 0.9
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

#fun <- ll_fun  #log-likelihood function
beta <- seq(.2, .6, by=.01)   # beta (x)
alpha <- seq(.2, .6, by=.01)   # alpha (y)

# Generate data for the plot

likelihood <- outer(alpha, beta, lik_fun)
VV <- V2 <- likelihood


#windows(width = 23,height=13)

#par(mfrow = c(1, 2))
par(mar = c(4, 4, 1, 2))        # Reduce space around plots
par(mfrow = c(1, 2),oma=c(.5,0,0,0))
## 2D Plot
image2D(x=alpha,y=beta,z=likelihood,type = "contour",#color=I("blue"),
        colkey = list(side = 4, length = 0.5), clab = "z",
        colorbar = list(title = "z", tickfont = list(size = 20))) 

## 3D Plot
persp3D(x=beta, y=alpha, z=VV,
        colvar=V2, phi = 40, theta = 295, nticks=5, 
        zlab="z", ylab = "y", xlab="x", colkey = FALSE, #colkey = list(side = 3, length = 0.5),
        box = TRUE, border = NA, shade = .1,ticktype="simple",
        facets = TRUE, scale= TRUE, expand = 0.8,
        colorbar = list(title = "z", tickfont = list(size = 20)))
