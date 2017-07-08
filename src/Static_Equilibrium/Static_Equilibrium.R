#Solve for static equilibrium

# Libraries ---------------------------------------------------------------

library(stargazer, quietly = TRUE)
library(lmtest, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(dplyr)
library(reshape2)
library(ggfortify)
library(data.table)
library(RMySQL)
library(plotly)
library(latex2exp)
library(truncnorm)
Sys.setenv("plotly_username"="tlarrouc")
Sys.setenv("plotly_api_key"="bdqUsPxMwTYk9lI7KmY4")
dir = '~/Dropbox/Technology/Codes/Robots/'
setwd(dir)


# Parameters --------------------------------------------------------------
#Household 
xi = 0.455          #Utility function
phi = 1             #Cost of labor effort  
lambda_g = 0.5      #Measure healthy workers
theta_L = 1.001     #Domain for theta
theta_H = 4
m_L = 0             #Domain for medical expenditure shocks
m_F = Inf
mu_g = 2            #Mean of theta for healthy workers
mu_b = 1.001        #Mean of theta for healthy workers
sd_g = 1            #Standard deviation of theta for healthy workers
sd_b = 1            #Standard deviation of theta for unhealthy workers
rate_g = 1          #Rate for exponential distribution for medical exp. healthy workers  
rate_b = 0.5        #Rate for exponential distribution for medical exp. unhealthy workers
util_min = 0.001    #Minimum consumption minus labor effort a household can have (can't be 0 or blows up)
#Firm
N = 1               #Range of tasks (upper limit)
eta = 0.5           #Distribution parameter of the CES
rho = 0.9           #Relative labor productivity of unhealthy workers
psi = 1             #Price of intermediates
sigma = 2           #Elasticity of substitution between tasks
zeta = 2            #Elasticity of substitution between factors
C_IN = 10           #Health Insurance Fixed Cost
A = 1               #Parameter in labor productivity
lambda_d = 1        #Parameter in sorting function
alpha_d = 1         #Parameter in sorting function
D = 1               #Parameter in Automation Cost function

# Primitive Functions ---------------------------------------------------------------
#Distributions
#Conditional cdf of Households' types
F_g  = function(theta){
  aux =  ptruncnorm(theta, a = theta_L, b = theta_H, mean = mu_g, sd = sd_g) #Good health
  return(aux)
} 
F_b  = function(theta){
  aux =  ptruncnorm(theta, a = theta_L, b = theta_H, mean = mu_b, sd = sd_b) #Bad health
  return(aux)
} 
#Conditional pdf of Medical Expenditure
h_g  = function(m){
  aux =  dexp(m, rate = rate_g) #Good health
  return(aux)
}     
h_b  = function(m){
  aux =  dexp(m, rate = rate_b) #Bad health
  return(aux)
}
#Conditional cdf of Medical Expenditure
H_g  = function(m){
  aux =  pexp(m, rate = rate_g) #Good health
  return(aux)
}     
H_b  = function(m){
  aux =  pexp(m, rate = rate_b) #Bad health
  return(aux)
}
#Parametrized functions
#Labor productivity 
gamma_prod = function(i){   
  aux = exp(A*i)
  return(aux)
}   
#Sorting of workers
delta_sort = function(i){
  exp(lambda_d*i - alpha_d)/(1+exp(lambda_d*i - alpha_d))
}
#Capital productivity
z_prod = function(i){
  aux = 1
  return(aux)
}
#Function for elasticity of subst
zeta_elas = function(i){
  aux = zeta
  return(aux)
}
#Automation Fixed Cost
C_A = function(i){
  aux = exp(D*i)
  return(aux)
}

# Equilibrium Equations ---------------------------------------------------
#Households' optimal choices
#Individual labor suply for No insurance
l0_s = function(w0){
  aux = (w0/phi)^(1/xi)
  return(aux)
}
#Individual labor suply for insurance
l1_s = function(w1){
  aux = (w1/phi)^(1/xi)
  return(aux)
}
#Reservation wage
w_res = function(w0){
  aux = w0/(1+xi)
  return(aux)
}
#Utility function under  insurance
u1 = function(theta,h,w1){
  l1 = l1_s(w1)
  aux = (1/(1-theta))*((w1*l1 - phi*((l1^(1+xi))/(1+xi)))^(1-theta)-1)
  return(aux)
}
#Utility function under No insurance
u0 = function(theta,h,w0,m){
  l0 = l0_s(w0)
  aux = vector(length = length(m)) #Adapting the function to retrieve a vector (to use integrate function later)
  for(i in 1:length(aux)){ #Here Im taking the max between the argument of the utility and util_min
    aux[i] = (1/(1-theta))*((max(w0*l0 - m[i] - phi*((l0^(1+xi))/(1+xi)),util_min))^(1-theta)-1)  
  }
  return(aux)
}
#Expected utility under no insurance
E_u0 = function(theta,h,w0){
  l0 = l0_s(w0)
  integrand_u0 = function(m){ #Creating a function that a returns a vectorized integrand
    aux = vector(length = length(m))
    if(h == 'g'){ #If healthy worker
      for(i in 1:length(m)){
        aux[i] = u0(theta,h,w0,m = m[i])*h_g(m[i]) #Utility times pdf of healthy, for each medical shock
      }
    }
    else{
      for(i in 1:length(m)){
        aux[i] = u0(theta,h,w0,m = m[i])*h_b(m[i]) #Utility times pdf of unhealthy, for each medical shock
      }
    }
    return(aux)
  }
  integral = integrate(integrand_u0, lower = m_L, upper = m_F)
  return(integral$value) #return just the value of the inetgral
}
#Threshold for household insurance decision
theta_ins = function(h,w0,w1){
  initial = theta_L  #Search over
  final = theta_H
  fun = function (theta) E_u0(theta,h,w0) - u1(theta,h,w1)   #As a function of just theta, this is a decreasing function of theta
  if(E_u0(theta_L,h,w0) - u1(theta_L,h,w1) < 0){aux = theta_L} #If at lower bound is negative
  else if(E_u0(theta_H,h,w0) - u1(theta_H,h,w1) > 0){aux = theta_H}
  else{aux = uniroot(fun, c(initial,final), tol = 1e-13, extendInt = "downX")$root}  #Get the root, downX is to tell that is decresing on theta
    return(aux) 
}


# Plots -------------------------------------------------------------------

#Plot to show that FOSD in Assumption 1 holds for this case
ggplot(data.frame(x=c(0, 15)), aes(x=x)) + 
  stat_function(fun=H_g, geom="line", aes(colour = "H_g")) + xlab("x") + 
  ylab("y") + stat_function(fun=H_b, geom="line",aes(colour = "H_b")) 
#Plot to show that FOSD in Assumption 2 holds for this case
ggplot(data.frame(x=c(1.001, 4)), aes(x=x)) + 
  stat_function(fun=F_g, geom="line", aes(colour = "F_g")) + xlab("x") + 
  ylab("y") + stat_function(fun=F_b, geom="line",aes(colour = "F_b")) 
#Plot gamma_prod
ggplot(data.frame(x=c(N-1,N)), aes(x=x)) + 
  stat_function(fun=gamma_prod, geom="line") + xlab("x") + ylab("y") 
#Plot delta_sort
ggplot(data.frame(x=c(N-1,N)), aes(x=x)) + 
  stat_function(fun=delta_sort, geom="line") + xlab("x") + ylab("y") 
#Plot C_A
ggplot(data.frame(x=c(N-1,N)), aes(x=x)) + 
  stat_function(fun=C_A , geom="line") + xlab("x") + ylab("y") 











