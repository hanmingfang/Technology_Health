#Solve for static equilibrium 

# Libraries ---------------------------------------------------------------

library(ggplot2, quietly = TRUE)
library(latex2exp)
library(truncnorm)
library(profvis)
library(numDeriv)
library(rootSolve)
library(microbenchmark)
library(nloptr)
#library(benchmarkme)
library(distrEx)
library(memoise)
library(nleqslv)

#Set my directory
dir = '~/Technology_Health/'
setwd(dir)

# Parameters --------------------------------------------------------------
#Household 
sH = 5              #Skill High type 
sL = 1              #Skill Low type
lambda_gH = 0.25    #Measure healthy workers High skill
lambda_gL = 0.25    #Measure healthy workers Low skill
lambda_bH = 0.25    #Measure unhealthy workers High skill
lambda_bL = 0.25    #Measure unhealthy workers Low skill
theta_L = 0         #Domain for theta
theta_H = Inf
m_L = 0             #Domain for medical expenditure shocks
m_F = Inf
shape_gH = 1.6         #Shape parameter of theta (Gamma distribution) High skill
shape_gL = 1         #Shape parameter of theta (Gamma distribution) Low skill
shape_bH = 0.6       #Shape parameter of theta (Gamma distribution) High skill
shape_bL = 0.1       #Shape parameter of theta (Gamma distribution) Low skill
                     #For a given scale parameter, higher shape parameter means more risk averse households
scale_gH = 1         #Scale parameter of theta (Gamma distribution) High skill
scale_gL = 1         #Scale parameter of theta (Gamma distribution) Low skill
scale_bH = 1         #Scale parameter of theta (Gamma distribution) High skill
scale_bL = 1         #Scale parameter of theta (Gamma distribution) Low skill
rate_g = 1           #Rate for exponential distribution for medical exp. healthy workers  (Mean=1/rate)
rate_b = 0.25        #Rate for exponential distribution for medical exp. unhealthy workers
P_0g =  0.7          #Probability of 0 medical expenditure for healthy worker 
P_0b =  0.5          #Probability of 0 medical expenditure for unhealthy worker
theta_ins_final = 10 #As we can not evaluate f(Inf) I use an upper bound number, but allowing uniroot to extend it in theta_ins
#Firm
N = 1               #Range of tasks (upper limit)
eta = 0.5           #Distribution parameter of the CES
rho = 0.8           #Relative labor productivity of unhealthy workers
psi = 1             #Price of intermediates
sigma = 2           #Elasticity of substitution between tasks
zeta = 2            #Elasticity of substitution between factors (if fixed), just to define zeta_elas
#Change this to a small positive number if equilibrium is not found
C_IN = 0.01          #Health Insurance Fixed Cost (we can start with a very low one)
A = 1               #Parameter in labor productivity
A_0 = 1             #Parameter in labor productivity
delta_H =  1.5       #Parameter in labor productivity of High skill type
lambda_d = 10       #Parameter in sorting function
alpha_d = 5         #Parameter in sorting function
D = 1               #Parameter in Automation Cost function
s_ccp = 1           #Parameter to scale ccps 
tol = 1e-8          #Tolerance for unitroot, affects computation time
K = 1               #Capital stock in the economy

# Primitive Functions ---------------------------------------------------------------
#Distribution objects
#Distribution for Positive part of Medical expenditure
D_mg = Exp(rate_g)
D_mb = Exp(rate_b)
#Distribution for Risk aversion parameter
D_theta_gH = Gammad(scale = scale_gH,shape = shape_gH)
D_theta_gL = Gammad(scale = scale_gL,shape = shape_gL)
D_theta_bH = Gammad(scale = scale_bH,shape = shape_bH)
D_theta_bL = Gammad(scale = scale_bL,shape = shape_bL)
#Conditional cdf of Households' risk aversion parameter
F_gH  = function(theta){ 
  aux = p(D_theta_gH)(theta) #Expectation is shape*rate, Good health
  return(aux)
}
F_gL  = function(theta){ 
  aux = p(D_theta_gL)(theta) #Expectation is shape*rate, Good health
  return(aux)
} 
F_bH  = function(theta){
  aux = p(D_theta_bH)(theta) #Bad health
  return(aux)
}
F_bL  = function(theta){
  aux = p(D_theta_bL)(theta) #Bad health
  return(aux)
}
#Conditional cdf of Medical Expenditure (Positive part!)
H_g  = function(m){
  aux =  p(D_mg)(m) #Good health
  return(aux)
}     
H_b  = function(m){
  aux =  p(D_mb)(m) #Bad health
  return(aux)
}
#Parametrized functions
#Labor productivity (per skill unit)
gamma_prod = function(s,i){   
  if(s == sL){aux = A_0*exp(A*i)}
  else{aux = A_0*exp(A*i)*exp(delta_H*i)}
  return(aux)
}
#Sorting of workers
#Primitive function
sorting_f = function(i){
  #aux = 1
  #aux = i
  aux = exp(lambda_d*i - alpha_d)/(1+exp(lambda_d*i - alpha_d))
  return(aux)
}
#Normalizing constant
norm_const_delta_sort = integrate(Vectorize(sorting_f),lower = N-1, upper = N)$value 
#Sorting function
delta_sort = function(i){
  aux = sorting_f(i)/norm_const_delta_sort
  return(aux)
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
#Normalizing constant in production function
B = function(i){
  zetai = zeta_elas(i)
  aux = (1-eta)^(zetai/(1-zetai))
  return(aux)
}
#Automation Fixed Cost
C_A = function(i){
  aux = 0.1*exp(D*i)
  return(aux)
}

# Equilibrium solutions -------------------------------------------------------------------

#Sourcing equations from the model
source("src/Static_Equilibrium/Household_fcns.R")
source("src/Static_Equilibrium/Firm_fcns.R")
source("src/Static_Equilibrium/Equilibrium_fcns.R")

#Non linear equation solver (is very fast but a little more sensitive to initial conditions)
#Using exponentials in the equations to ensure positivity
#Change method to "Newton" if "Broyden" doesn't converge.
#If the algorithm doesn't find a better point, try decreasing C_IN
nles_sol = nleqslv(x = c(log(14),log(2),log(12),log(0.5),log(1.6),log(40)), 
                   fn = F_zeros, jac=NULL, method = "Broyden", jacobian=FALSE,
                   control = list("allowSingular"=TRUE), global="dbldog")
nles_sol

w0H = exp(nles_sol$x[1])
w0L = exp(nles_sol$x[2])
w1H = exp(nles_sol$x[3])
w1L = exp(nles_sol$x[4])
R = exp(nles_sol$x[5])
Y = exp(nles_sol$x[6])
val = c(w0H,w0L,w1H,w1L,R,Y)
val
#Checking for wieghted wages 
L1_H = L1_s('g',sH,w0H,w1H)+L1_s('b',sH,w0H,w1H)
L1_L = L1_s('g',sL,w0L,w1L)+L1_s('b',sL,w0L,w1L)
w1 = (w1H*L1_H + w1L*L1_L)/(L1_H+L1_L) 

L0_H = L0_s('g',sH,w0H,w1H)+L0_s('b',sH,w0H,w1H)
L0_L = L0_s('g',sL,w0L,w1L)+L0_s('b',sL,w0L,w1L)
w0 = (w0H*L0_H + w0L*L0_L)/(L0_H+L0_L) 
c(w0,w1)
#Checking if w1>w0
w0-w1






