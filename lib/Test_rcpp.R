#Solve for static equilibrium (Less risk averse households)

# Libraries ---------------------------------------------------------------

library(stargazer, quietly = TRUE)
library(lmtest, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(dplyr)
library(reshape2)
library(ggfortify)
library(data.table)
library(plotly)
library(latex2exp)
library(truncnorm)
library(profvis)
library(numDeriv)
library(rootSolve)
library(microbenchmark)
#Change directory
dir = '~/Dropbox/Technology/Codes/Technology_Health/'
setwd(dir)

# Parameters --------------------------------------------------------------
#Household 
xi = 0.455          #Utility function
phi = 1             #Cost of labor effort  
lambda_g = 0.5      #Measure healthy workers
theta_L = 0         #Domain for theta
theta_H = 10
m_L = 0             #Domain for medical expenditure shocks
m_F = Inf
mu_g = 2            #Mean of theta for healthy workers
mu_b = 0            #Mean of theta for healthy workers
sd_g = 1            #Standard deviation of theta for healthy workers
sd_b = 1            #Standard deviation of theta for unhealthy workers
rate_g = 2          #Rate for exponential distribution for medical exp. healthy workers  
rate_b = 1          #Rate for exponential distribution for medical exp. unhealthy workers
util_min = 0.001    #Minimum consumption minus labor effort a household can have (can't be 0 or blows up)
#This is not in the same way in the document!
#Firm
N = 1               #Range of tasks (upper limit)
eta = 0.5           #Distribution parameter of the CES
rho = 0.8           #Relative labor productivity of unhealthy workers
psi = 1             #Price of intermediates
sigma = 2           #Elasticity of substitution between tasks
zeta = 2            #Elasticity of substitution between factors (if fixed), just to define zeta_elas
C_IN = 0.5          #Health Insurance Fixed Cost
A = 1               #Parameter in labor productivity
A_0 = 1             #Parameter in labor productivity
lambda_d = 10        #Parameter in sorting function
alpha_d = 5         #Parameter in sorting function
D = 1               #Parameter in Automation Cost function
tol = 1e-3          #Tolerance for unitroot, affects computation time
K = 0.2             #Capital stock in the economy
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
  aux = A_0*exp(A*i)
  return(aux)
}   
#Sorting of workers
#TODO: edit this function
delta_sort = function(i){
  aux = exp(lambda_d*i - alpha_d)/(1+exp(lambda_d*i - alpha_d))
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

# Equilibrium Equations ---------------------------------------------------
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
  fun = function (theta) E_u0(theta,h,w0) - u1(theta,h,w1)   #This is a decreasing function of theta
  if(E_u0(theta_L,h,w0) - u1(theta_L,h,w1) < 0){aux = theta_L} #If at lower bound is negative
  else if(E_u0(theta_H,h,w0) - u1(theta_H,h,w1) > 0){aux = theta_H}
  else{aux = uniroot(fun, c(initial,final), tol = tol, extendInt = "downX")$root}  #Get the root, 
  #downX is to tell that is decresing on theta, so can look further than the specified range, 
  #although Im not using this given the boudary cases
  return(aux) 
}
#Aggregate labor supply for no insurance
L0_s = function(h,w0,w1){
  if(h == 'g'){aux = lambda_g*l0_s(w0)*F_g(theta_ins(h,w0,w1))} # L^0_g
  else{aux = (1-lambda_g)*l0_s(w0)*F_b(theta_ins(h,w0,w1))} # L^0_b
  return(aux)
}
#Aggregate labor supply for insurance
L1_s = function(h,w0,w1){
  if(h == 'g'){aux = lambda_g*l1_s(w1)*(1-F_g(theta_ins(h,w0,w1)))} # L^1_g
  else{aux = (1-lambda_g)*l1_s(w1)*(1-F_b(theta_ins(h,w0,w1)))} # L^1_b
  return(aux)
}
#Endogenous proportion of healthy workers for no health insurance
#TODO: Edit this function
Chi_0g = function(w0,w1,i){
  aux = 1*(L0_s('g',w0,w1)/(L0_s('g',w0,w1) + L0_s('b',w0,w1)))
  return(aux)
}
#Endogenous proportion of healthy workers for health insurance
#TODO: edit this function
Chi_1g = function(w0,w1,i){
  aux =  delta_sort(i)*(L1_s('g',w0,w1)/(L1_s('g',w0,w1) + L1_s('b',w0,w1)))
  return(aux)
}
#Expected expenditure shock
#TODO: edit this function, just for Exponentialnow
E_m = function(h){
  if(h == 'g'){aux = 1/rate_g}
  else{aux = 1/rate_b}
  return(aux)
}
#Expected firm's medical expenditure
M = function(w0,w1,i){
  aux = (E_m('g')*Chi_1g(w0,w1,i)+E_m('b')*(1-Chi_1g(w0,w1,i)))/l1_s(w1)
  return(aux)
}
#Average labor productivity for contract without health insurance
gamma_prod_bar_0 = function(w0,w1,i){
  aux = gamma_prod(i)*((1-rho)*Chi_0g(w0,w1,i)+rho)
  return(aux)
}
#Average labor productivity for contract with health insurance
gamma_prod_bar_1 = function(w0,w1,i){
  aux = gamma_prod(i)*((1-rho)*Chi_1g(w0,w1,i)+rho)
  return(aux)
}
#Effective wage without health insurance
w_hat0 = function(w0,w1,i){
  aux = w0/gamma_prod_bar_0(w0,w1,i)
  return(aux)
}
#Effective wage with health insurance
w_hat1 = function(w0,w1,i){
  aux = (w1+M(w0,w1,i))/gamma_prod_bar_1(w0,w1,i)
  return(aux)
}
#Effective price of capital
R_hat = function(R,i){
  aux = R/z_prod(i)
  return(aux)
}
#Sufficient statistic for effective prices without insurance
p_0 = function(w0,w1,i){
  aux = (eta*w_hat0(w0,w1,i))/((1-eta)*psi)
  return(aux)
}
#Sufficient statistic for effective prices with insurance
p_1 = function(w0,w1,i){
  aux = (eta*w_hat1(w0,w1,i))/((1-eta)*psi)
  return(aux)
}
#Sufficient statistic for effective prices with capital
p_k = function(R,i){
  aux = (eta*R_hat(R,i))/((1-eta)*psi)
  return(aux)
}
#Conditional output function without health insurance
y_0 = function(w0,w1,i,Y){
  p0i = p_0(w0,w1,i)
  zetai = zeta_elas(i)
  w_hat0i = w_hat0(w0,w1,i)
  aux = (Y*((sigma-1)/sigma)^sigma)*((B(i)*(eta*p0i^(zetai-1)+(1-eta))^(zetai/(zetai-1)))/(w_hat0i + psi*p0i^zetai))^sigma
  return(aux)
}
#Conditional output function with health insurance
y_1 = function(w0,w1,i,Y){
  p1i = p_1(w0,w1,i)
  zetai = zeta_elas(i)
  w_hat1i = w_hat1(w0,w1,i)
  aux = (Y*((sigma-1)/sigma)^sigma)*((B(i)*(eta*p1i^(zetai-1)+(1-eta))^(zetai/(zetai-1)))/(w_hat1i + psi*p1i^zetai))^sigma
  return(aux)
}
#Conditional output function with capital
y_k = function(R,i,Y){
  pki = p_k(R,i)
  zetai = zeta_elas(i)
  R_hati = R_hat(R,i)
  aux = (Y*((sigma-1)/sigma)^sigma)*((B(i)*(eta*pki^(zetai-1)+(1-eta))^(zetai/(zetai-1)))/(R_hati + psi*pki^zetai))^sigma
  return(aux)
}
#Conditional effective labor wihout health insurance 
l_hat0 = function(w0,w1,i,Y){
  zetai = zeta_elas(i)
  p0i = p_0(w0,w1,i)
  aux = y_0(w0,w1,i,Y)/(B(i)*(eta*p0i^(zetai-1)+(1-eta))^(zetai/(zetai-1)))
  return(aux)
}
#Conditional effective labor with health insurance 
l_hat1 = function(w0,w1,i,Y){
  zetai = zeta_elas(i)
  p1i = p_1(w0,w1,i)
  aux = y_1(w0,w1,i,Y)/(B(i)*(eta*p1i^(zetai-1)+(1-eta))^(zetai/(zetai-1)))
  return(aux)
}
#Conditional capital 
k = function(R,i,Y){
  zetai = zeta_elas(i)
  pki = p_k(R,i)
  aux = y_k(R,i,Y)/(B(i)*(eta*pki^(zetai-1)+(1-eta))^(zetai/(zetai-1)))
  return(aux)
} 
#Conditional intermediates without health insurance
q_0 = function(w0,w1,i,Y){
  zetai = zeta_elas(i)
  p0i = p_0(w0,w1,i)
  aux = (y_0(w0,w1,i,Y)*p0i^zetai)/(B(i)*(eta*p0i^(zetai-1)+(1-eta))^(zetai/(zetai-1)))
  return(aux)
} 
#Conditional intermediates with health insurance
q_1 = function(w0,w1,i,Y){
  zetai = zeta_elas(i)
  p1i = p_1(w0,w1,i)
  aux = (y_1(w0,w1,i,Y)*p1i^zetai)/(B(i)*(eta*p1i^(zetai-1)+(1-eta))^(zetai/(zetai-1)))
  return(aux)
} 
#Conditional intermediates with capital
q_k = function(R,i,Y){
  zetai = zeta_elas(i)
  pki = p_k(R,i)
  aux = (y_k(R,i,Y)*pki^zetai)/(B(i)*(eta*pki^(zetai-1)+(1-eta))^(zetai/(zetai-1)))
  return(aux)
} 
#Indifference threshold between insurance and not insurance
X_tilde = function(w0,w1,Y){
  LHS = function(i) C_IN/(((B(i)*(sigma-1)/sigma)^(sigma-1))*Y/sigma)
  #Do not call the same function more than one time if is not neccesary
  L0_s_g_var = L0_s('g',w0,w1)
  L0_s_b_var = L0_s('b',w0,w1)
  Chi_0gi = function(i) 1*(L0_s_g_var)/(L0_s_g_var + L0_s_b_var)
  gamma_prod_bar_0i = function(i) gamma_prod(i)*((1-rho)*Chi_0gi(i)+rho)
  w_hat0i = function(i) w0/gamma_prod_bar_0i(i)
  p_0i = function(i) (eta*w_hat0i(i))/((1-eta)*psi)
  ###
  L1_s_g_var = L1_s('g',w0,w1)
  L1_s_b_var = L1_s('b',w0,w1)
  Chi_1gi = function(i) delta_sort(i)*((L1_s_g_var)/(L1_s_g_var + L1_s_b_var))
  gamma_prod_bar_1i = function(i) gamma_prod(i)*((1-rho)*Chi_1gi(i)+rho)
  E_mg = E_m('g')
  E_mb = E_m('b')
  Mi = function(i) (E_mg*Chi_1gi(i)+E_mb*(1-Chi_1gi(i)))/l1_s(w1)
  w_hat1i = function(i) (w1+Mi(i))/gamma_prod_bar_1i(i)
  p_1i = function(i) (eta*w_hat1i(i))/((1-eta)*psi)
  RHS1 = function(i) (((eta*p_1i(i)^(zeta_elas(i)-1)+(1-eta))^(zeta_elas(i)/(zeta_elas(i)-1)))/(w_hat1i(i)+psi*p_1i(i)^zeta_elas(i)))^(sigma-1)
  RHS0 = function(i) (((eta*p_0i(i)^(zeta_elas(i)-1)+(1-eta))^(zeta_elas(i)/(zeta_elas(i)-1)))/(w_hat0i(i)+psi*p_0i(i)^zeta_elas(i)))^(sigma-1)
  
  fun = function (i) LHS(i)-RHS1(i)+RHS0(i)   #This should be a DECREASING   function of i (but depends on the me
  #dical expenditure)
  initial = N-1 #Lower bound for i
  final = N #Upper boud for i
  #Boundary cases
  if(is.na(w_hat0i(N)) & !is.na(w_hat1i(N))){aux = N-1} #If w_hat0 is NA is because every household works for contract
  #with insurance, thus X_tilde has to be N-1
  else if(!is.na(w_hat0i(N)) & is.na(w_hat1i(N))){aux = N}#The opposite than before
  else if(fun(N)>0){aux = N} #This should be correct given that the function is DECREASING
  else if(fun(N-1)<0){aux = N-1}
  else{aux = uniroot(fun, c(initial,final), tol = tol)$root}
  return(aux) 
}
#Indifference threshold between capital and not insurance
I_tilde0 = function(w0,w1,R,Y){
  LHS = function(i) C_A(i)/(((B(i)*(sigma-1)/sigma)^(sigma-1))*Y/sigma)
  #Do not call the same function more than one time if is not neccesary
  L0_s_g_var = L0_s('g',w0,w1)
  L0_s_b_var = L0_s('b',w0,w1)
  Chi_0gi = function(i) 1*(L0_s_g_var)/(L0_s_g_var + L0_s_b_var)
  gamma_prod_bar_0i = function(i) gamma_prod(i)*((1-rho)*Chi_0gi(i)+rho)
  w_hat0i = function(i) w0/gamma_prod_bar_0i(i)
  p_0i = function(i) (eta*w_hat0i(i))/((1-eta)*psi)
  ###
  R_hati = function(i) R/z_prod(i)
  p_ki = function(i) (eta*R_hati(i))/((1-eta)*psi)
  RHSk = function(i) (((eta*p_ki(i)^(zeta_elas(i)-1)+(1-eta))^(zeta_elas(i)/(zeta_elas(i)-1)))/(R_hati(i)+psi*p_ki(i)^zeta_elas(i)))^(sigma-1)
  RHS0 = function(i) (((eta*p_0i(i)^(zeta_elas(i)-1)+(1-eta))^(zeta_elas(i)/(zeta_elas(i)-1)))/(w_hat0i(i)+psi*p_0i(i)^zeta_elas(i)))^(sigma-1)
  
  fun = function (i) LHS(i)-RHSk(i)+RHS0(i) # Under the assumption that profits conditional on capital are non increasing on i
  #this function should be INCREASING
  initial = N-1 #Lower bound for i
  final = N #Upper boud for i
  #boundary cases
  if(is.na(w_hat0i(N))){aux = N} #If w_hat0 is NA is because every household works for contract
  #with insurance, thus I_tilde0 has to be N (just use capital in this conditional case)
  else if(fun(N)<0){aux = N} #This should be correct given that the function is INCREASING in this case
  else if(fun(N-1)>0){aux = N-1}
  else{aux = uniroot(fun, c(initial,final), tol = tol)$root}
  return(aux) 
}
#Indifference threshold between capital and insurance
I_tilde1 = function(w0,w1,R,Y){
  LHS = function(i) (C_IN-C_A(i))/(((B(i)*(sigma-1)/sigma)^(sigma-1))*Y/sigma)
  ###
  L1_s_g_var = L1_s('g',w0,w1)
  L1_s_b_var = L1_s('b',w0,w1)
  Chi_1gi = function(i) delta_sort(i)*((L1_s_g_var)/(L1_s_g_var + L1_s_b_var))
  gamma_prod_bar_1i = function(i) gamma_prod(i)*((1-rho)*Chi_1gi(i)+rho)
  E_mg = E_m('g')
  E_mb = E_m('b')
  Mi = function(i) (E_mg*Chi_1gi(i)+E_mb*(1-Chi_1gi(i)))/l1_s(w1)
  w_hat1i = function(i) (w1+Mi(i))/gamma_prod_bar_1i(i)
  p_1i = function(i) (eta*w_hat1i(i))/((1-eta)*psi)
  ###
  R_hati = function(i) R/z_prod(i)
  p_ki = function(i) (eta*R_hati(i))/((1-eta)*psi)
  RHSk = function(i) (((eta*p_ki(i)^(zeta_elas(i)-1)+(1-eta))^(zeta_elas(i)/(zeta_elas(i)-1)))/(R_hati(i)+psi*p_ki(i)^zeta_elas(i)))^(sigma-1)
  RHS1 = function(i) (((eta*p_1i(i)^(zeta_elas(i)-1)+(1-eta))^(zeta_elas(i)/(zeta_elas(i)-1)))/(w_hat1i(i)+psi*p_1i(i)^zeta_elas(i)))^(sigma-1)
  
  fun = function (i) LHS(i)-RHS1(i)+RHSk(i) #Different sign than before (consistent with the document) This function is 
  #DECREASING
  initial = N-1 #Lower bound for i
  final = N #Upper bound for i
  #boundary cases
  if(is.na(w_hat1i(N))){aux = N} #If w_hat1 is NA is because every household works for contract
  #without insurance, thus I_tilde1 has to be N (just use capital in this conditional case)
  else if(fun(N)>0){aux = N} #This should be correct given that the function is DECREASING in this case
  else if(fun(N-1)<0){aux = N-1}
  else{aux = uniroot(fun, c(initial,final), tol = tol)$root}
  return(aux) 
}
#Excess demand for Capital
k_excess_d = function(w0,w1,R,Y){
  I0 = I_tilde0(w0,w1,R,Y)
  I1 = I_tilde1(w0,w1,R,Y)
  integrand_k = function(r){ #Creating a function that  returns a vectorized integrand (recieves a vector with 
    #values to be evaluated at
    aux = vector(length = length(r))
    for(i in 1:length(r)){
      aux[i] = k(R,r[i],Y) #integrand capital demanded, Remember that is r[i]
    }
    return(aux)
  }
  integral = integrate(integrand_k, lower = N-1, upper = min(I0,I1)) #Bounds specified in the document
  aux = integral$value-K
  return(aux)  
}
#Excess demand for labor without insurance
#TODO: CHECK THIS FOR NON INTERIOR CASES
l0_excess_d = function(w0,w1,R,Y){
  X = X_tilde(w0,w1,Y)
  I0 = I_tilde0(w0,w1,R,Y)
  I1 = I_tilde1(w0,w1,R,Y)
  if(I0<I1){ #Checking for interior cases
    integrand_l0 = function(r){ #Creating a function that a returns a vectorized integrand
      aux = vector(length = length(r))
      for(i in 1:length(r)){
        aux[i] = l_hat0(w0,w1,r[i],Y)/gamma_prod_bar_0(w0,w1,r[i]) #integrand of l0 demanded
      }
      return(aux)
    }
    integral = integrate(integrand_l0, lower = I0, upper = X) #Bounds specified in the document
    aux = integral$value - (L0_s('g',w0,w1)+L0_s('b',w0,w1)) #subtracting Labor supplied for no insurance
  }
  else{aux = 0 - (L0_s('g',w0,w1)+L0_s('b',w0,w1))}#If there is no labor without insurance in equilibrium
  return(aux)  
}
#Excess demand fo labor with insurance
#TODO: CHECK THIS FOR NON INTERIOR CASES
l1_excess_d = function(w0,w1,R,Y){
  X = X_tilde(w0,w1,Y)
  I0 = I_tilde0(w0,w1,R,Y)
  I1 = I_tilde1(w0,w1,R,Y)
  integrand_l1 = function(r){ #Creating a function that a returns a vectorized integrand
    aux = vector(length = length(r))
    for(i in 1:length(r)){
      aux[i] = l_hat1(w0,w1,r[i],Y)/gamma_prod_bar_1(w0,w1,r[i]) #integrand of l1 demanded
    }
    return(aux)
  }
  integral = integrate(integrand_l1, lower = max(I1,X), upper = N) #Bounds specified in the document
  aux = integral$value - (L1_s('g',w0,w1)+L1_s('b',w0,w1)) #subtracting Labor supplied for insurance
  return(aux)  
}
#Total consumption good
#TODO: Include indicator function, boudarie cases and Review
Y_excess_s = function(w0,w1,R,Y){
  X = X_tilde(w0,w1,Y)
  I0 = I_tilde0(w0,w1,R,Y)
  I1 = I_tilde1(w0,w1,R,Y)  
  Y_k = function(Y){
    integrand_yk = function(r){ #Creating a function that returns a vectorized integrand
      aux = vector(length = length(r))
      for(i in 1:length(r)){
        aux[i] = (y_k(R,r[i],Y))^((sigma-1)/sigma) #integrand of y_k in the consumption good production
      }
      return(aux)
    }
    integral = integrate(integrand_yk, lower = N-1, upper = min(I1,I0))
    aux = integral$value
    return(aux)
  }
  Y_0 = function(Y){
    integrand_y0 = function(r){ #Creating a function that returns a vectorized integrand
      aux = vector(length = length(r))
      for(i in 1:length(r)){
        aux[i] = (y_0(w0,w1,r[i],Y))^((sigma-1)/sigma) #integrand of y_0 in the consumption good production
      }
      return(aux)
    }
    integral = integrate(integrand_y0, lower = I0, upper = X)
    if(I0<I1){aux = integral$value}
    else{aux = 0}
    return(aux)
  }
  Y_1 = function(Y){
    integrand_y1 = function(r){ #Creating a function that returns a vectorized integrand
      aux = vector(length = length(r))
      for(i in 1:length(r)){
        aux[i] = (y_1(w0,w1,r[i],Y))^((sigma-1)/sigma) #integrand of y_1 in the consumption good production
      }
      return(aux)
    }
    integral = integrate(integrand_y1, lower = max(I1,X), upper = N)
    aux = integral$value
    return(aux)
  }
  aux = Y - (Y_k(Y)+Y_0(Y)+Y_1(Y))^(sigma/(sigma-1))
  return(aux) 
}


# Testing Rcpp ------------------------------------------------------------

Rcpp::sourceCpp('~/Dropbox/Technology/Codes/Technology_Health/lib/l0_s_cpp.cpp')

trial_seq = seq(from=0, to=1000, by=1)
microbenchmark(times=1000, unit="ms", l0_s_vector_R(trial_seq), l0_s_cpp(trial_seq))

#So just changing the function to Rcpp becomes 4 times faster.
#Now testing vectoriazed function u0
Rcpp::sourceCpp('~/Dropbox/Technology/Codes/Technology_Health/lib/u0_cpp.cpp')


microbenchmark(times=1000, unit="ms", u0(theta = 0.5, h='g', w0 = 2.4, m = trial_seq),
               u0_cpp(theta = 0.5, w0 = 2.4, m = trial_seq))
#In this case is 8 times faster!


