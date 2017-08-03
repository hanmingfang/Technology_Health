#Solve for static equilibrium 

# Libraries ---------------------------------------------------------------

#library(stargazer, quietly = TRUE)
#library(lmtest, quietly = TRUE)
library(ggplot2, quietly = TRUE)
#library(dplyr)
#library(reshape2)
#library(ggfortify)
#library(data.table)
#library(plotly)
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
#Change directory
dir = '~/Dropbox/Technology/Codes/Technology_Health/'
setwd(dir)

# Parameters --------------------------------------------------------------
#Household 
xi = 1              #Utility function
phi = 1             #Cost of labor effort  
lambda_g = 0.5      #Measure healthy workers
theta_L = 0         #Domain for theta
theta_H = Inf
m_L = 0             #Domain for medical expenditure shocks
m_F = Inf
mu_g = 2            #Mean of theta for healthy workers (Truncated Normal)
mu_b = 0            #Mean of theta for healthy workers (Truncated Normal)
sd_g = 1            #Standard deviation of theta for healthy workers (Truncated Normal)
sd_b = 2            #Standard deviation of theta for unhealthy workers (Truncated Normal)
shape_g = 2         #Shape parameter of theta (Gamma distribution) 
shape_b = 0.2       #Shape parameter of theta (Gamma distribution) 
scale_g = 1         #Scale parameter of theta (Gamma distribution) 
scale_b = 1         #Scale parameter of theta (Gamma distribution) 
rate_g = 1          #Rate for exponential distribution for medical exp. healthy workers  (Mean=1/rate)
rate_b = 0.25       #Rate for exponential distribution for medical exp. unhealthy workers
util_min = 0.0001   #Minimum consumption minus labor effort a household can have (can't be 0 or blows up)
#Check papaers of Erick French and Christina Dinadi (Medicaid)
#This is not in the same way in the document!
P_0g =  0.3         #Probability of 0 medical expenditure for healthy worker 
P_0b =  0.1         #Probability of 0 medical expenditure for unhealthy worker
theta_ins_final = 10 #As we can not evaluate f(Inf) I use an upper bound number, but allowing uniroot to extend it in theta_ins
#Firm
N = 1               #Range of tasks (upper limit)
eta = 0.5           #Distribution parameter of the CES
rho = 0.8           #Relative labor productivity of unhealthy workers
psi = 1             #Price of intermediates
sigma = 2           #Elasticity of substitution between tasks
zeta = 2            #Elasticity of substitution between factors (if fixed), just to define zeta_elas
#Change this to a small positive number
C_IN = 0.0          #Health Insurance Fixed Cost (we can start with a very low one)
A = 1               #Parameter in labor productivity
A_0 = 1             #Parameter in labor productivity
lambda_d = 10       #Parameter in sorting function
alpha_d = 5         #Parameter in sorting function
D = 1               #Parameter in Automation Cost function
tol = 1e-8          #Tolerance for unitroot, affects computation time
K = 3               #Capital stock in the economy
# Primitive Functions ---------------------------------------------------------------
#Distribution objects
#Distribution for Positive part of Medical expenditure
D_mg = Exp(rate_g)
D_mb = Exp(rate_b)
#Distribution for Risk aversion parameter
D_theta_g = Gammad(scale = scale_g,shape = shape_g)
D_theta_b = Gammad(scale = scale_b,shape = shape_b)
#Conditional cdf of Households' risk aversion parameter
F_g  = function(theta){
  #aux =  ptruncnorm(theta, a = theta_L, b = theta_H, mean = mu_g, sd = sd_g) #Good health
  aux = p(D_theta_g)(theta) #Expectation is shape*rate
  return(aux)
} 
F_b  = function(theta){
  #aux =  ptruncnorm(theta, a = theta_L, b = theta_H, mean = mu_b, sd = sd_b) #Bad health
  aux = p(D_theta_b)(theta)
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
#Labor productivity 
gamma_prod = function(i){   
  aux = A_0*exp(A*i)
  return(aux)
}   
#Sorting of workers 
#TODO: edit this function
#Making that integrates to 1
norm_const_delta_sort = integrate(Vectorize(function(i) exp(lambda_d*i - alpha_d)/(1+exp(lambda_d*i - alpha_d))),
                                  lower = N-1, upper = N)$value 
delta_sort = function(i){
  aux = exp(lambda_d*i - alpha_d)/(1+exp(lambda_d*i - alpha_d))/norm_const_delta_sort
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
  aux = 1
  return(aux)
}
#Individual labor suply for insurance
l1_s = function(w1){
  aux = 1
  return(aux)
}
#Utility function under  insurance
u1 = function(theta,h,w1){
  l1 = l1_s(w1)
  aux = -exp(-theta*(w1*l1))
  return(aux)
}
#Utility function under No insurance
u0 = function(theta,h,w0,m){
  l0 = l0_s(w0)
  aux =  -exp(-theta*(w0*l0 - m))
  return(aux)
}
#Expected utility under no insurance (simulated version)
E_u0 = function(theta,h,w0){
  #require("distrEx")
  #Assigning distribution
  if(h == 'g'){ #If healthy worker
    aux = P_0g*u0(theta,h,w0,m=0) + (1-P_0g)*E(D_mg, u0, theta = theta, h = h, w0 = w0)
  }
  else{
    aux = P_0b*u0(theta,h,w0,m=0) + (1-P_0b)*E(D_mb, u0, theta = theta, h = h, w0 = w0)
  }
  return(aux)
}
#Threshold for household insurance decision
#TODO: include \n in the message
theta_ins = function(h,w0,w1){
  initial = theta_L  #Search over
  #final = theta_H
  #As we can not evaluate f(Inf) I use a lower number, but allowing uniroot to extend it
  final = theta_ins_final 
  fun = function (theta) E_u0(theta,h,w0) - u1(theta,h,w1)   #This is a decreasing function of theta
  if(E_u0(theta_L,h,w0) - u1(theta_L,h,w1) < 0){
    aux = theta_L
    } #If at lower bound is negative
  #TODO: Uncomment the next line if we have bouded support for Theta
  #else if(E_u0(theta_H,h,w0) - u1(theta_H,h,w1) > 0){aux = theta_H} #This line works only with bounded support for theta
  else{
    aux = tryCatch(
      {
       return(uniroot(fun, c(initial,final), tol = tol, extendInt = "downX")$root) #Get the root, 
        #downX is to tell that is decresing on theta, so can look further than the specified range
      },
      error = function(cond){
        message(paste(cond,". Taking theta_ins(",round(w0,2),",",round(w1,2),") = +Inf"))
        return(Inf)
      }
    )
  }
  return(aux) 
}
#Aggregate labor supply for no insurance
L0_s = function(h,w0,w1){
  if(h == 'g'){aux = lambda_g*l0_s(w0)*F_g(theta_ins(h,w0,w1))} # L^0_g
  else{aux = (1-lambda_g)*l0_s(w0)*F_b(theta_ins(h,w0,w1))} # L^0_b
  return(aux)
}
#Aggregate labor supply for no insurance with memoise
L0_s_memo = memoise(L0_s)
#Aggregate labor supply for insurance
L1_s = function(h,w0,w1){
  if(h == 'g'){aux = lambda_g*l1_s(w1)*(1-F_g(theta_ins(h,w0,w1)))} # L^1_g
  else{aux = (1-lambda_g)*l1_s(w1)*(1-F_b(theta_ins(h,w0,w1)))} # L^1_b
  return(aux)
}
#Aggregate labor supply for insurance with memoise
L1_s_memo = memoise(L1_s) 
#Endogenous proportion of healthy workers for no health insurance
#TODO: Change beliefs out of path when appropiate
#Same delta_sort() for the moment
Chi_0g = function(w0,w1,i){
  #Assign the same beliefs for the moment if there is no Labor supply
  if(is.na(L0_s('g',w0,w1)/(L0_s('g',w0,w1) + L0_s('b',w0,w1)))){
    aux = Chi_1g(w0,w1,i)
  }
  else{
    aux = (delta_sort(i)*L0_s('g',w0,w1)/(delta_sort(i)*L0_s('g',w0,w1) + L0_s('b',w0,w1))) 
  }
  return(aux)
}
#Endogenous proportion of healthy workers for health insurance
#TODO: edit this function
Chi_1g = function(w0,w1,i){
  #Assign the same beliefs for the moment if there is no Labor supply
  if(is.na(L1_s('g',w0,w1)/(L1_s('g',w0,w1) + L1_s('b',w0,w1)))){
    aux = Chi_0g(w0,w1,i)
  }
  else{
    aux = (delta_sort(i)*L1_s('g',w0,w1)/(delta_sort(i)*L1_s('g',w0,w1) + L1_s('b',w0,w1))) 
  }
  return(aux)
}
#Assign beliefs out of path
#Expected expenditure shock
E_m = function(h){
  if(h == 'g'){aux = (1-P_0g)*E(D_mg)}
  else{aux = (1-P_0b)*E(D_mb)}
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
#TODO: Check monotonicity of the function
#Check monotonicity of the threshold itself (on wages)
#There is a discontinuity too
#TODO: Adapt the function for Adverse selection too
X_tilde = function(w0,w1,Y){
  LHS = function(i) C_IN/(((B(i)*(sigma-1)/sigma)^(sigma-1))*Y/sigma)
  #Do not call the same function more than one time if is not neccesary
  L0_s_g_var = L0_s_memo('g',w0,w1)
  L0_s_b_var = L0_s_memo('b',w0,w1)
  L1_s_g_var = L1_s_memo('g',w0,w1)
  L1_s_b_var = L1_s_memo('b',w0,w1)
  #Assigning the same beliefs if no labor supply
  if(is.na((L0_s_g_var)/(L0_s_g_var + L0_s_b_var))){
    Chi_0gi = function(i) (delta_sort(i)*L1_s_g_var)/(delta_sort(i)*L1_s_g_var + L1_s_b_var)
    Chi_1gi = function(i) ((delta_sort(i)*L1_s_g_var)/(delta_sort(i)*L1_s_g_var + L1_s_b_var))
  }
  else if(is.na((L1_s_g_var)/(L1_s_g_var + L1_s_b_var))){
    Chi_0gi = function(i) (delta_sort(i)*L0_s_g_var)/(delta_sort(i)*L0_s_g_var + L0_s_b_var)
    Chi_1gi = function(i) ((delta_sort(i)*L0_s_g_var)/(delta_sort(i)*L0_s_g_var + L0_s_b_var))
  }
  else{
    Chi_0gi = function(i) (delta_sort(i)*L0_s_g_var)/(delta_sort(i)*L0_s_g_var + L0_s_b_var)
    Chi_1gi = function(i) ((delta_sort(i)*L1_s_g_var)/(delta_sort(i)*L1_s_g_var + L1_s_b_var))
  }
  gamma_prod_bar_0i = function(i) gamma_prod(i)*((1-rho)*Chi_0gi(i)+rho)
  w_hat0i = function(i) w0/gamma_prod_bar_0i(i)
  p_0i = function(i) (eta*w_hat0i(i))/((1-eta)*psi)
  ###
  gamma_prod_bar_1i = function(i) gamma_prod(i)*((1-rho)*Chi_1gi(i)+rho)
  E_mg = E_m('g')
  E_mb = E_m('b')
  Mi = function(i) (E_mg*Chi_1gi(i)+E_mb*(1-Chi_1gi(i)))/l1_s(w1)
  w_hat1i = function(i) (w1+Mi(i))/gamma_prod_bar_1i(i)
  p_1i = function(i) (eta*w_hat1i(i))/((1-eta)*psi)
  RHS1 = function(i) (((eta*p_1i(i)^(zeta_elas(i)-1)+(1-eta))^(zeta_elas(i)/(zeta_elas(i)-1)))/(w_hat1i(i)+psi*p_1i(i)^zeta_elas(i)))^(sigma-1)
  RHS0 = function(i) (((eta*p_0i(i)^(zeta_elas(i)-1)+(1-eta))^(zeta_elas(i)/(zeta_elas(i)-1)))/(w_hat0i(i)+psi*p_0i(i)^zeta_elas(i)))^(sigma-1)
  
  fun = function (i) LHS(i)-RHS1(i)+RHS0(i)   #This is NOT a monotonic function!. Although if C_IN =0 and we have advantageous selection, 
  #then if the function crosses the X axis, then we know it will cross downward.
  initial = N-1 #Lower bound for i
  final = N #Upper bound for i
  #Boundary cases
  #if(is.na(w_hat0i(N)) & !is.na(w_hat1i(N))){aux = N-1} #If w_hat0 is NA is because every household works for contract
  #with insurance, thus X_tilde has to be N-1
  #else if(!is.na(w_hat0i(N)) & is.na(w_hat1i(N))){aux = N}#The opposite than before
  if(fun(N)>0){aux = N} #This should be correct given that the function is DECREASING?
  else if(fun(N-1)<0){aux = N-1}
  else{aux = uniroot(fun, c(initial,final), tol = tol)$root}
  return(aux) 
}
#Indifference threshold between insurance and not insurance with memoise
X_tilde_memo = memoise(X_tilde)
#Indifference threshold between capital and not insurance
I_tilde0 = function(w0,w1,R,Y){
  LHS = function(i) C_A(i)/(((B(i)*(sigma-1)/sigma)^(sigma-1))*Y/sigma)
  #Do not call the same function more than one time if is not neccesary
  L0_s_g_var = L0_s_memo('g',w0,w1)
  L0_s_b_var = L0_s_memo('b',w0,w1)
  L1_s_g_var = L1_s_memo('g',w0,w1)
  L1_s_b_var = L1_s_memo('b',w0,w1)
  #Assigning the same beliefs if no labor supply
  if(is.na((L0_s_g_var)/(L0_s_g_var + L0_s_b_var))){
    Chi_0gi = function(i) (delta_sort(i)*L1_s_g_var)/(delta_sort(i)*L1_s_g_var + L1_s_b_var)
    Chi_1gi = function(i) ((delta_sort(i)*L1_s_g_var)/(delta_sort(i)*L1_s_g_var + L1_s_b_var))
  }
  else if(is.na((L1_s_g_var)/(L1_s_g_var + L1_s_b_var))){
    Chi_0gi = function(i) (delta_sort(i)*L0_s_g_var)/(delta_sort(i)*L0_s_g_var + L0_s_b_var)
    Chi_1gi = function(i) ((delta_sort(i)*L0_s_g_var)/(delta_sort(i)*L0_s_g_var + L0_s_b_var))
  }
  else{
    Chi_0gi = function(i) (delta_sort(i)*L0_s_g_var)/(delta_sort(i)*L0_s_g_var + L0_s_b_var)
    Chi_1gi = function(i) ((delta_sort(i)*L1_s_g_var)/(delta_sort(i)*L1_s_g_var + L1_s_b_var))
  }
  gamma_prod_bar_0i = function(i) gamma_prod(i)*((1-rho)*Chi_0gi(i)+rho)
  w_hat0i = function(i) w0/gamma_prod_bar_0i(i)
  p_0i = function(i) (eta*w_hat0i(i))/((1-eta)*psi)
  ###
  R_hati = function(i) R/z_prod(i)
  p_ki = function(i) (eta*R_hati(i))/((1-eta)*psi)
  RHSk = function(i) (((eta*p_ki(i)^(zeta_elas(i)-1)+(1-eta))^(zeta_elas(i)/(zeta_elas(i)-1)))/(R_hati(i)+psi*p_ki(i)^zeta_elas(i)))^(sigma-1)
  RHS0 = function(i) (((eta*p_0i(i)^(zeta_elas(i)-1)+(1-eta))^(zeta_elas(i)/(zeta_elas(i)-1)))/(w_hat0i(i)+psi*p_0i(i)^zeta_elas(i)))^(sigma-1)
  
  fun = function (i) LHS(i)-RHSk(i)+RHS0(i) # Under the assumption that profits conditional on capital are non increasing on i
  #this function should be INCREASING (regardless of advantageous selection, due to the increasing sorting of delta_sort)
  initial = N-1 #Lower bound for i
  final = N #Upper boud for i
  #boundary cases
  if(fun(N)<0){aux = N} #This should be correct given that the function is INCREASING in this case
  else if(fun(N-1)>0){aux = N-1}
  else{aux = uniroot(fun, c(initial,final), tol = tol)$root}
  return(aux) 
}
#Indifference threshold between capital and not insurance with memoise
I_tilde0_memo = memoise(I_tilde0)
#Indifference threshold between capital and insurance
I_tilde1 = function(w0,w1,R,Y){
  LHS = function(i) (C_IN-C_A(i))/(((B(i)*(sigma-1)/sigma)^(sigma-1))*Y/sigma)
  ###
  L0_s_g_var = L0_s_memo('g',w0,w1)
  L0_s_b_var = L0_s_memo('b',w0,w1)
  L1_s_g_var = L1_s_memo('g',w0,w1)
  L1_s_b_var = L1_s_memo('b',w0,w1)
  #Assigning the same beliefs if no labor supply
  if(is.na((L0_s_g_var)/(L0_s_g_var + L0_s_b_var))){
    Chi_0gi = function(i) (delta_sort(i)*L1_s_g_var)/(delta_sort(i)*L1_s_g_var + L1_s_b_var)
    Chi_1gi = function(i) ((delta_sort(i)*L1_s_g_var)/(delta_sort(i)*L1_s_g_var + L1_s_b_var))
  }
  else if(is.na((L1_s_g_var)/(L1_s_g_var + L1_s_b_var))){
    Chi_0gi = function(i) (delta_sort(i)*L0_s_g_var)/(delta_sort(i)*L0_s_g_var + L0_s_b_var)
    Chi_1gi = function(i) ((delta_sort(i)*L0_s_g_var)/(delta_sort(i)*L0_s_g_var + L0_s_b_var))
  }
  else{
    Chi_0gi = function(i) (delta_sort(i)*L0_s_g_var)/(delta_sort(i)*L0_s_g_var + L0_s_b_var)
    Chi_1gi = function(i) ((delta_sort(i)*L1_s_g_var)/(delta_sort(i)*L1_s_g_var + L1_s_b_var))
  }
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
  if(fun(N)>0){aux = N} #This should be correct given that the function is DECREASING in this case
  else if(fun(N-1)<0){aux = N-1}
  else{aux = uniroot(fun, c(initial,final), tol = tol)$root}
  return(aux) 
}
#Indifference threshold between capital and insurance with memoise
I_tilde1_memo = memoise(I_tilde1) 
#Excess demand for Capital Fast and Vectorized (the fastest)
k_excess_d_fast_vec = function(p){
  w0 = p[1]
  w1 = p[2]
  R = p[3]
  Y = p[4]
  I0 = I_tilde0_memo(w0,w1,R,Y)
  I1 = I_tilde1_memo(w0,w1,R,Y)
  ###
  R_hati = function(i) R/z_prod(i)
  p_ki = function(i) (eta*R_hati(i))/((1-eta)*psi)
  ###
  yki = function(i) {
    aux = (Y*((sigma-1)/sigma)^sigma)*
      ((B(i)*(eta*p_ki(i)^(zeta_elas(i)-1)+(1-eta))^
          (zeta_elas(i)/(zeta_elas(i)-1)))/(R_hati(i) + psi*p_ki(i)^zeta_elas(i)))^sigma
    return(aux)
  }
  integrand_k = function(i) yki(i)/(B(i)*(eta*p_ki(i)^(zeta_elas(i)-1)+(1-eta))^(zeta_elas(i)/(zeta_elas(i)-1)))
  integral = integrate(Vectorize(integrand_k), lower = N-1, upper = min(I0,I1)) #Bounds specified in the document
  #We need to vectorize the function to use integrate
  aux = integral$value-K
  return(aux)
}
#Excess demand for labor without insurance fast and vectorized
#TODO: CHECK THIS FOR NON INTERIOR CASES
l0_excess_d_fast_vec = function(p){
  w0 = p[1]
  w1 = p[2]
  R = p[3]
  Y = p[4]
  X = X_tilde_memo(w0,w1,Y)
  I0 = I_tilde0_memo(w0,w1,R,Y)
  I1 = I_tilde1_memo(w0,w1,R,Y)
  L0_s_g_var = L0_s_memo('g',w0,w1)
  L0_s_b_var = L0_s_memo('b',w0,w1)
  L1_s_g_var = L1_s_memo('g',w0,w1)
  L1_s_b_var = L1_s_memo('b',w0,w1)
  #Assigning the same beliefs if no labor supply
  if(is.na((L0_s_g_var)/(L0_s_g_var + L0_s_b_var))){
    Chi_0gi = function(i) (delta_sort(i)*L1_s_g_var)/(delta_sort(i)*L1_s_g_var + L1_s_b_var)
    Chi_1gi = function(i) ((delta_sort(i)*L1_s_g_var)/(delta_sort(i)*L1_s_g_var + L1_s_b_var))
  }
  else if(is.na((L1_s_g_var)/(L1_s_g_var + L1_s_b_var))){
    Chi_0gi = function(i) (delta_sort(i)*L0_s_g_var)/(delta_sort(i)*L0_s_g_var + L0_s_b_var)
    Chi_1gi = function(i) ((delta_sort(i)*L0_s_g_var)/(delta_sort(i)*L0_s_g_var + L0_s_b_var))
  }
  else{
    Chi_0gi = function(i) (delta_sort(i)*L0_s_g_var)/(delta_sort(i)*L0_s_g_var + L0_s_b_var)
    Chi_1gi = function(i) ((delta_sort(i)*L1_s_g_var)/(delta_sort(i)*L1_s_g_var + L1_s_b_var))
  }
  gamma_prod_bar_0i = function(i) gamma_prod(i)*((1-rho)*Chi_0gi(i)+rho)
  w_hat0i = function(i) w0/gamma_prod_bar_0i(i)
  p_0i = function(i) (eta*w_hat0i(i))/((1-eta)*psi)
  if(I0<I1 & I0<X){ #Checking for interior cases
    y0i = function(i){
      aux = (Y*((sigma-1)/sigma)^sigma)*
        ((B(i)*(eta*p_0i(i)^(zeta_elas(i)-1)+(1-eta))^
            (zeta_elas(i)/(zeta_elas(i)-1)))/(w_hat0i(i) + psi*p_0i(i)^zeta_elas(i)))^sigma
      return(aux)
    }  
    #The integrand of l_0i should be l_hat_0i divided by gamma gamma_prod_bar_0i
    integrand_l0 = function(i) (y0i(i)/(B(i)*(eta*p_0i(i)^(zeta_elas(i)-1)+(1-eta))^(zeta_elas(i)/(zeta_elas(i)-1))))/gamma_prod_bar_0i(i)
    integral = integrate(Vectorize(integrand_l0), lower = I0, upper = X) #Bounds specified in the document, Vectorizing
    aux = integral$value - (L0_s_g_var + L0_s_b_var) #subtracting Labor supplied for no insurance
  }
  else{aux = 0 - (L0_s_g_var + L0_s_b_var)}#If there is no labor without insurance in equilibrium
  return(aux)  
}
#Excess demand fo labor with insurance
#TODO: CHECK THIS FOR NON INTERIOR CASES
l1_excess_d_fast_vec = function(p){
  w0 = p[1]
  w1 = p[2]
  R = p[3]
  Y = p[4]
  X = X_tilde_memo(w0,w1,Y)
  I0 = I_tilde0_memo(w0,w1,R,Y)
  I1 = I_tilde1_memo(w0,w1,R,Y)
  ###
  L0_s_g_var = L0_s_memo('g',w0,w1)
  L0_s_b_var = L0_s_memo('b',w0,w1)
  L1_s_g_var = L1_s_memo('g',w0,w1)
  L1_s_b_var = L1_s_memo('b',w0,w1)
  #Assigning the same beliefs if no labor supply
  if(is.na((L0_s_g_var)/(L0_s_g_var + L0_s_b_var))){
    Chi_0gi = function(i) (delta_sort(i)*L1_s_g_var)/(delta_sort(i)*L1_s_g_var + L1_s_b_var)
    Chi_1gi = function(i) ((delta_sort(i)*L1_s_g_var)/(delta_sort(i)*L1_s_g_var + L1_s_b_var))
  }
  else if(is.na((L1_s_g_var)/(L1_s_g_var + L1_s_b_var))){
    Chi_0gi = function(i) (delta_sort(i)*L0_s_g_var)/(delta_sort(i)*L0_s_g_var + L0_s_b_var)
    Chi_1gi = function(i) ((delta_sort(i)*L0_s_g_var)/(delta_sort(i)*L0_s_g_var + L0_s_b_var))
  }
  else{
    Chi_0gi = function(i) (delta_sort(i)*L0_s_g_var)/(delta_sort(i)*L0_s_g_var + L0_s_b_var)
    Chi_1gi = function(i) ((delta_sort(i)*L1_s_g_var)/(delta_sort(i)*L1_s_g_var + L1_s_b_var))
  }
  gamma_prod_bar_1i = function(i) gamma_prod(i)*((1-rho)*Chi_1gi(i)+rho)
  E_mg = E_m('g')
  E_mb = E_m('b')
  Mi = function(i) (E_mg*Chi_1gi(i)+E_mb*(1-Chi_1gi(i)))/l1_s(w1)
  w_hat1i = function(i) (w1+Mi(i))/gamma_prod_bar_1i(i)
  p_1i = function(i) (eta*w_hat1i(i))/((1-eta)*psi)
  if(max(I1,X)<N){
    y1i = function(i){
      aux = (Y*((sigma-1)/sigma)^sigma)*
        ((B(i)*(eta*p_1i(i)^(zeta_elas(i)-1)+(1-eta))^
            (zeta_elas(i)/(zeta_elas(i)-1)))/(w_hat1i(i) + psi*p_1i(i)^zeta_elas(i)))^sigma
      return(aux)
    }
    integrand_l1 = function(i) (y1i(i)/(B(i)*(eta*p_1i(i)^(zeta_elas(i)-1)+(1-eta))^(zeta_elas(i)/(zeta_elas(i)-1))))/gamma_prod_bar_1i(i)
    integral = integrate(Vectorize(integrand_l1), lower = max(I1,X), upper = N) #Bounds specified in the document, vectorize
    aux = integral$value - (L1_s_g_var + L1_s_b_var) #subtracting Labor supplied for insurance
  }
  else{aux = 0 - (L1_s_g_var + L1_s_b_var)}
  return(aux)  
}
#Total consumption good
#TODO: Include indicator function, boudarie cases and Review
Y_excess_s_fast_vec = function(p){
  w0 = p[1]
  w1 = p[2]
  R = p[3]
  Y = p[4]
  X = X_tilde_memo(w0,w1,Y)
  I0 = I_tilde0_memo(w0,w1,R,Y)
  I1 = I_tilde1_memo(w0,w1,R,Y)
  L0_s_g_var = L0_s_memo('g',w0,w1)
  L0_s_b_var = L0_s_memo('b',w0,w1)
  L1_s_g_var = L1_s_memo('g',w0,w1)
  L1_s_b_var = L1_s_memo('b',w0,w1)
  #Assigning the same beliefs if no labor supply
  if(is.na((L0_s_g_var)/(L0_s_g_var + L0_s_b_var))){
    Chi_0gi = function(i) (delta_sort(i)*L1_s_g_var)/(delta_sort(i)*L1_s_g_var + L1_s_b_var)
    Chi_1gi = function(i) ((delta_sort(i)*L1_s_g_var)/(delta_sort(i)*L1_s_g_var + L1_s_b_var))
  }
  else if(is.na((L1_s_g_var)/(L1_s_g_var + L1_s_b_var))){
    Chi_0gi = function(i) (delta_sort(i)*L0_s_g_var)/(delta_sort(i)*L0_s_g_var + L0_s_b_var)
    Chi_1gi = function(i) ((delta_sort(i)*L0_s_g_var)/(delta_sort(i)*L0_s_g_var + L0_s_b_var))
  }
  else{
    Chi_0gi = function(i) (delta_sort(i)*L0_s_g_var)/(delta_sort(i)*L0_s_g_var + L0_s_b_var)
    Chi_1gi = function(i) ((delta_sort(i)*L1_s_g_var)/(delta_sort(i)*L1_s_g_var + L1_s_b_var))
  }
  Y_k = function(Y){
    R_hati = function(i) R/z_prod(i)
    p_ki = function(i) (eta*R_hati(i))/((1-eta)*psi)
    yki = function(i) {
        aux = (Y*((sigma-1)/sigma)^sigma)*
          ((B(i)*(eta*p_ki(i)^(zeta_elas(i)-1)+(1-eta))^
              (zeta_elas(i)/(zeta_elas(i)-1)))/(R_hati(i) + psi*p_ki(i)^zeta_elas(i)))^sigma
        return(aux)
      }
    integrand_yk = function(i) (yki(i))^((sigma-1)/sigma) #integrand of y_k in the consumption good production
    integral = integrate(Vectorize(integrand_yk), lower = N-1, upper = min(I1,I0)) #Vectorizing
    aux = integral$value
    return(aux)
  }
  Y_0 = function(Y){
    gamma_prod_bar_0i = function(i) gamma_prod(i)*((1-rho)*Chi_0gi(i)+rho)
    w_hat0i = function(i) w0/gamma_prod_bar_0i(i)
    p_0i = function(i) (eta*w_hat0i(i))/((1-eta)*psi)
    y0i = function(i){
      aux = (Y*((sigma-1)/sigma)^sigma)*
        ((B(i)*(eta*p_0i(i)^(zeta_elas(i)-1)+(1-eta))^
            (zeta_elas(i)/(zeta_elas(i)-1)))/(w_hat0i(i) + psi*p_0i(i)^zeta_elas(i)))^sigma
      return(aux)
    }  
    integrand_y0 = function(i) (y0i(i))^((sigma-1)/sigma) #integrand of y_0 in the consumption good production
    if(I0<I1 & I0<X){
      integral = integrate(Vectorize(integrand_y0), lower = I0, upper = X)
      aux = integral$value
    }
    else{aux = 0}
    return(aux)
  }
  Y_1 = function(Y){
    gamma_prod_bar_1i = function(i) gamma_prod(i)*((1-rho)*Chi_1gi(i)+rho)
    E_mg = E_m('g')
    E_mb = E_m('b')
    Mi = function(i) (E_mg*Chi_1gi(i)+E_mb*(1-Chi_1gi(i)))/l1_s(w1)
    w_hat1i = function(i) (w1+Mi(i))/gamma_prod_bar_1i(i)
    p_1i = function(i) (eta*w_hat1i(i))/((1-eta)*psi)
    y1i = function(i){
      aux = (Y*((sigma-1)/sigma)^sigma)*
        ((B(i)*(eta*p_1i(i)^(zeta_elas(i)-1)+(1-eta))^
            (zeta_elas(i)/(zeta_elas(i)-1)))/(w_hat1i(i) + psi*p_1i(i)^zeta_elas(i)))^sigma
      return(aux)
    }
    integrand_y1 = function(i) (y1i(i))^((sigma-1)/sigma) #integrand of y_1 in the consumption good production
    if(max(I1,X)<N){
      integral = integrate(Vectorize(integrand_y1), lower = max(I1,X), upper = N) #Vectorizing
      aux = integral$value
    }
    else{aux = 0}
    return(aux)
  }
  aux = Y - (Y_k(Y)+Y_0(Y)+Y_1(Y))^(sigma/(sigma-1))
  return(aux) 
}
#Objective function to minimize
obj_fun = function(p){  
  Y = p[4]
  aux  = sqrt((k_excess_d_fast_vec(p)/K)^2 + (l0_excess_d_fast_vec(p))^2 +
                (l1_excess_d_fast_vec(p))^2 + (Y_excess_s_fast_vec(p)/Y)^2)
  return(aux)#Here I need to normalize in some way the Excess demands
}

#TODO: Solve the corner cases for the thresholds and excess demands
#Fix the inconsistency with the beliefs out of path (Assign some beliefs there)

# Optim -------------------------------------------------------------------
# Start the clock!
ptm = proc.time()
#Optimize and get the parameters
#Is stopping  "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
#This one works if tol= 1e-8 for the unitroot
optim_object = optim(par=c(2,1,1,30), fn = obj_fun, method = "L-BFGS-B", 
                  lower = c(0.001,0.001,0.001,0.001),  
                  control = list(trace = 1, pgtol = 1e-12),
                  upper = c(Inf,Inf,Inf,Inf))
# Stop the clock
proc.time() - ptm
optim_object
#Checking results
p = optim_object$par
w0 = p[1]
w1 = p[2]
R = p[3]
Y = p[4]

# nloptr ------------------------------------------------------------------
#The next ones are quite far from the result
# Start the clock!
ptm = proc.time()
#Optimize and get the parameters
nloptim_object = nloptr(x0=c(2,1,1,10), eval_f = obj_fun,
                      lb = c(0.001,0.001,0.001,0.001), 
                      opts= list(algorithm ="NLOPT_GN_ORIG_DIRECT", maxeval = 10000, xtol_abs = 1.0e-10),
                      ub = c(10,10,10,100))
print(nloptim_object)  
# Stop the clock
proc.time() - ptm
p = nloptim_object$solution
w0 = p[1]
w1 = p[2]
R = p[3]
Y = p[4]
#Local search with global solution as initial guess
# Start the clock!
ptm = proc.time()
#Optimize and get the parameters
nloptim_object = nloptr(x0= nloptim_object$solution, eval_f = obj_fun,
                        lb = c(0.001,0.001,0.001,0.001), 
                        opts= list(algorithm ="NLOPT_LN_COBYLA",xtol_absl =1.0e-10, maxeval = 10000),
                        ub = c(10,10,10,100))
print(nloptim_object)  
# Stop the clock
proc.time() - ptm

#Global optimizer with CRS2 (is the best working now)
ptm = proc.time()
crs2_sol = crs2lm(x0=c(2,1,1,10), fn = obj_fun,
             lower = c(0.001,0.001,0.001,0.001),
             upper = c(50,50,50,500),
             maxeval = 10000,
             xtol_rel = 1e-6)
proc.time() - ptm
crs2_sol
p = crs2_sol$par
w0 = p[1]
w1 = p[2]
R = p[3]
Y = p[4]


# Multiroot ---------------------------------------------------------------
#Trying to use Multiroot
#Works, but finds the corner
model = function(p) c(F1 = k_excess_d_fast_vec(c(exp(p[1]),exp(p[2]),exp(p[3]),exp(p[4]))),
                       F2 = l0_excess_d_fast_vec(c(exp(p[1]),exp(p[2]),exp(p[3]),exp(p[4]))),
                       F3 = l1_excess_d_fast_vec(c(exp(p[1]),exp(p[2]),exp(p[3]),exp(p[4]))),
                       F4 = Y_excess_s_fast_vec(c(exp(p[1]),exp(p[2]),exp(p[3]),exp(p[4]))))

#Using exponential to ensure strictly positive values
ss = multiroot(model, c(2,1,1,20), useFortran = TRUE)

#Changing to get a different root
ss = multiroot(model, c(3,2,1.7,30), useFortran = TRUE)

p = exp(ss$root)#remeber to take the exponential after
w0 = p[1]
w1 = p[2]
R = p[3]
Y = p[4]

# Plots -------------------------------------------------------------------
dir = '~/Dropbox/Technology/Codes/Technology_Health/plots/'
setwd(dir)
#Wages to play
w0 = 2.5
w1 = 1
R = 1.7
#Need to specify also total output to play
Y = 10
#New plots
#Check for Proportions of Helathy workers
ggplot(data.frame(x=c(N-1,N)), aes(x=x)) + 
  stat_function(fun=Chi_0gi, geom="line", aes(colour = "Chi_0g")) + xlab("i") + 
  ylab("") + stat_function(fun=Chi_1gi, geom="line",aes(colour = "Chi_1g"))
#Check crossing of effective wages 
ggplot(data.frame(x=c(N-1,N)), aes(x=x)) + 
  stat_function(fun=w_hat0i, geom="line", aes(colour = "what0")) + xlab("i") + 
  ylab("") + stat_function(fun=w_hat1i, geom="line",aes(colour = "what1"))
#Check crossing of function for thresholds
ggplot(data.frame(x=c(N-1,N)), aes(x=x)) + 
  stat_function(fun=fun, geom="line", aes(colour = "Function for threshold")) + xlab("i") + 
  ylab("")
#Check for expected Medical expenditure
ggplot(data.frame(x=c(N-1,N)), aes(x=x)) + 
  stat_function(fun=Mi, geom="line", aes(colour = "Mi")) + xlab("i") + 
  ylab("")
#Graphing Advantageous selection
#Is very non-monotonic, like oscilating
Adv_sel = function(w0) F_b(theta_ins('b',w0,1)) - F_g(theta_ins('g',w0,1))
ggplot(data.frame(x=c(0,10)), aes(x=x)) + 
  stat_function(fun = Vectorize(Adv_sel), geom="line", aes(colour = "Advan Sel")) + xlab("w0") + ylab("")
#Graphing X_tilde as a function of w0
#Good! It seems that X_tilde is continuous and monotonic with respect to changes in wages 
X_tilde_plot = function(w0) X_tilde(w0,1,10)
ggplot(data.frame(x=c(2.75,4.75)), aes(x=x)) + 
  stat_function(fun = Vectorize(X_tilde_plot), geom="line", aes(colour = "X_tilde")) + xlab("w0") + ylab("")
###
#Plot to show that FOSD in Assumption 1 holds for this case
ggplot(data.frame(x=c(0, 15)), aes(x=x)) + 
  stat_function(fun=H_g, geom="line", aes(colour = "H_g")) + xlab("x") + 
  ylab("y") + stat_function(fun=H_b, geom="line",aes(colour = "H_b")) 
#Plot to show that FOSD in Assumption 2 holds for this case
ggplot(data.frame(x=c(0, 4)), aes(x=x)) + 
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
#Plot Chi_0g and Chi_1g (change wages to get advantageous selection) 
Chi_0g_plot = function(i) Chi_0g(w0=w0, w1=w1,i)
Chi_1g_plot = function(i) Chi_1g(w0=w0, w1=w1,i)
ggplot(data.frame(x=c(N-1,N)), aes(x=x)) + 
  stat_function(fun=Chi_0g_plot, geom="line", aes(colour = "Chi_0g")) + xlab("i") + 
  ylab("") + stat_function(fun=Chi_1g_plot, geom="line",aes(colour = "Chi_1g"))
ggsave(file="endogenous_proportion_healthy_experiment.pdf", width=8, height=5)
#PlotEvolution of Expected medical expenditure across i under Advantageous selection
M_plot = function(i) M(w0=w0, w1=w1,i)
ggplot(data.frame(x=c(N-1,N)), aes(x=x)) + 
  stat_function(fun=M_plot, geom="line", aes(colour = "M")) + 
  xlab("i") +  ylab("")
ggsave(file="expected_medical_expenditure_experiment.pdf", width=8, height=5)
#Plot Labor average productivity
gamma_prod_bar_0_plot = function(i) gamma_prod_bar_0(w0=w0, w1=w1,i)
gamma_prod_bar_1_plot = function(i) gamma_prod_bar_1(w0=w0, w1=w1,i)
ggplot(data.frame(x=c(N-1,N)), aes(x=x)) + 
  stat_function(fun=gamma_prod_bar_0_plot, geom="line",  aes(colour = "gamma_bar0")) + 
  xlab("i") +  ylab("") + stat_function(fun=gamma_prod_bar_1_plot, geom="line",
                aes(colour = "gamma_bar1"))
ggsave(file="average_labor_productivity_experiment.pdf", width=8, height=5)
#Plot effective wages and prices
#Be careful here, for some wages the effective wages wont be well defined, 
#because the endogenous proportion is computed to be the equilibrium one,
#and depends on aggregate Labor suply, that can be corner, thus, the proportion 
#is not well defined if supply is 0 for example.
w_hat0_plot = function(i) w_hat0(w0=w0, w1=w1,i)
w_hat1_plot = function(i) w_hat1(w0=w0, w1=w1,i)
R_hat_plot = function(i) R_hat(R=R,i)
X = X_tilde(w0,w1,Y)
I0 = I_tilde0(w0,w1,R,Y)
I1 = I_tilde1(w0,w1,R,Y)
ggplot(data.frame(x=c(N-1,N)), aes(x=x)) + 
  stat_function(fun=w_hat0_plot, geom="line",  aes(colour = "w_hat0")) + 
  xlab("i") +  ylab("") + stat_function(fun=w_hat1_plot, geom="line",
                                         aes(colour = "w_hat1")) +
  stat_function(fun=R_hat_plot, geom="line", aes(colour = "R_hat")) +
  geom_vline(xintercept = X,linetype=4, colour="black") +
  geom_vline(xintercept = I0,linetype=3, colour="black") +
  geom_vline(xintercept = I1,linetype=2, colour="black") +
  geom_text(mapping = aes(label = "X", y = 0, x = X+0.02),colour="blue") +
  geom_text(mapping = aes(label = "I0", y = 0, x = I0-0.02),colour="blue") +
  geom_text(mapping = aes(label = "I1", y = 0, x = I1+0.02),colour="blue") +
  ggtitle(paste("(w0,w1,R,Y) = (",round(w0,2),",",round(w1,2),",",round(R,2),",",round(Y,2),")"))
ggsave(file="effective_wages_experiment.pdf", width=8, height=5)
#Plot conditional labor demanded and capital
#If the medical expenditure is too big, then for health insurance, seems almost 
#like flat, although it is increasing, showing that Proposition 8 and 9 hold
#k should be flat if z_prod(i)=constant.
l_0d_plot = function(i) l_hat0(w0=w0,w1=w1,i,Y=10)/gamma_prod_bar_0(w0=w0,w1=w1,i)
l_1d_plot = function(i) l_hat1(w0=w0,w1=w1,i,Y=10)/gamma_prod_bar_1(w0=w0,w1=w1,i)
k_plot = function(i) k(R=1,i,Y=10)
ggplot(data.frame(x=c(N-1,N)), aes(x=x)) +  xlab("x") +  ylab("y") + 
  stat_function(fun=l_0d_plot, geom="line",  aes(colour = "l_0d")) + 
  stat_function(fun=l_1d_plot, geom="line",  aes(colour = "l_1d")) +
  stat_function(fun=k_plot, geom="line",  aes(colour = "k"))
#Plot market clearing for Y
#Just testing
Y_seq = seq(from=0, to=20, by=1)
excess_Y_vec = vector(length=length(Y_seq))
for(j in Y_seq){
  excess_Y_vec[j] = fun(j)
}
plot(x=Y_seq, y=excess_Y_vec)
  

