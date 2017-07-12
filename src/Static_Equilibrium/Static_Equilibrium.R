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
lambda_d = 5        #Parameter in sorting function
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
  aux = A_0*exp(A*i)
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
#Normalizing constant in production function
B = function(i){
  zetai = zeta_elas(i)
  aux = (1-eta)^(zetai/(1-zetai))
  return(aux)
}
#Automation Fixed Cost
C_A = function(i){
  aux = 0.1*exp(D*i)
  #aux = 0
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
  else{aux = uniroot(fun, c(initial,final), tol = 1e-13, extendInt = "downX")$root}  #Get the root, 
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
  aux =  exp(2*lambda_d*i - 5*alpha_d)/(1+exp(2*lambda_d*i - 5*alpha_d))*(L1_s('g',w0,w1)/(L1_s('g',w0,w1) + L1_s('b',w0,w1)))
  return(aux)
}
#Expected expenditure shock
E_m = function(h){
  integrand_exp = function(m){ #Creating a function that a returns a vectorized integrand
    aux = vector(length = length(m))
    if(h == 'g'){ #If healthy worker
      for(i in 1:length(m)){
        aux[i] = m[i]*h_g(m[i]) #integrand of expectation healty, for each medical shock
      }
    }
    else{
      for(i in 1:length(m)){
        aux[i] = m[i]*h_b(m[i]) #integrand of expectation unhealty, for each medical shock
      }
    }
    return(aux)
  }
  integral = integrate(integrand_exp, lower = m_L, upper = m_F)
  return(integral$value) #return just the value of the inetgral
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
#TODO: Edit this function later (has to be consistent with the document)
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
  p_0i = function(i) p_0(w0=w0,w1=w1,i)
  p_1i = function(i) p_1(w0=w0,w1=w1,i)
  w_hat0i = function(i) w_hat0(w0=w0,w1=w1,i)
  w_hat1i = function(i) w_hat1(w0=w0,w1=w1,i)
  RHS1 = function(i) (((eta*p_1i(i)^(zeta_elas(i)-1)+(1-eta))^(zeta_elas(i)/(zeta_elas(i)-1)))/(w_hat1i(i)+psi*p_1i(i)^zeta_elas(i)))^(sigma-1)
  RHS0 = function(i) (((eta*p_0i(i)^(zeta_elas(i)-1)+(1-eta))^(zeta_elas(i)/(zeta_elas(i)-1)))/(w_hat0i(i)+psi*p_0i(i)^zeta_elas(i)))^(sigma-1)
  
  fun = function (i) LHS(i)-RHS1(i)+RHS0(i)   #This should be a  decreasing   function of i (but depends on the me
  #dical expenditure)
  initial = N-1 #Lower bound for i
  final = N #Upper boud for i
  #TODO: include boundary cases later
  aux = uniroot(fun, c(initial,final), tol = 1e-13)$root
  return(aux) 
}
#Indifference threshold between capital and not insurance
I_tilde0 = function(w0,w1,R,Y){
  LHS = function(i) C_A(i)/(((B(i)*(sigma-1)/sigma)^(sigma-1))*Y/sigma)
  p_0i = function(i) p_0(w0=w0,w1=w1,i)
  p_ki = function(i) p_k(R=R,i)
  w_hat0i = function(i) w_hat0(w0=w0,w1=w1,i)
  R_hati = function(i) R_hat(R=R,i)
  RHSk = function(i) (((eta*p_ki(i)^(zeta_elas(i)-1)+(1-eta))^(zeta_elas(i)/(zeta_elas(i)-1)))/(R_hati(i)+psi*p_ki(i)^zeta_elas(i)))^(sigma-1)
  RHS0 = function(i) (((eta*p_0i(i)^(zeta_elas(i)-1)+(1-eta))^(zeta_elas(i)/(zeta_elas(i)-1)))/(w_hat0i(i)+psi*p_0i(i)^zeta_elas(i)))^(sigma-1)
  
  fun = function (i) LHS(i)-RHSk(i)+RHS0(i) 
  initial = N-1 #Lower bound for i
  final = N #Upper boud for i
  #TODO: include boundary cases later
  aux = uniroot(fun, c(initial,final), tol = 1e-5)$root
  return(aux) 
}
#Indifference threshold between capital and insurance
I_tilde1 = function(w0,w1,R,Y){
  LHS = function(i) (C_IN-C_A(i))/(((B(i)*(sigma-1)/sigma)^(sigma-1))*Y/sigma)
  p_1i = function(i) p_1(w0=w0,w1=w1,i)
  p_ki = function(i) p_k(R=R,i)
  w_hat1i = function(i) w_hat1(w0=w0,w1=w1,i)
  R_hati = function(i) R_hat(R=R,i)
  RHSk = function(i) (((eta*p_ki(i)^(zeta_elas(i)-1)+(1-eta))^(zeta_elas(i)/(zeta_elas(i)-1)))/(R_hati(i)+psi*p_ki(i)^zeta_elas(i)))^(sigma-1)
  RHS1 = function(i) (((eta*p_1i(i)^(zeta_elas(i)-1)+(1-eta))^(zeta_elas(i)/(zeta_elas(i)-1)))/(w_hat1i(i)+psi*p_1i(i)^zeta_elas(i)))^(sigma-1)
  
  fun = function (i) LHS(i)-RHS1(i)+RHSk(i) #Different sign than before (consistent with the document)
  initial = N-1 #Lower bound for i
  final = N #Upper boud for i
  #TODO: include boundary cases later
  aux = uniroot(fun, c(initial,final), tol = 1e-5)$root
  return(aux) 
}
# NOTES --------------------------------------------------------------------
#How do we specify the reservation utility problem? IN the notse gives 0, but 
#has to be bigger than that.
#We can get advantageous selection for some wages and Adverse selection for others
#if we take (w0,w1)=(2,0.25) we get Advantageous Sel, and if (w0,w1)=(2,1) we get 
#Adverse Sel The question is what happens in equilibrium!
#We can calculate the expectations analitically too to be faster
#Should the Normalizing constant B depend on i too ? It is in the code nut not in the document
#One problem of the code is that whenever wages are such that there is no Labor supply
#for one of the contracts, then the problem of the firm for that contract is not 
#well defined, because Chi is not well defined.
#Effective wages do not cross! It seems that regardless of rho, effective wage approaches
#assymptotically to 0 and they never cross, just get closer and closer (across i).


# Plots -------------------------------------------------------------------
dir = '~/Dropbox/Technology/Codes/Technology_Health/plots/'
setwd(dir)
#Wages to play
w0 = 2
w1 = 1
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
R_hat_plot = function(i) R_hat(R = 1.5,i)
ggplot(data.frame(x=c(N-1,N)), aes(x=x)) + 
  stat_function(fun=w_hat0_plot, geom="line",  aes(colour = "w_hat0")) + 
  xlab("i") +  ylab("") + stat_function(fun=w_hat1_plot, geom="line",
                                         aes(colour = "w_hat1")) +
  ggtitle("(w0,w1)=(2,1)")
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


