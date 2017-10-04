#Solve for static equilibrium 
rm(list = ls())
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
library(statmod)

#Set my directory
dir = '~/Technology_Health/'
setwd(dir)

# Parameters --------------------------------------------------------------
#Household 
sH = 4              #Skill High type 
sL = 1              #Skill Low type
lambda_gH = 0.25    #Measure healthy workers High skill
lambda_gL = 0.25    #Measure healthy workers Low skill
lambda_bH = 0.25    #Measure unhealthy workers High skill
lambda_bL = 0.25    #Measure unhealthy workers Low skill
theta_L = 0         #Domain for theta
theta_H = Inf
#m_L = 0             #Domain for medical expenditure shocks
#m_F = Inf
shape_gH = 1         #Shape parameter of theta (Gamma distribution) High skill
shape_gL = 1         #Shape parameter of theta (Gamma distribution) Low skill
shape_bH = 1       #Shape parameter of theta (Gamma distribution) High skill
shape_bL = 1       #Shape parameter of theta (Gamma distribution) Low skill
#For a given scale parameter, higher shape parameter means more risk averse households
scale_gH = 1         #Scale parameter of theta (Gamma distribution) High skill
scale_gL = 1         #Scale parameter of theta (Gamma distribution) Low skill
scale_bH = 1         #Scale parameter of theta (Gamma distribution) High skill
scale_bL = 1         #Scale parameter of theta (Gamma distribution) Low skill
rate_g = 1.5           #Rate for exponential distribution for medical exp. healthy workers  (Mean=1/rate)
rate_b = 0.25        #Rate for exponential distribution for medical exp. unhealthy workers
M_trunc = 150        #Truncation for TEXP medical expenditure
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
delta_H =  4        #Parameter in labor productivity of High skill type
lambda_d = 10       #Parameter in sorting function
#alpha_d = 5         #Parameter in sorting function
D = 1               #Parameter in Automation Cost function
D_0 = 1               #Parameter in Automation Cost function
s_ccp = 1           #Parameter to scale ccps 
tol = 1e-12          #Tolerance for unitroot, affects computation time
K = 1               #Capital stock in the economy
n_nodes = 20        #Number of nodes to evaluate numerical integrals

# Equilibrium solutions -------------------------------------------------------------------

#Sourcing equations from the model
source("src/Static_Equilibrium/Primitive_fcns.R")
source("src/Static_Equilibrium/Household_fcns.R")
source("src/Static_Equilibrium/Firm_fcns.R")
source("src/Static_Equilibrium/Equilibrium_fcns.R")

#Non linear equation solver (is very fast but a little more sensitive to initial conditions)
#Using exponentials in the equations to ensure positivity
#Change method to "Newton" if "Broyden" doesn't converge. However Broyden is way faster
#If the algorithm doesn't find a better point, try decreasing C_IN
ptm = proc.time()
nles_sol = nleqslv(x = c(log(120),log(5),log(109),log(2.7),log(4.2),log(258)), 
                   fn = F_zeros, jac=NULL, method = "Broyden", jacobian=FALSE, xscalm = "fixed",
                   control = list("allowSingular"=TRUE, scalex = c(0.01,1,0.01,1,1,0.01), trace = 1, btol=.001),
                   global="dbldog")
proc.time() - ptm
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
L1_gH = L1_s('g',sH,w0H,w1H)
L1_bH = L1_s('b',sH,w0H,w1H)
L1_H = L1_gH + L1_bH
L1_gL = L1_s('g',sL,w0L,w1L)
L1_bL = L1_s('b',sL,w0L,w1L)
L1_L = L1_gL + L1_bL
w1 = (w1H*L1_H + w1L*L1_L)/(L1_H+L1_L) 

L0_gH = L0_s('g',sH,w0H,w1H)
L0_bH = L0_s('b',sH,w0H,w1H)
L0_H = L0_gH + L0_bH
L0_gL = L0_s('g',sL,w0L,w1L)
L0_bL = L0_s('b',sL,w0L,w1L)
L0_L = L0_gL + L0_bL
w0 = (w0H*L0_H + w0L*L0_L)/(L0_H+L0_L) 
c(w0,w1)
#Checking if w1>w0
w0-w1
#Checking if the compensation package s higher (w1+P)>w0
#PremiumH = 
#PremiumL =   
#w1P = ((w1H+PremiumH)*L1_H + (w1L+PremiumL)*L1_L)/(L1_H+L1_L)
#To test different algorithms uncomment next line
#ptm = proc.time()
#testnslv(x =c(log(10),log(2),log(8),log(0.2),log(1.3),log(28)), fn = F_zeros)
#proc.time() - ptm


