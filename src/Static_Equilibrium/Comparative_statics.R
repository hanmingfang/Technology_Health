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
s_ccp = 1           #Parameter to scale ccps 
tol = 1e-12          #Tolerance for unitroot, affects computation time
K = 1               #Capital stock in the economy
n_nodes = 40        #Number of nodes to evaluate numerical integrals
# Equilibrium solutions -------------------------------------------------------------------

#Sourcing equations from the model
source("src/Static_Equilibrium/Primitive_fcns.R")
source("src/Static_Equilibrium/Household_fcns.R")
source("src/Static_Equilibrium/Firm_fcns.R")
source("src/Static_Equilibrium/Equilibrium_fcns.R")

#Creating list of objects
n_grid = 20
n_grid = length(grid_par)
w0H = vector(length = n_grid)
w0L = vector(length = n_grid)
w1H = vector(length = n_grid)
w1L = vector(length = n_grid)
R = vector(length = n_grid)
Y = vector(length = n_grid)
L1_gH = vector(length = n_grid)
L1_bH = vector(length = n_grid)
L1_H = vector(length = n_grid)
L1_gL = vector(length = n_grid)
L1_bL = vector(length = n_grid)
L1_L = vector(length = n_grid)
w1 = vector(length = n_grid)
L0_gH = vector(length = n_grid)
L0_bH = vector(length = n_grid)
L0_H = vector(length = n_grid)
L0_gL = vector(length = n_grid)
L0_bL = vector(length = n_grid)
L0_L = vector(length = n_grid)
w0 = vector(length = n_grid)
m_k = vector(length = n_grid)
m_l0H = vector(length = n_grid)
m_l0L = vector(length = n_grid)
m_l1H = vector(length = n_grid)
m_l1L = vector(length = n_grid)
nles_sol_vec = vector(length = n_grid)
data_eq = list()
#Param changes
param_changes =  c('D', 'lambda_d','rho', 'delta_H','lambda_H','lambda_g')
#grid parameter values
grid_par = list()
grid_par[['D']] = seq(from = 1, by = 1, length.out = n_grid)
grid_par[['lambda_d']] = seq(from = 10, by = 1, length.out = n_grid)
grid_par[['rho']] = seq(from = 0.8, by = 0.01, length.out = n_grid)
grid_par[['delta_H']] = seq(from = 4, by = 0.1, length.out = n_grid)
grid_par[['lambda_H']] = seq(from = 0.25, by = 0.005, length.out = n_grid)
grid_par[['lambda_g']] = seq(from = 0.25, by = 0.005, length.out = n_grid)
#Computing the equilibrium for different values
#Start with a very close initial condition, and then use 
#the previous equilibrium solution as initial condition
x_initial = c(log(120),log(5.1),log(109),log(2.8),log(4.3),log(258))
#Start clock
ptm = proc.time()
for(j in param_changes){
  lambda_gH = 0.25    #Measure healthy workers High skill
  lambda_gL = 0.25    #Measure healthy workers Low skill
  lambda_bH = 0.25    #Measure unhealthy workers High skill
  lambda_bL = 0.25    #Measure unhealthy workers Low skill
  rho = 0.8           #Relative labor productivity of unhealthy workers
  delta_H =  4        #Parameter in labor productivity of High skill type
  lambda_d = 10       #Parameter in sorting function
  D = 1               #Parameter in Automation Cost function
  for(i in 1:n_grid){
    #Changing parameter value depending on the simulation
    if(j =='D'){ D = grid_par[[j]][i]}
    else if(j =='lambda_d'){ lambda_d = grid_par[[j]][i]}
    else if(j =='delta_H'){ delta_H = grid_par[[j]][i]}
    else if(j =='lambda_H'){ 
      lambda_gH = grid_par[[j]][i]
      lambda_bH = grid_par[[j]][i]
      lambda_gL =  0.25 - (grid_par[[j]][i]-0.25)
      lambda_bL =  0.25 - (grid_par[[j]][i]-0.25)
      print(paste("lambda_H = ", lambda_gH + lambda_bH, "lambda_L = ", lambda_gL + lambda_bL))
    }
    else if(j =='lambda_g'){ 
      lambda_gH = grid_par[[j]][i]
      lambda_bH = 0.25 - (grid_par[[j]][i]-0.25)
      lambda_gL = grid_par[[j]][i]
      lambda_bL = 0.25 - (grid_par[[j]][i]-0.25)
      print(paste("lambda_g = ", lambda_gH + lambda_gL,"lambda_b = ", lambda_bH + lambda_bL))
    }
    else{ rho = grid_par[[j]][i]}
    #Sourcing functions again
    source("src/Static_Equilibrium/Primitive_fcns.R")
    source("src/Static_Equilibrium/Household_fcns.R")
    source("src/Static_Equilibrium/Firm_fcns.R")
    source("src/Static_Equilibrium/Equilibrium_fcns.R")
    #TODO: Check why fails the try at Even iterations, there is a pattern (mistake in the code)
    nles_sol = tryCatch(
      { if(i>1){
          x_initial = c(log(w0H[i-1]),log(w0L[i-1]),log(w1H[i-1]),log(w1L[i-1]),log(R[i-1]),log(Y[i-1]))
        }
        print(paste("Iteration ",i))
        #The try will return the last assignement part
        nles_sol = nleqslv(x = x_initial, 
                           fn = F_zeros, jac=NULL, 
                           method = "Broyden", 
                           jacobian=FALSE, 
                           xscalm = "fixed",
                           control = list("allowSingular"=TRUE, scalex = c(0.01,1,0.01,1,1,0.01),
                                          trace = 1, btol=.001, xtol = 1e-10, maxit = 100),
                           global="dbldog")
      },
      error = function(cond){
        x_initial = c(log(120),log(5.1),log(109),log(2.8),log(4.3),log(258))
        message(paste(cond,"Taking starting initial condition"))
        nles_sol = nleqslv(x = x_initial, 
                           fn = F_zeros, jac=NULL, 
                           method = "Broyden", 
                           jacobian=FALSE, 
                           xscalm = "fixed",
                           control = list("allowSingular"=TRUE, scalex = c(0.01,1,0.01,1,1,0.01),
                                          trace = 1, btol=.001, xtol = 1e-10, maxit = 100),
                           global="dbldog")
        return(nles_sol) 
      }
    )
    w0H[i] = exp(nles_sol$x[1])
    w0L[i] = exp(nles_sol$x[2])
    w1H[i] = exp(nles_sol$x[3])
    w1L[i] = exp(nles_sol$x[4])
    R[i] = exp(nles_sol$x[5])
    Y[i] = exp(nles_sol$x[6])
    nles_sol_vec[i] = nles_sol$message
    #Checking for wieghted wages 
    L1_gH[i] = L1_s('g',sH,w0H[i],w1H[i])
    L1_bH[i] = L1_s('b',sH,w0H[i],w1H[i])
    L1_H[i] = L1_gH[i] + L1_bH[i]
    L1_gL[i] = L1_s('g',sL,w0L[i],w1L[i])
    L1_bL[i] = L1_s('b',sL,w0L[i],w1L[i])
    L1_L[i] = L1_gL[i] + L1_bL[i]
    w1[i] = (w1H[i]*L1_H[i] + w1L[i]*L1_L[i])/(L1_H[i]+L1_L[i]) 
    L0_gH[i] = L0_s('g',sH,w0H[i],w1H[i])
    L0_bH[i] = L0_s('b',sH,w0H[i],w1H[i])
    L0_H[i] = L0_gH[i] + L0_bH[i]
    L0_gL[i] = L0_s('g',sL,w0L[i],w1L[i])
    L0_bL[i] = L0_s('b',sL,w0L[i],w1L[i])
    L0_L[i] = L0_gL[i] + L0_bL[i]
    w0[i] = (w0H[i]*L0_H[i] + w0L[i]*L0_L[i])/(L0_H[i]+L0_L[i])
    #Calculating measure of firms that demand each input
    m_k[i] = measure_k(w0H[i],w0L[i],w1H[i],w1L[i],R[i],Y[i])
    m_l0H[i] = measure_l0H(w0H[i],w0L[i],w1H[i],w1L[i],R[i],Y[i])
    m_l0L[i] = measure_l0L(w0H[i],w0L[i],w1H[i],w1L[i],R[i],Y[i])
    m_l1H[i] = measure_l1H(w0H[i],w0L[i],w1H[i],w1L[i],R[i],Y[i])
    m_l1L[i] = measure_l1L(w0H[i],w0L[i],w1H[i],w1L[i],R[i],Y[i])
  }
#Saving results in a data frame
  data_eq[[j]] = data.frame(grid_par[[j]], w0H,w0L,w1H,w1L, R,Y, nles_sol_vec, L1_gH, L1_bH, L1_H, L1_gL, L1_bL,
                        L1_L, w1, L0_gH, L0_bH, L0_H, L0_gL, L0_bL, L0_L, w0, m_k, m_l0H, m_l0L, m_l1H, m_l1L)
}
#Stop clock
proc.time() - ptm
#Saving R data
save(data_eq,file="data/data.Rda")














