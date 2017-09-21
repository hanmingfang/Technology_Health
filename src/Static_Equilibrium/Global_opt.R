
#Objective function for Genetic algorithm
#Take away the exponentials in the excess demand functions to use it
obj_fun_prob = function(p){  
  Y = p[6]
  aux  = sqrt((k_excess_d_prob(p)/K)^2 + (l0H_excess_d_prob(p))^2  + (l0L_excess_d_prob(p))^2+
                (l1H_excess_d_prob(p))^2 + (l1L_excess_d_prob(p))^2 + (Y_excess_s_prob(p)/Y)^2)
  return(aux)#Here I need to normalize in some way the Excess demands
}


#Global optimizer with CRS2 (is slow but quite robust)
#Take away exponentials in the excess demand functions if you are using this method
ptm = proc.time()
crs2_sol = crs2lm(x0 = c(log(300),log(5.4),log(300),log(3.3),log(4.6),log(600)),
                  fn = obj_fun_prob,
                  lower = c(log(0.001),log(0.001),log(0.001),log(0.001),log(0.001),log(0.001)),
                  upper = c(log(1000),log(1000),log(1000),log(1000),log(1000),log(10000)),
                  maxeval = 10000,
                  xtol_rel = 1e-6)
proc.time() - ptm
crs2_sol
crs2_sol$par
w0H = exp(crs2_sol$par[1])
w0L = exp(crs2_sol$par[2])
w1H = exp(crs2_sol$par[3])
w1L = exp(crs2_sol$par[4])
R = exp(crs2_sol$par[5])
Y = exp(crs2_sol$par[6])
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

