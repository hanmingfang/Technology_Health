# Household Equations ---------------------------------------------------
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
#New Theta bar threshold using MGF
theta_ins_MGF = function(h,w0,w1){
  Delta_w = w0-w1
  if(h == 'g'){
    rate_h = rate_g
    P_0h = P_0g
  }
  else{
    rate_h = rate_b
    P_0h = P_0b
  }
  fun = function(theta) rate_h/(rate_h-theta) - (exp(theta*(Delta_w))-P_0h)/(1-P_0h)
  initial = theta_L + 1e-10            #strictly greater than 0 (this number is arbitrary though)
  final = rate_h                      #Due to the MGF, theta < rate_h
  if(Delta_w <= (1-P_0h)/rate_h | fun(initial)>0){
    aux = 0
  }
  else{
    aux = uniroot(fun, c(initial,final), tol = tol, f.upper = Inf)$root
  }
  return(aux)
}
theta_ins = theta_ins_MGF
#Aggregate labor supply for no insurance
L0_s = function(h,s,w0,w1){
  if(h == 'g' & s == sH){aux = lambda_gH*l0_s(w0)*F_gH(theta_ins(h,w0,w1))} # L^0_gH
  else if(h == 'g' & s == sL){aux = lambda_gL*l0_s(w0)*F_gL(theta_ins(h,w0,w1))} # L^0_gL
  else if(h == 'b' & s == sH){aux = lambda_bH*l0_s(w0)*F_bH(theta_ins(h,w0,w1))} # L^0_bH
  else{aux = lambda_bL*l0_s(w0)*F_bL(theta_ins(h,w0,w1))} # L^0_bL
  return(aux)
}
#Aggregate labor supply for no insurance with memoise, this improves the speed if the same value is computed
L0_s_memo = memoise(L0_s)
#Aggregate labor supply for insurance
L1_s = function(h,s,w0,w1){
  if(h == 'g' & s == sH){aux = lambda_gH*l1_s(w1)*(1-F_gH(theta_ins(h,w0,w1)))} # L^1_gH
  else if(h == 'g' & s == sL){aux = lambda_gL*l1_s(w1)*(1-F_gL(theta_ins(h,w0,w1)))} # L^1_gL
  else if(h == 'b' & s == sH){aux = lambda_bH*l1_s(w1)*(1-F_bH(theta_ins(h,w0,w1)))} # L^1_bH
  else{aux = lambda_bL*l1_s(w1)*(1-F_bL(theta_ins(h,w0,w1)))} # L^1_bL
  return(aux)
}
#Aggregate labor supply for insurance with memoise, this improves the speed if the same value is computed
L1_s_memo = memoise(L1_s) 



