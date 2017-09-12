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
#Threshold for household insurance decision
#TODO: include \n in the message
#Notice that as we have assumed a common domain for theta regardless of the skill type
#We do not need to index this function by skill, and will depend just on wages and health
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


